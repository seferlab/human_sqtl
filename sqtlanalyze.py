import numpy as np
import scipy as sp
import random
import ensembl_rest
from Ensembl_converter import EnsemblConverter
import pandas as pd
from scipy.optimize import linear_sum_assignment # For bipartite matching (min cost)
from scipy.stats import wilcoxon # For Wilcoxon signed-rank test
import math
    
def construct_sample_space(gene_lengths, distances, c, Q_tilde, n_p):
    Omega_Q_tilde = []

    while len(Omega_Q_tilde) <= len(Q_tilde) * n_p:
        # Step 2: Choose gene length l
        l = random.choice(gene_lengths)

        # Step 3: Choose a distance d
        d = random.choice(distances)

        # Step 4: Choose direction (True = downstream, False = upstream)
        downstream = random.choice([True, False])

        # Step 5–9: Choose SNP position s
        if downstream:
            s = random.randint(l + d, c)
        else:
            s = random.randint(1, c - l - d)

        # Step 10: Construct a random sQTL q' (dummy representation for now)
        # This part depends on the specific format of an sQTL; let's assume it's a tuple (l, d, s, direction)
        q_prime = (l, d, s, 'downstream' if downstream else 'upstream')

        # Step 11: Append to Omega_Q_tilde
        Omega_Q_tilde.append(q_prime)

    return Omega_Q_tilde


def construct_tad_sample_space(tagged_intervals):
    # Step 1–2: Create O as list of (tag, length) pairs
    O = list(tagged_intervals)
    
    # Step 3: Separate into domain and non-domain collections
    O_domain = [length for tag, length in O if tag == 'domain']
    O_non = [length for tag, length in O if tag == 'non-domain']
    
    # Step 4: Shuffle both collections
    random.shuffle(O_domain)
    random.shuffle(O_non)
    
    # Step 5: Reconstruct a randomized O with same structure (tags preserved, lengths shuffled)
    L = []
    for tag, _ in O:
        if tag == 'domain':
            length = O_domain.pop()  # Get next shuffled domain length
        else:
            length = O_non.pop()  # Get next shuffled non-domain length
        L.append((tag, length))
    
    # Step 6: Assign starting positions based on cumulative sum of lengths
    tad_positions = []
    start = 0
    for tag, length in L:
        if tag == 'domain':
            tad_positions.append((start, start + length))  # Domain start and end
        start += length

    return tad_positions  # List of (start, end) positions for each domain


  
# --- 1. Data Structures (Conceptual) ---
# In a real scenario, these would be loaded from files (e.g., BED, TSV)
# and processed into appropriate data structures (e.g., pandas DataFrames, custom objects).

class HiCFragment:
    """Represents a 40kb Hi-C genomic fragment."""
    def __init__(self, chrom, start, end, fragment_id):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.fragment_id = fragment_id # Unique ID for the fragment

class Gene:
    """Represents a gene with its genomic location and associated Hi-C fragments."""
    def __init__(self, gene_name, chrom, start, end, fragments):
        self.gene_name = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.fragments = fragments # List of HiCFragment objects overlapping the gene

class SNP:
    """Represents a SNP with its genomic location and associated Hi-C fragment."""
    def __init__(self, snp_id, chrom, position, fragment):
        self.snp_id = snp_id
        self.chrom = chrom
        self.position = position
        self.fragment = fragment # HiCFragment object where the SNP maps

class EQTLGenePair:
    """Represents an eQTL-gene pair before equivalence class formation."""
    def __init__(self, snp, gene):
        self.snp = snp
        self.gene = gene

class EQTLGeneEquivalenceClass:
    """
    Represents an eQTL-gene equivalence class.
    All pairs in this class map to the same SNP Hi-C fragment and target gene.
    """
    def __init__(self, snp_fragment, gene_fragments, original_pairs):
        self.snp_fragment = snp_fragment # HiCFragment object for the SNP
        self.gene_fragments = gene_fragments # List of HiCFragment objects for the gene
        self.original_pairs = original_pairs # List of original EQTLGenePair objects
        self.chrom = snp_fragment.chrom # Chromosome of this equivalence class

        # Properties (will be calculated)
        self.snp_gene_distance = None
        self.gene_length = None
        self.total_snp_frequency = None
        self.maximum_total_gene_frequency = None
        self.spatial_proximity = None # p(q)

    def calculate_properties(self, hic_interaction_data):
        """
        Calculates the five properties for this eQTL-gene equivalence class.
        hic_interaction_data: A dictionary or DataFrame representing Hi-C interactions.
                              Key: tuple (fragment_id1, fragment_id2), Value: interaction_frequency
        """
        # d(q): SNP-gene distance
        # Midpoint of SNP fragment
        snp_mid = (self.snp_fragment.start + self.snp_fragment.end) / 2

        # Find midpoint of closest gene fragment to SNP
        min_dist = float('inf')
        closest_gene_frag_mid = None
        for gf in self.gene_fragments:
            gene_frag_mid = (gf.start + gf.end) / 2
            dist = abs(snp_mid - gene_frag_mid)
            if dist < min_dist:
                min_dist = dist
                closest_gene_frag_mid = gene_frag_mid
        self.snp_gene_distance = abs(snp_mid - closest_gene_frag_mid) if closest_gene_frag_mid is not None else 0

        # l(q): Gene length
        if self.gene_fragments:
            gene_starts = [f.start for f in self.gene_fragments]
            gene_ends = [f.end for f in self.gene_fragments]
            self.gene_length = max(gene_ends) - min(gene_starts)
        else:
            self.gene_length = 0

        # t_s(q): Total SNP frequency
        # This requires summing interactions for the SNP fragment with all its neighbors.
        # In a real implementation, you'd query hic_interaction_data for all interactions
        # involving self.snp_fragment.fragment_id.
        self.total_snp_frequency = self._get_total_fragment_frequency(self.snp_fragment.fragment_id, hic_interaction_data)

        # t_g(q): Maximum total gene frequency
        # Max t(v) for any fragment v overlapping the gene
        self.maximum_total_gene_frequency = 0
        for gf in self.gene_fragments:
            self.maximum_total_gene_frequency = max(self.maximum_total_gene_frequency,
                                                    self._get_total_fragment_frequency(gf.fragment_id, hic_interaction_data))

        # p(q): Spatial proximity
        # Maximum frequency f'(e) between SNP fragment and gene fragments
        self.spatial_proximity = 0
        for gf in self.gene_fragments:
            interaction_key = tuple(sorted((self.snp_fragment.fragment_id, gf.fragment_id)))
            if interaction_key in hic_interaction_data:
                self.spatial_proximity = max(self.spatial_proximity, hic_interaction_data[interaction_key])

    def _get_total_fragment_frequency(self, fragment_id, hic_interaction_data):
        """Helper to sum interaction frequencies for a given fragment."""
        total_freq = 0
        for (f1, f2), freq in hic_interaction_data.items():
            if f1 == fragment_id or f2 == fragment_id:
                total_freq += freq
        return total_freq

    def get_feature_vector(self):
        """Returns the feature vector for optimal matching."""
        return np.array([
            self.snp_gene_distance,
            self.gene_length,
            self.total_snp_frequency,
            self.maximum_total_gene_frequency
        ])

class TopologicalDomain:
    """Represents a topological domain."""
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = end - start
       
def load_data(tissue, resol):
    """
    Loads data for Hi-C fragments, genes, SNPs, eQTLs, and domains.
    """
    #Hi-C fragments (40kb resolution)
    fragments = []
    chro2len = {}
    df_temp = pd.read_csv("human_chromosome_lengths.csv", sep="\t")
    for item in df_temp["Chromosome,Length (bp)"]:
        chro, length = item.split(",")
        if chro in ['chrX', 'chrY','chrM']:
            continue
        
        length = int(length)
        chro2len[chro] = length
    chro2fragments = {}
    for chro, length in chro2len.items():
        block = int(math.ceil(length/resol))
        chro2fragments[chro] = []
        for i in range(block): # 100 fragments on chromosome 1
            fragments.append(HiCFragment(chro, i * resol, (i + 1) * resol, f'frag_{i}'))
            chro2fragments[chro].append(HiCFragment(chro, i * resol, (i + 1) * resol, f'frag_{i}'))
    
    # Hi-C interaction data (fragment_id1, fragment_id2): frequency
    #hic_interactions = {}
    #hicpath = "{0}.mat".format(tissue)
    #row = 0
    #with open(hicpath, "r") as infile:
    #    for line in infile:
    #        splitted = [int(item) for item in line.split("\t")]
    #        for column,item2 in enumerate(splitted):
    #            hic_interactions[tuple(sorted((f'frag_{row}', f'frag_{column}')))] = random.randint(10, 100)

    # Genes (overlapping fragments)
    outfname = "dataset/sqtls/{0}.csv".format(tissue)
    df = pd.read_csv(outfname)
    seen = set()
    genes = []
    snps = []
    sqtl_gene_pairs = []
    for rind,row in df.iterrows():
        chro = row["Gene chromosome"]
        start = int(row["Gene start"] / resol)
        end = int(row["Gene end"] / resol)
        if chro in ["X", "Y"]:
            continue
        sub = [chro2fragments["chr{0}".format(chro)][ind] for ind in range(start,end+1)]
        gene = Gene(row["Gene"], 'chr{0}'.format(chro), row["Gene start"], row["Gene end"], sub)

        name = "snp{0}".format(rind+1)
        chro = row["Sqtl chromosome"]
        if chro in ["X", "Y"]:
            continue
        pos = row["Sqtl position"]
        frag = int(pos / resol)
        snps.append(SNP(name, 'chr{0}'.format(chro), pos, chro2fragments["chr{0}".format(chro)][frag]))
        sqtl_gene_pairs.append(EQTLGenePair(snps[-1], gene))
        
        if row["Gene"] not in seen:
            genes.append(Gene(row["Gene"], 'chr{0}'.format(chro), row["Gene start"], row["Gene end"], sub))

    # Topological Domains
    #domains = []
    #domainpath = "{0}.tad".format(tissue)
    #with open(domainpath, "r") as infile:
    #    for line in infile:
    #        chro,start,end=line.split("\t")
    #        domains.append(TopologicalDomain('chr{0}'.format(chro), start, end))

    hic_interactions = {}
    domains = []
    return fragments, hic_interactions, genes, snps, sqtl_gene_pairs, domains


# --- 3. Core Analysis Functions ---

def create_eqtl_equivalence_classes(eqtl_gene_pairs):
    """
    Aggregates eQTL-gene pairs into equivalence classes.
    An equivalence class is defined by the same SNP Hi-C fragment and target gene.
    """
    eq_classes = {}
    for pair in eqtl_gene_pairs:
        snp_frag_id = pair.snp.fragment.fragment_id
        gene_name = pair.gene.gene_name
        key = (snp_frag_id, gene_name)

        if key not in eq_classes:
            eq_classes[key] = EQTLGeneEquivalenceClass(
                pair.snp.fragment,
                pair.gene.fragments,
                [] # Initialize with empty list of original pairs
            )
        eq_classes[key].original_pairs.append(pair)
    return list(eq_classes.values())

def generate_random_eqtl_equivalence_classes(
    observed_eq_classes,
    all_hic_fragments,
    all_genes,
    hic_interaction_data,
    num_random_multiplier=100, # n_p from paper
    chromosome_length=4000000 # Example chromosome length
):
    """
    Generates a sample space of random eQTL equivalence classes.
    This is a simplified version of the procedure described in the paper.
    """
    random_eq_classes = []
    target_size = len(observed_eq_classes) * num_random_multiplier

    # Collect observed gene lengths and SNP-gene distances
    observed_gene_lengths = [eq.gene_length for eq in observed_eq_classes if eq.gene_length is not None]
    observed_snp_gene_distances = [eq.snp_gene_distance for eq in observed_eq_classes if eq.snp_gene_distance is not None]

    # Calculate observed probability of SNP upstream/downstream
    upstream_count = 0
    downstream_count = 0
    for eq in observed_eq_classes:
        snp_mid = (eq.snp_fragment.start + eq.snp_fragment.end) / 2
        gene_mid = (eq.gene_fragments[0].start + eq.gene_fragments[-1].end) / 2 if eq.gene_fragments else snp_mid
        if snp_mid < gene_mid:
            upstream_count += 1
        elif snp_mid > gene_mid:
            downstream_count += 1
    total_oriented = upstream_count + downstream_count
    prob_upstream = upstream_count / total_oriented if total_oriented > 0 else 0

    iteration_count = 0
    while len(random_eq_classes) < target_size and iteration_count < 10000000: # Max iterations
        iteration_count += 1

        # 1. Select a random gene length
        if not observed_gene_lengths: continue
        rand_gene_length = random.choice(observed_gene_lengths)

        # 2. Select a random SNP-gene distance
        if not observed_snp_gene_distances: continue
        rand_snp_gene_distance = random.choice(observed_snp_gene_distances)

        # 3. Randomly decide upstream/downstream
        is_upstream = random.random() < prob_upstream

        # 4. Build a random eQTL q' and determine its equivalence class
        # Simplified gene and SNP fragment selection for mock data
        # In a real scenario, you'd need to ensure valid fragment ranges for gene/snp.
        if not all_genes or not all_hic_fragments: continue

        rand_gene = random.choice(all_genes)
        # Ensure gene has fragments for properties calculation
        if not rand_gene.fragments: continue

        # Determine SNP location relative to gene
        gene_start = rand_gene.start
        gene_end = rand_gene.end

        # Calculate potential SNP fragment range
        if is_upstream:
            # SNP is upstream: [1, gene_start - (rand_gene_length + rand_snp_gene_distance)]
            # Simplified: just pick a fragment before the gene
            possible_snp_fragments = [f for f in all_hic_fragments if f.end < gene_start]
        else:
            # SNP is downstream: [gene_end + (rand_gene_length + rand_snp_gene_distance), chromosome_length]
            # Simplified: just pick a fragment after the gene
            possible_snp_fragments = [f for f in all_hic_fragments if f.start > gene_end]

        if not possible_snp_fragments: continue
        rand_snp_fragment = random.choice(possible_snp_fragments)

        # Create a mock SNP and EQTLGenePair for the equivalence class
        mock_snp = SNP(f'rand_snp_{iteration_count}', rand_snp_fragment.chrom,
                       (rand_snp_fragment.start + rand_snp_fragment.end) / 2, rand_snp_fragment)
        mock_pair = EQTLGenePair(mock_snp, rand_gene)

        # Create a temporary equivalence class to calculate properties
        temp_eq_class = EQTLGeneEquivalenceClass(mock_snp.fragment, rand_gene.fragments, [mock_pair])
        temp_eq_class.calculate_properties(hic_interaction_data)

        # Ensure properties are not None before adding
        if all(p is not None for p in temp_eq_class.get_feature_vector()):
            random_eq_classes.append(temp_eq_class)

    print(f"Generated {len(random_eq_classes)} random eQTL equivalence classes in {iteration_count} iterations.")
    return random_eq_classes


def perform_optimal_matching(observed_eq_classes, random_eq_classes, k=1):
    """
    Performs optimal matching using a bipartite matching algorithm.
    Here, we use scipy's linear_sum_assignment for minimum cost.
    For maximum weight, you'd typically convert weights to costs (max_weight - weight).
    """
    if not observed_eq_classes or not random_eq_classes:
        print("Cannot perform matching: observed or random classes are empty.")
        return [], []

    # Create cost matrix (Euclidean distance between feature vectors)
    # Feature vector: [d(q), l(q), ts(q), tg(q)]
    observed_features = np.array([eq.get_feature_vector() for eq in observed_eq_classes])
    random_features = np.array([eq.get_feature_vector() for eq in random_eq_classes])

    num_observed = len(observed_features)
    num_random = len(random_features)

    # Calculate Euclidean distances
    # This creates a matrix where element (i, j) is distance between observed_features[i] and random_features[j]
    cost_matrix = np.zeros((num_observed, num_random))
    for i in range(num_observed):
        for j in range(num_random):
            cost_matrix[i, j] = np.linalg.norm(observed_features[i] - random_features[j])

    # Perform bipartite matching (linear_sum_assignment finds min cost)
    # This will match each observed element to one random element with minimum total cost.
    # For k matches, you'd need a more advanced algorithm or iterative approach.
    # For simplicity, we'll do k=1 here.
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    matched_observed = [observed_eq_classes[i] for i in row_ind]
    matched_random = [random_eq_classes[j] for j in col_ind]

    return matched_observed, matched_random

def statistical_analysis_eqtls(matched_observed, matched_random):
    """
    Performs Wilcoxon signed-rank test on spatial proximity (p(q)).
    """
    if not matched_observed or not matched_random:
        print("No matched pairs for statistical analysis.")
        return

    observed_spatial_proximity = np.array([eq.spatial_proximity for eq in matched_observed])
    random_spatial_proximity = np.array([eq.spatial_proximity for eq in matched_random])

    # Ensure arrays are of the same length for Wilcoxon test
    min_len = min(len(observed_spatial_proximity), len(random_spatial_proximity))
    observed_spatial_proximity = observed_spatial_proximity[:min_len]
    random_spatial_proximity = random_spatial_proximity[:min_len]

    if min_len == 0:
        print("Insufficient data for Wilcoxon test.")
        return

    # Wilcoxon signed-rank test: tests if two related paired samples come from the same distribution.
    # Here, we test if observed spatial proximity is significantly different from random.
    # A one-sided test could be used if we hypothesize 'greater' spatial proximity.
    stat, p_value = wilcoxon(observed_spatial_proximity, random_spatial_proximity, alternative='greater')

    print("\n--- Statistical Analysis of eQTLs ---")
    print(f"Wilcoxon Signed-Rank Test (Observed vs. Random Spatial Proximity):")
    print(f"  Test Statistic: {stat:.4f}")
    print(f"  P-value: {p_value:.4e}") # Scientific notation for small p-values

    if p_value < 0.05:
        print("  Conclusion: Observed eQTLs show significantly greater spatial proximity than random counterparts (p < 0.05).")
    else:
        print("  Conclusion: No significant evidence that observed eQTLs have greater spatial proximity than random counterparts (p >= 0.05).")

def generate_shuffled_domains(observed_domains, chromosome_length, num_shuffles=10000):
    """
    Generates shuffled topological domain sequences.
    This is a simplified conceptual implementation.
    """
    shuffled_domain_sequences = []

    # 1. Label intervals as 'domain' or 'non-domain' and get lengths
    intervals = []
    current_pos = 0
    domain_starts = sorted([d.start for d in observed_domains])
    domain_ends = sorted([d.end for d in observed_domains])

    all_positions = sorted(list(set(domain_starts + domain_ends + [0, chromosome_length])))

    for i in range(len(all_positions) - 1):
        segment_start = all_positions[i]
        segment_end = all_positions[i+1]
        if segment_start == segment_end:
            continue

        # Check if this segment overlaps with any observed domain
        is_domain = False
        for d in observed_domains:
            if max(segment_start, d.start) < min(segment_end, d.end):
                is_domain = True
                break
        intervals.append({'type': 'domain' if is_domain else 'non-domain', 'length': segment_end - segment_start})

    # 2. Create lists of domain and non-domain lengths
    domain_lengths = [i['length'] for i in intervals if i['type'] == 'domain']
    non_domain_lengths = [i['length'] for i in intervals if i['type'] == 'non-domain']
    original_types_order = [i['type'] for i in intervals]

    print(f"Original intervals: {len(intervals)} (Domains: {len(domain_lengths)}, Non-domains: {len(non_domain_lengths)})")

    for _ in range(num_shuffles):
        shuffled_domain_lengths = list(domain_lengths)
        shuffled_non_domain_lengths = list(non_domain_lengths)
        random.shuffle(shuffled_domain_lengths)
        random.shuffle(shuffled_non_domain_lengths)

        current_shuffled_domains = []
        current_pos = 0
        domain_idx = 0
        non_domain_idx = 0

        for segment_type in original_types_order:
            if segment_type == 'domain':
                if domain_idx < len(shuffled_domain_lengths):
                    length = shuffled_domain_lengths[domain_idx]
                    current_shuffled_domains.append(TopologicalDomain(observed_domains[0].chrom, current_pos, current_pos + length))
                    current_pos += length
                    domain_idx += 1
                else:
                    # Handle cases where shuffled lists might run out if original_types_order has more domains/non-domains than available lengths
                    # This implies an issue with interval generation or shuffling logic for very complex cases.
                    # For simplicity, we'll just break or continue, but in a real scenario, this needs robust handling.
                    pass
            else: # non-domain
                if non_domain_idx < len(shuffled_non_domain_lengths):
                    length = shuffled_non_domain_lengths[non_domain_idx]
                    current_pos += length # Non-domains are gaps, just advance position
                    non_domain_idx += 1
                else:
                    pass
        shuffled_domain_sequences.append(current_shuffled_domains)

    print(f"Generated {len(shuffled_domain_sequences)} shuffled domain sequences.")
    return shuffled_domain_sequences

def analyze_domain_crossings(observed_eq_classes, observed_domains, shuffled_domain_sequences):
    """
    Analyzes how often eQTL-gene pairs cross topological domains.
    This is a conceptual outline.
    """
    def count_crossings(eq_class, domains):
        """Counts if an eQTL-gene interval crosses a domain."""
        snp_mid = (eq_class.snp_fragment.start + eq_class.snp_fragment.end) / 2
        gene_mid = (eq_class.gene_fragments[0].start + eq_class.gene_fragments[-1].end) / 2 if eq_class.gene_fragments else snp_mid
        interval_start = min(snp_mid, gene_mid)
        interval_end = max(snp_mid, gene_mid)

        for d in domains:
            # Check for overlap
            overlap_start = max(interval_start, d.start)
            overlap_end = min(interval_end, d.end)

            if overlap_start < overlap_end: # There is an overlap
                # If the interval is NOT completely contained within the domain, it crosses
                if not (interval_start >= d.start and interval_end <= d.end):
                    return 1 # Crosses at least one domain
        return 0

    observed_crossings = sum(count_crossings(eq, observed_domains) for eq in observed_eq_classes)
    print(f"\nObserved eQTL-gene pairs crossing domains: {observed_crossings}")

    random_crossings_counts = []
    for shuffled_domains in shuffled_domain_sequences:
        random_crossings_counts.append(sum(count_crossings(eq, shuffled_domains) for eq in observed_eq_classes))

    # Calculate P-value: proportion of random permutations with equal or more crossings
    p_value = sum(1 for c in random_crossings_counts if c >= observed_crossings) / len(random_crossings_counts)

    print(f"P-value for domain crossings (Observed >= Random): {p_value:.4e}")
    if p_value < 0.05:
        print("  Conclusion: Observed eQTL-gene pairs cross domains significantly more often than expected by chance (p < 0.05).")
    else:
        print("  Conclusion: No significant evidence that observed eQTL-gene pairs cross domains more often than expected by chance (p >= 0.05).")


def main(tissue):
    # 1. Load Data
    resol = 40000
    fragments, hic_interactions, genes, snps, eqtl_pairs, domains = load_data(tissue, resol)
    print(f"Loaded {len(fragments)} fragments, {len(hic_interactions)} Hi-C interactions, "
          f"{len(genes)} genes, {len(snps)} SNPs, {len(eqtl_pairs)} eQTL-gene pairs, "
          f"and {len(domains)} topological domains (mock data).")

    # 2. Create eQTL Equivalence Classes and Calculate Properties
    observed_eq_classes = create_eqtl_equivalence_classes(eqtl_pairs)
    print(f"Created {len(observed_eq_classes)} eQTL equivalence classes.")

    for eq_class in observed_eq_classes:
        eq_class.calculate_properties(hic_interactions)
        # print(f"EQ Class SNP: {eq_class.snp_fragment.fragment_id}, Gene: {eq_class.original_pairs[0].gene.gene_name}")
        # print(f"  d: {eq_class.snp_gene_distance}, l: {eq_class.gene_length}, ts: {eq_class.total_snp_frequency}, tg: {eq_class.maximum_total_gene_frequency}, p: {eq_class.spatial_proximity}")

    # 3. Generate Random eQTL Equivalence Classes
    random_eq_classes = generate_random_eqtl_equivalence_classes(
        observed_eq_classes,
        fragments, # Pass all fragments for random SNP/gene placement
        genes,     # Pass all genes for random gene selection
        hic_interactions,
        num_random_multiplier=10 # Reduced for quicker mock run
    )

    # 4. Perform Optimal Matching
    # Note: For k > 1, a more sophisticated matching algorithm or iterative approach is needed.
    # scipy.optimize.linear_sum_assignment is for k=1 (minimum total cost).
    matched_observed, matched_random = perform_optimal_matching(observed_eq_classes, random_eq_classes)
    print(f"\nPerformed optimal matching, found {len(matched_observed)} matched pairs.")

    # 5. Statistical Analysis of eQTLs (Spatial Proximity)
    statistical_analysis_eqtls(matched_observed, matched_random)

    # 6. Generate Shuffled Topological Domains
    shuffled_domain_sequences = generate_shuffled_domains(domains, chromosome_length=4000000, num_shuffles=100) # Reduced for mock run

    # 7. Analyze Domain Crossings
    analyze_domain_crossings(observed_eq_classes, domains, shuffled_domain_sequences)


    
tissues = {
"Adrenal Gland": "Adrenal Gland",
"Aorta": "Artery - Aorta",
"Dorsolateral Prefrontal Cortex": "Brain - Cortex",
"Hippocampus": "Brain - Hippocampus",
"Heart Left Ventricle": "Heart - Left Ventricle",
"Liver": "Liver",
"Lung": "Lung",
"Ovary": "Ovary",
"Pancreas": "Pancreas",
"Small Bowel (Intestine)": "Small Intestine - Terminal Ileum",
"Spleen": "Spleen"
}

for tissue in tissues.keys():
    main(tissue)

    
