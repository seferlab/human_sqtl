     
def load_mock_data():
    """
    Loads mock data for Hi-C fragments, genes, SNPs, eQTLs, and domains.
    This is highly simplified and for conceptual illustration only.
    """
    # Mock Hi-C fragments (40kb resolution)
    mock_fragments = []
    for i in range(100): # 100 fragments on chromosome 1
        mock_fragments.append(HiCFragment('chr1', i * 40000, (i + 1) * 40000, f'frag_{i}'))

    # Mock Hi-C interaction data (fragment_id1, fragment_id2): frequency
    mock_hic_interactions = {}
    for i in range(99):
        # Local interactions
        mock_hic_interactions[tuple(sorted((f'frag_{i}', f'frag_{i+1}')))] = random.randint(10, 100)
    # Some 'long-range' interactions
    mock_hic_interactions[tuple(sorted(('frag_10', 'frag_50')))] = random.randint(5, 50)
    mock_hic_interactions[tuple(sorted(('frag_25', 'frag_75')))] = random.randint(5, 50)


    # Mock Genes (overlapping fragments)
    mock_genes = [
        Gene('geneA', 'chr1', 100000, 150000, [mock_fragments[2], mock_fragments[3]]),
        Gene('geneB', 'chr1', 500000, 580000, [mock_fragments[12], mock_fragments[13], mock_fragments[14]]),
        Gene('geneC', 'chr1', 1000000, 1050000, [mock_fragments[25], mock_fragments[26]]),
        Gene('geneD', 'chr1', 2000000, 2080000, [mock_fragments[50], mock_fragments[51]]),
    ]

    # Mock SNPs (mapping to fragments)
    mock_snps = [
        SNP('snp1', 'chr1', 110000, mock_fragments[2]), # Near geneA
        SNP('snp2', 'chr1', 520000, mock_fragments[13]), # Within geneB
        SNP('snp3', 'chr1', 900000, mock_fragments[22]), # Distal to geneC
        SNP('snp4', 'chr1', 1900000, mock_fragments[47]), # Distal to geneD
    ]

    # Mock eQTL-gene pairs
    mock_eqtl_gene_pairs = [
        EQTLGenePair(mock_snps[0], mock_genes[0]), # snp1 -> geneA (within-domain likely)
        EQTLGenePair(mock_snps[1], mock_genes[1]), # snp2 -> geneB (within-domain likely)
        EQTLGenePair(mock_snps[2], mock_genes[2]), # snp3 -> geneC (potentially across-domain)
        EQTLGenePair(mock_snps[3], mock_genes[3]), # snp4 -> geneD (potentially across-domain)
        # Add a duplicate for equivalence class testing
        EQTLGenePair(SNP('snp1_dup', 'chr1', 115000, mock_fragments[2]), mock_genes[0]),
    ]

    # Mock Topological Domains
    mock_domains = [
        TopologicalDomain('chr1', 0, 400000), # frag_0 to frag_9
        TopologicalDomain('chr1', 800000, 1200000), # frag_20 to frag_29
        TopologicalDomain('chr1', 1600000, 2000000), # frag_40 to frag_49
    ]

    return mock_fragments, mock_hic_interactions, mock_genes, mock_snps, mock_eqtl_gene_pairs, mock_domains
