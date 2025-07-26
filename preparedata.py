import numpy as np
import scipy as sp
import random
import ensembl_rest
from Ensembl_converter import EnsemblConverter
import pandas as pd
import pickle

def convert_sqtl_location(fname, outfname):
    # Create an instance of EnsemblConverter
    converter = EnsemblConverter()
    
    df = pd.read_csv(fname,sep="\t")
    values = set(df["geneId"])
    #values = list(values)[0:500]
    symbol2info = {}
    for value in values:
        symbol_raw = value.split(".")[0]
        symbol = converter.convert_ids([symbol_raw])["Symbol"][0]
        symbol = symbol.split("-")[0]
    
        try:
            result = ensembl_rest.symbol_lookup(
                species='homo sapiens',
                symbol=symbol
            )
            print(f"Symbol: {symbol}, Gene ID: {result['id']}, Chromosome: {result['seq_region_name']}, Start: {result['start']}, End: {result['end']}")
            symbol2info[symbol_raw] = result
        except ensembl_rest.HTTPError as err:
            error_code = err.response.status_code
            error_message = err.response.json()['error']
            if error_code == 400:
                if 'No valid lookup found for symbol' in error_message:
                    print(f"No valid lookup found for symbol: {symbol}")
                else:
                    print(f"HTTP Error {error_code}: {error_message}")
            elif error_code == 503:
                continue
            else:
                raise

    output = []
    for row in df.iterrows():
        gene = row[1]["geneId"]
        snp = row[1]["snpId"]
        chro, pos = snp.split("_")[0:2]
        pos = int(pos)

        symbol_raw = gene.split(".")[0]
        if symbol_raw in symbol2info:
            output.append({
                "Sqtl chromosome": chro,
                "Sqtl position": pos,
                "Gene chromosome": symbol2info[symbol_raw]["seq_region_name"],
                "Gene start": symbol2info[symbol_raw]["start"],
                "Gene end": symbol2info[symbol_raw]["end"],
                "Gene": symbol_raw,
            })

    output = pd.DataFrame(output)
    output.to_csv(outfname, index=False)

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

# Process SQTL Dataset
for tissue, mapped in tissues.items():
    fname = "sqtls-0.05fdr.permuted.tsv.gz"
    fpath = "dataset/sqtls/{0}/{1}".format(mapped,fname)
    outfname = "dataset/sqtls/{0}.csv".format(mapped)
    print("Processing {0}".format(fpath))
    convert_sqtl_location(fpath, outfname)
  
