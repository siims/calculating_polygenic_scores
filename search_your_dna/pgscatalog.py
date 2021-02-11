import json
import time
import traceback
from enum import Enum

import numpy as np
import pandas as pd
import shutil
from pathlib import Path
from typing import Tuple, Dict, Optional, List

import requests
from tqdm import tqdm

from search_your_dna.snp_store import search_for_rsids
from search_your_dna.util import read_raw_zipped_polygenic_score_file


class MethodCategories(Enum):
    UNKNOWN = 0,
    GWAS_SNPS = 1,
    CLUMPING_THRESHOLDING = 2,
    PRUNING_THRESHOLDING = 3,


PGS_METHOD_MAPPING_TO_METHOD_CATEGORIES = {
    'GWAS-significant SNPs': MethodCategories.GWAS_SNPS,
    'GWAS significant SNPs': MethodCategories.GWAS_SNPS,
    'GWAS-sgnificant SNPs': MethodCategories.GWAS_SNPS,
    'GWAS significant variants': MethodCategories.GWAS_SNPS,
    'Genome-wide significant SNPs': MethodCategories.GWAS_SNPS,
    'Genome-wide significant variants.': MethodCategories.GWAS_SNPS,
    'Genome-wide significant variants': MethodCategories.GWAS_SNPS,
    'GWAS': MethodCategories.GWAS_SNPS,
    'GWAS Catalog SNPs (unweighted allele count)': "GWAS_SNPS_UNWEIGHTED",
    'GWAS Catalog SNPs': MethodCategories.GWAS_SNPS,
    'GWAS-significant (lead and secondary) variants': MethodCategories.GWAS_SNPS,
    'GWAS-significant variants': MethodCategories.GWAS_SNPS,
    'Genome-wide significant associations': MethodCategories.GWAS_SNPS,
    'Likelihood ratios for genomewide-significant SNPs': MethodCategories.GWAS_SNPS,
    'Genomewide-significant variants (sourced from PMID:28509669)': MethodCategories.GWAS_SNPS,
    'SNPs passing genome-wide significance': MethodCategories.GWAS_SNPS,
    'Genomewide-significant SNPs': MethodCategories.GWAS_SNPS,
    'Genomewide-significant SNPs, filtered to be specific to Chinese ancestry individuals': MethodCategories.GWAS_SNPS,
    'Genome-wide significant associations and interaction modelling': MethodCategories.GWAS_SNPS,
    'Known susceptibility variants (genome-wide significant SNPs)': MethodCategories.GWAS_SNPS,
    'GWAS-significant variants, HLA-specific significant variants.': MethodCategories.GWAS_SNPS,
    'Known susceptibility loci (genome-wide significant SNPs)': MethodCategories.GWAS_SNPS,
    'Clumping + Thresholding': MethodCategories.CLUMPING_THRESHOLDING,
    'LD-clumping and  p-value thresholding (Ricopilli)': MethodCategories.CLUMPING_THRESHOLDING,
    'Pruning and Thresholding (P+T)': MethodCategories.PRUNING_THRESHOLDING,
    'P+T': MethodCategories.PRUNING_THRESHOLDING,
    'Pruning + Thresholding': MethodCategories.PRUNING_THRESHOLDING,
    'Clumping and Thresholding (C+T)': MethodCategories.CLUMPING_THRESHOLDING,
    'Pruning + Clumping ': MethodCategories.PRUNING_THRESHOLDING,
    'LD thinning': MethodCategories.UNKNOWN,
    'PRS-CS': MethodCategories.UNKNOWN,
    'Log‐additive GRS': MethodCategories.UNKNOWN,
    'Polygenic Hazard Score': MethodCategories.UNKNOWN,
    'Weighted sum of risk alleles (with random forest selection)': MethodCategories.UNKNOWN,
    'Hazard model with stepwise selection of SNP inclusion': MethodCategories.UNKNOWN,
    'Lassosum': MethodCategories.UNKNOWN,
    '313 variants from Mavaddatt et al (PGS000005)': MethodCategories.UNKNOWN,
    'Curated variant associations': MethodCategories.UNKNOWN,
    'LASSO (biglasso)': MethodCategories.UNKNOWN,
    'Hard thresholding': MethodCategories.UNKNOWN,
    'Composite PRS (component scores combined by log(HR))': MethodCategories.UNKNOWN,
    'log-OR weighted sum of risk allele dosages': MethodCategories.UNKNOWN,
    'metaGRS': MethodCategories.UNKNOWN,
    '313 variants from Mavaddatt et al (PGS000006)': MethodCategories.UNKNOWN,
    'Weighted sum of risk alleles': MethodCategories.UNKNOWN,
    'LD Clumping, PRSice-2': MethodCategories.UNKNOWN,
    'SNP associations curated from the literature': MethodCategories.UNKNOWN,
    '313 variants from Mavaddatt et al (PGS000004)': MethodCategories.UNKNOWN,
    'LDpred': MethodCategories.UNKNOWN,
    'Fixed-effects two-stage polytomous model': MethodCategories.UNKNOWN,
    'snpnet': MethodCategories.UNKNOWN,
    'Established Lipid Loci; independent genome-wide significant variants': MethodCategories.UNKNOWN,
    'PRSice': MethodCategories.UNKNOWN,
    'SparSNP': MethodCategories.UNKNOWN,
    'composite likelihood ratio': MethodCategories.UNKNOWN,
    'Established lipid loci': MethodCategories.UNKNOWN,
    'Literature-derived SNP selection': MethodCategories.UNKNOWN,
    'metaGRS of 19 component PGSs': MethodCategories.UNKNOWN
}


def get_all_pgs_api_data(api_endpoint: str):
    cache_file = f"data/pgs/pgs_catalog_{api_endpoint.replace('/', '-')}.json"
    if Path(cache_file).exists():
        print(f"Found cache file {cache_file}. Loading data from cache.")
        with open(cache_file, "r") as f:
            return json.load(f)
    limit = 50
    offset = 0
    traits = []
    while True:
        url = f"https://www.pgscatalog.org/rest/{api_endpoint}?limit={limit}&offset={offset}"
        print(f"Requesting pgs data from {url}")
        traits_response = requests.get(url=url)
        data = traits_response.json()
        traits.extend(data["results"])
        if data["next"] == None:
            break
        offset += limit
    with open(cache_file, "w") as f:
        json.dump(traits, f)
    return traits


def download_file(url: str, local_filename: str) -> None:
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def read_or_download_pgs_scoring_file(pgs_id: str) -> Tuple[Dict, pd.DataFrame]:
    cache_file_score_file = f"data/pgs/{pgs_id}.txt.gz"
    cache_file_response = f"data/pgs/{pgs_id}.json"

    if Path(cache_file_response).exists():
        with open(cache_file_response) as f:
            response_data = json.load(f)
    else:
        url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
        print(f"Requesting pgs {pgs_id} data from {url}. Cache file {cache_file_response} not found.")
        time.sleep(0.5)  # Not to overload api with requests
        data = requests.get(url)
        response_data = data.json()
        with open(cache_file_response, "w") as f:
            json.dump(response_data, f)

    if not Path(cache_file_score_file).exists():
        print(
            f"Downloading pgs {pgs_id} scoring file from {response_data['ftp_scoring_file']}. Cache file {cache_file_score_file} not found."
        )
        scoring_file_url = response_data["ftp_scoring_file"]
        download_file(scoring_file_url, cache_file_score_file)
    return response_data, read_raw_zipped_polygenic_score_file(cache_file_score_file)


def to_gene_dosage_df(variance_str_list: List[str]) -> pd.DataFrame:
    gene_dosages = []
    for v in variance_str_list:
        chrom, pos, rsid, ref, alt, qual, filter, info, format, _ = tuple(v.split())
        gene_dosage = int(info[info.find("AC") + 3:info.find("AC") + 4])
        gene_dosages.append({"rsid": rsid, "gene_dosage": gene_dosage})
    res = pd.DataFrame(gene_dosages)
    res["gene_dosage"] = pd.to_numeric(res["gene_dosage"])
    return res


def clean_rsids(rsids: pd.Series, pgs_name: str) -> List[str]:
    if np.any(rsids.isna()):
        print(f"PGS {pgs_name} has {np.count_nonzero(rsids.isna())} missing rsids")
        rsids = rsids.dropna()
    start_with_rs = rsids.str.startswith("rs")
    if np.any(~start_with_rs):
        print(f"PGS {pgs_name} has {np.count_nonzero(~start_with_rs)} non rsid values")
        rsids = rsids[start_with_rs]
    values_with_commas_or_underscores = rsids.str.contains(",") | rsids.str.contains("_")
    if np.any(values_with_commas_or_underscores):
        print(f"PGS {pgs_name} has {np.count_nonzero(values_with_commas_or_underscores)} rsids containing multiple values")
        rsids = rsids[~values_with_commas_or_underscores]
    return rsids.to_list()


def calc_polygenic_score(
        max_pgs_alleles: Optional[int] = None,
        pgs_df: Optional[pd.DataFrame] = None,
        pgs_file: Optional[str] = None,
        file_my_vcf: str = "data/GFX0237425.GRCh38.p7.annotated.hg38_multianno.updated.vcf.gz",
) -> Tuple[float, pd.DataFrame]:
    if pgs_df is None:
        pgs_df = read_raw_zipped_polygenic_score_file(pgs_file)
    if max_pgs_alleles is not None and len(pgs_df.index) > max_pgs_alleles:
        raise Exception(f"Too many snps for {pgs_file}. Total {len(pgs_df.index)}")
    pgs_rsids = clean_rsids(pgs_df['rsid'], Path(pgs_file).stem)
    my_variance = search_for_rsids(pgs_rsids, file_my_vcf=file_my_vcf)
    gene_dosage_df = to_gene_dosage_df(my_variance)
    merged_df = pgs_df.merge(gene_dosage_df, on="rsid", how="outer")
    merged_df["gene_dosage"] = merged_df["gene_dosage"].fillna(0) # assuming that if rsid wasn't present in the vcf I don't have it
    merged_df["effect"] = merged_df["gene_dosage"] * merged_df["effect_weight"]
    return merged_df["effect"].sum(), merged_df


def calc_all_polygenic_scores(files: List[str], file_my_vcf: str) -> Tuple[pd.DataFrame, Dict]:
    def get_pgs_id_from_filename(filename: str) -> str:
        return Path(filename).stem[:-4]

    errors = {}
    all_pgs_scores = pd.DataFrame(columns=["pgs_id", "trait", "score", "method_categorized", "method", "file"])
    for pgs_file in tqdm(sorted(files)):
        try:
            pgs_id = get_pgs_id_from_filename(pgs_file)
            pgs, _ = calc_polygenic_score(max_pgs_alleles=200, pgs_file=pgs_file, file_my_vcf=file_my_vcf)

            metadata_json, _ = read_or_download_pgs_scoring_file(pgs_id)

            all_pgs_scores = all_pgs_scores.append(
                {
                    "pgs_id": pgs_id,
                    "trait": metadata_json["trait_reported"],
                    "score": pgs,
                    "method_categorized": PGS_METHOD_MAPPING_TO_METHOD_CATEGORIES.get(metadata_json["method_name"]),
                    "method": metadata_json["method_name"],
                    "file": pgs_file
                }, ignore_index=True)
        except Exception as e:
            errors[pgs_file] = [str(e), ''.join(traceback.format_exception(None, e, e.__traceback__))]
    return all_pgs_scores, errors
