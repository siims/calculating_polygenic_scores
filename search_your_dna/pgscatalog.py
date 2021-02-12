import json
import time
import traceback
from enum import Enum

import numpy as np
import pandas as pd
import shutil
from pathlib import Path
from typing import Tuple, Dict, Optional, List

import pysam
import requests
from tqdm import tqdm

from search_your_dna.util import read_raw_zipped_polygenic_score_file, search_for_rsids, \
    read_raw_zipped_polygenic_score_file_with_chrom_pos


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
        print(
            f"PGS {pgs_name} has {np.count_nonzero(values_with_commas_or_underscores)} rsids containing multiple values")
        rsids = rsids[~values_with_commas_or_underscores]
    return rsids.to_list()


def fetch_hg19_rsids_based_on_chrom_pos(
        df: pd.DataFrame,
        hg19_rsid_chrom_pos_mapping_file: str
) -> Tuple[pd.DataFrame, List[str]]:
    """
    :param df: data frame with expected columns: "chrom", "pos", "effect_allele", "reference_allele"
    :return: data frame with columns "rsid", "chrom", "pos" that were present in the db
    """
    assert {"chrom", "pos", "effect_allele", "reference_allele"} - set(df.columns) != {}, "Missing columns"
    hg19_vcf = pysam.Tabixfile(hg19_rsid_chrom_pos_mapping_file)

    pgs_locations = pd.DataFrame(columns=["rsid", "chrom", "pos"])
    warnings = []
    for _, pgs_row in df.iterrows():
        counter = 0
        for rsid_mapping_file_row in hg19_vcf.fetch(pgs_row["chrom"], pgs_row["pos"] - 1, pgs_row["pos"]):
            values = rsid_mapping_file_row.split("\t")
            chrom, pos, _, reference_allele, effect_allele, rsid = tuple(values)
            pos = int(pos)
            if pgs_row["pos"] == pos and pgs_row["effect_allele"] == effect_allele and pgs_row[
                "reference_allele"] == reference_allele:
                pgs_locations = pgs_locations.append({"rsid": rsid, "chrom": chrom, "pos": pos}, ignore_index=True)
                if counter != 0:
                    warnings.append(f"multipe rsids for {pgs_row['chrom']}, {pgs_row['pos']}. {rsid_mapping_file_row}")
                counter += 1
    return pgs_locations, warnings


def do_polygenic_score_calculation(gene_dosage: pd.Series, effect_weight: pd.Series) -> Tuple[float, pd.Series]:
    pgs_effect = gene_dosage * effect_weight
    return pgs_effect.sum(), pgs_effect


def read_polygenic_score_file(pgs_file):
    try:
        return read_raw_zipped_polygenic_score_file(pgs_file)
    except Exception as e:
        first_error = [str(e), ''.join(traceback.format_exception(None, e, e.__traceback__))]
    try:
        return read_raw_zipped_polygenic_score_file_with_chrom_pos(pgs_file)
    except Exception as e:
        second_error = [str(e), ''.join(traceback.format_exception(None, e, e.__traceback__))]

    if first_error is not None and second_error is not None:
        raise Exception(f"ERROR: read_polygenic_score_file failed. Exceptions: {second_error + first_error}")


def calc_polygenic_score(my_vcf_file: str, pgs_file: str, hg19_rsid_chrom_pos_mapping_file: str, max_pgs_alleles=200):
    pgs_df = read_polygenic_score_file(pgs_file)
    assert pgs_df.shape[0] < max_pgs_alleles, f"Too many snps for {pgs_file}. Total {pgs_df.shape[0]}"

    is_pgs_file_with_rsids = "rsid" in pgs_df.columns
    if is_pgs_file_with_rsids:
        print("calc pgs based on rsid")
        pgs_locations = pd.DataFrame(pgs_df["rsid"], columns=["rsid"])
        pgs_rsids = clean_rsids(pgs_locations["rsid"], Path(pgs_file).stem)
        my_variance = search_for_rsids(pgs_rsids, file_my_vcf=my_vcf_file)
    else:
        print("calc pgs based on chr-pos")
        assert pgs_df.attrs["metadata"]["hg_build"] == "hg19", (
            f"Can handle only pgs files with hg19 builds. Cannot handle {pgs_df.attrs['metadata']['hg_build']}"
        )

        pgs_locations, _ = fetch_hg19_rsids_based_on_chrom_pos(pgs_df, hg19_rsid_chrom_pos_mapping_file)
        my_variance = search_for_rsids(pgs_locations["rsid"], file_my_vcf=my_vcf_file)

    gene_dosage_df = to_gene_dosage_df(my_variance)

    if is_pgs_file_with_rsids:
        merged_df = pgs_df.merge(gene_dosage_df, on="rsid", how="inner")
    else:
        gene_dosage_df = gene_dosage_df.merge(pgs_locations, how="outer", on="rsid")
        merged_df = pgs_df.merge(gene_dosage_df, on=["chrom", "pos"], how="inner")
        merged_df["pos"] = merged_df["pos"].astype("Int64")

    # assuming that if gene_dosage is nan then it wasn't found in my genome
    merged_df["gene_dosage"] = merged_df["gene_dosage"].fillna(0)
    pgs_score, merged_df["effect"] = do_polygenic_score_calculation(merged_df["gene_dosage"],
                                                                    merged_df["effect_weight"])
    return pgs_score, merged_df


def calc_all_polygenic_scores(
        files: List[str],
        file_my_vcf: str,
        hg19_rsid_chrom_pos_mapping_file: str = "data/humandb/hg19_avsnp150.txt.gz"
) -> Tuple[pd.DataFrame, Dict]:
    def get_pgs_id_from_filename(filename: str) -> str:
        return Path(filename).stem[:-4]

    errors = {}
    all_pgs_scores = pd.DataFrame(columns=["pgs_id", "trait", "score", "method_categorized", "method", "file"])
    for pgs_file in tqdm(sorted(files)):
        try:
            pgs_id = get_pgs_id_from_filename(pgs_file)
            # TODO: add caching
            pgs, _ = calc_polygenic_score(my_vcf_file=file_my_vcf,
                                          pgs_file=pgs_file,
                                          hg19_rsid_chrom_pos_mapping_file=hg19_rsid_chrom_pos_mapping_file,
                                          max_pgs_alleles=200)

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
