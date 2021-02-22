import json
import time
import traceback
from enum import Enum
from multiprocessing import Pool

import numpy as np
import pandas as pd
import shutil
from pathlib import Path
from typing import Tuple, Dict, List, Union

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


def get_all_pgs_api_data(api_endpoint: str, cache_dir: str):
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True, parents=True)
    cache_file = cache_dir / f"api_results_{api_endpoint.replace('/', '-')}.json"
    if Path(cache_file).exists():
        print(f"Found cache file {cache_file}. Loading data from cache.")
        with open(cache_file, "r") as f:
            return json.load(f)
    limit = 50
    offset = 0
    results = []
    while True:
        url = f"https://www.pgscatalog.org/rest/{api_endpoint}?limit={limit}&offset={offset}"
        print(f"Requesting pgs data from {url}")
        traits_response = requests.get(url=url)
        data = traits_response.json()
        results.extend(data["results"])
        if data["next"] == None:
            break
        offset += limit
    with open(cache_file, "w") as f:
        json.dump(results, f)
    return results


def download_file(url: str, local_filename: str) -> None:
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def read_or_download_pgs_scoring_file(pgs_id: str, cache_dir: Union[str, Path]) -> Tuple[str, Dict]:
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True, parents=True)
    cache_score_file = get_pgs_score_file_from_id(pgs_id, cache_dir)
    cache_file_response = Path(cache_dir) / f"{pgs_id}.json"

    if Path(cache_file_response).exists():
        with open(cache_file_response) as f:
            response_data = json.load(f)
    else:
        url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
        print(f"Requesting pgs {pgs_id} data from {url}. Cache file {cache_file_response} not found.")
        time.sleep(0.5)  # Not to overload api with requests - there is a 100 requests/minute limit
        data = requests.get(url)
        response_data = data.json()
        with open(cache_file_response, "w") as f:
            json.dump(response_data, f)

    if not Path(cache_score_file).exists():
        print(
            f"Downloading pgs {pgs_id} scoring file from {response_data['ftp_scoring_file']}. Cache file {cache_score_file} not found."
        )
        scoring_file_url = response_data["ftp_scoring_file"]
        download_file(scoring_file_url, cache_score_file)
    return cache_score_file, response_data


def to_gene_dosage_df(variance_str_list: List[str]) -> pd.DataFrame:
    if len(variance_str_list) == 0:
        return pd.DataFrame(columns=["rsid", "gene_dosage"])
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


def _calc_polygenic_score(vcf_file: str, pgs_file: str, hg19_rsid_chrom_pos_mapping_file: str, max_pgs_alleles=200):
    pgs_df = read_polygenic_score_file(pgs_file)
    assert pgs_df.shape[0] < max_pgs_alleles, f"Too many snps for {pgs_file}. Total {pgs_df.shape[0]}"

    is_pgs_file_with_rsids = "rsid" in pgs_df.columns
    if is_pgs_file_with_rsids:
        print(f"calc pgs based on rsid for {get_pgs_id_from_filename(pgs_file)}")
        pgs_locations = pd.DataFrame(pgs_df["rsid"], columns=["rsid"])
        pgs_rsids = clean_rsids(pgs_locations["rsid"], Path(pgs_file).stem)
        my_variance = search_for_rsids(pgs_rsids, my_vcf_file=vcf_file)
    else:
        print(f"calc pgs based on chr-pos {get_pgs_id_from_filename(pgs_file)}")
        assert pgs_df.attrs["metadata"]["hg_build"] == "hg19", (
            f"Can handle only pgs files with hg19 builds. Cannot handle {pgs_df.attrs['metadata']['hg_build']}"
        )

        pgs_locations, _ = fetch_hg19_rsids_based_on_chrom_pos(pgs_df, hg19_rsid_chrom_pos_mapping_file)
        my_variance = search_for_rsids(pgs_locations["rsid"], my_vcf_file=vcf_file)

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


def calc_all_polygenic_scores_parallel(
        pgs_ids: List[str],
        vcf_file: str,
        hg19_rsid_chrom_pos_mapping_file: str,
        max_pgs_alleles: int = 20_000,
        num_parallel_processes: int = 6,
        pgs_catalog_input_cache_dir: str = "data/pgs",
        pgs_result_cache_dir: str = "data/pgs_results",
):
    with Pool(num_parallel_processes) as p:
        pgs_series_list = p.map(do_calc_polygenic_score_single_input_arg,
                                [(pgs_id,
                                  vcf_file,
                                  hg19_rsid_chrom_pos_mapping_file,
                                  max_pgs_alleles,
                                  pgs_catalog_input_cache_dir,
                                  pgs_result_cache_dir)
                                 for pgs_id in pgs_ids])
        if len(pgs_series_list) == 0:
            return pd.DataFrame()
        else:
            return pd.concat(pgs_series_list, axis=1).T


def get_pgs_id_from_filename(filename: str) -> str:
    return Path(filename).stem[:-4]


def get_pgs_score_file_from_id(pgs_id: str, path_to_pgs_files: Union[str, Path]) -> str:
    return Path(path_to_pgs_files) / f"{pgs_id}.txt.gz"


def do_calc_polygenic_score(
        pgs_id: str,
        vcf_file: str,
        hg19_rsid_chrom_pos_mapping_file: str,
        max_pgs_alleles: int,
        pgs_catalog_input_cache_dir: str,
        pgs_result_cache_dir: str,
) -> pd.Series:
    cache_dir = Path(pgs_result_cache_dir)
    cache_dir.mkdir(exist_ok=True, parents=True)
    pgs_result_cache_file = cache_dir / f"{pgs_id}.tsv"
    if Path(pgs_result_cache_file).exists():
        results = pd.read_csv(pgs_result_cache_file, sep="\t", squeeze=True, header=None, index_col=0,
                              dtype={"error": str})
        return results

    downloaded_pgs_score_file, _ = read_or_download_pgs_scoring_file(pgs_id, pgs_catalog_input_cache_dir)
    result_df = pd.Series(
        [
            pgs_id,
            np.nan,
            np.nan
        ],
        index=["pgs_id", "score", "error"]
    )
    try:
        pgs, _ = _calc_polygenic_score(
            vcf_file=vcf_file,
            pgs_file=downloaded_pgs_score_file,
            hg19_rsid_chrom_pos_mapping_file=hg19_rsid_chrom_pos_mapping_file,
            max_pgs_alleles=max_pgs_alleles
        )
        result_df["score"] = pgs
        result_df.to_csv(pgs_result_cache_file, sep="\t", header=False)
        return result_df
    except Exception as e:
        # store exceptions in the cache as well
        errors = [str(e), ''.join(traceback.format_exception(None, e, e.__traceback__))]
        result_df["error"] = str(errors)
        result_df.to_csv(pgs_result_cache_file, sep="\t", header=False)
        return result_df


def do_calc_polygenic_score_single_input_arg(
        input_args: Tuple[str, str, str, int, str, str]
) -> pd.Series:
    return do_calc_polygenic_score(*input_args)


def get_pgs_metadata(pgs_id: str, pgs_catalog_input_cache_dir: str) -> pd.Series:
    downloaded_pgs_score_file, metadata_json = read_or_download_pgs_scoring_file(pgs_id, pgs_catalog_input_cache_dir)
    return pd.Series(
        [
            pgs_id,
            metadata_json["trait_reported"],
            PGS_METHOD_MAPPING_TO_METHOD_CATEGORIES.get(metadata_json["method_name"]),
            metadata_json["method_name"],
            [v["ancestry_broad"] for v in metadata_json["samples_variants"]]
        ],
        index=["pgs_id", "trait", "method_categorized", "method", "ancestry"]
    )
