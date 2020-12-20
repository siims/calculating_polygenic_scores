import json
import pandas as pd
import shutil
from pathlib import Path
from typing import Tuple

import requests

from search_your_dna.snp_store import query_my_genotypes_for_rsids
from search_your_dna.util import read_raw_zipped_polygenic_score_file


def get_all_pgs_api_data(api_endpoint: str):
    cache_file = f"data/pgs_catalog_{api_endpoint.replace('/', '-')}.json"
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


def read_or_download_pgs_scoring_file(pgs_id: str):
    cache_file = f"data/{pgs_id}.txt.gz"
    if Path(cache_file).exists():
        print(f"Found cache file {cache_file}. Loading data from cache.")
        return read_raw_zipped_polygenic_score_file(cache_file)
    url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
    print(f"Requesting pgs data from {url}")
    data = requests.get(url)
    response_data = data.json()
    scoring_file_url = response_data["ftp_scoring_file"]
    download_file(scoring_file_url, cache_file)
    return read_raw_zipped_polygenic_score_file(cache_file)


def calc_polygenic_score(snp_db_file: str, pgs_file: str, max_pgs_alleles: int) -> Tuple[float, pd.DataFrame]:
    pgs_df = read_raw_zipped_polygenic_score_file(pgs_file)
    if len(pgs_df.index) > max_pgs_alleles:
        raise Exception(f"Too many snps for {pgs_file}. Total {len(pgs_df.index)}")
    genotype = query_my_genotypes_for_rsids(snp_db_file, pgs_df["rsid"].to_list())
    merged_df = pgs_df.merge(genotype, on="rsid")
    merged_df["effect_allele_1"] = merged_df["genotype"].map(lambda x: x[0]) == merged_df["effect_allele"]
    merged_df["effect_allele_2"] = merged_df["genotype"].map(lambda x: x[1]) == merged_df["effect_allele"]
    merged_df["effect_allele_1"] = merged_df["effect_allele_1"].astype(int)
    merged_df["effect_allele_2"] = merged_df["effect_allele_2"].astype(int)
    merged_df["gene_dosage"] = merged_df["effect_allele_1"] + merged_df["effect_allele_2"]
    merged_df["effect"] = merged_df["gene_dosage"] * merged_df["effect_weight"]
    return merged_df["effect"].sum(), merged_df
