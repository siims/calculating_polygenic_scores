import json
import shutil
from pathlib import Path

import requests

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
