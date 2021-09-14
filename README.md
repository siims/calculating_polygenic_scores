# Calculating polygenic scores

[PGS Catalog](https://www.pgscatalog.org/) is `An open database of polygenic scores and the relevant metadata required for accurate application and evaluation.`

> A polygenic score (PGS) aggregates the effects of many genetic variants into a single number which predicts genetic predisposition for a phenotype. 

Repo has a set of notebooks and utility functions for calculating polygenic scores for a genome from raw reads.

It includes notebooks on:
1) how to run sequence alignment
2) how to run variant calling
3) how to pick alternate contigs
4) how to annotate variants
5) how to filter variants [TODO]
6) how to calculate individual's polygenic scores quickly

PGS Catalog doesn't have any special API keys, anyone can query it. So there isn't really anything except a raw human genome file (fastq) you need before starting (vcf file would allow to skip first 2-4 notebooks).

*Disclaimer: Genome analysis is computation heavy, some steps here might take a whole day to run depending on the hardware (ex: sequence alignment is especially heavy).*  

*TODO:*
* [ ] how to search the dna for a specific variance
* [ ] polygenic score interpretation. There are many different ways methods how to calculate PGS and score numbers results vary a lot. 

### Installation

```
# create a virtual environment with your favorite venv tool
sudo apt-get install libsnappy-dev # for pandas to_parquet and read_parquet
pip install -r requirements.txt
pre-commit install
```
