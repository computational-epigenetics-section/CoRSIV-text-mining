# CoRSIV-text-mining

Repository containing code and analyses for the paper "Human systemic epigenetic variants are implicated in neurodevelopmental and metabolic disorders". 


## Repository Structure

- `probe_collection/` - Scripts to search PubMed articles and extract probes for disease categories
  - `pubmed_search_pipeline.ipynb` - Main script to collect, scrape, and extract probes from PubMed articles
  - `handle_zip.py` - Helper script to handle zip files
  - `get_probe_supplementary.py` - Helper script to extract probes from articles' supplementary files

- `permutation_testing/` - Suite for permutation testing to evaluate the statistical significance of CoRSIV enrichment, refer to `permutation_testing/README.md` for more details.

- `corsiv_regions/` - Scripts and data to merge previously identified SIV regions into a unified list of CoRSIVs.
  - `CoRSIV_annotation.ipynb` - Main script to annotate and merge SIV regions.
  - `SIV.hg38.bed` / `ME.hg38.bed` / `ESS.hg38.bed` / `corsiv2019.txt` - input SIV, ME, ESS, and CoRSIV regions.

- `controls/` - Scripts to generate control regions / probes.
  - `generate_lookup.ipynb` - Generate lookup table on a chromosome-by-chromosome basis that are later used to sample control regions from.
  - `sample_from_lookup.ipynb` - Sample control regions from lookup table based on CoRSIV metrics.
  - `process_control.ipynb` - Post-processing of control regions.

- `GSEA.ipynb` - Scripts to perform gene-set enrichment analysis on CoRSIVs.
- `median_iir_icc.ipynb` - Scripts to calculate median IRR and ICC values for different region types and categories.
- `main_figure_code.ipynb` - Code to generate main figures.
- `supplementary_figure_code.ipynb` - Code to generate supplementary figures.
- `util.py` - Utility functions used in the above scripts.

## Requirements

- Python 3.7+
- High-performance computing environment recommended for probe extraction and permutation testing.

