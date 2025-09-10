# protein-foundations-ukb-olink-50k

# protein‑foundations‑ukb‑olink‑50k

This repository contains the analysis code accompanying the paper **“Machine learning‑guided deconvolution of plasma protein levels”** (medRxiv preprint) ([medrxiv.org](https://www.medrxiv.org/content/10.1101/2025.01.09.25320257v2)). The study applies machine learning methods to disentangle the biological determinants of plasma protein variation, leveraging Olink proteomic data from \~43,240 UK Biobank participants.

## Important Notice

The code provided here is **not intended to work out of the box**. Adaptation of file paths, auxiliary files, and other settings is required to successfully run the scripts. Please refer to the paper for details on the necessary input files, data structure, and preprocessing steps to replicate the analysis.

## Overview

Population-scale plasma proteomics offers a powerful window into the biological pathways that underlie variation in protein levels, aging, and disease risk. This work focuses on identifying the contributions of genetics, demographic factors, and other biological determinants through ML‑based deconvolution of ∼3,000 plasma protein measurements in a cohort of over 43,000 participants ([medrxiv.org](https://www.medrxiv.org/content/10.1101/2025.01.09.25320257v2)).

---

## Repository Structure

The repository is organized into five numbered folders, each corresponding to key steps of the analysis pipeline:

* `01_data_preparation/`
  Scripts to preprocess raw Olink proteomic data and covariates, including normalization, filtering, and integration.

* `02_feature_selection/`
  Code to identify relevant features (e.g., genetic variants, demographic variables) for downstream modeling.

* `03_regression_analysis/`
  Implements the machine learning and regression frameworks to quantify the variance in protein levels explained by selected features.

* `04_drug_disease_network/`
  Constructs networks linking proteins with drug targets and disease associations based on analysis results.

* `05_prot_dis_assoc/`
  Investigates protein–disease associations, such as predictive protein biomarkers.

Additional files:

* `README.md` — (this document)
* `LICENSE` — distributed under the **GPL‑3.0 License** ([github.com](https://github.com/comp-med/protein-foundations-ukb-olink-50k))

---

## Getting Started

### Prerequisites

* **R** (version ≥ 4.0 recommended)
* Essential R packages (e.g., `data.table`, `tidyverse`, `glmnet`, `randomForest`, plus domain-specific libs)
* Access to processed Olink data and covariates via UK Biobank / UKB‑PPP
* Optional: high-performance compute infrastructure for large-scale modeling

---

## How to Cite

Please cite the following preprint when using this repository:

> *Machine learning‑guided deconvolution of plasma protein levels*, medRxiv (preprint), DOI: 10.1101/2025.01.09.25320257v1 ([medrxiv.org](https://www.medrxiv.org/content/10.1101/2025.01.09.25320257v2))

---

## Licensing

This project is licensed under the **GPL‑3.0 License**, which permits reuse and modification provided derivative work maintains the same license ([github.com](https://github.com/comp-med/protein-foundations-ukb-olink-50k)).

---

## Acknowledgements & Related Tools

* The work is part of the **UK Biobank Pharma Proteomics Project (UKB‑PPP)**, involving high-throughput Olink assays of \~3,000 plasma proteins in UK Biobank participants ([pmc.ncbi.nlm.nih.gov](https://www.nature.com/articles/s41586-023-06592-6), [omicscience.org](https://omicscience.org/apps/prot_foundation/)).
* A related **web-based visualization tool** enables exploration of protein variance and explained variance from the study ([omicscience.org](https://omicscience.org/apps/prot_foundation/)).

---

## Contact & Support

For questions regarding the analysis pipeline, data access, or usage:

* **GitHub issues** — report bugs or raise questions
* **Email** — reach the corresponding authors as listed in the preprint

---

## Summary Table

| Feature                 | Description                                                                                  |
| ----------------------- | -------------------------------------------------------------------------------------------- |
| **Purpose**             | ML-driven deconvolution of protein-level variation                                           |
| **Data Source**         | \~3,000 proteins measured via Olink in \~43,240 UK Biobank participants                      |
| **Repository Language** | R (primary)                                                                                  |
| **Steps Covered**       | Data prep → Feature selection → Regression → Drug/disease network → Protein–disease analysis |
| **License**             | GPL‑3.0                                                                                      |
| **Out-of-the-box**      | Code requires adaptation of paths and auxiliary files as described in the paper              |

