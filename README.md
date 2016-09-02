# entrez_qiime

A utility to generate input files for taxonomy assignment in QIIME from the NCBI database

This repository contains Python code (entrez_qiime.py) required for a workflow that takes an input FASTA file generated from the NCBI database (e.g. through an Entrez search) and generates the files needed to BLAST your metabarcode data against those sequences using the QIIME script assign_taxonomy.py. These and other output files may also be useful for programs such as MEGAN. The repository includes a PDF (entrez_qiime.pdf) with guidelines for carrying out the full workflow.

[![DOI](https://zenodo.org/badge/24221/bakerccm/entrez_qiime.svg)](https://zenodo.org/badge/latestdoi/24221/bakerccm/entrez_qiime)
