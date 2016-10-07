# entrez_qiime

A utility to generate input files for taxonomy assignment in QIIME from the NCBI database

This repository contains Python code (entrez_qiime.py) required for a workflow that takes an input FASTA file generated from the NCBI database (e.g. through an Entrez/gquery search) and generates the id-to-taxonomy file needed to BLAST your metabarcode data against those sequences using the QIIME script assign_taxonomy.py. The repository includes a PDF (entrez_qiime.pdf) with guidelines for carrying out the full workflow.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.159607.svg)](https://doi.org/10.5281/zenodo.159607)
