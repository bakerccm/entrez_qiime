This repository contains Python code (entrez_qiime.py) required for a workflow that takes an input FASTA file generated from the NCBI database (e.g. through an Entrez/gquery search) and generates the id-to-taxonomy file needed to BLAST your metabarcode data against those sequences using the QIIME script assign_taxonomy.py. The repository includes a PDF (entrez_qiime.pdf) with guidelines for carrying out the full workflow.

- GitHub repository: [https://github.com/bakerccm/entrez_qiime](https://github.com/bakerccm/entrez_qiime)

- Latest release (v2.0): [https://github.com/bakerccm/entrez_qiime/releases/tag/v2.0](https://github.com/bakerccm/entrez_qiime/releases/tag/v2.0)

    Version v2.0 involves a major update to entrez_qiime.py and its documentation. This script now takes an input FASTA file with NCBI accession numbers as sequence identifiers, in preparation for the NCBI phasing out GI numbers in sequence downloads. See the release documentation for additional changes.

- Initial release (v1.0): [https://github.com/bakerccm/entrez_qiime/releases/tag/v1.0](https://github.com/bakerccm/entrez_qiime/releases/tag/v1.0)

    Version v1.0 is the formal GitHub release of the original utility, which has been circulating outside GitHub since 2012. In v1.0, the Python script entrez_qiime.py takes an input FASTA file with GI numbers as sequence identifiers.
