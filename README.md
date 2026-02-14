# Nanocomporev2_manuscript
Analysis scripts and manuscript organization for the Nanocompore v2 manuscript


Content:
- `src/`: The directory contains the Julia scripts used for generating the plots in the manuscript.
    - `src/scripts`: Various scripts used in the study.
        - `src/scripts/test_single_molecule_comodifications_contingency.jl`: used for testing if pairs of predicted modified sites are associated an the molecule level.
        - `src/scripts/v1`: Scripts for running Nanocompore v1 on to prepare and compare the RNA002 native and m6A-depleted samples.
        - `src/scripts/RNA002`: Scripts used for preparing and comparing RNA002 samples with Nanocompore v2. Includes running aligning of the samples, running f5c eventalign + collapsing the output, and scripts to compare samples. Benchmarking scripts that compare v1 and v2 are also found here.
        - `src/scripts/RNA004`: Similar as above, but for the samples sequenced with the RNA004 chemistry. The scripts to run Nanocompore comparisons here also include different resquigglers (f5c eventalign, Uncalled4, and Remora) and evaluations at different coverage thresholds.
    - `src/notebooks`: This contains **very** messy Pluto notebooks used for exploration and experimentation.

