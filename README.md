# Skipjack-bigeye lateral separation analysis

The code in this repository reproduces figures in the following paper:
<br>Leung, S., Mislan, K. A. S., & Thompson, L. ENSO drives lateral separation of FAD-associated skipjack and bigeye tuna in the Western Tropical Pacific. <i>Submitted to PLOS ONE.</i>

Please cite the above paper and the code itself (here: https://doi.org/10.5281/zenodo.3904134) if you use any of it.

This code was written using Python 3.7.3.

How to use this code:
1. Download this repository.
2. Download data from https://doi.org/10.5281/zenodo.3904157.
3. Build a docker container with all the required packages using Dockerfile inside the docker dir. (Or create a virtual environment with the packages listed in the Dockerfile.)
4. To recreate the figures from the above paper, simply run bet_skj_sep_main.ipynb in the python dir. Change dpath (path to the data dir you downloaded in step 2) and figpath (path to dir where you want to save figures to) as needed at the beginning of bet_skj_sep_main.ipynb.
