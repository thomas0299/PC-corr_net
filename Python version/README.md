# PC-corr algorithm associated to PCA analysis

## Released under MIT License
Copyright (c) 16 Dec 2017 Sara Ciucci, Yan Ge, Claudio Durán and Carlo Vittorio Cannistraci

## Please cite:
Enlightening discriminative network functional modules behind Principal Component Analysis separation in differential-omic science studies.
Sara Ciucci, Yan Ge, Claudio Durán?n, Alessandra Palladini, Víctor Jiménez Jiménez, Luisa María Martínez Sánchez, 
Yuting Wang, Susanne Sales, Andrej Shevchenko, Steven W. Poser, Maik Herbig, Oliver Otto, Andreas Androutsellis-Theotokis, 
Jochen Guck, Mathias J. Gerl and Carlo Vittorio Cannistraci 
Scientific Reports, 2017

## Setup

1. Create virtual environment

`python -m venv pc_corr_env`

2. Download dependencies

`pip install requirements.txt`

3. Activate environment

`source pc_corr_env/bin/activate`

## Run

`python PC_corr.py path_to_csv.csv --sample_labels="after,before,after..." --default_colors="green,red"`


## Note

Not all features from the R/MATLAB code have been translated to Python. Following are missing features which will be added in the future:
- More than 2 groups in sample labels.
- Non-centered PCA (Currently non-centered PCA results are the same as centered PCA).
- Trustworthiness measure of segregation.
        