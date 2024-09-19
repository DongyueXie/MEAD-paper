
This folder contains all the functions to re-produce the results in the real data analysis section.


The generate_ref_p.R script pre-generates necessary analysis parameters and save them in desired folder.
The read_data_func.R script contains all the functions to perform the analysis. The function reads in the analysis parameters generated above.
To get the analysis results, run real_data_run.R script.

The dataset that used in the analysis are included in the data folder submission. Other dataset for MuSiC, CIBERSORT and RNA-Sieve can be downloaded here:  https://drive.google.com/drive/folders/1JznNX1Sh8Ty5hJlJfv6RFfX4KAYwc6KU?usp=sharing ; data for running RNA-Sieve is under  'data\neuron\rnasieve'.

RNA-Sieve is implemented in Python, and the code for running RNA-Sieve is included in the Jupyter notebook file, rnasieve_real.ipynb.
