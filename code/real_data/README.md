
This folder contains all the functions to re-produce the results in the real data analysis section.

The `generate_ref_p.R` script pre-generates necessary analysis parameters and save them in desired folder. Please change the directory accordingly.
The `read_data_func.R` script contains all the functions to perform the analysis. The function reads in the analysis parameters generated above.
To get the analysis results, run `real_data_run.R` script.

RNA-Sieve is implemented in Python, and the code for running RNA-Sieve is included in the Jupyter notebook, `rnasieve_real.ipynb`.
