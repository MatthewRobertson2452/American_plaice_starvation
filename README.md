American plaice starvation
================
Matthew Robertson
2025-02-14

# Citation

Robertson, M.D. 2025. Impact of starvation-induced mortality on the
collapse and lack of recovery of American plaice on the Newfoundland
Grand Banks: Zenodo code release.
<https://doi.org/10.5281/zenodo.14872897>.

# Contact

Author - [Matthew D. Robertson](matthew.robertson@mi.mun.ca), Marine
Institute of Memorial University

# Summary

This is the code accompanying Robertson et al. (Accepted) where we
investigate if starvation mortality is an important part of the recent,
elevated natural mortality rates for the American plaice population on
the Newfoundland Grand Banks.

# Code purpose and order

This code is intended to allow readers/reviewers the ability to repeat
analyses in the paper.

1.  The first script to run is `run_SpTem_WL.R`. This script will run
    the spatiotemporal weight-length model (`SpTemp_WL.cpp`) for
    American plaice and calculate the starvation-induced mortality
    index.

2.  The second script to run is `run_state_space_popdy.R` then uses the
    starvation-induced mortality index from the previously run model as
    an input in a state-space population dynamics model
    (`state_space_pop_dy.cpp`) to assess the contribution of starvation
    on total American plaice natural mortality.

# File descriptions

`R\run_SpTem_WL.R` - R script for loading data and running the
spatiotemporal weight-length model.

`R\run_state_space_popdy.R` - R script for loading data and running the
state-space population dynamics model.

`src\SpTemp_WL.cpp` - c++ code for the spatiotemporal weight-length
model.

`src\state_space_pop_dy.cpp` - c++ code for the state-space population
dynamics model.

`Data\distance_mat_strat` - data file containing distance matrix between
survey strata, used in the spatiotemporal weight-length model.

`Data\strat_area` - data file containing survey strata spatial areas,
used in the spatiotemporal weight-length model.

`Data\lw_mod_inputs.RData` - Processed data inputs for the
spatiotemporal weight-length model.

`Data\pop_dy_inputs.RData` - Processed data inputs for the state-space
population dynamics model.

# Disclaimer

Data were provided by Fisheries and Oceans Canada. Processed data files
are available in this repository. Raw data can be accessed through a
data request by contacting the co-author ([Laura
Wheeland](laura.wheeland@dfo-mpo.gc.ca)) from Fisheries and Oceans
Canada.

# Acknowledgements

This research associated with this data release was supported by an
NSERC Discovery Grant and Marine Institute of Memorial University
start-up funding.
