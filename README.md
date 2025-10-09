# Overall-and-Serotype-Specific-Dengue-Transmission-Trends-in-Mexico

This repository includes .rds files containing the data used to run the catalytic models used in the paper 
"Overall and Serotype-Specific Dengue Transmission Trends in Mexico, 2016-2023", 
an R script ("Code to run models.r") showing how to run the models for given states in Mexico, 
and .STAN files of the three models which appear in theResults section of the paper. These models are:

*Model 1, the non-serotype-specific time-constant model from Imai et al. 2016
*Model 5, an adapted version of the non-serotype-specific time-varying model from O'Driscoll et al. 2019, modified to allow differential
reporting rates for children under 15 years old
*Model 6, a serotype-specific time-varying model developed for this paper

The models are run using the "rstan" package, and functions from the "tidyverse" package are needed to perform some of the data 
manipulation required to format the case and population data to then run the models.
