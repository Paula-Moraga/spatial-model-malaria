# spatial-modeling-malaria

The files in this repository correspond to the paper "Bayesian Spatial Modelling of Geostatistical Data using INLA and SPDE methods: A Case Study Predicting Malaria Risk in Mozambique" published in Spatial and Spatio-temporal Epidemiology.

Data `d.csv` contains prevalence survey data for Mozambique and selected covariates in surveyed locations (altitude `alt`, maximum temperature `temp`, precipitation `prec`, humidity `hum`, population density `pop` and distance to nearest inland water bodies `dist_aqua`).

Data `dp.csv` specifies the locations where we wish to predict the prevalence together with values of covariates in these locations.

Code `code.R` contains the R code to run the analysis and visualize the results.
