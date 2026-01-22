# Social-Contact-Smoothing---LPS

This repository contains the code and data supporting the paper *“Estimating age patterns and grouped temporal trends in human contact patterns with Bayesian P-splines”* by Bryan Sumalinab, Oswaldo Gressani, Niel Hens, and Christel Faes.

It implements a Laplacian-P-splines (LPS) framework for estimating smooth age-specific contact rates and temporally grouped trends in social contact data. The core estimation procedure is provided in `contactLPS.R` which performs LPS estimation of age-specific contact rates and relies on `cubicbs.R` for the construction of the cubic B-spline basis.

The scripts `Holiday Groupings.R` and `Stringency Groupings.R` reproduce the analyses and figures for the holiday-based and stringency-based temporal groupings presented in the paper, respectively, using the datasets `data_holiday.RDS` and `data_stringency.RDS`.
