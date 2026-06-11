# StratRecRes_supp

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20644847.svg)](https://doi.org/10.5281/zenodo.20644847)

Code and supplementary data for "Sedimentation rates determine temporal resolution and disorder of the stratigraphic record".

## Authors

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

## Requirements

R version >= 4.2, optionally Matlab 2006a or more recent.

## Reproduction

For details on reproduction of the study see REPRODUCEME.md.

## Repository structure

* code : directory with matlab and r code
* data : directory with data
  * matlab_outputs
    * ADDs_from_matlab.mat : Matlab simulation results
    * codebook.md : description of contents of ADDs_from_matlab.mat
  * r_outputs : initially empty, filled with intermediate data after running of `analysis.R`
  * res : initially empty, filled with results after running `analysis.R`
  * smlbase : SMLBase v1.04, Hohmann (2022)
  * tomasovych_et_al_2018: Raw data from Tomasovych et al. (2018)
* figs : initially empty, filled with figures after running `analysis.R`
* renv : directory for the `renv` package
* .gitignore : untracked files
* .Rprofile : R session setup
* LICENSE : Apache 2.0 license text
* README : Readme file
* renv.lock : lock file for `renv` package
* REPRODUCEME.md : Instructions for reproduction of the results
* StratRecRes_supp.Rproj : Rproject file

## References

This repository contains data from

* Hohmann, Niklas. 2022. “Global Compilation of Surface Mixed Layer Parameters (Sedimentation Rate, Bioturbation Depth, Mixing Intensity) from Marine Environments: The SMLBase v1.0.” Frontiers in Earth Science 10 (December):1013174. [DOI: 10.3389/feart.2022.1013174](https://doi.org/10.3389/feart.2022.1013174).
* Tomašových, Adam, Ivo Gallmetzer, Alexandra Haselmair, Darrell S. Kaufman, Martina Kralj, Daniele Cassin, Roberto Zonta, and Martin Zuschin. 2018. “Tracing the Effects of Eutrophication on Molluscan Communities in Sediment Cores: Outbreaks of an Opportunistic Species Coincide with Reduced Bioturbation and High Frequency of Hypoxia in the Adriatic Sea.” Paleobiology 44 (4): 575–602. [DOI: 10.1017/pab.2018.22](https://doi.org/10.1017/pab.2018.22).

## Citation

To cite this repository, please use

* Hohmann, N. (2026). Supplementary Code for "Sediment accumulation controls the depositional resolution of the stratigraphic record" (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.20644847

## Copyright

Copyright 2023-2026 Netherlands eScience Center and Utrecht University

## License

Apache 2.0 License, see LICENSE file for license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.
![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
