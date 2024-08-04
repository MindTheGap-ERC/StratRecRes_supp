# Codebook

Describes the contents of the file `tavg_from_matlab.mat`.

## Contents

_.mat_ file containing the following variables (using matlab notation):

* `tavg_below_sml` : 2d array of size `length(t_dimless)`,  `length(peclet_numbers)`. `tavg_below_sml(i,k)` is the density of particles at time `t_dimless(i)` below the SML of the ADD with Peclet number `peclet_numbers(k)`
* `t_timless` : vector, dimensionless times at which the ADDs are evaluated
* `epsilon`: positive scalar, parameter characterizing the initial condition
* `odeoptions` : struct, settings for ode solver called by the pde solver `pdepe`
* `peclet_numbers` : vector, Peclet numbers at which the ADDs are calculated

## Usage

Can be opened in Matlab, is loaded into R via _code/r/import_adds.R_

## Origin

Either downloaded from Zenodo (INSERT DOI) via FILENAME HERE, or generated from scratch by running _code/matlab/create_adds.m_. See REPRODUCEME.md for details.
