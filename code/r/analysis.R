library(dplyr)

## load helper functions
source(file = "code/R/utils.R")

## import data from matlab
source(file = "code/R/import_from_matlab.R")

## load data
load(file = "data/r_outputs/tavg_below_sml.Rdata")

## determine dimensionless time averaging below SML
tavg_dimless = rep(NA, length(peclet_numbers))
for (i in seq_along(tavg_list)){
  cdf = cumul_int_trap(x = t_dimless, y = tavg_list[[i]]$den, normalize = TRUE)
  tavg_dimless[i] = iqr_from_cdf(cdf, t_dimless)
}

plot(log10(peclet_numbers), tavg_dimless, type = "l", ylim = c(0, max(tavg_dimless)))

## functions for tavg and disorder

tavg_fun = function(Pe, S, L){
  tavg = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L/S
  return(tavg)
}
disorder_fun = function(Pe, L){
  tavg = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L
  return(tavg)
}


## load empirical data
smlbase = read.csv(
  file = "data/smlbase_raw/SMLBase v1_03.csv",
  sep = ","
)

df = data.frame(id = smlbase$boxModelID)

df$S_min = pmin(smlbase$SMin, smlbase$SMax, smlbase$S, na.rm = TRUE)
df$S_max = pmax(smlbase$SMin, smlbase$SMax, smlbase$S, na.rm = TRUE)

df$L_min = pmin(smlbase$L1, smlbase$L1Max, smlbase$L1Min, smlbase$L2, na.rm = TRUE)
df$L_max = pmax(smlbase$L1, smlbase$L1Max, smlbase$L1Min, smlbase$L2, na.rm = TRUE)

df$M_min = pmin(smlbase$M1, smlbase$M1Max, smlbase$M1Min, smlbase$M2, na.rm = TRUE)
df$M_max = pmax(smlbase$M1, smlbase$M1Max, smlbase$M1Min, smlbase$M2, na.rm = TRUE)

df$Pe_min = df$L_min * df$S_min / df$M_max
df$Pe_max = df$L_max * df$S_max / df$M_min

df$S = 0.5 * (df$S_max + df$S_min)
df$M = 0.5 * (df$M_max + df$M_min)
df$L = 0.5 * (df$L_max + df$L_min)
df$Pe = df$L * df$S / df$M

df = dplyr::filter(df, !is.nan(df$Pe_min))
df = dplyr::filter(df, df$S != 0)

df$tavg_max = tavg_fun(Pe = df$Pe_min, S = df$S_min, L = df$L_max)
df$tavg_min = tavg_fun(Pe = df$Pe_max,  S = df$S_max, L = df$L_min)

df$tavg = tavg_fun(Pe = df$Pe, S = df$S, L = df$L)
df$disorder = disorder_fun(Pe = df$Pe, L = df$L)

df$disorder_max = disorder_fun(Pe = df$Pe_min, L = df$L_max)
df$disorder_min = disorder_fun(Pe = df$Pe_max, L = df$L_min)
hist(log10(df$tavg))
hist(log10(df$disorder))
hist(df$disorder)
