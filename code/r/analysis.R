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


#### load empirical data from SMLBase ####
smlbase = read.csv(
  file = "data/smlbase/SMLBase v1_03.csv.csv",
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

df = dplyr::filter(df, !(df$S_min <= 0 | df$S_max <= 0 | df$M_min <= 0 | df$M_max <= 0))

df$S = 0.5 * (df$S_max + df$S_min)
df$M = 0.5 * (df$M_max + df$M_min)
df$L = 0.5 * (df$L_max + df$L_min)
df$Pe = df$L * df$S / df$M


signif(quantile(log10(df$Pe), c(0.025, 0.975)), 3)

signif(median(log10(df$Pe)), 3)
signif(quantile(log10(df$Pe), c(0.25, 0.75)), 3)

df$tavg_max = tavg_fun(Pe = df$Pe_min, S = df$S_min, L = df$L_max)
df$tavg_min = tavg_fun(Pe = df$Pe_max,  S = df$S_max, L = df$L_min)

df$tavg = tavg_fun(Pe = df$Pe, S = df$S, L = df$L)
df$disorder = disorder_fun(Pe = df$Pe, L = df$L)

df$disorder_max = disorder_fun(Pe = df$Pe_min, L = df$L_max)
df$disorder_min = disorder_fun(Pe = df$Pe_max, L = df$L_min)
hist(log10(df$tavg))
hist(log10(df$disorder))
hist(df$disorder)

#### Load data from Matlab ####

Pe_matlab = peclet_numbers
tavg_dless_matlab = tavg_dimless

Pe_min_log = min(log10(Pe_matlab))
Pe_max_log = max(log10(Pe_matlab))
Pe_step = 1/3

a = hist(log10(df$Pe), plot = FALSE, freq = FALSE, breaks = seq(Pe_min_log, Pe_max_log, Pe_step))

plot(a, add = TRUE)

plot(a)
lines(log10(pe_matlab), tavg_dless_matlab)
inc_density = 2
# 
# plot(NULL,
#      xlim = range(log10(Pe_matlab)),
#      ylim = c(0, max(tavg_dless_matlab)) * 1.1,
#      xlab = "Log10(Peclet number)",
#      ylab = "")
tavg_lwd = 3
tavg_col = "black"
png("figs/fig2.png")
plot(log10(pe_matlab), tavg_dless_matlab, type = "l",
     xlab = "log10(Peclet number)",
     ylab = "Dimensionless time-averaging",
     lwd = tavg_lwd,
     col = tavg_col)
par(new = TRUE)

plot(a, freq = FALSE, axes = FALSE,
     xlab = "",
     ylab = "",
     main = "")
axis(4)
mtext("Frequency", side = 4)

par(new = FALSE)
dev.off()

#### distribution of F_mix ####

F_mix = approxfun(log10(Pe_matlab), tavg_dless_matlab)
png("figs/supp_fig_1.png")
hist(F_mix(df$Pe))
dev.off()

#### Data from Tomasovych et al 2018

to = read.csv("data/tomasovych_et_al_2018/S0094837318000222sup002.csv", sep = ";")
unique(to$Station)
to2 = dplyr::filter(to, Station == "Po3 M13")

plot(to2$Age, to2$max..depth..cm.)
 df_te = data.frame("age" = to2$Age[!is.na(to2$Age)], "depth" = as.numeric(to2$max..depth..cm.)[!is.na(to2$Age)])
 library(ggplot2)
ggplot(df_te, aes(x = age, y = depth)) + 
  geom_point() +
  scale_y_reverse() +
  geom_density2d_filled(alpha = 0.5)
ggsave("figs/tomasovych_et_al_2018_density.png")

df_te$depth
IQR(df_te$age[df_te$depth == 100])
