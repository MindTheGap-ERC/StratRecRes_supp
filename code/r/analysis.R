#### Load required packages ####
library(ggplot2)

#### Load SMLBase and clean it ####
smlbase = read.csv(file = "data/smlbase/SMLBase_v1_04.csv", sep = ",")

df = data.frame(id = smlbase$boxModelID)

df$S_min = pmin(smlbase$SMin, smlbase$SMax, smlbase$S, na.rm = TRUE)
df$S_max = pmax(smlbase$SMin, smlbase$SMax, smlbase$S, na.rm = TRUE)

df$L_min = pmin(smlbase$L1, smlbase$L1Max, smlbase$L1Min, smlbase$L2, na.rm = TRUE)
df$L_max = pmax(smlbase$L1, smlbase$L1Max, smlbase$L1Min, smlbase$L2, na.rm = TRUE)

df$M_min = pmin(smlbase$M1, smlbase$M1Max, smlbase$M1Min, smlbase$M2, na.rm = TRUE)
df$M_max = pmax(smlbase$M1, smlbase$M1Max, smlbase$M1Min, smlbase$M2, na.rm = TRUE)

df$wd_min = pmin(smlbase$WDMin, smlbase$WDMax, smlbase$WD, na.rm = TRUE)
df$wd_max = pmax(smlbase$WDMin, smlbase$WDMax, smlbase$WD, na.rm = TRUE)
df$wd = 0.5 * (df$wd_max + df$wd_min)

df = dplyr::filter(df, !(df$S_min <= 0 | df$S_max <= 0 | df$M_min <= 0 | df$M_max <= 0| df$wd_max <= 0| df$wd_min <= 0))
df$S = 0.5 * (df$S_max + df$S_min)
df$M = 0.5 * (df$M_max + df$M_min)
df$L = 0.5 * (df$L_max + df$L_min)
df$Pe_min = df$L_min * df$S_min / df$M_max
df$Pe_max = df$L_max * df$S_max / df$M_min
df$Pe = df$L * df$S / df$M


dd = data.frame(df$S, df$L, df$M)
cor(dd)

#### Import data from matlab and determine dimless tavg ####
## load helper functions
source(file = "code/R/utils.R")
## import data from matlab
source(file = "code/R/import_from_matlab.R")
## load matlab data
load(file = "data/r_outputs/tavg_below_sml.Rdata")

## determine dimensionless time averaging below SML
tavg_dimless = rep(NA, length(peclet_numbers))
for (i in seq_along(tavg_list)){
  cdf = cumul_int_trap(x = t_dimless, y = tavg_list[[i]]$den, normalize = TRUE)
  tavg_dimless[i] = iqr_from_cdf(cdf, t_dimless)
}

plot(log10(peclet_numbers), tavg_dimless, type = "l", ylim = c(0, max(tavg_dimless)))

## reverse scaling for tavg and disorder

tavg_fun = function(Pe, S, L){
  tavg = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L/S
  return(tavg)
}
disorder_fun = function(Pe, L){
  disorder = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L
  return(disorder)
}

#### Determine tavg and disorder for SMLBase ####
df$tavg = tavg_fun(Pe = df$Pe, S = df$S, L = df$L)
df$disorder = disorder_fun(Pe = df$Pe, L = df$L)

df$tavg_max = tavg_fun(Pe = df$Pe_min, S = df$S_min, L = df$L_max)
df$tavg_min = tavg_fun(Pe = df$Pe_max,  S = df$S_max, L = df$L_min)
df$disorder_max = disorder_fun(Pe = df$Pe_min, L = df$L_max)
df$disorder_min = disorder_fun(Pe = df$Pe_max, L = df$L_min)

df$F_mix = approx(log10(peclet_numbers), tavg_dimless, xout = df$Pe, rule = 2)$y

hist(log10(df$tavg))
hist(df$disorder)

#### Empirical results from SMLBase ####
# negative sign to convert to mixing intensity instead of peclet number

signif(quantile(-log10(df$Pe), c(0.025, 0.975)), 3)
signif(median(-log10(df$Pe)), 3)
signif(quantile(-log10(df$Pe), c(0.25, 0.75)), 3)

signif(quantile(df$F_mix, c(0.25, 0.75), 3))
df$F_mix |> median() |> signif(3)
df$F_mix |> quantile() |> signif(3)

#### GLM ####

m1 = log10(df$M)
s1 = log10(df$S)
l1 = log10(df$L)
ta = log10(df$tavg)
di = log10(df$disorder)

tavg_glm = glm(ta ~ m1 + s1 + l1)
car::vif(tavg_glm)
disorder_glm = stats::glm(di ~ m1 + s1 + l1)
car::vif(disorder_glm)

tavg_glm
summary(tavg_glm)
step(tavg_glm)
step(disorder_glm)
#fl = visreg::visreg(tavg_glm,ylim=c(0,4.4))
rsq::rsq.partial(disorder_glm,adj=TRUE)


#### tavg and disorder vs. water depth ####
logtavg = log10(df$tavg)
wd_l = log10(df$wd)
lm_tavg_wd = lm( logtavg ~ wd_l)

di = log10(df$disorder)
lm_disorder_wd = lm(di ~ wd_l)

summary(lm_tavg_wd)
step(lm_tavg_wd)
step(lm_disorder_wd)
summary(lm_disorder_wd)

#### Figure 1 ####
# Data from Tomasovych et al 2018
to = read.csv("data/tomasovych_et_al_2018/S0094837318000222sup002.csv", sep = ";")
to2 = dplyr::filter(to, Station == "Po4 M21" & Outlier == "")

df_te = data.frame("age" = to2$Age[!is.na(to2$Age)], "depth" = as.numeric(to2$max..depth..cm.)[!is.na(to2$Age)])

age_lim = 150
depth_lim = 120


#ggsave("figs/tomasovych_et_al_2018_density.png", p1)

#df_te$depth
IQR(df_te$age[df_te$depth == 95])

##
# read data from matlab
da = R.matlab::readMat("code/matlab/Po4_res.mat")

ages = da$ages[1,]
depths = da$depths[1,]
u = da$u

#dimnames(u) <- list( "depths" = depths, "ages" = ages)
mydf = reshape2::melt(u)

p1 = ggplot(df_te, aes(x = age, y = depth)) + 
  geom_point() +
  xlim(0, age_lim) +
  scale_y_reverse(lim = c(depth_lim, 0)) +
  geom_density2d_filled(alpha = 0.5, show.legend = FALSE) +
  ggtitle("Core Po4 Particle distribution") +
  xlab("Age [years]") +
  ylab("Depth [cm]") + 
  ggtitle("Empirical") +
  scale_color_viridis_c()
p1
br = c(seq(0, 0.01, by = 0.001), seq(0.01, 0.1, by = 0.02))
 p2 = ggplot(mydf, aes(Var2, Var1, z = value)) + 
   geom_contour_filled(breaks = br, show.legend = FALSE) + 
   geom_contour(breaks = br)+
   xlim(range(mydf$Var2)) + 
   ylim(range(mydf$Var1)) +
   xlim(0, age_lim) +
   scale_y_reverse(lim = c(depth_lim, 0)) +
   xlab("Age [years]") +
   ylab("Depth [cm]") +
   ggtitle("Model") +
   scale_color_viridis_c()
p2 
ggsave("figs/Po4_modeled.png", p2)

jo = egg::ggarrange(p1, p2, ncol = 2, nrow = 1)

ggsave("figs/joint_Po4.png", jo)


#### GLM plot ####
sedr_label = expression( log[10]*"("*Sed*"."*~rate*")"*  ~  group("[",cm/a,"]") )
mix_label =  expression( log[10]* group("(", Biodiffusion*".",")")  ~  group("[",cm^2/a,"]") ) #"log10(Mixing int.) [cm^2/a]"
mix_depth_label =  expression( log[10]* group("(", Mixing~depth,")")  ~  group("[",cm,"]") )  #"log10(Mixing depth) [cm]"
tavg_label =  expression( log[10]* group("(", Time ~averaging,")")  ~  group("[",a,"]") ) # "log10(Time averaging) [a]"
disorder_label = expression( log[10]* group("(", Stratigraphic~disorder,")")  ~  group("[",cm,"]") )  #"log10(Stratigraphic disorder) [cm]"

tavg_y_axis_lims = range(log10(c(df$tavg* 1.3, df$tavg * 0.8)))
disorder_y_axis_lims = range(log10(c(1.1 * df$disorder, 0.9 * df$disorder)))

two_col_width_cm = 12.28
x_axis_text_size_pts = 7
y_axis_text_size_pts = 7
glm_plot_height_cm = 14
annot_size_pts = 7

make_glm_plot = function(fig_name){
  rsq = rsq::rsq.partial(tavg_glm)
  # sed rate
  vv = paste(": ", signif(rsq$partial.rsq[2], 3))
  annot = bquote(partial ~ R^2* .(vv))
  slope_annot = paste("Slope: ", signif(coefficients(tavg_glm)["s1"], 3))
  tavg_glm_plot_s = visreg::visreg(tavg_glm,xvar = "s1",
                                   gg = TRUE) +
    xlab(sedr_label) +
    ylab(tavg_label) +
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$tavg)) - 0.5, label = slope_annot, size = annot_size_pts / .pt)
  tavg_glm_plot_s
  
  # biodiffusion
  vv = paste(": ", signif(rsq$partial.rsq[1], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  tavg_glm_plot_m = visreg::visreg(tavg_glm,xvar = "m1",
                                   gg = TRUE) +
    xlab(mix_label) +
    ylab(tavg_label)  + 
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) 
  
  # mixing depth
  vv = paste(": ", signif(rsq$partial.rsq[3], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  tavg_glm_plot_l= visreg::visreg(tavg_glm,xvar = "l1",
                                  gg = TRUE) +
    xlab(mix_depth_label) +
    ylab(tavg_label)  + 
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) 
  
  rsq = rsq::rsq.partial(disorder_glm)
  # sedimentation rate
  vv = paste(": ", signif(rsq$partial.rsq[2], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  disorder_glm_plot_s = visreg::visreg(disorder_glm,xvar = "s1",
                                       gg = TRUE) +
    xlab(sedr_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) 
  
  # biodiffusion
  vv = paste(": ", signif(rsq$partial.rsq[1], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  disorder_glm_plot_m = visreg::visreg(disorder_glm,xvar = "m1",
                                       gg = TRUE) +
    xlab(mix_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) 
  
  # mixing depth
  vv = paste(": ", signif(rsq$partial.rsq[3], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  disorder_glm_plot_l= visreg::visreg(disorder_glm,xvar = "l1",
                                      gg = TRUE) +
    xlab(mix_depth_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) 
  
  figure = egg::ggarrange(tavg_glm_plot_s, tavg_glm_plot_l, tavg_glm_plot_m,
                          disorder_glm_plot_s, disorder_glm_plot_l, disorder_glm_plot_m,
                          nrow = 2, ncol = 3,
                          labels = LETTERS[1:6])
  
  ggsave(paste0("figs/", fig_name, ".png"), figure,
         width = two_col_width_cm, units = c("cm"),
         height = glm_plot_height_cm)
}

make_glm_plot("glm_res")


#### Figure 2 ####
# plots histogram of Peclet numbers and F_mix
plot_hist_and_fmix = function(file_name){
  aa = par(no.readonly = TRUE)
  
  pe_min_log = min(-log10(peclet_numbers))
  pe_max_log = max(-log10(peclet_numbers))
  pe_step = 1/3
  tavg_lwd = 3
  tavg_col = "black"
  a = hist(-log10(df$Pe), plot = FALSE, breaks = seq(pe_min_log, pe_max_log, pe_step))
  png(file_name,
      width = )
  par(mar = c(5.1, 4.1, 4.1, 4.1))
  plot(-log10(peclet_numbers), tavg_dimless, type = "l",
       xlab = "log10(Mixing Intensity)",
       ylab = "Dimensionless time-averaging",
       lwd = tavg_lwd,
       col = tavg_col,
       mar = c(5.1, 4.1, 4.1, 4.1))
  par(new = TRUE)
  
  plot(a, freq = FALSE, axes = FALSE,
       xlab = "",
       ylab = "",
       main = "")
  axis(4)
  mtext("Frequency", side = 4, line = 2.8)
  
  par(new = FALSE)
  dev.off()
}
plot_hist_and_fmix("figs/fig_2.png")


#### Supplementary Figure 1 ####
plot_fmix_dist = function(fig_name){
  ggplot(df, aes(x = F_mix)) + 
    geom_histogram() +
    xlab(expression(Distribution ~ of ~ F[mix]))
  ggsave(fig_name)
}
plot_fmix_dist("figs/supp_distribution_Fmix.png")

#### Supplementary figure 2 ####
plot_sml_histograms = function(fig_name){
  mix_label =  expression(  Biodiffusion*""  ~  group("[",cm^2/a,"]") )
  p_m = ggplot(df, aes(x = M)) + 
    geom_histogram() + 
    scale_x_log10() +
    xlab(mix_label)
  p_l = ggplot(df, aes(x = L)) + 
    geom_histogram() +
    xlab("Mixing depth [cm]")
  p_s =ggplot(df, aes(x = S)) +
    geom_histogram() +
    scale_x_log10() +
    xlab("log10(Sedimentation rate) [cm/a]")
  
  figure = egg::ggarrange(p_m, p_l, p_s, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
  ggsave(fig_name, figure)
}

plot_sml_histograms("figs/suppl_summary_SML.png")

#### tavg and disorder vs. water depth plot ####

wd_plot = function(filepath){
  lm1 = visreg::visreg(lm_tavg_wd, gg = TRUE) +
    xlab("log10(water depth) [m]") +
    ylab("log10(time averaging) [a]")
  lm2 = visreg::visreg(lm_disorder_wd, gg = TRUE) +
    xlab("log10(water depth) [m]") +
    ylab("Stratigraphic disorder [cm]")
  fig = egg::ggarrange(lm1, lm2, ncol = 2)
  ggsave(filepath, fig)
}

wd_plot("figs/Fig_4.png")
