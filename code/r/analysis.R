#### Load required packages ####
library(ggplot2)

#### Set seed ####
set.seed(42)

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
df$Pe = df$L * df$S / df$M # peclet number
df$G = 1/df$Pe # mixing intensity

#### Import data from matlab and determine dimless tavg ####
## load helper functions
source(file = "code/R/utils.R")
## import data from matlab, saves intermediate file under data/r_outputs/matlab_res.Rdata
source(file = "code/R/import_from_matlab.R")
## load matlab data
load(file = "data/r_outputs/matlab_res.Rdata")

## determine dimensionless time averaging below SML
tavg_dimless = rep(NA, length(peclet_numbers))
for (i in seq_along(tavg_list)){
  # cumulative distribution function of particle ages 
  cdf = cumul_int_trap(x = t_dimless, y = tavg_list[[i]]$den, normalize = TRUE)
  # determine IQR 
  tavg_dimless[i] = iqr_from_cdf(cdf, t_dimless) # corresponds to F_mix at peclet numbers
}

#plot(log10(peclet_numbers), tavg_dimless, type = "l", ylim = c(0, max(tavg_dimless)))

## reverse scaling for tavg and disorder

tavg_fun = function(Pe, S, L){
  tavg = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L/S
  return(tavg)
}
disorder_fun = function(Pe, L){
  disorder = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L
  return(disorder)
}

#### Determine F_mix, tavg and disorder for SMLBase ####
df$tavg = tavg_fun(Pe = df$Pe, S = df$S, L = df$L)
df$disorder = disorder_fun(Pe = df$Pe, L = df$L)
df$tavg_max = tavg_fun(Pe = df$Pe_min, S = df$S_min, L = df$L_max)
df$tavg_min = tavg_fun(Pe = df$Pe_max,  S = df$S_max, L = df$L_min)
df$disorder_max = disorder_fun(Pe = df$Pe_min, L = df$L_max)
df$disorder_min = disorder_fun(Pe = df$Pe_max, L = df$L_min)
df$F_mix = approx(log10(peclet_numbers), tavg_dimless, xout = df$Pe, rule = 2)$y

#hist(log10(df$tavg))
#hist(df$disorder)

#### Empirical results from SMLBase ####
results = c()
round_res = 3
results["median log10 G"] = round(median(log10(df$G)), round_res)
results["1st quartile log10 G"] = round(quantile(log10(df$G), p = 0.25), round_res)
results["3rd quartile log10 G"] = round(quantile(log10(df$G), p = 0.75), round_res)
results["1st quartile log10 G"] = round(quantile(log10(df$G), p = 0.25), round_res)
results["lower bound HDR log10 G"] = round(quantile(log10(df$G), p = 0.025), round_res)
results["upper bound HDR log10 G"] = round(quantile(log10(df$G), p = 1 - 0.025), round_res)
results["median F_mix"] = round(median(df$F_mix), round_res)

#### GLMs ####
m1 = log10(df$M)
s1 = log10(df$S)
l1 = log10(df$L)
ta = log10(df$tavg)
di = log10(df$disorder)

tavg_glm = glm(ta ~ m1 + s1 + l1)
disorder_glm = stats::glm(di ~ m1 + s1 + l1)

## Check for multicolinearity
dd = data.frame(m1, s1, l1)
cor_mat = cor(dd)
# check variance inflation factors
vif_tavg = car::vif(tavg_glm)
vif_disorder = car::vif(disorder_glm)
if (any(vif_disorder> 5) | any(vif_tavg > 5)){ 
  stop("Multicolinearity detected")  
}

## Partial R square
rsq_part_disorder = rsq::rsq.partial(disorder_glm,adj=TRUE)
rsq_part_tavg = rsq::rsq.partial(tavg_glm, adj = TRUE)

## Save results for Rsq and coefficients
results[c("partial R2 disorder D_b", "partial R2 disorder S", "partial R2 disorder L")] = round(rsq_part_disorder$partial.rsq, round_res)
results[c("coeff disorder D_b", "coeff disorder S", "coeff disorder L")] = round(disorder_glm$coefficients[-1], round_res)
results[c("partial R2 tavg D_b", "partial R2 tavg S", "partial R2 tavg L")] = round(rsq_part_tavg$partial.rsq, round_res)
results[c("coeff tavg D_b", "coeff tavg S", "coeff tavg L")] = round(tavg_glm$coefficients[-1], round_res)

#### LM water depth vs tavg and disorder ####
logtavg = log10(df$tavg)
wd_l = log10(df$wd)
logdi = log10(df$disorder)
lm_tavg_wd = lm( logtavg ~ wd_l)
lm_disorder_wd = lm(logdi ~ wd_l)

#### Figure 1: ADD in Po4 core ####
# plotting limits
age_lim = 150
depth_lim = 120

# Data from Tomasovych et al 2018
to = read.csv("data/tomasovych_et_al_2018/S0094837318000222sup002.csv", sep = ";")
to2 = dplyr::filter(to, Station == "Po4 M21" & Outlier == "")
df_te = data.frame("age" = to2$Age[!is.na(to2$Age)], "depth" = as.numeric(to2$max..depth..cm.)[!is.na(to2$Age)])
results["IQR 110 cm depth Po4"] = IQR(df_te$age[df_te$depth == 110])
results["deepest C. gibba shell younger than 25 years in Po4"] = max(df_te$depth[df_te$age <= 25])

make_fig_1 = function(){
  # remove samples deeper and older than plot boundaries
  df_te = dplyr::filter(df_te, depth <= depth_lim & age <= age_lim)
  p1 = ggplot(df_te, aes(x = age, y = depth)) + 
    geom_point() +
    xlim(0, age_lim) +
    scale_y_reverse(lim = c(depth_lim, 0)) +
    geom_density2d_filled(alpha = 0.5, show.legend = FALSE) +
    ggtitle("Particle Distribution Core Po4") +
    xlab("Age [a]") +
    ylab("Depth [cm]") +
    scale_color_viridis_c()
  #p1
  u2 = u[ 1:depth_lim +1, 1:age_lim + 1]
  mydf = reshape2::melt(u2)
  br = c(seq(0, 0.01, by = 0.001), seq(0.01, 0.1, by = 0.02), Inf) # breaks for density
  p2 = ggplot(mydf, aes(Var2, Var1, z = value)) + 
    geom_contour_filled(breaks = br, show.legend = FALSE, alpha = 0.8) + 
    geom_contour(breaks = br)+
    xlim(0, age_lim) +
    scale_y_reverse(lim = c(depth_lim, 0)) +
    xlab("Age [a]") +
    ylab("Depth [cm]") +
    ggtitle("Modelled Age-Depth Distribution") +
    scale_color_viridis_c()
  #p2 
  #ggsave("figs/Po4_modeled.png", p2)
  jo = egg::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("A", "B"), draw = FALSE)
  ggsave("figs/Figure1.png", jo)
}
make_fig_1()

#### Figure 2: histogram of mixing intensity and F_mix ####
fmix_label =  expression( F[mix] )  #"log10(Mixing depth) [cm]"
g_label =  expression( log[10]* group("(", G,")")  ~  group("[","-","]") )  #"log10(Mixing depth) [cm]"


make_fig_2 = function(){
  aa = par(no.readonly = TRUE)
  pe_min_log = min(-log10(peclet_numbers))
  pe_max_log = max(-log10(peclet_numbers))
  pe_step = 1/3
  tavg_lwd = 3
  tavg_col = "black"
  a = hist(-log10(df$Pe), plot = FALSE, breaks = seq(pe_min_log, pe_max_log, pe_step))
  png("figs/Figure2.png",
      width = )
  par(mar = c(5.1, 4.1, 4.1, 4.1))
  plot(-log10(peclet_numbers), tavg_dimless, type = "l",
       xlab = g_label,
       ylab = fmix_label,
       lwd = tavg_lwd,
       col = tavg_col,
       ylim = c(0, 1.1 * max(tavg_dimless)),
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
make_fig_2()

#### Figure 3: GLM plot ####
sedr_label = expression( log[10]*"("*Sed*"."*~rate*")"*  ~  group("[",cm/a,"]") )
mix_label =  expression( log[10]* group("(", Biodiffusion*"",")")  ~  group("[",cm^2/a,"]") ) #"log10(Mixing int.) [cm^2/a]"
mix_depth_label =  expression( log[10]* group("(", Mixing~Depth,")")  ~  group("[",cm,"]") )  #"log10(Mixing depth) [cm]"
tavg_label =  expression( log[10]* group("(", Time-averaging,")")  ~  group("[",a,"]") ) # "log10(Time averaging) [a]"
disorder_label = expression( log[10]* group("(", Stratigraphic~Disorder,")")  ~  group("[",cm,"]") )  #"log10(Stratigraphic disorder) [cm]"
wd_label =  expression( log[10]* group("(", Water~Depth,")")  ~  group("[",m,"]") )  #"log10(Mixing depth) [cm]"

tavg_y_axis_lims = range(log10(c(df$tavg* 1.3, df$tavg * 0.8)))
disorder_y_axis_lims = range(log10(c(1.1 * df$disorder, 0.9 * df$disorder)))

two_col_width_cm = 12.28
x_axis_text_size_pts = 7
y_axis_text_size_pts = 7
glm_plot_height_cm = 14
annot_size_pts = 7

make_fig_3 = function(){
  ## time averaging
  rsq = rsq_part_tavg
  # sed rate
  vv = paste(": ", round(rsq$partial.rsq[2], 3))
  annot = bquote(partial ~ R^2* .(vv))
  slope_annot = paste("Slope: ", round(coefficients(tavg_glm)["s1"], 4))
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
  vv = paste(": ", round(rsq$partial.rsq[1], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  slope_annot = paste("Slope: ", round(coefficients(tavg_glm)["m1"], 3))
  tavg_glm_plot_m = visreg::visreg(tavg_glm,xvar = "m1",
                                   gg = TRUE) +
    xlab(mix_label) +
    ylab(tavg_label)  + 
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$tavg)) - 0.5, label = slope_annot, size = annot_size_pts / .pt)
  
  # mixing depth
  vv = paste(": ", round(rsq$partial.rsq[3], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  slope_annot = paste("Slope: ", round(coefficients(tavg_glm)["l1"], 3))
  tavg_glm_plot_l= visreg::visreg(tavg_glm,xvar = "l1",
                                  gg = TRUE) +
    xlab(mix_depth_label) +
    ylab(tavg_label)  + 
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$tavg)) - 0.5, label = slope_annot, size = annot_size_pts / .pt)
  
  ## disorder
  rsq = rsq_part_disorder
  # sedimentation rate
  vv = paste(": ", round(rsq$partial.rsq[2], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  slope_annot = paste("Slope:", round(coefficients(disorder_glm)["s1"], 3))
  disorder_glm_plot_s = visreg::visreg(disorder_glm,xvar = "s1",
                                       gg = TRUE) +
    xlab(sedr_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$disorder)) - 0.2, label = slope_annot, size = annot_size_pts / .pt)
  
  # biodiffusion
  vv = paste(": ", round(rsq$partial.rsq[1], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  slope_annot = paste("Slope:", round(coefficients(disorder_glm)["m1"], 3))
  disorder_glm_plot_m = visreg::visreg(disorder_glm,xvar = "m1",
                                       gg = TRUE) +
    xlab(mix_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$disorder)) - 0.2, label = slope_annot, size = annot_size_pts / .pt)
  
  # mixing depth
  vv = paste(": ", round(rsq$partial.rsq[3], 3))
  annot = bquote(partial ~ R^2*  .(vv))
  slope_annot = paste("Slope:", round(coefficients(disorder_glm)["l1"], 3))
  disorder_glm_plot_l= visreg::visreg(disorder_glm,xvar = "l1",
                                      gg = TRUE) +
    xlab(mix_depth_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$disorder)) - 0.2, label = slope_annot, size = annot_size_pts / .pt)
  
  figure = egg::ggarrange(tavg_glm_plot_s, tavg_glm_plot_l, tavg_glm_plot_m,
                          disorder_glm_plot_s, disorder_glm_plot_l, disorder_glm_plot_m,
                          nrow = 2, ncol = 3,
                          labels = LETTERS[1:6],
                          draw = FALSE)
  ggsave("figs/Figure3.png", figure,
         width = two_col_width_cm, units = c("cm"),
         height = glm_plot_height_cm)
}

make_fig_3()

#### Figure 4: water depth vs tavg and disorder ####
make_fig_4 = function(){
  lm1 = visreg::visreg(lm_tavg_wd, gg = TRUE) +
    xlab(wd_label) +
    ylab(tavg_label)
  lm2 = visreg::visreg(lm_disorder_wd, gg = TRUE) +
    xlab(wd_label) +
    ylab(disorder_label)
  fig = egg::ggarrange(lm1, lm2, ncol = 2, draw = FALSE, labels = c("A", "B"))
  ggsave("figs/Figure4.png", fig)
}
make_fig_4()

#### Supplementary Figure 1 ####
make_supp_fig_1 = function(){
  ggplot(df, aes(x = F_mix)) + 
    geom_histogram(bins = 30) +
    xlab(expression(Distribution ~ of ~ F[mix])) +
    ylab("Count")
  ggsave("figs/supp_Figure1.png")
}
make_supp_fig_1()

#### Supplementary figure 2 ####
make_supp_fig_2 = function(){
  mix_label =  expression( log[10] (  Biodiffusion)*""  ~  group("[",cm^2/a,"]") )
  p_m = ggplot(df, aes(x = M)) + 
    geom_histogram(bins = 30) + 
    scale_x_log10() +
    xlab(mix_label) +
    ylab("Count")
  p_l = ggplot(df, aes(x = L)) + 
    geom_histogram(bins = 30) +
    xlab("Mixing Depth [cm]") +
    ylab("Count")
  p_s =ggplot(df, aes(x = S)) +
    geom_histogram(bins = 30) +
    scale_x_log10() +
    xlab("log10(Sedimentation Rate) [cm/a]") +
    ylab("Count")
  
  figure = egg::ggarrange(p_m, p_l, p_s, labels = c("A", "B", "C"), ncol = 2, nrow = 2,
                          draw = FALSE)
  ggsave("figs/supp_Figure2.png", figure)
}
make_supp_fig_2()

#### Save results ####
re.df = data.frame(description = names(results), values = unname(results))
write.csv(re.df, file = "data/res/results.csv")

cat("Done. Results are under \"figs/\"and \"data/res/results.csv\"\n")
