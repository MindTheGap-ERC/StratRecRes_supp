#### Load required packages ####
library(ggplot2)
library(dotwhisker)
library(gapminder)

#### Set seed ####
set.seed(296573) # generated via random.org

# threshold for multicollinearity 
vif_threshold = 3

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

tavg_glm = glm(ta ~  s1 +  m1 +l1)
disorder_glm = stats::glm(di ~ s1 +  m1 +l1)

step(tavg_glm)
step(disorder_glm)

# Write summary to csv file
x = summary(tavg_glm)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/tavg_glm_res.csv")

x = summary(disorder_glm)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/disorder_glm_res.csv")

## Check for multicolinearity
dd = data.frame( s1, m1, l1)
cor_mat = cor(dd)
# check variance inflation factors
vif_tavg = car::vif(tavg_glm)
vif_disorder = car::vif(disorder_glm)
if (any(vif_disorder> vif_threshold) | any(vif_tavg > vif_threshold)){ 
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

#### Normalized regression ####
data = data.frame(tavg_st = scale(ta),
                  disorder_st = scale(di),
                  m_st = scale(m1),
                  l_st = scale(l1),
                  s_st = scale(s1))

tavg_glm_std = glm(tavg_st ~ s_st + m_st  +  l_st, data = data)
disorder_glm_std = glm(disorder_st ~ s_st + m_st  +  l_st, data = data)

vif_tavg_std = car::vif(tavg_glm_std)
vif_disorder_std = car::vif(disorder_glm_std)
if (any(vif_disorder_std> vif_threshold) | any(vif_tavg_std > vif_threshold)){ 
  stop("Multicolinearity detected")  
}

step(tavg_glm_std)
step(disorder_glm_std)

# write summary to csv
x = summary(tavg_glm_std)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/tavg_glm_std_res.csv")

x = summary(disorder_glm_std)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/disorder_glm_std_res.csv")

rsq_part_disorder_std = rsq::rsq.partial(disorder_glm_std,adj=TRUE)
rsq_part_tavg_std = rsq::rsq.partial(tavg_glm_std, adj = TRUE)


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

plot_adms = function(){
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
  p = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("A", "B"))
  return(p)

}
ggsave(filename = "figs/adms_po_and_modeled.png",
       plot = plot_adms(),
       bg = "white")

#### Figure 2: histogram of mixing intensity and F_mix ####
fmix_label =  expression( F[mix] )  #"log10(Mixing depth) [cm]"
g_label =  expression( log[10]* group("(", G,")")  ~  group("[","-","]") )  #"log10(Mixing depth) [cm]"

plot_dimless_mixing_and_tavg = function(){
  aa = par(no.readonly = TRUE)
  pe_min_log = min(-log10(peclet_numbers))
  pe_max_log = max(-log10(peclet_numbers))
  pe_step = 1/3
  tavg_lwd = 3
  tavg_col = "black"
  a = hist(-log10(df$Pe), plot = FALSE, breaks = seq(pe_min_log, pe_max_log, pe_step))
  png("figs/dimless_mixing_and_tavg.png",
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
plot_dimless_mixing_and_tavg()

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

plot_glm_summary_fig = function(){
  ## time averaging
  rsq = rsq_part_tavg
  # sed rate
  vv = paste(": ", round(rsq$partial.rsq[2], 3))
  annot = as.expression(bquote(partial ~ R^2* .(vv)))
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
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
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
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
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
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
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
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
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
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
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
  
  figure = ggpubr::ggarrange(tavg_glm_plot_s,  tavg_glm_plot_m, tavg_glm_plot_l,
                          disorder_glm_plot_s,  disorder_glm_plot_m, disorder_glm_plot_l,
                          nrow = 2, ncol = 3,
                          labels = LETTERS[1:6])
  return(figure)
}


ggsave("figs/glm_summary_fig.png", 
       plot = plot_glm_summary_fig(),
       width = two_col_width_cm,
       units = c("cm"),
       height = glm_plot_height_cm,
       bg = "white")

#### figure: water depth vs tavg and disorder ####
signif_dig = 3
plot_wd_vs_params = function(){
  tavg_r2 = summary(lm_tavg_wd)$adj.r.squared |> signif(digits = signif_dig)
  tavg_slope = lm_tavg_wd$coefficients["wd_l"] |> unname() |> signif(digits = signif_dig)
  disorder_r2 = summary(lm_disorder_wd)$adj.r.squared |> signif(digits = signif_dig)
  disorder_slope = lm_disorder_wd$coefficients["wd_l"] |> unname() |> signif(digits = signif_dig)
  
  
  yrange = range(log10(df$tavg))
  vv = paste0(": ", tavg_r2)
  annot = as.expression(bquote(adj. ~ R^2*  .(vv)))
  slope_annot = paste("Slope:", tavg_slope)
  lm1 = visreg::visreg(lm_tavg_wd, gg = TRUE) +
    xlab(wd_label) +
    ylab(tavg_label) +
    annotate("text",
             y = 0.9 * max(yrange),
             x = mean(range(log10(df$wd))),
             label = annot) +
    annotate("text",
             x = mean(range(log10(df$wd))),
             y = 1 * max(yrange),
             label = slope_annot)
  
  yrange = range(log10(df$disorder))
  vv = paste0(": ", disorder_r2)
  annot = as.expression(bquote(adj. ~ R^2*  .(vv)))
  slope_annot = paste("Slope:", disorder_slope)
  lm2 = visreg::visreg(lm_disorder_wd, gg = TRUE) +
    xlab(wd_label) +
    ylab(disorder_label)+
    annotate("text",
             y = 0.9 * max(yrange),
             x = mean(range(log10(df$wd))),
             label = annot) +
    annotate("text",
             x = mean(range(log10(df$wd))),
             y = 1 * max(yrange),
             label = slope_annot)
  p = ggpubr::ggarrange(lm1, lm2, ncol = 2, nrow = 1,
                          labels = LETTERS[1:2])
  return(p)
}

ggsave(filename = "figs/water_depth_vs_params.png",
       plot = plot_wd_vs_params(),
       bg = "white")


#### Figure: Overview of SML parameters ####
plot_sml_overview_fig = function(){
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
  p_fmix = ggplot(df, aes(x = F_mix)) + 
    geom_histogram(bins = 30) +
    xlab(expression(F[mix] * " [-]")) +
    ylab("Count")
  
  p = ggpubr::ggarrange(p_s, p_m, p_l, p_fmix,
                          labels = LETTERS[1:4],
                          ncol = 2, nrow = 2)
  return(p)
}
ggsave(filename = "figs/sml_overview.png",
       plot = plot_sml_overview_fig(),
       bg = "white")

#### Coefficient plot ####

make_coeff_plots = function(){
  dot_args = list(size = 3)
  whisker_args = list(size = 1)
  y_axis_labels = labels = rev(c(
    expression(log[10](S)),
    expression(log[10](D[b])),
    expression(log[10](L))
  ))
  a1 = dwplot(tavg_glm,
              dot_args = dot_args,
              whisker_args = whisker_args) +
    ggtitle("Time-averaging") + 
    labs(x = "Regression coefficient", y = "Independent variables")   +
    scale_y_discrete(labels = y_axis_labels) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_text(angle = 45))
  a2 = dwplot(tavg_glm_std,
              dot_args = dot_args,
              whisker_args = whisker_args) +
    ggtitle("") +
    labs(x = "Beta coefficient", y = "Independent variables") +
    scale_y_discrete(labels = y_axis_labels)+
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_text(angle = 45))
  
  a3 = dwplot(disorder_glm,
              dot_args = dot_args,
              whisker_args = whisker_args) +
    ggtitle("Stratigraphic disorder") + 
    labs(x = "Regression coefficient", y = "Independent variables") +
    scale_y_discrete(labels = y_axis_labels)+
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_text(angle = 45))
  a4 = dwplot(disorder_glm_std,
              dot_args = dot_args,
              whisker_args = whisker_args) +
    ggtitle("") +
    labs(x = "Beta coefficient", y = "Independent variables") +
    scale_y_discrete(labels = y_axis_labels)+
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_text(angle = 45))
  p = ggpubr:: ggarrange(a1, a2, a3, a4,
                         ncol = 2,
                         nrow = 2,
                         labels = LETTERS[1:4])
  return(p)
}

ggsave(filename = "figs/coefficient_plots.png",
       plot = make_coeff_plots(),
       bg = "white")

#### Plot dimensionless ADM ####
# serves as basis for Fig 2?
plot_adm_dimless_sketch = function(){
  i = 25 # select one density of the many calculated in Matlab
  x = t_dimless
  y = tavg_list[[i]]$den
  #plot(x, y, type = "l")
  # specifies ADM dimensionless (assuming the SML is well mixed)
  # below the sml (d < 1) adm is specified by eq 3 main text
  f = function(a, d){
    z = approx(x = d - 1 + t_dimless- 0.3, y = y, xout = a, rule = 2)$y
    z[d < 1] = approx(x = t_dimless - 0.3, y = y, xout = a, rule = 2)$y
    return(z)
  }
  
  a = seq(0, 3, by = 0.05)
  d = seq(0, 3, by = 0.05)
  df = expand.grid(a, d)
  df$z = mapply(f, df$Var1, df$Var2)
  p = ggplot(df, aes(x = Var1, y = Var2, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient(low = "black", # linear white scale from min to max
                        high = "white",
                        trans = "sqrt") + # nonlinear trans to buff low values
    scale_y_reverse() +
    labs(x = "Dimensionless Age [-]",
         y = "Dimensionless Depth [-]") +
    theme(legend.position = "none")
  return(p)
}

ggsave(filename = "figs/adm_dimless_whitescale.tiff",
       plot = plot_adm_dimless_sketch())

plot_dimless_add = function(){
  d_select = 1.2
  a_select = 2.8
  a = seq(0, 4, by = 0.05)
  d = seq(0, 4, by = 0.05)
  
  i = 25 # select one density of the many calculated in Matlab
  x = t_dimless
  y = tavg_list[[i]]$den
  #plot(x, y, type = "l")
  # specifies ADM dimensionless (assuming the SML is well mixed)
  # below the sml (d < 1) adm is specified by eq 3 main text
  f = function(a, d){
    z = approx(x = d - 1 + t_dimless- 0.3, y = y, xout = a, rule = 2)$y
    z[d < 1] = approx(x = t_dimless - 0.3, y = y, xout = a, rule = 2)$y
    return(z)
  }
  
  fixed_depth_color = "black"
  fixed_age_color = "red"
  transect_width = 2
  df = expand.grid(a, d)
  df$z = mapply(f, df$Var1, df$Var2)
  p1 = ggplot(df, aes(x = Var1, y = Var2, z = z)) +
    geom_contour_filled(breaks = sort(seq(0,1.1, by = 0.1))) +
    geom_contour(aes(z = z), color = "black", linewidth = 0.4,
                 breaks = sort(seq(0,1.1, by = 0.1))) +
    labs(x = "Age [-]",
         y = "Depth [-]") +
    annotate("segment",
             x = min(a), xend = max(a),
             y = d_select, yend = d_select,
             linewidth =transect_width,
             color = fixed_depth_color) +  
    annotate("segment",
             x = a_select, xend = a_select,
             y = min(d), yend = max(d),
             linewidth = transect_width,
             color = fixed_age_color) +
    coord_cartesian(clip = "on") +
    scale_y_reverse() +
    theme(legend.position = "none")
  
  
  df1 = data.frame(d = d, den_vals =  mapply(f, a_select, d))
  p2 = df1 |> ggplot(aes(x =  d, y = den_vals)) +
    geom_line(color = fixed_age_color,
              linewidth = transect_width) +
    coord_flip() +
    scale_x_reverse() +
    labs(x = "Depth [-]",
         y = "Density")
  
  df2 = data.frame(a = a, den_vals = mapply(f, a, d_select))
  p3 = df2 |> ggplot(aes(x = a, y = den_vals)) +
    geom_line(color = fixed_depth_color,
              linewidth = transect_width) +
    labs(x = "Age [-]",
         y = "Density")
  
  p4 = NULL
  
  p = ggpubr::ggarrange(p1, p2, p3, p4,
                        nrow = 2, 
                        ncol = 2,
                        widths = c(2,1),
                        heights = c(2,1))
  return(p)
}



ggsave(filename = "figs/dimless_adm.png",
       plot = plot_dimless_add(),
       height = 7,
       width = 7,
       bg = "white")


#### Save results ####
re.df = data.frame(description = names(results), values = unname(results))
write.csv(re.df, file = "data/res/results.csv")

cat("Done. Results are under \"figs/\"and \"data/res/results.csv\"\n")
