#### Load required packages ####
library(ggplot2)
library(dotwhisker)
library(gapminder)
library(patchwork)
library(metR)
library(dplyr)
library(visreg)

#### Plotting constants ####
## Colors
col_not_saturated = "goldenrod"
col_saturated = "mediumpurple"
sat_cols = c("Not saturated" = col_not_saturated, "Saturated" = col_saturated)
col_ribbon = "grey85"
col_hist_fill = "grey80"
col_hist_border = "grey40"
col_contour = "grey30"
col_line = "black"
col_low_risk = "aquamarine"
col_high_risk = "tomato"
fixed_depth_color = "black"
fixed_age_color = "red"

## Point shapes
sat_shapes = c("Not saturated" = 17, "Saturated" = 16) # filled triangle, filled circle

## Transparency
region_alpha = 0.25
box_alpha = 0.6
sat_jitter_alpha = 0.9
logit_point_alpha = 0.8
logit_ribbon_alpha = 0.2
density_alpha = 0.5
contour_fill_alpha = 0.8
risk_fill_alpha = 0.3

## Text sizes (pts)
title_size_pts = 10
axis_title_size_pts = 8
axis_text_size_pts = 7
x_axis_text_size_pts = 7
y_axis_text_size_pts = 7
annot_size_pts = 7
annot_size_large_pts = 8

## Line widths and point sizes
fit_linewidth = 1
contour_linewidth = 0.4
highlight_contour_linewidth = 1.2
transect_width = 2
sat_jitter_size = 0.8
logit_point_size = 1
coef_dot_size = 3
coef_whisker_size = 1

## Figure widths (cm)
single_col_width_cm = 8.5
double_col_width_cm = 19

## Axis labels
fmix_label = expression( "Dimless time-averaging/disorder" ~ F[mix] ~ " [-]" )
g_label = expression( Mixing~intensity~G  ~  group("[","-","]") )
# G and all S, D_b, L, tavg, disorder, SSD axes are plotted on log10 coordinates; these place
# breaks at powers of ten (or at 1, 2, 5 x 10^k for axes spanning less than ~2 decades)
# and label them with ordinary numbers
log10_breaks = function(lims) {
  in_range = function(b) b[b >= lims[1] - 1e-9 & b <= lims[2] + 1e-9]
  k = seq(floor(lims[1]), ceiling(lims[2]))
  b = in_range(k)
  if (length(b) < 3) b = in_range(as.vector(outer(log10(c(1, 2, 5)), k, "+")))
  b
}
# ordinary numbers, except compact powers of ten (e.g. 10^3) beyond 0.01 to 100 to avoid cramped labels
log10_labels = function(b) {
  lab = vapply(b, function(x) {
    if (abs(x - round(x)) < 1e-9 && abs(x) > 2) paste0("10^", round(x))
    else format(signif(10^x, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  }, character(1))
  parse(text = lab)
}
log10_x_scale = function(...) scale_x_continuous(breaks = log10_breaks, labels = log10_labels, ...)
log10_y_scale = function(...) scale_y_continuous(breaks = log10_breaks, labels = log10_labels, ...)
sedr_label = expression( "Sed. acc." ~ group("[", cm/a, "]") )
mix_label = expression( Biodiffusion ~ group("[",cm^2/a,"]") )
mix_depth_label = expression( Mixing~depth ~ group("[",cm,"]") )
tavg_label = expression( "Time-averaging" ~ group("[",a,"]") )
disorder_label = expression( Strat.~disorder ~ group("[",cm,"]") )
wd_lin_label = "Water depth [m]"
ssd_log_label = expression( Safe~sampling~dist. ~ group("[",cm,"]") )
ssd_log_label_short = expression( Safe~sampling~dist. ~ group("[",cm,"]") )

#### Set seed ####
set.seed(296573) # generated via random.org

# threshold for multicollinearity 
vif_threshold = 3
signif_dig = 3

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
df$F_mix = approx(peclet_numbers, tavg_dimless, xout = df$Pe, rule = 2)$y

#hist(log10(df$tavg))
#hist(df$disorder)

#### Empirical results from SMLBase ####

F_mix_threshold = 0.9 * max(tavg_dimless)
ff = function(x) approx(x = -log10(peclet_numbers), y = tavg_dimless - F_mix_threshold, xout = x)$y
G_threshold = 10^uniroot(ff, interval = c(-2,2))$root 
# determine which locations are (not) saturated
df$saturated = df$G > G_threshold
df$saturated[df$saturated] = "Saturated"
df$saturated[df$saturated == "FALSE"] = "Not saturated"
df <- df |> mutate(sat = as.integer(saturated == "Saturated"))
df = df |> mutate(S_log = log10(S))
df = df |> mutate(M_log = log10(M))
df = df |> mutate(L_log = log10(L))
df = df |> mutate(wd_log = log10(wd))
df = df |> mutate(tavg_log = log10(tavg))
df = df |> mutate(disorder_log = log10(disorder))
## Determine minimum safe sampling distance

get_min_safe_sampling_distance = function(){
  lim = c()
  for (i in seq_along(peclet_numbers)){
    fval = tavg_list[[i]]$den
    pe = peclet_numbers[i]
    
    aa = cumul_int_trap(t_dimless, fval, normalize = TRUE)
    den = c(0, diff(aa)/diff(t_dimless))
    CDF = approxfun(x = t_dimless, y = aa, yleft = 0, yright = 1)
    
    PDF = approxfun(x = t_dimless, y = den, yleft = 0, yright = 0)
    
    comb_f = function(y){
      comb_fun = function(x){PDF(x)* CDF(x-y)}
      return(integrate(comb_fun, lower = - Inf, upper = Inf)[[1]] - 0.1)
    }
    lim = c(lim, uniroot(comb_f, lower = 0, upper = 3)$root)
    
  }
  plot(-log10(peclet_numbers), lim, type = "l",
       main = "Distance to ensure orderedness",
       xlab = "log10 G",
       ylab = "Proportion of mixing depth")
  
  lim_fun = approxfun(-log10(peclet_numbers), lim)
  safe_sampling_distance = lim_fun(log10(df$G)) * df$L
  return(safe_sampling_distance)
}
df$ssd = get_min_safe_sampling_distance()
df = df |> mutate(ssd_log = log10(ssd))

results = c()
round_res = 3
signif_digits = 3
results["median log10 G"] = round(median(log10(df$G)), round_res)
results["1st quartile log10 G"] = round(quantile(log10(df$G), p = 0.25), round_res)
results["3rd quartile log10 G"] = round(quantile(log10(df$G), p = 0.75), round_res)
results["1st quartile log10 G"] = round(quantile(log10(df$G), p = 0.25), round_res)
results["lower bound HDR log10 G"] = round(quantile(log10(df$G), p = 0.025), round_res)
results["upper bound HDR log10 G"] = round(quantile(log10(df$G), p = 1 - 0.025), round_res)
results["median F_mix"] = round(median(df$F_mix), round_res)
results["prop of saturated SMLs"] = c(sum(df$G > G_threshold) / length(df$id)) |>
  (\(x) x * 100)() |> signif(digits = 3)



results["F_mix saturation threshold"] = F_mix_threshold |> signif(digits = signif_digits)
results["G saturation threshold"] = G_threshold |> signif(digits = signif_digits)
#### GLMs ####

p_from_glm = function(x){
  d = drop1(x, test = "F")
  return(setNames(d$`Pr(>F)`, rownames(d))[-1])   # drop the "<none>" row
}
m1 = log10(df$M)
s1 = log10(df$S)
l1 = log10(df$L)
ta = log10(df$tavg)
di = log10(df$disorder)

tavg_glm = glm(ta ~  s1 +  m1 +l1)
disorder_glm = stats::glm(di ~ s1 +  m1 +l1)
# rows with NA/Inf SSD (e.g. G outside the simulated Peclet range) are dropped;
# group labels are carried along so they stay aligned with the residuals in the plots
ssd_dat = data.frame(ssd = df$ssd_log,
                     s1 = s1, m1 = m1, l1 = l1,
                     group = df$saturated) |>
  dplyr::filter(is.finite(ssd))
ssd_glm = glm(ssd ~ s1 + m1 + l1, data = ssd_dat)

step(tavg_glm)
step(disorder_glm)
step(ssd_glm)

# Write summary to csv file
x = summary(tavg_glm)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/tavg_glm_res.csv")

x = summary(disorder_glm)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/disorder_glm_res.csv")

x = summary(ssd_glm)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/ssd_glm_res.csv")

## Check for multicolinearity
dd = data.frame( s1, m1, l1)
cor_mat = cor(dd)
# check variance inflation factors
vif_tavg = car::vif(tavg_glm)
vif_disorder = car::vif(disorder_glm)
vif_ssd = car::vif(ssd_glm)
if (any(vif_disorder> vif_threshold) | any(vif_tavg > vif_threshold) | any(vif_ssd > vif_threshold)){ 
  stop("Multicolinearity detected")  
}

## Partial R square
rsq_part_disorder = rsq::rsq.partial(disorder_glm,adj=TRUE)
rsq_part_tavg = rsq::rsq.partial(tavg_glm, adj = TRUE)
rsq_part_ssd = rsq::rsq.partial(ssd_glm, adj = TRUE)

## Save results for Rsq and coefficients, matched by variable name rather than position
save_glm_results = function(name, model, rsq_obj){
  pr2 = setNames(rsq_obj$partial.rsq, rsq_obj$variable)
  results[paste("partial R2", name, c("S", "D_b", "L"))] <<- round(pr2[c("s1", "m1", "l1")], round_res)
  results[paste("coeff",      name, c("S", "D_b", "L"))] <<- round(coef(model)[c("s1", "m1", "l1")], round_res)
}
save_glm_results("disorder", disorder_glm, rsq_part_disorder)
save_glm_results("tavg",     tavg_glm,     rsq_part_tavg)
save_glm_results("ssd",      ssd_glm,      rsq_part_ssd)

#### Normalized GLM regression ####
data = data.frame(tavg_st = scale(ta),
                  disorder_st = scale(di),
                  m_st = scale(m1),
                  l_st = scale(l1),
                  s_st = scale(s1))

tavg_glm_std = glm(tavg_st ~ s_st + m_st  +  l_st, data = data)
disorder_glm_std = glm(disorder_st ~ s_st + m_st  +  l_st, data = data)
ssd_glm_std = glm(scale(log10(ssd)) ~ scale(S_log) + scale(M_log) + scale(L_log), data = df)

vif_tavg_std = car::vif(tavg_glm_std)
vif_disorder_std = car::vif(disorder_glm_std)
vif_ssd_std = car::vif(ssd_glm_std)
if (any(vif_disorder_std> vif_threshold) | any(vif_tavg_std > vif_threshold) | any(vif_ssd_std > vif_threshold)){ 
  stop("Multicolinearity detected")  
}

step(tavg_glm_std)
step(disorder_glm_std)
step(ssd_glm_std)

# write summary to csv
x = summary(tavg_glm_std)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/tavg_glm_std_res.csv")

x = summary(disorder_glm_std)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/disorder_glm_std_res.csv")

x = summary(ssd_glm_std)
coeff = as.data.frame(x$coefficients)
coeff[,1:3] = round(coeff[,1:3], 3)
write.csv(coeff, file = "data/res/ssd_glm_std_res.csv")

rsq_part_disorder_std = rsq::rsq.partial(disorder_glm_std,adj=TRUE)
rsq_part_tavg_std = rsq::rsq.partial(tavg_glm_std, adj = TRUE)
rsq_part_ssd_std = rsq::rsq.partial(ssd_glm_std, adj = TRUE)

#### LM water depth vs tavg, disorder and tavg ####
lm_ssd_wd = lm(ssd_log ~ wd_log, data = df)
lm_tavg_wd = lm(tavg_log ~ wd_log, data = df)
lm_disorder_wd = lm(disorder_log ~ wd_log, data = df)

rsq_adj_from_lm = function(x){
  return(summary(x)$adj.r.squared)
}

p_from_lm = function(x){
  f = summary(x)$fstatistic        # value, numdf, dendf
  return(unname(pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE)))
}

slope_from_lm = function(x){
  return(summary(x)$coefficients[2,1])
}

rsq_adj_from_lm(lm_ssd_wd)


#### Logistic regression saturation vs. water depth ####
fit_wd <- glm(sat ~ wd_log, family = binomial, data = df)

summary(fit_wd)                           # slope b, SE, Wald z and p
exp(coef(fit_wd)["wd_log"])               # odds ratio per tenfold increase in water depth
exp(confint(fit_wd)["wd_log", ])          # 95% CI (profile likelihood)
anova(fit_wd, test = "Chisq")             # likelihood-ratio test vs. intercept-only model
predict(fit_wd, data.frame(wd_log = log10(c(100, 4000))), type = "response")  # P(sat) at two depths

#### S, M, L, SSD vs saturation ####
wilcox_stats = function(var, data = df){
  x   = data[[var]]
  res = wilcox.test(x ~ data$sat)
  med = tapply(x, data$sat, median) |> signif(digits = signif_dig)
  iqr = tapply(x, data$sat, IQR) |> signif(digits = signif_dig)
  p   = res$p.value |> signif(digits = signif_dig)
  A   = unname(res$statistic) / (sum(data$sat) * sum(!data$sat))
  list(var = var, W = unname(res$statistic), med = med, iqr = iqr, p = p, A = A |> signif(digits = signif_dig))
}
w_stats_s = wilcox_stats("S")
w_stats_m = wilcox_stats("M")
w_stats_l = wilcox_stats("L")
w_stats_tavg = wilcox_stats("tavg")
w_stats_disorder = wilcox_stats("disorder")
w_stats_ssd = wilcox_stats("ssd")

#### Save results ####
re.df = data.frame(description = names(results), values = unname(results))
write.csv(re.df, file = "data/res/results.csv")

cat("Done. Results are under \"figs/\"and \"data/res/results.csv\"\n")

########################################################################
#### Plotting ####
########################################################################

#### ADD in Po4 core & modeled ####

# Data from Tomasovych et al 2018
to = read.csv("data/tomasovych_et_al_2018/S0094837318000222sup002.csv", sep = ";")
to2 = dplyr::filter(to, Station == "Po4 M21" & Outlier == "")
df_te = data.frame("age" = to2$Age[!is.na(to2$Age)], "depth" = as.numeric(to2$max..depth..cm.)[!is.na(to2$Age)])
results["IQR 110 cm depth Po4"] = IQR(df_te$age[df_te$depth == 110])
results["deepest C. gibba shell younger than 25 years in Po4"] = max(df_te$depth[df_te$age <= 25])

plot_adms = function(){
  ## Axis limits for the Po4 figure
  age_lim = 150
  depth_lim = 115
  # remove samples deeper and older than plot boundaries
  df_te = dplyr::filter(df_te, depth <= depth_lim & age <= age_lim)
  p1 = ggplot(df_te, aes(x = age, y = depth)) + 
    geom_point() +
    xlim(0, age_lim) +
    scale_y_reverse(lim = c(depth_lim, 0)) +
    geom_density2d_filled(alpha = density_alpha, show.legend = FALSE) +
    ggtitle("Particle distribution core Po4") +
    xlab("Age [a]") +
    ylab("Depth [cm]") +
    scale_color_viridis_c()  +
    theme(
      plot.title = element_text(size = title_size_pts),
      axis.title.x = element_text(size = axis_title_size_pts),
      axis.title.y = element_text(size = axis_title_size_pts)
    )
  #p1
  u2 = u[ 1:depth_lim +1, 1:age_lim + 1]
  mydf = reshape2::melt(u2)
  br = c(seq(0, 0.01, by = 0.001), seq(0.01, 0.1, by = 0.02), Inf) # breaks for density
  p2 = ggplot(mydf, aes(Var2, Var1, z = value)) + 
    geom_contour_filled(breaks = br, show.legend = FALSE, alpha = contour_fill_alpha) + 
    geom_contour(breaks = br)+
    xlim(0, age_lim) +
    scale_y_reverse(lim = c(depth_lim, 0)) +
    xlab("Age [a]") +
    ylab("Depth [cm]") +
    ggtitle("Modeled distribution") +
    scale_color_viridis_c()  +
    theme(
      plot.title = element_text(size = title_size_pts),
      axis.title.x = element_text(size = axis_title_size_pts),
      axis.title.y = element_text(size = axis_title_size_pts)
    )
  #p2 
  #ggsave("figs/Po4_modeled.png", p2)
  p = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"))
  return(p)
  
}
adms_height_cm = 9
ggsave(filename = "figs/adms_po_and_modeled.pdf",
       plot = plot_adms(),
       bg = "white",
       width = double_col_width_cm,
       height = adms_height_cm,
       units = "cm")

#### dimensionless ADD ####
plot_dimless_add = function(){
  d_select = 1.0
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
  
  df = expand.grid(a, d)
  df$z = mapply(f, df$Var1, df$Var2)
  p1 = ggplot(df, aes(x = Var1, y = Var2, z = z)) +
    geom_contour_filled(breaks = sort(seq(0,1.1, by = 0.1))) +
    geom_contour(aes(z = z), color = col_line, linewidth = contour_linewidth,
                 breaks = sort(seq(0,1.1, by = 0.2))) +
    labs(x = "Age [-]",
         y = "Depth [-]",
         title = "Age-depth distribution") +
    annotate("segment",
             x = min(a), xend = max(a),
             y = d_select, yend = d_select,
             linewidth =transect_width/.pt,
             color = fixed_depth_color) +  
    annotate("segment",
             x = a_select, xend = a_select,
             y = min(d), yend = max(d),
             linewidth = transect_width/.pt,
             color = fixed_age_color) +
    coord_cartesian(clip = "on") +
    scale_y_reverse() +
    theme(legend.position = "none",
          plot.title = element_text(size = title_size_pts),
          axis.title.x = element_text(size = axis_title_size_pts),
          axis.title.y = element_text(size = axis_title_size_pts),
          axis.text = element_text(size = axis_text_size_pts)
    )
  
  
  df1 = data.frame(d = d, den_vals =  mapply(f, a_select, d))
  p2 = df1 |> ggplot(aes(x =  d, y = den_vals)) +
    geom_line(color = fixed_age_color,
              linewidth = transect_width/.pt) +
    geom_vline(xintercept = d_select,
               linetype = "dashed",
               color = fixed_depth_color) +
    scale_y_continuous(n.breaks = 3) + 
    coord_flip() +
    scale_x_reverse() +
    labs(x = "Depth [-]",
         y = "Density",
         title = expression("Position " ~ f[dis])) +
    theme(
      plot.title = element_text(size = title_size_pts),
      axis.title.x = element_text(size = axis_title_size_pts),
      axis.title.y = element_text(size = axis_title_size_pts),
      axis.text = element_text(size = axis_text_size_pts)
    )
  
  df2 = data.frame(a = a, den_vals = mapply(f, a, d_select))
  p3 = df2 |> ggplot(aes(x = a, y = den_vals)) +
    geom_line(color = fixed_depth_color,
              linewidth = transect_width/.pt) +
    geom_vline(xintercept = a_select,
               linetype = "dashed",
               color = fixed_age_color) +
    labs(x = "Age [-]",
         y = "Density",
         title = expression("Age " ~ f[tavg])) + 
    scale_y_continuous(n.breaks = 3) +
    theme(
      plot.title = element_text(size = title_size_pts),
      axis.title.x = element_text(size = axis_title_size_pts),
      axis.title.y = element_text(size = axis_title_size_pts),
      axis.text = element_text(size = axis_text_size_pts)
    )
  
  p4 = NULL
  
  p = ggpubr::ggarrange(p1, p2, p3, p4,
                        nrow = 2, 
                        ncol = 2,
                        widths = c(2,1),
                        heights = c(2,1),
                        labels = c("a", "b", "c"))
  return(p)
}

dimless_add_height_cm = 8.5

ggsave(filename = "figs/dimless_add.pdf",
       plot = plot_dimless_add(),
       height = dimless_add_height_cm,
       width = single_col_width_cm,
       units = "cm",
       bg = "white")

#### F_mix and reversal probability in dimensionless system ####
plot_fmix_rev_prob = function(){
# "Burial dominated" / "Mixing dominated" below the left/right ends of the horizontal G axis,
# on the line of the axis title; needs clip = "off" in the coord. Under coord_flip() the
# horizontal axis is the y aesthetic, so `flipped` swaps the coordinates.
g_axis_end_labels = function(flipped = FALSE, vjust = 4.6){
  lab = function(pos, hjust, label){
    if (flipped) annotate("text", x = -Inf, y = pos, hjust = hjust, vjust = vjust, label = label,
                          fontface = "italic", size = annot_size_pts / .pt)
    else         annotate("text", x = pos, y = -Inf, hjust = hjust, vjust = vjust, label = label,
                          fontface = "italic", size = annot_size_pts / .pt)
  }
  list(lab(-Inf, 0, "Burial dominated"), lab(Inf, 1, "Mixing dominated"))
}

plot_dimless_mixing_and_tavg = function(){
  pe_min_log = min(-log10(peclet_numbers))
  pe_max_log = max(-log10(peclet_numbers))
  pe_step = 1/3
  h = hist(-log10(df$Pe),
           plot = FALSE,
           breaks = seq(pe_min_log, pe_max_log, pe_step))
  
  hist_df = data.frame(x = h$mids,
                       den = h$density)
  
  line_df = data.frame(x = -log10(peclet_numbers),
                       y = tavg_dimless)
  scale_fact = max(line_df$y)/ max(hist_df$den)
  
  p = ggplot() +
    annotate("rect", xmin = -Inf, xmax = log10(G_threshold), ymin = -Inf, ymax = Inf,
             fill = col_not_saturated, alpha = region_alpha) +
    annotate("rect", xmin = log10(G_threshold), xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = col_saturated, alpha = region_alpha) +
    geom_col(data = hist_df,
             aes(x = x, y = den * scale_fact),
             width = pe_step,
             fill = col_hist_fill,
             color = col_hist_border) +
    geom_line(data = line_df,
              aes(x = x, y = y),
              linewidth = fit_linewidth,
              color = col_line) +
    geom_segment(aes(
      x = log10(G_threshold),
      xend = log10(G_threshold),
      y = 0,
      yend = Inf),
      linetype = "dashed") +
    scale_y_continuous(
      name = fmix_label,
      sec.axis = sec_axis(
        ~ ./ scale_fact,
        name = "Frequency of observations in SMLBase"
      )
    ) +
    scale_x_continuous(breaks = log10_breaks, labels = log10_labels) +
    labs(x = g_label,
         title = "Dimensionless time-averaging / disorder") +
    annotate("text",
             x = -1.5,
             y = 1,
             label = "Not saturated",
             size = annot_size_large_pts / .pt,
             fontface = "bold") +
    annotate("text",
             x = 3,
             y = 1,
             label = "Saturated",
             size = annot_size_large_pts / .pt,
             fontface = "bold")  +
    g_axis_end_labels() +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(size = title_size_pts),
      axis.title.x = element_text(size = axis_title_size_pts),
      axis.title.y = element_text(size = axis_title_size_pts)
    )
  return(p)
}

plot_rev_prob = function(){
  depths = seq(0, 3, by = 0.02)
  names = c("Pe", "P", "L_dimless")
  df_p = data.frame(matrix(nrow = 0, ncol = length(names)))
  
  for (i in seq_along(peclet_numbers)){
    fval = tavg_list[[i]]$den
    pe = peclet_numbers[i]
    
    aa = cumul_int_trap(t_dimless, fval, normalize = TRUE)
    den = c(0, diff(aa)/diff(t_dimless))
    CDF = approxfun(x = t_dimless, y = aa, yleft = 0, yright = 1)
    
    PDF = approxfun(x = t_dimless, y = den, yleft = 0, yright = 0)
    
    comb_f = function(y){
      comb_fun = function(x){PDF(x)* CDF(x-y)}
      return(integrate(comb_fun, lower = - Inf, upper = Inf)[[1]])
    }
    
    depths = seq(0, 3, by = 0.02)
    df_p = rbind(df_p, data.frame("Pe" = rep(pe, length(depths)),
                                  "P" = sapply(depths, comb_f),
                                  "L_dimless" = depths))
  }
  
  df_p$logG = -log10(df_p$Pe)
  
  levels = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01)
  p = df_p |>
    ggplot(aes(x = L_dimless, y = logG, z = P)) +
    geom_contour_filled(breaks = c(0, 0.1), fill = col_low_risk, alpha = risk_fill_alpha) +
    geom_contour_filled(breaks = c(0.1, 0.5), fill = col_high_risk, alpha = risk_fill_alpha) +
    geom_contour(breaks = setdiff(levels, 0.1), colour = col_contour) +
    geom_contour(breaks = 0.1, colour = col_line, linewidth = highlight_contour_linewidth) +   # highlighted contour
    geom_hline(yintercept = log10(G_threshold), linetype = "dashed") +   # saturated vs. not saturated
    geom_text_contour(breaks = levels, skip = 0, stroke = 0.2, size = annot_size_large_pts / .pt) +
    annotate("text", x = 1.9, y = 2.5, label = "Low reversal risk",
             colour = col_line, fontface = "bold", size = annot_size_large_pts / .pt) +
    annotate("text", x = 1.3, y = 2.5, label = "High reversal risk",
             colour = col_line, fontface = "bold", size = annot_size_large_pts / .pt) +
    annotate("text", x = 2.7, y = -1.5, label = "Not saturated", angle = 45,
             colour = col_line, fontface = "bold", size = annot_size_large_pts / .pt) +
    annotate("text", x = 2.7, y = 3, label = "Saturated", angle = 45,
             colour = col_line, fontface = "bold", size = annot_size_large_pts / .pt) +
    labs(title = "Probability of age reversals",
         x = "Sampling distance [prop. SML depth]",
         y = g_label) +
    scale_y_continuous(breaks = log10_breaks, labels = log10_labels) +
    g_axis_end_labels(flipped = TRUE) +
    coord_flip(clip = "off") +
    theme(
      plot.title = element_text(size = title_size_pts),
      axis.title.x = element_text(size = axis_title_size_pts),
      axis.title.y = element_text(size = axis_title_size_pts)
    )
  return(p)
}

p = plot_rev_prob()

p_comb = ggpubr::ggarrange(plot_dimless_mixing_and_tavg(), plot_rev_prob(),
                           nrow = 1, ncol = 2,
                           labels = c("a", "b"))
return(p_comb)
}

dimless_mixing_height_cm = 10
ggsave(filename = "figs/f_mix_and_rev_prob.pdf",
       plot = plot_fmix_rev_prob(),
       bg = "white",
       width = double_col_width_cm,
       height = dimless_mixing_height_cm,
       units = "cm")
#### Saturation state plot (log regression and boxplots) ####

saturation_state_plot = function(){
  # Mann-Whitney U test annotation (same definitions as wilcox_stats()):
  # W is the rank-sum statistic of the "Not saturated" group, A = W / (n_not_sat * n_sat)
  # is the probability of superiority P(Not saturated > Saturated).
  # Only finite values (= the points shown in the panel) enter the test; the test is
  # rank-based, so testing on the log scale gives the same result as on the linear scale.
  mwu_annot = function(yvar){
    d   = df[is.finite(df[[yvar]]), ]
    res = suppressWarnings(wilcox.test(d[[yvar]] ~ d$sat))
    W   = unname(res$statistic)
    A   = W / (sum(d$sat == 1) * sum(d$sat == 0))
    p_lab = if (res$p.value < 0.001) "p < 0.001" else paste("p =", signif(res$p.value, signif_digits))
    sprintf("W = %s, %s\nA = %s", signif(W, signif_digits), p_lab, signif(A, signif_digits))
  }

  sat_boxplot = function(yvar, ylab, title){
    ggplot(df, aes(x = saturated, y = .data[[yvar]], fill = saturated)) +
      geom_boxplot(alpha = box_alpha, outlier.shape = NA) +
      geom_jitter(aes(colour = saturated, shape = saturated),
                  width = 0.15, height = 0, size = sat_jitter_size, alpha = sat_jitter_alpha) +
      annotate("text", x = 1.5, y = Inf, label = mwu_annot(yvar),
               vjust = 1.2, lineheight = 0.9, size = annot_size_pts / .pt) +
      log10_y_scale(expand = expansion(mult = c(0.05, 0.4))) +   # headroom for the annotation
      scale_fill_manual(values = sat_cols, guide = "none") +
      scale_colour_manual(values = sat_cols, guide = "none") +
      scale_shape_manual(values = sat_shapes, guide = "none") +
      labs(x = NULL, y = ylab, title = title) +
      theme_bw() +
      theme(plot.title   = element_text(size = title_size_pts),
            axis.title.y = element_text(size = y_axis_text_size_pts),
            axis.text.x  = element_text(size = x_axis_text_size_pts))
  }
  
  p1 = sat_boxplot("S_log",   sedr_label,      "Sed. acc. rate")
  p2 = sat_boxplot("M_log",   mix_label,       "Biodiffusion")
  p3 = sat_boxplot("L_log",   mix_depth_label, "Mixing depth")
  p4 = sat_boxplot("tavg_log",     tavg_label,          "Time-averaging")
  p5 = sat_boxplot("disorder_log", disorder_label,      "Strat. disorder")
  p6 = sat_boxplot("ssd_log",      ssd_log_label_short, "Safe sampling dist.")

  # prediction grid for the logistic regression: fitted P(saturated) with 95% CI,
  # computed on the link scale and back-transformed so bounds stay within [0, 1]
  nd = data.frame(wd = 10^seq(min(df$wd_log), max(df$wd_log), length.out = 200))
  nd$wd_log = log10(nd$wd)
  pr = predict(fit_wd, newdata = nd, type = "link", se.fit = TRUE)
  nd$p   = plogis(pr$fit)
  nd$lwr = plogis(pr$fit - qnorm(0.975) * pr$se.fit)
  nd$upr = plogis(pr$fit + qnorm(0.975) * pr$se.fit)
  
  # annotations: odds ratio per tenfold increase in water depth, and likelihood-ratio test
  or_wd = exp(coef(fit_wd)["wd_log"])
  or_ci = exp(confint(fit_wd)["wd_log", ])
  lrt = anova(fit_wd, test = "Chisq")
  lrt_df = lrt["wd_log", "Df"]
  lrt_chi = signif(lrt["wd_log", "Deviance"], signif_digits)
  lrt_p = lrt["wd_log", "Pr(>Chi)"]
  p_lab = if (lrt_p < 0.001) "< 0.001" else paste("=", signif(lrt_p, signif_digits))
  annot1 = sprintf("Odds ratio per 10x depth: %s (95%% CI: %s to %s)",
                   signif(or_wd, signif_digits),
                   signif(or_ci[1], signif_digits),
                   signif(or_ci[2], signif_digits))
  annot2 = as.expression(bquote(LRT ~ chi^2 * "(" * .(lrt_df) * ")" == .(lrt_chi) * "," ~ italic(p) ~ .(p_lab)))
  
  # colour the logistic-regression points the same way
  p0 = ggplot() +
    geom_jitter(data = df, aes(wd, sat, colour = saturated, shape = saturated),
                height = 0.03, width = 0, alpha = logit_point_alpha, size = logit_point_size) +
    geom_ribbon(data = nd, aes(wd, ymin = lwr, ymax = upr), alpha = logit_ribbon_alpha) +
    geom_line(data = nd, aes(wd, p), linewidth = fit_linewidth) +
    scale_colour_manual(values = sat_cols, name = NULL) +
    scale_shape_manual(values = sat_shapes, name = NULL) +
    scale_x_log10() +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Water depth [m]", y = "P(Saturated)",
         title = "Probability of mixing saturation vs. water depth") +
    theme_bw() +
    theme(plot.title = element_text(size = title_size_pts)) +
    annotate("text", x = min(df$wd), y = 0.25, label = annot1, hjust = 0, size = annot_size_pts / .pt) +
    annotate("text", x = min(df$wd), y = 0.15, label = annot2, hjust = 0, size = annot_size_pts / .pt)
  
  fig = p0 / (p1 | p2 | p3) / (p4 | p5 | p6) +
    plot_layout(guides = "collect",
                heights = c(1.2, 1, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom")
  return(fig)
}

saturation_state_height_cm = 20
ggsave(filename = "figs/saturation_state.pdf",
       plot = saturation_state_plot(),
       width = double_col_width_cm,
       height = saturation_state_height_cm,
       units = "cm",
       device = cairo_pdf)

#### GLM S, M, L vs tavg, disorder, ssd ####
tavg_y_axis_lims = range(log10(c(df$tavg* 1.3, df$tavg * 0.8)))
disorder_y_axis_lims = range(log10(c(1.1 * df$disorder, 0.9 * df$disorder)))


# ssd_glm, ssd_dat and rsq_part_ssd are defined in the GLM section above
ssd_y_axis_lims = range(c(ssd_dat$ssd + log10(1.2), ssd_dat$ssd + log10(0.85)))

## Generic partial-residual panel used for all rows of the GLM summary figure
plot_partial_panel = function(model, xvar, xlab, ylab, ylims, rsq_obj, group){
  v = visreg::visreg(model, xvar = xvar, plot = FALSE)
  res = v$res
  res$group = group
  
  rsq_val = rsq_obj$partial.rsq[rsq_obj$variable == xvar]
  vv = paste(": ", round(rsq_val, 3))
  p_val = p_from_glm(model)[[xvar]]
  p_lab = if (p_val < 0.001) "< 0.001" else paste("=", signif(p_val, signif_digits))
  annot = as.expression(bquote(partial ~ R^2 * .(vv) * "," ~ italic(p) ~ .(p_lab)))
  slope_annot = paste("Slope:", round(coef(model)[xvar], 3))
  
  x_mid = mean(range(res[[xvar]]))
  y_top = max(ylims) - 0.05 * diff(ylims)
  dy = 0.1 * diff(ylims)
  
  ggplot() +
    geom_ribbon(data = v$fit, aes(x = .data[[xvar]], ymin = visreg_lwr, ymax = visreg_upr),
                fill = col_ribbon) +
    geom_line(data = v$fit, aes(x = .data[[xvar]], y = visreg_fit), linewidth = fit_linewidth) +
    geom_point(data = res, aes(x = .data[[xvar]], y = visreg_res, colour = group, shape = group)) +
    scale_colour_manual(name = NULL,
                        values = sat_cols) +
    scale_shape_manual(name = NULL,
                       values = sat_shapes) +
    log10_x_scale() +
    log10_y_scale() +
    xlab(xlab) +
    ylab(ylab) +
    coord_cartesian(ylim = ylims) +   # zooms without dropping points (unlike ylim())
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = x_mid, y = y_top,      label = annot,       size = annot_size_pts / .pt) +
    annotate("text", x = x_mid, y = y_top - dy, label = slope_annot, size = annot_size_pts / .pt)
}

plot_glm_summary_fig = function(){

  # one row of panels: S, D_b and L as predictors of one response
  # the row title sits on the middle panel; blank titles keep the panels aligned
  panel_row = function(row_title, model, ylab, ylims, rsq_obj, group){
    row_panels = list(plot_partial_panel(model, "s1", sedr_label,      ylab, ylims, rsq_obj, group),
                      plot_partial_panel(model, "m1", mix_label,       ylab, ylims, rsq_obj, group),
                      plot_partial_panel(model, "l1", mix_depth_label, ylab, ylims, rsq_obj, group))
    Map(function(p, title){
      p + ggtitle(title) + theme(plot.title = element_text(size = title_size_pts, hjust = 0.5))
    }, row_panels, c(" ", row_title, " "))
  }

  panels = c(panel_row("Time-averaging",         tavg_glm,     tavg_label,     tavg_y_axis_lims,     rsq_part_tavg,     df$saturated),
             panel_row("Disorder",               disorder_glm, disorder_label, disorder_y_axis_lims, rsq_part_disorder, df$saturated),
             panel_row("Safe sampling distance", ssd_glm,      ssd_log_label,  ssd_y_axis_lims,      rsq_part_ssd,      ssd_dat$group))

  ggpubr::ggarrange(plotlist = panels,
                    nrow = 3, ncol = 3,
                    labels = letters[1:9],
                    common.legend = TRUE,
                    legend = "bottom")
}
glm_summary_height_cm = 21
ggsave("figs/glm_summary_fig.pdf",
       plot = plot_glm_summary_fig(),
       width = double_col_width_cm,
       height = glm_summary_height_cm,
       units = "cm",
       bg = "white")

#### wd vs tavg, disorder, and ssd (plot ssd histogram) ####


aa = summary(lm_ssd_wd)
aa$coefficients |> signif(digits = signif_dig)


# x axis: data stay on log10 scale, ticks are labelled in metres
wd_x_scale = list(
  scale_x_continuous(breaks = log10(c(10, 100, 1000, 10000)),
                     labels = c("10", "100", "1000", "10000"),
                     minor_breaks = log10(as.vector(outer(1:9, 10^(1:4))))),
  annotation_logticks(sides = "b", outside = FALSE,
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"))
)

signif_dig = 3
## Regression panel against water depth, same styling as the GLM summary figure
plot_wd_panel = function(model, title, ylab, ylims, group){
  v = visreg::visreg(model, xvar = "wd_log", plot = FALSE)
  res = v$res
  res$group = group

  vv = paste(": ", round(rsq_adj_from_lm(model), 3))
  p_val = p_from_lm(model)
  p_lab = if (p_val < 0.001) "< 0.001" else paste("=", signif(p_val, signif_digits))
  annot = as.expression(bquote(adj. ~ R^2 * .(vv) * "," ~ italic(p) ~ .(p_lab)))
  slope_annot = paste("Slope:", round(slope_from_lm(model), 3))

  x_mid = mean(range(res$wd_log))
  y_top = max(ylims) - 0.05 * diff(ylims)
  dy = 0.1 * diff(ylims)

  ggplot() +
    geom_ribbon(data = v$fit, aes(x = wd_log, ymin = visreg_lwr, ymax = visreg_upr),
                fill = col_ribbon) +
    geom_line(data = v$fit, aes(x = wd_log, y = visreg_fit), linewidth = fit_linewidth) +
    geom_point(data = res, aes(x = wd_log, y = visreg_res, colour = group, shape = group)) +
    scale_colour_manual(name = NULL,
                        values = sat_cols) +
    scale_shape_manual(name = NULL,
                       values = sat_shapes) +
    wd_x_scale +
    log10_y_scale() +
    ggtitle(title) +
    xlab(wd_lin_label) +
    ylab(ylab) +
    coord_cartesian(ylim = ylims) +
    theme(plot.title = element_text(size = title_size_pts, hjust = 0.5),
          axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = x_mid, y = y_top,      label = annot,       size = annot_size_pts / .pt) +
    annotate("text", x = x_mid, y = y_top - dy, label = slope_annot, size = annot_size_pts / .pt)
}

plot_wd_vs_params = function(){
  # same limits as the GLM summary figure, plus room above the data for the annotations
  with_headroom = function(ylims) c(ylims[1], ylims[2] + 0.25 * diff(ylims))
  ssd_lims = with_headroom(ssd_y_axis_lims)

  lm1 = plot_wd_panel(lm_tavg_wd,     "Time-averaging",         tavg_label,     with_headroom(tavg_y_axis_lims),     df$saturated)
  lm2 = plot_wd_panel(lm_disorder_wd, "Disorder",               disorder_label, with_headroom(disorder_y_axis_lims), df$saturated)
  lm3 = plot_wd_panel(lm_ssd_wd,      "Safe sampling distance", ssd_log_label,  ssd_lims,                            df$saturated)

  hist1 = df |>
    ggplot(aes(y = ssd_log)) +   # y instead of x + coord_flip()
    geom_histogram(binwidth = 0.1, fill = col_hist_fill, color = col_hist_border) +
    geom_hline(yintercept = median(df$ssd_log),   # = log10 of the median SSD
               linetype = "dashed") +
    ggtitle("Distribution of safe sampling distance") +
    labs(y = ssd_log_label,
         x = "Count in SMLBase") +
    log10_y_scale() +
    coord_cartesian(ylim = ssd_lims) +
    theme(plot.title = element_text(size = title_size_pts, hjust = 0.5),
          axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))

  ggpubr::ggarrange(lm1, lm2, lm3, hist1, ncol = 2, nrow = 2,
                    labels = letters[1:4],
                    common.legend = TRUE,
                    legend = "bottom")
}

wd_params_height_cm = 15
ggsave(filename = "figs/water_depth_vs_params.pdf",
       plot = plot_wd_vs_params(),
       bg = "white",
       width = double_col_width_cm,
       height = wd_params_height_cm,
       unit = "cm")

#### Figure: Overview of SML parameters ####
plot_sml_overview_fig = function(){
  hist_theme = theme(plot.title = element_text(size = title_size_pts, hjust = 0.5),
                     axis.title.x = element_text(size = x_axis_text_size_pts),
                     axis.title.y = element_text(size = y_axis_text_size_pts))

  # histogram of one SML parameter, on the same log10 scale as in the other figures
  param_hist = function(xvar, xlab, title){
    ggplot(df, aes(x = .data[[xvar]])) +
      geom_histogram(bins = 30, fill = col_hist_fill, color = col_hist_border) +
      ggtitle(title) +
      log10_x_scale() +
      labs(x = xlab, y = "Count in SMLBase") +
      hist_theme
  }

  p_s = param_hist("S_log", sedr_label,      "Sediment accumulation rate")
  p_m = param_hist("M_log", mix_label,       "Biodiffusion")
  p_l = param_hist("L_log", mix_depth_label, "Mixing depth")

  p_fmix = ggplot(df, aes(x = F_mix)) +
    annotate("rect", xmin = -Inf, xmax = F_mix_threshold, ymin = -Inf, ymax = Inf,
             fill = col_not_saturated, alpha = region_alpha) +
    annotate("rect", xmin = F_mix_threshold, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = col_saturated, alpha = region_alpha) +
    geom_histogram(bins = 30, fill = col_hist_fill, color = col_hist_border) +
    geom_vline(xintercept = F_mix_threshold, linetype = "dashed") +
    annotate("text", x = 0.5, y = Inf, vjust = 1.5,
             label = "Not saturated",
             fontface = "bold", size = annot_size_large_pts / .pt) +
    annotate("text", x = 1.03, y = Inf, hjust = 1.1, angle = 90,
             label = "Saturated",
             fontface = "bold", size = annot_size_large_pts / .pt) +
    ggtitle("Dimensionless time-averaging/disorder") +
    labs(x = fmix_label, y = "Count in SMLBase") +
    hist_theme

  ggpubr::ggarrange(p_s, p_m, p_l, p_fmix,
                    labels = letters[1:4],
                    ncol = 2, nrow = 2)
}
sml_overview_height_cm = 18
ggsave(filename = "figs/sml_overview.pdf",
       plot = plot_sml_overview_fig(),
       bg = "white",
       width = double_col_width_cm,
       height = sml_overview_height_cm,
       units = "cm")

#### Coefficient plot ####
aa = summary(ssd_glm_std)
aa$coefficients |> signif(digits = signif_dig)
make_coeff_plots = function(){
  # predictor labels in model order (S, D_b, L)
  predictor_labels = c(expression(log[10](S)),
                       expression(log[10](D[b])),
                       expression(log[10](L)))

  # coefficient estimates with 95% CI for one model
  coef_panel = function(model, xlab){
    terms = names(coef(model))[-1]   # model terms, matched to the labels by name
    dwplot(model,
           dot_args = list(size = coef_dot_size, colour = col_line),
           whisker_args = list(size = coef_whisker_size, colour = col_line)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_y_discrete(labels = function(b) predictor_labels[match(b, terms)]) +
      labs(x = xlab, y = "Independent variables") +
      theme(legend.position = "none",
            axis.text.y = element_text(angle = 45),
            axis.title.x = element_text(size = x_axis_text_size_pts),
            axis.title.y = element_text(size = y_axis_text_size_pts))
  }

  # one row per response: raw coefficients (left) and normalized coefficients (right),
  # with the row title centred above both panels as in the GLM summary figure
  coef_row = function(row_title, model, model_std, panel_labels){
    ggpubr::annotate_figure(
      ggpubr::ggarrange(coef_panel(model,     "Regression coefficient"),
                        coef_panel(model_std, "Normalized regression coefficient"),
                        ncol = 2, labels = panel_labels),
      top = ggpubr::text_grob(row_title, size = title_size_pts))
  }

  ggpubr::ggarrange(coef_row("Time-averaging",         tavg_glm,     tavg_glm_std,     c("a", "b")),
                    coef_row("Disorder",               disorder_glm, disorder_glm_std, c("c", "d")),
                    coef_row("Safe sampling distance", ssd_glm,      ssd_glm_std,      c("e", "f")),
                    nrow = 3)
}

coeff_plots_height_cm = 21
ggsave(filename = "figs/coefficient_plots.pdf",
       plot = make_coeff_plots(),
       width = double_col_width_cm,
       height = coeff_plots_height_cm,
       units = "cm",
       bg = "white")

#### Sediment accumulation rate vs. water depth ####
lm_wd_S = lm(S_log ~ wd_log, data = df)
lm_wd_M = lm(M_log ~ wd_log, data = df)
lm_wd_L = lm(L_log ~ wd_log, data = df)

summary_lm_wd_S = summary(lm_wd_S)
summary_lm_wd_S$coefficients |> signif(digits = signif_dig)
summary_lm_wd_S$adj.r.squared |> signif(digits = signif_dig)
f = summary_lm_wd_S$fstatistic          # named vector: value, numdf, dendf
p_r2 = pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE) |> unname()

# x axis: data stay on log10 scale, ticks are labelled in metres
wd_x_scale = list(
  scale_x_continuous(breaks = log10(c(10, 100, 1000, 10000)),
                     labels = c("10", "100", "1000", "10000"),
                     minor_breaks = log10(as.vector(outer(1:9, 10^(1:4))))),
  annotation_logticks(sides = "b", outside = FALSE,
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"))
)

plot_wd_vs_s = function(){
  # padded like the GLM summary limits, plus room above the data for the annotations
  s_lims = range(log10(c(df$S * 1.3, df$S * 0.8)))
  s_lims = c(s_lims[1], s_lims[2] + 0.25 * diff(s_lims))

  plot_wd_panel(lm_wd_S, "Sediment accumulation rate", sedr_label, s_lims, df$saturated) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal")
}

wd_vs_sed_height_cm = 9
ggsave("figs/wd_vs_sed_acc.pdf",
       plot = plot_wd_vs_s(),
       width = single_col_width_cm,
       height = wd_vs_sed_height_cm,
       units = "cm",
       bg = "white")

#### S, M, L vs. water depth ####
plot_wd_vs_sml = function(){
  # padded like the GLM summary limits, plus room above the data for the annotations
  lims = function(x){
    l = range(log10(c(x * 1.3, x * 0.8)))
    c(l[1], l[2] + 0.35 * diff(l))
  }

  p_s = plot_wd_panel(lm_wd_S, "Sediment accumulation rate", sedr_label,      lims(df$S), df$saturated)
  p_m = plot_wd_panel(lm_wd_M, "Biodiffusion",               mix_label,       lims(df$M), df$saturated)
  p_l = plot_wd_panel(lm_wd_L, "Mixing depth",               mix_depth_label, lims(df$L), df$saturated)

  ggpubr::ggarrange(p_s, p_m, p_l, ncol = 3, nrow = 1,
                    labels = letters[1:3],
                    common.legend = TRUE,
                    legend = "bottom")
}

wd_vs_sml_height_cm = 8
ggsave("figs/wd_vs_sml.pdf",
       plot = plot_wd_vs_sml(),
       width = double_col_width_cm,
       height = wd_vs_sml_height_cm,
       units = "cm",
       bg = "white")

#### Pairwise scatterplots of S, M, L ####
plot_sml_pairs = function(){
  pair_panel = function(xvar, yvar, xlab, ylab, title){
    d = df[complete.cases(df[, c(xvar, yvar)]), ]
    ct = cor.test(d[[xvar]], d[[yvar]], method = "pearson")
    p_lab = if (ct$p.value < 0.001) "< 0.001" else paste("=", signif(ct$p.value, signif_digits))
    annot = as.expression(bquote(italic(r) == .(round(unname(ct$estimate), 3)) * "," ~ italic(p) ~ .(p_lab)))

    # room above the data for the annotation
    ylims = range(d[[yvar]])
    ylims = c(ylims[1], ylims[2] + 0.2 * diff(ylims))

    ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]], colour = saturated, shape = saturated)) +
      geom_point() +
      scale_colour_manual(name = NULL, values = sat_cols) +
      scale_shape_manual(name = NULL, values = sat_shapes) +
      log10_x_scale() +
      log10_y_scale() +
      ggtitle(title) +
      xlab(xlab) +
      ylab(ylab) +
      coord_cartesian(ylim = ylims) +
      theme(plot.title = element_text(size = title_size_pts, hjust = 0.5),
            axis.title.x = element_text(size = x_axis_text_size_pts),
            axis.title.y = element_text(size = y_axis_text_size_pts)) +
      annotate("text", x = mean(range(d[[xvar]])), y = max(ylims) - 0.05 * diff(ylims),
               label = annot, size = annot_size_pts / .pt)
  }

  p_sm = pair_panel("S_log", "M_log", sedr_label, mix_label,       "Sed. acc. rate vs. biodiffusion")
  p_sl = pair_panel("S_log", "L_log", sedr_label, mix_depth_label, "Sed. acc. rate vs. mixing depth")
  p_ml = pair_panel("M_log", "L_log", mix_label,  mix_depth_label, "Biodiffusion vs. mixing depth")

  ggpubr::ggarrange(p_sm, p_sl, p_ml, ncol = 3, nrow = 1,
                    labels = letters[1:3],
                    common.legend = TRUE,
                    legend = "bottom")
}

sml_pairs_height_cm = 8
ggsave("figs/sml_pairs.pdf",
       plot = plot_sml_pairs(),
       width = double_col_width_cm,
       height = sml_pairs_height_cm,
       units = "cm",
       bg = "white")
