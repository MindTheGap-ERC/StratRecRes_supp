#### Load required packages ####
library(ggplot2)
library(dotwhisker)
library(gapminder)

#### Plotting constants
ax_label_ssd = 

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
df = df |> mutate(ssd_log = log10(ssd))
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
ssd_glm = glm(ssd_log ~ S_log + M_log + L_log, data = df)

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

## Save results for Rsq and coefficients
results[c("partial R2 disorder D_b", "partial R2 disorder S", "partial R2 disorder L")] = round(rsq_part_disorder$partial.rsq, round_res)
results[c("coeff disorder D_b", "coeff disorder S", "coeff disorder L")] = round(disorder_glm$coefficients[-1], round_res)
results[c("partial R2 tavg D_b", "partial R2 tavg S", "partial R2 tavg L")] = round(rsq_part_tavg$partial.rsq, round_res)
results[c("coeff tavg D_b", "coeff tavg S", "coeff tavg L")] = round(tavg_glm$coefficients[-1], round_res)

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
vif_ssd_std = car::vif(disorder_glm_std)
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

summary(fit)                              # slope b, SE, Wald z and p
exp(coef(fit)["wd_log"])                  # odds ratio per tenfold increase in water depth
exp(confint(fit)["wd_log", ])             # 95% CI (profile likelihood)
anova(fit, test = "Chisq")                # likelihood-ratio test vs. intercept-only model
predict(fit, data.frame(wd_log = log10(c(100, 4000))), type = "response")  # P(sat) at two depths

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
w_stats_ssd = wilcox_stats("ssd")



sat_cols = c("Not saturated" = "goldenrod", "Saturated" = "mediumpurple")

title_size_pts = 10

#### Plot: saturation state vs. S, M, L, and SSD ####
saturation_state_plot = function(){
  sat_boxplot = function(yvar, ylab, title){
  ggplot(df, aes(x = saturated, y = .data[[yvar]], fill = saturated)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(colour = saturated, shape = saturated),
                width = 0.15, height = 0, size = 0.8, alpha = 0.9) +
    scale_fill_manual(values = sat_cols, guide = "none") +
    scale_colour_manual(values = sat_cols, guide = "none") +
    scale_shape_manual(values = c("Not saturated" = 17, "Saturated" = 16), guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_bw() +
    theme(plot.title   = element_text(size = title_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts),
          axis.text.x  = element_text(size = x_axis_text_size_pts))
}

ssd_log_label = expression( log[10]* group("(", Safe~sampling~dis.,")") ~ group("[",cm,"]") )
df$ssd_log = log10(df$safe_sampling_distance)

p1 = sat_boxplot("S_log",   sedr_label,      "Sed. acc. rate")
p2 = sat_boxplot("M_log",   mix_label,       "Biodiffusion")
p3 = sat_boxplot("L_log",   mix_depth_label, "Mixing depth")
p4 = sat_boxplot("ssd_log", ssd_log_label,   "Safe sampling dist.")

# colour the logistic-regression points the same way
p0 = ggplot() +
  geom_jitter(data = df, aes(wd, sat, colour = saturated, shape = saturated),
              height = 0.03, width = 0, alpha = 0.8, size = 1) +
  geom_ribbon(data = nd, aes(wd, ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(data = nd, aes(wd, p), linewidth = 1) +
  scale_colour_manual(values = sat_cols, name = NULL) +
  scale_shape_manual(values = c("Not saturated" = 17, "Saturated" = 16), name = NULL) +
  scale_x_log10() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Water depth [m]", y = "P(Saturated)",
       title = "Probability of mixing saturation vs. water depth") +
  theme_bw() +
  theme(plot.title = element_text(size = title_size_pts)) +
  annotate("text", x = min(df$wd), y = 0.25, label = annot1, hjust = 0, size = annot_size_pts / .pt) +
  annotate("text", x = min(df$wd), y = 0.15, label = annot2, hjust = 0, size = annot_size_pts / .pt)

library(patchwork)
fig = p0 / (p1 | p2 | p3 | p4) +
  plot_layout(guides = "collect",
              heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
  return(fig)
}

ggsave(filename = "figs/saturation_state.pdf",
       plot = saturation_state_plot(),
       width = 19,
       height = 12,
       units = "cm",
       device = cairo_pdf)
#### ####



#### Safe sampling distance
glm_ssd = glm(log10(safe_sampling_distance) ~ S_log + L_log + M_log, data = df)
summary(glm_ssd)
visreg(glm_ssd)

lm_ssd = lm(safe_sampling_distance ~ wd_log, data = df)
visreg(lm_ssd)


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
    scale_color_viridis_c()  +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
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
    ggtitle("Modelled Distribution") +
    scale_color_viridis_c()  +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  #p2 
  #ggsave("figs/Po4_modeled.png", p2)
  p = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"))
  return(p)

}
ggsave(filename = "figs/adms_po_and_modeled.pdf",
       plot = plot_adms(),
       bg = "white",
       width = 18,
       height = 8,
       units = "cm")

#### Figure 2: histogram of mixing intensity and F_mix ####
fmix_label =  expression( F[mix] ~ " [-]" )  #"log10(Mixing depth) [cm]"
g_label =  expression( log[10]* group("(", G,")")  ~  group("[","-","]") )  #"log10(Mixing depth) [cm]"

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
             fill = "goldenrod", alpha = 0.25) +
    annotate("rect", xmin = log10(G_threshold), xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "mediumpurple", alpha = 0.25) +
    geom_col(data = hist_df,
             aes(x = x, y = den * scale_fact),
             width = pe_step,
             fill = "grey80",
             color = "grey40") +
    geom_line(data = line_df,
              aes(x = x, y = y),
              linewidth = 1,
              color = "black") +
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
    labs(x = g_label,
         title = "Dimensionless time-averaging / disorder") +
    annotate("text",
             x = -1.5,
             y = 1,
             label = "Not saturated",
             size = 8/.pt,
             fontface = "bold") +
    annotate("text",
             x = 3,
             y = 1,
             label = "Saturated",
             size = 8/.pt,
             fontface = "bold")  +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
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
  
  library(ggplot2)
  library(metR)
  levels = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01)
  p = df_p |>
    ggplot(aes(x = L_dimless, y = logG, z = P)) +
    geom_contour_filled(breaks = c(0, 0.1), fill = "aquamarine", alpha = 0.3) +
    geom_contour_filled(breaks = c(0.1, 0.5), fill = "tomato", alpha = 0.3) +
    geom_contour(breaks = setdiff(levels, 0.1), colour = "grey30") +
    geom_contour(breaks = 0.1, colour = "black", linewidth = 1.2) +   # highlighted contour
    geom_text_contour(breaks = levels, skip = 0, stroke = 0.2, size = 8/.pt) +
    annotate("text", x = 1.9, y = 2.5, label = "Low reversal risk",
             colour = "black", fontface = "bold", size = 8/.pt) +
    annotate("text", x = 1.3, y = 2.5, label = "High reversal risk",
             colour = "black", fontface = "bold", size = 8/.pt) +
    labs(title = "Probability of age reversals",
         x = "Safe sampling distance [prop. SML depth]",
         y = g_label) +
    coord_flip() +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  return(p)
}

p = plot_rev_prob()

p_comb = ggpubr::ggarrange(plot_dimless_mixing_and_tavg(), plot_rev_prob(),
                           nrow = 1, ncol = 2,
                           labels = c("a", "b"))

ggsave(filename = "figs/dimless_mixing_and_tavg.pdf",
       plot = p_comb,
       bg = "white",
       width = 18,
       height = 10,
       units = "cm")


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

ggplot() +
  geom_col(data = hist_df,
           aes(x = x, y = den * scale_fact),
           width = pe_step,
           fill = "grey80",
           color = "grey40") +
  geom_line(data = line_df,
            aes(x = x, y = y),
            linewidth = 1,
            color = "black") +
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
      name = "Frequency"
    )
  ) +
  labs(x = g_label) +
  annotate("text",
           x = -1.5,
           y = 1,
           label = "Not saturated") +
  annotate("text",
           x = 3,
           y = 1,
           label = "Saturated")

#### Figure 3: GLM plot ####
sedr_label = expression( log[10] * "(" * "Sed. acc." * ")" ~ group("[", cm/a, "]") )
mix_label =  expression( log[10]* group("(", Biodiffusion*"",")")  ~  group("[",cm^2/a,"]") ) #"log10(Mixing int.) [cm^2/a]"
mix_depth_label =  expression( log[10]* group("(", Mixing~depth,")")  ~  group("[",cm,"]") )  #"log10(Mixing depth) [cm]"
tavg_label =  expression( log[10]* group("(", Time-averaging,")")  ~  group("[",a,"]") ) # "log10(Time averaging) [a]"
disorder_label = expression( log[10]* group("(", Strat.~disorder,")")  ~  group("[",cm,"]") )  #"log10(Stratigraphic disorder) [cm]"
wd_label =  expression( log[10]* group("(", Water~depth,")")  ~  group("[",m,"]") )  #"log10(Mixing depth) [cm]"

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
  vv = paste(": ", round(rsq$partial.rsq[1], 3))
  annot = as.expression(bquote(partial ~ R^2* .(vv)))
  slope_annot = paste("Slope: ", round(coefficients(tavg_glm)["s1"], 4))
  v = visreg::visreg(tavg_glm, xvar = "s1", gg = TRUE, plot = FALSE)
  
  
  res <- v$res
  res$group <- df$saturated   # see the NA note below
  
  tavg_glm_plot_s = ggplot() +
    geom_ribbon(data = v$fit, aes(x = s1, ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = s1, y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = s1, y = visreg_res, colour = group, shape = group)) +
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")) +
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16)) + # filled circle
    xlab(sedr_label) +
    ylab(tavg_label) +
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$tavg)) - 0.5, label = slope_annot, size = annot_size_pts / .pt)
  tavg_glm_plot_s
  
  # biodiffusion
  vv = paste(": ", round(rsq$partial.rsq[2], 3))
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
  slope_annot = paste("Slope: ", round(coefficients(tavg_glm)["m1"], 3))
  v = visreg::visreg(tavg_glm, xvar = "m1", plot = FALSE)
  
  res <- v$res
  res$group <- df$saturated
  
  tavg_glm_plot_m <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = m1, ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = m1, y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = m1, y = visreg_res, colour = group, shape = group)) +
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")) +
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16)) + # filled circle
    xlab(mix_label) +
    ylab(tavg_label) +
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$tavg)) - 0.5, label = slope_annot, size = annot_size_pts / .pt)
  tavg_glm_plot_m
  # mixing depth
  vv = paste(": ", round(rsq$partial.rsq[3], 3))
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
  slope_annot = paste("Slope: ", round(coefficients(tavg_glm)["l1"], 3))
  v = visreg::visreg(tavg_glm, xvar = "l1", plot = FALSE)
  
  res <- v$res
  res$group <- df$saturated
  
  tavg_glm_plot_l <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = l1, ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = l1, y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = l1, y = visreg_res, colour = group, shape = group)) +
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")) +
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16)) + # filled circle
    xlab(mix_depth_label) +
    ylab(tavg_label) +
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$tavg)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$tavg)) - 0.5, label = slope_annot, size = annot_size_pts / .pt)
  tavg_glm_plot_l
  
  ## shared point styling
  sat_scales <- list(
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")),
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16))   # filled circle
  )
  
  ## disorder
  rsq = rsq_part_disorder
  
  # sedimentation rate
  vv = paste(": ", round(rsq$partial.rsq[1], 3))
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
  slope_annot = paste("Slope:", round(coefficients(disorder_glm)["s1"], 3))
  v = visreg::visreg(disorder_glm, xvar = "s1", plot = FALSE)
  res <- v$res
  res$group <- df$saturated
  
  disorder_glm_plot_s <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = s1, ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = s1, y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = s1, y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(sedr_label) +
    ylab(disorder_label) +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$S))), y = max(log10(df$disorder)) - 0.2, label = slope_annot, size = annot_size_pts / .pt)
  disorder_glm_plot_s
  # biodiffusion
  vv = paste(": ", round(rsq$partial.rsq[2], 3))
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
  slope_annot = paste("Slope:", round(coefficients(disorder_glm)["m1"], 3))
  v = visreg::visreg(disorder_glm, xvar = "m1", plot = FALSE)
  res <- v$res
  res$group <- df$saturated
  
  disorder_glm_plot_m <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = m1, ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = m1, y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = m1, y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(mix_label) +
    ylab(disorder_label) +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$M))), y = max(log10(df$disorder)) - 0.2, label = slope_annot, size = annot_size_pts / .pt)
  disorder_glm_plot_m
  # mixing depth
  vv = paste(": ", round(rsq$partial.rsq[3], 3))
  annot = as.expression(bquote(partial ~ R^2*  .(vv)))
  slope_annot = paste("Slope:", round(coefficients(disorder_glm)["l1"], 3))
  v = visreg::visreg(disorder_glm, xvar = "l1", plot = FALSE)
  res <- v$res
  res$group <- df$saturated
  
  disorder_glm_plot_l <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = l1, ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = l1, y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = l1, y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(mix_depth_label) +
    ylab(disorder_label) +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$disorder)), label = annot, size = annot_size_pts / .pt) + 
    annotate("text", x = mean(range(log10(df$L))), y = max(log10(df$disorder)) - 0.2, label = slope_annot, size = annot_size_pts / .pt)
  disorder_glm_plot_l
  figure = ggpubr::ggarrange(tavg_glm_plot_s,  tavg_glm_plot_m, tavg_glm_plot_l,
                          disorder_glm_plot_s,  disorder_glm_plot_m, disorder_glm_plot_l,
                          nrow = 2, ncol = 3,
                          labels = LETTERS[1:6],
                          common.legend = TRUE,
                          legend = "bottom")
  return(figure)
}


ggsave("figs/glm_summary_fig.pdf", 
       plot = plot_glm_summary_fig(),
       width = 19,
       units = c("cm"),
       height = 12,
       bg = "white")

#### wd vs tavg, disorder, and ssd (plot ssd histogram) ####

df$ssd_log = log10(df$safe_sampling_distance)

aa = summary(lm_ssd_log)
aa$coefficients |> signif(digits = signif_dig)

ssd_log_label = expression( log[10]* group("(", Safe~sampling~dist.,")") ~ group("[",cm,"]") )

# x axis: data stay on log10 scale, ticks are labelled in metres
wd_lin_label = "Water depth [m]"
wd_x_scale = list(
  scale_x_continuous(breaks = log10(c(10, 100, 1000, 10000)),
                     labels = c("10", "100", "1000", "10000"),
                     minor_breaks = log10(as.vector(outer(1:9, 10^(1:4))))),
  annotation_logticks(sides = "b", outside = FALSE,
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"))
)

signif_dig = 3
# Annotation "adj. R^2: x, p = y" for an lm, with both values to `digits` significant digits.
# The p-value is the overall F-test of the model (= slope p-value for one predictor).
r2_p_annot = function(model, digits = signif_dig){
  s  = summary(model)
  f  = s$fstatistic
  r2 = signif(s$adj.r.squared, digits)
  p  = pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE) |> unname()
  p_txt = format(signif(p, digits), digits = digits)
  as.expression(bquote(adj. ~ R^2 * ": " * .(r2) * "," ~ italic(p) == .(p_txt)))
}

plot_wd_vs_params = function(){
  tavg_slope = lm_tavg_wd$coefficients["wd_l"] |> unname() |> signif(digits = signif_dig)
  disorder_slope = lm_disorder_wd$coefficients["wd_l"] |> unname() |> signif(digits = signif_dig)
  ssd_slope = lm_ssd_log$coefficients["wd_log"] |> unname() |> signif(digits = signif_dig)
  
  sat_scales <- list(
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")),
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16))   # filled circle
  )
  
  yrange = range(log10(df$tavg))
  annot = r2_p_annot(lm_tavg_wd)
  slope_annot = paste("Slope:", tavg_slope)
  
  v = visreg::visreg(lm_tavg_wd, plot = FALSE)
  xname <- v$meta$x            # name of the predictor, as used in v$fit / v$res
  res <- v$res
  res$group <- df$saturated
  
  lm1 <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = .data[[xname]], ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = .data[[xname]], y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = .data[[xname]], y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(wd_lin_label) +
    wd_x_scale +
    ylab(tavg_label) +
    annotate("text",
             y = 0.9 * max(yrange),
             x = mean(range(log10(df$wd))),
             label = annot,
             size = 8/.pt) +
    annotate("text",
             x = mean(range(log10(df$wd))),
             y = 1 * max(yrange),
             label = slope_annot,
             size = 8/.pt) +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  lm1
  
  yrange = range(log10(df$disorder))
  annot = r2_p_annot(lm_disorder_wd)
  slope_annot = paste("Slope:", disorder_slope)
  
  v = visreg::visreg(lm_disorder_wd, plot = FALSE)
  xname <- v$meta$x
  res <- v$res
  res$group <- df$saturated
  
  lm2 <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = .data[[xname]], ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = .data[[xname]], y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = .data[[xname]], y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(wd_lin_label) +
    wd_x_scale +
    ylab(disorder_label) +
    annotate("text",
             y = 0.9 * max(yrange),
             x = mean(range(log10(df$wd))),
             label = annot,
             size = 8/.pt) +
    annotate("text",
             x = mean(range(log10(df$wd))),
             y = 1 * max(yrange),
             label = slope_annot,
             size = 8/.pt) +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  lm2
  
  ## safe sampling distance, log10 scale
  annot = r2_p_annot(lm_ssd_log)
  slope_annot = paste("Slope:", ssd_slope)
  
  v = visreg::visreg(lm_ssd_log, plot = FALSE)
  xname <- v$meta$x
  res <- v$res
  res$group <- df$saturated
  
  # shared y-limits for the regression panel and the histogram
  ylims <- range(c(res$visreg_res, v$fit$visreg_lwr, v$fit$visreg_upr,
                   df$ssd_log), na.rm = TRUE)
  # annotation positions as fractions of the axis range
  # (multiplying a log value, e.g. 0.9 * max, is unreliable when values can be <= 0)
  y_slope = ylims[2]
  y_r2    = ylims[2] - 0.08 * diff(ylims)
  
  lm3 <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = .data[[xname]], ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = .data[[xname]], y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = .data[[xname]], y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(wd_lin_label) +
    wd_x_scale +
    ylab(ssd_log_label) +
    coord_cartesian(ylim = ylims) +
    annotate("text",
             y = y_r2,
             x = mean(range(log10(df$wd))),
             label = annot,
             size = 8/.pt) +
    annotate("text",
             x = mean(range(log10(df$wd))),
             y = y_slope,
             label = slope_annot,
             size = 8/.pt) +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  lm3
  
  hist1 <- df |>
    ggplot(aes(y = ssd_log)) +   # y instead of x + coord_flip()
    geom_histogram(binwidth = 0.1) +
    labs(y = ssd_log_label,
         x = "Count in SMLBase") +
    coord_cartesian(ylim = ylims) +
    geom_hline(yintercept = median(df$ssd_log),   # = log10 of the median SSD
               linetype = "dashed") +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  hist1
  p = ggpubr::ggarrange(lm1, lm2, lm3, hist1, ncol = 2, nrow = 2,
                        labels = LETTERS[1:4],
                        common.legend = TRUE,
                        legend = "bottom")
  return(p)
}

ggsave(filename = "figs/water_depth_vs_params.pdf",
       plot = plot_wd_vs_params(),
       bg = "white",
       width = 19,
       height = 15,
       unit = "cm")

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
    xlab("log10(Sediment Acc. Rate) [cm/a]") +
    ylab("Count")
  p_fmix = ggplot(df, aes(x = F_mix)) + 
    annotate("rect", xmin = -Inf, xmax = F_mix_threshold, ymin = -Inf, ymax = Inf,
             fill = "goldenrod", alpha = 0.25) +
    annotate("rect", xmin = F_mix_threshold, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "mediumpurple", alpha = 0.25) +
    geom_histogram(bins = 30) +
    xlab(expression(F[mix] * " [-]")) +
    ylab("Count") +
    geom_vline(xintercept = F_mix_threshold,
               linetype = "dashed") +
    annotate("text",
             x = 0.5,
             y = 90,
             hjust = 0.5,
             vjust = 0.5,
             label = "Not\nsaturated") +
    annotate("text",
             x = 1.03,
             y = 90,
             angle = 90,
             hjust = 0.5,
             vjust = 0.5,
             label = "Saturated")
  p_fmix
  
  p = ggpubr::ggarrange(p_s, p_m, p_l, p_fmix,
                          labels = LETTERS[1:4],
                          ncol = 2, nrow = 2)
  return(p)
}
ggsave(filename = "figs/sml_overview.pdf",
       plot = plot_sml_overview_fig(),
       bg = "white",
       width = 19,
       height = 18,
       units = "cm")

#### Coefficient plot ####


aa = summary(ssd_glm_std)
aa$coefficients |> signif(digits = signif_dig)
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
  a5 = dwplot(ssd_glm,
              dot_args = dot_args,
              whisker_args = whisker_args) +
    ggtitle("Safe sampling distance") + 
    labs(x = "Regression coefficient", y = "Independent variables") +
    scale_y_discrete(labels = y_axis_labels)+
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_text(angle = 45))
  a6 = dwplot(ssd_glm_std,
              dot_args = dot_args,
              whisker_args = whisker_args) +
    ggtitle("") +
    labs(x = "Beta coefficient", y = "Independent variables") +
    scale_y_discrete(labels = y_axis_labels)+
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_text(angle = 45))
  
  p = ggpubr:: ggarrange(a1, a2, a3, a4, a5, a6,
                         ncol = 2,
                         nrow = 3,
                         labels = letters[1:6])
  return(p)
}

ggsave(filename = "figs/coefficient_plots.pdf",
       plot = make_coeff_plots(),
       width = 19,
       height = 21,
       units = "cm")

#### Plot dimensionless ADD ####
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
  
  fixed_depth_color = "black"
  fixed_age_color = "red"
  transect_width = 2
  df = expand.grid(a, d)
  df$z = mapply(f, df$Var1, df$Var2)
  p1 = ggplot(df, aes(x = Var1, y = Var2, z = z)) +
    geom_contour_filled(breaks = sort(seq(0,1.1, by = 0.1))) +
    geom_contour(aes(z = z), color = "black", linewidth = 0.4,
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
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text = element_text(size = 7)
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
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text = element_text(size = 7)
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
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
  
  p4 = NULL
  
  p = ggpubr::ggarrange(p1, p2, p3, p4,
                        nrow = 2, 
                        ncol = 2,
                        widths = c(2,1),
                        heights = c(2,1),
                        labels = c("a", "b", "b"))
  return(p)
}



ggsave(filename = "figs/dimless_add.pdf",
       plot = plot_dimless_add(),
       height = 8.5,
       width = 8.5,
       units = "cm",
       bg = "white")

## water depth vs sedimentation rate

#### Sediment accumulation rate vs. water depth ####
lm_wd = lm(S_log ~ wd_log, data = df)

summary_lm_wd = summary(lm_wd)
summary_lm_wd$coefficients |> signif(digits = signif_dig)
summary_lm_wd$adj.r.squared |> signif(digits = signif_dig)
f = summary_lm_wd$fstatistic          # named vector: value, numdf, dendf
p_r2 = pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE) |> unname()

# x axis: data stay on log10 scale, ticks are labelled in metres
wd_lin_label = "Water depth [m]"
wd_x_scale = list(
  scale_x_continuous(breaks = log10(c(10, 100, 1000, 10000)),
                     labels = c("10", "100", "1000", "10000"),
                     minor_breaks = log10(as.vector(outer(1:9, 10^(1:4))))),
  annotation_logticks(sides = "b", outside = FALSE,
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"))
)

plot_wd_vs_s = function(){
  s_lm      = summary(lm_wd)
  wd_r2     = s_lm$adj.r.squared |> signif(digits = signif_dig)
  wd_slope  = lm_wd$coefficients["wd_log"] |> unname() |> signif(digits = signif_dig)
  f         = s_lm$fstatistic
  wd_p      = pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE) |> unname()
  p_lab     = if (wd_p < 0.001) "< 0.001" else paste("=", signif(wd_p, signif_dig))
  
  sat_scales <- list(
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")),
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16))   # filled circle
  )
  
  annot = as.expression(bquote(adj. ~ R^2 * ": " * .(wd_r2) * "," ~ italic(p) ~ .(p_lab)))
  slope_annot = paste("Slope:", wd_slope)
  
  v = visreg::visreg(lm_wd, plot = FALSE)
  xname <- v$meta$x
  res <- v$res
  res$group <- df$saturated
  
  ylims <- range(c(res$visreg_res, v$fit$visreg_lwr, v$fit$visreg_upr,
                   df$S_log), na.rm = TRUE)
  
  # annotation positions: below the data, with extra room added at the bottom
  ann_x    <- mean(range(v$fit[[xname]]))
  ann_step <- 0.08 * diff(ylims)
  ann_bot  <- ylims[1] - 0.05 * diff(ylims)
  
  plot_ylims <- c(ylims[1] - 0.22 * diff(ylims), ylims[2])
  
  p <- ggplot() +
    geom_ribbon(data = v$fit, aes(x = .data[[xname]], ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = .data[[xname]], y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = .data[[xname]], y = visreg_res, colour = group, shape = group)) +
    sat_scales +
    xlab(wd_lin_label) +
    wd_x_scale +
    ylab(sedr_label) +
    coord_cartesian(ylim = plot_ylims) +
    annotate("text", x = ann_x, y = ann_bot,
             label = slope_annot, size = 8/.pt) +
    annotate("text", x = ann_x, y = ann_bot - ann_step,
             label = annot, size = 8/.pt) +
    theme(
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
  return(p)
}

ggsave("figs/wd_vs_sed_acc.pdf",
       plot = plot_wd_vs_s(),
       width = 9, height = 8, units = "cm",
       bg = "white", device = cairo_pdf)

ggsave("figs/wd_vs_sed_acc.pdf",
       plot = plot_wd_vs_s(),
       width = 8.5,
       height = 9,
       units = "cm")

#### Save results ####
re.df = data.frame(description = names(results), values = unname(results))
write.csv(re.df, file = "data/res/results.csv")

cat("Done. Results are under \"figs/\"and \"data/res/results.csv\"\n")

#### Figure X: GLM plot, S / M / L vs. safe sampling distance ####
# Paste after the "Figure 3: GLM plot" section (needs s1, m1, l1, df, labels,
# vif_threshold, round_res, results and the plotting constants defined there).

ssd_label = expression( log[10]* group("(", Safe~sampling~distance,")") ~ group("[",cm,"]") )

## Fit on log scale with the same predictors as Fig. 3
# rows with NA/Inf SSD (e.g. G outside the simulated Peclet range) are dropped;
# group labels are carried along so they stay aligned with the residuals
ssd_dat = data.frame(ssd = log10(df$safe_sampling_distance),
                     s1 = s1, m1 = m1, l1 = l1,
                     group = df$saturated) |>
  dplyr::filter(is.finite(ssd))

ssd_glm = glm(ssd ~ s1 + m1 + l1, data = ssd_dat)

vif_ssd = car::vif(ssd_glm)
if (any(vif_ssd > vif_threshold)) stop("Multicolinearity detected")

rsq_part_ssd = rsq::rsq.partial(ssd_glm, adj = TRUE)

# store results, matched by variable name rather than position
pr2 = setNames(rsq_part_ssd$partial.rsq, rsq_part_ssd$variable)
results[c("partial R2 ssd S", "partial R2 ssd D_b", "partial R2 ssd L")] =
  round(pr2[c("s1", "m1", "l1")], round_res)
results[c("coeff ssd S", "coeff ssd D_b", "coeff ssd L")] =
  round(coef(ssd_glm)[c("s1", "m1", "l1")], round_res)

## Generic partial-residual panel (same styling as Fig. 3)
plot_partial_panel = function(model, xvar, xlab, ylab, ylims, rsq_obj, group){
  v = visreg::visreg(model, xvar = xvar, plot = FALSE)
  res = v$res
  res$group = group
  
  rsq_val = rsq_obj$partial.rsq[rsq_obj$variable == xvar]
  vv = paste(": ", round(rsq_val, 3))
  annot = as.expression(bquote(partial ~ R^2 * .(vv)))
  slope_annot = paste("Slope:", round(coef(model)[xvar], 3))
  
  x_mid = mean(range(res[[xvar]]))
  y_top = max(ylims) - 0.05 * diff(ylims)
  dy = 0.1 * diff(ylims)
  
  ggplot() +
    geom_ribbon(data = v$fit, aes(x = .data[[xvar]], ymin = visreg_lwr, ymax = visreg_upr),
                fill = "grey85") +
    geom_line(data = v$fit, aes(x = .data[[xvar]], y = visreg_fit), linewidth = 1) +
    geom_point(data = res, aes(x = .data[[xvar]], y = visreg_res, colour = group, shape = group)) +
    scale_colour_manual(name = NULL,
                        values = c("Not saturated" = "goldenrod",
                                   "Saturated"     = "mediumpurple")) +
    scale_shape_manual(name = NULL,
                       values = c("Not saturated" = 17,    # filled triangle
                                  "Saturated"     = 16)) + # filled circle
    xlab(xlab) +
    ylab(ylab) +
    coord_cartesian(ylim = ylims) +   # zooms without dropping points (unlike ylim())
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts)) +
    annotate("text", x = x_mid, y = y_top,      label = annot,       size = annot_size_pts / .pt) +
    annotate("text", x = x_mid, y = y_top - dy, label = slope_annot, size = annot_size_pts / .pt)
}

plot_ssd_glm_fig = function(){
  ssd_y_axis_lims = range(c(ssd_dat$ssd + log10(1.2), ssd_dat$ssd + log10(0.85)))
  
  p_s = plot_partial_panel(ssd_glm, "s1", sedr_label,      ssd_label, ssd_y_axis_lims, rsq_part_ssd, ssd_dat$group)
  p_m = plot_partial_panel(ssd_glm, "m1", mix_label,       ssd_label, ssd_y_axis_lims, rsq_part_ssd, ssd_dat$group)
  p_l = plot_partial_panel(ssd_glm, "l1", mix_depth_label, ssd_label, ssd_y_axis_lims, rsq_part_ssd, ssd_dat$group)
  
  ggpubr::ggarrange(p_s, p_m, p_l,
                    nrow = 1, ncol = 3,
                    labels = LETTERS[1:3],
                    common.legend = TRUE,
                    legend = "bottom")
}

ggsave("figs/ssd_glm_summary_fig.pdf",
       plot = plot_ssd_glm_fig(),
       width = 19,
       height = 7,
       units = "cm",
       bg = "white")


