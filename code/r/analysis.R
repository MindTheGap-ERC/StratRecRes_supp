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

df = dplyr::filter(df, !(df$S_min <= 0 | df$S_max <= 0 | df$M_min <= 0 | df$M_max <= 0))
df$S = 0.5 * (df$S_max + df$S_min)
df$M = 0.5 * (df$M_max + df$M_min)
df$L = 0.5 * (df$L_max + df$L_min)
df$Pe_min = df$L_min * df$S_min / df$M_max
df$Pe_max = df$L_max * df$S_max / df$M_min
df$Pe = df$L * df$S / df$M

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
  tavg = approx(x = peclet_numbers, y = tavg_dimless, xout = Pe, rule = 2)$y * L
  return(tavg)
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

signif(quantile(log10(df$Pe), c(0.025, 0.975)), 3)
signif(median(log10(df$Pe)), 3)
signif(quantile(log10(df$Pe), c(0.25, 0.75)), 3)

signif(quantile(df$F_mix, c(0.25, 0.75), 3))
df$F_mix |> median() |> signif(3)

#### GLM ####

m1 = log10(df$M)
s1 = log10(df$S)
l1 = log10(df$L)
ta = log10(df$tavg)
di = log10(df$disorder)

tavg_glm = glm(ta ~ m1 + s1 + l1)
disorder_glm = stats::glm(di ~ m1 + s1 + l1)

tavg_glm
summary(tavg_glm)
step(tavg_glm)
step(disorder_glm)
#fl = visreg::visreg(tavg_glm,ylim=c(0,4.4))




tavg_glm_plot_s

rsq::rsq.partial(g,adj=TRUE)

g = stats::glm(di ~ m1 + s1 + l1)

g
summary(g)
stats::step(g)
visreg::visreg(g)

visreg::visreg(g,"s1")

rsq::rsq.partial(g,adj=TRUE)

#### Figure 2 ####
# plots histogram of Peclet numbers and F_mix
plot_hist_and_fmix = function(){
  pe_min_log = min(log10(peclet_numbers))
  pe_max_log = max(log10(peclet_numbers))
  pe_step = 1/3
  tavg_lwd = 3
  tavg_col = "black"
  a = hist(log10(df$Pe), plot = FALSE, freq = FALSE, breaks = seq(pe_min_log, pe_max_log, pe_step))
  png("figs/fig2.png")
  plot(log10(peclet_numbers), tavg_dimless, type = "l",
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
}


#### Supplementary Figure 1 ####
plot_fmix_dist = function(){
  ggplot(df, aes(x = F_mix)) + geom_histogram()
  ggsave("figs/supp_fig_1.png")

}

plot_fmix_dist()

#### Supplementary figure 2 ####

plot_sml_histograms = function(fig_name){
  p_m = ggplot(df, aes(x = M)) + geom_histogram() + scale_x_log10()
  p_l = ggplot(df, aes(x = L)) + geom_histogram()
  p_s =ggplot(df, aes(x = S)) + geom_histogram() + scale_x_log10()
  
  figure = egg::ggarrange(p_m, p_l, p_s, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
  ggsave(paste0("figs/", fig_name, ".png"), figure)
}

plot_sml_histograms("suppl_fig_2")


#### Figure 1 ####
# Data from Tomasovych et al 2018
to = read.csv("data/tomasovych_et_al_2018/S0094837318000222sup002.csv", sep = ";")
to2 = dplyr::filter(to, Station == "Po4 M21" & Outlier == "")

df_te = data.frame("age" = to2$Age[!is.na(to2$Age)], "depth" = as.numeric(to2$max..depth..cm.)[!is.na(to2$Age)])

age_lim = 150
depth_lim = 120

ggplot(df_te, aes(x = age, y = depth)) + 
  geom_point() +
  xlim(0, age_lim) +
  scale_y_reverse(lim = c(depth_lim, 0)) +
  geom_density2d_filled(alpha = 0.5) +
  ggtitle("Core Po4 Particle distribution") +
  xlab("Age [years]") +
  ylab("Depth [cm]")
ggsave("figs/tomasovych_et_al_2018_density.png")

#df_te$depth
#IQR(df_te$age[df_te$depth == 95])

##
# read data from matlab
da = R.matlab::readMat("code/matlab/Po4_res.mat")

ages = da$ages[1,]
depths = da$depths[1,]
u = da$u

#dimnames(u) <- list( "depths" = depths, "ages" = ages)
mydf = reshape2::melt(u)


 p = ggplot(mydf, aes(Var2, Var1, z = value)) + geom_contour(binwidth = 0.001, show.legend = FALSE) + xlim(range(mydf$Var2)) + 
   ylim(range(mydf$Var1)) +
   xlim(0, age_lim) +
   scale_y_reverse(lim = c(depth_lim, 0)) +
   xlab("Age [years]") +
   ylab("Depth [cm]")
p 
ggsave("figs/Po4_modeled.png")


#### GLM plot ####
sedr_label = expression( log[10]*"("*Sed*"."*~rate*")"*  ~  group("[",cm/a,"]") )
mix_label =  expression( log[10]* group("(", Mixing~int*".",")")  ~  group("[",cm/a,"]") ) #"log10(Mixing int.) [cm^2/a]"
mix_depth_label =  expression( log[10]* group("(", Mixing~depth,")")  ~  group("[",cm,"]") )  #"log10(Mixing depth) [cm]"
tavg_label =  expression( log[10]* group("(", Time ~averaging,")")  ~  group("[",a,"]") ) # "log10(Time averaging) [a]"
disorder_label = expression( log[10]* group("(", Stratigraphic~disorder,")")  ~  group("[",cm,"]") )  #"log10(Stratigraphic disorder) [cm]"

tavg_y_axis_lims = range(log10(c(df$tavg* 1.1, df$tavg * 0.9)))
disorder_y_axis_lims = range(log10(c(1.1 * df$disorder, 0.9 * df$disorder)))

two_col_width_cm = 12.28
x_axis_text_size_pts = 7
y_axis_text_size_pts = 7
glm_plot_height_cm = 14

make_glm_plot = function(fig_name){
  tavg_glm_plot_s = visreg::visreg(tavg_glm,xvar = "s1",
                                   gg = TRUE) +
    xlab(sedr_label) +
    ylab(tavg_label) +
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))
  
  tavg_glm_plot_m = visreg::visreg(tavg_glm,xvar = "m1",
                                   gg = TRUE) +
    xlab(mix_label) +
    ylab(tavg_label)  + 
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))
  
  tavg_glm_plot_l= visreg::visreg(tavg_glm,xvar = "l1",
                                  gg = TRUE) +
    xlab(mix_depth_label) +
    ylab(tavg_label)  + 
    ylim(tavg_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))
  
  
  disorder_glm_plot_s = visreg::visreg(disorder_glm,xvar = "s1",
                                       gg = TRUE) +
    xlab(sedr_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))
  
  disorder_glm_plot_m = visreg::visreg(disorder_glm,xvar = "m1",
                                       gg = TRUE) +
    xlab(mix_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))
  
  disorder_glm_plot_l= visreg::visreg(disorder_glm,xvar = "l1",
                                      gg = TRUE) +
    xlab(mix_depth_label) +
    ylab(disorder_label)  +
    ylim(disorder_y_axis_lims) +
    theme(axis.title.x = element_text(size = x_axis_text_size_pts),
          axis.title.y = element_text(size = y_axis_text_size_pts))
  
  figure = egg::ggarrange(tavg_glm_plot_s, tavg_glm_plot_l, tavg_glm_plot_m,
                          disorder_glm_plot_s, disorder_glm_plot_l, disorder_glm_plot_m,
                          nrow = 2, ncol = 3,
                          labels = LETTERS[1:6])
  
  ggsave(paste0("figs/", fig_name, ".png"), figure,
         width = two_col_width_cm, units = c("cm"),
         height = glm_plot_height_cm)
}

make_glm_plot("glm_res")
