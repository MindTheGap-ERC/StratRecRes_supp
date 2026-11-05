df2 = df[df$G < 1.56, ]

df3 = df2[,c("id", "S", "M", "L", "wd", "Pe", "G")]
summary(df3)
df
smlbase$long2[, df3$]
write.csv(df3, file = "data/non_mixed.csv")

load(file = "data/r_outputs/matlab_res.Rdata")

i = 30

fval = tavg_list[[i]]
pe = peclet_numbers[i]
source(file = "code/R/utils.R")

plot(t_dimless, fval$den)

aa = cumul_int_trap(t_dimless, fval$den, normalize = TRUE)
max(aa)
plot(t_dimless, aa)
aa[1
  ]
den = c(0, diff(aa)/diff(t_dimless))
plot(t_dimless, den)
CDF = approxfun(x = t_dimless, y = aa, yleft = 0, yright = 1)
plot(t_dimless, CDF(t_dimless), type = "l")

PDF = approxfun(x = t_dimless, y = den, yleft = 0, yright = 0)

comb_f = function(y){
  comb_fun = function(x){PDF(x)* CDF(x-y)}
  return(integrate(comb_fun, lower = - Inf, upper = Inf)[[1]])
}

comb_f(1)

depths = seq(0, 3, by = 0.02)
plot(depths, sapply(depths, comb_f), ylim = c(0, 0.5), type = "l")

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
df_p |>
  ggplot(aes(x = L_dimless, y = logG, z = P)) +
  geom_contour_filled(breaks = c(0, 0.1), fill = "aquamarine", alpha = 0.3) +
  geom_contour_filled(breaks = c(0.1, 0.5), fill = "tomato", alpha = 0.3) +
  geom_contour(breaks = setdiff(levels, 0.1), colour = "grey30") +
  geom_contour(breaks = 0.1, colour = "black", linewidth = 1.2) +   # highlighted contour
  geom_text_contour(breaks = levels, skip = 0, stroke = 0.2) +
  annotate("text", x = 2.2, y = 3, label = "Low reversal risk",
           colour = "black", fontface = "bold", size = 4) +
  annotate("text", x = 0.8, y = 3, label = "High reversal risk",
           colour = "black", fontface = "bold", size = 4) +
  labs(title = "Probability of age reversal",
       x = "Sampling distance [proportion SML depth]",
       y = "log10 G") +
  coord_flip()

# find sampling distance where prob of reversal is below 0.1
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
df$safe_sampling_distance = lim_fun(log10(df$G)) * df$L
df$safe_sampling_distance |> hist()

plot(df$safe_sampling_distance, df$wd)
wd = df$wd
wd1 = log10(df$wd)
ssd = df$safe_sampling_distance
ttt = lm(ssd ~ wd1)
summary(ttt)
oo = visreg::visreg(ttt)
typeof(oo)
