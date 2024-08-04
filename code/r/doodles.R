source(file = "code/R/utils.R")

import_ADDs_from_matlab = FALSE
if (import_ADDs_from_matlab) source(file = "code/R/import_ADD_to_R.R")

load_ADDs = TRUE
if(load_ADDs) load(file = "data/R_outputs/ADDs.RData")

#### extract data from adds ####

df_dimless = data.frame(peclet_no = peclet_numbers)
below_sml = min(which(d_dimless > 1))
tavgs = rep(NA, length(peclet_numbers))
for (i in seq_along(ADD_list)){
  add = ADD_list[[i]]$ADD
  cdf = cumul_int_trap(x = t_dimless, y = add[ , below_sml])
  tavg = iqr_from_cdf(cdf, t_dimless)
  tavgs[i] = tavg
  df_dimless$tavg = tavgs
}



#### Make contour plot ####

add_no = 12
pe = ADD_list[[add_no]]$peclet_number

maxval = max(add)

contour(y = d_dimless,
        x = t_dimless,
        z = add,
        levels = seq(0,1,length.out = 10),
        xlim = c(0,1))

melt(add)

a = as.data.frame(add)

ggplot(a, aes = )

# some plots
i = 1
below_sml = min(which(d_dimless > 1))
plot(
  x = t_dimless,
  y = ADD_list[[i]]$ADD[,below_sml],
  type = "l",
  main = paste(c("Peclet number: ", ADD_list[[i]]$peclet_number))
)

plot(NULL,
     xlim = range(t_dimless),
     ylim = c(0,9))
for ( i in seq_along(peclet_numbers)){
  lines(
    x = t_dimless,
    y = ADD_list[[i]]$ADD[,below_sml]
  )
}




plot(x = t_dimless,
     y = cumul_int_trap(
       x = t_dimless,
       y = ADD_list[[i]]$ADD[,below_sml]
    ),
     type = "l"
)

for (i in seq_along(peclet_numbers)){
  ADD = ADD_list[[i]]$ADD
  
  ADD_list[[i]]$cdf_bottom = cumul_int_trap(
    x = t_dimless,
    y = pmax(ADD[,below_sml],rep(0,length.out = length(t_dimless))),
    normalize = TRUE
    )
  ADD_list[[i]]$mat_sml = sapply(
    X = seq_along(t_dimless),
    FUN = function(j) int_trap(x = d_dimless, y = ADD[j,])
  )
  ADD_list[[i]]$tavg = iqr_from_cdf(
    cdf = ADD_list[[i]]$cdf_bottom,
    x = t_dimless
    )
  ADD_list[[i]]$disorder = iqr_from_cdf(
    cdf = 1-ADD_list[[i]]$mat_sml,
    x = t_dimless
    )
}

i = 20
plot(
  x= t_dimless,
  y = ADD_list[[i]]$cdf_bottom,
  type = "l"
  )

plot(
  x = t_dimless,
  y = ADD_list[[i]]$mat_sml,
  type = "l"
)


tavg_dimless = sapply(seq_along(peclet_numbers), function(i) ADD_list[[i]]$tavg)
disorder_dimless = sapply(seq_along(peclet_numbers), function(i) ADD_list[[i]]$disorder)

plot(log10(peclet_numbers),tavg,type = "l",ylim = c(0,max(tavg)))
plot(log10(peclet_numbers),disorder,type = "l",ylim = c(0,max(disorder)))

# numerics!
plot(log10(peclet_numbers),tavg / disorder,type = "l",ylim = pmax(range(tavg / disorder),c(0,0)))


#### read data from smlbase ####

SMLBase = read.csv(
  file = "data/smlbase_raw/SMLBase v1_03_temp.csv",
  sep = ","
)

S_max = pmax(SMLBase$SMin, SMLBase$SMax, SMLBase$S, na.rm = TRUE)
S_min = pmin(SMLBase$SMin, SMLBase$SMax, SMLBase$S, na.rm = TRUE)

L_max = pmax(SMLBase$L1Max, SMLBase$L1Min, SMLBase$L1, SMLBase$L2, na.rm = TRUE)
L_min = pmin(SMLBase$L1Max, SMLBase$L1Min, SMLBase$L1, SMLBase$L2, na.rm = TRUE)

M_max = pmax(SMLBase$M1, SMLBase$M1Max, SMLBase$M1Min, SMLBase$M2, na.rm = TRUE)
M_min = pmin(SMLBase$M1, SMLBase$M1Max, SMLBase$M1Min, SMLBase$M2, na.rm = TRUE)

M_use = runif(length(M_min),M_min,M_max)
L_use = runif(length(L_min),L_min,L_max)
S_use = runif(length(S_min),S_min,S_max)

Pe_max = L_max * S_max / M_min
max(Pe_max, na.rm = TRUE)

Pe_min = L_min * S_min/ M_max
hist(log10(Pe_min))

#### Check SMLBase ####
## Check SMLBase at these locatinos:
# S is 0
SMLBase$boxModelID[SMLBase$S <= 0]

SMLBase$boxModelID[S_max <= 0]

SMLBase$boxModelID[L_max - L_min > 10]

SMLBase$boxModelID[M_max - M_min > 10]



tavg_fun = function(M,S,L){
  peclet_no = L * S / M
  tavg = approx(x = peclet_numbers, y = tavg_dimless, xout = peclet_no)$y * L/S
  return(tavg)
}
disorder_fun = function(M,S,L){
  peclet_no = L * S / M
  tavg = approx(x = peclet_numbers, y = disorder_dimless, xout = peclet_no)$y * L
  return(tavg)
}

tavg_fun(
  M = 1,
  S = 2,
  L = 3
     )

disorder_fun(
  M = 1,
  S = 2,
  L = 3
)

ind = S_min > 0



tavg_a = tavg_fun(M,S,L)
disorder_cm = disorder_fun(M,S,L)
hist(log10(tavg_a))
hist(disorder_cm)

hist(log10(L*S/M))

m1 = log10(M)
s1 = log10(S)
l1 = log10(L)
ta = log10(tavg_a)
di = log10(disorder_cm)

g = glm(ta ~ m1 + s1 + l1)

g
summary(g)
step(g)
visreg::visreg(g,ylim=c(0,4.4))

visreg(g,"S")

rsq.partial(g,adj=TRUE)

g = glm(di ~ m1 + s1 + l1)

g
summary(g)
step(g)
visreg::visreg(g)

visreg::visreg(g,"s1")

rsq.partial(g,adj=TRUE)



tavg_a = tavg_fun(M_use,S_use,L_use)
disorder_cm = disorder_fun(M,S,L)
hist(log10(tavg_a))
hist(disorder_cm)

hist(log10(L*S/M))

m1 = log10(M)
s1 = log10(S)
l1 = log10(L)
ta = log10(tavg_a)
di = log10(disorder_cm)

g = glm(ta ~ m1 + s1 + l1)
