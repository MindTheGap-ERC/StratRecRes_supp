# read data from matlab outputs.
cat("Importing matlab data into R, please wait.\n")
matlab_outputs = R.matlab::readMat("data/matlab_outputs/tavg_from_matlab.mat")

## remove additional dimensions
peclet_numbers = matlab_outputs$peclet.numbers[1,]
t_dimless = matlab_outputs$t.dimless[1,]
epsilon = matlab_outputs$epsilon[1,1]
ode_solver_options = matlab_outputs$odeoptions

tavg_list = vector(mode = "list", length = length(peclet_numbers))
for (i in seq_along(peclet_numbers)){
  tavg_list[[i]] = list(
    den = matlab_outputs$tavg.below.sml[, i],
    peclet_number = matlab_outputs$peclet.numbers[1,i]
  )
}

remove(i,matlab_outputs)
var_names = c("peclet_numbers",
              "t_dimless",
              "epsilon",
              "ode_solver_options",
              "tavg_list")
save(
  list = var_names,
  file = "data/r_outputs/tavg_below_sml.Rdata"
)

remove(list = c(var_names,"var_names"))

cat("Done. Data is saved in /data/r_outputs/tavg_below_sml.Rdata\n")
