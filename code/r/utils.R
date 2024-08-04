int_trap = function(x,y){
  #' @title integral via trapezoidal rule
  #' 
  #' @description integrates a piecewise linear function using
  #' 
  #' @param x abscissa, strictly increasing
  #' @param y ordinate
  #' 
  #' @returns numeric, integral over approxfun(x,y)
  #' 
  xdiff =diff(x)
  ymean = 0.5 * (y[1:(length(y)-1)] + y[2:length(y)])
  return(sum(xdiff * ymean))               
}

cumul_int_trap = function(x,y, normalize = FALSE){
  #' @title cumulative integral via trapezoidal rule
  #' 
  #' @description calculated the cumulative integral over piecewise linear function
  #' 
  #' @param x abscissa, strictly increasing
  #' @param y ordinate
  #' @param normalize logical, should the integral be normalized to 1? (for usage with cdfs)
  #' 
  #' @returns numeric vector of same length as x
  xdiff =diff(x)
  ymean = 0.5 * (y[1:(length(y)-1)] + y[2:length(y)])
  r = cumsum(c(0,xdiff * ymean))
  if (normalize) {r = r/max(r)} 
  return(r)   
  
}

iqr_from_cdf = function(cdf,x){
  #' @title interquartile range from cumulative distribution function
  #' 
  #' @param cdf vector, values of cdf
  #' @param x abscissa of cdf
  #' 
  #' @returns numeric, interquartile range of the cdf
  f2 = approxfun(x = x,
                 y = cdf - 0.25)
  f1 = approxfun(x = x,
                 y = cdf - 0.75)
  iqr = uniroot(f1, range(x))$root - uniroot(f2,range(x))$root
  return(iqr)
}
