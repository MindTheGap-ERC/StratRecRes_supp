load(file = "data/r_outputs/matlab_res.Rdata")
i = 25
x = t_dimless
y = tavg_list[[i]]$den
plot(x, y, type = "l")


f = function(a, d){
  z = approx(x = d - 1 + t_dimless- 0.3, y = y, xout = a, rule = 2)$y
  z[d < 1] = approx(x = t_dimless - 0.3, y = y, xout = a, rule = 2)$y
  return(z)
}

a = seq(0, 4, by = 0.05)
d = seq(0, 4, by = 0.05)
df = expand.grid(a, d)
df$z = mapply(f, df$Var1, df$Var2)
p = ggplot(df, aes(x = Var1, y = Var2, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(low = "black", high = "white", trans = "sqrt") +
  scale_y_reverse() +
  labs(x = "Age",
       y = "Depth") +
  theme(legend.position = "none")
ggsave(filename = "figs/whitescale_heights.tiff",
       plot = p)
