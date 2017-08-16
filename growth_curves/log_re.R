library(data.table)
library(dplyr)
library(drc)
library(ggplot2)

name <- "short"
time <- c(0,30,60,90,120,150,180)
growth_w <- c(6000000.0,16000000.0,18000000.0,410000000.0,4600000000.0,6200000000.0,4900000000.0)
growth_a <- c(6000000.00,16000000.00,18000000.00,410000000.00,4600000000.00,2200000000.00,4300000000.00)
#
# name <- "long"
# time <- c(0,30,60,90,120,150, 180, 210,240,270,300)
# growth_w <- c(6.00e+06, 1.60e+07, 1.80e+07, 4.10e+08, 4.60e+09, 6.20e+09, 4.90e+09, 5.60e+10, 1.56e+10, 2.00e+10, 4.20e+10)
# growth_a <- c(6.0e+06, 1.6e+07, 1.8e+07, 4.1e+08, 4.6e+09, 2.2e+09, 4.3e+09, 9.2e+09, 2.8e+10, 4.8e+09, 5.7e+09)

growth = list(growth_w, growth_a)
data.long = list()
data.fits = list()
legend = c("With Induction", "After Induction")

for (i in 1:2) {
  d = data.table(X=time, Y=growth[[i]])
  d.long <- reshape2::melt(d,id.vars = "X") #reshapes the demo dataset to long format
  d.long.LL.5 <- drm(data = d.long, value~X, fct=LL.5(c(NA, NA, NA, NA, NA)), na.action = na.omit, type="continuous")
  d.fits <- expand.grid(conc=seq(min(time), max(time), length=100))
  pm <- predict(d.long.LL.5, newdata=d.fits, interval="confidence") 
  d.fits$p <- pm[,1]
  d.fits$pmin <- pm[,2]
  d.fits$pmax <- pm[,3]
  d.long$XX <- d.long$X
  d.long$XX[d.long$XX == 0] <- 1.00e-09
  data.long[[i]] = cbind(d.long, "bg" = legend[i])
  data.fits[[i]] = cbind(d.fits, "bg" = legend[i])
}
data.long <- do.call(rbind,data.long)
data.fits <- do.call(rbind,data.fits)

weird <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log(abs(x)),
                           inverse=function(x) sign(x)*exp(abs(x)))

g <- ggplot(data.long, aes(x = XX, y = value, color=bg, shape=bg)) +
  geom_point(size=2, show.legend = FALSE) +
  geom_line(data=data.fits, aes(x=conc, y=p, color=bg), size = 1) +  
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(breaks = round(seq(min(time), max(time), by = 30),1)) +
  labs(title = paste0("Bacterial Growth Curve (",min(time), "-",max(time),"s)"), 
       x="Time (min)", 
       y="Millions of Colony Forming Units per mL",
       color="Bacterial Growth") +
  scale_colour_discrete("Bacterial Growth") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0("gc_", name, ".png"), g, width = 9, height = 4)

# p1 <- ggplot(data.long, aes(x = XX, y = value, color=bg)) + 
#   geom_point() + 
#   scale_x_continuous(breaks = round(seq(min(time), max(time), by = 30),1)) +
#   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
#   labs(title = "Bacterial Growth Curves", x="Time (min)", y="Millions of Colony Forming Units per mL") + 
#   scale_colour_discrete(name="Bacterial Growth")
# 
# p2 <- ggplot(data.fits, aes(x =conc, y = p, color=bg)) +
#   geom_line(size = 1.2) + 
#   scale_x_continuous(breaks = round(seq(min(time), max(time), by = 30),1)) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
# 
# g1 <- ggplotGrob(p1)
# g2 <- ggplotGrob(p2)
#scale_y_continuous(trans=weird) +
#scale_y_log10(limits=c(10^5, max(c(d.fits$p, d.long$value)))) +



#plot(drm(growth ~ time, data = d, fct = LL.5(c(NA, NA, NA, NA, NA))),log="y",ylim=c(1, 10^10))


