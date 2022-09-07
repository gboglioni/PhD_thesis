# Individual Risk Model with three PIBD variables

# Load packages
library(ggplot2)
library(gridExtra)
# Set "rate" of the Exponentials
l <- 1

# Independent case
df  <- function(x){
  return(3*dgamma(x,shape=1, rate=l)/8 + 3*dgamma(x,shape=2, rate=l)/8
         + dgamma(x,shape=3, rate=l)/8)}
cdf <- function(x){
  return(1/8*pbinom(x+1,1,1) + 3*pgamma(x,shape=1, rate=l)/8 
         +3*pgamma(x,shape=2, rate=l)/8 + pgamma(x,shape=3, rate=l)/8)}

# PIBD case
df.p  <- function(x){3*dgamma(x,shape=1, rate=l)/4 + dgamma(x,shape=3, rate=l)/4}
cdf.p <- function(x){3*pgamma(x,shape=1, rate=l)/4 + pgamma(x,shape=3, rate=l)/4}
x      <- seq(0.01, 7.5, by = 0.01)
pdf.x  <- df(x)
cdf.x  <- cdf(x)
pdf.p.x  <- df.p(x)
cdf.p.x  <- cdf.p(x)

#CDF comparison
df <- data.frame(x, pdf.x, cdf.x, pdf.p.x, cdf.p.x)

#PDF comparison
mp=data.frame(x=c(0), y=c(0), vx=c(0), vy=c(1/8))

df.plot <- ggplot()+
  geom_line(data=df,aes(x, y=pdf.x, colour="darkblue"),size=1.25)+
  geom_line(data=df,aes(x, y=pdf.p.x, colour="red"),size=1.25)+
  geom_segment(data=mp, mapping=aes(x=x, y=y, xend=vx, yend=vy), 
  size=1.5, color="#F8766D")+  geom_point()+
  geom_point(aes(x=0, y=1/8), size = 3, colour="#F8766D")+
  labs(y= "f(x)", x = "x") + theme_grey(base_size = 16)+
  theme(legend.position = "none")+ theme(legend.title = element_blank())


#CDF
cdf.plot <- ggplot()+
  geom_line(data=df,aes(x, y=cdf.x, colour="darkblue"),size=1.25)+
  geom_line(data=df,aes(x, y=cdf.p.x, colour="red"),size=1.25)+
  geom_point()+
  geom_point(aes(x=0, y=1/8), size = 3, colour="#F8766D")+
  geom_point(aes(x=0, y=0), size = 3, colour="#00BFC4")+
  scale_color_discrete(labels = c("mutual indep.", "PIBD"))+ 
  labs(y= "F(x)", x = "x") + theme_grey(base_size = 16)+ 
  theme(legend.position = c(0.75, 0.1)) + theme(legend.title = element_blank())

# Put two plots side by side
grid.arrange(df.plot, cdf.plot, ncol=2)
