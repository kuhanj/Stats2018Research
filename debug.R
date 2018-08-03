## tests for rma.exact tensor and optimization versions
########################################

source("metafunction.R")

K=10
c0=4
tau20=12.5
varhat=(seq(1, 5, length=K))^2
mu0=0
c0 <- 1
ystar=rnorm(K)*sqrt(varhat+tau20)+mu0
Z <- matrix(rnorm(K*1e3),nrow=K)
yi <- ystar
vi <- varhat

## 1. try new contouring function at a point (mu,tau2)
Z <- matrix(rnorm(K*1e3),nrow=K)
mu <- -5
tau2 <- 1
c0 <- 3
rma.exact2(mu,tau2,ystar,varhat,Z,c0=c0)
rma.test(mu,tau2,ystar,varhat,Z,c0=c0)

cov.new=mean(apply(result.exact-mu0,1,prod)<0)
cov.tau=mean(apply(tau.range-tau20, 1, prod)<0)
length.new=mean(result.exact[,2]-result.exact[,1])
print(c(cov.new, cov.dl, cov.tau, length.new, length.dl))




## tau.range=matrix(0, 1000, 2)
## result.exact=result.dl=matrix(0, 1000, 2)
## for(b in 1:1000)
##    {set.seed(b)
##     ystar=rnorm(K)*sqrt(varhat+tau20)+mu0
##     fit.dl=rma.uni(ystar, varhat, method="DL")
##     result.dl[b,]=c(fit.dl$ci.lb, fit.dl$ci.ub)
##     }

## cov.dl=mean(apply(result.dl-mu0,1,prod)<0)
## length.dl=mean(result.dl[,2]-result.dl[,1])
## print(c(cov.dl, length.dl))


## for(b in 1:1000)
##    {set.seed(b)
##     ystar=rnorm(K)*sqrt(varhat+tau20)+mu0
##     fit.dl=rma.uni(ystar, varhat, method="DL")
##     result.dl[b,]=c(fit.dl$ci.lb, fit.dl$ci.ub)

##     fit=rma.exact(ystar, varhat, B=5000, num.mu=40, num.tau=200, c0=c0)
##     result.exact[b,]=fit$ci
##     tau.range[b,]=range(fit$tau.grd)
##     print(b)
##     print(fit$precision)
##     print(mean(apply(result.exact[1:b,,drop=F]-mu0,1,prod)<0))
##     plot(c(min(fit$mu.min), max(fit$mu.max)), range(fit$tau.grd), type="n", xlab="mu",  ylab="tau^2")
##     lines(fit$mu.min, fit$tau.grd)
##     lines(fit$mu.max, fit$tau.grd)
##     lines(rep(mu0, 2), c(0, 1000), col=2, lwd=2)
##     lines(c(-100,100), rep(tau20,2), col=3, lwd=2)
##     lines(c(-100,100), rep(fit.dl$tau2, 2), col=4)
##     }


## try new contouring function

tau2 <- seq(1,100,length.out=1e2)
mu <- c(seq(-5,-3,length.out=1e2),seq(-1,1,length.out=5e1))
mu <- seq(-5,-3,length.out=1e2)
mu <- seq(-5,1,length.out=1e2)
p.vals <- rma.exact2(mu,tau2,yi=ystar,vi=varhat,Z=Z)
p.vals2 <- rma.exact.tensor(mu,tau2,yi=ystar,vi=varhat,Z=Z)
system.time(rma.exact2(mu,tau2,yi=ystar,vi=varhat,Z=Z))
system.time(rma.exact.tensor(mu,tau2,yi=ystar,vi=varhat,Z=Z))
system.time(rma.exact.fast(mu,tau2,yi=ystar,vi=varhat,Z=Z))
p.vals <- rma.exact.fast(mu,tau2,yi=ystar,vi=varhat,Z=Z,c0=c(1,3))
plot(p.vals[[2]])
system.time(rma.exact(thetahat=yi,varhat=vi,B=1e3,num.tau=100,num.mu=100))
p.vals['-4.7979798','18']
plot(p.vals)
confint(p.vals,.05)$length

mu <- rnorm(1); tau2 <- abs(rnorm(1))
## rma.exact(thetahat=yi,varhat=vi,)


ggplot(df,aes(x=mu,y=tau2,z=p))+stat_contour()
p.vals <- rma.exact2(mu=seq(-5,-3,length.out=1e2),tau2,yi=ystar,vi=varhat,Z=matrix(rnorm(K*1e3),nrow=K)) 
p.vals2 <- rma.exact2(mu=seq(-1,1,length.out=1e2),tau2,yi=ystar,vi=varhat,Z=matrix(rnorm(K*1e3),nrow=K))

plot(p.vals)
contour(mu,tau2,p.vals,levels=(1:20)/100)
contour(mu,tau2,p.vals,levels=.05)
contour(mu,tau2,p.vals)
image(mu,tau2,p.vals,col=topo.colors(12))
contour(mu,tau2,p.vals,add=T,levels=(1:10)/100)
contour(mu,tau2,p.vals,add=T,levels=.05)

p.vals2 <- p.vals

uu <- expand.grid(as.numeric(rownames(p.vals)),as.numeric(colnames(p.vals)))
pp <- data.frame(p.val=as.numeric(p.vals),mu=as.numeric(uu$Var1),tau2=as.numeric(uu$Var2))
uu2 <- expand.grid(as.numeric(rownames(p.vals2)),as.numeric(colnames(p.vals2)))
pp2 <- data.frame(p.val=as.numeric(p.vals2),mu=as.numeric(uu2$Var1),tau2=as.numeric(uu2$Var2))
pp3 <- rbind(pp,pp2)
pp3 <- summarize(group_by(pp3,mu,tau2),p.val=mean(p.val))
pp3 <- as.data.frame(arrange(pp3,tau2,mu))
mu <- unique(pp3$mu)
tau2 <- unique(pp3$tau2)
p.vals3 <- matrix(pp3$p.val,nrow=length(mu),ncol=length(tau2))
p.vals3 <- p.vals+p.vals2
rownames(p.vals3) <- mu
colnames(p.vals3) <- tau2
contour(mu,tau2,p.vals3)


## 2. get bounds for mu, tau2
K <- 10
c0 <- 2
tau2 <- 12.5
vi <- (seq(1, 5, length=K))^2
mu0=0
yi=rnorm(K)*sqrt(vi+tau2)+mu0
Z <- matrix(rnorm(K*1e3),nrow=K)


fit=rma.uni(yi=yi, vi=vi, method="DL")
mu.bounds <- fit$b+c(-1,1)*fit$se*qnorm(.999)
tau2.bounds <- tau.ci(thetahat=yi,varhat=vi)
mu <- seq(mu.bounds[1],mu.bounds[2],length.out=5e2)
tau2 <- seq(tau2.bounds[1],tau2.bounds[2],length.out=5e2)
rma0 <- rma.exact.fast(mu,tau2,yi=yi,vi=vi,Z=Z)
## using defaults
rma0 <- rma.exact.fast(yi=yi,vi=vi,Z=Z)
plot(rma0,levels=c(.01,.05,.1,.15,.2))
plot(rma0)

## 3. adding rma objects to improve resolution, avoiding mem overflow
source('metafunction.R')
K <- 10
c0 <- 2
tau2 <- 12.5
vi <- (seq(1, 5, length=K))^2
mu0=0
yi=rnorm(K)*sqrt(vi+tau2)+mu0
Z <- matrix(rnorm(K*1e3),nrow=K)

rma0 <- NULL
rma0 <- rma.exact.fast(yi=yi,vi=vi,Z=Z)
rma1 <- rma.exact.fast(mu = attr(rma0,'mu'),yi=yi,vi=vi,Z=Z)
plot( rma0+rma1)
for(i in 1:3) {
    rma0 <- rma0+rma.exact.fast(mu = attr(rma0,'mu'),yi=yi,vi=vi,Z=Z)
    rma0 <- rma0+rma.exact.fast(tau2 = attr(rma0,'tau2'),yi=yi,vi=vi,Z=Z)
}
plot(rma0)



## 4. arbitrary grids
## plotting
source('rma.exact.R')
K <- 10
c0 <- 2
tau2 <- 12.5
vi <- (seq(1, 5, length=K))^2
mu0=0
yi=rnorm(K)*sqrt(vi+tau2)+mu0

Z <- matrix(rnorm(K*1e3),nrow=K)

alpha <- .05
mu <- NULL
tau2 <- NULL
grid.points <- NULL
resolution <- 1e2
rma0 <- rma.exact(grid.points=grid.points,yi=yi,vi=vi,Z=Z,resolution=2e1)
plot(rma0)


mu <- runif(1e3)*5-1
tau2 <- runif(1e3)*60
grid.points <- data.frame(mu=mu,tau2=tau2)
rma0 <- rma.exact(grid.points=grid.points,yi=yi,vi=vi,Z=Z)
plot(rma0)

## adding
source('rma.exact.R')
mu <- runif(1e3)*5-1
tau2 <- runif(1e3)*60
grid.points <- data.frame(mu=mu,tau2=tau2)
rma0 <- rma.exact(grid.points=grid.points,yi=yi,vi=vi,Z=Z)
mu <- runif(1e3)*5-1
tau2 <- runif(1e3)*60
grid.points <- data.frame(mu=mu,tau2=tau2)
rma1 <- rma.exact(grid.points=grid.points,yi=yi,vi=vi,Z=Z)
plot(rma0+rma1)

## 5. Lu's zooming suggestion
source('rma.exact.R')
K <- 10
c0 <- 2
tau2 <- 12.5
vi <- (seq(1, 5, length=K))^2
mu0=0
yi=rnorm(K)*sqrt(vi+tau2)+mu0

Z <- matrix(rnorm(K*1e3),nrow=K)

alpha <- .05
mu <- runif(1e3)*5-1
tau2 <- runif(1e3)*60
grid.points <- data.frame(mu=mu,tau2=tau2)
system.time({
rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=7e1)
plt <- plot(rma0)
plt
}) ## 1.9sec
## ggsave('041017b.png',plt)
system.time({
rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=4e1)
plt <- plot(rma0)
plt
}) ## 3.3sec
## ggsave('041017b.png',plt)
contour <- getContourLines(rma0$mu,rma0$tau2,rma0$p.val,levels=.05)
sd.tau2 <- sd(loess(y ~ x, data=contour)$resid)
sd.mu <- sd(loess(x ~ y, data=contour)$resid)
tau2.upper <- contour$y+sd.tau2
tau2.lower <- contour$y-sd.tau2
## png('041017c.png')
plot(contour$x,contour$y,col='red',type='l')
N <- 1e4
scale.factor <- .5
mu.noise <- sd.mu*rnorm(N)*scale.factor
tau2.noise <- sd.tau2*rnorm(N)*scale.factor
points(contour$x+mu.noise,contour$y+tau2.noise,col=rgb(0,0,0,.1))
 ## dev.off()
grid.points <- data.frame(mu=contour$x+mu.noise,tau2=contour$y+tau2.noise)
grid.points$tau2 <- pmax(0,grid.points$tau2)
rma1 <- rma.exact(grid.points=grid.points,yi=yi,vi=vi,Z=Z)
plot(rma0+rma1,levels=alpha)
## ggsave('041017d.png')
## plot(contours$x,contours$y)
## sd.tau2 <- sd(loess(y ~ x, data=contours)$resid)
## sd.mu <- sd(loess(x ~ y, data=contours)$resid)
## tau2.upper <- contours$y+sd.tau2
## tau2.lower <- contours$y-sd.tau2
## lines(contours$x,tau2.upper)

## lines(contours$x,loess0$fitted)
## plot(loess0$fitted,loess0$resid)
## N <- 1e3
## plot(contours$x,contours$y)
## lines(contours$x,contours$y)
## sd.y <- sd(loess(y ~ x, data=contours)$resid)
## sd.x <- sd(loess(x ~ y, data=contours)$resid)
## x.noise <- sd.x*rnorm(N)/5
## y.noise <- sd.y*rnorm(N)/5
## points(contours$x+x.noise,contours$y+y.noise)
## require(contoureR)


## 6. Lu's second zooming suggestion, grid around widest points

source('rma.exact.R')
K <- 10
c0 <- 2
tau2 <- 12.5
vi <- (seq(1, 5, length=K))^2
mu0=0
yi=rnorm(K)*sqrt(vi+tau2)+mu0

Z <- matrix(rnorm(K*1e3),nrow=K)

alpha <- .05
mu <- runif(1e3)*5-1
tau2 <- runif(1e3)*60
grid.points <- data.frame(mu=mu,tau2=tau2)
rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=3e1)
plt <- plot(rma0)
plt
rma1 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=3e1)
plot(rma0 +rma1)
for(i in 1:10) rma1 <- rma1 + rma.exact(yi=yi,vi=vi,Z=Z,resolution=3e1)
for(i in 1:5) rma1 <- rma1 + rma.exact(yi=yi,vi=vi,Z=matrix(rnorm(K*1e3),nrow=K),resolution=3e1)
plot(rma1)
plot(rma.exact(yi=yi,vi=vi,Z=Z,resolution=round(sqrt(5*30^2))))
with(rma.exact(yi=yi,vi=vi,Z=Z,resolution=round(sqrt(5*30^2))), plot(mu,tau2,pch='.',cex=2))

## using base graphics
cc <- contourLines(unique(rma0$mu),unique(rma0$tau2),matrix(rma0$p.val,nrow=length(unique(rma0$mu))),levels=.05)
lapply(cc,function(cc.i)plot(y ~ x, data=cc.i,type='l'))
sd.mu <- sd(rma0$mu)
sd.tau2 <- sd(rma0$tau2)
CI <- confint(rma0,alpha=.05)
scale <- 1.3
tau2 <- pmax(0,runif(5e1)*scale*sd.tau2 + (CI[,'tau2']-scale*sd.tau2/2)) %>% unique %>% sort
mu <- (runif(5e1)*scale*sd.mu + (CI[,'min.mu']-scale*sd.mu/2)) %>% sort
## mu <- seq(-6,-4,length.out=6e1)
## tau2 <- seq(0,60,length.out=6e1)
rma2 <- rma.exact(grid.points=expand.grid(mu=mu,tau2=tau2),yi=yi,vi=vi,B=1e3)
rma2.mtx <- matrix(rma2$p.val,nrow=length(unique(rma2$mu)))
new.lines <- contourLines(mu,tau2,rma2.mtx,levels=.05)
lapply(new.lines,function(cc)lines(cc$x,cc$y,col='red',lwd=3))
mu.upper <- (runif(5e1)*scale*sd.mu + (CI[,'max.mu']-scale*sd.mu/2)) %>% sort
rma.upper <- rma.exact(grid.points=expand.grid(mu=mu.upper,tau2=tau2),yi=yi,vi=vi,B=1e3)
upper.lines <- contourLines(mu.upper,tau2,matrix(rma.upper$p.val,nrow=length(unique(rma.upper$mu))),levels=.05)
invisible(lapply(upper.lines,function(cc)lines(cc$x,cc$y,col='red',lwd=3)))

## ggplot
cc <- contourLines(unique(rma0$mu),unique(rma0$tau2),matrix(rma0$p.val,nrow=length(unique(rma0$mu))),levels=.05)
lapply(cc,function(cc.i)plot(y ~ x, data=cc.i,type='l'))
plt <- plot(rma0)
plt
sd.mu <- sd(rma0$mu)
sd.tau2 <- sd(rma0$tau2)
CI <- confint(rma0,alpha=.05)
scale <- 1.3
tau2 <- pmax(0,runif(5e1)*scale*sd.tau2 + (CI[,'tau2']-scale*sd.tau2/2)) %>% unique %>% sort
mu <- (runif(5e1)*scale*sd.mu + (CI[,'min.mu']-scale*sd.mu/2)) %>% sort
## mu <- seq(-6,-4,length.out=6e1)
## tau2 <- seq(0,60,length.out=6e1)
rma.lower <- rma.exact(grid.points=expand.grid(mu=mu,tau2=tau2),yi=yi,vi=vi,B=1e3)
lower.lines <- contourLines(mu,tau2,matrix(rma.lower$p.val,nrow=length(mu)),levels=.05)
lower.lines <- mutate(as.data.frame(lower.lines),p.val=paste0(alpha,'new'),group='upper')
## plt <- plt+geom_path(data=lower.lines,aes(x=x,y=y,group=p.val,color=p.val))+scale_color_manual(values=c(cm.colors(6),'red'))
mu.upper <- (runif(5e1)*scale*sd.mu + (CI[,'max.mu']-scale*sd.mu/2)) %>% sort
rma.upper <- rma.exact(grid.points=expand.grid(mu=mu.upper,tau2=tau2),yi=yi,vi=vi,B=1e3)
upper.lines <- contourLines(mu.upper,tau2,matrix(rma.upper$p.val,nrow=length(unique(rma.upper$mu))),levels=.05)
upper.lines <- mutate(as.data.frame(upper.lines),p.val=paste0(alpha,'new'),group='lower')
new.lines <- rbind(upper.lines,lower.lines)
plt+geom_path(data=new.lines,aes(x=x,y=y,group=group,color=p.val))#+scale_color_manual(values=c(cm.colors(6),'red'))

## test zoom

source('rma.exact.R')
K <- 4
c0 <- 20
tau2 <- 12.5
vi <- (seq(1, 5, length=K))^2
mu0=0
yi=rnorm(K)*sqrt(vi+tau2)+mu0

Z <- matrix(rnorm(K*1e3),nrow=K)

alpha <- .05
mu <- runif(1e3)*5-1
tau2 <- runif(1e3)*60
grid.points <- data.frame(mu=mu,tau2=tau2)
rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=3e1,c0=c0)
rma0 <- zoom(rma0,alpha=.05,resolution=7e1,plot=TRUE)
confint(rma0)
plot(dd,levels=.05)
rma.uni(yi=yi,vi=vi,method='SJ')
## lines(cc2[[2]]$x,cc2[[2]]$y,col='red')

## plot(rma0+rma2)
## confint(rma0,alpha=.05)

## plot(rma0+rma1,levels=.05)
## cc <- contourLines(matrix(rma0$p.val,nrow=50),levels=.05)
## plot(cc[[1]],type='l')
## system.time({
## rma0 <- rma.exact(yi=yi,vi=vi,Z=Z,resolution=4e1)
## plt <- plot(rma0)
## plt
## }) ## 3.3sec
## ## ggsave('041017b.png',plt)
## contour <- getContourLines(rma0$mu,rma0$tau2,rma0$p.val,levels=.05)
## sd.tau2 <- sd(loess(y ~ x, data=contour)$resid)
## sd.mu <- sd(loess(x ~ y, data=contour)$resid)
## tau2.upper <- contour$y+sd.tau2
## tau2.lower <- contour$y-sd.tau2
## ## png('041017c.png')
## plot(contour$x,contour$y,col='red',type='l')
## N <- 1e4
## scale.factor <- .5
## mu.noise <- sd.mu*rnorm(N)*scale.factor
## tau2.noise <- sd.tau2*rnorm(N)*scale.factor
## points(contour$x+mu.noise,contour$y+tau2.noise,col=rgb(0,0,0,.1))
##  ## dev.off()
## grid.points <- data.frame(mu=contour$x+mu.noise,tau2=contour$y+tau2.noise)
## grid.points$tau2 <- pmax(0,grid.points$tau2)
## rma1 <- rma.exact(grid.points=grid.points,yi=yi,vi=vi,Z=Z)
## plot(rma0+rma1,levels=alpha)
## ## ggsave('041017d.png')
## ## plot(contours$x,contours$y)
## ## sd.tau2 <- sd(loess(y ~ x, data=contours)$resid)
## ## sd.mu <- sd(loess(x ~ y, data=contours)$resid)
## ## tau2.upper <- contours$y+sd.tau2
## ## tau2.lower <- contours$y-sd.tau2
## ## lines(contours$x,tau2.upper



## 7. compare scales of T1,T2
require(tensorA)

test.stats <- function(mu,tau2,yi,vi) {
    K <- length(yi)
    inv.vi <- as.tensor(1/vi,dims=c(K=K))
    sum.inv.vi <- sum(inv.vi)
    mu.fixed.obs <- yi%*%inv.vi/sum.inv.vi
    tau2.DL <- max(0, ((yi-mu.fixed.obs)^2%*%inv.vi - (K-1)) / (sum.inv.vi-sum(1/vi^2)/sum(1/vi)))
    mu.DL <- sum(yi/(vi+tau2.DL))/sum(1/(vi+tau2.DL))

    T1 <- (mu.DL-mu)^2 * sum(1/(tau2.DL+vi))
    T2 <- sum( (yi-mu)^2/(tau2+vi)+log(tau2+vi) - (yi-mu.DL)^2/(tau2.DL+vi) - log(tau2.DL+vi))/2
    return(c(T1=T1,T2=T2))
}

K <- 3
tau2 <- 2
vi <- .01*(seq(1, 5, length=K))^2
mu=0
reps <- 1e3
ans <- replicate(reps,{
    yi=rnorm(K)*sqrt(vi+tau2)+mu
    test.stats(mu=mu,tau2=tau2,yi=yi,vi=vi)
})
quantile(ans[2,]/ans[1,],probs=c(.2,.4,.6,.5,.8,.9,.95,.99))


## fix c0 vectorization
## find recommended c0, vary tau2 and K
## do  K=8,12,16,20 or larger simulation with DL, SJ, etc