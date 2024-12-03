library(mixtools)
library(ggplot2)
library(gridExtra)

### Read the data
load("data_main.RData")
### Covariates matrix
t         <- rep(seq(1, 96), 4)/96
pr3$inter <- as.numeric(pr3$age)*as.numeric(pr3$sexe)
sint      <- sin(2*pi*t/3)
cost      <- cos(2*pi*t/3)
covars    <- cbind(t, pr3$age, pr3$sexe, pr3$inter, sint, cost)

### Direct estimation via mixtools (used as initial values for MLE)
w0     <- 0.7
q0     <- 0.5
prova  <- regmixEM(pr3$incid, covars, lambda=c(w0, (1-w0)), 
                   beta=matrix(c(mean(pr3$incid), 0, 0, 0, 0, 0, 0, mean(pr3$incid)/q0, 0, 0, 0, 0, 0 ,0), ncol=2, nrow=ncol(covars)+1), 
                   sigma=c(sd(pr3$incid), sd(pr3$incid)/q0),
                   k=2, addintercept=TRUE, epsilon=1e-12, maxit=10000)

### Initial values for covariates from linear regression model
linmod <- lm(pr3$incid~covars[, 1]+covars[, 2]+covars[, 3]+covars[, 4]+covars[, 5]+covars[, 6])

### Estimates, standard errors and confidence intervals for w and q
max.llh <- nlm(f=llh, p=c(log(prova$lambda[which(prova$beta[1,]==min(prova$beta[1,]))]/(1-prova$lambda[which(prova$beta[1,]==min(prova$beta[1,]))])), 
                          -0.5, linmod$coefficients, prova$sigma[which(prova$beta[1,]==min(prova$beta[1,]))],
                          log((prova$beta[1,which(prova$beta[1,]==min(prova$beta[1,]))]/
                                 prova$beta[1,which(prova$beta[1,]==max(prova$beta[1,]))])/(1-prova$beta[1,which(prova$beta[1,]==min(prova$beta[1,]))]/
                                                                                              prova$beta[1,which(prova$beta[1,]==max(prova$beta[1,]))]))),
               data=pr3$incid, covars=covars, hessian=TRUE)

q <- exp(max.llh$estimate[11])/(1+exp(max.llh$estimate[11]))

sigma        <- solve(max.llh$hessian)
lim.inf95_1  <- max.llh$estimate[1] - qnorm(0.975)*sqrt(diag(sigma))[1] #intercept w
lim.sup95_1  <- max.llh$estimate[1] + qnorm(0.975)*sqrt(diag(sigma))[1] #intercept w
lim.inf95_2  <- max.llh$estimate[2] - qnorm(0.975)*sqrt(diag(sigma))[2] #time w
lim.sup95_2  <- max.llh$estimate[2] + qnorm(0.975)*sqrt(diag(sigma))[2] #time w
lim.inf95_3  <- max.llh$estimate[3] - qnorm(0.975)*sqrt(diag(sigma))[3] #intercept mean
lim.sup95_3  <- max.llh$estimate[3] + qnorm(0.975)*sqrt(diag(sigma))[3] #intercept mean
lim.inf95_4  <- max.llh$estimate[4] - qnorm(0.975)*sqrt(diag(sigma))[4] #t
lim.sup95_4  <- max.llh$estimate[4] + qnorm(0.975)*sqrt(diag(sigma))[4] #t
lim.inf95_5  <- max.llh$estimate[5] - qnorm(0.975)*sqrt(diag(sigma))[5] #age
lim.sup95_5  <- max.llh$estimate[5] + qnorm(0.975)*sqrt(diag(sigma))[5] #age
lim.inf95_6  <- max.llh$estimate[6] - qnorm(0.975)*sqrt(diag(sigma))[6] #sex
lim.sup95_6  <- max.llh$estimate[6] + qnorm(0.975)*sqrt(diag(sigma))[6] #sex
lim.inf95_7  <- max.llh$estimate[7] - qnorm(0.975)*sqrt(diag(sigma))[7] #interaction age*sex
lim.sup95_7  <- max.llh$estimate[7] + qnorm(0.975)*sqrt(diag(sigma))[7] #interaction age*sex
lim.inf95_8  <- max.llh$estimate[8] - qnorm(0.975)*sqrt(diag(sigma))[8] #sin
lim.sup95_8  <- max.llh$estimate[8] + qnorm(0.975)*sqrt(diag(sigma))[8] #sin
lim.inf95_9  <- max.llh$estimate[9] - qnorm(0.975)*sqrt(diag(sigma))[9] #cos
lim.sup95_9  <- max.llh$estimate[9] + qnorm(0.975)*sqrt(diag(sigma))[9] #cos
lim.inf95_10 <- max.llh$estimate[10] - qnorm(0.975)*sqrt(diag(sigma))[10] #sd
lim.sup95_10 <- max.llh$estimate[10] + qnorm(0.975)*sqrt(diag(sigma))[10] #sd
lim.inf95_11 <- max.llh$estimate[11] - qnorm(0.975)*sqrt(diag(sigma))[11] #logit(q)
lim.sup95_11 <- max.llh$estimate[11] + qnorm(0.975)*sqrt(diag(sigma))[11] #logit(q)

### Confidence interval for q
exp(lim.inf95_11)/(1+exp(lim.inf95_11)); exp(lim.sup95_11)/(1+exp(lim.sup95_11))

### Reconstruction of the hidden processes
# a. Females, 15-29 years old
gw.women1 <- pr3[pr3$sexe==0 & pr3$age==0, ]
### Calculation of the posterior probabilities
post <- matrix(nrow=96, ncol=2)
w    <- vector()
for (i in 1:96)
{
  w[i] <- exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96)/(1+exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96))
  post[i, 1] <- w[i]*dnorm(gw.women1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[8]*covars[i, 5]+
                                                       max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])/
    (w[i]*dnorm(gw.women1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[8]*covars[i, 5]+
                                            max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.women1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[8]*covars[i, 5]+
                                                  max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
  post[i, 2] <- (1-w[i])*dnorm(gw.women1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[8]*covars[i, 5]+
                                                           max.llh$estimate[9]*covars[i, 6])/q, sd=max.llh$estimate[10]/q)/
    (w[i]*dnorm(gw.women1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[8]*covars[i, 5]+
                                            max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.women1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[8]*covars[i, 5]+
                                                  max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
}
xrec1 <- ifelse(post[, 2] > 0.5, gw.women1$incid, gw.women1$incid/q)

mean(gw.women1$incid); mean(xrec1)
(mean(xrec1)-mean(gw.women1$incid))/mean(gw.women1$incid)*100

par(mfrow=c(2, 2))
gw.dones1.ts <- ts(gw.women1$incid, start=c(2009, 1), end=c(2016, 12), freq=12)
ts.plot(gw.dones1.ts, ylim=c(9, 32), ylab="Incidence x 100,000", main="Women 15-29 years old")
lines(seq(2009, 2016.99, 1/12), xrec1, col="red", lty=2)

# b. Females, over 30 years old
gw.women2 <- pr3[pr3$sexe==0 & pr3$age==1, ]
### Calculation of the posterior probabilities
post <- matrix(nrow=96, ncol=2)
w    <- vector()
for (i in 1:96)
{
  w[i] <- exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96)/(1+exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96))
  post[i, 1] <- w[i]*dnorm(gw.women2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+
                                                       max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])/
    (w[i]*dnorm(gw.women2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[8]*covars[i, 5]+
                                            max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.women2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[8]*covars[i, 5]+
                                                  max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
  post[i, 2] <- (1-w[i])*dnorm(gw.women2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[8]*covars[i, 5]+
                                                           max.llh$estimate[9]*covars[i, 6])/q, sd=max.llh$estimate[10]/q)/
    (w[i]*dnorm(gw.women2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[8]*covars[i, 5]+
                                            max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.women2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[8]*covars[i, 5]+
                                                  max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
}
xrec2 <- ifelse(post[,2] > 0.5, gw.women2$incid, gw.women2$incid/q)

mean(gw.women2$incid); mean(xrec2)
(mean(xrec2)-mean(gw.women2$incid))/mean(gw.women2$incid)*100

gw.dones2.ts <- ts(gw.women2$incid, start=c(2009, 1), end=c(2016, 12), freq=12)
ts.plot(gw.dones2.ts, ylim=c(1, 8), ylab="Incidence x 100,000", main="Women 30-94 years old")
lines(seq(2009, 2016.99, 1/12), xrec2, col="red", lty=2)

# c. Males, 15-29 years old
gw.men1 <- pr3[pr3$sexe==1 & pr3$age==0, ]
### Calculation of the posterior probabilities
post <- matrix(nrow=96, ncol=2)
w    <- vector()
for (i in 1:96)
{
  w[i] <- exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96)/(1+exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96))
  post[i, 1] <- w[i]*dnorm(gw.men1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[6]+
                                                     max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])/
    (w[i]*dnorm(gw.men1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[6]+max.llh$estimate[8]*covars[i, 5]+
                                          max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.men1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[6]+max.llh$estimate[8]*covars[i, 5]+
                                                max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
  post[i, 2] <- (1-w[i])*dnorm(gw.men1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[6]+max.llh$estimate[8]*covars[i, 5]+
                                                         max.llh$estimate[9]*covars[i, 6])/q, sd=max.llh$estimate[10]/q)/
    (w[i]*dnorm(gw.men1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[6]+max.llh$estimate[8]*covars[i, 5]+
                                          max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.men1$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[6]+max.llh$estimate[8]*covars[i, 5]+
                                                max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
}
xrec3 <- ifelse(post[,2] > 0.5, gw.men1$incid, gw.men1$incid/q)

mean(gw.men1$incid); mean(xrec3)
(mean(xrec3)-mean(gw.men1$incid))/mean(gw.men1$incid)*100

gw.homes1.ts <- ts(gw.men1$incid, start=c(2009, 1), end=c(2016, 12), freq=12)
ts.plot(gw.homes1.ts, ylim=c(4, 32), ylab="Incidence x 100,000", main="Men 15-29 years old")
lines(seq(2009, 2016.99, 1/12), xrec3, col="red", lty=2)

# d. Males, 30-94 years old
gw.men2 <- pr3[pr3$sexe==1 & pr3$age==1, ]
### Calculation of the posterior probabilities
post <- matrix(nrow=96, ncol=2)
w    <- vector()
for (i in 1:96)
{
  w[i] <- exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96)/(1+exp(max.llh$estimate[1]+max.llh$estimate[2]*i/96))
  post[i, 1] <- w[i]*dnorm(gw.men2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[6]+
                                                     max.llh$estimate[7]+max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])/
    (w[i]*dnorm(gw.men2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[6]+
                                          max.llh$estimate[7]+max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.men2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[6]+
                                                max.llh$estimate[7]+max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
  post[i, 2] <- (1-w[i])*dnorm(gw.men2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[6]+
                                                         max.llh$estimate[7]+max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6])/q, sd=max.llh$estimate[10]/q)/
    (w[i]*dnorm(gw.men2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[6]+
                                          max.llh$estimate[7]+max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6]), sd=max.llh$estimate[10])+
       (1-w[i])*dnorm(gw.men2$incid[i], mean=(max.llh$estimate[3]+max.llh$estimate[4]*i/96+max.llh$estimate[5]+max.llh$estimate[6]+
                                                max.llh$estimate[7]+max.llh$estimate[8]*covars[i, 5]+max.llh$estimate[9]*covars[i, 6])/q, 
                      sd=max.llh$estimate[10]/q))
}
xrec4 <- ifelse(post[,2] > 0.5, gw.men2$incid, gw.men2$incid/q)

mean(gw.men2$incid); mean(xrec4)
(mean(xrec4)-mean(gw.men2$incid))/mean(gw.men2$incid)*100

gw.homes2.ts <- ts(gw.men2$incid, start=c(2009, 1), end=c(2016, 12), freq=12)
ts.plot(gw.homes2.ts, ylim=c(1, 12), ylab="Incidence x 100,000", main="Men 30-94 years old")
lines(seq(2009, 2016.99, 1/12), xrec4, col="red", lty=2)


postscript("/Users/dmorina/Desktop/fig3.eps", width=12, height=8)
par(mfrow=c(2, 2))
ts.plot(gw.dones1.ts, ylim=c(9, 32), ylab="Incidence x 100,000", main="Women 15-29 years old")
lines(seq(2009, 2016.99, 1/12), xrec1, col="red", lty=2)
ts.plot(gw.dones2.ts, ylim=c(1, 8), ylab="Incidence x 100,000", main="Women 30-94 years old")
lines(seq(2009, 2016.99, 1/12), xrec2, col="red", lty=2)
ts.plot(gw.homes1.ts, ylim=c(4, 32), ylab="Incidence x 100,000", main="Men 15-29 years old")
lines(seq(2009, 2016.99, 1/12), xrec3, col="red", lty=2)
ts.plot(gw.homes2.ts, ylim=c(1, 12), ylab="Incidence x 100,000", main="Men 30-94 years old")
lines(seq(2009, 2016.99, 1/12), xrec4, col="red", lty=2)
dev.off()

### Global validation (residuals analysis)
y_est <- vector()
w     <- vector()
for (i in 1:384)
{
  j <- (i %% 96)/96
  if (j == 0) j <- 1
  w[i] <- exp(max.llh$estimate[1]+max.llh$estimate[2]*j)/(1+exp(max.llh$estimate[1]+max.llh$estimate[2]*j))
  m    <- max.llh$estimate[3]+max.llh$estimate[4]*j+max.llh$estimate[5]*pr3$age[i]+
    max.llh$estimate[6]*pr3$sexe[i]+max.llh$estimate[7]*pr3$inter[i]+max.llh$estimate[8]*covars[i, 5]+
    max.llh$estimate[9]*covars[i, 6]
  y_est[i] <- w[i]*m+(1-w[i])*m/q
}

y_est_agg_temp <- data.frame(t=rep(seq(1:96), 4), sexe=c(rep(0, 192), rep(1, 192)), 
                             edat=c(rep(0, 96), rep(1, 96), rep(0, 96), rep(1, 96)), y_est)
tw             <- sum(unique(pr3$Pob)) 
y_est_agg      <- aggregate(y_est_agg_temp$y_est, by=list(y_est_agg_temp$t), FUN=sum)
y_agg          <- aggregate(pr3$incid, by=list(pr3$mes_any_problema), FUN=sum)
y_est_agg$x[1:12]  <- y_est_agg$x[1:12] *sum(unique(pr3$Pob[pr3$Year==2009]))/tw
y_est_agg$x[13:24] <- y_est_agg$x[13:24]*sum(unique(pr3$Pob[pr3$Year==2010]))/tw
y_est_agg$x[25:36] <- y_est_agg$x[25:36]*sum(unique(pr3$Pob[pr3$Year==2011]))/tw
y_est_agg$x[37:48] <- y_est_agg$x[37:48]*sum(unique(pr3$Pob[pr3$Year==2012]))/tw
y_est_agg$x[49:60] <- y_est_agg$x[49:60]*sum(unique(pr3$Pob[pr3$Year==2013]))/tw
y_est_agg$x[61:72] <- y_est_agg$x[61:72]*sum(unique(pr3$Pob[pr3$Year==2014]))/tw
y_est_agg$x[73:84] <- y_est_agg$x[73:84]*sum(unique(pr3$Pob[pr3$Year==2015]))/tw
y_est_agg$x[85:96] <- y_est_agg$x[85:96]*sum(unique(pr3$Pob[pr3$Year==2016]))/tw

y_agg$x[1:12]  <- y_agg$x[1:12] *sum(unique(pr3$Pob[pr3$Year==2009]))/tw
y_agg$x[13:24] <- y_agg$x[13:24]*sum(unique(pr3$Pob[pr3$Year==2010]))/tw
y_agg$x[25:36] <- y_agg$x[25:36]*sum(unique(pr3$Pob[pr3$Year==2011]))/tw
y_agg$x[37:48] <- y_agg$x[37:48]*sum(unique(pr3$Pob[pr3$Year==2012]))/tw
y_agg$x[49:60] <- y_agg$x[49:60]*sum(unique(pr3$Pob[pr3$Year==2013]))/tw
y_agg$x[61:72] <- y_agg$x[61:72]*sum(unique(pr3$Pob[pr3$Year==2014]))/tw
y_agg$x[73:84] <- y_agg$x[73:84]*sum(unique(pr3$Pob[pr3$Year==2015]))/tw
y_agg$x[85:96] <- y_agg$x[85:96]*sum(unique(pr3$Pob[pr3$Year==2016]))/tw

### Residuals
resid <- y_est_agg$x-y_agg$x

### ACF and PACF
bacf   <- acf(resid, lag.max = 10, plot = FALSE)
bacfdf <- with(bacf[1:10], data.frame(lag, acf))
conf.level <- 0.95
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(y_agg$x))
q1 <- ggplot(data = bacfdf, mapping = aes(x = as.integer(lag), y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = ciline), linetype=2) + geom_hline(aes(yintercept = -ciline), linetype=2)+
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ylab("") + xlab("Lag") + ggtitle("ACF") + theme(plot.title = element_text(hjust = 0.5))
bacf   <- pacf(resid, lag.max = 10, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(y_agg$x))
q2 <- ggplot(data = bacfdf, mapping = aes(x = as.integer(lag), y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = ciline), linetype=2) + geom_hline(aes(yintercept = -ciline), linetype=2)+
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ylab("") + xlab("Lag") + ggtitle("PACF") + theme(plot.title = element_text(hjust = 0.5))
postscript("/Users/dmorina/Desktop/fig1.eps", width=10, height=5)
grid.arrange(q1, q2, ncol=2)
dev.off()

### Residuals qqplot
df <- data.frame(resid)
p <- ggplot(df, aes(sample=resid))
postscript("/Users/dmorina/Desktop/fig2.eps", width=10, height=5)
p+stat_qq()+stat_qq_line()+xlab("Teoretical quantiles")+ylab("Sample quantiles")
dev.off()

