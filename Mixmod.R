library(mixtools)
library(ggplot2)
library(gridExtra)

#source("R/llh_covars.R") ### (-) log-likelihood function (try to write it!!)
llh<-function(params,data,covars){
  loglik<-0
  for( i in 1:length(data)){
    w<-exp(params[1]+params[2]*covars[i,1])/(1+exp(params[1]+params[2]*covars[i,1]))
    q<-exp(params[11])/(1+exp(params[11]))
    mu<-params[3]+params[4]*covars[i,1]+params[5]*covars[i,2]+params[6]*covars[i,3]+params[7]*covars[i,4]+params[8]*covars[i,5]
    +params[9]*covars[i,6]
    sigma<-params[10]
    loglik<-loglik+log((1-w)*dnorm(data[i],mu/q,sigma/q)+w*dnorm(data[i],mu,sigma))
  }
  
  return (-1*loglik)
  
}

### Load the data
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
print(prova)
### Initial values for covariates from linear regression model
linmod <- lm(pr3$incid~covars[, 1]+covars[, 2]+covars[, 3]+covars[, 4]+covars[, 5]+covars[, 6])

### Estimates, standard errors and confidence intervals for w and q
max.llh <- nlm(f=llh, p=c(log(prova$lambda[which(prova$beta[1,]==min(prova$beta[1,]))]/(1-prova$lambda[which(prova$beta[1,]==min(prova$beta[1,]))])), 
                          -0.5, linmod$coefficients, prova$sigma[which(prova$beta[1,]==min(prova$beta[1,]))],
                          log((prova$beta[1,which(prova$beta[1,]==min(prova$beta[1,]))]/
                                 prova$beta[1,which(prova$beta[1,]==max(prova$beta[1,]))])/(1-prova$beta[1,which(prova$beta[1,]==min(prova$beta[1,]))]/
                                                                                              prova$beta[1,which(prova$beta[1,]==max(prova$beta[1,]))]))),
               data=pr3$incid, covars=covars, hessian=TRUE)

q<-exp(max.llh$estimate[11])/(1+exp(max.llh$estimate[11]))
varcovar<-solve(max.llh$hessian)
sds<-sqrt(diag(varcovar))
# confidence interval for alpha0
max.llh$estimate[1]-qnorm(0.975)*sds[1]
max.llh$estimate[1]+qnorm(0.975)*sds[1]

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
bacf   <- acf(resid, lag.max = 10, plot = TRUE)
bacfdf <- with(bacf[1:10], data.frame(lag, acf))
conf.level <- 0.95
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(y_agg$x))
q1 <- ggplot(data = bacfdf, mapping = aes(x = as.integer(lag), y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = ciline), linetype=2) + geom_hline(aes(yintercept = -ciline), linetype=2)+
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ylab("") + xlab("Lag") + ggtitle("ACF") + theme(plot.title = element_text(hjust = 0.5))
bacf   <- pacf(resid, lag.max = 10, plot = TRUE)
bacfdf <- with(bacf, data.frame(lag, acf))
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(y_agg$x))
q2 <- ggplot(data = bacfdf, mapping = aes(x = as.integer(lag), y = acf)) +
  geom_hline(aes(yintercept = 0)) + geom_hline(aes(yintercept = ciline), linetype=2) + geom_hline(aes(yintercept = -ciline), linetype=2)+
  geom_segment(mapping = aes(xend = lag, yend = 0)) + ylab("") + xlab("Lag") + ggtitle("PACF") + theme(plot.title = element_text(hjust = 0.5))
postscript("/Users/vishalnair/Downloads/Modelling/Mixture models/fig1.eps", width=10, height=5)
grid.arrange(q1, q2, ncol=2)
dev.off()

### Residuals qqplot
df <- data.frame(resid)
p <- ggplot(df, aes(sample=resid))
postscript("/Users/vishalnair/Downloads/Modelling/Mixture models/fig2.eps", width=10, height=5)
p+stat_qq()+stat_qq_line()+xlab("Teoretical quantiles")+ylab("Sample quantiles")
dev.off()

### Reconstruction of the hidden processes
# a. Females, 15-29 years old
gw.women1 <- pr3[pr3$sexe==0 & pr3$age==0, ]
post<-matrix(nrow=dim(gw.women1)[1],ncol=2)
w<-vector()
### Calculation of the posterior probabilities
for(i in 1:dim(gw.women1)[1]){
  w[i]<-exp(2.9743256 -4.2817967*i)/(1+exp(2.9743256-4.2817967*i))
  post[i,1]<-w[i]*dnorm(gw.women1$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                          max.llh$estimate[5]*covars[i,2]+
                          max.llh$estimate[6]*covars[i,3]+
                          max.llh$estimate[7]*covars[i,4]+
                          max.llh$estimate[8]*covars[i,5]+
                          max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])/
    (w[i]*dnorm(gw.women1$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
     max.llh$estimate[5]*covars[i,2]+
     max.llh$estimate[6]*covars[i,3]+
     max.llh$estimate[7]*covars[i,4]+
     max.llh$estimate[8]*covars[i,5]+
     max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])+(1-w[i])*dnorm(gw.women1$incid[i],mean=(max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                                                                       max.llh$estimate[5]*covars[i,2]+
                                                                       max.llh$estimate[6]*covars[i,3]+
                                                                       max.llh$estimate[7]*covars[i,4]+
                                                                       max.llh$estimate[8]*covars[i,5]+
                                                                       max.llh$estimate[9]*covars[i,6])/q,sd=max.llh$estimate[10]/q))
  post[i,2]<- 1-post[i,1] 
}
# b. Females, over 30 years old
gw.women2 <- pr3[pr3$sexe==0 & pr3$age==1, ]
### Calculation of the posterior probabilities
post2<-matrix(nrow=dim(gw.women2)[1],ncol=2)
w2<-vector()
### Calculation of the posterior probabilities
for(i in 1:dim(gw.women2)[1]){
  w2[i]<-exp(2.9743256 -4.2817967*i)/(1+exp(2.9743256-4.2817967*i))
  post[i,1]<-w2[i]*dnorm(gw.women2$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                          max.llh$estimate[5]*covars[i,2]+
                          max.llh$estimate[6]*covars[i,3]+
                          max.llh$estimate[7]*covars[i,4]+
                          max.llh$estimate[8]*covars[i,5]+
                          max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])/
    (w2[i]*dnorm(gw.women2$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                  max.llh$estimate[5]*covars[i,2]+
                  max.llh$estimate[6]*covars[i,3]+
                  max.llh$estimate[7]*covars[i,4]+
                  max.llh$estimate[8]*covars[i,5]+
                  max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])+(1-w2[i])*dnorm(gw.women2$incid[i],mean=(max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                                                                                                                     max.llh$estimate[5]*covars[i,2]+
                                                                                                                     max.llh$estimate[6]*covars[i,3]+
                                                                                                                     max.llh$estimate[7]*covars[i,4]+
                                                                                                                     max.llh$estimate[8]*covars[i,5]+
                                                                                                                     max.llh$estimate[9]*covars[i,6])/q,sd=max.llh$estimate[10]/q))
  post2[i,2]<- 1-post2[i,1] 
}
# c. Males, 15-29 years old
gw.men1 <- pr3[pr3$sexe==1 & pr3$age==0, ]
### Calculation of the posterior probabilities
post3<-matrix(nrow=dim(gw.men1)[1],ncol=2)
w3<-vector()
### Calculation of the posterior probabilities
for(i in 1:dim(gw.men1)[1]){
  w3[i]<-exp(2.9743256 -4.2817967*i)/(1+exp(2.9743256-4.2817967*i))
  post3[i,1]<-w3[i]*dnorm(gw.men1$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                          max.llh$estimate[5]*covars[i,2]+
                          max.llh$estimate[6]*covars[i,3]+
                          max.llh$estimate[7]*covars[i,4]+
                          max.llh$estimate[8]*covars[i,5]+
                          max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])/
    (w3[i]*dnorm(gw.men1$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                  max.llh$estimate[5]*covars[i,2]+
                  max.llh$estimate[6]*covars[i,3]+
                  max.llh$estimate[7]*covars[i,4]+
                  max.llh$estimate[8]*covars[i,5]+
                  max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])+(1-w3[i])*dnorm(gw.men1$incid[i],mean=(max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                                                                                                                     max.llh$estimate[5]*covars[i,2]+
                                                                                                                     max.llh$estimate[6]*covars[i,3]+
                                                                                                                     max.llh$estimate[7]*covars[i,4]+
                                                                                                                     max.llh$estimate[8]*covars[i,5]+
                                                                                                                     max.llh$estimate[9]*covars[i,6])/q,sd=max.llh$estimate[10]/q))
  post3[i,2]<- 1-post3[i,1] 
}
# d. Males, 30-94 years old
gw.men2 <- pr3[pr3$sexe==1 & pr3$age==1, ]
### Calculation of the posterior probabilities
post4<-matrix(nrow=dim(gw.men2)[1],ncol=2)
w4<-vector()
### Calculation of the posterior probabilities
for(i in 1:dim(gw.men2)[1]){
  w4[i]<-exp(2.9743256 -4.2817967*i)/(1+exp(2.9743256-4.2817967*i))
  post[i,1]<-w4[i]*dnorm(gw.men2$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                          max.llh$estimate[5]*covars[i,2]+
                          max.llh$estimate[6]*covars[i,3]+
                          max.llh$estimate[7]*covars[i,4]+
                          max.llh$estimate[8]*covars[i,5]+
                          max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])/
    (w4[i]*dnorm(gw.men2$incid[i],mean=max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                  max.llh$estimate[5]*covars[i,2]+
                  max.llh$estimate[6]*covars[i,3]+
                  max.llh$estimate[7]*covars[i,4]+
                  max.llh$estimate[8]*covars[i,5]+
                  max.llh$estimate[9]*covars[i,6],sd=max.llh$estimate[10])+(1-w4[i])*dnorm(gw.men2$incid[i],mean=(max.llh$estimate[3]+max.llh$estimate[4]*covars[i,1]+
                                                                                                                     max.llh$estimate[5]*covars[i,2]+
                                                                                                                     max.llh$estimate[6]*covars[i,3]+
                                                                                                                     max.llh$estimate[7]*covars[i,4]+
                                                                                                                     max.llh$estimate[8]*covars[i,5]+
                                                                                                                     max.llh$estimate[9]*covars[i,6])/q,sd=max.llh$estimate[10]/q))
  post4[i,2]<- 1-post4[i,1] 
}