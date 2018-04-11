#################################################################
#### For 1:5 case-control setting
# Using the artifical dataset

simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta1","Beta2")
ests.mat.p0<-unpoolclogit(r,m=5, simdata)
names(ests.mat.p0)<-c("Beta1","Beta2")
ests.mat.p2<-multi.cc.pool(s=2,r,m=5,simdata)
names(ests.mat.p2)<-c("Beta1","Beta2")
ests.mat.p4<-multi.cc.pool4(s=4,r,m=5,simdata)
names(ests.mat.p4)<-c("Beta1","Beta2")
ests.mat.p6<-multi.cc.pool6(s=6,r,m=5,simdata)
names(ests.mat.p6)<-c("Beta1","Beta2")

# compiling the beta 1 estimates
set.seed(2334)
par(mar=c(5,5,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta1`,ests.mat.p0$`Beta1`,
                       ests.mat.p2$`Beta1`,ests.mat.p4$`Beta1`,
                       ests.mat.p6$`Beta1`
)
boxplot(est.mat.b0,main='',ylab ="Hazard Ratio",
        xlab ="Type of model",
        names = c('Cox','P0-1:5M','P2-1:5M','P4-1:5M','P6-1:5M'),
        #col = "gray",
        ylim=c(-.2,3.6),las=1);
abline(h=c(1.5),lty=2,col='red');

## The beta 2 estimates for the different pool sizes
set.seed(2556)
par(mar=c(5,5,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta2`,ests.mat.p0$`Beta2`,
                       ests.mat.p2$`Beta2`,ests.mat.p4$`Beta2`,
                       ests.mat.p6$`Beta2`
)
boxplot(est.mat.b1,
        #main='Hazard Ratio of X2 in 1-5 matched sets',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c('Cox','P0-1:5M','P2-1:5M','P4-1:5M','P6-1:5M'),
        #col = "gray",
        ylim=c(-3.1,-1.1),las=1);
abline(h=c(-2.01),lty=2,col='red')

######### SMART data   ######################
## bootstrap via repititions

simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  

r=100

cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m=10,simdata)
p2mod<-multi.cc.pool.SE(s=2,r,m=10,simdata)
p4mod<-multi.cc.pool4.SE(s=4,r,m=10,simdata)
p6mod<-multi.cc.pool6.SE(s=6,r,m=10,simdata)
p6mod10<-multi.cc.pool6.SE(s=6,r,m=20,simdata)


beta1est<-c(cox$estimate[1],unpool$estimate[1],p2mod$estimate[1],p4mod$estimate[1],
            p6mod$estimate[1],
            p6mod10$estimate[1])
beta2est<-c(cox$estimate[2],unpool$estimate[2],p2mod$estimate[2],p4mod$estimate[2],
            p6mod$estimate[2],
            p6mod10$estimate[2])
SEbeta1est<-c(cox$SE[1],unpool$SE[1],p2mod$SE[1],p4mod$SE[1],
            p6mod$SE[1],
            p6mod10$SE[1])
SEbeta2est<-c(cox$SE[2],unpool$SE[2],p2mod$SE[2],p4mod$SE[2],
            p6mod$SE[2],
            p6mod10$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:6))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:6))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox","P0-1:10M","P2-1:10M","P4-1:10M","P6-1:10M","P6-1:20M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.005), xlim = c(0.5,7),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0002),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.6,1.4),xlim = c(0.5,6.6),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.04),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)


### Need all the possible combination of the result
##











