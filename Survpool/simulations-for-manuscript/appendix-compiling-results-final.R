###############################################################
### (1)       Simulations for pools of size s/g=2
###############################################################
###                                                       ###
###                 Artificial dataset                    ###
#############################################################
simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
s = 2

ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta0-cox","Beta1-cox")
ests.mat.p0<-unpoolclogit(r,m=1, simdata)
names(ests.mat.p0)<-c("Beta0-P0","Beta1-P0")
ests.mat.p2<-multi.cc.pool(s,r,m=1,simdata)
names(ests.mat.p2)<-c("Beta0-P2","Beta1-P2")
ests.mat.m2<-multi.cc.pool(s,r,m=2,simdata)
names(ests.mat.m2)<-c("Beta0-m2","Beta1-m2")
ests.mat.m4<-multi.cc.pool(s,r,m=4,simdata)
names(ests.mat.m4)<-c("Beta0-m4","Beta1-m4")
ests.mat.m6<-multi.cc.pool(s,r,m=6,simdata)
names(ests.mat.m6)<-c("Beta0-m6","Beta1-m6")

# compiling the beta 0 estimates
set.seed(23)
par(mar=c(5,5,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta0-cox`,ests.mat.p0$`Beta0-P0`,
                       ests.mat.p2$`Beta0-P2`,ests.mat.m2$`Beta0-m2`,
                       ests.mat.m4$`Beta0-m4`,ests.mat.m6$`Beta0-m6`)
boxplot(est.mat.b0,main='Hazard Ratio of X by model',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c("Cox","Unpooled","P2-1:1","P2-1:2","P2-1:4","P2-1:6"),
        col = "gray",ylim=c(-.5,3.5),las=1);
abline(h=c(1.5),lty=2,col='red');

## The beta 1 estimates for the different pool sizes
set.seed(25)
par(mar=c(5,5,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta1-cox`,ests.mat.p0$`Beta1-P0`,
                       ests.mat.p2$`Beta1-P2`,ests.mat.m2$`Beta1-m2`,ests.mat.m4$`Beta1-m4`,
                       ests.mat.m6$`Beta1-m6`)
boxplot(est.mat.b1,main='Hazard Ratio of Z by model',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c("Cox","Unpooled","P2-1:1","P2-1:2","P2-1:4","P2-1:6"),
        col = "gray",ylim=c(-3.0,-1.2),las=1);
abline(h=c(-2),lty=2,col='red')
## the results seem to be consistent 


###                      SMART dataset    poolsize g=2    ###
#############################################################
simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]])) 
dim(simdata[[1]])  

s = 2    # poolsize
r = r

ucox<-unpoolcox.SE(r, simdata)
#uclogit<-unpoolclogit.SE(r,m=1,simdata)
p2clogit<-multi.cc.pool.SE(s,r,m=1,simdata)
m212<-multi.cc.pool.SE(s,r,m=2,simdata)
m214<-multi.cc.pool.SE(s,r,m=5,simdata)
m216<-multi.cc.pool.SE(s,r,m=10,simdata)
m2110<-multi.cc.pool.SE(s,r,m=20,simdata)


beta1est<-c(ucox$estimate[1],
            #uclogit$estimate[1],
            p2clogit$estimate[1],m212$estimate[1],
            m214$estimate[1],m216$estimate[1],m2110$estimate[1])
beta2est<-c(ucox$estimate[2],
            #uclogit$estimate[2],
            p2clogit$estimate[2],m212$estimate[2],
            m214$estimate[2],m216$estimate[2],m2110$estimate[2])
SEbeta1est<-c(ucox$SE[1],
             # uclogit$SE[1],
              p2clogit$SE[1],m212$SE[1],
              m214$SE[1],m216$SE[1],m2110$SE[1])
SEbeta2est<-c(ucox$SE[2],
              #uclogit$SE[2],
              p2clogit$SE[2],m212$SE[2],
              m214$SE[2],m216$SE[2],m2110$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:6))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:6))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox",
        #"P0-1:1M",
        "P2-1:1M","P2-1:2M","P2-1:5M","P2-1:10M","P2-1:20M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.007), xlim = c(0.5,6.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0008),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.5,1.7),xlim = c(0.5,6.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)
  
##########################  END    ###########################


#############################################################
### (2)       Simulations for pools of size s/g=4
#############################################################
###                                                       ###
###                 Artificial dataset                    ###
#############################################################
simdata<-simdata.artificial
summary(simdata[[1]])
r = 100  ## 100 repititions specified above
s = 4  # the poolsize

ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta0-cox","Beta1-cox")
ests.mat.p0<-unpoolclogit(r,m=1, simdata)
names(ests.mat.p0)<-c("Beta0-P0","Beta1-P0")
ests.mat.p2<-multi.cc.pool4(s,r,m=1,simdata)
names(ests.mat.p2)<-c("Beta0-P2","Beta1-P2")
ests.mat.m2<-multi.cc.pool4(s,r,m=2,simdata)
names(ests.mat.m2)<-c("Beta0-m2","Beta1-m2")
ests.mat.m4<-multi.cc.pool4(s,r,m=4,simdata)
names(ests.mat.m4)<-c("Beta0-m4","Beta1-m4")
ests.mat.m6<-multi.cc.pool4(s,r,m=6,simdata)
names(ests.mat.m6)<-c("Beta0-m6","Beta1-m6")

# compiling the beta estimates
set.seed(2323)
par(mar=c(5,5,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta0-cox`,ests.mat.p0$`Beta0-P0`,
                       ests.mat.p2$`Beta0-P2`,ests.mat.m2$`Beta0-m2`,
                       ests.mat.m4$`Beta0-m4`,ests.mat.m6$`Beta0-m6`)
boxplot(est.mat.b0,main='Hazard Ratio of X by model',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c("Cox","Unpooled","P4-1:1","P4-1:2","P4-1:4","P4-1:6"),
        col = "gray",ylim=c(-.5,4.5),las=1);
abline(h=c(1.5),lty=2,col='red');

## The gamma estimates for the different pool sizes
set.seed(25)
par(mar=c(5,5,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta1-cox`,ests.mat.p0$`Beta1-P0`,
                       ests.mat.p2$`Beta1-P2`,ests.mat.m2$`Beta1-m2`,ests.mat.m4$`Beta1-m4`,
                       ests.mat.m6$`Beta1-m6`)
boxplot(est.mat.b1,main='Hazard Ratio of Z by model',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c("Cox","Unpooled","P4-1:1","P4-1:2","P4-1:4","P4-1:6"),
        col = "gray",ylim=c(-3.2,-0.9),las=1);
abline(h=c(-2),lty=2,col='red')
## the results seem to be consistent and in general unbiased

############################################################


###                                                       ###
###                      SMART dataset                    ###
#############################################################
simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  

####
# lets compile the results
r = r
s = 4    # poolsize g=4

cox<-unpoolcox.SE(r, simdata)
#unpool<-unpoolclogit.SE(r,m=1,simdata)
model1<-multi.cc.pool.SE(s,r,m=1,simdata)
model2<-multi.cc.pool.SE(s,r,m=2,simdata)
model3<-multi.cc.pool.SE(s,r,m=5,simdata)
model4<-multi.cc.pool.SE(s,r,m=10,simdata)
model5<-multi.cc.pool.SE(s,r,m=20,simdata)


beta1est<-c(cox$estimate[1],
            #unpool$estimate[1],
            model1$estimate[1],model2$estimate[1],
            model3$estimate[1],model4$estimate[1],model5$estimate[1])
beta2est<-c(cox$estimate[2],
            #unpool$estimate[2],
            model1$estimate[2],model2$estimate[2],
            model3$estimate[2],model4$estimate[2],model5$estimate[2])
SEbeta1est<-c(cox$SE[1],
              #unpool$SE[1],
              model1$SE[1],model2$SE[1],
              model3$SE[1],model4$SE[1],model5$SE[1])
SEbeta2est<-c(cox$SE[2],
              #unpool$SE[2],
              model1$SE[2],model2$SE[2],
              model3$SE[2],model4$SE[2],model5$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:6))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:6))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox",
        #"P0-1:1M",
        "P4-1:1M","P4-1:2M","P4-1:5M","P4-1:10M","P4-1:20M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.007), xlim = c(0.5,6.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0008),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of gamma
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.4,1.8),xlim = c(0.5,6.6),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)



#############################################################
### (3)       Simulations for poolsize s/g=6
#############################################################
###                                                       ###
###                 Artificial dataset                    ###
#############################################################
simdata<-simdata.artificial
summary(simdata[[1]])
r = r  ## 100 repititions specified above
s = 6  # the poolsize

ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta0-cox","Beta1-cox")
ests.mat.p0<-unpoolclogit(r,m=1, simdata)  # matching is set at m=1 but this can vary
names(ests.mat.p0)<-c("Beta0-P0","Beta1-P0")
ests.mat.p2<-multi.cc.pool6(s,r,m=1,simdata)
names(ests.mat.p2)<-c("Beta0-P2","Beta1-P2")
ests.mat.m2<-multi.cc.pool6(s,r,m=2,simdata)
names(ests.mat.m2)<-c("Beta0-m2","Beta1-m2")
ests.mat.m4<-multi.cc.pool6(s,r,m=4,simdata)
names(ests.mat.m4)<-c("Beta0-m4","Beta1-m4")
ests.mat.m6<-multi.cc.pool6(s,r,m=6,simdata)
names(ests.mat.m6)<-c("Beta0-m6","Beta1-m6")

# compiling the beta estimates
set.seed(23)
par(mar=c(5,5,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta0-cox`,ests.mat.p0$`Beta0-P0`,
                       #ests.mat.p2$`Beta0-P2`,ests.mat.m2$`Beta0-m2`,
                       ests.mat.m4$`Beta0-m4`,ests.mat.m6$`Beta0-m6`)

boxplot(est.mat.b0,main='Hazard Ratio of X by model',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c("Cox","Unpooled",
                                         #"P6-1:1","P6-1:2",
                                         "P6-1:4","P6-1:6"),
        col = "gray",ylim=c(-1,5),las=1);
abline(h=c(1.5),lty=2,col='red');

## The beta 1 estimates for the different pool sizes
set.seed(2535)
par(mar=c(5,5,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta1-cox`,ests.mat.p0$`Beta1-P0`,
                       #  ests.mat.p2$`Beta1-P2`,ests.mat.m2$`Beta1-m2`,
                       ests.mat.m4$`Beta1-m4`,
                       ests.mat.m6$`Beta1-m6`)

boxplot(est.mat.b1,main='Hazard Ratio of Z by model',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c("Cox","Unpooled",
                  #"P6-1:1","P6-1:2",
                  "P6-1:4","P6-1:6"),
        col = "gray",ylim=c(-3.5,-1),las=1);
abline(h=c(-2),lty=2,col='red')

############################################################


###                                                       ###
###                      SMART dataset                    ###
#############################################################
simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  
#
# lets compile the results
s = 6    # poolsize g=6
r = r
cox<-unpoolcox.SE(r, simdata)
#unpool<-unpoolclogit.SE(r,m=1,simdata)
model1<-multi.cc.pool.SE(s,r,m=1,simdata)
model2<-multi.cc.pool.SE(s,r,m=2,simdata)
model3<-multi.cc.pool.SE(s,r,m=5,simdata)
model4<-multi.cc.pool.SE(s,r,m=10,simdata)
model5<-multi.cc.pool.SE(s,r,m=30,simdata)


beta1est<-c(cox$estimate[1],
            #unpool$estimate[1],
            model1$estimate[1],model2$estimate[1],
            model3$estimate[1],model4$estimate[1],model5$estimate[1])
beta2est<-c(cox$estimate[2],
            #unpool$estimate[2],
            model1$estimate[2],model2$estimate[2],
            model3$estimate[2],model4$estimate[2],model5$estimate[2])
SEbeta1est<-c(cox$SE[1],
              #unpool$SE[1],
              model1$SE[1],model2$SE[1],
              model3$SE[1],model4$SE[1],model5$SE[1])
SEbeta2est<-c(cox$SE[2],
              #unpool$SE[2],
              model1$SE[2],model2$SE[2],
              model3$SE[2],model4$SE[2],model5$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:6))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:6))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox",
        #"P0-1:1M",
        "P6-1:1M","P6-1:2M","P6-1:5M","P6-1:10M","P6-1:20M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.007), xlim = c(0.5,6.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0008),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.4,1.7),xlim = c(0.5,6.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)



#############################################################
##            Compiling the results by pool sizes
#############################################################
###                   1:1 matched sets

#### For a fixed case-control setting i.e. 1:m case-control setting
# Using the artifical dataset

simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta1","Beta2")
ests.mat.p0<-unpoolclogit(r,m=1, simdata)
names(ests.mat.p0)<-c("Beta1","Beta2")
ests.mat.p2<-multi.cc.pool(s=2,r,m=1,simdata)
names(ests.mat.p2)<-c("Beta1","Beta2")
ests.mat.p4<-multi.cc.pool4(s=4,r,m=1,simdata)
names(ests.mat.p4)<-c("Beta1","Beta2")
ests.mat.p6<-multi.cc.pool6(s=6,r,m=1,simdata)
names(ests.mat.p6)<-c("Beta1","Beta2")

# compiling the beta estimates
set.seed(2334)
par(mar=c(5,5,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta1`,ests.mat.p0$`Beta1`,
                       ests.mat.p2$`Beta1`,ests.mat.p4$`Beta1`
                       #ests.mat.p6$`Beta1`
                       )
boxplot(est.mat.b0,main='HR of X in 1-1 matched control per case',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c('Cox','Pool0','Pool2','Pool4'
                                        # 'Pool6'
                                         ),
        col = "gray",ylim=c(-.5,4.0),las=1);
abline(h=c(1.54),lty=2,col='red');

## The beta 1 estimates for the different pool sizes
set.seed(2556)
par(mar=c(5,5,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta2`,ests.mat.p0$`Beta2`,
                       ests.mat.p2$`Beta2`,ests.mat.p4$`Beta2`
                      # ests.mat.p6$`Beta2`
                       )
boxplot(est.mat.b1,main='HR of Z in 1-1 matched control per case',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c('Cox','Pool0','Pool2','Pool4'
                 # 'Pool6'
                  ),
        col = "gray",ylim=c(-3.2,-0.9),las=1);
abline(h=c(-2),lty=2,col='red')

######### lets try the same for the real  ######################

simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  

m =1 
r =r

cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m,simdata)
model1<-multi.cc.pool.SE(s=2,r,m,simdata)
model2<-multi.cc.pool4.SE(s=4,r,m,simdata)
model3<-multi.cc.pool6.SE(s=6,r,m,simdata)


beta1est<-c(cox$estimate[1],
            #unpool$estimate[1],
            model1$estimate[1],model2$estimate[1],
            model3$estimate[1])
beta2est<-c(cox$estimate[2],
            #unpool$estimate[2],
            model1$estimate[2],model2$estimate[2],
            model3$estimate[2])
SEbeta1est<-c(cox$SE[1],
              #unpool$SE[1],
              model1$SE[1],model2$SE[1],
              model3$SE[1])
SEbeta2est<-c(cox$SE[2],
              #uclogit$SE[2],
              model1$SE[2],model2$SE[2],
              model3$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:4))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:4))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox",
        #"P0-1:1M",
        "P2-1:1M","P4-1:1M","P6-1:1M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.007), xlim = c(0.5,4.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0008),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.4,1.7),xlim = c(0.5,4.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)

#################### sim data 1:2 #############################
#### For 1:2 case-control setting
# Using the artifical dataset

simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta1","Beta2")
ests.mat.p0<-unpoolclogit(r,m=2, simdata)
names(ests.mat.p0)<-c("Beta1","Beta2")
ests.mat.p2<-multi.cc.pool(s=2,r,m=2,simdata)
names(ests.mat.p2)<-c("Beta1","Beta2")
ests.mat.p4<-multi.cc.pool4(s=4,r,m=2,simdata)
names(ests.mat.p4)<-c("Beta1","Beta2")
ests.mat.p6<-multi.cc.pool6(s=6,r,m=2,simdata)
names(ests.mat.p6)<-c("Beta1","Beta2")

# compiling the beta 0 estimates
set.seed(2334)
par(mar=c(4,4,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta1`,ests.mat.p0$`Beta1`,
                       ests.mat.p2$`Beta1`,ests.mat.p4$`Beta1`
                       #ests.mat.p6$`Beta1`
)
boxplot(est.mat.b0,main='Hazard Ratio of X in 1-2 matched sets',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c('Cox','Pool0','Pool2','Pool4'
                                          #'Pool6'
        ),
        col = "gray",ylim=c(-.5,4.5),las=1);
abline(h=c(1.55),lty=2,col='red');

## The beta 1 estimates for the different pool sizes
set.seed(2556)
par(mar=c(4,4,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta2`,ests.mat.p0$`Beta2`,
                       ests.mat.p2$`Beta2`,ests.mat.p4$`Beta2`
                      # ests.mat.p6$`Beta2`
)
boxplot(est.mat.b1,main='Hazard Ratio of Z in 1-2 matched sets',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c('Cox','Pool0','Pool2','Pool4'
                 #  'Pool6'
        ),
        col = "gray",ylim=c(-3.9,-0.9),las=1);
abline(h=c(-2),lty=2,col='red')

######### lets try the same for the real data   ######################

simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  

r=r
m=2

cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m,simdata)
model1<-multi.cc.pool.SE(s=2,r,m,simdata)
model2<-multi.cc.pool4.SE(s=4,r,m,simdata)
model3<-multi.cc.pool6.SE(s=6,r,m,simdata)


beta1est<-c(cox$estimate[1],unpool$estimate[1],model1$estimate[1],model2$estimate[1],
            model3$estimate[1])
beta2est<-c(cox$estimate[2],unpool$estimate[2],model1$estimate[2],model2$estimate[2],
            model3$estimate[2])
SEbeta1est<-c(cox$SE[1],unpool$SE[1],model1$SE[1],model2$SE[1],
              model3$SE[1])
SEbeta2est<-c(cox$SE[2],uclogit$SE[2],model1$SE[2],model2$SE[2],
              model3$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:5))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:5))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox","P0-1:2M","P2-1:2M","P4-1:2M","P6-1:2M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.008), xlim = c(0.5,5.8),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0008),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.4,1.8),xlim = c(0.5,5.6),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)

#####
######################sim data##############################
#### For 1:4 case-control setting
# Using the artifical dataset

simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta1","Beta2")
ests.mat.p0<-unpoolclogit(r,m=4, simdata)
names(ests.mat.p0)<-c("Beta1","Beta2")
ests.mat.p2<-multi.cc.pool(s=2,r,m=4,simdata)
names(ests.mat.p2)<-c("Beta1","Beta2")
ests.mat.p4<-multi.cc.pool4(s=4,r,m=4,simdata)
names(ests.mat.p4)<-c("Beta1","Beta2")
ests.mat.p6<-multi.cc.pool6(s=6,r,m=4,simdata)
names(ests.mat.p6)<-c("Beta1","Beta2")

# compiling the beta 0 estimates
set.seed(2334)
par(mar=c(4,4,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta1`,ests.mat.p0$`Beta1`,
                       ests.mat.p2$`Beta1`,ests.mat.p4$`Beta1`
                       #ests.mat.p6$`Beta1`
)
boxplot(est.mat.b0,main='Hazard Ratio of X in 1-4 matched sets',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c('Cox','Pool0','Pool2','Pool4'
                                         #'Pool6'
        ),
        col = "gray",ylim=c(-.2,4.0),las=1);
abline(h=c(1.55),lty=2,col='red');

## The beta 1 estimates for the different pool sizes
set.seed(2556)
par(mar=c(4,4,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta2`,ests.mat.p0$`Beta2`,
                       ests.mat.p2$`Beta2`,ests.mat.p4$`Beta2`
                       # ests.mat.p6$`Beta2`
)
boxplot(est.mat.b1,main='Hazard Ratio of Z in 1-4 matched sets',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c('Cox','Pool0','Pool2','Pool4'
                  #  'Pool6'
        ),
        col = "gray",ylim=c(-3.1,-1.1),las=1);
abline(h=c(-2),lty=2,col='red')

######### lets try the same for the real data   ######################

simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  


m=4
r=r
cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m,simdata)
model1<-multi.cc.pool.SE(s=2,r,m,simdata)
model2<-multi.cc.pool4.SE(s=4,r,m,simdata)
model3<-multi.cc.pool6.SE(s=6,r,m,simdata)


beta1est<-c(cox$estimate[1],unpool$estimate[1],model1$estimate[1],model2$estimate[1],
            model3$estimate[1])
beta2est<-c(cox$estimate[2],unpool$estimate[2],model1$estimate[2],model2$estimate[2],
            model3$estimate[2])
SEbeta1est<-c(cox$SE[1],unpool$SE[1],model1$SE[1],model2$SE[1],
              model3$SE[1])
SEbeta2est<-c(cox$SE[2],uclogit$SE[2],model1$SE[2],model2$SE[2],
              model3$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:5))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:5))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox","P0-1:4M","P2-1:4M","P4-1:4M","P6-1:4M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.008), xlim = c(0.5,5.8),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0008),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.4,1.8),xlim = c(0.5,5.6),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[2],lty=2,col=2)

#################################################################
#### For 1:10 case-control setting
# Using the artifical dataset

simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
m=10
ests.mat.uncox<-unpoolcox(r, simdata)
names(ests.mat.uncox)<-c("Beta1","Beta2")
ests.mat.p0<-unpoolclogit(r,m=10, simdata)
names(ests.mat.p0)<-c("Beta1","Beta2")
ests.mat.p2<-multi.cc.pool(s=2,r,m=10,simdata)
names(ests.mat.p2)<-c("Beta1","Beta2")
ests.mat.p4<-multi.cc.pool4(s=4,r,m=10,simdata)
names(ests.mat.p4)<-c("Beta1","Beta2")
ests.mat.p6<-multi.cc.pool6(s=6,r,m=10,simdata)
names(ests.mat.p6)<-c("Beta1","Beta2")

# compiling the beta 0 estimates
set.seed(2334)
par(mar=c(5,5,3,2))
est.mat.b0<-data.frame(ests.mat.uncox$`Beta1`,ests.mat.p0$`Beta1`,
                       ests.mat.p2$`Beta1`,ests.mat.p4$`Beta1`,
                       ests.mat.p6$`Beta1`
)
boxplot(est.mat.b0,main='HR of X in 1-10 matched controls per case',ylab ="Hazard Ratio",
        xlab ="Type of model", names = c('Cox','Pool0','Pool2','Pool4',
                                         'Pool6'
        ),
        col = "gray",ylim=c(-.0,3.2),las=1);
abline(h=c(1.5),lty=2,col='red');

## The beta 1 estimates for the different pool sizes
set.seed(2556)
par(mar=c(5,5,3,2))
est.mat.b1<-data.frame(ests.mat.uncox$`Beta2`,ests.mat.p0$`Beta2`,
                       ests.mat.p2$`Beta2`,ests.mat.p4$`Beta2`,
                       ests.mat.p6$`Beta2`
)
boxplot(est.mat.b1,main='Hazard Ratio of Z in 1-10 matched sets',
        ylab ="Hazard Ratio", xlab ="Type of model",
        names = c('Cox','Pool0','Pool2','Pool4',
                  'Pool6'
        ),
        col = "gray",ylim=c(-2.6,-1.5),las=1);
abline(h=c(-2),lty=2,col='red')

######### lets try the same for the real data   ######################

simdata<-simdata.smart
head(simdata[[1]])
tail(head(simdata[[1]]))
dim(simdata[[1]])  


r=r  # the number of bootstraps
m=10
cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m,simdata)
model1<-multi.cc.pool.SE(s=2,r,m,simdata)
model2<-multi.cc.pool4.SE(s=4,r,m,simdata)
model3<-multi.cc.pool6.SE(s=6,r,m,simdata)


beta1est<-c(cox$estimate[1],unpool$estimate[1],model1$estimate[1],model2$estimate[1],
            model3$estimate[1])
beta2est<-c(cox$estimate[2],unpool$estimate[2],model1$estimate[2],model2$estimate[2],
            model3$estimate[2])
SEbeta1est<-c(cox$SE[1],unpool$SE[1],model1$SE[1],model2$SE[1],
              model3$SE[1])
SEbeta2est<-c(cox$SE[2],uclogit$SE[2],model1$SE[2],model2$SE[2],
              model3$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:5))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:5))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox","P0-1:10M","P2-1:10M","P4-1:10M","P6-1:10M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.001,0.005), xlim = c(0.5,5.8),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0004),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.5,1.55),xlim = c(0.5,5.6),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)


########################################
### lets try 20 controls per case matching for the different pools
r=r  # the number of bootstraps
m=25
cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m,simdata)
model1<-multi.cc.pool.SE(s=2,r,m,simdata)
model2<-multi.cc.pool4.SE(s=4,r,m,simdata)
model3<-multi.cc.pool6.SE(s=6,r,m,simdata)


beta1est<-c(cox$estimate[1],unpool$estimate[1],model1$estimate[1],model2$estimate[1],
            model3$estimate[1])
beta2est<-c(cox$estimate[2],unpool$estimate[2],model1$estimate[2],model2$estimate[2],
            model3$estimate[2])
SEbeta1est<-c(cox$SE[1],unpool$SE[1],model1$SE[1],model2$SE[1],
              model3$SE[1])
SEbeta2est<-c(cox$SE[2],uclogit$SE[2],model1$SE[2],model2$SE[2],
              model3$SE[2])

dfbeta1<-data.frame(b1=beta1est, b1SE=SEbeta1est, x=rep(1:5))
dfbeta2<-data.frame(b2=beta2est, b2SE=SEbeta2est, x=rep(1:5))

# 
### using the plot in R base
set.seed(3433)
par(mar=c(5,5,3,3))
Text<-c("Cox","P0-1:20M","P2-1:20M","P4-1:20M","P6-1:20M")
plot(dfbeta1$x,dfbeta1$b1,
     ylim = c(0.0015,0.004), xlim = c(0.5,5.5),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of Creatinine by model'
)
arrows(dfbeta1$x, dfbeta1$b1-dfbeta1$b1SE, dfbeta1$x, dfbeta1$b1+dfbeta1$b1SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta1$x,y = dfbeta1$b1-(dfbeta1$b1SE+0.0004),labels = Text, cex = .8)
abline(h = dfbeta1$b1[1],col=2,lty=2)

# plot of lthe estimate of beta2
set.seed(3433)
par(mar=c(5,5,3,3))
plot(dfbeta2$x,dfbeta2$b2,
     ylim = c(.5,1.55),xlim = c(0.5,5.6),
     pch=19, xlab="Type of model", ylab='Hazard Ratio',
     main='Hazard Ratio of IMT by model')
arrows(dfbeta2$x, dfbeta2$b2-dfbeta2$b2SE, dfbeta2$x, dfbeta2$b2+dfbeta2$b2SE,
       length=0.05, angle=90, code=3)
text(x = dfbeta2$x,y = dfbeta2$b2-(dfbeta2$b2SE+0.08),labels = Text, cex = .8)
abline(h = dfbeta2$b2[1],lty=2,col=2)
