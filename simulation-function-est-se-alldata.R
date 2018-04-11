###########################################################################
############     Privacy preserving analysis for survival data   ##########
rm(list = ls())
library(survival)

n = 5000        # lets fix some constants
beta1 = 1.5
beta2 = -2

### 20% prevalence, lT=.007; lC=.014 ### 90%, lT=.007; lC=.984 
lambdaT = .18  # baseline hazard   
lambdaC = .28  # hazard of censoring

r<-100   # the number of repititions
simdata<-list()
for(i in 1:r){
  x1 = rnorm(n,mean = 3,sd = .2)
  x2 = rnorm(n,mean = 3.5, sd = .6)
  # baseline hazard and censoring hazard are controlled above
  t = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1-beta2*x2)) # true event time
  C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
  time = pmin(t,C)  #observed time is min of censored and true time
  event<-ifelse(time==t,1,0)
  survdata<-data.frame(time=time, event=event, x1 = x1, x2 = x2)
  simdata[[i]]<-survdata 
  }

head(simdata[[1]]) # preview of the first dataset | 1000 datasets in the list
head(simdata[[r]]) # the last dataset

# check the prevalence of the event
(prevalence<-sum(simdata[[1]]$event)/length(simdata[[1]]$event)*100)
sum(simdata[[1]]$event)

#
simdata.artificial<-simdata

##############################################################################
### We need to modify the functions to produce the estimates and the SDs

####  Unpooled cox model

unpoolcox.SE<-function(r,simdata){
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  for(i in 1:r){
    survdata<-simdata[[i]]
    ress<-summary(coxph(Surv(time, event)~ x1 + x2, method="breslow", data = survdata))
    beta0vec[i]<-ress$coefficients[1,1]
    SEbeta0vec[i]<-ress$coefficients[1,3]
    beta1vec[i]<-ress$coefficients[2,1]
    SEbeta1vec[i]<-ress$coefficients[2,3]
  }
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec), median(beta1vec) ),
                     SE = c(median(SEbeta0vec), median(SEbeta1vec) ) )
  #result<-data.frame(beta0vec ,beta1vec)
  return(result)
}
#unpoolcox(r,simdata)


## equivalent reduced form: Conditional unpooled logistic regression
unpoolclogit.SE<-function(r,m,simdata){
  # m is the the number of controls in the  1:m pairing; the could vary accordingky
  
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    
    #  lets set the number of controls to m
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    
    # create a list of the event ids starting with the first event that took place
    # We have multiple events occuringt at the same time
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    ## we need to recover the eventids in the list of all the events
    eventid<-which(event==1)[tempfoo]
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-m matched set
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    
    ress<-summary(clogit(formula = event~x1+x2+strata(risksetid), data = survdata))
    beta0vec[z]<-ress$coefficients[1,1]
    SEbeta0vec[z]<-ress$coefficients[1,3]
    beta1vec[z]<-ress$coefficients[2,1]
    SEbeta1vec[z]<-ress$coefficients[2,3]
    
  }
  result<-data.frame(parameter=c("beta1", "beta2"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  #result<-data.frame(beta0vec ,beta1vec)
  return(result)
}
#unpoolclogit(r, simdata) 


############### All the multiple c-c controls settings  #########################

## model for pools of size 2
multi.cc.pool.SE<-function(s,r,m,simdata){  # this is for a 1:m case control pairing
  # s is the poolsize
  # r is the number of repititions
  # m is the the number of controls in the  1:m pairing
  
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  poolsize = s
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    
    #  lets set the number of controls to m
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    
    # create a list of the event ids starting with the first event that took place
    # We have multiple events occuringt at the same time
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    ## we need to recover the eventids in the list of all the events
    eventid<-which(event==1)[tempfoo]
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-m matched set
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    # lets determine the number of pools we could create | this is were the pooling starts
    poolcase.count<-sum(survdata$event, na.rm = T)
    
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count-1),2) # identifier for picking a new riskset to create a pool e.g. case 1and2 form pool 1, 3and4 form pool 2, etc
    
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|survdata$risksetid==(poolmarker[i]+1)]<-i
    }
    #head(survdata)
    ## We now pool case-control matched sets | recall how ids were assigned by time order,
    #  so close closer events are pooled together
    
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    
    ### lets try to transform the dataframe to a wide format eg reshape in R
    Rpoolsurvdata<-poolsurvdata[,2:6]
    cclen<-1+m
    wide.pooldata<-data.frame(event=rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1B = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2B= rep(NA,max(poolsurvdata$poolid)*cclen),
                              poolid= rep(NA,max(poolsurvdata$poolid)*cclen))
    
    Rpoolsurvdata<-Rpoolsurvdata[order(Rpoolsurvdata$risksetid),] #ensure the index are fixed
    
    oddmarker<-seq(1,max(Rpoolsurvdata$risksetid)-1,2)
    evenmarker<-seq(2,max(Rpoolsurvdata$risksetid),2)
    wide.pooldata[,c("event","x1A",
                     "x2A","poolid")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%oddmarker),
                                                     c("event","x1","x2","poolid")]
    wide.pooldata[,c("x1B","x2B")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%evenmarker),
                                                  c("x1","x2")]
    
    #head(wide.pooldata)
    wide.pooldata$sumx1<-wide.pooldata$x1A+wide.pooldata$x1B
    wide.pooldata$sumx2<-wide.pooldata$x2A+wide.pooldata$x2B
    
    #
    
    survdata<-wide.pooldata[,c("event","poolid","sumx1","sumx2")]
    survdata$id<-survdata$poolid
    
    ## conditional logistic for the pools
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
    
  }
  
  result<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  #result<-data.frame(beta0vec ,beta1vec)
  return(result)
}


## model for pools of size 4
multi.cc.pool4.SE<-function(s,r,m,simdata){  # this is for a 1:m case control pairing
  # s is the poolsize
  # r is the number of repititions
  # m is the the number of controls in the  1:m pairing
  
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  poolsize = s
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    
    #  lets set the number of controls to m
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    
    # create a list of the event ids starting with the first event that took place
    # We have multiple events occuringt at the same time
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    ## we need to recover the eventids in the list of all the events
    eventid<-which(event==1)[tempfoo]
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-m matched set
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    # lets determine the number of pools we could create | this is were the pooling starts
    poolcase.count<-sum(survdata$event, na.rm = T)
    
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count),s) # identifier for picking a new riskset to create a pool
    
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|
                        survdata$risksetid==(poolmarker[i]+1)|
                        survdata$risksetid==(poolmarker[i]+2)|
                        survdata$risksetid==(poolmarker[i]+3)]<-i
    }
    #head(survdata)
    ## We now pool case-control matched sets | recall how ids were assigned by time order,
    #  so close closer events are pooled together
    
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    ### lets try to transform the dataframe to a wide format eg reshape in R
    Rpoolsurvdata<-poolsurvdata[,2:6]
    cclen<-1+m
    wide.pooldata<-data.frame(event=rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1B = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1D = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2B= rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2D= rep(NA,max(poolsurvdata$poolid)*cclen),
                              id= rep(NA,max(poolsurvdata$poolid)*cclen),
                              test= rep(NA,max(poolsurvdata$poolid)*cclen))
    Rpoolsurvdata<-Rpoolsurvdata[order(Rpoolsurvdata$risksetid),]
    #Rpoolsurvdata<-Rpoolsurvdata[,1:5]
    
    marker1<-seq(1,max(Rpoolsurvdata$risksetid),4)
    marker2<-seq(2,max(Rpoolsurvdata$risksetid),4)
    marker3<-seq(3,max(Rpoolsurvdata$risksetid),4)
    marker4<-seq(4,max(Rpoolsurvdata$risksetid),4)
    
    wide.pooldata[,c("event","x1A",
                     "x2A","id")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker1),
                                                 c("event","x1","x2","poolid")]
    wide.pooldata[,c("x1B","x2B")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker2),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1C","x2C")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker3),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1D","x2D")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker4),
                                                  c("x1","x2")]
    
    #head(wide.pooldata)
    wide.pooldata$sumx1<-wide.pooldata$x1A+wide.pooldata$x1B+wide.pooldata$x1C+wide.pooldata$x1D
    wide.pooldata$sumx2<-wide.pooldata$x2A+wide.pooldata$x2B+wide.pooldata$x2C+wide.pooldata$x2D
    
    #
    survdata<-wide.pooldata[,c("event","id","sumx1","sumx2")]
    #survdata$id<-survdata$poolid
    
    ## conditional logistic for the pools
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
    
  }
  result<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  #result<-data.frame(beta0vec ,beta1vec)
  return(result)
}



## model for pools of size 6
multi.cc.pool6.SE<-function(s,r,m,simdata){  # this is for a 1:m case control pairing
  # s is the poolsize == g in this case
  # r is the number of repititions
  # m is the the number of controls in the  1:m pairing
  
  beta0vec<-NA
  SEbeta0vec<-NA
  beta1vec<- NA
  SEbeta1vec<-NA
  poolsize = s
  r = r
  for(z in 1:r){
    survdata<-simdata[[z]]
    survdata$risksetid<-NA
    eventid<-NA
    time=survdata$time
    event=survdata$event
    rowlength<-sum(event)
    
    #  lets set the number of controls to m
    resurvdata<-NA
    resurvdata<-data.frame(time = c(survdata$time,rep(NA,rowlength*m)),
                           event = c(survdata$event,rep(NA,rowlength*m)),
                           x1 = c(survdata$x1,rep(NA,rowlength*m)),
                           x2 = c(survdata$x2,rep(NA,rowlength*m)),
                           risksetid = c(survdata$risksetid,rep(NA,rowlength*m)) )
    
    # create a list of the event ids starting with the first event that took place
    # We have multiple events occuringt at the same time
    lfoo<-sort.int(time[event==1], decreasing = FALSE, index.return = TRUE)
    tempfoo<-lfoo$ix    # a list of the event index
    ## we need to recover the eventids in the list of all the events
    eventid<-which(event==1)[tempfoo]
    
    foo<-c()       # an emptly list to keep the ids that are matched on progressively
    l1riskset<-NA  # list for riskset excluding future events
    
    # A loop starting with the first event that is taking place
    # create the pooled ids for all the events
    # allocate a control for each case/event untill all the cases are used up
    
    for(i in 1:length(eventid)){
      l1<-which(time>time[eventid[i]]) # all the eligible controls
      
      if(length(l1)>=m ){
        resurvdata$risksetid[eventid[i]]<-i    # risksetid for the case
        index<-sample(x = l1,size = m, replace = T)  # randomly select one of the controls to match on 
        survdata$risksetid[index]<-i         # risksetid for the control in a 1-m case-control setting
        temp<-(which(is.na(resurvdata$time))[1]):(which(is.na(resurvdata$time))[1]+m-1)
        resurvdata[temp,]<-survdata[index,]
        resurvdata$event[temp]<-0        # the selected person could be a future case. make sure it is a control
        
      }
    }
    #head(resurvdata)  # A dataset of 1-m matched set
    survdata<-resurvdata[!is.na(resurvdata$risksetid),]
    
    # lets determine the number of pools we could create | this is were the pooling starts
    poolcase.count<-sum(survdata$event, na.rm = T)
    
    case.count<-floor(poolcase.count/poolsize)*poolsize #we need even case-control count to pool
    survdata$poolid<-NA
    poolmarker<-seq(1,(case.count),s) # identifier for picking a new riskset to create a pool
    
    for(i in 1:length(poolmarker)){ # poolids for risksets; 1,2==id 1, 3,4= id 2, etc
      survdata$poolid[survdata$risksetid==poolmarker[i]|
                        survdata$risksetid==(poolmarker[i]+1)|
                        survdata$risksetid==(poolmarker[i]+2)|
                        survdata$risksetid==(poolmarker[i]+3)|
                        survdata$risksetid==(poolmarker[i]+4)|
                        survdata$risksetid==(poolmarker[i]+5)]<-i
    }
    #head(survdata)
    ## We now pool case-control matched sets | recall how ids were assigned by time order,
    #  so close closer events are pooled together
    
    poolsurvdata<-survdata[!is.na(survdata$poolid),]   
    survdata=NA
    ### lets try to transform the dataframe to a wide format eg reshape in R
    Rpoolsurvdata<-poolsurvdata[,2:6]
    cclen<-1+m
    wide.pooldata<-data.frame(event=rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1B = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1D = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1E = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x1F = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2A = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2B= rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2C = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2D= rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2E = rep(NA,max(poolsurvdata$poolid)*cclen),
                              x2F= rep(NA,max(poolsurvdata$poolid)*cclen),
                              id= rep(NA,max(poolsurvdata$poolid)*cclen))
    Rpoolsurvdata<-Rpoolsurvdata[order(Rpoolsurvdata$risksetid),]
    #Rpoolsurvdata<-Rpoolsurvdata[,1:5]
    
    marker1<-seq(1,max(Rpoolsurvdata$risksetid),s)
    marker2<-seq(2,max(Rpoolsurvdata$risksetid),s)
    marker3<-seq(3,max(Rpoolsurvdata$risksetid),s)
    marker4<-seq(4,max(Rpoolsurvdata$risksetid),s)
    marker5<-seq(5,max(Rpoolsurvdata$risksetid),s)
    marker6<-seq(6,max(Rpoolsurvdata$risksetid),s)
    
    wide.pooldata[,c("event","x1A",
                     "x2A","id")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker1),
                                                 c("event","x1","x2","poolid")]
    wide.pooldata[,c("x1B","x2B")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker2),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1C","x2C")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker3),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1D","x2D")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker4),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1E","x2E")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker5),
                                                  c("x1","x2")]
    wide.pooldata[,c("x1F","x2F")]<-Rpoolsurvdata[which(Rpoolsurvdata$risksetid%in%marker6),
                                                  c("x1","x2")]
    
    #head(wide.pooldata)
    wide.pooldata$sumx1<-wide.pooldata$x1A+wide.pooldata$x1B+wide.pooldata$x1C+
      wide.pooldata$x1D+wide.pooldata$x1E+wide.pooldata$x1F
    wide.pooldata$sumx2<-wide.pooldata$x2A+wide.pooldata$x2B+wide.pooldata$x2C+
      wide.pooldata$x2D+wide.pooldata$x2E+wide.pooldata$x2F
    
    #
    survdata<-wide.pooldata[,c("event","id","sumx1","sumx2")]
    #survdata$id<-survdata$poolid
    
    ## conditional logistic for the pools
    ssp2<-summary(clogit(formula = event~sumx1+sumx2+strata(id), data = survdata))
    beta0vec[z]<-ssp2$coefficients[1,1]
    SEbeta0vec[z]<-ssp2$coefficients[1,3]
    beta1vec[z]<-ssp2$coefficients[2,1]
    SEbeta1vec[z]<-ssp2$coefficients[2,3]
    
  }
  result<-data.frame(parameter=c("beta0", "beta1"),
                     estimate=c(median(beta0vec) ,median(beta1vec) ),
                     SE = c(median(SEbeta0vec),median(SEbeta1vec) )  )
  #result<-data.frame(beta0vec ,beta1vec)
  return(result)
}



############### 2) Simulation using the SMART dataset ########################
##############################################################################
require(hdnom)
smartdata<-(smarto)
head(smartdata)

# The biomakers to consider in this analysis 
# Markers of atherosclerosis (disease)

# HOMOC - Homocyste; 463 missing values
# GLUT - Glutamine; 19 missing values
# CREAT - Creatinine clearance, in mL/min; 17 missing values
# ALBUMIN - Albumin (no, micro, macro); 28 missing values
# IMT - Intima media thickness, in mm; 98 missing values
# STENOSIS - Carotid artery stenosis > 50%; 93 missing values

# first, lets select a subset that is useful for our analysis
subsmart<- smartdata[,c("TEVENT","EVENT","SEX","AGE","DIABETES",
                        "STENOSIS","HOMOC","GLUT","CREAT","IMT")]

#### data prep: same format as the input format
data1<-subsmart[,c("TEVENT","EVENT","CREAT","IMT")]
names(data1)<-c("time","event","x1","x2")
dataf<-list()
dataf[[1]]<-data1[!is.na(data1$x1&data1$x2),]
length(dataf[1])

simdata.smart<-dataf
r<-length(dataf[1])    # the number of datasets
#sum(dataf[[1]]$event)  # this is the number of events


###########################################################

simdata<-simdata.artificial
summary(simdata[[1]])
r = 100
m = 5
cox<-unpoolcox.SE(r, simdata)
unpool<-unpoolclogit.SE(r,m,simdata)
model1<-multi.cc.pool.SE(s=2,r,m,simdata)
model2<-multi.cc.pool4.SE(s=4,r,m,simdata)
model3<-multi.cc.pool6.SE(s=6,r,m,simdata)

### use the functions in the file "all sim functions" taylored for
### the SMART dataset


