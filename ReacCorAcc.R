#### Reactivation Correlated with Recognition Accuracy Stats ####
# assumes all required files are in the working directory
#
# instructions:
# - run code in "Dependances" (download any missing packages)
# - run code in "Get Data"
# - run code in "Reactivation Residuals"
# - run code in "Z-score Variables"
# - the data in Figures 2-4 can be recreated by running the code below with the appropriate headings 
#   - the models are provided under the main figure headings (e.g. "Figure 2: ...")
#   - sub-figure code (e.g. "Figure 2a: ...") uses precalculated model results
#   - model/bootstrap results are precalculated so it is not necessary to run the models (can take >3 hours)



#### Dependances ####

library(lme4)
library(ggplot2)



#### Get Data ####
# load the desired reactivation (rank format) results:
#
# 'ReacRankRecall.RData' contains the recall reactivation data used in the paper.
# trials are along the rows (sorted alphabetically by cue label, not in temporal order).
# columns:
# Sub = subject id
# Label = descriptive label for image pair
# ProbeCResp = the correct response, i.e. the recognition condition (1 if old, 2 if lure)
# ProbeAcc = 1 if old/lure recognition task response is correct, 0 otherwise
# RateVivid = vividness rating (1-4, NA means no response), wasn't used in current study
# Confidence = confidence rating (1-4), wasn't used in current study
# AccSub = subjects' average accuracy
# AccSubOld = subjects' average old condition accuracy
# AccSubLure = subjects' average lure condition accuracy
# ReacROI<HipAnt, HipPst, EV>Lvl<1-4, i.e. low, mid, high, semantic> = feature-specific reactivation 
#            rank measure for the specified roi and feature level during recall
# ReacROI<HipAnt, HipPst, EV> = item-specific reactivation rank measure for the specified roi during recall
#
# ROIs:
#       - HipAnt: bilateral anterior hippocampus (reactivation rank averaged over hemispheres), ROIs: 17_4, 17_5, 53_4, 53_5
#       - HipPst: bilateral posterior hippocampus (reactivation rank averaged over hemispheres), ROIs: 17_1, 17_2, 53_1, 53_2
#       - EV: bilateral calcarine sulcus (reactivation rank averaged over hemispheres), ROIs: 11145, 12145
#
# Feature Levels:
#       - 1 (low visual): reactivation rank averaged over vgg16 layers 1-4
#       - 2 (mid visual): reactivation rank averaged over vgg16 layers 5-9
#       - 3 (high visual): reactivation rank averaged over vgg16 layers 10-13
#       - 4 (semantic): reactivation rank averaged over vgg16 layers 14-16

## data to use ##
load(file="ReacRankRecall.RData") # reactivation
##



#### Reactivation Residuals ####
## replace reactivation rank values with residuals to remove shared varience between feature levels
## and the anterior and posterior hippocampus

dat2 = ReacRankRecall
dat3 = dat2

## neocortex ##
tempNames = names(dat2)[18:21]

roi = 1
for (roi in 1:1) {
  for (lvl in 1:4) {
    ivix = (1:4)[-lvl]
    dv = dat2[,tempNames[4*(roi-1)+lvl]]
    iv1 = dat2[,tempNames[4*(roi-1)+ivix[1]]]
    iv2 = dat2[,tempNames[4*(roi-1)+ivix[2]]]
    iv3 = dat2[,tempNames[4*(roi-1)+ivix[3]]]
    
    dat3[,tempNames[4*(roi-1)+lvl]] = resid(lm(dv~iv1+iv2+iv3-1))
  }
}

## hippo ant pst ##
tempNames = names(dat2)[10:17]

roi = 1
for (roi in 1:1) {
  # lvl includes ant/pst
  for (lvl in 1:8) {
    ivix = (1:8)[-lvl]
    dv = dat2[,tempNames[8*(roi-1)+lvl]]
    iv1 = dat2[,tempNames[8*(roi-1)+ivix[1]]]
    iv2 = dat2[,tempNames[8*(roi-1)+ivix[2]]]
    iv3 = dat2[,tempNames[8*(roi-1)+ivix[3]]]
    iv4 = dat2[,tempNames[8*(roi-1)+ivix[4]]]
    iv5 = dat2[,tempNames[8*(roi-1)+ivix[5]]]
    iv6 = dat2[,tempNames[8*(roi-1)+ivix[6]]]
    iv7 = dat2[,tempNames[8*(roi-1)+ivix[7]]]
    
    
    dat3[,tempNames[8*(roi-1)+lvl]] = resid(lm(dv~iv1+iv2+iv3+iv4+iv5+iv6+iv7-1))
  }
}

## hippo ant pst (item spec) ##
tempNames = names(dat2)[22:23]

roi = 1
for (roi in 1:2) {
  dv = dat2[,tempNames[roi]]
  iv1 = dat2[,tempNames[-roi]]
  
  dat3[,tempNames[roi]] = resid(lm(dv~iv1-1))
}



#### Z-score Variables ####
dat3$AccSub = scale(dat3$AccSub)
dat3$AccSubNew = scale(dat3$AccSubNew)
dat3$AccSubOld = scale(dat3$AccSubOld)
for (col in names(dat3)[10:24]) {
  dat3[,col] = scale(dat3[,col])
}



#### Figure 2: Hippo Reactivation COR Recognition Accuracy (Within-Subject) ####
# generate stats in Figure 2 using "ReacRankRecall.RData"

## Feature Specific Model ##
# run model and bootstrap (takes about 2 hours - you can skip this because the results are pre-calculated)
BootN = 1000
fmla = formula(paste0("ProbeACC ~ ProbeCResp + ",
                      "AccSubNew*ReacROIHipPstLvl1 + AccSubNew*ReacROIHipPstLvl2 + AccSubNew*ReacROIHipPstLvl3 + AccSubNew*ReacROIHipPstLvl4 - AccSubNew + ",
                      "AccSubNew*ReacROIHipAntLvl1 + AccSubNew*ReacROIHipAntLvl2 + AccSubNew*ReacROIHipAntLvl3 + AccSubNew*ReacROIHipAntLvl4 - AccSubNew + ",
                      "(1|Sub) + (1|Label)"))
nlm = glmer(fmla,data=dat3,family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun = 100000)))
save(nlm,file='lmerModel_AccCORReac_HipAntPst_InterNewAcc.RData')
nlmBoot = bootMer(nlm,function (fit) {fixef(fit)},nsim=BootN,verbose=T)$t
save(nlmBoot,file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_Boot.RData')

## Item Specific Model ##
# run model and bootstrap (takes about 1 hour - you can skip this because the results are pre-calculated)
BootN = 1000
fmla = formula(paste0("ProbeACC ~ ProbeCResp + ",
                      "AccSubNew*ReacROIHipPst + AccSubNew*ReacROIHipAnt - AccSubNew + ",
                      "(1|Sub) + (1|Label)"))
nlm = glmer(fmla,data=dat3,family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun = 100000)))
save(nlm,file='lmerModel_AccCORReacItem_HipAntPst_InterNewAcc.RData')
nlmBoot = bootMer(nlm,function (fit) {fixef(fit)},nsim=BootN,verbose=T)$t
save(nlmBoot,file='lmerModel_AccCORReacItem_HipAntPst_InterNewAcc_Boot.RData')


#### Figure 2a: interaction coef ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean(nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']>0)*2
mean(nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']>0)*2
mean(nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']>0)*2
mean(nlmBoot[,'AccSubNew:ReacROIHipAntLvl4']>0)*2
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl4']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipPstLvl4']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl4']<0)

# graph
statData = data.frame('coef'=rep(NA,8),
                      'ub'=rep(NA,8),
                      'lb'=rep(NA,8),
                      'lvl'=rep(1:4,2),
                      'roi'=rep(c('HippoAnt','HippoPst'),each=4))
statData[1:4,'coef'] = summary(nlm)$coefficients[15:18,1]
statData[5:8,'coef'] = summary(nlm)$coefficients[11:14,1]
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,14+x])[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,10+x])[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,14+x])[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,10+x])[.05*BootN]})
pdf(file=paste0('AccCORReac_HipAntPst_InterNewAcc.pdf'),8,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=lvl, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Feature Level") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_Inter')))
dev.off()


#### Figure 2b: interaction coef feat average ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean((nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']+
        nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']+
        nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']+
        nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])>0)*2
mean((nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']+
        nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']+
        nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']+
        nlmBoot[,'AccSubNew:ReacROIHipPstLvl4'])<0)
mean((nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']+
        nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']+
        nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']+
        nlmBoot[,'AccSubNew:ReacROIHipPstLvl4'])-
       (nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']+
          nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']+
          nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']+
          nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])<0)

# graph
statData = data.frame('coef'=rep(NA,2),
                      'ub'=rep(NA,2),
                      'lb'=rep(NA,2),
                      'roi'=rep(c('HippoAnt','HippoPst'),each=1))
statData[1,'coef'] = mean(summary(nlm)$coefficients[15:18,1])
statData[2,'coef'] = mean(summary(nlm)$coefficients[11:14,1])
statData[1,'ub'] = sort(rowMeans(nlmBoot[,15:18]))[.95*BootN]
statData[2,'ub'] = sort(rowMeans(nlmBoot[,11:14]))[.95*BootN]
statData[1,'lb'] = sort(rowMeans(nlmBoot[,15:18]))[.05*BootN]
statData[2,'lb'] = sort(rowMeans(nlmBoot[,11:14]))[.05*BootN]
pdf(file=paste0('AccCORReac_HipAntPst_InterNewAcc_LvlAv.pdf'),5,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=roi, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("ROI") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_Inter_FeatAv')))
dev.off()


#### Figure 2c: interaction coef item specific reac ####
BootN = 1000
load(file='lmerModel_AccCORReacItem_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReacItem_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean(nlmBoot[,'AccSubNew:ReacROIHipPst']<0)
mean(nlmBoot[,'AccSubNew:ReacROIHipAnt']>0)*2
mean(nlmBoot[,'AccSubNew:ReacROIHipPst']-nlmBoot[,'AccSubNew:ReacROIHipAnt']<0)

# graph
statData = data.frame('coef'=rep(NA,2),
                      'ub'=rep(NA,2),
                      'lb'=rep(NA,2),
                      'roi'=rep(c('HippoAnt','HippoPst'),each=1))
statData[1,'coef'] = summary(nlm)$coefficients[6,1]
statData[2,'coef'] = summary(nlm)$coefficients[5,1]
statData[1,'ub'] = sort(nlmBoot[,6])[.95*BootN]
statData[2,'ub'] = sort(nlmBoot[,5])[.95*BootN]
statData[1,'lb'] = sort(nlmBoot[,6])[.05*BootN]
statData[2,'lb'] = sort(nlmBoot[,5])[.05*BootN]
pdf(file=paste0('AccCORReac_HipAntPst_InterNewAcc_ItemSpec.pdf'),5,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=roi, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("ROI") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_Inter_ImSpecif')))
dev.off()


#### Figure 2d: coef high lure accuracy ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIHipAntLvl1']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']>0)*2
mean(nlmBoot[,'ReacROIHipAntLvl2']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']>0)*2
mean(nlmBoot[,'ReacROIHipAntLvl3']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']>0)*2
mean(nlmBoot[,'ReacROIHipAntLvl4']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl4']>0)*2
mean(nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']<0)
mean(nlmBoot[,'ReacROIHipPstLvl2']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']<0)
mean(nlmBoot[,'ReacROIHipPstLvl3']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']<0)
mean(nlmBoot[,'ReacROIHipPstLvl4']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl4']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']-(nlmBoot[,'ReacROIHipAntLvl1']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl1'])<0)
mean(nlmBoot[,'ReacROIHipPstLvl2']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']-(nlmBoot[,'ReacROIHipAntLvl2']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl2'])<0)
mean(nlmBoot[,'ReacROIHipPstLvl3']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']-(nlmBoot[,'ReacROIHipAntLvl3']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl3'])<0)
mean(nlmBoot[,'ReacROIHipPstLvl4']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl4']-(nlmBoot[,'ReacROIHipAntLvl4']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])<0)

# graph
statData = data.frame('coef'=rep(NA,16),
                      'ub'=rep(NA,16),
                      'lb'=rep(NA,16),
                      'lvl'=rep(1:4,4),
                      'roi'=rep(rep(c('HippoAnt','HippoPst'),each=4),2),
                      'group'=rep(c('lowAcc','highAcc'),each=8))
statData[1:4,'coef'] = summary(nlm)$coefficients[7:10,1]-summary(nlm)$coefficients[15:18,1]
statData[5:8,'coef'] = summary(nlm)$coefficients[3:6,1]-summary(nlm)$coefficients[11:14,1]
statData[9:12,'coef'] = summary(nlm)$coefficients[7:10,1]+summary(nlm)$coefficients[15:18,1]
statData[13:16,'coef'] = summary(nlm)$coefficients[3:6,1]+summary(nlm)$coefficients[11:14,1]
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]-nlmBoot[,14+x])[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]-nlmBoot[,10+x])[.95*BootN]})
statData[9:12,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]+nlmBoot[,14+x])[.95*BootN]})
statData[13:16,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]+nlmBoot[,10+x])[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]-nlmBoot[,14+x])[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]-nlmBoot[,10+x])[.05*BootN]})
statData[9:12,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]+nlmBoot[,14+x])[.05*BootN]})
statData[13:16,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]+nlmBoot[,10+x])[.05*BootN]})
statData2 = statData[which(statData$group=='highAcc'),]
pdf(file=paste0('AccCORReac_HipAntPst_HighNewAcc.pdf'),8,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData2, aes(x=lvl, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Feature Level") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_HighNewAcc')))
dev.off()


#### Figure 2e: coef low lure accuracy ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIHipAntLvl1']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']<0)
mean(nlmBoot[,'ReacROIHipAntLvl2']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']<0)
mean(nlmBoot[,'ReacROIHipAntLvl3']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']<0)
mean(nlmBoot[,'ReacROIHipAntLvl4']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl4']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']>0)*2
mean(nlmBoot[,'ReacROIHipPstLvl2']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']<0)
mean(nlmBoot[,'ReacROIHipPstLvl3']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']<0)
mean(nlmBoot[,'ReacROIHipPstLvl4']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl4']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']-(nlmBoot[,'ReacROIHipAntLvl1']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl1'])>0)*2
mean(nlmBoot[,'ReacROIHipPstLvl2']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']-(nlmBoot[,'ReacROIHipAntLvl2']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl2'])<0)
mean(nlmBoot[,'ReacROIHipPstLvl3']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']-(nlmBoot[,'ReacROIHipAntLvl3']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl3'])<0)
mean(nlmBoot[,'ReacROIHipPstLvl4']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl4']-(nlmBoot[,'ReacROIHipAntLvl4']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])<0)

# graph
statData = data.frame('coef'=rep(NA,16),
                      'ub'=rep(NA,16),
                      'lb'=rep(NA,16),
                      'lvl'=rep(1:4,4),
                      'roi'=rep(rep(c('HippoAnt','HippoPst'),each=4),2),
                      'group'=rep(c('lowAcc','highAcc'),each=8))
statData[1:4,'coef'] = summary(nlm)$coefficients[7:10,1]-summary(nlm)$coefficients[15:18,1]
statData[5:8,'coef'] = summary(nlm)$coefficients[3:6,1]-summary(nlm)$coefficients[11:14,1]
statData[9:12,'coef'] = summary(nlm)$coefficients[7:10,1]+summary(nlm)$coefficients[15:18,1]
statData[13:16,'coef'] = summary(nlm)$coefficients[3:6,1]+summary(nlm)$coefficients[11:14,1]
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]-nlmBoot[,14+x])[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]-nlmBoot[,10+x])[.95*BootN]})
statData[9:12,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]+nlmBoot[,14+x])[.95*BootN]})
statData[13:16,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]+nlmBoot[,10+x])[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]-nlmBoot[,14+x])[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]-nlmBoot[,10+x])[.05*BootN]})
statData[9:12,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,6+x]+nlmBoot[,14+x])[.05*BootN]})
statData[13:16,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,2+x]+nlmBoot[,10+x])[.05*BootN]})
statData2 = statData[which(statData$group=='lowAcc'),]
pdf(file=paste0('AccCORReac_HipAntPst_LowNewAcc.pdf'),8,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData2, aes(x=lvl, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Feature Level") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_LowNewAcc')))
dev.off()


#### Figure 2f: coef high low lure accuracy feat average ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean((nlmBoot[,'ReacROIHipAntLvl1']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']+
        nlmBoot[,'ReacROIHipAntLvl2']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']+
        nlmBoot[,'ReacROIHipAntLvl3']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']+
        nlmBoot[,'ReacROIHipAntLvl4']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])>0)*2
mean((nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']+
        nlmBoot[,'ReacROIHipPstLvl2']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']+
        nlmBoot[,'ReacROIHipPstLvl3']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']+
        nlmBoot[,'ReacROIHipPstLvl4']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl4'])<0)
mean((nlmBoot[,'ReacROIHipAntLvl1']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']+
        nlmBoot[,'ReacROIHipAntLvl2']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']+
        nlmBoot[,'ReacROIHipAntLvl3']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']+
        nlmBoot[,'ReacROIHipAntLvl4']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])<0)
mean((nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']+
        nlmBoot[,'ReacROIHipPstLvl2']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']+
        nlmBoot[,'ReacROIHipPstLvl3']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']+
        nlmBoot[,'ReacROIHipPstLvl4']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl4'])<0)
mean((nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']+
        nlmBoot[,'ReacROIHipPstLvl2']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']+
        nlmBoot[,'ReacROIHipPstLvl3']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']+
        nlmBoot[,'ReacROIHipPstLvl4']+nlmBoot[,'AccSubNew:ReacROIHipPstLvl4'])-
       (nlmBoot[,'ReacROIHipAntLvl1']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']+
          nlmBoot[,'ReacROIHipAntLvl2']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']+
          nlmBoot[,'ReacROIHipAntLvl3']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']+
          nlmBoot[,'ReacROIHipAntLvl4']+nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])<0)
mean((nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl1']+
        nlmBoot[,'ReacROIHipPstLvl2']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl2']+
        nlmBoot[,'ReacROIHipPstLvl3']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl3']+
        nlmBoot[,'ReacROIHipPstLvl4']-nlmBoot[,'AccSubNew:ReacROIHipPstLvl4'])-
       (nlmBoot[,'ReacROIHipAntLvl1']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl1']+
          nlmBoot[,'ReacROIHipAntLvl2']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl2']+
          nlmBoot[,'ReacROIHipAntLvl3']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl3']+
          nlmBoot[,'ReacROIHipAntLvl4']-nlmBoot[,'AccSubNew:ReacROIHipAntLvl4'])<0)

# graph
statData = data.frame('coef'=rep(NA,4),
                      'ub'=rep(NA,4),
                      'lb'=rep(NA,4),
                      'roi'=rep(c('HippoAnt','HippoPst'),2),
                      'group'=rep(c('lowAcc','highAcc'),each=2))
statData[1,'coef'] = mean(summary(nlm)$coefficients[7:10,1]-summary(nlm)$coefficients[15:18,1])
statData[2,'coef'] = mean(summary(nlm)$coefficients[3:6,1]-summary(nlm)$coefficients[11:14,1])
statData[3,'coef'] = mean(summary(nlm)$coefficients[7:10,1]+summary(nlm)$coefficients[15:18,1])
statData[4,'coef'] = mean(summary(nlm)$coefficients[3:6,1]+summary(nlm)$coefficients[11:14,1])
statData[1,'ub'] = sort(rowMeans(nlmBoot[,7:10]-nlmBoot[,15:18]))[.95*BootN]
statData[2,'ub'] = sort(rowMeans(nlmBoot[,3:6]-nlmBoot[,11:14]))[.95*BootN]
statData[3,'ub'] = sort(rowMeans(nlmBoot[,7:10]+nlmBoot[,15:18]))[.95*BootN]
statData[4,'ub'] = sort(rowMeans(nlmBoot[,3:6]+nlmBoot[,11:14]))[.95*BootN]
statData[1,'lb'] = sort(rowMeans(nlmBoot[,7:10]-nlmBoot[,15:18]))[.05*BootN]
statData[2,'lb'] = sort(rowMeans(nlmBoot[,3:6]-nlmBoot[,11:14]))[.05*BootN]
statData[3,'lb'] = sort(rowMeans(nlmBoot[,7:10]+nlmBoot[,15:18]))[.05*BootN]
statData[4,'lb'] = sort(rowMeans(nlmBoot[,3:6]+nlmBoot[,11:14]))[.05*BootN]
pdf(file=paste0('AccCORReac_HipAntPst_HighLowNewAcc_LvlAv.pdf'),6,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=group, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Lure Acc") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_FeatAv')))
dev.off()


#### Figure 2g: coef high low lure accuracy item specific reac ####
BootN = 1000
load(file='lmerModel_AccCORReacItem_HipAntPst_InterNewAcc.RData')
summary(nlm)
load(file='lmerModel_AccCORReacItem_HipAntPst_InterNewAcc_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIHipAnt']+nlmBoot[,'AccSubNew:ReacROIHipAnt']>0)*2
mean(nlmBoot[,'ReacROIHipPst']+nlmBoot[,'AccSubNew:ReacROIHipPst']<0)
mean((nlmBoot[,'ReacROIHipPst']+nlmBoot[,'AccSubNew:ReacROIHipPst'])-
       (nlmBoot[,'ReacROIHipAnt']+nlmBoot[,'AccSubNew:ReacROIHipAnt'])<0)
mean(nlmBoot[,'ReacROIHipAnt']-nlmBoot[,'AccSubNew:ReacROIHipAnt']<0)
mean(nlmBoot[,'ReacROIHipPst']-nlmBoot[,'AccSubNew:ReacROIHipPst']<0)
mean((nlmBoot[,'ReacROIHipPst']-nlmBoot[,'AccSubNew:ReacROIHipPst'])-
       (nlmBoot[,'ReacROIHipAnt']-nlmBoot[,'AccSubNew:ReacROIHipAnt'])<0)

# graph
statData = data.frame('coef'=rep(NA,2),
                      'ub'=rep(NA,2),
                      'lb'=rep(NA,2),
                      'roi'=rep(c('HippoAnt','HippoPst'),2),
                      'group'=rep(c('lowAcc','highAcc'),each=2))
statData[1,'coef'] = summary(nlm)$coefficients[4,1]-summary(nlm)$coefficients[6,1]
statData[2,'coef'] = summary(nlm)$coefficients[3,1]-summary(nlm)$coefficients[5,1]
statData[3,'coef'] = summary(nlm)$coefficients[4,1]+summary(nlm)$coefficients[6,1]
statData[4,'coef'] = summary(nlm)$coefficients[3,1]+summary(nlm)$coefficients[5,1]
statData[1,'ub'] = sort(nlmBoot[,4]-nlmBoot[,6])[.95*BootN]
statData[2,'ub'] = sort(nlmBoot[,3]-nlmBoot[,5])[.95*BootN]
statData[3,'ub'] = sort(nlmBoot[,4]+nlmBoot[,6])[.95*BootN]
statData[4,'ub'] = sort(nlmBoot[,3]+nlmBoot[,5])[.95*BootN]
statData[1,'lb'] = sort(nlmBoot[,4]-nlmBoot[,6])[.05*BootN]
statData[2,'lb'] = sort(nlmBoot[,3]-nlmBoot[,5])[.05*BootN]
statData[3,'lb'] = sort(nlmBoot[,4]+nlmBoot[,6])[.05*BootN]
statData[4,'lb'] = sort(nlmBoot[,3]+nlmBoot[,5])[.05*BootN]
pdf(file=paste0('AccCORReac_HipAntPst_HighLowNewAcc_ItemSpec.pdf'),6,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=group, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Lure Acc") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_ImSpecif')))
dev.off()



#### Figure 3: Hippo Reactivation COR Recognition Lure Accuracy (Between-Subject) ####
# generate stats in Figure 3 using "ReacRankRecall.RData"

## Feature Specific Model ##
# run model and bootstrap (you can skip this because the results are pre-calculated)
BootN = 1000
fmla = formula(paste0("AccSubNew ~ ",
                      "ReacROIHipPstLvl1 + ReacROIHipPstLvl2 + ReacROIHipPstLvl3 + ReacROIHipPstLvl4 + ",
                      "ReacROIHipAntLvl1 + ReacROIHipAntLvl2 + ReacROIHipAntLvl3 + ReacROIHipAntLvl4"))
nlm = lm(fmla,data=dat3)
save(nlm,file='lmModel_LureAccCORReac_HipAntPst_Between.RData')
subN = length(unique(dat3$Sub))
nlmBoot = matrix(NA,nrow=BootN,ncol=nrow(summary(nlm)$coefficients))
colnames(nlmBoot) = rownames(summary(nlm)$coefficients)
for (n in 1:BootN) {
  ixBoot = rep((sample(1:subN,replace=T)-1)*90,each=90) + rep(1:90,subN)
  nlm = lm(fmla,data=dat3[ixBoot,])
  nlmBoot[n,] = summary(nlm)$coefficients[,1]
}
save(nlmBoot,file='lmModel_LureAccCORReac_HipAntPst_Between_Boot.RData')

## Item Specific Model ##
# run model and bootstrap (you can skip this because the results are pre-calculated)
BootN = 1000
fmla = formula(paste0("AccSubNew ~ ",
                      "ReacROIHipPst + ReacROIHipAnt"))
nlm = lm(fmla,data=dat3)
save(nlm,file='lmModel_LureAccCORReacItem_HipAntPst_Between.RData')
subN = length(unique(dat3$Sub))
nlmBoot = matrix(NA,nrow=BootN,ncol=nrow(summary(nlm)$coefficients))
colnames(nlmBoot) = rownames(summary(nlm)$coefficients)
for (n in 1:BootN) {
  ixBoot = rep((sample(1:subN,replace=T)-1)*90,each=90) + rep(1:90,subN)
  nlm = lm(fmla,data=dat3[ixBoot,])
  nlmBoot[n,] = summary(nlm)$coefficients[,1]
}
save(nlmBoot,file='lmModel_LureAccCORReacItem_HipAntPst_Between_Boot.RData')


#### Figure 3a: coef ####
BootN = 1000
load(file='lmModel_LureAccCORReac_HipAntPst_Between.RData')
summary(nlm)
load(file='lmModel_LureAccCORReac_HipAntPst_Between_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIHipAntLvl1']>0)*2
mean(nlmBoot[,'ReacROIHipAntLvl2']<0)
mean(nlmBoot[,'ReacROIHipAntLvl3']<0)
mean(nlmBoot[,'ReacROIHipAntLvl4']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']<0)
mean(nlmBoot[,'ReacROIHipPstLvl2']<0)
mean(nlmBoot[,'ReacROIHipPstLvl3']<0)
mean(nlmBoot[,'ReacROIHipPstLvl4']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'ReacROIHipAntLvl1']<0)
mean(nlmBoot[,'ReacROIHipPstLvl2']-nlmBoot[,'ReacROIHipAntLvl2']<0)
mean(nlmBoot[,'ReacROIHipPstLvl3']-nlmBoot[,'ReacROIHipAntLvl3']<0)
mean(nlmBoot[,'ReacROIHipPstLvl4']-nlmBoot[,'ReacROIHipAntLvl4']<0)

# graph
statData = data.frame('coef'=rep(NA,8),
                      'ub'=rep(NA,8),
                      'lb'=rep(NA,8),
                      'lvl'=rep(1:4,2),
                      'roi'=rep(c('HippoAnt','HippoPst'),each=4))
statData[1:4,'coef'] = summary(nlm)$coefficients[6:9,1]
statData[5:8,'coef'] = summary(nlm)$coefficients[2:5,1]
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,5+x])[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,1+x])[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,5+x])[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,1+x])[.05*BootN]})
pdf(file=paste0('NewAccCORReac_HipAntPst_Between.pdf'),8,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=lvl, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Feature Level") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Between')))
dev.off()


#### Figure 3b: coef feat average ####
BootN = 1000
load(file='lmModel_LureAccCORReac_HipAntPst_Between.RData')
summary(nlm)
load(file='lmModel_LureAccCORReac_HipAntPst_Between_Boot.RData')

# p-values
mean((nlmBoot[,'ReacROIHipAntLvl1']+nlmBoot[,'ReacROIHipAntLvl2']+nlmBoot[,'ReacROIHipAntLvl3']+nlmBoot[,'ReacROIHipAntLvl4'])<0)
mean((nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'ReacROIHipPstLvl2']+nlmBoot[,'ReacROIHipPstLvl3']+nlmBoot[,'ReacROIHipPstLvl4'])<0)
mean((nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'ReacROIHipPstLvl2']+nlmBoot[,'ReacROIHipPstLvl3']+nlmBoot[,'ReacROIHipPstLvl4'])-
       (nlmBoot[,'ReacROIHipAntLvl1']+nlmBoot[,'ReacROIHipAntLvl2']+nlmBoot[,'ReacROIHipAntLvl3']+nlmBoot[,'ReacROIHipAntLvl4'])<0)

# graph
statData = data.frame('coef'=rep(NA,2),
                      'ub'=rep(NA,2),
                      'lb'=rep(NA,2),
                      'roi'=c('HippoAnt','HippoPst'))
statData[1,'coef'] = mean(summary(nlm)$coefficients[6:9,1])
statData[2,'coef'] = mean(summary(nlm)$coefficients[2:5,1])
statData[1,'ub'] = sort(rowMeans(nlmBoot[,6:9]))[.95*BootN]
statData[2,'ub'] = sort(rowMeans(nlmBoot[,2:5]))[.95*BootN]
statData[1,'lb'] = sort(rowMeans(nlmBoot[,6:9]))[.05*BootN]
statData[2,'lb'] = sort(rowMeans(nlmBoot[,2:5]))[.05*BootN]
pdf(file=paste0('NewAccCORReac_HipAntPst_Between_LvlAv.pdf'),5,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=roi, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("ROI") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_ImSpecif')))
dev.off()


#### Figure 3c: coef item specific reac ####
BootN = 1000
load(file='lmModel_LureAccCORReacItem_HipAntPst_Between.RData')
summary(nlm)
load(file='lmModel_LureAccCORReacItem_HipAntPst_Between_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIHipAnt']<0)
mean(nlmBoot[,'ReacROIHipPst']<0)
mean(nlmBoot[,'ReacROIHipPst']-nlmBoot[,'ReacROIHipAnt']<0)

# graph
statData = data.frame('coef'=rep(NA,2),
                      'ub'=rep(NA,2),
                      'lb'=rep(NA,2),
                      'roi'=c('HippoAnt','HippoPst'))
statData[1,'coef'] = summary(nlm)$coefficients[3,1]
statData[2,'coef'] = summary(nlm)$coefficients[2,1]
statData[1,'ub'] = sort(nlmBoot[,3])[.95*BootN]
statData[2,'ub'] = sort(nlmBoot[,2])[.95*BootN]
statData[1,'lb'] = sort(nlmBoot[,3])[.05*BootN]
statData[2,'lb'] = sort(nlmBoot[,2])[.05*BootN]
pdf(file=paste0('NewAccCORReac_HipAntPst_Between_ItemSpec.pdf'),5,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData, aes(x=roi, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("ROI") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_ImSpecif')))
dev.off()



#### Figure 4: Hippo Reactivation COR Recognition Accuracy - EV interaction (Within-Subject) ####
# generate stats in Figure 4 using "ReacRankRecall.RData"

## Feature Specific Model ##
# run model and bootstrap (takes about 3 hours - you can skip this because the results are pre-calculated)
BootN = 1000
fmla = formula(paste0("ProbeACC ~ ProbeCResp + ",
                      "ReacROIEVLvl1*ReacROIHipPstLvl1 + ReacROIEVLvl1*ReacROIHipPstLvl2 + ReacROIEVLvl1*ReacROIHipPstLvl3 + ReacROIEVLvl1*ReacROIHipPstLvl4 + ",
                      "ReacROIEVLvl1*ReacROIHipAntLvl1 + ReacROIEVLvl1*ReacROIHipAntLvl2 + ReacROIEVLvl1*ReacROIHipAntLvl3 + ReacROIEVLvl1*ReacROIHipAntLvl4 + ",
                      "AccSubNew*ReacROIEVLvl1*ReacROIHipPstLvl1 + AccSubNew*ReacROIEVLvl1*ReacROIHipPstLvl2 + AccSubNew*ReacROIEVLvl1*ReacROIHipPstLvl3 + AccSubNew*ReacROIEVLvl1*ReacROIHipPstLvl4 - AccSubNew + ",
                      "AccSubNew*ReacROIEVLvl1*ReacROIHipAntLvl1 + AccSubNew*ReacROIEVLvl1*ReacROIHipAntLvl2 + AccSubNew*ReacROIEVLvl1*ReacROIHipAntLvl3 + AccSubNew*ReacROIEVLvl1*ReacROIHipAntLvl4 - AccSubNew + ",
                      "(1|Sub) + (1|Label)"))
nlm = glmer(fmla,data=dat3,family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun = 100000)))
save(nlm,file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV.RData')
nlmBoot = bootMer(nlm,function (fit) {fixef(fit)},nsim=BootN,verbose=T)$t
save(nlmBoot,file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV_Boot.RData')


#### Figure 4b: hippo-EV interaction coef low-level ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIEVLvl1:ReacROIHipAntLvl1']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipAntLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1:ReacROIHipAntLvl1']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipAntLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1:ReacROIHipAntLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)

# graph
statData = data.frame('coef'=rep(NA,16),
                      'ub'=rep(NA,16),
                      'lb'=rep(NA,16),
                      'lvl'=rep(1:4,4),
                      'roi'=rep(rep(c('HippoAnt','HippoPst'),each=4),2),
                      'group'=rep(c('lowAcc','highAcc'),each=8))
statData[1:4,'coef'] = summary(nlm)$coefficients[16:19,1]-summary(nlm)$coefficients[33:36,1]
statData[5:8,'coef'] = summary(nlm)$coefficients[12:15,1]-summary(nlm)$coefficients[29:32,1]
statData[9:12,'coef'] = summary(nlm)$coefficients[16:19,1]+summary(nlm)$coefficients[33:36,1]
statData[13:16,'coef'] = summary(nlm)$coefficients[12:15,1]+summary(nlm)$coefficients[29:32,1]
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,15+x]-nlmBoot[,32+x])[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,11+x]-nlmBoot[,28+x])[.95*BootN]})
statData[9:12,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,15+x]+nlmBoot[,32+x])[.95*BootN]})
statData[13:16,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,11+x]+nlmBoot[,28+x])[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,15+x]-nlmBoot[,32+x])[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,11+x]-nlmBoot[,28+x])[.05*BootN]})
statData[9:12,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,15+x]+nlmBoot[,32+x])[.05*BootN]})
statData[13:16,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,11+x]+nlmBoot[,28+x])[.05*BootN]})
statData2 = statData[which(statData$lvl=='1'),]
pdf(file=paste0('AccCORReac_HipAntPstLvl1_InterEVLvl1_HighLowNewAcc.pdf'),6,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData2, aes(x=group, y=coef, fill=roi)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("HippoAnt" = "#00BFC4", "HippoPst" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Lure Acc") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within')))
dev.off()


#### Figure 4c: pHC coef low-level EV high-low ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'ReacROIHipPstLvl1:AccSubNew']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']+nlmBoot[,'ReacROIHipPstLvl1:AccSubNew']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'ReacROIHipPstLvl1:AccSubNew']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']>0)*2
mean(nlmBoot[,'ReacROIHipPstLvl1']-nlmBoot[,'ReacROIHipPstLvl1:AccSubNew']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']>0)*2

# graph
statData = data.frame('coef'=rep(NA,32),
                      'ub'=rep(NA,32),
                      'lb'=rep(NA,32),
                      'lvl'=rep(1:4,8),
                      'roi'=rep(rep(c('HippoAnt','HippoPst'),each=4),4),
                      'groupEV'=rep(rep(c('lowEV','highEV'),each=8),2),
                      'group'=rep(c('lowAcc','highAcc'),each=16))
statData[1:4,'coef'] = summary(nlm)$coefficients[8:11,1]-summary(nlm)$coefficients[25:28,1]-
  (summary(nlm)$coefficients[16:19,1]-summary(nlm)$coefficients[33:36,1])
statData[5:8,'coef'] = summary(nlm)$coefficients[4:7,1]-summary(nlm)$coefficients[21:24,1]-
  (summary(nlm)$coefficients[12:15,1]-summary(nlm)$coefficients[29:32,1])
statData[9:12,'coef'] = summary(nlm)$coefficients[8:11,1]-summary(nlm)$coefficients[25:28,1]+
  (summary(nlm)$coefficients[16:19,1]-summary(nlm)$coefficients[33:36,1])
statData[13:16,'coef'] = summary(nlm)$coefficients[4:7,1]-summary(nlm)$coefficients[21:24,1]+
  (summary(nlm)$coefficients[12:15,1]-summary(nlm)$coefficients[29:32,1])
statData[17:20,'coef'] = summary(nlm)$coefficients[8:11,1]+summary(nlm)$coefficients[25:28,1]-
  (summary(nlm)$coefficients[16:19,1]+summary(nlm)$coefficients[33:36,1])
statData[21:24,'coef'] = summary(nlm)$coefficients[4:7,1]+summary(nlm)$coefficients[21:24,1]-
  (summary(nlm)$coefficients[12:15,1]+summary(nlm)$coefficients[29:32,1])
statData[25:28,'coef'] = summary(nlm)$coefficients[8:11,1]+summary(nlm)$coefficients[25:28,1]+
  (summary(nlm)$coefficients[16:19,1]+summary(nlm)$coefficients[33:36,1])
statData[29:32,'coef'] = summary(nlm)$coefficients[4:7,1]+summary(nlm)$coefficients[21:24,1]+
  (summary(nlm)$coefficients[12:15,1]+summary(nlm)$coefficients[29:32,1])
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]-nlmBoot[,24+x]-(nlmBoot[,15+x]-nlmBoot[,32+x]))[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]-nlmBoot[,20+x]-(nlmBoot[,11+x]-nlmBoot[,28+x]))[.95*BootN]})
statData[9:12,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]-nlmBoot[,24+x]+(nlmBoot[,15+x]-nlmBoot[,32+x]))[.95*BootN]})
statData[13:16,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]-nlmBoot[,20+x]+(nlmBoot[,11+x]-nlmBoot[,28+x]))[.95*BootN]})
statData[17:20,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]+nlmBoot[,24+x]-(nlmBoot[,15+x]+nlmBoot[,32+x]))[.95*BootN]})
statData[21:24,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]+nlmBoot[,20+x]-(nlmBoot[,11+x]+nlmBoot[,28+x]))[.95*BootN]})
statData[25:28,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]+nlmBoot[,24+x]+(nlmBoot[,15+x]+nlmBoot[,32+x]))[.95*BootN]})
statData[29:32,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]+nlmBoot[,20+x]+(nlmBoot[,11+x]+nlmBoot[,28+x]))[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]-nlmBoot[,24+x]-(nlmBoot[,15+x]-nlmBoot[,32+x]))[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]-nlmBoot[,20+x]-(nlmBoot[,11+x]-nlmBoot[,28+x]))[.05*BootN]})
statData[9:12,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]-nlmBoot[,24+x]+(nlmBoot[,15+x]-nlmBoot[,32+x]))[.05*BootN]})
statData[13:16,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]-nlmBoot[,20+x]+(nlmBoot[,11+x]-nlmBoot[,28+x]))[.05*BootN]})
statData[17:20,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]+nlmBoot[,24+x]-(nlmBoot[,15+x]+nlmBoot[,32+x]))[.05*BootN]})
statData[21:24,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]+nlmBoot[,20+x]-(nlmBoot[,11+x]+nlmBoot[,28+x]))[.05*BootN]})
statData[25:28,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,7+x]+nlmBoot[,24+x]+(nlmBoot[,15+x]+nlmBoot[,32+x]))[.05*BootN]})
statData[29:32,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3+x]+nlmBoot[,20+x]+(nlmBoot[,11+x]+nlmBoot[,28+x]))[.05*BootN]})
statData2 = statData[which(statData$lvl=='1' & statData$roi=='HippoPst'),]
pdf(file=paste0('AccCORReac_HipPstLvl1_HighLowEVLvl1_HighLowNewAcc.pdf'),6,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData2, aes(x=group, y=coef, fill=groupEV)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("lowEV" = "#00BFC4", "highEV" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Lure Acc") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_ImSpecif')))
dev.off()


#### Figure 4d: EV coef low-level pHC high-low ####
BootN = 1000
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV.RData')
summary(nlm)
load(file='lmerModel_AccCORReac_HipAntPst_InterNewAcc_InterEV_Boot.RData')

# p-values
mean(nlmBoot[,'ReacROIEVLvl1']+nlmBoot[,'ReacROIEVLvl1:AccSubNew']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1']+nlmBoot[,'ReacROIEVLvl1:AccSubNew']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1']-nlmBoot[,'ReacROIEVLvl1:AccSubNew']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)
mean(nlmBoot[,'ReacROIEVLvl1']-nlmBoot[,'ReacROIEVLvl1:AccSubNew']-nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1']+nlmBoot[,'ReacROIEVLvl1:ReacROIHipPstLvl1:AccSubNew']<0)

# graph
statData = data.frame('coef'=rep(NA,32),
                      'ub'=rep(NA,32),
                      'lb'=rep(NA,32),
                      'lvl'=rep(1:4,8),
                      'roi'=rep(rep(c('HippoAnt','HippoPst'),each=4),4),
                      'groupEV'=rep(rep(c('lowHip','highHip'),each=8),2),
                      'group'=rep(c('lowAcc','highAcc'),each=16))
statData[1:4,'coef'] = summary(nlm)$coefficients[3,1]-summary(nlm)$coefficients[20,1]-
  (summary(nlm)$coefficients[16:19,1]-summary(nlm)$coefficients[33:36,1])
statData[5:8,'coef'] = summary(nlm)$coefficients[3,1]-summary(nlm)$coefficients[20,1]-
  (summary(nlm)$coefficients[12:15,1]-summary(nlm)$coefficients[29:32,1])
statData[9:12,'coef'] = summary(nlm)$coefficients[3,1]-summary(nlm)$coefficients[20,1]+
  (summary(nlm)$coefficients[16:19,1]-summary(nlm)$coefficients[33:36,1])
statData[13:16,'coef'] = summary(nlm)$coefficients[3,1]-summary(nlm)$coefficients[20,1]+
  (summary(nlm)$coefficients[12:15,1]-summary(nlm)$coefficients[29:32,1])
statData[17:20,'coef'] = summary(nlm)$coefficients[3,1]+summary(nlm)$coefficients[20,1]-
  (summary(nlm)$coefficients[16:19,1]+summary(nlm)$coefficients[33:36,1])
statData[21:24,'coef'] = summary(nlm)$coefficients[3,1]+summary(nlm)$coefficients[20,1]-
  (summary(nlm)$coefficients[12:15,1]+summary(nlm)$coefficients[29:32,1])
statData[25:28,'coef'] = summary(nlm)$coefficients[3,1]+summary(nlm)$coefficients[20,1]+
  (summary(nlm)$coefficients[16:19,1]+summary(nlm)$coefficients[33:36,1])
statData[29:32,'coef'] = summary(nlm)$coefficients[3,1]+summary(nlm)$coefficients[20,1]+
  (summary(nlm)$coefficients[12:15,1]+summary(nlm)$coefficients[29:32,1])
statData[1:4,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]-(nlmBoot[,15+x]-nlmBoot[,32+x]))[.95*BootN]})
statData[5:8,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]-(nlmBoot[,11+x]-nlmBoot[,28+x]))[.95*BootN]})
statData[9:12,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]+(nlmBoot[,15+x]-nlmBoot[,32+x]))[.95*BootN]})
statData[13:16,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]+(nlmBoot[,11+x]-nlmBoot[,28+x]))[.95*BootN]})
statData[17:20,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]-(nlmBoot[,15+x]+nlmBoot[,32+x]))[.95*BootN]})
statData[21:24,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]-(nlmBoot[,11+x]+nlmBoot[,28+x]))[.95*BootN]})
statData[25:28,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]+(nlmBoot[,15+x]+nlmBoot[,32+x]))[.95*BootN]})
statData[29:32,'ub'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]+(nlmBoot[,11+x]+nlmBoot[,28+x]))[.95*BootN]})
statData[1:4,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]-(nlmBoot[,15+x]-nlmBoot[,32+x]))[.05*BootN]})
statData[5:8,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]-(nlmBoot[,11+x]-nlmBoot[,28+x]))[.05*BootN]})
statData[9:12,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]+(nlmBoot[,15+x]-nlmBoot[,32+x]))[.05*BootN]})
statData[13:16,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]-nlmBoot[,20]+(nlmBoot[,11+x]-nlmBoot[,28+x]))[.05*BootN]})
statData[17:20,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]-(nlmBoot[,15+x]+nlmBoot[,32+x]))[.05*BootN]})
statData[21:24,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]-(nlmBoot[,11+x]+nlmBoot[,28+x]))[.05*BootN]})
statData[25:28,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]+(nlmBoot[,15+x]+nlmBoot[,32+x]))[.05*BootN]})
statData[29:32,'lb'] = sapply(1:4, function (x) {sort(nlmBoot[,3]+nlmBoot[,20]+(nlmBoot[,11+x]+nlmBoot[,28+x]))[.05*BootN]})
statData2 = statData[which(statData$lvl=='1' & statData$roi=='HippoPst'),]
pdf(file=paste0('AccCORReac_EVLvl1_HighLowHipPstLvl1_HighLowNewAcc.pdf'),6,6)
limits = aes(ymax = ub, ymin = lb)
dodge = position_dodge(width=0.9)
print(ggplot(statData2, aes(x=group, y=coef, fill=groupEV)) + 
        theme(panel.background = element_rect(fill = 'white', colour = 'grey'),text = element_text(size=20),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
        scale_fill_manual("group", values = c("lowHip" = "#00BFC4", "highHip" = "#F8766D")) +
        geom_bar(position=dodge, stat="identity") + 
        geom_errorbar(limits, position=dodge, width=0.20, size=.5) + 
        xlab("Lure Acc") +
        ylab("coefficient") +
        ggtitle(paste0('Reac_Within_ImSpecif')))
dev.off()

