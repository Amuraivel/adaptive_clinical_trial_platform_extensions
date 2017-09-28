setwd() #directory where simulation results are located

############################################################
############################################################
###Bivariate scatterplot of power and average type-1 error rate for RR=0.7, PP thresholds identified under constant case, adaptive randomization, and interim monitoring
###Values taken from tables included in manuscript

###Constant mortality
plot(x=-100,y=-100,ylim=c(.4,.8), xlim=c(.02,.03), ylab='Power', xlab='Average Type-1 Error Rate Across Scenarios')
title('Constant Mortality')
cex.val <- 1.5
plot.legend <- c('PREVAIL II',expression(paste('MEM ',pi[EB[10]])),expression(paste('MEM ',pi[e])),'Naive Pooling')
legend('topleft',pch=c(6,15,17,18,15,15,15,15),col=c(rep('gray65',4),'black','gray65','orangered2','blue'),legend=c(plot.legend,'Segment 2','Segment 3','Segment 4','Segment 5'))

#PREVAIL-II
points( y=c(.432,.431,.434,.441), x=c(.028,.02875,.030,.02875), pch=6, col=c('black','gray65','orangered2','blue'), cex=cex.val)
#MEM with EB10 prior
points( y=c(.470,.509,.548,.560), x=c(.026,.02575,.02975,.02575), pch=15, col=c('black','gray65','orangered2','blue'), cex=cex.val)
#MEM with uniform/equal weight prior
points( y=c(.556,.611,.694,.745), x=c(.027,.02125,.02675,.0285), pch=17, col=c('black','gray65','orangered2','blue'), cex=cex.val)
#MEM with pooled prior
points( y=c(.581,.682,.731,.778), x=c(.022,.02525,.02675,.02925), pch=18, col=c('black','gray65','orangered2','blue'), cex=cex.val)

###Varying mortality
plot(x=-100,y=-100,ylim=c(.2,1), xlim=c(.02,.6), ylab='Power', xlab='Average Type-1 Error Rate Across Scenarios')
title('Varying Mortality')
cex.val <- 1.5
plot.legend <- c('PREVAIL II',expression(paste('MEM ',pi[EB[10]])),expression(paste('MEM ',pi[e])),'Naive Pooling')
legend('bottomright',pch=c(6,15,17,18,15,15,15,15),col=c(rep('gray65',4),'black','gray65','orangered2','blue'),legend=c(plot.legend,'Segment 2','Segment 3','Segment 4','Segment 5'))

#PREVAIL-II
points( y=c(.764,.555,.386,.235), x=c(.030,.03075,.02625,.026), pch=6, col=c('black','gray65','orangered2','blue'), cex=cex.val)
#MEM with EB10 prior
points( y=c(.781,.637,.511,.354), x=c(.040,.0465,.05775,.06425), pch=15, col=c('black','gray65','orangered2','blue'), cex=cex.val)
#MEM with uniform/equal weight prior
points( y=c(.854,.768,.717,.637), x=c(.102,.14725,.199,.27325), pch=17, col=c('black','gray65','orangered2','blue'), cex=cex.val)
#MEM with pooled prior
points( y=c(.997,.986,.947,.955), x=c(.342,.5895,.45175,.6085), pch=18, col=c('black','gray65','orangered2','blue'), cex=cex.val)


############################################################
############################################################
###Boxplots

###Read in files for constant mortality and varying mortality cases
#Column 1: combination of priors/other parameters from simulation (1=EB10, 3=uniform, 5=PREVAIL-II, 6=pooled: all with PP thresholds set under constant mortality)
#Columns 2-6: N per each scenario of null case/2/3/4/5
#Columns 7-11: proportion assigned to treatment in null/2/3/4/5 scenarios from segments 2-5
#Column 12: proportion surviving from segments 2-5 in null scenario
#Columns 13-16: proportion surviving in non-null segment for scenario 2/3/4/5

datc <- read.table('paper2_simresults_constant.txt', header=F) #load results for constant mortality scenario
datv <- read.table('paper2_simresults_varying.txt', header=F) #load results for varying mortality scenario

###Select 
dat <- datc[which(datc[,1] %in% c(1,3,5,6)),] #only select simulations used in main manuscript for plotting
dat <- datv[which(datv[,1] %in% c(1,3,5,6)),] #only select simulations used in main manuscript for plotting

colnames(dat) <- c('simnum','n1','n2','n3','n4','n5','t1','t2','t3','t4','t5','s1','s2','s3','s4','s5')
dat$AK_simrenum <- NA
dat$AK_simrenum[which(dat$simnum==5)] <- 1 #relist PREVAIL-II as 1
dat$AK_simrenum[which(dat$simnum==1)] <- 2 #relist EB10 as 2
dat$AK_simrenum[which(dat$simnum==3)] <- 3 #relist uniform prior as 3
dat$AK_simrenum[which(dat$simnum==6)] <- 4 #relist pooled prior as 4

datT <- reshape(dat, varying = c('t1','t2','t3','t4','t5'), v.names = 'proptrt', timevar = 'segment', times = c(1,2,3,4,5), direction = 'long')
datS <- reshape(dat, varying = c('s1','s2','s3','s4','s5'), v.names = 'propsurv', timevar = 'segment', times = c(1,2,3,4,5), direction = 'long')

boxplot(proptrt ~ AK_simrenum*segment, data=datT, at=c(1,1.8,2.6,3.4, 5,5.8,6.6,7.4, 9,9.8,10.6,11.4, 13,13.8,14.6,15.4, 17,17.8,18.6,19.4), width=rep(2,20), xaxt='n', ylab='Proportion Randomized to Treatment', col=c('white','gray75','orangered2','blue'))
axis(1, at=c(2.2, 6.2, 10.2, 14.2, 18.2), c('Null','Scenario 2','Scenario 3','Scenario 4','Scenario 5'), tick=F)
abline(v=c(4.2,8.2,12.2,16.2), lty=2, col='gray65')
plot.legend <- c('PREVAIL II',expression(paste('MEM ',pi[EB[10]])),expression(paste('MEM ',pi[e])),'Naive Pooling')
legend('top',bty='n',fill=c('white','gray75','orangered2','blue'), legend=plot.legend, horiz=T, inset=-.07, xpd=T)
title('(a)') #add title to help reference in manuscript, choose a/c depending on constant/varying mortality

boxplot(propsurv ~ AK_simrenum*segment, data=datS, at=c(1,1.8,2.6,3.4, 5,5.8,6.6,7.4, 9,9.8,10.6,11.4, 13,13.8,14.6,15.4, 17,17.8,18.6,19.4), width=rep(2,20), xaxt='n', ylab='Proportion Surviving', col=c('white','gray75','orangered2','blue'))
axis(1, at=c(2.2, 6.2, 10.2, 14.2, 18.2), c('Null','Scenario 2','Scenario 3','Scenario 4','Scenario 5'), tick=F)
abline(v=c(4.2,8.2,12.2,16.2), lty=2, col='gray65')
plot.legend <- c('PREVAIL II',expression(paste('MEM ',pi[EB[10]])),expression(paste('MEM ',pi[e])),'Naive Pooling')
legend('top',bty='n',fill=c('white','gray75','orangered2','blue'), legend=plot.legend, horiz=T, inset=-.07, xpd=T)
title('(b)') #add title to help reference in manuscript, choose b/d depending on constant/varying mortality

