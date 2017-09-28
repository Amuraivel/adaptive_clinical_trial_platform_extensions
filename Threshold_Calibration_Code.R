#######################################################################
###
### This file includes the code used to identify the posterior probability thresholds
###
#######################################################################
setwd() #directory where results are stored from simulations

source('GitHub_Functions_AdaptivePlatformDesign.R')


##########################################################################################
##########################################################################################
### Thresholds identified for "Constant" scenario

############################
### MEMs with pi_EB10 prior

pp.vec <- c(.975,.975,.975,.975,.975)
descent.vec <- rev(seq(from=.8,to=.975,by=.0025))


for(i in c(2,3,4,5)){ #loop through each segment

	mean.err <- rep(0,5)
	dw <- 1 #initialize values from descent.vec

	while(mean.err[i] <= 0.025){

		pp.vec[i] <- descent.vec[dw]
		write.table(t(pp.vec), file ='pi_eb10_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

		cn <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c2 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c3 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c4 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c5 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
	
		asd <- rbind(cn,c2,c3,c4,c5)
		diag(asd) <- 0 #set off-diagonals equal to 0 in order to calculate average type-1 error per segment
		mean.err <- colSums(asd)/4 #average type-1 error per segment, aiming for 0.025

		if(mean.err[i] >= 0.025){next}
		dw <- dw + 1
	}

	#Identify the values where mean type-1 error was just below and just above 0.025 to further hone in value
	low <- descent.vec[dw]
	hi <- descent.vec[dw-1]
	pp.vec[i] <- (low+hi)/2 #take average of low and high to begin honing in of value		

	while(mean.err[i] > 0.024 & mean.err[i] < 0.026){

		cn <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c2 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c3 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c4 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c5 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
	
		asd <- rbind(cn,c2,c3,c4,c5)
		diag(asd) <- 0 #set off-diagonals equal to 0 in order to calculate average type-1 error per segment
		mean.err <- colSums(asd)/4 #average type-1 error per segment, aiming for 0.025

		if(mean.err[i] > 0.024 & mean.err[i] < 0.026){break}
		if(mean.err[i] <= 0.024){
			hi <- pp.vec[i]
			pp.vec[i] <- (low+hi)/2
			write.table(t(pp.vec), file ='pi_eb10_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
		}
		if(mean.err[i] >= 0.026){
			low <- pp.vec[i]
			pp.vec[i] <- (low+hi)/2
			write.table(t(pp.vec), file ='pi_eb10_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
		}
	}

}

############################
### MEMs with pi_e prior

pp.vec <- c(.975,.975,.975,.975,.975)
descent.vec <- rev(seq(from=.8,to=.975,by=.0025))


for(i in c(2,3,4,5)){ #loop through each segment

	mean.err <- rep(0,5)
	dw <- 1 #initialize values from descent.vec

	while(mean.err[i] <= 0.025){

		pp.vec[i] <- descent.vec[dw]
		write.table(t(pp.vec), file ='pi_e_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

		cn <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c2 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c3 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c4 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c5 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
	
		asd <- rbind(cn,c2,c3,c4,c5)
		diag(asd) <- 0 #set off-diagonals equal to 0 in order to calculate average type-1 error per segment
		mean.err <- colSums(asd)/4 #average type-1 error per segment, aiming for 0.025

		if(mean.err[i] >= 0.025){next}
		dw <- dw + 1
	}

	#Identify the values where mean type-1 error was just below and just above 0.025 to further hone in value
	low <- descent.vec[dw]
	hi <- descent.vec[dw-1]
	pp.vec[i] <- (low+hi)/2 #take average of low and high to begin honing in of value		

	while(mean.err[i] > 0.024 & mean.err[i] < 0.026){

		cn <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c2 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c3 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c4 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c5 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
	
		asd <- rbind(cn,c2,c3,c4,c5)
		diag(asd) <- 0 #set off-diagonals equal to 0 in order to calculate average type-1 error per segment
		mean.err <- colSums(asd)/4 #average type-1 error per segment, aiming for 0.025

		if(mean.err[i] > 0.024 & mean.err[i] < 0.026){break}
		if(mean.err[i] <= 0.024){
			hi <- pp.vec[i]
			pp.vec[i] <- (low+hi)/2
			write.table(t(pp.vec), file ='pi_e_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
		}
		if(mean.err[i] >= 0.026){
			low <- pp.vec[i]
			pp.vec[i] <- (low+hi)/2
			write.table(t(pp.vec), file ='pi_e_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
		}
	}

}

############################
### Naive pooling

pp.vec <- c(.975,.975,.975,.975,.975)
descent.vec <- rev(seq(from=.8,to=.975,by=.0025))


for(i in c(2,3,4,5)){ #loop through each segment

	mean.err <- rep(0,5)
	dw <- 1 #initialize values from descent.vec

	while(mean.err[i] <= 0.025){

		pp.vec[i] <- descent.vec[dw]
		write.table(t(pp.vec), file ='pooled_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)

		cn <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c2 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c3 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c4 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c5 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
	
		asd <- rbind(cn,c2,c3,c4,c5)
		diag(asd) <- 0 #set off-diagonals equal to 0 in order to calculate average type-1 error per segment
		mean.err <- colSums(asd)/4 #average type-1 error per segment, aiming for 0.025

		if(mean.err[i] >= 0.025){next}
		dw <- dw + 1
	}

	#Identify the values where mean type-1 error was just below and just above 0.025 to further hone in value
	low <- descent.vec[dw]
	hi <- descent.vec[dw-1]
	pp.vec[i] <- (low+hi)/2 #take average of low and high to begin honing in of value		

	while(mean.err[i] > 0.024 & mean.err[i] < 0.026){

		cn <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c2 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c3 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c4 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
		c5 <- colMeans(t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg.res)))
	
		asd <- rbind(cn,c2,c3,c4,c5)
		diag(asd) <- 0 #set off-diagonals equal to 0 in order to calculate average type-1 error per segment
		mean.err <- colSums(asd)/4 #average type-1 error per segment, aiming for 0.025

		if(mean.err[i] > 0.024 & mean.err[i] < 0.026){break}
		if(mean.err[i] <= 0.024){
			hi <- pp.vec[i]
			pp.vec[i] <- (low+hi)/2
			write.table(t(pp.vec), file ='pooled_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
		}
		if(mean.err[i] >= 0.026){
			low <- pp.vec[i]
			pp.vec[i] <- (low+hi)/2
			write.table(t(pp.vec), file ='pooled_constant.txt', append=TRUE, col.names=FALSE, row.names=FALSE)
		}
	}

}


##########################################################################################
##########################################################################################
### Thresholds identified for "Both" scenarios (constant and varying underlying mortality, i.e., control both varying type-1 error and constant power)
#NOTE: pp.vec is final thresholds identified, approach is manual in that you start with pp.vec=0.975 for all segments and adjust the value segment-by-segment to identify desirable trade-off of power and type-1 error (will probably differ for each case)

p.vary <- c(.74,.61,.48,.36,.23)

tab.boththresh <- function(vn,v2,v3,v4,v5,cn,c2,c3,c4,c5){
	asdf <- rbind( c( format(round(colMeans(vn)[1:5],3), nsmall=3), paste0(round( median(vn[,6]),0)," (",round(quantile(vn[,6],.25),0),",",round(quantile(vn[,6],.75),0),")")),
		c( format(round(colMeans(v2)[1:5],3), nsmall=3), paste0(round( median(v2[,6]),0)," (",round(quantile(v2[,6],.25),0),",",round(quantile(v2[,6],.75),0),")")),
		c( format(round(colMeans(v3)[1:5],3), nsmall=3), paste0(round( median(v3[,6]),0)," (",round(quantile(v3[,6],.25),0),",",round(quantile(v3[,6],.75),0),")")),
		c( format(round(colMeans(v4)[1:5],3), nsmall=3), paste0(round( median(v4[,6]),0)," (",round(quantile(v4[,6],.25),0),",",round(quantile(v4[,6],.75),0),")")),
		c( format(round(colMeans(v5)[1:5],3), nsmall=3), paste0(round( median(v5[,6]),0)," (",round(quantile(v5[,6],.25),0),",",round(quantile(v5[,6],.75),0),")")),
		c( format(round(colMeans(cn)[1:5],3), nsmall=3), paste0(round( median(cn[,6]),0)," (",round(quantile(cn[,6],.25),0),",",round(quantile(cn[,6],.75),0),")")),
		c( format(round(colMeans(c2)[1:5],3), nsmall=3), paste0(round( median(c2[,6]),0)," (",round(quantile(c2[,6],.25),0),",",round(quantile(c2[,6],.75),0),")")),
		c( format(round(colMeans(c3)[1:5],3), nsmall=3), paste0(round( median(c3[,6]),0)," (",round(quantile(c3[,6],.25),0),",",round(quantile(c3[,6],.75),0),")")),
		c( format(round(colMeans(c4)[1:5],3), nsmall=3), paste0(round( median(c4[,6]),0)," (",round(quantile(c4[,6],.25),0),",",round(quantile(c4[,6],.75),0),")")),
		c( format(round(colMeans(c5)[1:5],3), nsmall=3), paste0(round( median(c5[,6]),0)," (",round(quantile(c5[,6],.25),0),",",round(quantile(c5[,6],.75),0),")")))
	return(asdf)
}

###EB with 0.10 constraint, burn-in period of 60 subjects (30 per arm) and 5 blocks of 28 subjects to adaptively randomize
pp.vec <- c(.975,.975,.9775,.975,.975) 

vn <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=rep(1,5), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v2 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v3 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v4 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v5 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))

cn <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c2 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c3 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c4 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c5 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='opt.source_constrain',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))

xtable(tab.boththresh(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5))

####Thresholds to control both varying type-1 error and constant power
###Uniform prior, burn-in period of 60 subjects (30 per arm) and 5 blocks of 28 subjects to adaptively randomize
pp.vec <- c(.975,.9815,.985,.9875,.9875) 

vn <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=rep(1,5), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v2 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v3 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v4 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v5 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))

cn <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c2 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c3 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c4 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c5 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='equal',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))

xtable(tab.boththresh(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5))

####Thresholds to control both varying type-1 error and constant power
###Pooled case, burn-in period of 60 subjects (30 per arm) and 5 blocks of 28 subjects to adaptively randomize
pp.vec <- c(.975,.995,.999,.9975,.999) 

vn <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=rep(1,5), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v2 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v3 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v4 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
v5 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=p.vary, arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))

cn <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c2 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,.7,1,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c3 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,.7,1,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c4 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,.7,1), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))
c5 <- t(sapply(1:2000, function(x) ebola.sim(pp=pp.vec, p=rep(.4,5), arm.RR=c(1,1,1,1,.7), nvec=rep(100,6), MEM=T, seed=x, prior='pool',constraint=.1,AR=T,n.burn=60,block.num=5)$seg_ntrt_res))

xtable(tab.boththresh(vn, v2, v3, v4, v5,cn, c2, c3, c4, c5))





