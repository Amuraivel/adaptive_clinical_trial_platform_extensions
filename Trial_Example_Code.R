setwd() #directory where functions file is located

source('GitHub_Functions_AdaptivePlatformDesign.R')

###Example trial (2 segments) with calculations for relevant aspects
##Use observed mortality in trial of 37% in control arm and 22% in treatment arm for first segment
##Assume 22% vs 11% in second segment

##Segment 1 (PREVAIL II and MEM approach)
round(ebola.posterior.prob(na=20, nb=20, xa=4, xb=7, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=40, nb=40, xa=9, xb=15, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=60, nb=60, xa=13, xb=22, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=80, nb=80, xa=18, xb=30, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=100, nb=100, xa=22, xb=37, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)

##Segment 2 (PREVAIL II, no borrowing, 22% vs 11%)
round(ebola.posterior.prob(na=20, nb=20, xa=2, xb=4, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=40, nb=40, xa=4, xb=9, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=60, nb=60, xa=7, xb=13, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=80, nb=80, xa=9, xb=18, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)
round(ebola.posterior.prob(na=100, nb=100, xa=11, xb=22, a=1, b=1, delta=0, mod.mat=NULL, mod.weight=NULL), 4)

##Segment 2 (MEM with pi_e, 22% vs 11%)
#Intialize
R.blck <- 160
n.a <- n.b <- 20 #40 observed before first interim analysis

#Interim Analysis 1 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(4, 22), nvec = c(20, 100), avec=c(1,1), bvec=c(1,1), prior='equal')$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(4, 22), nvec = c(20, 100), avec=c(1,1), bvec=c(1,1), prior='equal')
round(ebola.posterior.prob(nb=c(20,100), xb=c(4,22), na=20, xa=2, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Interim Analysis 2 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(7, 22), nvec = c(30, 100), avec=c(1,1), bvec=c(1,1), prior='equal')$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(7, 22), nvec = c(30, 100), avec=c(1,1), bvec=c(1,1), prior='equal')
round(ebola.posterior.prob(nb=c(30,100), xb=c(7,22), na=50, xa=6, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Interim Analysis 3 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(9, 22), nvec = c(39, 100), avec=c(1,1), bvec=c(1,1), prior='equal')$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(9, 22), nvec = c(39, 100), avec=c(1,1), bvec=c(1,1), prior='equal')
round(ebola.posterior.prob(nb=c(39,100), xb=c(9,22), na=81, xa=9, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Interim Analysis 4 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(11, 22), nvec = c(48, 100), avec=c(1,1), bvec=c(1,1), prior='equal')$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(11, 22), nvec = c(48, 100), avec=c(1,1), bvec=c(1,1), prior='equal')
round(ebola.posterior.prob(nb=c(48,100), xb=c(11,22), na=112, xa=12, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Final analysis
esss.blck <- summary_calc.betabin(xvec = c(13, 22), nvec = c(57, 100), avec=c(1,1), bvec=c(1,1), prior='equal')$esss

mem.est_im <- summary_calc.betabin(xvec = c(13, 22), nvec = c(57, 100), avec=c(1,1), bvec=c(1,1), prior='equal')
round(ebola.posterior.prob(nb=c(57,100), xb=c(13,22), na=143, xa=16, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)



##Segment 2 (MEM with pi_EB10, 22% vs 11%)
#Intialize
R.blck <- 160
n.a <- n.b <- 20 #40 observed before first interim analysis

#Interim Analysis 1 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(4, 22), nvec = c(20, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(4, 22), nvec = c(20, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)
round(ebola.posterior.prob(nb=c(20,100), xb=c(4,22), na=20, xa=2, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Interim Analysis 2 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(8, 22), nvec = c(36, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(8, 22), nvec = c(36, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)
round(ebola.posterior.prob(nb=c(36,100), xb=c(8,22), na=44, xa=5, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Interim Analysis 3 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(11, 22), nvec = c(51, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(11, 22), nvec = c(51, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)
round(ebola.posterior.prob(nb=c(51,100), xb=c(11,22), na=69, xa=7, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Interim Analysis 4 and AR for next block
esss.blck <- summary_calc.betabin(xvec = c(14, 22), nvec = c(65, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)$esss
pi.est <- 0.5 * ( ((esss.blck + n.b - n.a)/R.blck) + 1)
tau <- max(0, min(c(pi.est, 1)))

n.a_blck <- round(40*tau) #block.size_vec[blck]=40 for example
n.b_blck <- 40 - n.a_blck #block.size_vec[blck]=40 for example
R.blck <- R.blck - 40

n.a <- n.a + n.a_blck
n.b <- n.b + n.b_blck

mem.est_im <- summary_calc.betabin(xvec = c(14, 22), nvec = c(65, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)
round(ebola.posterior.prob(nb=c(65,100), xb=c(14,22), na=95, xa=10, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
round(tau,3)
n.a_blck; n.b_blck

#Final analysis
esss.blck <- summary_calc.betabin(xvec = c(17, 22), nvec = c(79, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)$esss

mem.est_im <- summary_calc.betabin(xvec = c(17, 22), nvec = c(79, 100), avec=c(1,1), bvec=c(1,1), prior='opt.source_constrain', constraint=.1)
round(ebola.posterior.prob(nb=c(79,100), xb=c(17,22), na=121, xa=13, mod.mat=mem.est_im$mod.mat, mod.weight=mem.est_im$mod.weight, a=1, b=1), 4)
round(esss.blck,1)
