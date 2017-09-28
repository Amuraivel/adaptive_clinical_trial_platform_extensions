setwd() #directory where functions file is located
source('GitHub_Functions_AdaptivePlatformDesign.R')library(xtable)

###TABLE 4 in manuscript
xvec.test <- c(50,52,45,65)
nvec.test <- rep(100,4)
res <- cbind(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=rep(1,4), bvec=rep(1,4), prior='equal')$mod.mat, 
round(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=rep(1,4), bvec=rep(1,4), prior='equal')$q,4),
round(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=c(1,1,1,1), bvec=c(1,1,1,1), prior='opt.source_constrain',constrain=1)$q,4),
round(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=c(1,1,1,1), bvec=c(1,1,1,1), prior='opt.source_constrain',constrain=.9)$q,4),
round(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=c(1,1,1,1), bvec=c(1,1,1,1), prior='opt.source_constrain',constrain=.5)$q,4),
round(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=c(1,1,1,1), bvec=c(1,1,1,1), prior='opt.source_constrain',constrain=.1)$q,4),
round(calc.MEM.betabin(xvec= xvec.test, nvec=nvec.test, avec=c(1,1,1,1), bvec=c(1,1,1,1), prior='opt.source_constrain',constrain=0)$q,4))
res

xtable(res, digits=c(rep(0,5),rep(3,6)))

