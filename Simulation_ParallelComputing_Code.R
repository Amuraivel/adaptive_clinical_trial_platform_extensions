###Code to run simulation scenarios with parallel computing via snowfall package
nsim=NA

library(snowfall)
sfInit(cpus = 12, type = 'SOCK', parallel = T)

sfLibrary(xtable)
sfExportAll()

sfSource("./GitHub_Functions_AdaptivePlatformDesign.R")

###Run simulation scenarios which have one effective drug, no futility monitoring
sfSapply(1:16, function(x) server.sim(sim.num=x) )

###Run simulation scenarios which have two effective drugs, no futility monitoring
sfSapply(1:4, function(x) server.sim_2nonnull(sim.num=x) )

###Run simulation scenarios which have one effective drug and include futility monitoring
sfSapply(1:16, function(x) server.sim_futility(sim.num=x) )

sfStop()