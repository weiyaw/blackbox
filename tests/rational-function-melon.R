## Proposal distribution and number of iterations depend on the cooling schedule and problem at hand.

######## fit a constrained model. ########
## Visualise the dataset
load_all()

## three temp: 15, 20, 30, use 15
melon15 <- melon[melon$temp == 15, 1:2]
scaled_melon <- as.data.frame(scale(melon15, FALSE, c(max(melon15$days), 10)))
## plot(height ~ days, data = scaled_melon,
##      main = "Cucumis melo at 15",
##      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

## Fit a rational model f(x) = (a + b*x + c*x^2) / (1 + d*x + e*x^2), subject to constraint
## Construct loss functions and an indicator
loss <- rational_rss_fac(scaled_melon$height, scaled_melon$days, 3L, 2L)
is_monotone <- is_monotones_fac(3L, 2L, min(scaled_melon$days),
                                max(scaled_melon$days), TRUE)

## Get starting values using the rearrange and regress method
rar_melon <- coef(lm(height ~ days + I(days^2) + I(-days*height) + I(-days^2*height),
                     scaled_melon))
names(rar_melon) <- letters[1:5]

library(doParallel)
registerDoParallel(cores = 8)

recip_time <- system.time(recip_melon <- foreach(i = 1:40) %dopar% {
    ## Randomly generate 1000 starting values arounf rar_hyper
    set.seed(100 + i)
    starting <- get_feasibles(rar_melon, is_monotone, 1000, dist_para = list(scale = 2))

    ## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 3000, 2 comp
    ## Use a more aggresive proposal
    rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 2)
    runtime <- system.time(rm1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 3000, 1000, FALSE, TRUE))
    rm1$runtime <- runtime
    rm1
})
cat("RECIP DONE.", recip_time[3], "\n")
## plot_rat(rm1$state[1:3], c(1, rm1$state[4:5]), col = 'red')
saveRDS(recip_melon, "~/Dropbox/honours/extras/melon/recip_melon.rds")

recip85_time <- system.time(recip85_melon <- foreach(i = 1:40) %dopar% {
  ## Randomly generate 1000 starting values arounf rar_hyper
  set.seed(100 + i)
  starting <- get_feasibles(rar_melon, is_monotone, 1000, dist_para = list(scale = 2))
  
  ## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 3000, 2 comp
  ## Use a more aggresive proposal
  rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 2)
  runtime <- system.time(rm1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule_0.85, 3000, 1000, FALSE, TRUE))
  rm1$runtime <- runtime
  rm1
})
cat("RECIP DONE.", recip85_time[3], "\n")
## plot_rat(rm1$state[1:3], c(1, rm1$state[4:5]), col = 'red')
saveRDS(recip85_melon, "~/Dropbox/honours/extras/melon/recip85_melon.rds")


log_time <- system.time(log_melon <- foreach(i = 1:40) %dopar% {
  ## Randomly generate 1000 starting values arounf rar_hyper
  set.seed(100 + i)
  starting <- get_feasibles(rar_melon, is_monotone, 1000, dist_para = list(scale = 2))

  ## SMCSA, log, fraction = 0.998, sd = 1, 3000 iter, N = 1000, 2 comp
  rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.998, n_pertub = 2)
  runtime <- system.time(rm2 <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 3000, FALSE, TRUE))
  rm2$runtime <- runtime
  rm2
})
cat("LOG DONE.", log_time[3], "\n")
saveRDS(log_melon, "~/Dropbox/honours/extras/melon/log_melon.rds")
## plot_rat(rm2$state[1:3], c(1, rm2$state[4:5]), col = 'blue')


melon_pilot <- coef(nls(height ~ (a + b*days + c*days^2) / (1 + d*days + e*I(days^2)),
                        scaled_melon, start = rar_melon))

cepso_time <- system.time(cepso_melon <- foreach(i = 1:40) %dopar% {
  ## Randomly generate 1000 starting values arounf rar_hyper
  set.seed(100 + i)
  starting <- get_feasibles(rar_melon, is_monotone, 1000, dist_para = list(scale = 2))

  ## CEPSO, 2000 iter, N = 3000, neighbour = 3000/5
  ## OLS as pilot estimate
  runtime <- system.time(rm3 <- CEPSO(starting, loss, melon_pilot, is_monotone, iter = 2000, N = 3000, verbose = TRUE))
  rm3$runtime <- runtime
  rm3
})
cat("CEPSO DONE.", cepso_time[3], "\n")
saveRDS(cepso_melon, "~/Dropbox/honours/extras/melon/cepso_melon.rds")
## plot_rat(rm3$theta[1:3], c(1, rm3$theta[4:5]), col = 'red')

