rm(list = ls())
source("./rational-function.R")
source("./generic-subor.R")
## get_feasibles
is_monotones <- is_monotones_fac(2, 1)
all(is_monotones(get_feasibles(matrix(1:3, 3, 1), is_monotones, 10)))
all(is_monotones(get_feasibles(matrix(1:12, 3, 4), is_monotones, 4)))
all(is_monotones(get_feasibles(matrix(1:12, 3, 4), is_monotones, 3)))
all(is_monotones(get_feasibles(matrix(1:12, 3, 4), is_monotones, 10000)))


## rnorm_rw_fac and rtnorm_rw_fac (with a dummy indicator)
## when k = 1, they all are equivalent
set.seed(250318)
rnorm_rw_fac(sd = 5)(matrix(1:12, 3, 4), 1)
set.seed(250318)
rtnorm_rw_fac(function(x) TRUE, sd = 5)(matrix(1:12, 3, 4), 1)
set.seed(250318)
matrix(1:12, 3, 4) + rnorm(12, sd = 5 * 1)

## when k = 2, they all are equivalent
set.seed(250318)
rnorm_rw_fac(sd = 5)(matrix(1:12, 3, 4), 2)
set.seed(250318)
rtnorm_rw_fac(function(x) TRUE, sd = 5)(matrix(1:12, 3, 4), 2)
set.seed(250318)
matrix(1:12, 3, 4) + rnorm(12, sd = 5 * 0.97)


## rtnorm_rw_fac (with a dummy indicator)
all(is_monotones(rtnorm_rw_fac(is_monotones, sd = 3)(matrix(1:12, 3, 4), 1)))
all(is_monotones(rtnorm_rw_fac(is_monotones, sd = 4)(matrix(1:12, 3, 4), 2)))
all(is_monotones(rtnorm_rw_fac(is_monotones, sd = 5)(matrix(1:12, 3, 4), 3)))


## recip_schedule and log_schedule
recip_schedule(1, 10)
recip_schedule(2, 10)
log_schedule(1, 10)
log_schedule(2, 10)
