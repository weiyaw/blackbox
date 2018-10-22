p <- function(x) { if (x > 0) {dnorm(x)}
                   else {0}}
theta.old <- 0.0
thetasamples <- rep(0,20000)
for (i in 1:20000) {
    theta.new <- rnorm(1,1,1)
    post.new <- p(theta.new)*dnorm(theta.old, 1, 1)
    post.old <- p(theta.old)*dnorm(theta.new, 1, 1)
    if (runif(1)<(post.new/post.old)) { theta.old <- theta.new }
    thetasamples[i] <- theta.old
}
mean(thetasamples[10001:20000])
hist(thetasamples[10001:20000])
acf(thetasamples[10001:20000])
plot(thetasamples[10001:20000], type = 'l')

simlidar <- read.csv("data/simlidar.csv")
x <- simlidar$s.range
y <- simlidar$s.y
plot(y ~ x)

nl <- nls(y ~ (a + b*x + c*x^2) / (1 + d*x + e*x^2), start = list(a = 1, b = 2, c = 3, d = 4, e = 5))

source("main.R")
source("subor.R")
source("generic.R")
source("generic-subor.R")
indicatorAll <- indicatorAllFac(3L,2L,xmax = 2)

starting <- matrix(nrow = 5, ncol = 1000)
for (i in 1:1000) {
    temp <- genMonoRatOnly(coef(nl)[1:3], c(1, coef(nl)[4:5]), c(T,T,T), c(F,T,T), 0, 2)
    starting[1:3, i] <- temp$f
    starting[4:5, i] <- temp$g[-1]
}
proposal <- proposalMonoFac(dist.para = list(mean = 0, sd = 1), indicator = indicatorAll)
rm1 <- SMCSA(rationalResidFac(y, x, 3L, 2L), proposal, starting, schedule, 5000, 100, T)
plotRat(coef(nl)[1:3], c(1, coef(nl)[4:5]), c(0, 5))
plot(y ~ x)
plotRat(rm1$theta[1:3], c(1, rm1$theta[4:5]), c(0, 5), col = 2)


rm2 <- monoRat(y, x, coef(nl)[1:3], c(1, coef(nl)[4:5]), limits = c(0, 2))
plotRat(rm2$f, rm2$g, c(0, 2), col = 4)

proposal <- proposalMonoFac(dist.para = list(mean = 0, sd = 1), indicator = indicatorAll)
rm3 <- multiSA(rationalResidFac(y, x, 3L, 2L), proposal, starting, schedule, 5000, 100, T)
plot(y ~ x)
plotRat(rm3$theta[1:3], c(1, rm3$theta[4:5]), c(0, 5), col = 2, T)


