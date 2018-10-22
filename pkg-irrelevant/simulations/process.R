rm(list = ls())
hyper <- readRDS("hyper.rds")
hyper.rob <- readRDS("hyper-rob.rds")

pdf("~/Dropbox/honours/talk/rss-smc-sa/tex/plots/hyper.pdf", width = 8, height = 5)
plot(y ~ x, data = hyper, main = "A simulated dataset with a sigmoidal trend",
     xlab = "time", ylab = "height of a plant", cex.lab = 1.3, cex.axis = 1.3)
dev.off()

pdf("~/Dropbox/honours/talk/rss-smc-sa/tex/plots/hyper-rob.pdf", width = 8, height = 5)
plot(y ~ x, data = hyper.rob,
     main = "A simulated dataset with outliers",
     xlab = "time", ylab = "height of a plant", cex.lab = 1.3, cex.axis = 1.3)
dev.off()

source("generic-subor.R")

## Gives the information for 'hyper' dataset. This function is model-specific.
## pic: gives a picture?
## out: supply a path ("string") if wanna export the pictures
info.hyper <- function(res, dat, pic = FALSE, out = FALSE) {
    n <- length(res)
    objv <- rep(NA, len = n)

    if (pic) {
        if (is.character(out)) {
            pdf(file = out, width = 8, height = 5)
        }
        ## Data
        plot(y ~ x, dat, main = "A simulated dataset with outliers",
             xlab = "time", ylab = "height of a plant", cex.lab = 1.3, cex.axis = 1.3)

        ## Actual function
        ## lines(1 + tanh(x - 3) ~ x, dat, col = 'blue')
    }
        ## Resulting curves
    for (i in seq_len(length(res))) {
        objv[i] <- res[[i]]$objv
        if (pic) {
            with(res[[i]], plot_rat(theta[1:2], c(1, theta[3:4]), range(dat$x),
                                   col = 'grey', lty = 1, lwd = 0.25))
        }
    }

    ## Unconstrained least squares
    rar.hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), dat))
    names(rar.hyper) <- letters[1:4]
    un.hyper <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), dat, rar.hyper)
    cat(crossprod(resid(un.hyper)), "\n")

    ## Gives the best curve (the minima)
    if (pic) {
        ## plot_rat(coef(un.hyper)[1:2], c(1, coef(un.hyper)[3:4]), range(dat$x),
        ##         col = 'black', lty = 1)
        best <- which.min(objv)
        with(res[[best]], plot_rat(theta[1:2], c(1, theta[3:4]), range(dat$x),
                                  col = 'black'))
    }

    info <- c(mean = mean(objv), sd = sd(objv),
              quantile(objv))
    ## print(round(info, digits = 4))
    cat(round(info, digits = 3), "\n", sep = " & ")
    best <- res[[which.min(objv)]]
    best$all.objv <- objv
    invisible(best)
    if (is.character(out) && pic) {
        dev.off()
    }
}



## Reciprocal hyper
sa.hyper <- readRDS("./recip/sa-hyper.rds")
smc.hyper <- readRDS("./recip/smc-hyper.rds")
info.hyper(sa.hyper, hyper, T, "../1-paper/tex/plots/recip/sa-hyper.pdf")
info.hyper(smc.hyper, hyper, T,
           "~/Dropbox/honours/talk/rss-smc-sa/tex/plots/hyper-rss-mono.pdf")

## CEPSO hyper
cep.hyper <- readRDS("./simulations/cepso/cepso-hyper.rds")
info.hyper(cep.hyper, hyper, T, "../1-paper/tex/plots/cepso/cepso-hyper.pdf")

## Logarithmic hyper
sa.log.hyper <- readRDS("./simulations/log/sa-log-hyper.rds")
smc.log.hyper <- readRDS("./simulations/log/smc-log-hyper.rds")
info.hyper(sa.log.hyper, hyper, T, "../1-paper/tex/plots/log/sa-hyper.pdf")
info.hyper(smc.log.hyper, hyper, T, "../1-paper/tex/plots/log/smc-hyper.pdf")


## Reciprocal hyper robust
sa.hyper.rob <- readRDS("./simulations/recip/sa-hyper-rob.rds")
smc.hyper.rob <- readRDS("./recip/smc-hyper-rob.rds")
info.hyper(sa.hyper.rob, hyper.rob)
info.hyper(smc.hyper.rob, hyper.rob, T, "~/Dropbox/honours/talk/rss-smc-sa/tex/plots/hyper-rob-mono.pdf")

## CEPSO hyper robust
cep.hyper.rob <- readRDS("./simulations/cepso/cepso-hyper-rob.rds")
info.hyper(cep.hyper.rob, hyper.rob, T, "../paper/tex/plots/cepso/cepso-hyper-rob.pdf")

## Logarithmic hyper robust
sa.log.hyper.rob <- readRDS("./simulations/log/sa-hyper-rob.rds")
smc.log.hyper.rob <- readRDS("./simulations/log/smc-hyper-rob.rds")
info.hyper(sa.log.hyper.rob, hyper.rob)
info.hyper(smc.log.hyper.rob, hyper.rob)

bad <- 0
for (i in 1:100) {
    if (cep.hyper.rob[[i]]$objv > 4) {bad <- bad + 1}
}



