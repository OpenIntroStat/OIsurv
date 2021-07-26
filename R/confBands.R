confBands <- function(
  x,
  confType = c('plain', 'log-log', 'asin-sqrt'),
  confLevel = c(0.90, 0.95, 0.99),
  type = c('ep', 'hall'),
  tL,
  tU
){
  confLevel <- match.arg(confLevel)
  confType <- match.arg(confType)
  type <- match.arg(type)
  fit  <- survfit(x ~ 1)
  n    <- length(x) / 2
  time <- summary(fit)$time
  std.err <- summary(fit)$std.err
  surv    <- summary(fit)$surv
  t.L <- ifelse(missing(tL), head(time, 1), tail(time[time >= tL], 1))
  t.U <- ifelse(missing(tU), max(x[x[, 2] == 1, 1]), tail(time[time <= tU], 1))
  t.L.pos <- which(time == t.L)
  t.U.pos <- which(time == t.U)
  sigma.2 <- (std.err / surv)^2
  a.L     <- n * sigma.2[time == t.L] / (1 + n * sigma.2[time == t.L])
  a.U     <- n * sigma.2[time == t.U] / (1 + n * sigma.2[time == t.U])
  
  #=====> pull value from table <====#
  aU <- format(c(round(50 * a.U) / 50, 0.01))[1]
  aL <- format(c(round(50 * a.L) / 50, 0.01))[1]
  if (type[1] == 'hall'){
    if (confLevel == 0.90){
      hwk <- OIsurv:::hw.k10
    } else if(confLevel == 0.95){
      hwk <- OIsurv:::hw.k05
    } else if(confLevel == 0.99){
      hwk <- OIsurv:::hw.k01
    } else {
      stop("Only confidence levels of 0.90, 0.95 and 0.99\nare allowed.")
    }
    m.fact <- hwk[aU, aL] * (1 + n * sigma.2) / sqrt(n)
  } else if(type[1] == 'ep'){
    if (confLevel == 0.90){
      epc <- OIsurv:::ep.c10
    } else if(confLevel == 0.95){
      epc <- OIsurv:::ep.c05
    } else if(confLevel == 0.99){
      epc <- OIsurv:::ep.c01
    } else {
      stop("Only confidence levels of 0.90, 0.95 and 0.99\nare allowed.")
    }
    m.fact <- epc[aU, aL] * sqrt(sigma.2)
  } else {
    stop('The "type" variable is not recognized.')
  }
  CI <- matrix(NA, length(surv), 2)
  if(confType[1] == 'log-log'){
    theta  <- exp(m.fact / log(surv))
    CI[, 1] <- surv^(1 / theta)
    CI[, 2] <- surv^theta
  } else if(confType[1] == 'asin-sqrt'){
    temp   <- asin(sqrt(surv)) - 0.5 * m.fact * sqrt(surv / (1 - surv))
    lower  <- apply(cbind(rep(0, length(temp)), temp), 1, max)
    CI[,1] <- sin(lower)^2
    temp   <- asin(sqrt(surv)) + 0.5 * m.fact * sqrt(surv / (1 - surv))
    upper  <- apply(cbind(rep(pi / 2, length(temp)), temp), 1, min)
    CI[,2] <- sin(upper)^2
  } else {
    CI[,1] <- surv * (1 - m.fact)
    CI[,2] <- surv * (1 + m.fact)
  }
  CI[CI[,1] < 0, 1] <- 0
  CI[CI[,2] > 1, 2] <- 1
  if (tail(time, 1) != max(x[, 1])) {
    time <- c(time, max(x[,1]))
    CI   <- rbind(CI, CI[nrow(CI), ])
  }
  tR        <- list(time = time, lower = CI[, 1], upper = CI[, 2])
  class(tR) <- "confBands"
  return(tR)
}

