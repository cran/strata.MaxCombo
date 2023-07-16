#' This function performs weighted log-rank test.
#' @importFrom stats aggregate pnorm
#' @param time time a vector of event or censored times.
#' @param event a vector with entries 0 or 1, indicating event occurrence (1) or time being censored (0).
#' @param group a vector indicating treatment groups.
#' @param rho a value of rho in Fleming-Harrington weight
#' @param gamma a value of gamma in Fleming-Harrington weight
#' @return A list with components:
#' \item{z}{z value of weighted log-rank test.}
#' \item{pvalue.two.sided}{p-value for two-sided test.}
#' \item{pvalue.less}{p-value for one-sided test (less).}
#' \item{pvalue.greater}{p-value for one-sided test (greater).}
#' \item{D}{a matrix containing components that can be used in stratified max-combo test.}
#'
#' @references
#' Ristl, R., Ballarini, N., Götte, H., Schüler, A., Posch, M. and König, F. (2021).
#' Delayed treatment effects, treatment switching and heterogeneous patient populations: How to design and analyze RCTs in oncology.
#' {\emph{Pharmaceutical Statistics.}, 20:129-145} \doi{10.1002/pst.2062}
#'
#' @examples
#' data(sim_data)
#' ##load survival data
#' time <- sim_data$event_time
#' event <- sim_data$event_status
#' group <- sim_data$group
#'
#' ##perform weighted log-rank test
#' WLRtest(time,event,group,rho=0,gamma=0)
#'
#' @seealso \code{\link{SMCtest}}
#'
#' @export

WLRtest <- function(time,event,group,rho=0,gamma=0) {
  n <- length(time)
  group <- factor(group)
  Ag <- aggregate(event,by=list(time=time,group=group),FUN=sum,drop=FALSE)
  Ag$x <- ifelse(is.na(Ag$x), 0, Ag$x)

  #distinct failure times | group 1 failures | group 2 failures
  tab <- data.frame(time=Ag$time[Ag$group==levels(group)[1]],event1=Ag$x[Ag$group==levels(group)[1]],event2=Ag$x[Ag$group==levels(group)[2]])

  #group 1 at risk | group 2 at risk
  tab$atrisk1 <- NA
  tab$atrisk2 <- NA
  tpoints <- dim(tab)[1]
  for(i in 1:tpoints) {
    tab$atrisk1[i] <- sum(time[group==levels(group)[1]]>=tab$time[i])
    tab$atrisk2[i] <- sum(time[group==levels(group)[2]]>=tab$time[i])
  }

  #total at risk | total failures
  tab$atrisk <- tab$atrisk1+tab$atrisk2
  tab$event <- tab$event1+tab$event2

  #unique event times
  D <- tab[tab$event>0,]
  #LR test components
  D$expected1 <- D$event*D$atrisk1/D$atrisk
  D$expected2 <- D$event*D$atrisk2/D$atrisk
  D$diff1 <- D$event1-D$expected1
  D$var <- D$atrisk1*D$atrisk2*D$event*(D$atrisk-D$event)/(D$atrisk^2*(D$atrisk-1))

  #FH weights
  D$S <- cumprod((D$atrisk-D$event)/D$atrisk)
  D$Sminus <- c(1,D$S[-length(D$S)])
  D$w <- D$Sminus^rho*(1-D$Sminus)^gamma

  #test statistics & pvalue
  D <- D[D$atrisk1>0 & D$atrisk2>0,]
  z <- sum(D$w*D$diff1)/sqrt(sum(D$w^2*D$var))
  pval.two.sided <- 2*(1-pnorm(abs(z)))
  pval.less <- pnorm(z)
  pval.greater <- 1-pnorm(z)

  return(list(z=z,pval.two.sided=pval.two.sided,pval.less=pval.less,pval.greater=pval.greater,D=D))
}
