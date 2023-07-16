#' This function performs stratified Max-Combo test.
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @importFrom stats cov2cor
#' @param time a vector of event or censored times.
#' @param event a vector with entries 0 or 1, indicating event occurrence (1) or time being censored (0).
#' @param group a vector indicating treatment groups.
#' @param stratum a vector of stratification factor.
#' @param alternative choose from \code{"two.sided"}, \code{"less"} or \code{"greater"} to determine which type of tests to conduct and calculate the corresponding p-values.
#' @param rho a vector indicating different values of rho in Fleming-Harrington weight.
#' @param gamma a vector indicating different values of gamma in Fleming-Harrington weight.
#' @return A list with components:
#' \item{z.list}{a vector of z values calculated from all startified weighted log-rank tests under stratification method 1.}
#' \item{z.max}{the z value that is the furthest away from 0 under stratification method 1.}
#' \item{cov}{the covariance matrix of different startified weighted log-rank tests under stratification method 1.}
#' \item{pval}{p-value of desired alternative test under stratification method 1.}
#' \item{z.list2}{a vector of z values calculated from all startified weighted log-rank tests under stratification method 2.}
#' \item{z.max2}{the z value that is the furthest away from 0 under stratification method 2.}
#' \item{cov2}{the covariance matrix of different startified weighted log-rank tests under stratification method 2.}
#' \item{pval2}{p-value of desired alternative test under stratification method 2.}
#' \item{z.list3}{a vector of z values calculated from all startified weighted log-rank tests under stratification method 3.}
#' \item{z.max3}{the z value that is the furthest away from 0 under stratification method 3.}
#' \item{cov3}{the covariance matrix of different startified weighted log-rank tests under stratification method 3.}
#' \item{pval3}{p-value of desired alternative test under stratification method 3.}
#'
#' @references
#' Magirr, D. and Jiménez, J. (2023).
#' Stratified modestly-weighted log-rank tests in settings with an anticipated delayed separation of survival curves.
#' {\emph{Biometrical Journal.}, 2023;65:2200126} \doi{10.1002/bimj.202200126}
#'
#' @examples
#' data(sim_data)
#' ##load survival data
#' time <- sim_data$event_time
#' event <- sim_data$event_status
#' group <- sim_data$group
#' stratum <- sim_data$strata
#'
#' ##perform stratified max-combo test
#' SMCtest(time,event,group,stratum,alternative="two.sided",rho=c(0,1,1,0),gamma=c(0,0,1,1))
#'
#' @seealso \code{\link{WLRtest}}
#'
#' @export


SMCtest <- function(time,event,group,stratum,alternative=c("two.sided","less","greater"),rho=c(0,1,1,0),gamma=c(0,0,1,1)){
  u <- unique(stratum)
  ns <- length(u)
  nt <- length(rho)
  #sample size of each stratum
  sz <- as.numeric(table(stratum))

  #z values for all tests
  z.list <- rep(0,nt)
  z.list2 <- rep(0,nt)
  z.list3 <- rep(0,nt)

  #covariance components
  n1.cov <- list()
  n2.cov <- list()
  denominator.cov <- rep(0,nt)

  n3.cov2 <- list()
  denominator.cov2 <- rep(0,nt)

  denominator.cov3 <- rep(0,nt)

  #z statistics and cov components
  for(i in 1:nt){
    #single z value components
    numerator <- 0
    denominator <- 0
    numerator2 <- 0
    denominator2 <- 0
    numerator3 <- 0
    denominator3 <- 0
    #covariance numerator components
    n1 <- matrix(0,length(time),ns)
    n2 <- matrix(0,length(time),ns)
    n3 <- matrix(0,length(time),ns)

    for(j in 1:ns){
      index <- which(stratum==u[j])
      T0 <- WLRtest(time[index],event[index],group[index],rho=rho[i],gamma=gamma[i])
      numerator <- numerator + sum((T0$D$w)*(T0$D$diff1))
      denominator <- denominator + sum(((T0$D$w)^2)*(T0$D$var))
      n1[(1:dim(T0$D)[1]),j] <- (T0$D$w)*(T0$D$var)
      n2[(1:dim(T0$D)[1]),j] <- T0$D$w

      numerator2 <- numerator2 + sqrt(sum(T0$D$var))*(sum((T0$D$w)*(T0$D$diff1))/sqrt(sum(((T0$D$w)^2)*(T0$D$var))))
      denominator2 <- denominator2 + sum(T0$D$var)
      n3[(1:dim(T0$D)[1]),j] <- T0$D$var

      numerator3 <- numerator3 + sz[j]*sum((T0$D$w)*(T0$D$diff1))/sum(((T0$D$w)^2)*(T0$D$var))
      denominator3 <- denominator3 + (sz[j]^2)/sum(((T0$D$w)^2)*(T0$D$var))
    }
    z.list[i] <- numerator/sqrt(denominator)
    n1.cov[[i]] <- n1
    n2.cov[[i]] <- n2
    denominator.cov[i] <- denominator

    z.list2[i] <- numerator2/sqrt(denominator2)
    n3.cov2[[i]] <- n3
    denominator.cov2[i] <- denominator2

    z.list3[i] <- numerator3/sqrt(denominator3)
    denominator.cov3[i] <- denominator3
  }

  #covariance matrix construction
  cov <- diag(1,nt)
  for(ic in 1:(nt-1)){
    for(jc in 2:nt){
      temp <- n1.cov[[ic]]*n2.cov[[jc]]
      n0 <- sum(temp)
      d0 <- sqrt(denominator.cov[ic]*denominator.cov[jc])
      v0 <- n0/d0
      cov[ic,jc]<-cov[jc,ic]<-v0
    }
  }
  cor <- cov

  cov2 <- diag(1,nt)
  for(ic in 1:nt){
    for(jc in 2:nt){
      temp1 <- colSums(n3.cov2[[ic]])/sqrt(colSums(n1.cov[[ic]]*n2.cov[[ic]])*colSums(n1.cov[[jc]]*n2.cov[[jc]]))
      temp2 <- colSums(n1.cov[[ic]]*n2.cov[[jc]])
      n0 <- sum(temp1*temp2)
      d0 <- denominator.cov2[ic]
      v0 <- n0/d0
      cov2[ic,jc]<-cov2[jc,ic]<-v0
    }
  }
  cor2 <- cov2cor(cov2)

  cov3 <- diag(1,nt)
  for(ic in 1:nt){
    for(jc in 2:nt){
      temp1 <- (sz^2)/(colSums(n1.cov[[ic]]*n2.cov[[ic]])*colSums(n1.cov[[jc]]*n2.cov[[jc]]))
      temp2 <- colSums(n1.cov[[ic]]*n2.cov[[jc]])
      n0 <- sum(temp1*temp2)
      d0 <- sqrt(denominator.cov3[ic]*denominator.cov3[jc])
      v0 <- n0/d0
      cov3[ic,jc]<-cov3[jc,ic]<-v0
    }
  }
  cor3 <- cov2cor(cov3)
  #pval
  algorithm = GenzBretz(maxpts = 50000, abseps = 0.00001, releps = 0)

  if(alternative=="two.sided") {
    z.max <- max(abs(z.list))
    pval <- 1-pmvnorm(lower=rep(-z.max,nt),upper=rep(z.max,nt),corr=cov,algorithm=algorithm)[1]

    z.max2 <- max(abs(z.list2))
    pval2 <- 1-pmvnorm(lower=rep(-z.max2,nt),upper=rep(z.max2,nt),corr=cov2,algorithm=algorithm)[1]

    z.max3 <- max(abs(z.list3))
    pval3 <- 1-pmvnorm(lower=rep(-z.max3,nt),upper=rep(z.max3,nt),corr=cov3,algorithm=algorithm)[1]
  }

  if(alternative=="less") {
    z.max <- min(z.list)
    pval <- 1-pmvnorm(lower=rep(z.max,nt),upper=rep(Inf,nt),corr=cov,algorithm=algorithm)[1]

    z.max2 <- min(z.list2)
    pval2 <- 1-pmvnorm(lower=rep(z.max2,nt),upper=rep(Inf,nt),corr=cov2,algorithm=algorithm)[1]

    z.max3 <- min(z.list3)
    pval3 <- 1-pmvnorm(lower=rep(z.max3,nt),upper=rep(Inf,nt),corr=cov3,algorithm=algorithm)[1]
  }

  if(alternative=="greater") {
    z.max <- max(z.list)
    pval <- 1-pmvnorm(lower=rep(-Inf,nt),upper=rep(z.max,nt),corr=cov,algorithm=algorithm)[1]

    z.max2 <- max(z.list2)
    pval2 <- 1-pmvnorm(lower=rep(-Inf,nt),upper=rep(z.max2,nt),corr=cov2,algorithm=algorithm)[1]

    z.max3 <- max(z.list3)
    pval3 <- 1-pmvnorm(lower=rep(-Inf,nt),upper=rep(z.max3,nt),corr=cov3,algorithm=algorithm)[1]
  }

  return(list(z.list=z.list,z.max=z.max,cor=cor,pval=pval,z.list2=z.list2,z.max2=z.max2,cor2=cor2,pval2=pval2,z.list3=z.list3,z.max3=z.max3,cor3=cor3,pval3=pval3))
}
