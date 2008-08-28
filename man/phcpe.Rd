\name{phcpe}
\alias{phcpe}
\title{Gonen \& Heller Concordance Probability Estimate for the Cox Proportional Harzards model}
\description{
   A function to calculate Gonen \& Heller concordance probability
   estimate (CPE) for the Cox proportional harzards model.
}
\usage{
phcpe(coxfit, CPE.SE=FALSE,out.ties=FALSE)
}
\arguments{
  \item{coxfit}{a coxph or cph object}
  \item{CPE.SE}{a logical value indicating whether the standard error of
    the CPE should be calculated}
  \item{out.ties}{If out.ties is set to FALSE,pairs of observations tied
    on covariates will be used to calculate the CPE.  Otherwise, they will not be used.
  }
}
\value{
  \item{CPE}{Concordance Probability Estimate}
  \item{CPE.SE}{the Standard Error of the Concordance Probability Estimate}
}

\examples{

### create a simple data set for testing
set.seed(199)
nn <- 1000
time <- rexp(nn)
status <- sample(0:1, nn, replace=TRUE)
covar <- matrix(rnorm(3*nn), ncol=3)
survd <- data.frame(time, status, covar)
names(survd) <- c("time","status","x1","x2","x3")

library(survival)
coxph.fit <- coxph(Surv(time,status)~x1+x2+x3,data=survd)


### Calculate CPE only (needs much less time).
library(CPE)
phcpe(coxph.fit)
phcpe(coxph.fit,out.ties=TRUE)
#result is identical because the covariates are not tied #

### Calculate CPE and CPE.SE
phcpe(coxph.fit, CPE.SE=TRUE)
phcpe(coxph.fit, CPE.SE=TRUE,out.ties=TRUE)

#*** For unknown reason, 'coxph.fit' may need to be removed before running cph()***
rm(coxph.fit)

library(Design)
cph.fit <- cph(Surv(time, status)~x1+x2+x3, data=survd,method="breslow")

### Calculate CPE only (needs much less time).
phcpe(cph.fit)
phcpe(cph.fit,out.ties=TRUE)


### Calculate CPE and CPE.SE
phcpe(cph.fit, CPE.SE=TRUE)
phcpe(cph.fit, CPE.SE=TRUE,out.ties=TRUE)

}

\author{Qianxing Mo, Mithat Gonen and Glenn Heller; \email{moq@mskcc.org}}

\references{
  Mithat Gonen and Glenn Heller. (2005). Concordance probability and
  discriminatory power in proportional hazards regression.
  \emph{Biometrika}, 92, 4, pp.965-970
}
\keyword{survival}
