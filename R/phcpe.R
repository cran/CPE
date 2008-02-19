#Qianxing Mo, moq@mskcc.org
#Department of Epidemiology and Biostatistics
#Memorial Sloan-Kettering Cancer Center, NY 10021

#This R function is used to call the c function that computes the concordance probability for
#the cox proportional hazard model (see Gonen and Heller, Biometrika(2005), 92, 4, pp. 965-970)
#The c function implements a fast algorithm to compute the probability primarily based on
#Dr. Ennapadam Venkatraman's R codes

#The input for phcpe fuction must be a 'coxph' or 'cph' object

phcpe <- function(coxfit, CPE.SE=TRUE) {
  if(class(coxfit)[1] != "coxph" && class(coxfit)[1] != "cph"){
    stop("Error! Input must be a coxph or cph object")
  }
  
  row <- as.integer(sum(coxfit$n))
  col <- as.integer(length(coxfit$coefficients))
  design <- model.matrix(coxfit)
  design <- as.double(as.vector(t(design[,-1])))
  xbeta <- as.double(as.vector(coxfit$linear.predictors))
  varbeta <- as.double(as.vector(t(coxfit$var)))
  bandwidth <- as.double(0.5*sd(coxfit$linear.predictors)*(row^(-1/3)))

  if(CPE.SE==TRUE){
    if(row >= 3000) {
      cat(c("It may take about n*n minutes to calculate 10000*n rows of data.\n"))
    }

    res <- .C("coxcpe",row, col, bandwidth, xbeta, design, varbeta, out=as.double(rep(0, 3)),PACKAGE="CPE")
    return(list(CPE = res$out[1], CPE.SE = res$out[3]))
  }else {
    res <- .C("coxcpeOnly",row,xbeta,out=as.double(0), PACKAGE="CPE")
    return(list(CPE=res$out))
  }
}

