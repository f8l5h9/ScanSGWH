#' @name LMScanSGWH_SAR
#' @rdname LMScanSGWH_SAR
#'
#' @title LM Scan test of SGWH in SAR Models
#'
#' @description Get the LM Scan Test for SGWH
#'
#' @param    y       : Data vector Rx1
#' @param    X       : Data vector Rxk include the constant
#' @param    W       : W matrix RxR
#' @param    mC      : minimúm cluster size.
#' @param    NNN     : Set of windows (circular or elliptic).
#' @param    it2     : The number of Monte Carlo replications, e.g., 999, 9999
#'
#' @return The value of the test and Most Likelihood Cluster.
#'
#' @references
#' Chasco C., Le Gallo, J and López F.A. A Scan test for spatial groupwise heteroscedasticity in cross-sectional models with an application on houses prices in Madrid
#' \cr
#' @seealso
#' \code{\link{nn_ellipse}}
#' @examples
#'
#' R=100
#' y<-runif(R,0,1)
#' XX<-as.matrix(cbind(matrix(1,R,1),runif(R,0,1)))
#' cx <- runif(R, 0, 1)
#' cy <- runif(R, 0, 1)
#' co <- cbind(cx,cy)
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co,k=5,longlat=F)))
#' mC <- 10
#' NN <- ScanSGWH::nn_ellipse(Cx=co[,1],Cy=co[,2],nn=50,p=30)
#' ## Circles
#' NNN<-NN$circles
#' it2 <- 999
#' LMScan <- LMScanSGWH_SAR(y=y,X=XX,W=W,mC=mC,it2=it2,NNN=NNN)
#' ## Ellipses
#' NNN<-NN$ellipses
#' it2 <- 100
#' LMScan <- LMScanSGWH_SAR(y=y,X=XX,W=W,mC=mC,it2=it2,NNN=NNN)
#' @export
#'
#'

LMScanSGWH_SAR <- function(y=y,X=XX,W=W,mC=mC,it2=it2,NNN=NNN){
  y <- as.matrix(y)
  R <- length(y)
  Ws<-W/matrix(rowSums(W),nrow=R,ncol=R)
  sar<-spdep::lagsarlm(y ~ XX[,-1], listw=mat2listw(W, style="W"))
  A <- diag(R)-sar$rho*Ws
  Beta <- Matrix::solve(t(XX)%*%XX,t(XX)%*%A%*%y)
  u <- (A%*%y-XX%*%Beta)
  sigma2 <- as.numeric(Matrix::crossprod(u)/R)
  LM0<-LMScanSGWH_SAR_Chi(y=y,X=XX,W=W,rho=sar$rho,mC=mC,NNN=NNN)
  ## Bootstrap
  LMs <- matrix(0,ncol = it2,nrow = 1)
  for (i in 1:it2){
    Rp<-sample(R)
    yp<-y[Rp]
    XXp <- XX[Rp,]
    Wp <- Ws
    Wp <- Wp[,Rp]
    Wp <- Wp[Rp,]
    LMp<-LMScanSGWH_SAR_Chi(y=yp,X=XXp,W=Wp,rho=sar$rho,mC=mC,NNN=NNN)
    LMs[i]<-LMp$testLM
  }

  pvalue <- sum(LMs>LM0$testLM)/(it2+1)
  a<- which(LM0$LM==max(LM0$LM),arr.ind = TRUE)
  fil <- a[1,1]
  col <- a[1,2]
  ZLM <- NNN[fil,1:(col+mC-1)]
  sigma2in <-as.numeric(t(u[ZLM])%*%u[ZLM]/length(ZLM))
  v <- u
  v<-v[-ZLM]
  sigma2out=as.numeric(t(v)%*%v/(R-length(ZLM)))
  cat("LMScan-SGWH test... ",LM0$testLM,"\n")
  cat("pvalue... ",pvalue,"\n")
  cat("sigma2... ",sigma2,"\n")
  cat("sigma2in... ",sigma2in,"\n")
  cat("sigma2out... ",sigma2out,"\n")

  results <- list(testLM=max(LM0$LM),pvalue=pvalue,rho=sar$rho,sigma2=sigma2,sigma2in=sigma2in,sigma2out=sigma2,ZLM=ZLM,size=length(ZLM),Nrep=it2)

}
