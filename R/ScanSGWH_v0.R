#' @name ScanSGWH_v0
#' @rdname ScanSGWH_v0
#'
#' @title Scan Spatial Groupwise Heteroscedasticity (SGWH) for Normal process
#'
#' @description Get the Scan Test for SGWH
#'
#' @param    y       : Data vector Rx1
#' @param    C       : Coordenate vector Rx2.
#' @param    nv      : maximun number of neighbors inside cluster.
#' @param    k       : Looking for only high incidence cluster (k=0); only low (k=1); both (k=2).
#' @param    Nrep    : The number of Monte Carlo replications, e.g., 999, 9999
#'
#' @return The value of the test and Most Likelihood Cluster.
#'
#' @references Chasco C., Le Gallo, J and López F.A. A Scan test for spatial groupwise heteroscedasticity in cross-sectional models with an application on houses prices in Madrid
#'
#' @examples
#'
#' y <- runif(100,0,1)
#' nv <- 25
#' k <- 2
#' Nrep <- 10
#' C <- cbind(runif(100,0,1),runif(100,0,1))
#' LRScan <- ScanSGWH_v0(y=y,C=C,nv=nv,k=k,Nrep=Nrep)
#' @export
#'

ScanSGWH_v0<-function(y=y,C=C,NNN=NNN,nv=nv,k=k,Nrep=Nrep){
  #### ---------------------------------------------------------------------
  #  Inputs
  #    Nrep    : The number of Monte Carlo replications, e.g., 999, 9999
  #    C       : Coordenate vector Rx2.
  #    y       : Data vector Rx1
  #    nv      : maximun number of neighbors inside cluster.
  #    e       : Cluster type (e=0) circular (e=1) elliptic.
  #    k       : Looking for only high incidence cluster (k=0); only low (k=1); both (k=2).
  #    info    : Optional argument
  #
  #  Outputs
  #    Lik     : Lik value of Scan_sigma
  #    Likm    : Lik value of Scan_mu. As in Papers in Regional Science
  #    Elem    : Observation in the windows Z for Scan_sigma
  #    size    : Number element inside Z
  #    Coor    : Coordenate of the center of Z
  #    Lik     : Lik value of Scan_sigma
  #    mElem   : Observation in the windows Z for Scan_mu
  #    size    : Number element inside Z

  #    p_value : Monte Carlo simulated (adjusted) p-value
  #    ScanSGWH: Value Scan statisitic
  #    size    : Number of element inside of most likely cluster (MLC)
  #    Coor    : Coordenates of the center of MLC
  #    angle   : Angle of ellipse of MLC
  #    shape   : Shape of ellipse of MLC
  #    Elem    : Element (Id's) inside MLC
  #    MLboots : Value of Lik in each boots

R <- length(y)
# NNN <- knearneigh(C, k=nv, longlat = NULL, RANN=TRUE)
# NNN <- NNN$nn
yNNN <- y[NNN[,1:nv]]
dim(yNNN) = dim(NNN)
Ra=dim(NNN)[1]
# stdY <- apply (yNNN,1,sd)
stdY <- sd(y)*sqrt((R-1)/R)
mY=mean(y);
X=sum(y);
LL0s <- c(0,0,0)
LL0m <- c(0,0,0)

for (i in 1:(nv-1)){
nz=i+1
yNN=yNNN[,1:nz]
dim(yNN)=c(Ra,nz)
xz <- rowSums(yNN)
landaz <- (X-xz)/(R-nz)
muz <- xz/nz
sz<-(rowSums(yNN^2)-2*xz*muz+nz*muz^2+(sum(y^2)-rowSums(yNN^2))-2*(X-xz)*landaz+(R-nz)*landaz^2)/R
Lm=R*log(stdY)+sum((y-mY)^2)/(2*stdY^2)-(R/2)-R*log(sqrt(sz))
stdZ=sqrt(rowSums((yNN-mY)^2)/nz)
stdZC=(sqrt((R*stdY^2-nz*stdZ^2)/(R-nz)))
Ls=-nz*log(stdZ)-(R-nz)*log(stdZC)
if (k==0) { # Look for only high incidence cluster
Ls=Ls*(stdZ>stdZC)
Lm=Lm*(muz>landaz)
} else if (k==1) {# Look for only low incidence cluster
Ls=Ls*(stdZ<stdZC)
Lm=Lm*(muz<landaz)
} else if (k==2) {# Look for cluster of high or low incidence
Ls=Ls
Lm=Lm
}
# Para obtener max Lik de cluster secundarios
LLs <- cbind(Ls,seq(1:Ra),nz*rep(1,Ra))
LLs <- rbind(LL0s,LLs)
LLs <- LLs[order(LLs[,1],LLs[,2],decreasing=F),]
LL0s <- LLs[(dim(LLs)[1]-min(R-1,dim(LLs)[1]-1)):dim(LLs)[1],] ## LL0=[valor de Lik ; fila de NNN en la que esta el cluster; n de elemetos de la fila que debo tomar]
# Para obtener max Lik de cluster secundarios
LLm <- cbind(Lm,seq(1:Ra),nz*rep(1,Ra))
LLm <- rbind(LL0m,LLm)
LLm <- LLs[order(LLm[,1],LLm[,2],decreasing=F),]
LL0m <- LLs[(dim(LLm)[1]-min(R-1,dim(LLm)[1]-1)):dim(LLm)[1],] ## LL0=[valor de Lik ; fila de NNN en la que esta el cluster; n de elemetos de la fila que debo tomar]
}
rm(stdZ,stdZC)
results <- list()
results$Lik <- LL0s[dim(LL0s)[1],1]  # valor del Log de Lik
results$Likm <- LL0m[dim(LL0m)[1],1] # valor del Log de Lik
results$Elem <- NNN[LL0s[dim(LL0s)[1],2],(1:LL0s[dim(LL0s)[1],3])] # Elementos que forman el cluster mas probable
results$size <- length(results$Elem) # tama?o del cluster
results$Coor <- C[results$Elem[1],] # coordenadas del centro del cluster
results$mElem <- NNN[LL0s[dim(LL0m)[1],2],(1:LL0s[dim(LL0m)[1],3])] # Elementos que forman el cluster mas probable
results$msize <- length(results$mElem);# tama?o del cluster
results$mCoor <- C[results$mElem[1],] # coordenadas del centro del cluster

#results$shape <- NNN[LL0s[dim(LL0s)[1],2],end] # excentricidad de la elipse
# results$angle <- NNN[LL0s[dim(LL0s)[1],2],(dim(LL0s)[1]-1)] # angulo de la elipse

results$m <- mY;

# results$m_in <- mean(y(results.Elem));
# results$m_out=(R*mY-results.size*results.m_in)/(R-results.size) ;
results$s=stdY;
results$s_in=sd(y[results$Elem],1);
# results$s_out=sqrt((R*stdY^2-results.size*results.s_in^2)/(R-results.size));
results$nw=Ra; # Total Number of differents windows
MLf <- results$Lik;
MLmf <- results$Likm;

## Permutational bootstrapping (Dwass 1957)
MLboots <- rep(0,Nrep)
mMLboots <- rep(0,Nrep)
paso=R*log(stdY)+sum((y-mY)^2)/(2*stdY^2)-(R/2)
I <- t(apply(matrix(1,Ra,nv),1,cumsum))
RstdY2 <- R*stdY^2

for (f in 1:Nrep){
Y=y[sample(R)]
YNNN=Y[NNN[,1:nv]]
dim(YNNN)<-c(R,nv)
cs <- t(apply(YNNN^2,1,cumsum))
cs2 <- t(apply(YNNN,1,cumsum))
cm2 <- cs2/I;
cm <- (X-cs2)/(R-I);
cm22 <- cm2^2;
LLm=R*log(((cs-2*cs2*cm2+I*cm22+(sum(Y^2)-cs)-2*(X-cs2)*cm+(R-I)*cm^2)/R))
LLm=paso-LLm/2

STDZ <- sqrt(t(apply((YNNN-mY)^2,1,cumsum))/I)
STDZC <- sqrt(((RstdY2-I*STDZ^2)/(R-I)))
STDZ <- STDZ[,(2:dim(STDZ)[2])]
STDZC=STDZC[,(2:dim(STDZC)[2])]
LLs=-I[,2:dim(I)[2]]*log(STDZ)-(R-I[,2:dim(I)[2]])*log(STDZC)

if (k==0){  # cluster de alta incidencia
LLs=LLs*(STDZ>STDZC)
LLm=LLm*(cm2>cm)
} else if (k==1) { # cluster de baja incidencia
LLs=LLs*(STDZ<STDZC);
LLm=LLm*(cm2<cm);
}
MLboots[f]=max(LLs)
mMLboots[f]=max(LLm)
}

## Salida de Resultados
MLboots <- MLboots[1:Nrep]
mMLboots <- mMLboots[1:Nrep]
results$pvalue <- sum(MLboots>MLf)/(Nrep+1)
results$mpvalue <- sum(mMLboots>MLmf)/(Nrep+1)
results$Nrep <- Nrep;
results$Coor <- C;
return(results)
}





