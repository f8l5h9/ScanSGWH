#' @name nn_ellipse
#' @rdname nn_ellipse
#'
#' @title obtain circular or ellipetic windows for a set of coordianates
#'
#' @description Obtain circular or ellipetic windows for a set of coordianates
#'
#' @param    Cx      : Coordenate (Longitude). A vector Rx1
#' @param    Cy      : Coordenate (Latitude). A vector Rx1
#' @param    nn      : maximun (number of observations) size of the windows.
#' @param    p       : rotation angle ellipse.
#'
#' @return ventanas \code{x} and \code{y}.
#'
#' @examples
#'
#' R <- 1000
#' nn<-50
#' Cx <- runif(R,0,1)
#' Cy <- runif(R,0,1)
#' A <- nn_ellipse(Cx=Cx,Cy=Cy,nn=nn,p=30)
#'
#' @export
#'


nn_ellipse = function(Cx=Cx,Cy=Cy,nn=nn,p=p){
  # PURPOSE: Obtiene todas las elipses formadas por los nn vecinos mas proximos
  #          a cada punto.
  #          Para cada punto se obtienen todas las elipses con angulo de rotacion de 10:10:180
  #          y excentricidades predeterminadas en el vector S
  # ---------------------------------------------------
  #  USAGE: results = nn_ellipse(Cx,Cy,nn)
  #  where: Cx  = Coordenate X - Latitude
  #         Cy  = Coordenate Y - Longitude
  #         nn  = max number of neinbourg
  # ---------------------------------------------------
  #  RETURNS: results.nn(p,s).eq = matrix Rxnn donde cada fila tiene los
  #  vecinos elipticos mas proximos
  #  --------------------------------------------------
  # written by:
  # Fernando A. Lopez Hernandez, 10/08/2017
  # Facultad de C.C. de la Empresa
  # Dpto. Metodos Cuantitativos e Inform?ticos
  # Universidad Polit?cnica de Cartagena
  # E-mail: Fernando.Lopez@upct.es

  repmat = function(X,m,n){
    ##R equivalent of repmat (matlab)
    X <- as.matrix(X)
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)}
  ##
  R=length(Cx)
  XX <- (repmat(t(Cx),R,1)-repmat(Cx,1,R))
  YY=(repmat(t(Cy),R,1)-repmat(Cy,1,R));
  S<- c(1,2,3,4) # S <- c(1,1.5, 2, 3, 4, 5, 10)
  results <- list()
  ## get the circle
  sh <- 1
  a <- 1
  b=a*sh
  # p <- 30; # el angulo de rotacion. No es importante para obtener los circulos
  angulo <- (pi/180)*p
  X=XX*cos(angulo)+YY*sin(angulo)
  Y <- -XX*sin(angulo)+YY*cos(angulo);
  F <- (X/a)^2+(Y/b)^2
  e2 <- apply(F,2,order)
  F1 <- t(e2[t(seq(1,R)),])
  F1 <- F1[,1:nn]
  results$circles <- F1
  ## Ellipses
  F2 <- as.numeric(matrix(, nrow = 0, ncol = nn))
  for (s in 2:length(S)){
    sh=S[s];
    a=1;
    b=a*sh;
    for(i in seq(10, 180, by = 30)){
      angulo=(pi/180)*p
      X <- XX*cos(angulo)+YY*sin(angulo)
      Y <- -XX*sin(angulo)+YY*cos(angulo)
      F <- (X/a)^2+(Y/b)^2
      e2 <- apply(F,2,order)
      F1 <- t(e2[t(seq(1,R)),])
      F1 <- F1[,1:nn]
      F2 <-rbind(F2,F1)
      # results.nne(p/10,s-1).eq=int16(F1(:,1:nn));
    }
  }
  results$ellipses <-F2
  return(results)
}
