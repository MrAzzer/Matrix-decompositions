Decomposition=function(A, b){
  #======================cholesky=================================
  if (ncol(A)==nrow(A) && all(A==t(A)) && all(eigen(A)$values>0)){
    print("Cholesky decomposition")
    L=chol(A)
    y=solve(t(L), b)
    x=solve(L, y)
    return (x)
  }
  else {
    warning("Can't be decomposed with cholesky")
  }
  
  #=======================LU============================#
  if (ncol(A)==nrow(A) && (det(A)!=0) && qr(A)$rank==ncol(A)){
    
    print("LU decomposition")
    
    LU=lu(A)
    L=expand(LU)$L
    U=expand(LU)$U
    y=solve(L, b)
    x=solve(U, y)
    return (x)
  }
  else {
    warning("Can't be decomposed with LU")
  }
  
  #====================QR===========================#
  if (qr(A)$rank==ncol(A)){
    print("QR decomposition")
    QR=qr(A)
    Q=qr.Q(QR)
    R=qr.R(QR)
    x=solve(R, t(Q) %*% b)
    return (x)
  }
  else {
    warning("Can't be decomposed with QR")
  }
  
  #=====================SVD=====================#
  print("SVD decomposition")
  SVD=svd(A)
  d=SVD$d
  v=SVD$v
  u=SVD$u
  d_inv=diag(1/d)
  x=v %*% d_inv %*% u %*% b
  return (x)
  }


A=matrix(c(4, 2, 2, 2, 3, 1, 2, 1, 3), nrow = 3, byrow = TRUE)
b=c(8, 7, 8)
Decomposition(A, b)

