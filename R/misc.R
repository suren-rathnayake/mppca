is.defined <- function (sym, env) {
  sym <- deparse(substitute(sym))
  exists(sym, env)
}

chol.inv <- function(x, ...){
  C <- chol(x)
  inv_x <- chol2inv(C)
  return(inv_x)
}
