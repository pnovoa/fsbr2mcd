roc_weights <- function(ncrit){
  
  m <- matrix(
    rep(1 / 1:ncrit, ncrit),
    nrow = ncrit,
    byrow = T
  )
  
  m[lower.tri(m)] <- 0.
  
  apply(m, MARGIN = 1, sum)/ncrit
  
}


rs_weights <- function(ncrit){
  
  sapply(1:ncrit, function(i) (2*(ncrit + 1 - i))/(ncrit*(ncrit + 1)))
}

rr_weights <- function(ncrit){
  
  sapply(1:ncrit, function(i) (1/i)/sum(1/(1:ncrit)))
}

ew_weights <- function(ncrit){
  rep(1/ncrit, ncrit)
}