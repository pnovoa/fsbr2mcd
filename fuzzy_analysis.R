library(FuzzyNumbers)


`[.myfuzzn` <- function(x, i, ...) structure(unclass(x)[i], class = "myfuzzn")

`==.myfuzzn` <- function(a, b) {identical(a, b)}

`>.myfuzzn` <- function(a, b) {
  
  a_ <- TriangularFuzzyNumber(a[[1]][1], a[[1]][2], a[[1]][3])
  b_ <- TriangularFuzzyNumber(b[[1]][1], b[[1]][2], b[[1]][3])
  
  a_ <- as.PiecewiseLinearFuzzyNumber(a_)
  b_ <- as.PiecewiseLinearFuzzyNumber(b_)
  
  return(
    possibilityStrictExceedance_comparison(a_, b_)
  )
  
}

is.na.myfuzzn <- function(a) FALSE


possibilityStrictExceedance_comparison <- function(fn_a, fn_b){
  poss_excee <- possibilityStrictExceedance(fn_a, fn_b)
  poss_under <- possibilityStrictUndervaluation(fn_a, fn_b)
  
  return(poss_excee > poss_under)
}

fuzzy_solution_sort <- function(fuzzy_num_def_matrix){
  
  
  # Defining the triangular fuzzy numbers
  lst_fns <- apply(fuzzy_num_def_matrix, MARGIN = 1, FUN = function(x) TriangularFuzzyNumber(x["LB"], x["CORE"], x["UB"]))
  lst_fns_ <- lapply(lst_fns, function(x) c(x["a1"], x["a2"], x["a4"]))
  
  lst_fns_cls <- structure(lst_fns_, class = "myfuzzn")
  
  # Sorting solutions
  sorted_fns <- sort(lst_fns_cls, decreasing = TRUE)
  
  return(
    names(sorted_fns)
  )
  
}