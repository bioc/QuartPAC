matrix_combine <-
function(A, B) {
  result <- rbind(cbind(A, array(0, dim = c(nrow(A), ncol(B)))),
              cbind(array(0, dim = c(nrow(B), ncol(A))), B))
  
  return(result)
}
