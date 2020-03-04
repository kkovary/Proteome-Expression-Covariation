library(mvtnorm)

fun = function(cor_matrix, list_distributions, L)
{
  n = length(list_distributions)
  # Correlated Gaussian variables
  Gauss = rmvnorm(n=L, mean = rep(0,n), sig=cor_matrix)
  # convert them to uniform distribution.
  Unif = pnorm(Gauss) 
  # Convert them to whatever I want
  vars = sapply(1:n, FUN = function(i) list_distributions[[i]](Unif[,i]))
  return(vars)
}

# Create a matrix of values with a given R, mean, and standard deviation
corProts <- function(nProt, nCell, R, mean, sd) {
  
  cor_matrix =  matrix(data = R,
                       ncol = nProt,
                       nrow = nProt)
  for (i in 1:nProt) {
    cor_matrix[i, i] = 1
  }
  
  list_distributions <- list()
  for (i in 1:nProt) {
    list_distributions <-
      c(list_distributions, function(nCell)
        round(qnorm(nCell, mean, sd)))
  }
  
  mat = fun(cor_matrix, list_distributions, nCell)
  return(mat)
}
