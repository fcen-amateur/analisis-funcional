# Estimador media
get_muhat <- function(X) {
    # X tiene N observaciones (columnas) p-variadas, y es de dimension p * n
    apply(X, ROWS, mean)
}
# Estimador covarianza
get_gammahat <- function(X, mu = NULL) {
    if (is.null(mu)) {
        muhat <- get_muhat(X)
    } else {
        muhat <- mu
    }
    N <- dim(X)[2]
    Xcent <- X - matrix(rep(muhat, N), ncol=N)
    ret <- (1 / N) * Xcent %*% t(Xcent)
    ret
} 

rmvnorm <- function(mu, sigma, n=1) {
    p <- length(mu)
    L <- chol(sigma)
    # M %*% t(M)
    Z <- matrix(rnorm(p * n),p, n) # 2 rows, N/2 columns
    ret <- t(L) %*% Z + matrix(rep(mu, n), nrow=p)
    return(ret)
}