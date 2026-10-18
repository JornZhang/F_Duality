set.seed(12345)

B <- 100000
n <- 8
p <- 3
q <- 2

## Projection matrix
Pmat <- function(A) {
  A %*% solve(crossprod(A)) %*% t(A)
}

############################################################
## Example 1: Theorem F-new
## fixed w, random Z
############################################################

w <- 1:n

F1 <- replicate(B, {
  
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  Z2 <- Z[, (q + 1):p, drop = FALSE]
  
  PZ  <- Pmat(Z)
  PZ2 <- Pmat(Z2)
  
  numerator <- drop(
    t(w) %*% (PZ - PZ2) %*% w
  ) / q
  
  denominator <- drop(
    t(w) %*% (diag(n) - PZ) %*% w
  ) / (n - p)
  
  numerator / denominator
})


############################################################
## Example 2: Theorem F-ols-new
## fixed epsilon, random X
############################################################

beta <- c(1, -1, 2)

C <- cbind(
  diag(q),
  matrix(0, nrow = q, ncol = p - q)
)

epsilon <- 1:8

F2 <- replicate(B, {
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  y <- drop(X %*% beta + epsilon)
  
  XtX_inv <- solve(crossprod(X))
  beta_hat <- drop(XtX_inv %*% crossprod(X, y))
  
  beta_diff <- beta_hat - beta
  
  numerator <- drop(
    t(beta_diff) %*%
      t(C) %*%
      solve(C %*% XtX_inv %*% t(C)) %*%
      C %*%
      beta_diff
  ) / q
  
  residual <- y - drop(X %*% beta_hat)
  
  denominator <- sum(residual^2) / (n - p)
  
  numerator / denominator
})


############################################################
## Q-Q plots
############################################################

prob <- seq(0.001, 0.999, length.out = 500)

F_theory <- qf(prob, df1 = q, df2 = n - p)

F1_emp <- quantile(F1, prob = prob, names = FALSE)
F2_emp <- quantile(F2, prob = prob, names = FALSE)

pdf("F_exact_simulation.pdf", width = 9, height = 4.2)

par(mfrow = c(1, 2))

lim1 <- range(c(F_theory, F1_emp))
plot(F_theory, F1_emp,
     xlim = lim1, ylim = lim1,
     xlab = "Theoretical quantiles of F(2,5)",
     ylab = "Empirical quantiles",
     main = expression("(a) Random " * Z * ", fixed " * w),
     pch = 16, cex = 0.6)
abline(0, 1, lty = 2)

lim2 <- range(c(F_theory, F2_emp))
plot(F_theory, F2_emp,
     xlim = lim2, ylim = lim2,
     xlab = "Theoretical quantiles of F(2,5)",
     ylab = "Empirical quantiles",
     main = expression("(b) Random " * X * ", fixed " * epsilon),
     pch = 16, cex = 0.6)
abline(0, 1, lty = 2)

dev.off()


############################################################
## Optional: empirical size at alpha = 0.05
############################################################

critical_value <- qf(0.95, df1 = q, df2 = n - p)

mean(F1 > critical_value)
mean(F2 > critical_value)
