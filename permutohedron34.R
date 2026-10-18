pdf("permutohedron34.pdf", height = 5.5, width = 10)

get_perms <- function(vec) {
  if (length(vec) == 1) return(matrix(vec, 1, 1))
  do.call(rbind, lapply(1:length(vec), function(i) {
    cbind(vec[i], get_perms(vec[-i]))
  }))
}

rotate_3d <- function(pts, theta = 35, phi = 25) {
  t <- theta * pi / 180
  p <- phi * pi / 180
  Rz <- matrix(c(cos(t), -sin(t), 0,
                 sin(t),  cos(t), 0,
                 0,       0,      1), 3, 3, byrow = TRUE)
  Rx <- matrix(c(1, 0,       0,
                 0, cos(p), -sin(p),
                 0, sin(p),  cos(p)), 3, 3, byrow = TRUE)
  return(pts %*% Rz %*% Rx)
}

# Set plot margins to accommodate axes and 1x2 layout
par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))

# ==============================================================================
# 1. Plot for n = 3 (Circle S^1)
# ==============================================================================
perms3 <- get_perms(c(1, 2, 3))
centered3 <- sweep(perms3, 2, colMeans(perms3), "-")
basis3 <- rbind(c(1, -1, 0) / sqrt(2), c(1, 1, -2) / sqrt(6))

pts2d_3 <- centered3 %*% t(basis3)
radius3 <- sqrt(2)

theta_seq <- seq(0, 2 * pi, length.out = 300)
plot(radius3 * cos(theta_seq), radius3 * sin(theta_seq), 
     type = "l", col = "#FFA500", lwd = 2, lty = 2, asp = 1,
     xlim = c(-2.3, 2.3), ylim = c(-2.3, 2.3),
     xlab = "", ylab = "",
     # main = expression("Permutohedron (" * n == 3 * ") on " * S^1))
     main = expression("Permutohedron (" * n == 3 * ")"))
grid()

# Draw edges
d_mat3 <- as.matrix(dist(pts2d_3))
for (i in 1:6) {
  for (j in i:6) {
    if (abs(d_mat3[i, j] - sqrt(2)) < 1e-5) {
      lines(pts2d_3[c(i, j), 1], pts2d_3[c(i, j), 2], col = "#4169E1", lwd = 2)
    }
  }
}

# Draw points and labels
points(pts2d_3[, 1], pts2d_3[, 2], col = "#DC143C", pch = 19, cex = 1.3)
text(pts2d_3[, 1] * 1.25, pts2d_3[, 2] * 1.25, apply(perms3, 1, paste, collapse = ""), font = 2, cex = 0.85)

# ==============================================================================
# 2. Plot for n = 4 (Sphere S^2)
# ==============================================================================
perms4 <- get_perms(c(1, 2, 3, 4))
centered4 <- sweep(perms4, 2, colMeans(perms4), "-")
basis4 <- rbind(c(1, -1, 0, 0) / sqrt(2), c(1, 1, -2, 0) / sqrt(6), c(1, 1, 1, -3) / sqrt(12))

pts3d_4 <- centered4 %*% t(basis4)
pts_rot4 <- rotate_3d(pts3d_4)
pts2d_4 <- pts_rot4[, 1:2]
z_depth <- pts_rot4[, 3]
radius4 <- sqrt(5)

plot(radius4 * cos(theta_seq), radius4 * sin(theta_seq), 
     type = "l", col = "#FFA500", lwd = 2, lty = 2, asp = 1,
     xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2),
     xlab = "", ylab = "",
     # main = expression("Permutohedron (" * n == 4 * ") on " * S^2)
     main = expression("Permutohedron (" * n == 4 * ")"))
grid()

d_mat4 <- as.matrix(dist(pts3d_4))

# Back edges (dashed)
for (i in 1:24) {
  for (j in i:24) {
    if (abs(d_mat4[i, j] - sqrt(2)) < 1e-5 && ((z_depth[i] + z_depth[j]) / 2 < 0)) {
      lines(pts2d_4[c(i, j), 1], pts2d_4[c(i, j), 2], col = "#B0B0B0", lty = 2, lwd = 1.2)
    }
  }
}

# Front edges (solid)
for (i in 1:24) {
  for (j in i:24) {
    if (abs(d_mat4[i, j] - sqrt(2)) < 1e-5 && ((z_depth[i] + z_depth[j]) / 2 >= 0)) {
      lines(pts2d_4[c(i, j), 1], pts2d_4[c(i, j), 2], col = "#4169E1", lty = 1, lwd = 1.8)
    }
  }
}

is_front <- z_depth >= 0
labels4 <- apply(perms4, 1, paste, collapse = "")

# Draw back points and labels
points(pts2d_4[!is_front, 1], pts2d_4[!is_front, 2], col = "#D3D3D3", pch = 19, cex = 0.9)
text(pts2d_4[!is_front, 1] * 1.12, pts2d_4[!is_front, 2] * 1.12, labels4[!is_front], col = "#A0A0A0", cex = 0.65)

# Draw front points and labels
points(pts2d_4[is_front, 1], pts2d_4[is_front, 2], col = "#DC143C", pch = 19, cex = 1.3)
text(pts2d_4[is_front, 1] * 1.15, pts2d_4[is_front, 2] * 1.15, labels4[is_front], col = "black", font = 2, cex = 0.75)

dev.off()
