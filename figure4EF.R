# Clear the workspace
rm(list = ls())

# Load required package
library(ggplot2)

# Define shared plotting constants
BLACK <- "black"
sqrt3 <- sqrt(3)

# ============================================================
# Helper functions for equilateral mixture triangle geometry
# ============================================================

# Convert barycentric coordinates (x1, x2, x3) to 2D Cartesian
# coordinates for an equilateral triangle
bary_to_xy <- function(x1, x2, x3) {
  x <- x2 + 0.5 * x3
  y <- (sqrt3 / 2) * x3
  data.frame(x, y)
}

# Define triangle vertices and centroid
v1 <- bary_to_xy(1, 0, 0)
v2 <- bary_to_xy(0, 1, 0)
v3 <- bary_to_xy(0, 0, 1)
tri2d <- rbind(v1, v2, v3, v1)
center <- bary_to_xy(1 / 3, 1 / 3, 1 / 3)

# Compute the outward unit normal vector for a triangle edge
edge_outward_normal <- function(vA, vB) {
  e <- c(vB$x - vA$x, vB$y - vA$y)
  n1 <- c(-e[2], e[1])
  n1 <- n1 / sqrt(sum(n1^2))
  n2 <- -n1
  
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  to_center <- c(center$x - mid[1], center$y - mid[2])
  
  if (sum(n1 * to_center) < 0) n1 else n2
}

# Build internal grid lines for specified mixture levels
build_grid_2d_levels <- function(levels = c(0.2, 0.4, 0.6, 0.8)) {
  mk <- function(p1, p2) {
    data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
  }
  
  out <- list()
  for (t in levels) {
    out[[length(out) + 1]] <- mk(bary_to_xy(t, 0, 1 - t), bary_to_xy(t, 1 - t, 0))
    out[[length(out) + 1]] <- mk(bary_to_xy(0, t, 1 - t), bary_to_xy(1 - t, t, 0))
    out[[length(out) + 1]] <- mk(bary_to_xy(0, 1 - t, t), bary_to_xy(1 - t, 0, t))
  }
  
  do.call(rbind, out)
}

# Build percentage tick labels outside a triangle edge
build_edge_tick_labels <- function(vA, vB, values, offset = 0.06) {
  n <- edge_outward_normal(vA, vB)
  
  do.call(rbind, lapply(values, function(t) {
    P <- c((1 - t) * vA$x + t * vB$x, (1 - t) * vA$y + t * vB$y)
    data.frame(
      x = P[1] + n[1] * offset,
      y = P[2] + n[2] * offset,
      txt = paste0(round(100 * t), "%")
    )
  }))
}

# Build an axis label centered outside a triangle edge
build_axis_label <- function(sym, vA, vB, offset = 0.10) {
  n <- edge_outward_normal(vA, vB)
  mid <- c((vA$x + vB$x) / 2, (vA$y + vB$y) / 2)
  
  data.frame(
    x = mid[1] + n[1] * offset,
    y = mid[2] + n[2] * offset,
    txt = sym
  )
}

# ============================================================
# Construct the base EMT plot
# ============================================================

grid2d <- build_grid_2d_levels(c(0.2, 0.4, 0.6, 0.8))
ticks <- c(0.2, 0.4, 0.6, 0.8)

# Map axes to triangle edges
edges <- list(
  z1 = list(A = v1, B = v2),
  z2 = list(A = v1, B = v3),
  z3 = list(A = v2, B = v3)
)

# Build edge tick labels and axis labels
labs_z1 <- build_edge_tick_labels(edges$z1$A, edges$z1$B, ticks, offset = 0.06)
labs_z2 <- build_edge_tick_labels(edges$z2$A, edges$z2$B, ticks, offset = 0.06)
labs_z3 <- build_edge_tick_labels(edges$z3$A, edges$z3$B, ticks, offset = 0.06)

axis_labs <- rbind(
  build_axis_label("r[1]", edges$z1$A, edges$z1$B, offset = 0.10),
  build_axis_label("r[2]", edges$z2$A, edges$z2$B, offset = 0.10),
  build_axis_label("r[3]", edges$z3$A, edges$z3$B, offset = 0.10)
)

# Base equilateral mixture triangle plot
p_base <- ggplot() +
  geom_path(data = tri2d, aes(x, y), linewidth = 1, colour = BLACK) +
  geom_segment(
    data = grid2d,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    alpha = 0.8,
    colour = BLACK
  ) +
  geom_text(data = labs_z1, aes(x, y, label = txt), size = 3.8, colour = BLACK) +
  geom_text(data = labs_z2, aes(x, y, label = txt), size = 3.8, colour = BLACK) +
  geom_text(data = labs_z3, aes(x, y, label = txt), size = 3.8, colour = BLACK) +
  geom_text(data = axis_labs, aes(x, y, label = txt), parse = TRUE, size = 5, colour = BLACK) +
  coord_fixed(
    xlim = c(-0.15, 1.15),
    ylim = c(-0.15, sqrt3 / 2 + 0.15),
    expand = FALSE
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    plot.title   = element_blank()
  )

# ============================================================
# Simplex-lattice design
# ============================================================

# Build a {3, m} simplex-lattice design consisting of all mixtures
# with components equal to integer multiples of 1/m
build_simplex_lattice <- function(m = 4) {
  pts <- lapply(0:m, function(i) {
    lapply(0:(m - i), function(j) {
      k <- m - i - j
      x1 <- i / m
      x2 <- j / m
      x3 <- k / m
      data.frame(x1, x2, x3)
    })
  })
  
  pts <- do.call(rbind, unlist(pts, recursive = FALSE))
  rownames(pts) <- NULL
  
  # Classify points by whether they fall at vertices, edges, or the interior
  zeros <- (pts$x1 == 0) + (pts$x2 == 0) + (pts$x3 == 0)
  type <- ifelse(
    zeros == 2, "vertex",
    ifelse(zeros == 1, "edge", "interior")
  )
  pts$type <- factor(type, levels = c("vertex", "edge", "interior"))
  
  xy <- bary_to_xy(pts$x1, pts$x2, pts$x3)
  cbind(pts, xy)
}

# Generate the simplex-lattice design
m <- 4
lattice_pts <- build_simplex_lattice(m)

# Define point aesthetics by design-point type
palette_types <- c(
  vertex = "#1b9e77",
  edge = "#d95f02",
  interior = "#7570b3"
)

sizes_types <- c(
  vertex = 3.8,
  edge = 3.2,
  interior = 2.8
)

# Overlay simplex-lattice design points on the EMT
p_lattice <- p_base +
  geom_point(
    data = lattice_pts,
    aes(x = x, y = y, fill = type, size = type),
    shape = 21,
    colour = "black",
    alpha = 0.95,
    stroke = 0.9
  ) +
  scale_fill_manual(values = palette_types, name = "Design point") +
  scale_size_manual(values = sizes_types, guide = "none")

# ============================================================
# Simplex-centroid design
# ============================================================

# Build the 3-component simplex-centroid design consisting of
# vertices, edge midpoints, and the overall centroid
build_simplex_centroid_k3 <- function() {
  pts <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1 / 2, 1 / 2, 0),
    c(1 / 2, 0, 1 / 2),
    c(0, 1 / 2, 1 / 2),
    c(1 / 3, 1 / 3, 1 / 3)
  )
  
  colnames(pts) <- c("x1", "x2", "x3")
  
  type <- c(
    rep("vertex", 3),
    rep("edge", 3),
    "centroid"
  )
  
  df <- data.frame(
    pts,
    type = factor(type, levels = c("vertex", "edge", "centroid"))
  )
  
  xy <- bary_to_xy(df$x1, df$x2, df$x3)
  cbind(df, xy)
}

# Generate centroid design points
centroid_pts <- build_simplex_centroid_k3()

# Define point aesthetics
palette_centroid <- c(
  vertex = "#1b9e77",
  edge = "#d95f02",
  centroid = "#7570b3"
)

sizes_centroid <- c(
  vertex = 3.8,
  edge = 3.2,
  centroid = 4.2
)

# Overlay centroid design points on the EMT
p_centroid <- p_base +
  geom_point(
    data = centroid_pts,
    aes(x = x, y = y, fill = type, size = type),
    shape = 21,
    colour = "black",
    alpha = 0.98,
    stroke = 0.9
  ) +
  scale_fill_manual(
    values = palette_centroid,
    name = "Centroid design",
    breaks = c("vertex", "edge", "centroid"),
    labels = c(vertex = "vertex", edge = "edge", centroid = "centroid")
  ) +
  scale_size_manual(values = sizes_centroid, guide = "none")

# ============================================================
# Constrained continuous candidate set and D-optimal design
# ============================================================

# Maximum allowable values for each simplex component
max_z1 <- 0.85
max_z2 <- 0.55
max_z3 <- 0.50

# Sample points uniformly from the simplex using exponential draws
sample_simplex <- function(N) {
  e <- matrix(rexp(3 * N, rate = 1), ncol = 3)
  s <- rowSums(e)
  x <- e / s
  colnames(x) <- c("x1", "x2", "x3")
  as.data.frame(x)
}

# Generate points along a line of constant x1, x2, or x3
points_on_constant <- function(axis = c("x1", "x2", "x3"), t, n = 101) {
  axis <- match.arg(axis)
  u <- seq(0, 1 - t, length.out = n)
  
  if (axis == "x1") {
    df <- data.frame(x1 = t, x2 = u, x3 = 1 - t - u)
  } else if (axis == "x2") {
    df <- data.frame(x1 = u, x2 = t, x3 = 1 - t - u)
  } else {
    df <- data.frame(x1 = u, x2 = 1 - t - u, x3 = t)
  }
  
  df[df$x1 >= 0 & df$x2 >= 0 & df$x3 >= 0, , drop = FALSE]
}

# Build a continuous candidate set and augment it with points
# along the constraint boundaries
set.seed(1)
N_cand <- 3000
cand <- sample_simplex(N_cand)

aug <- rbind(
  points_on_constant("x1", max_z1, n = 151),
  points_on_constant("x2", max_z2, n = 151),
  points_on_constant("x3", max_z3, n = 151)
)

cand <- unique(rbind(cand, aug))

# Apply component-wise upper-bound constraints
cand <- subset(cand, x1 <= max_z1 & x2 <= max_z2 & x3 <= max_z3)
stopifnot(nrow(cand) > 0)

# Project candidate points into EMT coordinates
cand_xy <- bary_to_xy(cand$x1, cand$x2, cand$x3)
cand <- cbind(cand, cand_xy)

# ============================================================
# Scheffé model matrix and D-optimal row-exchange search
# ============================================================

# Construct the model matrix for a linear or quadratic Scheffé model
model_matrix_scheffe <- function(x1, x2, x3, order = c("linear", "quadratic")) {
  order <- match.arg(order)
  
  if (order == "linear") {
    X <- cbind(x1, x2, x3)
    colnames(X) <- c("x1", "x2", "x3")
  } else {
    X <- cbind(x1, x2, x3, x1 * x2, x1 * x3, x2 * x3)
    colnames(X) <- c("x1", "x2", "x3", "x1x2", "x1x3", "x2x3")
  }
  
  X
}

order <- "quadratic"
Xcand <- with(cand, model_matrix_scheffe(x1, x2, x3, order = order))
p <- ncol(Xcand)

# Choose the number of design runs
n_runs <- max(p, min(nrow(cand), ceiling(2.0 * p)))

# Compute the log determinant of X'X as the D-optimality criterion
logdet_xtx <- function(X) {
  XtX <- crossprod(X) + diag(1e-12, ncol(X))
  R <- chol(XtX)
  2 * sum(log(diag(R)))
}

# Perform greedy one-row exchanges to improve D-optimality
one_row_exchange <- function(Xcand, idx, max_iter = 50) {
  chosen <- sort(unique(idx))
  pool <- setdiff(seq_len(nrow(Xcand)), chosen)
  cur_ld <- logdet_xtx(Xcand[chosen, , drop = FALSE])
  
  improved <- TRUE
  it <- 0
  
  while (improved && it < max_iter) {
    improved <- FALSE
    it <- it + 1
    
    for (pos in seq_along(chosen)) {
      best_ld <- cur_ld
      best_swap <- NA_integer_
      
      for (cnd in pool) {
        trial <- chosen
        trial[pos] <- cnd
        new_ld <- logdet_xtx(Xcand[trial, , drop = FALSE])
        
        if (is.finite(new_ld) && new_ld > best_ld + 1e-10) {
          best_ld <- new_ld
          best_swap <- cnd
        }
      }
      
      if (!is.na(best_swap)) {
        pool <- sort(setdiff(union(pool, chosen[pos]), best_swap))
        chosen[pos] <- best_swap
        cur_ld <- best_ld
        improved <- TRUE
      }
    }
  }
  
  list(idx = sort(chosen), logdet = cur_ld)
}

# Search over multiple random starts for a good D-optimal design
set.seed(2)
best <- list(idx = integer(0), logdet = -Inf)
n_starts <- min(40, 2 * nrow(Xcand))

for (s in 1:n_starts) {
  start_idx <- sample(seq_len(nrow(Xcand)), size = n_runs, replace = FALSE)
  res <- one_row_exchange(Xcand, start_idx, max_iter = 60)
  if (res$logdet > best$logdet) best <- res
}

opt_idx <- best$idx
opt_design <- cand[opt_idx, , drop = FALSE]

# ============================================================
# Plot constraint boundaries and selected optimal design points
# ============================================================

# Build the dashed boundary segment for a constant x1 constraint
constraint_segment_x1 <- function(t) {
  p1 <- bary_to_xy(t, 0, 1 - t)
  p2 <- bary_to_xy(t, 1 - t, 0)
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}

# Build the dashed boundary segment for a constant x2 constraint
constraint_segment_x2 <- function(t) {
  p1 <- bary_to_xy(0, t, 1 - t)
  p2 <- bary_to_xy(1 - t, t, 0)
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}

# Build the dashed boundary segment for a constant x3 constraint
constraint_segment_x3 <- function(t) {
  p1 <- bary_to_xy(0, 1 - t, t)
  p2 <- bary_to_xy(1 - t, 0, t)
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}

seg_z1 <- constraint_segment_x1(max_z1)
seg_z2 <- constraint_segment_x2(max_z2)
seg_z3 <- constraint_segment_x3(max_z3)

# Overlay the constrained optimal design on the EMT
p_opt <- p_base +
  geom_segment(
    data = seg_z1,
    aes(x = x, y = y, xend = xend, yend = yend),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "black"
  ) +
  geom_segment(
    data = seg_z2,
    aes(x = x, y = y, xend = xend, yend = yend),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "black"
  ) +
  geom_segment(
    data = seg_z3,
    aes(x = x, y = y, xend = xend, yend = yend),
    linetype = "dashed",
    linewidth = 0.8,
    colour = "black"
  ) +
  geom_point(
    data = opt_design,
    aes(x = x, y = y),
    shape = 16,
    colour = "black",
    size = 3.6
  )

# ============================================================
# Display plots
# ============================================================

p_lattice
p_centroid
p_opt