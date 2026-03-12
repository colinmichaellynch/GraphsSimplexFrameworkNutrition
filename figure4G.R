rm(list = ls())

#setwd("~/Documents/EMT VS RMT")

## ============================================================
## Packages
## ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(AlgDesign)
  library(fields)    # for Tps smoothing
  library(scales)
})

BLACK <- "black"
sqrt3 <- sqrt(3)

## ============================================================
## 1. EMT helpers (same coordinate system as yours)
## ============================================================
bary_to_xy <- function(r1, r2, r3){
  x <- r2 + 0.5 * r3
  y <- (sqrt(3)/2) * r3
  data.frame(x, y)
}

xy_to_bary <- function(x, y){
  r3 <- (2 / sqrt(3)) * y
  r2 <- x - 0.5 * r3
  r1 <- 1 - r2 - r3
  cbind(r1 = r1, r2 = r2, r3 = r3)
}

inside_simplex <- function(x, y, tol = 1e-9){
  b <- xy_to_bary(x, y)
  apply(b >= -tol, 1, all)
}

# triangle scaffold
v1 <- bary_to_xy(1,0,0)
v2 <- bary_to_xy(0,1,0)
v3 <- bary_to_xy(0,0,1)
tri2d  <- rbind(v1, v2, v3, v1)
center <- bary_to_xy(1/3,1/3,1/3)

# grid lines
build_grid_2d_levels <- function(levels = c(0.2,0.4,0.6,0.8)){
  mk <- function(p1, p2) data.frame(x=p1$x, y=p1$y, xend=p2$x, yend=p2$y)
  out <- list()
  for(t in levels){
    out[[length(out)+1]] <- mk(bary_to_xy(t,0,1-t), bary_to_xy(t,1-t,0)) # r1
    out[[length(out)+1]] <- mk(bary_to_xy(0,t,1-t), bary_to_xy(1-t,t,0)) # r2
    out[[length(out)+1]] <- mk(bary_to_xy(0,1-t,t), bary_to_xy(1-t,0,t)) # r3
  }
  do.call(rbind, out)
}
grid2d <- build_grid_2d_levels(c(0.2,0.4,0.6,0.8))

# outward normal for labels
edge_outward_normal <- function(vA, vB){
  e  <- c(vB$x - vA$x, vB$y - vA$y)
  n1 <- c(-e[2], e[1]); n1 <- n1 / sqrt(sum(n1^2))
  n2 <- -n1
  mid <- c((vA$x + vB$x)/2, (vA$y + vB$y)/2)
  to_center <- c(center$x - mid[1], center$y - mid[2])
  if (sum(n1 * to_center) < 0) n1 else n2
}

build_edge_tick_labels <- function(vA, vB, values, offset = 0.06){
  n <- edge_outward_normal(vA, vB)
  do.call(rbind, lapply(values, function(t){
    P <- c((1-t)*vA$x + t*vB$x, (1-t)*vA$y + t*vB$y)
    data.frame(x = P[1] + n[1]*offset,
               y = P[2] + n[2]*offset,
               txt = paste0(round(100*t), "%"))
  }))
}

build_axis_label <- function(sym, vA, vB, offset = 0.10){
  n <- edge_outward_normal(vA, vB)
  mid <- c((vA$x + vB$x)/2, (vA$y + vB$y)/2)
  data.frame(x = mid[1] + n[1]*offset,
             y = mid[2] + n[2]*offset,
             txt = sym)
}

ticks <- c(0.2,0.4,0.6,0.8)
labs_r1 <- build_edge_tick_labels(v1, v2, ticks, offset = 0.06)
labs_r2 <- build_edge_tick_labels(v1, v3, ticks, offset = 0.06)
labs_r3 <- build_edge_tick_labels(v2, v3, ticks, offset = 0.06)
axis_labs <- rbind(
  build_axis_label("r[1]", v1, v2, offset = 0.10),
  build_axis_label("r[2]", v1, v3, offset = 0.10),
  build_axis_label("r[3]", v2, v3, offset = 0.10)
)

base_xlim <- c(-0.15, 1.15)
base_ylim <- c(-0.15, sqrt3/2 + 0.15)

## ============================================================
## 2. constrained candidate set
## ============================================================
max_r1 <- 0.80
max_r2 <- 0.95
max_r3 <- 0.70

line_on_const <- function(axis = c("r1","r2","r3"), t, n = 81) {
  axis <- match.arg(axis)
  u <- seq(0, 1 - t, length.out = n)
  if (axis == "r1") {
    df <- data.frame(r1 = t, r2 = u, r3 = 1 - t - u)
  } else if (axis == "r2") {
    df <- data.frame(r1 = u, r2 = t, r3 = 1 - t - u)
  } else {
    df <- data.frame(r1 = u, r2 = 1 - t - u, r3 = t)
  }
  df[df$r1 >= 0 & df$r2 >= 0 & df$r3 >= 0, , drop = FALSE]
}

simplex_grid <- function(step = 0.03) {
  pts <- list()
  for (r1 in seq(0, 1, by = step)) {
    for (r2 in seq(0, 1 - r1, by = step)) {
      r3 <- 1 - r1 - r2
      pts[[length(pts) + 1]] <- c(r1 = r1, r2 = r2, r3 = r3)
    }
  }
  as.data.frame(do.call(rbind, pts))
}

cand <- simplex_grid(step = 0.03)
cand <- rbind(
  cand,
  line_on_const("r1", max_r1, n = 121),
  line_on_const("r2", max_r2, n = 121),
  line_on_const("r3", max_r3, n = 121)
)
cand <- subset(cand, r1 <= max_r1 & r2 <= max_r2 & r3 <= max_r3)
cand <- unique(cand)
cand_all <- cand

# dashed constraint segments
constraint_segment_r1 <- function(t) {
  p1 <- bary_to_xy(t, 0, 1 - t)
  p2 <- bary_to_xy(t, 1 - t, 0)
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}
constraint_segment_r2 <- function(t) {
  p1 <- bary_to_xy(0, t, 1 - t)
  p2 <- bary_to_xy(1 - t, t, 0)
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}
constraint_segment_r3 <- function(t) {
  p1 <- bary_to_xy(0, 1 - t, t)
  p2 <- bary_to_xy(1 - t, 0, t)
  data.frame(x = p1$x, y = p1$y, xend = p2$x, yend = p2$y)
}
seg_r1 <- constraint_segment_r1(max_r1)
seg_r2 <- constraint_segment_r2(max_r2)
seg_r3 <- constraint_segment_r3(max_r3)

## ============================================================
## 3. A-efficient design (AlgDesign)
## ============================================================
form_mix <- ~ r1 + r2 + r3 + r1:r2 + r1:r3 + r2:r3 - 1
nTrials  <- 18

set.seed(123)
des <- optFederov(
  form_mix,
  data     = cand[, c("r1","r2","r3")],
  nTrials  = nTrials,
  nRepeats = 30
)
design_pts <- des$design

## ============================================================
## 4. simulate true growth, fit model
## ============================================================
true_growth <- function(r1, r2, r3){
  peak <- exp(-((r1 - 0.45)^2/0.010 + (r2 - 0.35)^2/0.012 + (r3 - 0.20)^2/0.008))
  50 + 35 * peak + 2*r1 + 1.5*r2
}

set.seed(77)
sigma_eps <- 1.5

design_data <- design_pts
design_data$growth <- with(design_data,
                           true_growth(r1, r2, r3) + rnorm(nrow(design_data), 0, sigma_eps))

fit <- lm(growth ~ r1 + r2 + r3 + r1:r2 + r1:r3 + r2:r3 - 1,
          data = design_data)

## ============================================================
## 5. Predict on rectangular grid, but keep only simplex + constraints
## ============================================================
# rectangular grid like in your TPS script
make_rect_grid <- function(nx = 320, ny = 320){
  xlo <- -0.01; xhi <- 1.01
  ylo <- -0.01; yhi <- sqrt3/2 + 0.01
  xo <- seq(xlo, xhi, length.out = nx)
  yo <- seq(ylo, yhi, length.out = ny)
  expand.grid(x = xo, y = yo)
}
RECT <- make_rect_grid(320, 320)

# mask to simplex
RECT$in_simplex <- inside_simplex(RECT$x, RECT$y, tol = 1e-9)

# convert to barycentric to enforce constraints
bary_rect <- xy_to_bary(RECT$x, RECT$y)
RECT$r1 <- bary_rect[, "r1"]
RECT$r2 <- bary_rect[, "r2"]
RECT$r3 <- bary_rect[, "r3"]

RECT$in_constraints <- with(RECT,
                            in_simplex &
                              r1 >= 0 & r2 >= 0 & r3 >= 0 &
                              r1 <= max_r1 & r2 <= max_r2 & r3 <= max_r3
)

# predict model at design xy and TPS-smooth (style like your script)
design_xy <- bary_to_xy(design_pts$r1, design_pts$r2, design_pts$r3)
df_smooth <- data.frame(
  x = design_xy$x,
  y = design_xy$y,
  z = predict(fit, newdata = design_pts)
)

# de-duplicate
dup <- !duplicated(df_smooth[, c("x","y")])
df_smooth <- df_smooth[dup, , drop = FALSE]

# TPS fit (simple)
tps_fit <- fields::Tps(df_smooth[, c("x","y")], df_smooth$z)

RECT$pred <- NA_real_
RECT$pred[RECT$in_constraints] <- predict(tps_fit, RECT[RECT$in_constraints, c("x","y")])

# get predicted optimum from candidate set
cand_all$pred <- predict(fit, newdata = cand_all)
cand_all_xy   <- bary_to_xy(cand_all$r1, cand_all$r2, cand_all$r3)
cand_all      <- cbind(cand_all, cand_all_xy)
cand_all      <- cand_all[order(-cand_all$pred), ]
pred_opt      <- cand_all[1, ]

## ============================================================
## 6. 5 replicates at predicted optimum
## ============================================================
set.seed(101)
true_at_opt <- true_growth(pred_opt$r1, pred_opt$r2, pred_opt$r3)
replicates  <- rnorm(5, mean = true_at_opt, sd = sigma_eps)
rep_mean    <- mean(replicates)
rep_sd      <- sd(replicates)

cat("Predicted response at optimum:", round(pred_opt$pred, 2), "\n")
cat("Replicate mean:", round(rep_mean, 2), "SD:", round(rep_sd, 2), "\n")
cat("Bias (rep_mean - predicted):", round(rep_mean - pred_opt$pred, 2), "\n")

rep_xy <- bary_to_xy(pred_opt$r1, pred_opt$r2, pred_opt$r3)
rep_data <- data.frame(
  x = rep_xy$x + rnorm(5, 0, 0.004),
  y = rep_xy$y + rnorm(5, 0, 0.004),
  growth = replicates
)

## ============================================================
## 7. final plot: TPS surface (new palette), NO legend, green star
## ============================================================
# choose a different HCL palette from your earlier script
n_bins <- 22
pal <- hcl.colors(n_bins, palette = "Mint", rev = FALSE)

p_case <- ggplot() +
  # filled contours from TPS-smoothed surface, only where feasible
  stat_contour_filled(
    data = RECT[RECT$in_constraints, ],
    aes(x = x, y = y, z = pred, fill = after_stat(level)),
    bins = n_bins,
    alpha = 0.78
  ) +
  scale_fill_manual(values = pal, drop = FALSE) +
  # triangle outline and internal grid
  geom_path(data = tri2d, aes(x, y), linewidth = 1, colour = BLACK) +
  geom_segment(data = grid2d, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.45, alpha = 0.8, colour = BLACK) +
  # dashed constraints
  geom_segment(data = seg_r1, aes(x = x, y = y, xend = xend, yend = yend),
               linetype = "dashed", linewidth = 0.8, colour = "black") +
  geom_segment(data = seg_r2, aes(x = x, y = y, xend = xend, yend = yend),
               linetype = "dashed", linewidth = 0.8, colour = "black") +
  geom_segment(data = seg_r3, aes(x = x, y = y, xend = xend, yend = yend),
               linetype = "dashed", linewidth = 0.8, colour = "black") +
  # tick labels and axis labels
  geom_text(data = labs_r1, aes(x, y, label = txt), size = 3.5, colour = BLACK) +
  geom_text(data = labs_r2, aes(x, y, label = txt), size = 3.5, colour = BLACK) +
  geom_text(data = labs_r3, aes(x, y, label = txt), size = 3.5, colour = BLACK) +
  geom_text(data = axis_labs, aes(x, y, label = txt), parse = TRUE, size = 5, colour = BLACK) +
  # black design points
  geom_point(data = cbind(design_pts, bary_to_xy(design_pts$r1, design_pts$r2, design_pts$r3)),
             aes(x = x, y = y),
             shape = 16, colour = "black", size = 3.1) +
  # predicted optimum as green star
  geom_point(data = data.frame(x = pred_opt$x, y = pred_opt$y),
             aes(x = x, y = y),
             shape = 8, colour = "forestgreen", size = 4.7, stroke = 1.1) +
  # replicates at optimum
  #geom_point(data = rep_data, aes(x = x, y = y),
  #           shape = 21, fill = "forestgreen", colour = "black", size = 2.8) +
  coord_fixed(xlim = base_xlim, ylim = base_ylim, expand = FALSE) +
  theme_bw(base_size = 14) +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )

print(p_case)

## predicted value from the model at the optimum
pred_val <- pred_opt$pred    # numeric(1)

## actual 5 measurements at that same mix
y <- replicates              # numeric(5)

## basic summary
n          <- length(y)
obs_mean   <- mean(y)
obs_sd     <- sd(y)
bias       <- obs_mean - pred_val
abs_err    <- abs(y - pred_val)
sq_err     <- (y - pred_val)^2

stats <- list(
  n               = n,
  predicted_value = pred_val,
  observed_mean   = obs_mean,
  observed_sd     = obs_sd,
  bias_mean_minus_pred = bias,
  mean_abs_error  = mean(abs_err),
  rmse            = sqrt(mean(sq_err)),
  min_obs         = min(y),
  max_obs         = max(y),
  obs_quantiles   = quantile(y, probs = c(0.25, 0.5, 0.75))
)

str(stats)


