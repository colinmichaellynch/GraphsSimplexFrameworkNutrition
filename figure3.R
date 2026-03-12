rm(list = ls())

library(ggplot2)
library(patchwork)
library(dtw)

## ------------------------------------------------------------
## CONFIGURATION (top-level knobs)
## ------------------------------------------------------------
TRUE_N_STEPS       <- 25      # separate cap for the TRUE trajectory
TRUE_BETA          <- 5       # softness of the TRUE trajectory (lower = more stochastic)
N_STEPS            <- 40      # default max time steps for other trajs
N_STEPS_MACRO      <- 25      # generous cap for macro rules
ALPHA_COMMON       <- 0.1     # EMA step size
N_RANDOM_SIMS      <- 80      # number of random sims to average
TOL_NEAR_TARGET    <- 1e-3    # stopping tolerance
set.seed(123)

## ------------------------------------------------------------
## 1. EMT helpers
## ------------------------------------------------------------
sqrt3 <- sqrt(3)

bary_to_xy <- function(z1, z2, z3){
  x <- z2 + 0.5 * z3
  y <- (sqrt(3)/2) * z3
  data.frame(x, y)
}

# triangle border
v1 <- bary_to_xy(1,0,0)
v2 <- bary_to_xy(0,1,0)
v3 <- bary_to_xy(0,0,1)
tri2d <- rbind(v1, v2, v3, v1)

## ------------------------------------------------------------
## 2. Foods (5) and target INSIDE convex hull, target higher
## ------------------------------------------------------------
# z1 = carb, z2 = protein, z3 = fat
foods <- rbind(
  c(0.80, 0.05, 0.15),  # carb-ish
  c(0.10, 0.75, 0.15),  # protein-ish
  c(0.15, 0.10, 0.75),  # high fat (top)
  c(0.55, 0.35, 0.10),  # mixed lower
  c(0.20, 0.45, 0.35)   # mid, helps keep hull full
)
colnames(foods) <- c("z1","z2","z3")
M <- nrow(foods)

foods_xy <- bary_to_xy(foods[,"z1"], foods[,"z2"], foods[,"z3"])
foods_df <- cbind(as.data.frame(foods), foods_xy)

# center of mass (blue triangle)
p <- rep(1/M, M)
zbar <- colSums(foods * p)
zbar_xy <- bary_to_xy(zbar[1], zbar[2], zbar[3])

# intake target (kept interior but pushed "up")
w_target <- c(0.05, 0.05, 0.70, 0.05, 0.15)
w_target <- w_target / sum(w_target)
target   <- as.numeric(w_target %*% foods)
target_xy <- bary_to_xy(target[1], target[2], target[3])

# convex hull
h_idx   <- chull(foods_df$x, foods_df$y); h_idx <- c(h_idx, h_idx[1])
hull_df <- foods_df[h_idx, c("x","y")]

## ------------------------------------------------------------
## 3. Grid
## ------------------------------------------------------------
levs <- c(0.2, 0.4, 0.6, 0.8)
mkseg <- function(p1, p2) data.frame(x=p1$x, y=p1$y, xend=p2$x, yend=p2$y)

grid_z1 <- do.call(rbind, lapply(levs, function(t)
  mkseg(bary_to_xy(t,0,1-t), bary_to_xy(t,1-t,0)) ))
grid_z2 <- do.call(rbind, lapply(levs, function(t)
  mkseg(bary_to_xy(0,t,1-t), bary_to_xy(1-t,t,0)) ))
grid_z3 <- do.call(rbind, lapply(levs, function(t)
  mkseg(bary_to_xy(0,1-t,t), bary_to_xy(1-t,0,t)) ))

## ------------------------------------------------------------
## 4. Trajectory helpers
## ------------------------------------------------------------
ema_step <- function(cur, food, alpha=ALPHA_COMMON){
  nxt <- (1 - alpha)*cur + alpha*food
  nxt <- nxt / sum(nxt)
  nxt
}

# deterministic greedy to target (your original optimal idea)
policy_greedy_to <- function(goal, alpha=ALPHA_COMMON){
  function(cur, foods){
    dists <- apply(foods, 1, function(fj){
      cand <- (1 - alpha)*cur + alpha*fj
      sum((cand - goal)^2)
    })
    which.min(dists)
  }
}

# fixed-food policy
policy_fixed_food <- function(i){
  function(cur, foods){
    i
  }
}

# random policy
policy_random <- function(){
  function(cur, foods){
    sample(seq_len(nrow(foods)), 1)
  }
}

# macro-prioritizing policy (your previous version)
policy_macro <- function(target, priority_idx, alpha=ALPHA_COMMON, high_w = 0.95){
  W <- rep((1 - high_w)/2, 3)
  W[priority_idx] <- high_w
  function(cur, foods){
    scores <- apply(foods, 1, function(fj){
      cand <- (1 - alpha)*cur + alpha*fj
      deficit <- cand - target
      sum( W * deficit^2 )
    })
    which.min(scores)
  }
}

# ---- manuscript-style softmax policy (Section S1.3) ----
policy_softmax_to <- function(goal, beta = 4, w_food = NULL) {
  function(cur, foods) {
    v_cur  <- goal - cur
    norm_c <- sqrt(sum(v_cur^2))
    M      <- nrow(foods)
    scores <- numeric(M)
    
    for (i in seq_len(M)) {
      v_food  <- goal - foods[i, ]
      norm_f  <- sqrt(sum(v_food^2))
      if (norm_c == 0 || norm_f == 0) {
        align_score <- 0
      } else {
        cos_theta   <- sum(v_cur * v_food) / (norm_c * norm_f)
        cos_theta   <- max(min(cos_theta, 1), -1)
        # manuscript: d_{ti} = cos(theta - pi) = -cos(theta)
        align_score <- -cos_theta
      }
      scores[i] <- align_score
    }
    
    if (is.null(w_food)) {
      w_food <- rep(1, M)
    }
    
    probs <- w_food * exp(beta * scores)
    probs <- probs / sum(probs)
    
    sample.int(M, 1, prob = probs)
  }
}

# simulate until near goal
simulate_traj_to_goal <- function(start_bary,
                                  policy_fun,
                                  goal_bary,
                                  alpha = ALPHA_COMMON,
                                  max_steps = N_STEPS,
                                  tol = TOL_NEAR_TARGET,
                                  force_final = FALSE){
  cur <- start_bary
  out <- matrix(NA_real_, nrow = max_steps + 1, ncol = 3)
  out[1,] <- cur
  t <- 1
  dist2 <- sum((cur - goal_bary)^2)
  
  while (t <= max_steps && dist2 >= tol^2) {
    j <- policy_fun(cur, foods)
    cur <- ema_step(cur, foods[j,], alpha)
    out[t+1,] <- cur
    dist2 <- sum((cur - goal_bary)^2)
    t <- t + 1
  }
  
  out <- out[1:t, , drop = FALSE]
  
  if (force_final && dist2 >= tol^2) {
    out[t, ] <- goal_bary
  }
  
  out <- as.data.frame(out); names(out) <- c("z1","z2","z3")
  out
}

to_xy_traj <- function(tr_bary){
  xy <- bary_to_xy(tr_bary$z1, tr_bary$z2, tr_bary$z3)
  cbind(xy, t = seq_len(nrow(tr_bary)))
}

# trim a traj once it reaches target
trim_to_target <- function(df, target, tol = TOL_NEAR_TARGET){
  d2 <- (df$z1 - target[1])^2 + (df$z2 - target[2])^2 + (df$z3 - target[3])^2
  idx <- which(d2 < tol^2)
  if (length(idx) > 0) {
    df[1:idx[1], , drop = FALSE]
  } else {
    df
  }
}

## ------------------------------------------------------------
## 5. Start point and TRUE probabilistic trajectory
## ------------------------------------------------------------
start_bary <- c(0.975, 0.015, 0.010)

true_bary <- simulate_traj_to_goal(
  start_bary,
  policy_softmax_to(target, beta = TRUE_BETA),  # <-- new stochastic policy
  goal_bary = target,
  alpha     = ALPHA_COMMON,
  max_steps = TRUE_N_STEPS,
  tol       = TOL_NEAR_TARGET
)

true_bary <- trim_to_target(true_bary, target, TOL_NEAR_TARGET)
true_traj <- to_xy_traj(true_bary)

## ------------------------------------------------------------
## 6. Theoretical trajectories (with stop rules)
## ------------------------------------------------------------

## random average (goal = center of mass)
rand_paths <- array(NA_real_, dim = c(N_STEPS+1, 3, N_RANDOM_SIMS))
for(s in 1:N_RANDOM_SIMS){
  tr_b <- simulate_traj_to_goal(
    start_bary,
    policy_random(),
    goal_bary = zbar,
    alpha = ALPHA_COMMON,
    max_steps = N_STEPS,
    tol = TOL_NEAR_TARGET,
    force_final = FALSE
  )
  
  # pad to N_STEPS+1 by repeating last row
  pad_len <- (N_STEPS + 1) - nrow(tr_b)
  if (pad_len > 0) {
    last_row <- as.numeric(tr_b[nrow(tr_b), ])
    pad_mat  <- matrix(rep(last_row, pad_len),
                       ncol = 3, byrow = TRUE)
    colnames(pad_mat) <- names(tr_b)
    tr_b <- rbind(tr_b, as.data.frame(pad_mat))
  }
  
  rand_paths[,,s] <- as.matrix(tr_b)
}

rand_mean_bary <- apply(rand_paths, c(1,2), mean)
rand_mean_bary <- as.data.frame(rand_mean_bary); names(rand_mean_bary) <- c("z1","z2","z3")
rand_traj <- to_xy_traj(rand_mean_bary)

# closest distance (deterministic)
opt_traj <- to_xy_traj(
  simulate_traj_to_goal(
    start_bary,
    policy_greedy_to(target, alpha=ALPHA_COMMON),
    goal_bary = target,
    alpha = ALPHA_COMMON,
    max_steps = N_STEPS,
    tol = TOL_NEAR_TARGET,
    force_final = FALSE
  )
)

# 5 food trajectories: each ends at its food item
food_trajs <- lapply(1:M, function(i){
  to_xy_traj(
    simulate_traj_to_goal(
      start_bary,
      policy_fixed_food(i),
      goal_bary = foods[i,],
      alpha = ALPHA_COMMON,
      max_steps = N_STEPS,
      tol = TOL_NEAR_TARGET,
      force_final = FALSE
    )
  )
})

# 3 macronutrient rules (all end at intake target) with different routes
macro_trajs <- list(
  "Macronutrient: Carb" = to_xy_traj(
    simulate_traj_to_goal(
      start_bary,
      policy_macro(target, priority_idx = 1, alpha=ALPHA_COMMON, high_w = 0.95),
      goal_bary = target,
      alpha = ALPHA_COMMON,
      max_steps = N_STEPS_MACRO,
      tol = TOL_NEAR_TARGET,
      force_final = FALSE
    )
  ),
  "Macronutrient: Protein" = to_xy_traj(
    simulate_traj_to_goal(
      start_bary,
      policy_macro(target, priority_idx = 2, alpha=ALPHA_COMMON, high_w = 0.7),
      goal_bary = target,
      alpha = ALPHA_COMMON,
      max_steps = N_STEPS_MACRO,
      tol = TOL_NEAR_TARGET,
      force_final = FALSE
    )
  ),
  "Macronutrient: Fat" = to_xy_traj(
    simulate_traj_to_goal(
      start_bary,
      policy_macro(target, priority_idx = 3, alpha=ALPHA_COMMON, high_w = 0.7),
      goal_bary = target,
      alpha = ALPHA_COMMON,
      max_steps = N_STEPS_MACRO,
      tol = TOL_NEAR_TARGET,
      force_final = FALSE
    )
  )
)

## ------------------------------------------------------------
## 7. Distance / Similarity (DTW + simple path)
## ------------------------------------------------------------
path_dist <- function(obs_df, theo_df){
  L <- min(nrow(obs_df), nrow(theo_df))
  sum(sqrt((obs_df$x[1:L] - theo_df$x[1:L])^2 +
             (obs_df$y[1:L] - theo_df$y[1:L])^2))
}

dtw_dist_2d <- function(obs_df, theo_df){
  alignment <- dtw::dtw(
    as.matrix(obs_df[,c("x","y")]),
    as.matrix(theo_df[,c("x","y")]),
    dist.method = "Euclidean"
  )
  alignment$distance
}

sim_list <- list()

# random
sim_list[["Random Choice"]] <- list(
  euc = path_dist(true_traj, rand_traj),
  dtw = dtw_dist_2d(true_traj, rand_traj)
)

# closest distance
sim_list[["Closest Distance"]] <- list(
  euc = path_dist(true_traj, opt_traj),
  dtw = dtw_dist_2d(true_traj, opt_traj)
)

# food-based
for(i in 1:M){
  nm <- paste0("Food ", i)
  sim_list[[nm]] <- list(
    euc = path_dist(true_traj, food_trajs[[i]]),
    dtw = dtw_dist_2d(true_traj, food_trajs[[i]])
  )
}

# macronutrient
for(nm in names(macro_trajs)){
  sim_list[[nm]] <- list(
    euc = path_dist(true_traj, macro_trajs[[nm]]),
    dtw = dtw_dist_2d(true_traj, macro_trajs[[nm]])
  )
}

sim_df <- data.frame(
  name = names(sim_list),
  euc_dist = sapply(sim_list, function(x) x$euc),
  dtw_dist = sapply(sim_list, function(x) x$dtw),
  stringsAsFactors = FALSE
)

sim_df$type <- "Other"
sim_df$type[sim_df$name == "Random Choice"] <- "Random"
sim_df$type[sim_df$name == "Closest Distance"] <- "Closest Distance"
sim_df$type[grepl("^Food", sim_df$name)] <- "Food Choice"
sim_df$type[grepl("^Macronutrient", sim_df$name)] <- "Macronutrient"

alpha_s <- 4
sim_df$dtw_sim  <- exp(-alpha_s * sim_df$dtw_dist)
sim_df$dtw_sim_norm <- sim_df$dtw_sim / sum(sim_df$dtw_sim)
best_idx <- which.max(sim_df$dtw_sim_norm)
sim_df$highlight <- ifelse(seq_len(nrow(sim_df)) == best_idx, "best", "other")

## ------------------------------------------------------------
## 8. EMT plot (lines only, gray slightly transparent)
## ------------------------------------------------------------
theo_col <- "#555555"

p_emt <- ggplot() +
  geom_segment(data = grid_z1, aes(x=x, y=y, xend=xend, yend=yend),
               colour="grey90", linewidth=0.5) +
  geom_segment(data = grid_z2, aes(x=x, y=y, xend=xend, yend=yend),
               colour="grey90", linewidth=0.5) +
  geom_segment(data = grid_z3, aes(x=x, y=y, xend=xend, yend=yend),
               colour="grey90", linewidth=0.5) +
  geom_polygon(data = hull_df, aes(x, y), fill="grey85", alpha=0.45, colour=NA) +
  geom_path(data = hull_df, aes(x, y), colour="grey40", linewidth=0.8) +
  geom_path(data = rand_traj, aes(x, y), colour=theo_col, linewidth=1, alpha=0.45) +
  geom_path(data = opt_traj,  aes(x, y), colour=theo_col, linewidth=1, alpha=0.45)

for(i in 1:M){
  p_emt <- p_emt +
    geom_path(data = food_trajs[[i]], aes(x, y),
              colour=theo_col, linewidth=1, alpha=0.45)
}

p_emt <- p_emt +
  geom_path(data = macro_trajs[["Macronutrient: Carb"]],
            aes(x, y), colour=theo_col, linewidth=1, alpha=0.45) +
  geom_path(data = macro_trajs[["Macronutrient: Protein"]],
            aes(x, y), colour=theo_col, linewidth=1, alpha=0.45) +
  geom_path(data = macro_trajs[["Macronutrient: Fat"]],
            aes(x, y), colour=theo_col, linewidth=1, alpha=0.45) +
  geom_path(data = true_traj, aes(x, y), colour="#FF7F0E", linewidth=1.5) +
  geom_point(data = foods_df, aes(x, y), colour="black", size=3) +
  geom_point(aes(x=zbar_xy$x, y=zbar_xy$y), shape=24, fill="#0072B2", colour="black", size=4) +
  geom_point(aes(x=target_xy$x, y=target_xy$y), shape=21, fill="goldenrod2", colour="black", size=3.5) +
  geom_path(data = tri2d, aes(x, y), linewidth = 1.1, colour = "black") +
  coord_fixed(xlim=c(-0.15, 1.15), ylim=c(-0.15, sqrt3/2 + 0.15), expand=FALSE) +
  theme_void(base_size = 13)

# triangle vertices again for axis labels
v1 <- bary_to_xy(1,0,0)
v2 <- bary_to_xy(0,1,0)
v3 <- bary_to_xy(0,0,1)

# center of triangle
center <- bary_to_xy(1/3, 1/3, 1/3)

edge_outward_normal <- function(vA, vB, center){
  e  <- c(vB$x - vA$x, vB$y - vA$y)
  n1 <- c(-e[2], e[1]); n1 <- n1 / sqrt(sum(n1^2))
  n2 <- -n1
  mid <- c((vA$x + vB$x)/2, (vA$y + vB$y)/2)
  to_center <- c(center$x - mid[1], center$y - mid[2])
  if (sum(n1 * to_center) < 0) n1 else n2
}

build_axis_label <- function(label, vA, vB, center, offset=0.10){
  n   <- edge_outward_normal(vA, vB, center)
  mid <- c((vA$x + vB$x)/2, (vA$y + vB$y)/2)
  data.frame(
    x   = mid[1] + n[1]*offset,
    y   = mid[2] + n[2]*offset,
    txt = label
  )
}

axis_lab_right  <- build_axis_label("Fat %",     v2, v3, center, offset = 0.10)
axis_lab_left   <- build_axis_label("Protein %", v3, v1, center, offset = 0.10)
axis_lab_bottom <- build_axis_label("Carb %",    v1, v2, center, offset = 0.10)

p_emt <- p_emt +
  geom_text(data = axis_lab_right,  aes(x, y, label = txt),
            fontface = "bold", size = 5) +
  geom_text(data = axis_lab_left,   aes(x, y, label = txt),
            fontface = "bold", size = 5) +
  geom_text(data = axis_lab_bottom, aes(x, y, label = txt),
            fontface = "bold", size = 5)

## ------------------------------------------------------------
## 9. Plot 2: bar plot of DTW-based similarity
## ------------------------------------------------------------
facet_labs <- c(
  "Random"           = "Random",
  "Closest Distance" = "Closest\nDistance",
  "Food Choice"      = "Food Choice",
  "Macronutrient"    = "Macronutrient",
  "Other"            = "Other"
)

p_bar_dtw <- ggplot(sim_df, aes(x = name, y = dtw_sim_norm, fill = highlight)) +
  geom_col() +
  facet_grid(
    . ~ type,
    scales = "free_x",
    space  = "free_x",
    labeller = as_labeller(facet_labs)
  ) +
  scale_fill_manual(values = c(best = "goldenrod2", other = "grey70")) +
  labs(x = NULL, y = "Normalized DTW Similarity") +
  theme_bw(base_size = 16) +
  theme(
    legend.position   = "none",
    strip.background  = element_rect(fill = "white"),
    axis.text.x       = element_text(angle = 35, hjust = 1)
  )

## ------------------------------------------------------------
## 10. Plot 3: distance comparison
## ------------------------------------------------------------
metrics_long <- rbind(
  data.frame(
    name  = sim_df$name,
    type  = sim_df$type,
    metric = "Euclidean Dist",
    value = sim_df$euc_dist,
    stringsAsFactors = FALSE
  ),
  data.frame(
    name  = sim_df$name,
    type  = sim_df$type,
    metric = "DTW Dist",
    value = sim_df$dtw_dist,
    stringsAsFactors = FALSE
  )
)

p_metrics <- ggplot(metrics_long, aes(x = name, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  facet_grid(
    . ~ type,
    scales = "free_x",
    space  = "free_x",
    labeller = as_labeller(facet_labs)
  ) +
  labs(x = NULL, y = "Distance", fill = "Metric") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "white"),
    legend.position  = "top",
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.box.background = element_rect(fill = NA, colour = NA)
  )


## show
p_emt
p_bar_dtw
p_metrics

