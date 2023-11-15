### Decentralized optimization using ADMM
## DC Power Flow Assumption, ignoring reactive power flow

## A trivial case
wtp <- 10
wts <- 0
flow_limit <- 10

# f(x) = -10 * x
# L: f(S) + .5 * \rho * (-x - w + 10 + u)^2
# w = 0, u = 0
# 10 = \rho * (x - 10 + w - u)
# x = 10 + 10 / \rho + u - w
# w = max(0, -x + 10) = 0
# u = u - x - w + 10 = -x + 10 = -10 / \rho
# ......

rho <- .1
w <- 0
u <- 0

x <- 10 + 10 / rho + u - w
w <- max(0, -x + 10)
u <- u - x - w + 10

# x > -10 -> x - w_neg = -10
# x < 10 -> -x - w = -10 -> x + w_pos = 10
# x + w = 10























# ## Price finding func for residual load curve
# price_find <- function(rl, residual_load_curve){
#   if(residual_load_curve[1] > rl){
#     nodal_price <- 1
#   }else if(residual_load_curve[num_price] < rl){
#     nodal_price <- num_price
#   }else{
#     price_low <- 1
#     price_high <- num_price
#     
#     while(price_high - price_low > 1){
#       price_mid <- as.integer(.5 * (price_low + price_high))
#       
#       if((residual_load_curve[price_low] - rl) * (residual_load_curve[price_mid] - rl) <= 0){
#         price_high <- price_mid
#       }else{
#         price_low <- price_mid
#       }      
#     }
#     nodal_price <- price[price_mid] 
#   }
#   return(nodal_price)
# }
# 
# ### Network Parameters
# num_node <- 10
# num_line <- 9
# z_img <- .001 
# Cond = diag(rep(1 / z_img, num_line)) 
# V_limit <- pi / 6
# I_limit <- 100
# 
# ## Topology of the network
# Topology = matrix(nrow = num_line, ncol = 2)
# colnames(Topology) = c("From", "To")
# 
# for(i in 1:nrow(Topology)){
#   Topology[i, ] = c(i, i+1)
# }
# # Topology[nrow(Topology), ] = c(nrow(Topology), 1)
# 
# # Topology in form of Node * Line
# NL = matrix(0, ncol = num_node, nrow = num_line)
# for(i in 1:num_line){
#   NL[i, Topology[i, 1]] = 1
#   NL[i, Topology[i, 2]] = -1
# }
# 
# ## Equivalent Admittance Matrix
# Y_eq <- t(NL) %*% Cond %*% NL
# Z_eq <- solve(Y_eq[2:num_node, 2:num_node])
# 
# ## ADMM
# # Augmented Lagrangian:
# # L(x, w, u) = f(x) + g(w) + .5 * \rho * ||A %*% x + B %*% w - u||^2
# # x = voltage, power source / sink, and line current
# # w = box constraint indicator for voltage and line current
# # u = dual variables of the constraints
# # f(x): only S has cost function associated
# # g(w) = well function for the box constraint
# num_box <- num_node + num_line
# num_eq <- num_box
# num_constraints <- num_eq + num_box
# num_variable <- 2 * num_node + num_line
# A <- matrix(0, nrow = num_constraints, ncol = num_variable)
# B <- matrix(0, nrow = num_constraints, ncol = num_box)
# 
# # Power Flow Equation
# # Y_eq %*% V - S = 0
# A[1:num_node, 1:num_node] <- Y_eq
# A[1:num_node, num_node + 1:num_node] <- diag(rep(-1, num_node))
# 
# # Ohm's Law
# # Y_L %*% NL %*% V - I = 0
# A[num_node + 1:num_line, 1:num_node] <- Cond %*% NL
# A[num_node + 1:num_line, 2 * num_node + 1:num_line] <- diag(rep(-1, num_line))
# 
# # Box Constraint for Voltage
# # V - w_v = 0
# A[num_eq + 1:num_node, 1:num_node] <- diag(rep(1, num_node))
# B[num_eq + 1:num_node, 1:num_node] <- diag(rep(-1, num_node))
# 
# # Box Constraint for Line Current
# # I - w_i = 0
# A[num_eq + num_node + 1:num_line, 2 * num_node + 1:num_line] <- diag(rep(1, num_line))
# B[num_eq + num_node + 1:num_line, num_node + 1:num_line] <- diag(rep(-1, num_line))
# 
# ### Market Parameters
# num_price <- 101
# price <- seq(0, 100, length.out = num_price)
# moc_supply <- matrix(0, nrow = num_price, ncol = num_node)
# rownames(moc_supply) <- price
# colnames(moc_supply) <- 1:num_node
# for(node_iter in 1:num_node){
#   moc_supply[(10 * (node_iter - 1) + 1):num_price, node_iter] <- 1
# }
# moc_demand <- 1 - moc_supply
# moc_supply[, 1] <-  moc_supply[, 1] * 10
# cumsum_supply <- apply(moc_supply, 2, cumsum)
# cumsum_demand <- apply(apply(apply(moc_demand, 2, rev), 2, cumsum), 2, rev)
# residual_load_curve <- cumsum_supply - cumsum_demand
# 
# ### Main optimization process
# rho <- 10
# nodal_price <- rep(0, num_node)
# S_prev <- rep(0, num_node - 1)
# S_now <- rep(0, num_node - 1)
# w <- rep(0, num_box)
# u <- rep(0, num_constraints)
# 
# # Market clearing at each node
# # A %*% dx = (A_v %*% Z_eq + A_S + A_I %*% (Cond %*% NL)[, 2:num_node] %*% Z_eq) %*% dS
# A_eq <- A[, 2:num_node] %*% Z_eq + A[, num_node + 2:num_node] + A[, 2 * num_node + 1:num_line] %*% (Cond %*% NL)[, 2:num_node] %*% Z_eq
# ## Using inverse matrix is dangerous!!
# 
# ## Main loop
# penalty <- rho * t(A_eq  %*% S_prev + B %*% w + u) %*% A_eq
# for(node_iter in 2:num_node){
#   S_green_func <- rep(0, num_node - 1)
#   S_green_func[node_iter - 1] <- 1
#   penalty_green_func <- rho * t(A_eq %*% S_green_func) %*% A_eq
#   nodal_price[node_iter] <- nodal_price[1] - penalty[node_iter - 1] + S_prev[node_iter - 1] * penalty_green_func[node_iter - 1]
#   
#   # Bisection method
#   rl_l <- residual_load_curve[1, node_iter]
#   penalty_l <- nodal_price[node_iter] - rl_l * penalty_green_func[node_iter - 1]
#   price_l <- price_find(rl_l, residual_load_curve[, node_iter])
#   rl_h <- residual_load_curve[num_price, node_iter]
#   penalty_h <- nodal_price[node_iter] - rl_h * penalty_green_func[node_iter - 1]
#   price_h <- price_find(rl_h, residual_load_curve[, node_iter])
#   
#   while(rl_h - rl_l > .01){
#     rl_m <- .5 * (rl_l + rl_h)
#     penalty_m <- nodal_price[node_iter] - rl_m * penalty_green_func[node_iter - 1]
#     price_m <- price_find(rl_m, residual_load_curve[, node_iter])
#     if((price_l - penalty_l) * (price_m - penalty_m) <= 0){
#       rl_h <- rl_m
#       penalty_h <- penalty_m
#       price_h <- price_m
#     }else{
#       rl_l <- rl_m
#       penalty_l <- penalty_m
#       price_l <- price_m    
#     }    
#   }
#   
#   nodal_price[node_iter] <- penalty_m
#   S_now[node_iter - 1] <- rl_m
# }
# 
# # Update Nodal Price at Reference
# nodal_price[1] <- price_find(-sum(S_now), residual_load_curve[, 1])
# 
# # Update w variables
# V_now <- c(0, Z_eq %*% S_now)
# I_now <- Cond %*% NL %*% V_now
# w <- c(V_now * (abs(V_now) < V_limit), I_now * (abs(I_now) < I_limit))
# 
# # Update u variables
# u <- u + A_eq %*% S_now + B %*% w
# S_prev <- S_now
# 
# S_now
# nodal_price