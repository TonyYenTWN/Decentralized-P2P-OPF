### Decentralized optimization using ADMM, simplified variables
## DC Power Flow Assumption, ignoring reactive power flow
## Convergence guaranteed for sequential method; parallel method needs proof and validation of stabilization method

## merit order curve of residual load creation
moc_create <- function(bid_demand, bid_supply, M){
  moc_rl <- matrix(ncol = 2, nrow = 2 * (nrow(bid_demand) + nrow(bid_supply) - 3))
  price_demand_ID <- 1
  price_supply_ID <- 1
  
  price_rl_ID <- 1
  moc_rl[price_rl_ID, 1] <- -M
  moc_rl[price_rl_ID, 2] <- -sum(bid_demand[, 2])
  while(price_rl_ID < nrow(moc_rl) - 1){
    price_rl_ID <- price_rl_ID + 1
    if(bid_demand[price_demand_ID + 1, 1] <= bid_supply[price_supply_ID + 1, 1]){
      price_demand_ID <- price_demand_ID + (price_demand_ID < nrow(bid_demand) - 1)
      moc_rl[price_rl_ID, 1] <- bid_demand[price_demand_ID, 1]
      moc_rl[price_rl_ID, 2] <- moc_rl[price_rl_ID - 1, 2]
      
      price_rl_ID <- price_rl_ID + 1
      moc_rl[price_rl_ID, 1] <- bid_demand[price_demand_ID, 1]
      moc_rl[price_rl_ID, 2] <- moc_rl[price_rl_ID - 1, 2] + bid_demand[price_demand_ID, 2]
    }else{
      price_supply_ID <- price_supply_ID + (price_supply_ID < nrow(bid_supply) - 1)
      moc_rl[price_rl_ID, 1] <- bid_supply[price_supply_ID, 1]
      moc_rl[price_rl_ID, 2] <- moc_rl[price_rl_ID - 1, 2]
      
      price_rl_ID <- price_rl_ID + 1
      moc_rl[price_rl_ID, 1] <- bid_supply[price_supply_ID, 1]
      moc_rl[price_rl_ID, 2] <- moc_rl[price_rl_ID - 1, 2] + bid_supply[price_supply_ID, 2] 
    }
  }
  moc_rl[nrow(moc_rl), 1] <- M
  moc_rl[nrow(moc_rl), 2] <- sum(bid_supply[, 2])
  
  return(moc_rl)
}

## Price finding func for residual load curve
price_find <- function(rl, residual_load_curve){
  if(residual_load_curve[1] > rl){
    nodal_price <- price[1]
  }else if(residual_load_curve[num_price] < rl){
    nodal_price <- price[num_price + 1]
  }else{
    price_low <- 1
    price_high <- num_price + 1
    
    while(price_high - price_low > 1){
      price_mid <- round(.5 * (price_low + price_high))
      
      if((residual_load_curve[price_low] - rl) * (residual_load_curve[price_mid] - rl) <= 0){
        price_high <- price_mid
      }else{
        price_low <- price_mid
      }
    }
    nodal_price <- price[price_mid]
  }
  return(nodal_price)
}

### Network Parameters
num_node <- 4
num_line <- 3
z_img <- 1E-2
Cond = diag(rep(1 / z_img, num_line))
V_limit <- pi / 6
I_limit <- 100

## Topology of the network
Topology = matrix(nrow = num_line, ncol = 2)
colnames(Topology) = c("From", "To")

for(i in 1:nrow(Topology)){
  Topology[i, ] = c(i, i+1)
}
# Topology[nrow(Topology), ] = c(nrow(Topology), 1)

# Topology in form of Node * Line
NL = matrix(0, ncol = num_node, nrow = num_line)
for(i in 1:num_line){
  NL[i, Topology[i, 1]] = 1
  NL[i, Topology[i, 2]] = -1
}

## Equivalent Admittance Matrix
## Using inverse matrix is dangerous!!
Y_eq <- t(NL) %*% Cond %*% NL
Z_eq <- solve(Y_eq[2:num_node, 2:num_node])

## ADMM using iterative price guessing
# Augmented Lagrangian:
# L(x, w, u) = f(x) + g(w) + .5 * \rho * ||A %*% x + B %*% w - c + u||^2
# For symmetric box constraints, c will cancel out when optimizing x
# x = power source / sink
# w = box constraint indicator for power balance, voltage, power source / sink, and line current
# c = boundary of box constraints
# u = dual variables of the constraints
# f(x): cost function associated with S (cost function from reference node implicitly incorporated); assume constant at each iteration
# g(w) = indicator function for the box constraint
num_box <- 2 * num_node + num_line
num_constraints <- 2 * num_box
num_variable <- num_node
A <- matrix(0, nrow = num_constraints, ncol = num_variable)

# Power Balancing Equation
# {1} %*% S = -err_eq
# -{1} %*% S = -err_eq
A[1, 1:num_node] <- 1
A[2, 1:num_node] <- -1
constraint_now <- 2

# Box Constraint for Voltage
# Z_eq %*% S - w_{V, neg} = -V_limit
# -Z_eq %*% S - w_{V, pos} = -V_limit
A[constraint_now + 1:(num_node - 1), 2:num_node] <- Z_eq
constraint_now <- constraint_now + (num_node - 1)
A[constraint_now + 1:(num_node - 1), 2:num_node] <- -Z_eq
constraint_now <- constraint_now + (num_node - 1)

# Box Constraint for Power Source / Sink
# S - w_[S, neg] = S_min
# -S - w_{S, pos} = -S_max
A[constraint_now + 1:num_node, 1:num_node] <- diag(rep(1, num_node))
constraint_now <- constraint_now + num_node
A[constraint_now + 1:num_node, 1:num_node] <- diag(rep(-1, num_node))
constraint_now <- constraint_now + num_node

# Box Constraint for Line Current
# (Cond %*% NL)[, 2:num_node] %*% Z_eq %*% S - w_{I, neg} = -I_limit
# -(Cond %*% NL)[, 2:num_node] %*% Z_eq %*% S - w_{I, pos} = -I_limit
A[constraint_now + 1:num_line, 2:num_node] <- (Cond %*% NL)[, 2:num_node] %*% Z_eq
constraint_now <- constraint_now + num_line
A[constraint_now + 1:num_line, 2:num_node] <- -(Cond %*% NL)[, 2:num_node] %*% Z_eq

### Market Parameters
big_num <- 1E16
rl_boundary <- matrix(nrow = num_node, ncol = 2)
bid_demand <- list()
bid_supply <- list()
moc_rl <- list()
bid_demand[[1]] <- matrix(0, ncol = 2, nrow = 2)
bid_demand[[1]][1, ] <- c(-big_num, 0)
bid_demand[[1]][2, ] <- c(big_num, 0)
bid_supply[[1]] <- matrix(0, ncol = 2, nrow = 3)
bid_supply[[1]][1, ] <- c(-big_num, 0)
bid_supply[[1]][2, ] <- c(0, 100)
bid_supply[[1]][3, ] <- c(big_num, 0)
moc_rl[[1]] <- moc_create(bid_demand[[1]], bid_supply[[1]], big_num)
rl_boundary[1, ] <- range(moc_rl[[1]][, 2])

for(node_iter in 2:num_node){
  bid_demand[[node_iter]] <- matrix(0, ncol = 2, nrow = 3)
  bid_demand[[node_iter]][1, ] <- c(-big_num, 0)
  bid_demand[[node_iter]][2, ] <- c(100, 10)
  bid_demand[[node_iter]][3, ] <- c(big_num, 0)
  
  bid_supply[[node_iter]] <- matrix(0, ncol = 2, nrow = 2)
  bid_supply[[node_iter]][1, ] <- c(-big_num, 0)
  bid_supply[[node_iter]][2, ] <- c(big_num, 0)  
  
  moc_rl[[node_iter]] <- moc_create(bid_demand[[node_iter]], bid_supply[[node_iter]], big_num)
  rl_boundary[node_iter, ] <- range(moc_rl[[node_iter]][, 2])
}

### Main optimization process
# f(x) + g(w) + .5 * \rho * ||A %*% x - w - c + u||^2
# f(x) + g(w) + .5 * \rho * ||a_i %*% x_i + (A_!i %*% x_!i - w - c + u)||^2
# f(x) + g(w) + .5 * \rho * (a_i^2 x_i^2 + 2 * (b * a_i) %*% x_i + ...)
# KKT: p(x) + \rho * (a_i^2 x_i + (b * a_i)) = 0
tol <- 1E-6
rho <- 1
sensitivity_x <- c()
for(var_iter in 1:num_variable){
  sensitivity_x[var_iter] <- t(A[, var_iter]) %*% A[, var_iter]
}
x <- rep(0, num_variable)
x_best <- rep(0, num_variable)
w <- rep(0, num_constraints)
boundary <- c(0, 0, rep(-V_limit, 2 * (num_node - 1)), rl_boundary[, 1], -rl_boundary[, 2], rep(-I_limit, 2 * num_line))
u <- rep(0, num_constraints)
price_margin <- rep(0, num_node)
error_min <- Inf

loop <- 0
while(TRUE){
  # Update state variables
  for(node_iter in 1:num_node){
    constant <- A %*% x - w - boundary + u # Sequential if this is inside the loop for variables
    constant_temp <- constant - x[node_iter] * A[, node_iter]
    
    # Bisection method
    gap_price_ID <- c(1, nrow(moc_rl[[node_iter]]))
    mid_price_ID <- round(mean(gap_price_ID))
    while(diff(gap_price_ID) > 1){
      gap_rent <- moc_rl[[node_iter]][gap_price_ID, 1] + rho * (as.vector(sensitivity_x[node_iter] * moc_rl[[node_iter]][gap_price_ID, 2]) + c(t(constant_temp) %*% A[, node_iter]))
      mid_rent <- moc_rl[[node_iter]][mid_price_ID, 1] + rho * (as.vector(sensitivity_x[node_iter] * moc_rl[[node_iter]][mid_price_ID, 2]) + c(t(constant_temp) %*% A[, node_iter]))
      
      if(mid_rent * gap_rent[1] < 0){
        gap_price_ID[2] <- mid_price_ID
      }else{
        gap_price_ID[1] <- mid_price_ID
      }
      
      mid_price_ID <- round(mean(gap_price_ID))
    }
    
    # Update price and power
    if(gap_price_ID[1] %% 2 == 0){
      # price fixed case
      price_margin[node_iter] <- moc_rl[[node_iter]][mid_price_ID, 1]
      x[node_iter] <- (-price_margin[node_iter] / rho - t(constant_temp) %*% A[, node_iter]) / sensitivity_x[node_iter]
      
    }else{
      # power fixed case
      x[node_iter] <- moc_rl[[node_iter]][mid_price_ID, 2]
      price_margin[node_iter] <- -rho * (sensitivity_x[node_iter] * x[node_iter] + t(constant_temp) %*% A[, node_iter]) 
    }
  }
  
  # Update constraint indicators
  w <- A %*% x - boundary + u
  w <- w * (w > 0)
  
  # Update dual variables
  u <- u + A %*% x - w - boundary
  
  # u[5]
  # prime_error[5]
  # range(constant)
  
  # Check convergence
  prime_error <- A %*% x - w - boundary
  dual_error <- price_margin / rho + t(A) %*% u  
  
  if(max(prime_error^2, dual_error^2) < error_min){
    error_min <- max(prime_error^2, dual_error^2)
    x_best <- x
    
    if(max(prime_error^2, dual_error^2) < tol){
      break
    }
  }
  
  if(loop %% 10000 == 0){
    print(paste("Loop ", loop, " :"))
    print(paste("Prime Error: ", max(prime_error^2)))
    print(paste("Dual Error: ", max(dual_error^2)))
    print(paste("Multiplied: ", sum(prime_error^2) * sum(dual_error^2)))
    print(price_margin)
  }  
  loop <- loop + 1 
}

S <- x_best
V <- c(0, Z_eq %*% S[2:num_node])
I <- Cond %*% NL %*% V
print(S)
print(V)
print(t(I))