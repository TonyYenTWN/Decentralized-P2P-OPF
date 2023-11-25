### Decentralized optimization using ADMM, relaxed constraints
## DC Power Flow Assumption, ignoring reactive power flow
## Completely sequential might still not converge...

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

### Network Parameters
num_node <- 3
num_line <- num_node - 1
z_img <- 1 / num_line
Cond = diag(rep(1 / z_img, num_line))
V_limit <- pi / 6
I_limit <- .9

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

## ADMM
# Augmented Lagrangian:
# L(x, u) = f(x) + .5 * \rho * ||A %*% x + u||^2
# x = scaled voltage, power source / sink, and line current
# u = dual variables of the constraints
# f(x): cost function associated with V, S, and I (f(V) and f(I) are well functions)
num_variable <- 2 * num_node + num_line
num_constraints <- num_node + num_line
A <- matrix(0, nrow = num_constraints, ncol = num_variable)

# # Power Flow Equation
# # Y_eq %*% V - S = 0
# A[1:num_node, 1:num_node] <- Y_eq
# A[1:num_node, num_node + 1:num_node] <- diag(rep(-1, num_node))
# constraint_now <- num_node

# Power Balance Equation
# S - t(NL) %*% I = 0
A[1:num_node, num_node + 1:num_node] <- diag(rep(1, num_node))
A[1:num_node, 2 * num_node + 1:num_line] <- -t(NL)
constraint_now <- num_node

# Line Current Equation
# Y_l %*% NL %*% V - I = 0
# set U = mean(Y_l) V
A[constraint_now + 1:num_line, 1:num_node] <- Cond %*% NL
A[constraint_now + 1:num_line, 2 * num_node + 1:num_line] <- diag(rep(-1, num_line))
constraint_now <- constraint_now + num_line

D <- diag(diag(t(A) %*% A))
eigen(solve(D) %*% t(A) %*% A)

### Market Parameters and Cost FUnctions
# Cost functions for V and I
big_num <- 1E16
cost_V <- matrix(ncol = 2, nrow = 4)
cost_V[1, ] <- c(-big_num, -V_limit)
cost_V[2, ] <- c(0, -V_limit)
cost_V[3, ] <- c(0, V_limit)
cost_V[4, ] <- c(big_num, V_limit)
cost_V <- cost_V

cost_I <- matrix(ncol = 2, nrow = 4)
cost_I[1, ] <- c(-big_num, -I_limit)
cost_I[2, ] <- c(0, -I_limit)
cost_I[3, ] <- c(0, I_limit)
cost_I[4, ] <- c(big_num, I_limit)

# Cost functions for S
num_price <- 10
bid_demand <- list()
bid_supply <- list()
moc_rl <- list()
bid_demand[[1]] <- matrix(0, ncol = 2, nrow = 2)
bid_demand[[1]][1, ] <- c(-big_num, 0)
bid_demand[[1]][2, ] <- c(big_num, 0)
bid_supply[[1]] <- matrix(0, ncol = 2, nrow = 3)
bid_supply[[1]][1, ] <- c(-big_num, 0)
bid_supply[[1]][2, ] <- c(0, 2)
bid_supply[[1]][3, ] <- c(big_num, 0)
moc_rl[[1]] <- moc_create(bid_demand[[1]], bid_supply[[1]], big_num)

for(node_iter in 2:num_node){
  bid_demand[[node_iter]] <- matrix(0, ncol = 2, nrow = num_price + 2)
  bid_demand[[node_iter]][1, ] <- c(-big_num, 0)
  for(price_iter in 1:num_price){
    bid_demand[[node_iter]][price_iter + 1, ] <- c(price_iter, 2 / (num_node - 1) / num_price)
  }
  bid_demand[[node_iter]][num_price + 2, ] <- c(big_num, 0)
  
  bid_supply[[node_iter]] <- matrix(0, ncol = 2, nrow = 2)
  bid_supply[[node_iter]][1, ] <- c(-big_num, 0)
  bid_supply[[node_iter]][2, ] <- c(big_num, 0)  
  
  moc_rl[[node_iter]] <- moc_create(bid_demand[[node_iter]], bid_supply[[node_iter]], big_num)
}

# Group cost functions
cost_x <- list()
var_ID <- 1
for(node_iter in 1:num_node){
  cost_x[[var_ID]] <- cost_V
  var_ID <- var_ID + 1
}
for(node_iter in 1:num_node){
  cost_x[[var_ID]] <- moc_rl[[node_iter]]
  var_ID <- var_ID + 1
}
for(line_iter in 1:num_line){
  cost_x[[var_ID]] <- cost_I
  var_ID <- var_ID + 1
}

# Scaling of variables
# x' = [1 / m] %*% (x - x_0)
# m = (ub - lb) / 2
# x_0 = (ub + lb) / 2
# then A %*% x = c -> A %*% ([m] %*% x' + x_0) = c
# A %*% [m] %*% x' = c - A %*% x_0
# A' = A %*% [1 / m]
# c' = c - A %*% x_0
# also -1 <= x' <= 1
# grad(f(x')) = [m] %*% grad(f(x))
scale_var <- rep(1, num_variable)
shift_var <- rep(0, num_variable)
boundary <- rep(0, num_constraints)
A_trans <- A
for(var_iter in 1:num_variable){
  scale_var[var_iter] <- max(cost_x[[var_iter]][, 2]) - min(cost_x[[var_iter]][, 2])
  scale_var[var_iter] <- scale_var[var_iter] / 2
  shift_var[var_iter] <- max(cost_x[[var_iter]][, 2]) + min(cost_x[[var_iter]][, 2])
  shift_var[var_iter] <- shift_var[var_iter] / 2
  A_trans[, var_iter] <- A[, var_iter] * scale_var[var_iter]
  cost_x[[var_iter]][, 2] <- (cost_x[[var_iter]][, 2] - shift_var[var_iter]) / scale_var[var_iter]
  cost_x[[var_iter]][, 1] <- cost_x[[var_iter]][, 1] * scale_var[var_iter]
}
boundary <- -A%*%shift_var

# Mat <- t(A) %*% A
# Mat_trans <- t(A_trans) %*% A_trans
# Diag <- diag(1 / diag(Mat))
# Diag_trans <- diag(1 / diag(Mat_trans))
# Test_mat <- Diag %*% Mat
# Test_mat_trans <- Diag_trans %*% Mat_trans
# eigen(Test_mat)
# eigen(Test_mat_trans)

### Main optimization process
# f(x) + .5 * \rho * ||A %*% x - c + u||^2
# f(x) + .5 * \rho * ||a_i %*% x_i + (A_!i %*% x_!i - c + u)||^2
# f(x) + .5 * \rho * (a_i^2 x_i^2 + 2 * (b * a_i) %*% x_i + ...)
# KKT: p(x) + \rho * (a_i^2 x_i + (b * a_i)) = 0
tol <- 1E-6 / num_node
alpha <- 0
rho <- 100
sensitivity_x <- c()
for(var_iter in 1:num_variable){
  sensitivity_x[var_iter] <- t(A_trans[, var_iter]) %*% A_trans[, var_iter]
}
x_prev <- rep(0, num_variable)
x <- rep(0, num_variable)
u <- rep(0, num_constraints)
price_margin <- rep(0, num_variable)

start_time <- Sys.time()
loop <- 0
while(TRUE){
  # Update prime variables and associate prices
  # constant <- A %*% x - boundary + u # Parallel if this is outside the loop for variables
  for(sweep in 1:1){
    for(var_iter in 1:num_variable){
      constant <- A_trans %*% x - boundary + u # Sequential if this is inside the loop for variables
      constant_temp <- constant - x[var_iter] * A_trans[, var_iter]
      inner_prod <- c(t(constant_temp) %*% A_trans[, var_iter])
      
      # Update x
      gap_price_ID <- c(1, nrow(cost_x[[var_iter]]))
      mid_price_ID <- round(mean(gap_price_ID))
      while(diff(gap_price_ID) > 1){
        gap_rent <- cost_x[[var_iter]][gap_price_ID, 1] + rho * (as.vector(sensitivity_x[var_iter] * cost_x[[var_iter]][gap_price_ID, 2]) + inner_prod)
        mid_rent <- cost_x[[var_iter]][mid_price_ID, 1] + rho * (as.vector(sensitivity_x[var_iter] * cost_x[[var_iter]][mid_price_ID, 2]) + inner_prod)
        
        if(mid_rent * gap_rent[1] <= 0){
          gap_price_ID[2] <- mid_price_ID
        }else{
          gap_price_ID[1] <- mid_price_ID
        }
        
        mid_price_ID <- round(mean(gap_price_ID))
      }
      
      # Update price and power
      if(gap_price_ID[1] %% 2 == 0){
        # price fixed case
        price_margin[var_iter] <- cost_x[[var_iter]][mid_price_ID, 1]
        x[var_iter] <- (-price_margin[var_iter] / rho - inner_prod) / sensitivity_x[var_iter]
        
      }else{
        # power fixed case
        x[var_iter] <- cost_x[[var_iter]][mid_price_ID, 2]
        price_margin[var_iter] <- -rho * (sensitivity_x[var_iter] * x[var_iter] + inner_prod)
      }
      
      # Test use, should be commented when normal
      # inner_prod
      # var_iter <- var_iter + 1
      # Test use, should be commented when normal
    }
  }
  x <- alpha * x_prev + (1 - alpha) * x
  x_prev <- x
  
  # Update dual variables
  u <- u + A_trans %*% x - boundary
  
  # Check convergence
  prime_error <- A_trans %*% x - boundary
  dual_error <-  price_margin[num_node + 1:num_node] + rho * t(A_trans[, num_node + 1:num_node]) %*% u
  if(max(prime_error^2, dual_error^2) < tol){
    break
  }
  if(max(abs(dual_error)) != 0){
    if(abs(log(max(abs(prime_error)) / max(abs(dual_error)))) > log(1E3)){
      rho <- rho * sqrt(max(abs(prime_error)) / max(abs(dual_error)))
    }    
  }
  
  if(loop %% 10000 == 0){
    print(paste("Loop ", loop, " :"))
    print(paste("Prime Error: ", max(prime_error^2)))
    print(paste("Dual Error: ", max(dual_error^2)))
    print(paste("Multiplied: ", sum(prime_error^2) * sum(dual_error^2)))
    # print(price_margin)
  }  
  loop <- loop + 1 
}
end_time <- Sys.time()
end_time - start_time
(end_time - start_time) * 60
shift_var + x * scale_var