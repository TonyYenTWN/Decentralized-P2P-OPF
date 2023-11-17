### Decentralized optimization using ADMM
## DC Power Flow Assumption, ignoring reactive power flow
## Convergence guaranteed for sequential method; parallel method needs proof and validation of stabilization method

## Price finding func for residual load curve
price_find <- function(rl, residual_load_curve){
  if(residual_load_curve[1] > rl){
    nodal_price <- price[1]
  }else if(residual_load_curve[num_price] < rl){
    nodal_price <- price[num_price + 2]
  }else{
    price_low <- 1
    price_high <- num_price + 2
    
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
z_img <- .01
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
# x = voltage, power source / sink, line current
# w = box constraint indicator for voltage, power source / sink, and line current
# c = boundary of box constraints
# u = dual variables of the constraints
# f(x): cost function associated with S (cost function from reference node implicitly incorporated); assume constant at each iteration
# g(w) = indicator function for the box constraint
num_eq <- num_node + num_line
num_box <- 2 * num_node + num_line
num_constraints <- num_eq + 2 * num_box
num_variable <- 2 * num_node + num_line
A <- matrix(0, nrow = num_constraints, ncol = num_variable)

# Power flow equation
# Y_eq %*% V - S = 0
A[1:num_node, 1:num_node] <- Y_eq
A[1:num_node, num_node + 1:num_node] <- diag(rep(-1, num_node))

# Line current equation
# Y_l %*% NL %*% V - I = 0
A[num_node + 1:num_line, 1:num_node] <- Cond %*% NL
A[num_node + 1:num_line, 2 * num_node + 1:num_line] <- diag(rep(-1, num_line))

# Voltage box boundary
# V - w_{V, neg} = -V_limit
# -V - w_{V, pos} = -V_limit
A[num_eq + 1:num_node, 1:num_node] <- diag(rep(1, num_node))
A[num_eq + num_node + 1:num_node, 1:num_node] <- diag(rep(-1, num_node))

# Power source / sink box boundary
# S - w_{S, neg} = S_min
# -S - w_{S, pos} = -S_max
A[num_eq + 2 * num_node + 1:num_node, num_node + 1:num_node] <- diag(rep(1, num_node))
A[num_eq + 3 * num_node + 1:num_node, num_node + 1:num_node] <- diag(rep(-1, num_node))

# Line current box boundary
# I - w_{I, neg} = -I_limit
# -I - w_{I, pos} = -I_limit
A[num_eq + 4 * num_node + 1:num_line, 2 * num_node + 1:num_line] <- diag(rep(1, num_line))
A[num_eq + 4 * num_node + num_line + 1:num_line, 2 * num_node + 1:num_line] <- diag(rep(-1, num_line))

### Market Parameters
num_price <- 100
price <- seq(1, 100, length.out = num_price) - .5
price <- c(0, price, 100)
moc_supply <- matrix(0, nrow = num_price, ncol = num_node)
rownames(moc_supply) <- price[1:num_price + 1]
colnames(moc_supply) <- 1:num_node
for(node_iter in 1:num_node){
  moc_supply[(10 * (node_iter - 1) + 1):num_price, node_iter] <- 1
}
moc_demand <- 1 - moc_supply
moc_supply[, 1] <-  moc_supply[, 1] * 10
cumsum_supply <- apply(moc_supply, 2, cumsum)
cumsum_supply <- rbind(rep(0, num_node), cumsum_supply, cumsum_supply[num_price, ])
cumsum_demand <- apply(apply(apply(moc_demand, 2, rev), 2, cumsum), 2, rev)
cumsum_demand <- rbind(cumsum_demand[1, ], cumsum_demand, rep(0, num_node))
residual_load_curve <- cumsum_supply - cumsum_demand 

### Main optimization process
# f(x) + g(w) + .5 * \rho * ||A %*% x - w - c + u||^2
# f(x) + g(w) + .5 * \rho * ||a_i %*% x_i + (A_!i %*% x_!i - w - c + u)||^2
# f(x) + g(w) + .5 * \rho * (a_i^2 x_i^2 + 2 * (b * a_i) %*% x_i + ...)
# KKT: p(x) + \rho * (a_i^2 x_i + (b * a_i)) = 0
tol <- 1E-6
alpha <- 0
rho <- 1
sensitivity_x <- c()
for(var_iter in 1:num_variable){
  sensitivity_x[var_iter] <- t(A[, var_iter]) %*% A[, var_iter]
}
x_prev <- rep(0, num_variable)
x <- rep(0, num_variable)
w <- rep(0, num_constraints)
boundary <- c(rep(0, num_eq), rep(-V_limit, 2 * num_node), residual_load_curve[1, ], -residual_load_curve[num_price + 2, ], rep(-I_limit, 2 * num_line))
u <- rep(0, num_constraints)
price_margin <- rep(0, num_node)

loop <- 0
while(TRUE){
  # Update state variables
  # constant <- A %*% x - w - boundary + u # Parallel if this is outside the loop for variables
  for(var_iter in 1:num_variable){
    constant <- A %*% x - w - boundary + u # Sequential if this is inside the loop for variables
    constant_temp <- constant - x[var_iter] * A[, var_iter]
    if(var_iter <= num_node || var_iter > 2 * num_node){
      x[var_iter] <- -t(constant_temp) %*% A[, var_iter] / sensitivity_x[var_iter]
    }
    else{
      node_ID <- var_iter - num_node
      error <- price + rho * (as.vector(sensitivity_x[var_iter] * residual_load_curve[, node_ID]) + c(t(constant_temp) %*% A[, var_iter]))
      error <- abs(error)
      price_ID <- which.min(error)
      price_margin[node_ID] <- price[price_ID]
      if(price_ID  == 1){
        x_min <- residual_load_curve[price_ID, node_ID]
      }
      else{
        x_min <- residual_load_curve[price_ID - 1, node_ID]
      }
      x_max <- residual_load_curve[price_ID, node_ID]
      x[var_iter] <- (-price_margin[node_ID] / rho - t(constant_temp) %*% A[, var_iter]) / sensitivity_x[var_iter]
      
      # Check if bid is constrained
      x[var_iter] <- (x[var_iter] > x_min) * x[var_iter] + (x[var_iter] <= x_min) * x_min
      x[var_iter] <- (x[var_iter] < x_max) * x[var_iter] + (x[var_iter] >= x_max) * x_max
      price_margin[node_ID] <- -rho * (sensitivity_x[var_iter] * x[var_iter] + t(constant_temp) %*% A[, var_iter]) 
    }
  }
  x <- x_prev * alpha + x * (1 - alpha)
  x_prev <- x
  
  # Update constraint indicators
  w <- A %*% x - boundary + u
  w <- w * (w > 0)
  w[1:num_eq] <- 0
  
  # Update dual variables
  u <- u + A %*% x - w - boundary
  
  # Check if loop should continue
  prime_error <- A %*% x - w - boundary
  dual_error <- c(rep(0, num_node), price_margin, rep(0, num_line)) / rho + t(A) %*% u
  if(max(prime_error^2, dual_error^2) < tol){
    break
  }
  
  if(loop %% 1000 == 0){
    # print(x[num_node + 1:num_node])
    print(max(prime_error^2))
    print(max(dual_error^2))
    # print(range(x[1:num_node]))
    # print(range(x[2 * num_node + 1:num_line]))
  }  
  loop <- loop + 1
}

print(x[num_node + 1:num_node])