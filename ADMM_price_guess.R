### Decentralized optimization using ADMM
## DC Power Flow Assumption, ignoring reactive power flow

## Price finding func for residual load curve
price_min <- -500
price_max <- 3000
price_find <- function(rl, residual_load_curve){
  if(residual_load_curve[1] > rl){
    nodal_price <- price_min
  }else if(residual_load_curve[num_price] < rl){
    nodal_price <- price_max
  }else{
    price_low <- 1
    price_high <- price[num_price]

    while(price_high - price_low > 1){
      price_mid <- as.integer(.5 * (price_low + price_high))

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
num_node <- 3
num_line <- 2
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
# x = power source / sink
# w = box constraint indicator for voltage and line current
# c = boundary of box constraints (absolute values)
# u = dual variables of the constraints
# f(x): cost function associated with S (cost function from reference node implicitly incorporated); assume constant at each iteration
# g(w) = indicator function for the box constraint
num_box <- num_node - 1 + num_line
num_constraints <- 2 * num_box
num_variable <- num_node - 1
A <- matrix(0, nrow = num_constraints, ncol = num_variable)
B <- matrix(0, nrow = num_constraints, ncol = num_constraints)

# Box Constraint for Voltage
# Z_eq %*% S - w_{V, neg} = -V_limit
# -Z_eq %*% S - w_{V, pos} = -V_limit
A[1:(num_node - 1), 1:(num_node - 1)] <- Z_eq
A[(num_node - 1) + 1:(num_node - 1), 1:(num_node - 1)] <- -Z_eq
B[1:(2 * (num_node - 1)), 1:(2 * (num_node - 1))] <- diag(rep(-1, (2 * (num_node - 1))))

# Box Constraint for Line Current
# (Cond %*% NL)[, 2:num_node] %*% Z_eq %*% S - w_{I, neg} = -I_limit
# -(Cond %*% NL)[, 2:num_node] %*% Z_eq %*% S - w_{I, pos} = -I_limit
A[2 * (num_node - 1) + 1:num_line, 1:(num_node - 1)] <- (Cond %*% NL)[, 2:num_node] %*% Z_eq
A[2 * (num_node - 1) + num_line + 1:num_line, 1:(num_node - 1)] <- -(Cond %*% NL)[, 2:num_node] %*% Z_eq
B[2 * (num_node - 1) + 1:(2 * num_line), 2 * (num_node - 1) + 1:(2 * num_line)] <- diag(rep(-1, (2 * num_line)))

### Market Parameters
num_price <- 101
price <- seq(0, 100, length.out = num_price)
moc_supply <- matrix(0, nrow = num_price, ncol = num_node)
rownames(moc_supply) <- price
colnames(moc_supply) <- 1:num_node
for(node_iter in 1:num_node){
  moc_supply[(10 * (node_iter - 1) + 1):num_price, node_iter] <- 1
}
moc_demand <- 1 - moc_supply
moc_supply[, 1] <-  moc_supply[, 1] * 10
cumsum_supply <- apply(moc_supply, 2, cumsum)
cumsum_demand <- apply(apply(apply(moc_demand, 2, rev), 2, cumsum), 2, rev)
residual_load_curve <- cumsum_supply - cumsum_demand

### Main optimization process
# optimization
# p + \rho * t(A) %*% (A %*% x - w - c + u) = 0
# \rho * (t(A) %*% A) %*% x = -p + \rho * t(A) %*% (w + c - u)
# (t(A) %*% A) %*% x = -1 / rho * p + t(A) %*% (w + c - u)
tol <- 1E-12
alpha <- .5
rho <- 1
drho <- rho / 2
sensitivity_S <- t(A) %*% A
sensitivity_S_inv <- solve(sensitivity_S)
sensitivity_w <- t(A) %*% B
nodal_price <- rep(0, num_node)
S <- rep(0, num_node - 1)
error <- Inf
w <- rep(0, num_constraints)
c <- -c(rep(V_limit, 2 * (num_node - 1)), rep(I_limit, 2 * num_line))
u <- rep(0, num_constraints)
price_guess_prev <- c(price[1], rep(price[num_price], num_node - 1))
price_guess_now <- price_guess_prev

while(TRUE){
  while(TRUE){
    # Update x variables
    S <- sensitivity_S_inv %*% (-(price_guess_now[2:num_node] - price_guess_now[1]) / rho + t(A) %*% (w + c - u))
    
    # Update w variables
    V <- Z_eq %*% S
    I <- Cond %*% NL %*% c(0, V)
    w <- c(V, -V, I, -I) - c + u 
    w <- w * (w > 0)
    
    # Update u variables
    u <- u + A %*% S - w - c
    
    # Check if loop should continue
    prime_error <- A %*% S - w - c
    dual_error <- (price_guess_now[2:num_node] - price_guess_now[1]) / rho + t(A) %*% u
    if(max(prime_error^2, dual_error^2) < tol){
      break
    }
  }
  
  price_margin <- price_find(-sum(S), residual_load_curve[, 1])
  for(node_iter in 2:num_node){
    price_margin[node_iter] <- price_find(S[node_iter - 1], residual_load_curve[, node_iter])
  }
  price_guess_now <- alpha * price_guess_now + (1 - alpha) * price_margin
  for(node_iter in 1:num_node){
    if(price_margin[node_iter] < price[1]){
      price_margin[node_iter] <- price_min
    }else if(price_margin[node_iter] < price[num_price]){
      price_margin[node_iter] <- price_max
    }
  }
  
  if(max((price_guess_now - price_guess_prev)^2) < tol){
    break
  }
  price_guess_prev <- price_guess_now  
}
print(t(price_guess_now))
t(S)
t(V)
t(I)