### Decentralized optimization using ADMM, full constraint relaxation with appropriate prices, bi-section splitting method
## DC Power Flow Assumption, ignoring reactive power flow

## merit order curve of residual load creation
moc_create <- function(bid_demand, bid_supply, M){
  moc_rl <- matrix(ncol = 4, nrow = 2 * (nrow(bid_demand) + nrow(bid_supply) - 3))
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
  
  ### Calculate objective
  moc_rl[1, 3] <- -sum(bid_demand[, 1] * bid_demand[, 2])
  for(row_iter in 2:nrow(moc_rl)){
    moc_rl[row_iter, 3] <- moc_rl[row_iter - 1, 3] + moc_rl[row_iter - 1, 1] * (moc_rl[row_iter, 2] - moc_rl[row_iter - 1, 2])
  }
  
  return(moc_rl)
}

### Network Parameters
num_node <- 20
num_line <- num_node - 1
z_img <- 1 / num_line
Cond = diag(rep(1 / z_img, num_line))
V_limit <- pi / 6
I_limit <- 1.2

## Topology of the network
Topology = matrix(nrow = num_line, ncol = 2)
colnames(Topology) = c("From", "To")

for(i in 1:nrow(Topology)){
  Topology[i, ] = c(i, i+1)
}

# Topology in form of Node * Line
NL = matrix(0, ncol = num_node, nrow = num_line)
for(i in 1:num_line){
  NL[i, Topology[i, 1]] = 1
  NL[i, Topology[i, 2]] = -1
}

# ## Equivalent Admittance Matrix
# ## Using inverse matrix is dangerous!!
# Y_eq <- t(NL) %*% Cond %*% NL
# Z_eq <- solve(Y_eq[2:num_node, 2:num_node])

## Split Matrix
num_node_sub <- rep(num_node / 2, 2)
num_line_sub <- c((num_line - 1) / 2, (num_line - 1) / 2)
num_variable <- 2 * num_node_sub + num_line_sub
num_constraints <- num_node_sub + num_line_sub
start_node <- c(0, num_node / 2)
start_line <- c(0, (num_line - 1) / 2 + 1)
A <- list()
green_V <- list()
green_I <- list()
boundary_0 <- list()
for(matrix_iter in 1:2){
  A[[matrix_iter]] <- matrix(0, nrow = num_constraints[matrix_iter], ncol = num_variable[matrix_iter])
  
  # Power Balance Equation
  # S - t(NL) %*% I = 0 
  A[[matrix_iter]][1:num_node_sub[matrix_iter], num_node_sub[matrix_iter] + 1:num_node_sub[matrix_iter]] <- diag(rep(1, num_node_sub[matrix_iter]))
  for(line_iter in 1:num_line_sub[matrix_iter]){
    line_ID <- start_line[matrix_iter] + line_iter
    from_ID <- Topology[line_ID, 1] - start_node[matrix_iter]
    to_ID <- Topology[line_ID, 2] - start_node[matrix_iter]
    
    A[[matrix_iter]][from_ID, 2 * num_node_sub[matrix_iter] + line_iter] <- -1
    if(to_ID <= num_node_sub[matrix_iter]){
      A[[matrix_iter]][to_ID, 2 * num_node_sub[matrix_iter] + line_iter] <- 1
    }
  }
  constraint_now <- num_node_sub[matrix_iter]
  
  # Line Current Equation
  # Y_l %*% NL %*% V - I = 0
  for(line_iter in 1:num_line_sub[matrix_iter]){
    line_ID <- start_line[matrix_iter] + line_iter
    from_ID <- Topology[line_ID, 1] - start_node[matrix_iter]
    to_ID <- Topology[line_ID, 2] - start_node[matrix_iter]    
    
    A[[matrix_iter]][constraint_now + line_iter, from_ID] <- Cond[line_ID, line_ID]
    if(to_ID <= num_node_sub[matrix_iter]){
      A[[matrix_iter]][constraint_now + line_iter, to_ID] <- -Cond[line_ID, line_ID]
    }
    A[[matrix_iter]][constraint_now + line_iter, 2 * num_node_sub[matrix_iter] + line_iter] <- -1
  }
  constraint_now <- num_node_sub[matrix_iter] + num_line_sub[matrix_iter]
  
  # Dirchlet boundary conditions
  if(matrix_iter == 1){
    A[[matrix_iter]] <- A[[matrix_iter]][, -num_node_sub[matrix_iter]]
  }else{
    A[[matrix_iter]] <- A[[matrix_iter]][, 2:ncol(A[[matrix_iter]])]
  }
  
  # boundary vector
  boundary_0[[matrix_iter]] <- rep(0, constraint_now)  
  
  # Green function for voltage and line current
  free_node <- num_node_sub[matrix_iter] - 1
  green_V[[matrix_iter]] <- matrix(0, nrow = free_node, ncol = free_node)
  green_I[[matrix_iter]] <- matrix(0, nrow = free_node, ncol = free_node)
  z_cum <- cumsum(1 / diag(Cond)[start_line[matrix_iter] + 1:num_line_sub[matrix_iter]])
  z_cum <- c(0, z_cum)

  if(matrix_iter < 2){
    for(node_iter_1 in 1:free_node){
      z_L <- z_cum[node_iter_1]
      z_R <- max(z_cum) - z_L
      dV <- z_R

      for(node_iter_2 in 1:free_node){
        if(node_iter_2 <= node_iter_1){
          green_V[[matrix_iter]][node_iter_2, node_iter_1] <- dV
        }else{
          z_L_temp <- z_cum[node_iter_2]
          z_R_temp <- max(z_cum) - z_L_temp
          green_V[[matrix_iter]][node_iter_2, node_iter_1] <- dV * z_R_temp / z_R
        }
      }

      for(line_iter in node_iter_1:(start_line[matrix_iter] + num_line_sub[matrix_iter])){
        green_I[[matrix_iter]][line_iter, node_iter_1] <- 1
      }

      # Dirchlet BC case (keep for future)
      # z_L <- z_cum[node_iter_1]
      # z_R <- z_cum[num_line_sub[matrix_iter]] - z_L
      # dV <- 1 / (1 / z_L + 1 / z_R)
      #
      # for(node_iter_2 in 1:free_node){
      #   z_L_temp <- z_cum[node_iter_2]
      #   if(node_iter_2 <= node_iter_1){
      #     green_V[[matrix_iter]][node_iter_2, node_iter_1] <- dV * z_L_temp / z_L
      #   }else{
      #     z_R_temp <- z_cum[num_line_sub[matrix_iter]] - z_L_temp
      #     green_V[[matrix_iter]][node_iter_2, node_iter_1] <- dV * z_R_temp / z_R
      #   }
      # }
      #
      # for(line_iter in 1:nrow(green_I[[matrix_iter]])){
      #   if(line_iter <= node_iter_1){
      #     green_I[[matrix_iter]][line_iter, node_iter_1] <- -dV / z_L
      #   }else{
      #     green_I[[matrix_iter]][line_iter, node_iter_1] <- dV / z_R
      #   }
      # }
    }
  }else{
    for(node_iter_1 in 1:free_node){
      z_L <- z_cum[node_iter_1 + 1]
      dV <- z_L

      for(node_iter_2 in 1:free_node){
        if(node_iter_2 <= node_iter_1){
          z_L_temp <- z_cum[node_iter_2 + 1]
          green_V[[matrix_iter]][node_iter_2, node_iter_1] <- dV * z_L_temp / z_L
        }else{
          green_V[[matrix_iter]][node_iter_2, node_iter_1] <- dV
        }
      }

      for(line_iter in 1:node_iter_1){
        green_I[[matrix_iter]][line_iter, node_iter_1] <- -1
      }
    }
  }
}
num_variable <- num_variable - 1

### Market Parameters and Cost Functions
price_limit <- c(-500, 3000)

## Cost functions for V and I
# Assign appropriate price to constraints
big_num <- 1E16
cost_coeff_V <- rep(0, 2)
cost_coeff_I <- rep(0, 2)
for(matrix_iter in 1:2){
  for(node_iter_1 in 1:ncol(green_V[[matrix_iter]])){
    # Exchange with boundary
    green_I_temp <- green_I[[matrix_iter]][, node_iter_1]
    green_I_temp <- green_I_temp + big_num * (green_I_temp == 0)
    cost_coeff_V[matrix_iter] <- max(cost_coeff_V[matrix_iter], max(1 / abs(green_V[[matrix_iter]][, node_iter_1])))
    cost_coeff_I[matrix_iter] <- max(cost_coeff_I[matrix_iter], max(1 / abs(green_I_temp)))
    
    # Between nodes
    if(node_iter_1 < ncol(green_V[[matrix_iter]])){
      for(node_iter_2 in (node_iter_1 + 1):ncol(green_V[[matrix_iter]])){
        green_V_temp <- green_V[[matrix_iter]][, node_iter_1] - green_V[[matrix_iter]][, node_iter_2]
        green_V_temp <- green_V_temp + big_num * (green_V_temp == 0)
        green_I_temp <- green_I[[matrix_iter]][, node_iter_1] - green_I[[matrix_iter]][, node_iter_2]
        green_I_temp <- green_I_temp + big_num * (green_I_temp == 0)
        
        cost_coeff_V[matrix_iter] <- max(cost_coeff_V[matrix_iter], max(1 / abs(green_V_temp)))
        cost_coeff_I[matrix_iter] <- max(cost_coeff_I[matrix_iter], max(1 / abs(green_I_temp)))        
      }
    }
  }
}

bid_demand_V <- list()
bid_supply_V <- list()
bid_demand_I <- list()
bid_supply_I <- list()
cost_V <- list()
cost_I <- list()
for(matrix_iter in 1:2){
  bid_demand_V[[matrix_iter]] <- matrix(ncol = 2, nrow = 4)
  bid_demand_V[[matrix_iter]][1, ] <- c(-big_num, 0)
  bid_demand_V[[matrix_iter]][2, ] <- c(-diff(price_limit) * cost_coeff_V[matrix_iter], 1000)
  bid_demand_V[[matrix_iter]][3, ] <- c(0, V_limit)
  bid_demand_V[[matrix_iter]][4, ] <- c(big_num, 0)
  
  bid_supply_V[[matrix_iter]] <- matrix(ncol = 2, nrow = 4)
  bid_supply_V[[matrix_iter]][1, ] <- c(-big_num, 0)
  bid_supply_V[[matrix_iter]][2, ] <- c(0, V_limit)
  bid_supply_V[[matrix_iter]][3, ] <- c(diff(price_limit) * cost_coeff_V[matrix_iter], 1000)
  bid_supply_V[[matrix_iter]][4, ] <- c(big_num, 0)  
  
  cost_V[[matrix_iter]] <- moc_create(bid_demand_V[[matrix_iter]], bid_supply_V[[matrix_iter]], big_num)
  
  bid_demand_I[[matrix_iter]] <- matrix(ncol = 2, nrow = 4)
  bid_demand_I[[matrix_iter]][1, ] <- c(-big_num, 0)
  bid_demand_I[[matrix_iter]][2, ] <- c(-diff(price_limit) * cost_coeff_I[matrix_iter], 1000)
  bid_demand_I[[matrix_iter]][3, ] <- c(0, I_limit)
  bid_demand_I[[matrix_iter]][4, ] <- c(big_num, 0)
  
  bid_supply_I[[matrix_iter]] <- matrix(ncol = 2, nrow = 4)
  bid_supply_I[[matrix_iter]][1, ] <- c(-big_num, 0)
  bid_supply_I[[matrix_iter]][2, ] <- c(0, I_limit)
  bid_supply_I[[matrix_iter]][3, ] <- c(diff(price_limit) * cost_coeff_I[matrix_iter], 1000)
  bid_supply_I[[matrix_iter]][4, ] <- c(big_num, 0)    
  
  cost_I[[matrix_iter]] <- moc_create(bid_demand_I[[matrix_iter]], bid_supply_I[[matrix_iter]], big_num)
}

# Cost functions for S
num_price <- 10
bid_demand <- list()
bid_supply <- list()
moc_rl <- list()
bid_demand[[1]] <- matrix(0, ncol = 2, nrow = 3)
bid_demand[[1]][1, ] <- c(-big_num, 0)
bid_demand[[1]][2, ] <- c(price_limit[1], 1000)
bid_demand[[1]][3, ] <- c(big_num, 0)
bid_supply[[1]] <- matrix(0, ncol = 2, nrow = 4)
bid_supply[[1]][1, ] <- c(-big_num, 0)
bid_supply[[1]][2, ] <- c(0, 2)
bid_supply[[1]][3, ] <- c(price_limit[2], 1000)
bid_supply[[1]][4, ] <- c(big_num, 0)
moc_rl[[1]] <- moc_create(bid_demand[[1]], bid_supply[[1]], big_num)

for(node_iter in 2:num_node){
  bid_demand[[node_iter]] <- matrix(0, ncol = 2, nrow = num_price + 3)
  bid_demand[[node_iter]][1, ] <- c(-big_num, 0)
  bid_demand[[node_iter]][2, ] <- c(price_limit[1], 1000)
  for(price_iter in 1:num_price){
    bid_demand[[node_iter]][price_iter + 2, ] <- c(price_iter, 2 / (num_node - 1) / num_price)
  }
  bid_demand[[node_iter]][num_price + 3, ] <- c(big_num, 0)
  
  bid_supply[[node_iter]] <- matrix(0, ncol = 2, nrow = 3)
  bid_supply[[node_iter]][1, ] <- c(-big_num, 0)
  bid_supply[[node_iter]][2, ] <- c(price_limit[2], 1000)
  bid_supply[[node_iter]][3, ] <- c(big_num, 0)
  
  moc_rl[[node_iter]] <- moc_create(bid_demand[[node_iter]], bid_supply[[node_iter]], big_num)
}

# Group cost functions
cost_x <- list()
for(matrix_iter in 1:2){
  cost_x_temp <- list()
  
  var_ID <- 1
  for(node_iter in 2:num_node_sub[matrix_iter]){
    cost_x_temp[[var_ID]] <- cost_V[[matrix_iter]]
    var_ID <- var_ID + 1
  }
  for(node_iter in 1:num_node_sub[matrix_iter]){
    node_ID <- start_node[matrix_iter] + node_iter
    cost_x_temp[[var_ID]] <- moc_rl[[node_ID]]
    var_ID <- var_ID + 1
  }
  for(line_iter in 1:num_line_sub[matrix_iter]){
    cost_x_temp[[var_ID]] <- cost_I[[matrix_iter]]
    var_ID <- var_ID + 1
  }  
  
  cost_x[[matrix_iter]] <- cost_x_temp
}

### Main optimization process
control_node_ID <- num_node / 2 + 0:1
rho <- 100
tol_error <- 1E-8
tol_root <- 1E-4
sensitivity_x <- list()
for(matrix_iter in 1:2){
  sensitivity_x[[matrix_iter]] <- rep(0, num_variable[matrix_iter])
  for(var_iter in 1:num_variable[matrix_iter]){
    sensitivity_x[[matrix_iter]][var_iter] <- t(A[[matrix_iter]][, var_iter]) %*% A[[matrix_iter]][, var_iter]
  }
}
## Calculate default price gap
for(matrix_iter in 1:2){
  for(var_iter in 1:num_variable[matrix_iter]){
    cost_x[[matrix_iter]][[var_iter]][, 4] <- cost_x[[matrix_iter]][[var_iter]][, 1] / rho + as.vector(sensitivity_x[[matrix_iter]][var_iter] * cost_x[[matrix_iter]][[var_iter]][, 2])
  }  
}
x <- list()
u <- list()
price_margin <- list()
V_control_prev <- rep(0, 2)

V_control_1 <- c()
dV <- V_limit / 2
V_control_1[2] <- 0
V_control_1[1] <- V_control_1[2] - dV
V_control_1[3] <- V_control_1[2] + dV

start_time <- Sys.time()
while(TRUE){
  V_control_2 <- c()
  V_control_lb <- max(-V_limit, V_control_1[2] - I_limit / diag(Cond)[start_line[2]])
  V_control_ub <- min(V_limit, V_control_1[2] + I_limit / diag(Cond)[start_line[2]])
  dV <- V_control_ub - V_control_lb
  dV <- dV / 4
  V_control_2[2] <- V_control_1[2]
  V_control_2[1] <- V_control_2[2] - dV
  V_control_2[3] <- V_control_2[2] + dV
  obj <- rep(NA, 3)
  while(dV > tol_root){
    for(control_iter in 1:3){
      if(!is.na(obj[2])){
        if(control_iter == 2){
          next
        }      
      }
      
      boundary <- boundary_0
      boundary[[1]][num_node_sub[1]] <- V_control_1[2] - V_control_2[control_iter]
      boundary[[1]][num_node_sub[1]] <- boundary[[1]][num_node_sub[1]] * diag(Cond)[start_line[2]]
      boundary[[1]][num_constraints[1]] <- V_control_1[2] * diag(Cond)[start_line[2]]
      boundary[[2]][1] <- V_control_2[control_iter] - V_control_1[2]
      boundary[[2]][1] <- boundary[[2]][1] * diag(Cond)[start_line[2]]
      boundary[[2]][num_node_sub[2] + 1] <- -V_control_2[control_iter] * diag(Cond)[start_line[2]]
      obj[control_iter] <- 0
      
      for(matrix_iter in 1:2){
        # Initialization of prime and dual variables
        x[[matrix_iter]] <- rep(0, num_variable[matrix_iter])
        u[[matrix_iter]] <- rep(0, num_constraints[matrix_iter])
        price_margin[[matrix_iter]] <- rep(0, num_variable[matrix_iter])
        
        while(TRUE){
          obj_temp <- 0
          
          # Update prime variables and associate prices
          for(var_iter in 1:num_variable[matrix_iter]){
            constant <- A[[matrix_iter]] %*% x[[matrix_iter]] - boundary[[matrix_iter]] + u[[matrix_iter]] # Sequential if this is inside the loop for variables
            constant_temp <- constant - x[[matrix_iter]][var_iter] * A[[matrix_iter]][, var_iter]
            inner_prod <- sum(constant_temp * A[[matrix_iter]][, var_iter])
            
            # Update x
            gap_price_ID <- c()
            gap_price_ID[1] <- 1
            gap_price_ID[3] <- nrow(cost_x[[matrix_iter]][[var_iter]])
            gap_price_ID[2] <- (gap_price_ID[1] + gap_price_ID[3]) %/% 2
            while(gap_price_ID[3] - gap_price_ID[1] > 1){
              # gap_rent <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID, 1] + rho * (as.vector(sensitivity_x[[matrix_iter]][var_iter] * cost_x[[matrix_iter]][[var_iter]][gap_price_ID, 2]) + inner_prod)
              gap_rent <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID, 4] + inner_prod
              
              if(gap_rent[1] * gap_rent[2] <= 0){
                gap_price_ID[3] <- gap_price_ID[2]
              }else{
                gap_price_ID[1] <- gap_price_ID[2]
              }
              
              gap_price_ID[2] <- (gap_price_ID[1] + gap_price_ID[3]) %/% 2
            }
            
            # Update price and power
            if(gap_price_ID[1] %% 2 == 0){
              # price fixed case
              price_margin[[matrix_iter]][var_iter] <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID[2], 1]
              x[[matrix_iter]][var_iter] <- (-price_margin[[matrix_iter]][var_iter] / rho - inner_prod) / sensitivity_x[[matrix_iter]][var_iter]
              
              ratio <- x[[matrix_iter]][var_iter] - cost_x[[matrix_iter]][[var_iter]][gap_price_ID[1], 2]
              ratio <- ratio / (cost_x[[matrix_iter]][[var_iter]][gap_price_ID[3], 2] - cost_x[[matrix_iter]][[var_iter]][gap_price_ID[1], 2])
              obj_temp <- obj_temp + (1 - ratio) * cost_x[[matrix_iter]][[var_iter]][gap_price_ID[1], 3]
              obj_temp <- obj_temp + ratio * cost_x[[matrix_iter]][[var_iter]][gap_price_ID[3], 3]
              
            }else{
              # power fixed case
              x[[matrix_iter]][var_iter] <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID[2], 2]
              price_margin[[matrix_iter]][var_iter] <- -rho * (sensitivity_x[[matrix_iter]][var_iter] * x[[matrix_iter]][var_iter] + inner_prod)
              obj_temp <- obj_temp + cost_x[[matrix_iter]][[var_iter]][gap_price_ID[2], 3]
            }
          }
          
          # Update dual variables
          u[[matrix_iter]] <- u[[matrix_iter]] + A[[matrix_iter]] %*% x[[matrix_iter]] - boundary[[matrix_iter]]
          
          # Check convergence
          prime_error <- A[[matrix_iter]] %*% x[[matrix_iter]] - boundary[[matrix_iter]]
          dual_error <-  price_margin[[matrix_iter]] + rho * t(A[[matrix_iter]]) %*% u[[matrix_iter]]
          if(max(prime_error^2, dual_error^2) < tol_error){
            break
          }       
        }
        
        # Update objective value
        obj[control_iter] <- obj[control_iter] + obj_temp
      }
    }
    
    # Choose point with smallest obj value
    dV <- dV / 2
    V_control_2[2] <- V_control_2[which.min(obj)]
    obj[2] <- obj[which.min(obj)]
    V_control_2[1] <- V_control_2[2] - dV
    V_control_2[3] <- V_control_2[2] + dV
  }

  ## Check if solution is degenerate in the feasible voltage range
  # If yes, then the optimal solution is reached by optimizing seperating
  if(abs(max(x[[1]][1:(num_node_sub[1] - 1)], x[[2]][1:(num_node_sub[2] - 1)])) < V_limit){
    break
  }
  
  ## Shift voltage to center 
  V_upper_gap <- V_limit - max(x[[1]][1:(num_node_sub[1] - 1)], x[[2]][1:(num_node_sub[2] - 1)])
  V_lower_gap <- V_limit + min(x[[1]][1:(num_node_sub[1] - 1)], x[[2]][1:(num_node_sub[2] - 1)])
  V_control_1 <- V_control_1 + (V_upper_gap - V_lower_gap) / 2
  V_control_2 <- V_control_2 + (V_upper_gap - V_lower_gap) / 2
  if(abs(V_control_1[2] - V_control_prev[1]) + abs(V_control_2[2] - V_control_prev[2]) < tol_root){
    break
  }  
  V_control_prev[1] <- V_control_1[2]
  V_control_prev[2] <- V_control_2[2]  
  
  ### Optimize V_control_1
  V_control_1 <- c()
  V_control_lb <- max(-V_limit, V_control_2[2] - I_limit / diag(Cond)[start_line[2]])
  V_control_ub <- min(V_limit, V_control_2[2] + I_limit / diag(Cond)[start_line[2]])
  dV <- V_control_ub - V_control_lb
  dV <- dV / 4
  V_control_1[2] <- V_control_2[2]
  V_control_1[1] <- V_control_1[2] - dV
  V_control_1[3] <- V_control_1[2] + dV
  obj <- rep(NA, 3)
  while(dV > tol_root){
    for(control_iter in 1:3){
      if(!is.na(obj[2])){
        if(control_iter == 2){
          next
        }      
      }
      
      boundary <- boundary_0
      boundary[[1]][num_node_sub[1]] <- V_control_1[control_iter] - V_control_2[2]
      boundary[[1]][num_node_sub[1]] <- boundary[[1]][num_node_sub[1]] * diag(Cond)[start_line[2]]
      boundary[[1]][num_constraints[1]] <- V_control_1[control_iter] * diag(Cond)[start_line[2]]
      boundary[[2]][1] <- V_control_2[2] - V_control_1[control_iter]
      boundary[[2]][1] <- boundary[[2]][1] * diag(Cond)[start_line[2]]
      boundary[[2]][num_node_sub[2] + 1] <- -V_control_2[2] * diag(Cond)[start_line[2]]
      obj[control_iter] <- 0
      
      for(matrix_iter in 1:2){
        # Initialization of prime and dual variables
        x[[matrix_iter]] <- rep(0, num_variable[matrix_iter])
        u[[matrix_iter]] <- rep(0, num_constraints[matrix_iter])
        price_margin[[matrix_iter]] <- rep(0, num_variable[matrix_iter])
        
        while(TRUE){
          obj_temp <- 0
          
          # Update prime variables and associate prices
          for(var_iter in 1:num_variable[matrix_iter]){
            constant <- A[[matrix_iter]] %*% x[[matrix_iter]] - boundary[[matrix_iter]] + u[[matrix_iter]] # Sequential if this is inside the loop for variables
            constant_temp <- constant - x[[matrix_iter]][var_iter] * A[[matrix_iter]][, var_iter]
            inner_prod <- sum(constant_temp * A[[matrix_iter]][, var_iter])
            
            # Update x
            gap_price_ID <- c()
            gap_price_ID[1] <- 1
            gap_price_ID[3] <- nrow(cost_x[[matrix_iter]][[var_iter]])
            gap_price_ID[2] <- (gap_price_ID[1] + gap_price_ID[3]) %/% 2
            while(gap_price_ID[3] - gap_price_ID[1] > 1){
              # gap_rent <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID, 1] + rho * (as.vector(sensitivity_x[[matrix_iter]][var_iter] * cost_x[[matrix_iter]][[var_iter]][gap_price_ID, 2]) + inner_prod)
              gap_rent <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID, 4] + inner_prod
              
              if(gap_rent[1] * gap_rent[2] <= 0){
                gap_price_ID[3] <- gap_price_ID[2]
              }else{
                gap_price_ID[1] <- gap_price_ID[2]
              }
              
              gap_price_ID[2] <- (gap_price_ID[1] + gap_price_ID[3]) %/% 2
            }
            
            # Update price and power
            if(gap_price_ID[1] %% 2 == 0){
              # price fixed case
              price_margin[[matrix_iter]][var_iter] <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID[2], 1]
              x[[matrix_iter]][var_iter] <- (-price_margin[[matrix_iter]][var_iter] / rho - inner_prod) / sensitivity_x[[matrix_iter]][var_iter]
              
              ratio <- x[[matrix_iter]][var_iter] - cost_x[[matrix_iter]][[var_iter]][gap_price_ID[1], 2]
              ratio <- ratio / (cost_x[[matrix_iter]][[var_iter]][gap_price_ID[3], 2] - cost_x[[matrix_iter]][[var_iter]][gap_price_ID[1], 2])
              obj_temp <- obj_temp + (1 - ratio) * cost_x[[matrix_iter]][[var_iter]][gap_price_ID[1], 3]
              obj_temp <- obj_temp + ratio * cost_x[[matrix_iter]][[var_iter]][gap_price_ID[3], 3]
              
            }else{
              # power fixed case
              x[[matrix_iter]][var_iter] <- cost_x[[matrix_iter]][[var_iter]][gap_price_ID[2], 2]
              price_margin[[matrix_iter]][var_iter] <- -rho * (sensitivity_x[[matrix_iter]][var_iter] * x[[matrix_iter]][var_iter] + inner_prod)
              obj_temp <- obj_temp + cost_x[[matrix_iter]][[var_iter]][gap_price_ID[2], 3]
            }
          }
          
          # Update dual variables
          u[[matrix_iter]] <- u[[matrix_iter]] + A[[matrix_iter]] %*% x[[matrix_iter]] - boundary[[matrix_iter]]
          
          # Check convergence
          prime_error <- A[[matrix_iter]] %*% x[[matrix_iter]] - boundary[[matrix_iter]]
          dual_error <-  price_margin[[matrix_iter]] + rho * t(A[[matrix_iter]]) %*% u[[matrix_iter]]
          if(max(prime_error^2, dual_error^2) < tol_error){
            break
          }       
        }
        
        # Update objective value
        obj[control_iter] <- obj[control_iter] + obj_temp
      }
    }
    
    # Choose point with smallest obj value
    dV <- dV / 2
    V_control_1[2] <- V_control_1[which.min(obj)]
    obj[2] <- obj[which.min(obj)]
    V_control_1[1] <- V_control_1[2] - dV
    V_control_1[3] <- V_control_1[2] + dV
  }
  
  ## Check if solution is degenerate in the feasible voltage range
  # If yes, then the optimal solution is reached by optimizing seperating
  if(abs(max(x[[1]][1:(num_node_sub[1] - 1)], x[[2]][1:(num_node_sub[2] - 1)])) < V_limit){
    break
  }
  
  ## Shift voltage to center 
  V_upper_gap <- V_limit - max(x[[1]][1:(num_node_sub[1] - 1)], x[[2]][1:(num_node_sub[2] - 1)])
  V_lower_gap <- V_limit + min(x[[1]][1:(num_node_sub[1] - 1)], x[[2]][1:(num_node_sub[2] - 1)])
  V_control_1 <- V_control_1 + (V_upper_gap - V_lower_gap) / 2
  V_control_2 <- V_control_2 + (V_upper_gap - V_lower_gap) / 2
  if(abs(V_control_1[2] - V_control_prev[1]) + abs(V_control_2[2] - V_control_prev[2]) < tol_root){
    break
  }
  V_control_prev[1] <- V_control_1[2]
  V_control_prev[2] <- V_control_2[2]
}
end_time <- Sys.time()
end_time - start_time
(end_time - start_time) * 60

print(c(V_control_1[2], V_control_2[2]))
plot(1:num_node, c(x[[1]][1:(num_node_sub[1] - 1)], V_control_1[2], V_control_2[2], x[[2]][1:(num_node_sub[2] - 1)]), type = "l")