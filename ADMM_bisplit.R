### Decentralized optimization using ADMM, relaxed constraints, bi-section splitting method
## DC Power Flow Assumption, ignoring reactive power flow

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
num_node <- 6
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

## Equivalent Admittance Matrix
## Using inverse matrix is dangerous!!
Y_eq <- t(NL) %*% Cond %*% NL
Z_eq <- solve(Y_eq[2:num_node, 2:num_node])

## Split Matrix
num_node_sub <- rep(num_node / 2, 2)
num_line_sub <- rep((num_line - 1) / 2, 2)
num_variable <- 2 * num_node_sub + num_line_sub
num_constraints <- num_node_sub + num_line_sub
start_node <- c(0, num_node / 2)
start_line <- c(0, (num_line - 1) / 2 + 1)
A <- list()
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
    A[[matrix_iter]][to_ID, 2 * num_node_sub[matrix_iter] + line_iter] <- 1
  }
  constraint_now <- num_node_sub[matrix_iter]
  
  # Line Current Equation
  # Y_l %*% NL %*% V - I = 0
  for(line_iter in 1:num_line_sub[matrix_iter]){
    line_ID <- start_line[matrix_iter] + line_iter
    from_ID <- Topology[line_ID, 1] - start_node[matrix_iter]
    to_ID <- Topology[line_ID, 2] - start_node[matrix_iter]    
    
    A[[matrix_iter]][constraint_now + line_iter, from_ID] <- Cond[line_ID, line_ID]
    A[[matrix_iter]][constraint_now + line_iter, to_ID] <- -Cond[line_ID, line_ID]
    A[[matrix_iter]][constraint_now + line_iter, 2 * num_node_sub[matrix_iter] + line_iter] <- -1
  } 
}

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
bid_demand[[1]] <- matrix(0, ncol = 2, nrow = 3)
bid_demand[[1]][1, ] <- c(-big_num, 0)
bid_demand[[1]][1, ] <- c(-1, 2)
bid_demand[[1]][3, ] <- c(big_num, 0)
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
  
  bid_supply[[node_iter]] <- matrix(0, ncol = 2, nrow = 3)
  bid_supply[[node_iter]][1, ] <- c(-big_num, 0)
  bid_supply[[node_iter]][2, ] <- c(num_price + 1, 2 / (num_node - 1)) 
  bid_supply[[node_iter]][3, ] <- c(big_num, 0)  
  
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
cost_x_0 <- cost_x

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
scale_var <- list()
shift_var <- list()
boundary <- list()
A_trans <- A
for(matrix_iter in 1:2){
  scale_var[[matrix_iter]] <- rep(0, num_node_sub[matrix_iter])
  shift_var[[matrix_iter]] <- rep(0, num_node_sub[matrix_iter])
  for(node_iter in 1:num_node_sub[matrix_iter]){
    V_ID <- start_node[matrix_iter] + node_iter
    V_ID_sub <- node_iter
    S_ID <- num_node + start_node[matrix_iter] + node_iter
    S_ID_sub <- num_node_sub[matrix_iter] + node_iter
    
    scale_var[[matrix_iter]][V_ID_sub] <- max(cost_x[[V_ID]][, 2]) - min(cost_x[[V_ID]][, 2])
    scale_var[[matrix_iter]][V_ID_sub] <- scale_var[[matrix_iter]][V_ID_sub] / 2
    shift_var[[matrix_iter]][V_ID_sub] <- max(cost_x[[V_ID]][, 2]) + min(cost_x[[V_ID]][, 2])
    shift_var[[matrix_iter]][V_ID_sub] <- shift_var[[matrix_iter]][V_ID_sub] / 2
    A_trans[[matrix_iter]][, V_ID_sub] <- A[[matrix_iter]][, V_ID_sub] * scale_var[[matrix_iter]][V_ID_sub]
    cost_x[[V_ID]][, 2] <- (cost_x[[V_ID]][, 2] - shift_var[[matrix_iter]][V_ID_sub]) / scale_var[[matrix_iter]][V_ID_sub]
    cost_x[[V_ID]][, 1] <- cost_x[[V_ID]][, 1] * scale_var[[matrix_iter]][V_ID_sub]
    
    scale_var[[matrix_iter]][S_ID_sub] <- max(cost_x[[S_ID]][, 2]) - min(cost_x[[S_ID]][, 2])
    scale_var[[matrix_iter]][S_ID_sub] <- scale_var[[matrix_iter]][S_ID_sub] / 2
    shift_var[[matrix_iter]][S_ID_sub] <- max(cost_x[[S_ID]][, 2]) + min(cost_x[[S_ID]][, 2])
    shift_var[[matrix_iter]][S_ID_sub] <- shift_var[[matrix_iter]][S_ID_sub] / 2
    A_trans[[matrix_iter]][, S_ID_sub] <- A[[matrix_iter]][, S_ID_sub] * scale_var[[matrix_iter]][S_ID_sub]
    cost_x[[S_ID]][, 2] <- (cost_x[[S_ID]][, 2] - shift_var[[matrix_iter]][S_ID_sub]) / scale_var[[matrix_iter]][S_ID_sub]
    cost_x[[S_ID]][, 1] <- cost_x[[S_ID]][, 1] * scale_var[[matrix_iter]][S_ID_sub]    
  }
  
  for(line_iter in 1:num_line_sub[matrix_iter]){
    I_ID <- 2 * num_node + start_line[matrix_iter] + line_iter
    I_ID_sub <- 2 * num_node_sub[matrix_iter] + line_iter 
    
    scale_var[[matrix_iter]][I_ID_sub] <- max(cost_x[[I_ID]][, 2]) - min(cost_x[[I_ID]][, 2])
    scale_var[[matrix_iter]][I_ID_sub] <- scale_var[[matrix_iter]][I_ID_sub] / 2
    shift_var[[matrix_iter]][I_ID_sub] <- max(cost_x[[I_ID]][, 2]) + min(cost_x[[I_ID]][, 2])
    shift_var[[matrix_iter]][I_ID_sub] <- shift_var[[matrix_iter]][I_ID_sub] / 2
    A_trans[[matrix_iter]][, I_ID_sub] <- A[[matrix_iter]][, I_ID_sub] * scale_var[[matrix_iter]][I_ID_sub]
    cost_x[[I_ID]][, 2] <- (cost_x[[I_ID]][, 2] - shift_var[[matrix_iter]][I_ID_sub]) / scale_var[[matrix_iter]][I_ID_sub]
    cost_x[[I_ID]][, 1] <- cost_x[[I_ID]][, 1] * scale_var[[matrix_iter]][I_ID_sub]     
  }
  
  boundary[[matrix_iter]] <- -A[[matrix_iter]]%*%shift_var[[matrix_iter]]
}

# Voltage limits of the control node
V_limit_control <- c()
control_node <- num_node / 2 + 1
green_V <- Z_eq
green_I <- (Cond %*% NL)[, 2:num_node] %*% Z_eq
S_limit <- matrix(ncol = 2, nrow = num_node)
for(node_iter in 1:num_node){
  S_ID <- num_node + node_iter
  S_limit[node_iter, 1] <- min(cost_x_0[[S_ID]][, 2] )
  S_limit[node_iter, 2] <- max(cost_x_0[[S_ID]][, 2])
}

# Find voltage and current at extremes
tol_root <- 1E-12
S_lb <- S_limit[2:num_node, 1]
S_ub <- S_limit[2:num_node, 2]
# the case sum(S_lb) > 0 not in solvable set
if(sum(S_lb) + S_limit[1, 2] < 0){
  gap <- sum(S_lb) + S_limit[1, 2]
  
  # start from right because rhs affects voltage more (we want to know the most extreme case); might be different for other boundary conditions
  for(node_iter in rev(2:num_node - 1)){
    gap_fill <- min(-gap, -S_lb[node_iter])
    S_lb[node_iter] <- S_lb[node_iter] + gap_fill
    gap <- gap + gap_fill
    
    if(gap == 0){
      break
    }
  }
}
if(sum(S_ub) + S_limit[1, 1] > 0){
  gap <- sum(S_ub) + S_limit[1, 1]
  
  for(node_iter in rev(2:num_node - 1)){
    gap_fill <- min(gap, S_ub[node_iter])
    S_ub[node_iter] <- S_ub[node_iter] - gap_fill
    gap <- gap - gap_fill
    
    if(gap == 0){
      break
    }    
  }
}

## Lower voltage limit at control node
# Fix voltage violation first
S_lb_temp <- S_lb
S_lb_range <- c(0, length(S_lb))
V_lb_temp <- green_V %*% S_lb_temp
V_lb_exceed <- min(V_lb_temp) + V_limit < 0
if(V_lb_exceed){
  gap_ID <- c(0, length(S_lb))
  mid_ID <- round(mean(gap_ID))
  
  while(diff(gap_ID) > 1){
    gap_V_min <- c()
    if(gap_ID[1] == 0){
      gap_V_min[1] <- 0
    }else{
      gap_V_min[1] <- sum(green_V[num_node - 1, ] * c(S_lb[1:gap_ID[1]], rep(0, length(S_lb) - gap_ID[1])))
    }
    gap_V_min[2] <- sum(green_V[num_node - 1, ] * c(S_lb[1:gap_ID[2]], rep(0, length(S_lb) - gap_ID[2]))) 
    
    if(mid_ID[1] == 0){
      mid_V_min <- 0
    }else{
      mid_V_min <- sum(green_V[num_node - 1, ] * c(S_lb[1:mid_ID], rep(0, length(S_lb) - mid_ID)))
    }
    
    gap_rent <- gap_V_min + V_limit
    mid_rent <- mid_V_min + V_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_ID[2] <- mid_ID
    }else{
      gap_ID[1] <- mid_ID
    }
    mid_ID <- round(mean(gap_ID))
  }
  
  gap_S <- c(0, S_lb[gap_ID[2]])
  mid_S <- mean(gap_S)
  while(abs(diff(gap_S)) > tol_root){
    gap_V_min <- c()
    if(gap_S[1] == 0 && gap_ID[1] == 0){
      gap_V_min[1] <- 0
    }else{
      gap_V_min[1] <- sum(green_V[num_node - 1, ] * c(S_lb[1:gap_ID[1]], gap_S[1], rep(0, length(S_lb) - gap_ID[1] - 1)))
    }
    gap_V_min[2] <- sum(green_V[num_node - 1, ] * c(S_lb[1:gap_ID[1]], gap_S[2], rep(0, length(S_lb) - gap_ID[1] - 1)))   
    
    mid_V_min <- sum(green_V[num_node - 1, ] * c(S_lb[1:gap_ID[1]], mid_S, rep(0, length(S_lb) - gap_ID[1] - 1)))  
    
    gap_rent <- gap_V_min + V_limit
    mid_rent <- mid_V_min + V_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_S[2] <- mid_S
    }else{
      gap_S[1] <- mid_S
    }
    mid_S <- mean(gap_S)   
  }
  
  S_lb_temp <- c(S_lb[1:gap_ID[1]], mid_S, rep(0, length(S_lb) - gap_ID[1] - 1))
  S_lb_range <- c(0, gap_ID[1] + mid_S / S_lb[gap_ID[2]])
}

# Fix line current violation second
I_lb_temp <- green_I %*% S_lb_temp
I_lb_exceed <- max(I_lb_temp) - I_limit > 0
if(I_lb_exceed){
  top_ID <- as.integer(S_lb_range[2]) + (as.integer(S_lb_range[2]) != S_lb_range[2])
  gap_ID <- c(1, top_ID + 1)
  mid_ID <- round(mean(gap_ID))
  
  while(TRUE){
    gap_I_max <- c()
    gap_I_max[1] <- -sum(S_lb_temp[gap_ID[1]:top_ID]) 
    
    if(gap_ID[2] == top_ID + 1){
      gap_I_max[2] <- 0
    }else{
      gap_I_max[2] <- -sum(S_lb_temp[gap_ID[2]:top_ID])
    }
    
    if(mid_ID == top_ID + 1){
      mid_I_max <- 0
    }else{
      mid_I_max <- -sum(S_lb_temp[mid_ID:top_ID])
    }
    
    if(diff(gap_ID) == 1){
      break
    }    
    
    gap_rent <- gap_I_max - I_limit
    mid_rent <- mid_I_max - I_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_ID[2] <- mid_ID
    }else{
      gap_ID[1] <- mid_ID
    }
    mid_ID <- round(mean(gap_ID))    
  }
  
  gap_I_max_min <- gap_I_max[2]
  gap_S <- c(0, S_lb[gap_ID[1]])
  mid_S <- mean(gap_S)
  while(abs(diff(gap_S)) > tol_root){
    gap_I_max <- gap_I_max_min - gap_S
    mid_I_max <- gap_I_max_min - mid_S
    
    gap_rent <- gap_I_max - I_limit
    mid_rent <- mid_I_max - I_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_S[2] <- mid_S
    }else{
      gap_S[1] <- mid_S
    }
    mid_S <- mean(gap_S)      
  }
  
  for(node_iter in 1:length(S_lb_temp)){
    if(node_iter < gap_ID[1]){
      S_lb_temp[node_iter] <- 0
    }else if(node_iter == gap_ID[1]){
      S_lb_temp[node_iter] <- mid_S
      break
    }
  }
  
  S_lb_range[1] <- gap_ID[1] - mid_S / S_lb[gap_ID[1]]
}

# Check optimal position last
V_lb_temp <- green_V %*% S_lb_temp
V_lb_exceed <- min(V_lb_temp) + V_limit < 0
if(!V_lb_exceed){
  V_lb_gap <- min(V_lb_temp) + V_limit
  while(V_lb_gap > 0){
    low_ID <- as.integer(S_lb_range[1]) + (as.integer(S_lb_range[1]) != S_lb_range[1])
    high_ID <- as.integer(S_lb_range[2]) + (as.integer(S_lb_range[2]) != S_lb_range[2])
    
    exchange_diff <- green_V[num_node - 1, high_ID] - green_V[num_node - 1, low_ID]
    S_move <- -max(S_lb_temp[low_ID], S_lb[high_ID] - S_lb_temp[high_ID])
    S_move <- min(S_move, V_lb_gap / exchange_diff)
    
    S_lb_temp[low_ID] <- S_lb_temp[low_ID] + S_move
    S_lb_temp[high_ID] <- S_lb_temp[high_ID] - S_move
    V_lb_gap <- V_lb_gap - exchange_diff * S_move
    if(S_lb_temp[low_ID] == 0){
      S_lb_range[1] <- as.integer(S_lb_range[1]) + 1
      S_lb_range[2] <- high_ID + S_lb_temp[high_ID] / S_lb[high_ID]
    }else if(S_lb_temp[high_ID] == 0){
      S_lb_range[2] <- as.integer(S_lb_range[2])
      S_lb_range[1] <- low_ID + 1 - S_lb_temp[low_ID] / S_lb[low_ID]
    }
  }
}
S_lb <- S_lb_temp
V_limit_control[1] <- (green_V %*% S_lb)[control_node - 1]

## Upper voltage limit at control node
# Fix voltage violation first
S_ub_temp <- S_ub
S_ub_range <- c(0, length(S_ub))
V_ub_temp <- green_V %*% S_ub_temp
V_ub_exceed <- max(V_ub_temp) - V_limit > 0
if(V_ub_exceed){
  gap_ID <- c(0, length(S_ub))
  mid_ID <- round(mean(gap_ID))
  
  while(TRUE){
    gap_V_max <- c()
    if(gap_ID[1] == 0){
      gap_V_max[1] <- 0
    }else{
      gap_V_max[1] <- sum(green_V[num_node - 1, ] * c(S_ub[1:gap_ID[1]], rep(0, length(S_ub) - gap_ID[1])))
    }
    gap_V_max[2] <- sum(green_V[num_node - 1, ] * c(S_ub[1:gap_ID[2]], rep(0, length(S_ub) - gap_ID[2]))) 
    
    if(mid_ID[1] == 0){
      mid_V_max <- 0
    }else{
      mid_V_max <- sum(green_V[num_node - 1, ] * c(S_ub[1:mid_ID], rep(0, length(S_ub) - mid_ID)))
    }
    
    if(diff(gap_ID) == 1){
      break
    }
    
    gap_rent <- gap_V_max - V_limit
    mid_rent <- mid_V_max - V_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_ID[2] <- mid_ID
    }else{
      gap_ID[1] <- mid_ID
    }
    mid_ID <- round(mean(gap_ID))
  }  
  
  gap_V_max_min <- gap_V_max[1]
  gap_S <- c(0, S_ub[gap_ID[2]])
  mid_S <- mean(gap_S)
  while(abs(diff(gap_S)) > tol_root){
    if(gap_S[1] == 0 && gap_ID[1] == 0){
      gap_V_max[1] <- 0
    }else{
      gap_V_max[1] <- gap_V_max_min + green_V[num_node - 1, gap_ID[2]] * gap_S[1]
    }
    gap_V_max[2] <- gap_V_max_min + green_V[num_node - 1, gap_ID[2]] * gap_S[2]
    
    mid_V_max <- gap_V_max_min + green_V[num_node - 1, mid_ID] * mid_S 
    
    gap_rent <- gap_V_max - V_limit
    mid_rent <- mid_V_max - V_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_S[2] <- mid_S
    }else{
      gap_S[1] <- mid_S
    }
    mid_S <- mean(gap_S)     
  }
  
  S_ub_temp <- c(S_ub[1:gap_ID[1]], mid_S, rep(0, length(S_ub) - gap_ID[1] - 1))
  S_ub_range <- c(0, gap_ID[1] + mid_S / S_ub[gap_ID[2]])  
}

# Fix line current violation second
I_ub_temp <- green_I %*% S_ub_temp
I_ub_exceed <- min(I_ub_temp) + I_limit < 0
if(I_ub_exceed){
  top_ID <- as.integer(S_ub_range[2]) + (as.integer(S_ub_range[2]) != S_ub_range[2])
  gap_ID <- c(1, top_ID + 1)
  mid_ID <- round(mean(gap_ID))  
  
  while(TRUE){
    gap_I_min <- c()
    gap_I_min[1] <- -sum(S_ub_temp[gap_ID[1]:top_ID]) 
    
    if(gap_ID[2] == top_ID + 1){
      gap_I_min[2] <- 0
    }else{
      gap_I_min[2] <- -sum(S_ub_temp[gap_ID[2]:top_ID])
    }
    if(mid_ID == top_ID + 1){
      mid_I_min <- 0
    }else{
      mid_I_min <- -sum(S_ub_temp[mid_ID:top_ID])
    }    
    
    if(diff(gap_ID) == 1){
      break
    }    
    
    gap_rent <- gap_I_min + I_limit
    mid_rent <- mid_I_min + I_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_ID[2] <- mid_ID
    }else{
      gap_ID[1] <- mid_ID
    }
    mid_ID <- round(mean(gap_ID))      
  }
  
  gap_I_min_max <- gap_I_min[2]
  gap_S <- c(0, S_ub[gap_ID[1]])
  mid_S <- mean(gap_S)
  while(abs(diff(gap_S)) > tol_root){
    gap_I_min <- gap_I_min_max - gap_S
    mid_I_min <- gap_I_min_max - mid_S
    
    gap_rent <- gap_I_min + I_limit
    mid_rent <- mid_I_min + I_limit
    
    if(mid_rent * gap_rent[1] <= 0){
      gap_S[2] <- mid_S
    }else{
      gap_S[1] <- mid_S
    }
    mid_S <- mean(gap_S)      
  }
  
  for(node_iter in 1:length(S_ub_temp)){
    if(node_iter < gap_ID[1]){
      S_ub_temp[node_iter] <- 0
    }else if(node_iter == gap_ID[1]){
      S_ub_temp[node_iter] <- mid_S
      break
    }
  }
  
  S_ub_range[1] <- gap_ID[1] - mid_S / S_ub[gap_ID[1]]  
}

# Check optimal position last
V_ub_temp <- green_V %*% S_ub_temp
V_ub_exceed <- max(V_ub_temp) - V_limit > 0
if(!V_ub_exceed){
  V_ub_gap <- V_limit - max(V_ub_temp) 
  while(V_ub_gap > 0){
    low_ID <- as.integer(S_ub_range[1]) + (as.integer(S_ub_range[1]) != S_ub_range[1])
    high_ID <- as.integer(S_ub_range[2]) + (as.integer(S_ub_range[2]) != S_ub_range[2])
    
    exchange_diff <- green_V[num_node - 1, high_ID] - green_V[num_node - 1, low_ID]
    S_move <- min(S_ub_temp[low_ID], S_ub[high_ID] - S_ub_temp[high_ID])
    S_move <- min(S_move, V_ub_gap / exchange_diff)
    
    S_ub_temp[low_ID] <- S_ub_temp[low_ID] - S_move
    S_ub_temp[high_ID] <- S_ub_temp[high_ID] + S_move
    V_ub_gap <- V_ub_gap - exchange_diff * S_move
    if(S_ub_temp[low_ID] == 0){
      S_ub_range[1] <- as.integer(S_ub_range[1]) + 1
      S_ub_range[2] <- high_ID + S_ub_temp[high_ID] / S_ub[high_ID]
    }else if(S_ub_temp[high_ID] == 0){
      S_ub_range[2] <- as.integer(S_ub_range[2])
      S_ub_range[1] <- low_ID + 1 - S_ub_temp[low_ID] / S_ub[low_ID]
    }
  }
}
S_ub <- S_ub_temp
V_limit_control[2] <- (green_V %*% S_ub)[control_node - 1]