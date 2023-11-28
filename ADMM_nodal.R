### Decentralized optimization using ADMM method, relaxed constraints
## DC Power Flow Assumption, ignoring reactive power flow
## Nodal iteration

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
num_node <- 201
num_line <- num_node - 1
z_img <- 1 / num_line
Cond = rep(1 / z_img, num_line)
V_limit <- pi / 6
I_limit <- .9

## Topology of the network, lines
Topology_line = matrix(nrow = num_line, ncol = 2)
colnames(Topology_line) = c("From", "To")

for(i in 1:nrow(Topology_line)){
  Topology_line[i, ] = c(i, i+1)
}

## Topology of network, nodes
line_inform = list(line_ID = c(), neighbor = c())
node_inform = list(from = c(), to = c(), Y_ii = 0)
Topology_node = list()
for(node_iter in 1:num_node){
  Topology_node[[node_iter]] = node_inform
}
for(line_iter in 1:num_line){
  from_ID = Topology_line[line_iter, 1]
  to_ID = Topology_line[line_iter, 2]
  
  num_line <- length(Topology_node[[from_ID]]$from) + 1
  Topology_node[[from_ID]]$from[[num_line]] <- line_inform
  Topology_node[[from_ID]]$from[[num_line]]$line_ID <- line_iter
  Topology_node[[from_ID]]$from[[num_line]]$neighbor <- to_ID
  Topology_node[[from_ID]]$Y_ii = Topology_node[[from_ID]]$Y_ii + Cond[line_iter]
  
  num_line <- length(Topology_node[[to_ID]]$to) + 1
  Topology_node[[to_ID]]$to[[num_line]] <- line_inform
  Topology_node[[to_ID]]$to[[num_line]]$line_ID <- line_iter
  Topology_node[[to_ID]]$to[[num_line]]$neighbor <- from_ID
  Topology_node[[to_ID]]$Y_ii = Topology_node[[to_ID]]$Y_ii + Cond[line_iter]  
}

### Market Parameters and Cost FUnctions
# Cost functions for S
big_num <- 1E16
num_price <- 100
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

## Main iteration process
## For each node i
## given p_{N_i}, V_{N_i}
## min f(S_i) + (p_{N_i} * I_{N_i} + .5 * rho * (sum(S_!i) + S_i + u)) * S_i
tol <- 1E-6
rho = 10
alpha = .001
price_margin <- rep(0, num_node)
voltage <- rep(0, num_node)
S_prev <- rep(0, num_node)
S <- rep(0, num_node)
u <- 0
for(loop in 1:2000){
  
  # Update power
  constant <- sum(S)
  for(node_iter in 1:num_node){
    constant_temp <- constant - S[node_iter]
    
    # calculate terms from neighbor
    yV_from <- c()
    yp_from <- c()
    for(neighbor_iter in 1:length(Topology_node[[node_iter]]$from)){
      line_ID <- Topology_node[[node_iter]]$from[[neighbor_iter]]$line_ID
      neighbor_ID <- Topology_node[[node_iter]]$from[[neighbor_iter]]$neighbor
      yV_from[neighbor_iter] <- Cond[line_ID] * voltage[neighbor_ID]
      yp_from[neighbor_iter] <- Cond[line_ID] * price_margin[neighbor_ID]
    }
    
    yV_to <- c()
    yp_to <- c()
    for(neighbor_iter in 1:length(Topology_node[[node_iter]]$to)){
      line_ID <- Topology_node[[node_iter]]$to[[neighbor_iter]]$line_ID
      neighbor_ID <- Topology_node[[node_iter]]$to[[neighbor_iter]]$neighbor
      yV_to[neighbor_iter] <- Cond[line_ID] * voltage[neighbor_ID]
      yp_to[neighbor_iter] <- Cond[line_ID] * price_margin[neighbor_ID]
    }  
    
    weighted_price <- (sum(yp_from) + sum(yp_to)) / Topology_node[[node_iter]]$Y_ii
    weighted_voltage <- (sum(yV_from) + sum(yV_to)) / Topology_node[[node_iter]]$Y_ii
    
    lb_S <- -V_limit
    ub_S <- V_limit
    for(neighbor_iter in 1:length(Topology_node[[node_iter]]$from)){
      line_ID <- Topology_node[[node_iter]]$from[[neighbor_iter]]$line_ID
      neighbor_ID <- Topology_node[[node_iter]]$from[[neighbor_iter]]$neighbor
      lb_S <- max(lb_S, -I_limit / Cond[line_ID] + voltage[neighbor_ID])
      ub_S <- min(ub_S, I_limit / Cond[line_ID] + voltage[neighbor_ID])
    }
    
    for(neighbor_iter in 1:length(Topology_node[[node_iter]]$to)){
      line_ID <- Topology_node[[node_iter]]$to[[neighbor_iter]]$line_ID
      neighbor_ID <- Topology_node[[node_iter]]$to[[neighbor_iter]]$neighbor
      lb_S <- max(lb_S, -I_limit / Cond[line_ID] + voltage[neighbor_ID])
      ub_S <- min(ub_S, I_limit / Cond[line_ID] + voltage[neighbor_ID])
    }
    
    lb_S <- Topology_node[[node_iter]]$Y_ii * (lb_S - weighted_voltage)
    ub_S <- Topology_node[[node_iter]]$Y_ii * (ub_S - weighted_voltage)
    
    ### Unconstrained Optimization
    gap_price_ID <- c(1, nrow(moc_rl[[node_iter]]))
    mid_price_ID <- round(mean(gap_price_ID))
    while(diff(gap_price_ID) > 1){
      gap_rent <- moc_rl[[node_iter]][gap_price_ID, 1] - weighted_price + rho * (moc_rl[[node_iter]][gap_price_ID, 2] + constant_temp + u)
      mid_rent <- moc_rl[[node_iter]][mid_price_ID, 1] - weighted_price + rho * (moc_rl[[node_iter]][mid_price_ID, 2] + constant_temp + u)
      
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
      price_margin[node_iter] <- moc_rl[[node_iter]][mid_price_ID, 1]
      S_unconstrained <- -(price_margin[node_iter] -  weighted_price) / rho - (constant_temp + u)
      
    }else{
      # power fixed case
      S_unconstrained <- moc_rl[[node_iter]][mid_price_ID, 2]
    }    
    
    S_constrained <- (1 - alpha) * S[node_iter] + alpha * S_unconstrained
    S_constrained <- max(S_constrained, lb_S)
    S_constrained <- min(S_constrained, ub_S)
    S[node_iter] = S_constrained

    gap_price_ID <- c(1, nrow(moc_rl[[node_iter]]))
    mid_price_ID <- round(mean(gap_price_ID))
    while(diff(gap_price_ID) > 1){
      gap_rent <- moc_rl[[node_iter]][gap_price_ID, 2] - S[node_iter]
      mid_rent <- moc_rl[[node_iter]][mid_price_ID, 2] - S[node_iter]

      if(mid_rent * gap_rent[1] <= 0){
        gap_price_ID[2] <- mid_price_ID
      }else{
        gap_price_ID[1] <- mid_price_ID
      }

      mid_price_ID <- round(mean(gap_price_ID))
    }
    price_margin[node_iter] = moc_rl[[node_iter]][gap_price_ID[1], 1]
    voltage[node_iter] <- S[node_iter] / Topology_node[[node_iter]]$Y_ii + weighted_voltage
  }
  
  # Update dual
  u <- u + sum(S)
  
  # Check convergence
  error <- S - S_prev
  if(max(abs(error), abs(sum(S))) < tol){
    break
  }
  S_prev <- S
}
S
sum(S)
voltage
price_margin