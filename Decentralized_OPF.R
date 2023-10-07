### Decentralized Optimization 
## DC Power Flow Assumption, ignoring reactive power flow

### Network Parameters
num_node <- 10
num_line <- 9
z_img <- .001 
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
Y_eq <- t(NL) %*% Cond %*% NL
Z_eq <- solve(Y_eq[2:num_node, 2:num_node])

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

### Main optimization process
## Sensitivity matrix
V_sensitivity <- array(0, dim = c(num_node, num_node, num_node - 1))
I_sensitivity <- array(0, dim = c(num_node, num_node, num_line))
for(node_iter_1 in 1:(num_node - 1)){
  for(node_iter_2 in (node_iter_1 + 1):num_node){
    dS = rep(0, num_node)
    dS[node_iter_1] <- 1
    dS[node_iter_2] <- -1    
    V_sensitivity[node_iter_1, node_iter_2, ] <- Z_eq %*% dS[2:num_node]
    V_sensitivity[node_iter_2, node_iter_1, ] <- -V_sensitivity[node_iter_1, node_iter_2, ]
    I_sensitivity[node_iter_1, node_iter_2, ] <- Cond %*% NL %*% c(0, Z_eq %*% dS[2:num_node])
    I_sensitivity[node_iter_2, node_iter_1, ] <- -I_sensitivity[node_iter_1, node_iter_2, ]
  }
}

## Update grid rent between trade of each pair of nodes
# (first calculation) first node: injet (sell) power; 2nd node: draw (buy) power
# second value is just the opposite
# dual_voltage <- rep(0, num_node - 1)
dual_voltage <- c(rep(0, num_node - 2), -3670)
dual_current <- rep(0, num_line)
change_v <- rep(1, length(dual_voltage))
change_i <- rep(1, length(dual_current))

for(loop in 1:500){
  rent <- matrix(0, nrow = num_node, ncol = num_node)
  for(node_iter_1 in 1:(num_node - 1)){
    for(node_iter_2 in (node_iter_1 + 1):num_node){
      dS = rep(0, num_node)
      dS[node_iter_1] <- 1
      dS[node_iter_2] <- -1
      
      drent <- t(dual_voltage) %*% V_sensitivity[node_iter_1, node_iter_2, ]
      drent <- drent + t(dual_current) %*% I_sensitivity[node_iter_1, node_iter_2, ]
      rent[node_iter_1, node_iter_2] <- drent
      rent[node_iter_2, node_iter_1] <- -drent
    }
  }
  
  ## P2P Trading from highest price spread
  trade_confirmed <- matrix(0, nrow = num_node, ncol = num_node)
  moc_supply_temp <- moc_supply
  moc_demand_temp <- moc_demand
  
  # loop for P2P trade
  while(TRUE){
    margin_price_supply_ID <- rep(1, num_node)
    for(node_iter in 1:num_node){
      for(price_iter in margin_price_supply_ID[node_iter]:num_price){
        if(moc_supply_temp[price_iter, node_iter] != 0 || price_iter == num_price){
          margin_price_supply_ID[node_iter] <- price_iter
          break
        }
      }
    }
    
    margin_price_demand_ID <- rep(num_price, num_node)
    for(node_iter in 1:num_node){
      for(price_iter in margin_price_demand_ID[node_iter]:1){
        if(moc_demand_temp[price_iter, node_iter] != 0 || price_iter == 1){
          margin_price_demand_ID[node_iter] <- price_iter
          break
        }
      }
    }
    
    price_spread <- matrix(0, nrow = num_node, ncol = num_node)
    price_spread_max <- -Inf
    price_spread_max_pair <- c() # first = sell, second = buy
    for(node_iter_1 in 1:(num_node - 1)){
      for(node_iter_2 in (node_iter_1 + 1):num_node){
        price_spread[node_iter_1, node_iter_2] <- price[margin_price_demand_ID[node_iter_2]] - price[margin_price_supply_ID[node_iter_1]]
        price_spread[node_iter_1, node_iter_2] <- price_spread[node_iter_1, node_iter_2] - rent[node_iter_1, node_iter_2]
        if(price_spread[node_iter_1, node_iter_2] > price_spread_max){
          price_spread_max <- price_spread[node_iter_1, node_iter_2]
          price_spread_max_pair <- c(node_iter_1, node_iter_2)
        }
        
        price_spread[node_iter_2, node_iter_1] <- price[margin_price_demand_ID[node_iter_1]] - price[margin_price_supply_ID[node_iter_2]]
        price_spread[node_iter_2, node_iter_1] <- price_spread[node_iter_2, node_iter_1] - rent[node_iter_2, node_iter_1]
        if(price_spread[node_iter_2, node_iter_1] > price_spread_max){
          price_spread_max <- price_spread[node_iter_2, node_iter_1]
          price_spread_max_pair <- c(node_iter_2, node_iter_1)
        }        
      }
    }
    
    trade_quan <- moc_supply_temp[margin_price_supply_ID[price_spread_max_pair[1]], price_spread_max_pair[1]]
    trade_quan <- min(trade_quan, moc_demand_temp[margin_price_demand_ID[price_spread_max_pair[2]], price_spread_max_pair[2]])
    trade_confirmed[price_spread_max_pair[1], price_spread_max_pair[2]] <- trade_confirmed[price_spread_max_pair[1], price_spread_max_pair[2]] + trade_quan
    moc_supply_temp[margin_price_supply_ID[price_spread_max_pair[1]], price_spread_max_pair[1]] <- moc_supply_temp[margin_price_supply_ID[price_spread_max_pair[1]], price_spread_max_pair[1]] - trade_quan
    moc_demand_temp[margin_price_demand_ID[price_spread_max_pair[2]], price_spread_max_pair[2]] <- moc_demand_temp[margin_price_demand_ID[price_spread_max_pair[2]], price_spread_max_pair[2]] - trade_quan
    
    if(sum(price_spread > 0) == 0){
      break
    }
  }
  
  # Check network constraint violations
  S <- c()
  for(node_iter in 1:num_node){
    S[node_iter] <- sum(trade_confirmed[node_iter, ])
    S[node_iter] <- S[node_iter] - sum(trade_confirmed[, node_iter])
  }
  V <- Z_eq %*% S[2:num_node]
  I <- Cond %*% NL %*% c(0, V)
  
  # Update Lagrangian multipliers
  change_v <- rep(0, length(dual_voltage))
  change_i <- rep(0, length(dual_current))
  
  tol <- 1E-4
  ddV <- 1
  for(dual_v_iter in 1:length(dual_voltage)){
    # overvoltage correction
    if(V[dual_v_iter] > V_limit){
      dual_voltage[dual_v_iter] <- dual_voltage[dual_v_iter] + ddV
      change_v[dual_v_iter] <- 1
    }
    
    # undervoltage correction
    if(V[dual_v_iter] < -V_limit){
      dual_voltage[dual_v_iter] <- dual_voltage[dual_v_iter] - ddV
      change_v[dual_v_iter] <- 1
    }
    
    # dual infeasible correction
    if(abs(V[dual_v_iter]) <= V_limit - tol && abs(dual_voltage[dual_v_iter]) >= tol){
      dual_voltage[dual_v_iter] <- dual_voltage[dual_v_iter] * .99
      change_v[dual_v_iter] <- 1
    }
  }
  
  ddI <- 1
  for(dual_i_iter in 1:length(dual_current)){
    # overcurrent correction
    if(I[dual_i_iter] > I_limit){
      dual_current[dual_i_iter] <- dual_current[dual_i_iter] + ddI
      change_i[dual_i_iter] <- 1
    }
    
    # overcurrent correction (opposite direction)
    if(I[dual_i_iter] < -I_limit){
      dual_current[dual_i_iter] <- dual_current[dual_i_iter] - ddI
      change_i[dual_i_iter] <- 1
    }  
    
    # dual infeasible correction
    if(abs(I[dual_i_iter]) < I_limit - tol && abs(dual_current[dual_i_iter]) > tol){
      dual_current[dual_i_iter] <- dual_current[dual_i_iter] * .99
    }  
  }  
}
