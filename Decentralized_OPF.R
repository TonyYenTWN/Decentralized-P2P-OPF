### Decentralized Optimization 
## DC Power Flow Assumption, ignoring reactive power flow

### Network Parameters
num_node <- 10
num_line <- 10
z_img <- .1 
Cond = diag(rep(1 / z_img, num_line)) 

## Topology of the network
Topology = matrix(nrow = num_line, ncol = 2)
colnames(Topology) = c("From", "To")

for(i in 1:(nrow(Topology)-1)){
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
num_price <- 11
price <- seq(0, 100, length.out = num_price)
moc_supply <- matrix(0, nrow = num_price, ncol = num_node)
rownames(moc_supply) <- price
colnames(moc_supply) <- 1:num_node
for(node_iter in 1:num_node){
  moc_supply[node_iter:num_price, node_iter] <- 10
}
moc_demand <- 10 - moc_supply

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
dual_voltage <- rep(0, num_node - 1)
dual_current <- rep(0, num_line)
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
moc_supply_temp <- moc_supply
margin_price_supply_ID <- rep(1, num_node)
for(node_iter in 1:num_node){
  for(price_iter in margin_price_supply_ID[node_iter]:num_price){
    if(moc_supply_temp[price_iter, node_iter] != 0 || price_iter == num_price){
      margin_price_supply_ID[node_iter] <- price_iter
      break
    }
  }
}

moc_demand_temp <- moc_demand
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
    price_spread[node_iter_1, node_iter_2] <- price[margin_price_demand_ID[node_iter_2]] - price[margin_price_demand_ID[node_iter_1]]
    price_spread[node_iter_1, node_iter_2] <- price_spread[node_iter_1, node_iter_2] - drent[node_iter_1, node_iter_2]
    if(price_spread[node_iter_1, node_iter_2] > price_spread_max){
      price_spread_max <- price_spread[node_iter_1, node_iter_2]

    }

  }
}
