
###################################################################################################################################################################################
###EACH BRAND MAXIMIZES THE SUM OF POSTERIOR EXPECTED PROFITS OVER ALL PRODUCTS IT IS OFFERING AT THE MARKET!######################################################################
###IN THIS VERSION, BRAND i MAXIMIZES THE SUM OF POSTERIOR EXPECTED PROFITS WITH RESPECT TO p_j GIVEN p_{-j} FOR ALL j \in######################################################### 
###[1,...,J] PRODUCTS OFFERED BY BRAND i###########################################################################################################################################
###################################################################################################################################################################################
###IMPLEMENTATION:
###A grp VECTOR INDICATES THE SET OF PRODUCTS BELONGING TO THE SAME BRAND AND IS INPUT TO "nashEquilibrium"-FUNCTION
###THE grp-VECTOR IS GIVEN AS FOLLOWS FOR THE CHIPS EXAMPLE
#products = c("NikNaks 20g","Simba 36g","Crack-A Snack 50g","Lays 36g","Doritos 45g","Crack-A Snack 24g","Spookies 50g","Flyers Puffed Corn  20g","Unbranded Extrudes 20g","Bigga Naks 24g","Mamas Puffed Corn 50g","Fritos 25g","Snappers 22g","Stylos 24g","Snak Nak 20g","Jiggies 22g","Frimax Thunder 20g","Unbranded Extrudes 50g","Frimax Potato chips 30g","NikNaks 55g","Cheetos 16g","Fritos 65g","Simba 45g","Lays 40g","Doritos 55g","Simba 25g","Lays 23g","Doritos 30g")
#sort(products)
#brands = c("NikNaks","Simba","Crack-A Snack","Lays","Doritos","Spookies","Flyers Puffed Corn","Unbranded Extrudes","Bigga Naks","Mamas Puffed Corn","Fritos","Snappers","Stylos","Snak Nak","Jiggies","Frimax","Cheetos")
#grp = c(1,2,3,4,5,3,6,7,8,9,10,11,12,13,14,15,16,8,15,1,17,11,2,4,5,2,4,5)
#length(grp) == length(products)
#brand_ind = unique(grp)
#length(brand_ind) == length(brands)
#names(brand_ind) = brands
#brand_ind

### Begin: Functions ###
FOC <- function(prices, design, beta, costs, FC = FALSE) {
  #
  # Compute a vector of the derivatives of outsideGood Goods profits.
  # note: we assume outsideGood goods, with outsideGood-1 single product firms,
  #       outsideGood is the outside good
  #
  # prices is nProducts x 1 vector; the price of the outside good needs to be zero
  # costs is a nProducts x 1 vector; the costs of the outside good need to be zero
  # beta is nObs x nParams matrix of betadraws -- we assume last col is price coefficient
  #

  fullDesign <- cbind(design, prices)
  xbeta <- beta%*%t(fullDesign)

  # Correct for underflow and and overflow
  xbeta <- ifelse(xbeta > 500, 500, xbeta)
  xbeta <- ifelse(xbeta < -500, -500, xbeta)
  # Share of Preference
  probabilities <- exp(xbeta) / rowSums(exp(xbeta))
  # First Choice IF FC = TRUE
  if (FC) {
     probabilities <- t(apply(xbeta, 1, function(x) {(x == max(x))}) * 1)
    }
  insideGoods <- which(prices != 0)
  probsInside <- probabilities[, insideGoods]

  # probabilities is an nObs x nProducts matrix of probabilities
  #   for each good for each beta draw used
  # we only need first nProducts-1 columns of probabilities
  # The FOC are mean(Beta_p * pi_i * (1 - pi_i) * (p_i - c_i) + pi_i)

  # Handle Monopoly Case where only one inside good
  probsInside <- as.matrix(probsInside)
  ###FOC: DERIVATIVE OF EXPECTED PROFITS W.R.T. p_j
  delprofit = (colMeans(probsInside * (1-probsInside) * beta[,ncol(beta)])) * (rowSums(prices-costs)[insideGoods]) + colMeans(probsInside)

  return(delprofit)
}

Price_Elasticity <- function(design, prices, beta) {

  fullDesign <- cbind(prices,design)
  xbeta <- beta%*%t(fullDesign)
  #Stabilize
  max_xbeta = apply(xbeta,1,max)
  max_xbeta_big = array(max_xbeta,dim=c(dim(beta)[1],dim(fullDesign)[1])) #copy maximum p times
  max_xbeta_stab = xbeta - max_xbeta_big 
  #Share of Preference
  probabilities <- exp(max_xbeta_stab) / rowSums(exp(max_xbeta_stab))
  insideGoods <- which(prices != 0)
  probsInside <- probabilities[, insideGoods]
  #Monopolistic case where only one inside good
  probsInside <- as.matrix(probsInside)
  ###FOC: DERIVATIVE OF EXPECTED PROFITS W.R.T. p_j
  #elasticity = colMeans(((probsInside * (1-probsInside) * beta[,1])*prices[prices!=0])/probsInside)
  elasticity = colMeans((probsInside * (1-probsInside) * beta[,1]))
  return(elasticity)
}

computeShares <- function(prices, beta, design, pw = FALSE) {
	# Compute the share of preference for the specifified design
	# the price of the outside good needs to be zero
	# beta is a nObs X nParams matrix of draws of beta --
  #     The last column should be the price coefficient
  prices = prices
  if (pw == FALSE){
    fullDesign <- cbind(prices,design)
  }else{
    fullDesign = design
  }
  probabilities_log_cpp(beta,fullDesign)
}

computeShares_BC <- function(prices, beta, design, pr = 1) {
  # Compute the share of preference for the specifified design
  # the price of the outside good needs to be zero
  # beta is a nObs X nParams matrix of draws of beta --
  #     The last column should be the price coefficient
  
  fullDesign <- cbind(prices,design)
  probabilities_BC_log_cpp(beta,fullDesign,pr)
}


computeProfit <- function(design, prices, beta, costs, FC = FALSE, outside = FALSE) {
	# Compute the expected profit for a given design and price vector
	# prices should include the price of the outside good which should be 0.
	# costs should also include the outside good  as 0 in the last element.
  # beta is an ndraws by nparam array of draws for approximately integrals.
  # the outside good NEEDS TO BE LAST in the design matrix; thomas July 4, 2018 (see computeProfit...)

  fullDesign <- cbind(prices,design)
  probabilities = probabilities_log_cpp(beta,fullDesign)
	net_price = (prices - costs)
  if (outside) {
    # net_price = c(net_price[net_price!=0],0) ###Only guaranteed to work if costs are << min(price): dimension reduction otherwise
    net_price[length(net_price)] = 0  # a workaround the warnings generated when prices minus costs is zero for an inside alternative;  
    # need to rethink the "outside" flag, thomas, July 3, 2018 
  }else{
    net_price = net_price[net_price!=0]
  }
  profits <- as.vector(probabilities)  * net_price ###as.vector is NOT really relevant here
	list(profits = profits, demand = probabilities)
}


computeProfit_BC <- function(design, prices, beta, costs, FC = FALSE, outside = FALSE, pr = 1) {
  #pr indicates the column containing prices in the design; pr = 1 in Nash examples
  # the outside good NEEDS TO BE LAST in the design matrix; thomas July 4, 2018 (see computeProfit...)
  #print(prices) # for monitoring
  fullDesign <- cbind(prices,design)
  probabilities = probabilities_BC_log_cpp(beta,fullDesign,pr)
  net_price = (prices - costs)
  if (outside) {
    #net_price = c(net_price[net_price!=0],0)
    net_price[length(net_price)] = 0  # a workaround the warnings generated when prices minus costs is zero for an inside alternative;  
    # need to rethink the "outside" flag, thomas, July 3, 2018 
  }else{
    net_price = net_price[net_price!=0]
  }
  #if (length(as.vector(probabilities)) != length(net_price)) browser()
  profits <- as.vector(probabilities)  * net_price
  list(profits = profits, demand = probabilities)
}


nashEquilibrium <- function(design, startingPrices, costs, beta, priceRange,grp, FC = FALSE, outside = FALSE, freeze, BC = FALSE, pr = 1) {
	# Find the Nash equilibrium prices for a specific set of betas
	# This uses the best response method where each firm prices at the best
	# response given the other firms prices and offerings.
  # grp as indicator which product belongs to which brand
  # freeze is a vector of integer pointers to inside brands that do not price-compete
  # the outside good NEEDS TO BE LAST in the design matrix; thomas July 4, 2018 (see computeProfit...)
  # ALWAYS start with positive prices for inside goods
	prices <- as.matrix(startingPrices)
	numIters <- 0
	#path_eq_prices = NULL
	valid = FALSE
	# Convergence parameters
	tol = .01
	maxIter = 1000

	while(!valid){
    #if(numIters != 0){
      #cat("Completed Round", numIters, "\n")
    #}
		numIters <- numIters + 1
		oldPrices <- prices
		
    for(i in 1:(length(which(prices != 0)))) {
      test = (i == freeze)
      if(any(test)){
        prices[i, prices[i,] != 0] = oldPrices[i, prices[i,] != 0]
      }else{
        if(BC == TRUE){
		      cp <- function(p) {
		        prices[i, prices[i,] != 0] <- p
		        sum(computeProfit_BC(design, prices, beta, costs, FC, outside,pr)$profits[which(grp == grp[i])])  # added grouping list  rm 150227
		      }
        }else{
          cp <- function(p) {
            prices[i, prices[i,] != 0] <- p
            sum(computeProfit(design, prices, beta, costs, FC, outside)$profits[which(grp == grp[i])])  # added grouping list  rm 150227
          }
        }
		    prices[i, prices[i,] != 0] <- optimize(cp, interval=priceRange[i,], maximum=TRUE)$maximum
      }
    }
		print('iteration')
		print(numIters)
		print('prices and profits')
		print(cbind(prices,computeProfit_BC(design, prices, beta, costs, FC, outside,pr)$profits))
		if(sum(abs(prices - oldPrices)) < tol || numIters == maxIter) {
			  valid = TRUE
    }
	}
	return(list(eq_prices = prices[prices != 0], iterations = numIters))
}


computeNashDistribution <- function(draws, design, prices, costs, priceRange, grp, FC = FALSE) {
	# Compute the posterior distribution of equilibrium prices and associated shares
	# Draws is a list of length nResp of nDrawsPerIter by nParams matrix of draws per respondent
	# design is a nProducts by nParams mdesign matrix - The last row should be the ouside good
	# prices is a nProducts by nPriceAttributes array of starting prices -
  #   The outside good is assumed to have a 0 price
	# costs is a nProducts by nPriceAttributes array of costs - The outside good is assumed to have a 0 cost

	                ### nIters <- dim(draws)[1]     
  length_draws <- length(draws)                            # changed to list and nIter changed to nResp rm 150206

	                ### eqPrices <- array(double(dim(prices)[1]*dim(prices)[2]*nIters), dim=c(nIters, dim(prices)[1], dim(prices)[2]))
  ### eqPrices <- vector("list", length=nResp)          # changed to list rm 150206
  ### eqShares <- matrix(double(dim(prices)[1]*nResp), nrow=nResp)
  eqPrices <- array(0,dim=c(dim(draws[[1]])[1],dim(priceRange)[1])) 
  
	for(i in 1:length_draws) {
		              ### eqPrices[i,,] <- nashEquilibrium(design = design, startingPrices = prices, costs = costs,  beta = draws[i,,], priceRange = priceRange, FC)
    eqPrices[i,] <- nashEquilibrium(design = design, startingPrices = prices, costs = costs,  beta = as.matrix(draws[[i]]), priceRange = priceRange, FC)$eq_prices                        # beta = draws[[i]]; changed to list rm 150206
		### eqShares[i,] <- computeShares(eqPrices[[i]], beta = as.matrix(draws[[i]]), design = design, FC)  # eqPrices[[i]], beta = draws[[i]]; changed to list rm 150206
		if(i %% 10 == 0) {
			cat("Completed Draws ", i, "\n")
		}
	}
	### return(list(eqShares = eqShares, eqPrices = eqPrices))
  return(list(eqPrices = eqPrices))
}


computeOptimumNoCompetitiveReaction <- function(draws, design, prices, costs, targetProduct, FC = FALSE) {
	# Compute the posterior distribution of optimal price and associated shares for targetProduct
	# assuming that the remaining product prices remain fixed
	# Draws is a list of length nResp of nDrawsPerIter by nParams matrix of draws per respondent
	# design is a nProducts by nParams design matrix - The last rouw should be the outside good
	# prices is a nProducts by nPriceAttributes array of starting prices -
  #         The outside good is assumed to have a 0 price
	# costs is a nProduct by nPriceAttributes array of costs - The outside good is assumed to have a 0 cost
	# targetProduct is an interger specifying which product should be optimized

	oneGoodCriterion <- function(targetPrice, basePrices, costs, beta, design, targetProduct) {
		# Instead of minimizing the sum of First Order Conditions, we are simply minimize the target FOC
		# Returns a single value which is the Squared FOC for the targetBrand

		price = basePrices
		price[targetProduct] = targetPrice
		(FOC(c(price), design, beta, costs, FC)*10000)[targetProduct]
	}
  
	ldraws <- length(draws)
	finalPrice = double(ldraws)
	#finalShares = matrix(double(length(prices)*nResp), nrow=nResp)

	for(i in 1:ldraws) {
		finalPrice[i] = uniroot(oneGoodCriterion, interval = c(0, 1e4), basePrices = prices,
								costs = costs, beta = as.matrix(draws[[i]]), design = design,
								targetProduct = targetProduct)$root
		#finalShares[i,] = computeShares(prices = c(finalPrice[i], prices[-1]),
										#beta = draws[i,,], design = design, FC)
	}

	return(list(finalPrice = finalPrice))
}
#### End: Functions



