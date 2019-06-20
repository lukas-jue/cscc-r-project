
rm(list=ls())

# LOAD LIBRARIES REQUIRED TO CREATE THE SIMULATED DATA. YOU MAY NEED TO INSTALL THESE PACKAGES.
library(xtable)
library(devtools)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(bayesm)
library(ggplot2)
library(tikzDevice)
library(plyr)
library(latex2exp)  ## you may have to install "stringi" too
library(pracma) # required for approximating Jacobian

###Increase memory capacities
memory.limit(size=100000000)

load("Estimation_Data_Beer_20170423.Rdata")
products = c("Amstel Extra Lata 37,5 cl","Amstel Extra Lata 33 cl","Amstel Lata 37,5 cl","Amstel Lata 33 cl","Amstel Cl?sica Lata 33 cl",
             "Cruzcampo Lata 33 cl","Estrella Damm Lata 33 cl","Estrella Galicia Lata 33 cl","Heineken Lata 33 cl","Mahou 5 Estrellas Lata 33 cl",
             "Mahou Cl?sica Lata 33 cl","San Miguel Lata 33 cl","Voll Damm Lata 33 cl","Steinburg (Marca Blanca Mercadona) Lata 33 cl",
             "Marca Blanca Carrefour Lata 33 cl")

#reorder that price is first in the design matrix
for (i in 1:length(E_Data$lgtdata)){
  E_Data$lgtdata[[i]]$X <- E_Data$lgtdata[[i]]$X[,c(16,1:15)]
}

###Initialize lists to store results

#Number of players (in a pricing game)
nplayers = 15

####################################################################################################################
#####################TUNED BC SAMPLER: TRUE MODEL###################################################################
####################################################################################################################

#major difference to before: likelihood function
Rcpp::sourceCpp("rhierMnlRwMixture_rcpp_loop_sim_tuned_BC.cpp",showOutput = FALSE)
source('rhierMnlRwMixture_main_BC.R')
###Market share computations
Rcpp::sourceCpp("Speed++_MS_BC.cpp",showOutput = FALSE)
source('functions_Nash_v7_BC.R')

#number of constrained coefficients (budget & price)
nvar_c = 2
#position of price coefficient in design matrix
pr=1 #cpp notation

###Prior setting
Amu = diag(1/10, nrow = nvar_c, ncol = nvar_c)
mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
nu = 15 + nvar_c
V = nu * diag(nvar_c)*0.5

Prior = list(ncomp=1, Amu = Amu, mustarbarc = mustarbarc, nu = nu, V = V)
Mcmc = list(R=40000, keep=15)#, s=1.6)
#,s=c(0.1,0.5,0.5,0.5)
out_HB = rhierMnlRwMixture_SR(Data=E_Data,Prior=Prior,Mcmc=Mcmc,nvar_c=nvar_c,pr=pr,starting_budget = log(0.74))

betastar_id_HB = out_HB$betadraw
compdraw_HB = out_HB$nmix$compdraw
probdraw_HB = out_HB$nmix$probdraw
rejection = out_HB$rejection
loglike = out_HB$loglike

plot(out_HB$loglike, type="l")

###Get rid of burnin
burnin = 1000
R = dim(betastar_id_HB)[3]

### visually check if it's a stationary time series without burn-in
plot(out_HB$loglike[(burnin+1):R], type="l")

betastar_id_HB = betastar_id_HB[,,(burnin+1):R]
compdraw_HB = compdraw_HB[(burnin+1):R]
probdraw_HB = probdraw_HB[(burnin+1):R]
rejection = rejection[(burnin+1):R,]
loglike = loglike[(burnin+1):R]

R = dim(betastar_id_HB)[3]
N = dim(betastar_id_HB)[1]

###Compute rejection rate of sampler 
rej_rate_indi = apply(rejection,2,mean)
summary(rej_rate_indi)
rej_rate_agg = mean(rej_rate_indi)

###Heterogeneity distribution posterior means
beta_id_HB = betastar_id_HB
beta_id_HB[,1,] = exp(betastar_id_HB[,1,]) # budget
beta_id_HB[,2,] = -exp(betastar_id_HB[,2,]) # price

#Compute Posterior Means
beta_id_PM = apply(beta_id_HB,c(1,2),mean)
summary(beta_id_PM)

###Heterogeneity distribution hierarchical prior
l = 200
index = rep(rank(runif(R)),l)
betastar_id_HP <- array(0,dim=c(R*l,dim(betastar_id_HB)[2]))
#compute posterior predictive density of beta (hierarchical prior)
for(j in 1:(R*l)){
  betastar_id_HP[j,] = rmixture(1,probdraw_HB[index[j]],compdraw_HB[[index[j]]])$x    
}
#transform betastars to betas
beta_id_HP = betastar_id_HP
beta_id_HP[,1] = exp(betastar_id_HP[,1])
beta_id_HP[,2] = -exp(betastar_id_HP[,2])
summary(beta_id_HP)

###Heterogeneity distribution lower level model non-smoothed  
betastar_id_LLMns <- array(0,dim=c(R*l,dim(betastar_id_HB)[2]))
index_r = rep(rank(runif(R)),l)
index_n = rep(rank(runif(N)),round((R*l)/N)+1)

#data generation
for (i in 1:(R*l)){
  betastar_id_LLMns[i,] = betastar_id_HB[index_n[i],,index_r[i]]
}
#transform betastardraws to betadraws
beta_id_LLMns = betastar_id_LLMns
beta_id_LLMns[,1] = exp(betastar_id_LLMns[,1])
beta_id_LLMns[,2] = -exp(betastar_id_LLMns[,2])
summary(beta_id_LLMns)

###Compare posterior predictive densities
summary(beta_id_HP)
summary(beta_id_LLMns)
summary(beta_id_PM)

###Fraction of individuals being budget constrained 
length(which(beta_id_HP[,1]<3))/dim(beta_id_HP)[1]

###################################################
###Convergence of individual level coefficients####
###################################################
#Randomly Select 10 out of 1000 to directly compare recovery of individual level coefficients across T
x = 1:N
index_resp = sample(x, 10, replace = FALSE, prob = NULL)
index_resp

############
###Budget###
############
###Identify maximum price of choice & minimum price not chosen
price_chosen_perTask <- array(0,dim=c(N,length(E_Data$lgtdata[[1]]$y)))
p = dim(E_Data$lgtdata[[1]]$X)[1]/length(E_Data$lgtdata[[1]]$y)
for(indi in 1:N){
  for(task in 1:length(E_Data$lgtdata[[1]]$y)){
    price_chosen_perTask[indi,task] = E_Data$lgtdata[[indi]]$X[((task-1)*p+E_Data$lgtdata[[indi]]$y[task]),pr]
  }
}
max_price_chosen = apply(price_chosen_perTask,1,max)
### doesn't make sense in our case
# min_price_notchosen = array(0,dim=length(max_price_chosen))
# for(i in 1:indi){
#   diff = abs(max_price_chosen[i] - 3)
#   if(diff == 1){
#     min_price_notchosen[i] = 3
#   }else if(diff == 2){
#     min_price_notchosen[i] = 2
#   }else if(diff == 3){
#     min_price_notchosen[i] = 1
#   }else{
#     min_price_notchosen[i] = Inf
#   }
# }
summary(max_price_chosen)
#summary(min_price_notchosen)
###Identify respondents being budget constrained a posteriori###
budget_c_ind = which(max_price_chosen<0.74)
#budget_c_ind = budget_c_ind[105:115]

#Save as eps
#file_name = paste0("Saved_Results/PartII_Basic_Simulation_BC/Conv_betastar_id_HB_budget_TrueBC.eps")
#postscript(file_name, horizontal = FALSE, onefile = FALSE, width=800)
windows()
par(mfrow=c(3,4))  # multiple plots are filled by rows!!!
for (s in 1:length(budget_c_ind)){
  plot(betastar_id_HB[budget_c_ind[s],1,],xlab="R",ylab="betastar_budget",main=budget_c_ind[s],type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
  abline(h = log(max_price_chosen[budget_c_ind[s]]), lwd=5, col = "red")
  #abline(h = betastar[budget_c_ind[s],1], lwd = 5, col = "blue", lty = 2)
}
#dev.off()

###
###Illustrate with only two individuals
###

###Start with indi 10 who is not budget constrained but did not chose highest price
budget_c_ind[10]
E_Data$lgtdata[[budget_c_ind[10]]]
#exp(betastar[budget_c_ind[10],1])

#file_name = paste0("Saved_Results/PartII_Basic_Simulation_BC/Conv_betastar_budget_two_boundaries_not_constrained.eps")
#postscript(file_name, horizontal = FALSE, onefile = FALSE, width=800)
windows()
par(mfrow=c(1,1))  # multiple plots are filled by rows!!!
#for (s in 1:length(budget_c_ind)){
plot(betastar_id_HB[budget_c_ind[10],1,],xlab="R",ylab="betastar_budget",ylim=c(-2,4),main="Individual 10",type="l",cex.main=2,
     cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = log(max_price_chosen[budget_c_ind[10]]), lwd=5, col = "red")
#abline(h = log(min_price_notchosen[budget_c_ind[10]]), lwd=5, col = "red",lty=3)
#abline(h = betastar[budget_c_ind[10],1], lwd = 5, col = "blue", lty = 2)
legend("topright",expression("Log of true budget","Log of maximum price chosen","Log of minimum price not chosen"),cex=1.5,lty=c(2,1,3),lwd=c(5,5,5),
       col=c("blue","red","red"))
#}
#dev.off()

###Do the same for indi 10 who is budget constrained and chose p=1 as highest price
budget_c_ind[1]
E_Data$lgtdata[[budget_c_ind[1]]]
#exp(betastar[budget_c_ind[1],1])

#file_name = paste0("Saved_Results/PartII_Basic_Simulation_BC/Conv_betastar_budget_two_boundaries_constrained.eps")
#postscript(file_name, horizontal = FALSE, onefile = FALSE, width=800)
windows()
par(mfrow=c(1,1))  # multiple plots are filled by rows!!!
#for (s in 1:length(budget_c_ind)){
plot(betastar_id_HB[budget_c_ind[1],1,],xlab="R",ylim=c(-2,4),ylab="betastar_budget",main="Individual 200",type="l",cex.main=2,
     cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h = log(max_price_chosen[budget_c_ind[1]]), lwd=5, col = "red")
#abline(h = log(min_price_notchosen[budget_c_ind[1]]), lwd=5, col = "red",lty=3)
#abline(h = betastar[budget_c_ind[1],1], lwd = 5, col = "blue", lty = 2)
legend("topright",expression("Log of true budget","Log of maximum price chosen","Log of minimum price not chosen"),cex=1.5,lty=c(2,1,3),lwd=c(5,5,5),
       col=c("blue","red","red"))
#}
#dev.off()




#
###Intuition: market shares at highest prices - standard approach vs. population BC
#
designBase = rbind(diag(nplayers),rep(0,nplayers)) #choice set


######################################################################
################APPROXIMATE MARKET DEMAND#############################
######################################################################
#
###Standard model
#
price_grid = seq(0.01,2,0.08) ###increase size of grid for better approximation set 0.01
length(price_grid)
p_other = rep(0.6,15)
#
###Market demand true model (BC)
#
ms_grid_BC_esti <- array(0,dim=c(length(price_grid),16,1)) #second dim is 16 (15 brands + 1 price)

for(r in 1:length(price_grid)){
  price_temp = c(price_grid[r],p_other,0)
  ms_grid_BC_esti[r,,1] = computeShares_BC(price_temp,beta_id_LLMns,designBase,pr=2)
}

###AGGREGATE FIRST BEER BRAND###
windows()
par(mfrow=c(1,1))
plot(price_grid,ms_grid_BC_esti[,1,1],type="l",lty=1,lwd=4,pch=19,col="red",xlim=c(0,2),main=products[1],
     xlab="Price",ylab="Demand",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
lines(price_grid,ms_grid_BC_esti[,1,1],lwd=4,pch=19,col="blue")
#lines(price_grid,ms_grid_BC_true[,1,1],lwd=4,pch=19,col="black")
#lines(price_grid,ms_grid[,2,2],lwd=4,pch=19,col="blue",lty=2)
#lines(price_grid,ms_grid_BC_esti[,2,2],lwd=4,pch=19,col="green",lty=2)
abline(v=0,lwd=1,lty=1)
abline(v=p_other,lwd=4,lty=3,col="black")
legend("topright",expression("Standard model", "Inferred budget","Fixed price"),
       cex=1.2,lty=c(1,1,3),lwd=c(4,4,4),col=c("red","blue","black"))

# drawing price demand curve with ggplot2
# create data frame for ggplot()
market_demand_frame = data.frame(market_demand = ms_grid_BC_esti[,1,1], x = price_grid,
                                 model = rep("Inferred budgets",length(ms_grid_BC_esti[,1,1])))

ggplot(market_demand_frame, aes(x, market_demand)) +
  geom_line() +
  labs(title = products[1],
       subtitle = paste("Price-demand curve, price of other brands = ", "EUR ", p_other),
       x = "Price",
       y = "Demand") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  geom_vline(xintercept = p_other, color = "blue", linetype = "dashed")


#--------------------------------
#not applicable below
#--------------------------------


# ######################################################
# ###Do the same for low brand given fixed high brand###
# ######################################################
# 
# #
# ###Market demand true model (BC)
# #
# for(r in 1:length(price_grid)){
#   price_temp = c(p_other,price_grid[r],0)
#   ms_grid_BC_esti[r,,2] = computeShares_BC(price_temp,beta_id_LLMns,designBase,pr=1)
# }
# 
# ###AGGREGATE DEMAND HIGH BRAND###
# #file_name = paste0("Saved_Results/PartII_Basic_Simulation_BC/Empirical_Demand_Function_BC_PopStandardBCSampler_i_and_iv_h_b.eps")
# #postscript(file_name, horizontal = FALSE, onefile = FALSE, width=800)
# windows()
# par(mfrow=c(1,1))
# plot(price_grid,ms_grid[,1,1],type="l",lty=1,lwd=4,pch=19,col="red",xlim=c(0,10),main="High brand",
#      xlab="Price",ylab="Demand",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
# lines(price_grid,ms_grid_BC_esti[,1,1],lwd=4,pch=19,col="blue")
# #lines(price_grid,ms_grid_BC_true[,1,1],lwd=4,pch=19,col="black")
# #lines(price_grid,ms_grid[,2,2],lwd=4,pch=19,col="blue",lty=2)
# #lines(price_grid,ms_grid_BC_esti[,2,2],lwd=4,pch=19,col="green",lty=2)
# abline(v=0,lwd=1,lty=1)
# abline(v=p_other,lwd=4,lty=3,col="black")
# legend("topright",expression("Standard model", "Inferred budget","Fixed price"),
#        cex=1.2,lty=c(1,1,3),lwd=c(4,4,4),col=c("red","blue","black"))
# #dev.off()
# 
# #############
# ###GGPLOT2###
# #############
# market_demand_frame = data.frame(market_demand = c(ms_grid[,1,1],ms_grid_BC_esti[,1,1]), x = c(price_grid,price_grid),
#                                  model = c(rep("Standard model",length(ms_grid[,1,1])),
#                                            rep("Inferred budgets",length(ms_grid_BC_esti[,1,1]))))
# 
# market_demand_frame_lb = data.frame(market_demand = c(ms_grid[,2,2],ms_grid_BC_esti[,2,2]), x = c(price_grid,price_grid),
#                                     model = c(rep("Standard model",length(ms_grid[,1,1])),
#                                               rep("Inferred budgets",length(ms_grid_BC_esti[,1,1]))))
# 
# windows()
# ggplot(market_demand_frame, aes(x=x, y=market_demand, fill = model, colour = model)) +
#   ggtitle("High brand") +
#   geom_line(size=1.5) +
#   xlim(0,10) +
#   xlab("price") +
#   ylab("demand") +
#   theme(axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15),
#         axis.title.x = element_text(size=18),
#         axis.title.y = element_text(size=18),legend.title = element_blank(),legend.position="none",#legend.key.size = unit(1.2,"cm"),legend.text = element_text(size = 12)
#         plot.title = element_text(size = 18, hjust = 0.5))#,legend.position="top") 
# ###save as pdf & eps
# #ggsave("Saved_Results/PartII_Basic_Simulation_BC/Empirical_Demand_Function_BC_Standard_high_b.pdf",width = 20, height = 20, units = "cm")
# #ggsave("Saved_Results/PartII_Basic_Simulation_BC/Empirical_Demand_Function_BC_Standard_high_b.eps",width = 20, height = 20, units = "cm")
# #dev.off()
# 
# windows()
# ggplot(market_demand_frame_lb, aes(x=x, y=market_demand, fill = model, colour = model)) +
#   ggtitle("Low brand") +
#   geom_line(size=1.5) +
#   xlim(0,10) +
#   ylim(0,1) +
#   xlab("price") +
#   ylab("demand") +
#   theme(axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15),
#         axis.title.x = element_text(size=18),
#         axis.title.y = element_text(size=18),legend.title = element_blank(),legend.key.size = unit(1.2,"cm"),
#         legend.text = element_text(size = 12),plot.title = element_text(size = 18, hjust = 0.5))#,legend.position="top") 
# ###save as pdf & eps
# #ggsave("Saved_Results/PartII_Basic_Simulation_BC/Empirical_Demand_Function_BC_Standard_low_b.pdf",width = 20, height = 20, units = "cm")
# #ggsave("Saved_Results/PartII_Basic_Simulation_BC/Empirical_Demand_Function_BC_Standard_low_b.eps",width = 20, height = 20, units = "cm")
# #dev.off()

