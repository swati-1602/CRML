library(wooldridge)
library(foreach)
library(doParallel)
numCores<-detectCores()
cl<-makeCluster(numCores-1)
registerDoParallel(cl)
data<-mroz
index<-which(data$hours==0) 

# Count left-censored observations
left_censored_count <-length(index)
censoring<-(left_censored_count/753)*100


### first stage regression
data$res<-residuals(lm(nwifeinc~age+educ+exper+expersq+kidslt6+kidsge6+huseduc,data=data))

f<-function(tau,A){
  
  Ql<-function(Beta){
    
    ## Define quantile loss function
    qL <- function(x, tau) {
      # Calculate quantile loss
      loss <- ifelse(x < 0, (1 - tau) * -x, tau * x)
      
      # Return mean quantile loss
      return(loss)
    }
    
    ##
    L<-qL(A$hours-pmax(0,Beta[1]+Beta[2]*A$age+Beta[3]*A$educ+Beta[4]*A$exper+Beta[5]*A$expersq+Beta[6]*A$kidslt6+Beta[7]*A$kidsge6+Beta[8]*A$nwifeinc+Beta[9]*A$res),tau)
    return(mean(L))
    
  }
  est_m<-optim(c(rep(1,9)),Ql)$par
  return(est_m)
}


## function J
J <- function(tau, alpha=0.2) {
  if (tau > alpha && tau < 1 - alpha) {
    return(1 / (1 - 2 * alpha))
  } else {
    return(0)
  }
} 

### define fJ function in R code

fJ<-function(tau,A)
{
  return(f(tau,A)*J(tau))
}

## Monte Carlo estimate of the integral
## Generate 10^4 samples from uniform 0 to 1 

Ln<-function(A)
{
  U<-runif(1000,0.1,0.9)
  results<-foreach(j = 1:length(U), .combine = "rbind") %dopar%
    {
      J <- function(tau, alpha=0.1) {
        if (tau > alpha && tau < 1 - alpha) {
          return(1 / (1 - 2 * alpha))
        } else {
          return(0)
        }
      } 
      f<-function(tau,A){
        
        Ql<-function(Beta){
          
          ## Define quantile loss function
          qL <- function(x, tau) {
            # Calculate quantile loss
            loss <- ifelse(x < 0, (1 - tau) * -x, tau * x)
            
            # Return mean quantile loss
            return(loss)
          }
          
          ##
          L<-qL(A$hours-pmax(0,Beta[1]+Beta[2]*A$age+Beta[3]*A$educ+Beta[4]*A$exper+Beta[5]*A$expersq+Beta[6]*A$kidslt6+Beta[7]*A$kidsge6+Beta[8]*A$nwifeinc+Beta[9]*A$res),tau)
          return(mean(L))
          
        }
        est_m<-optim(c(rep(1,9)),Ql)$par
        return(est_m)
      }
      
      fJ<-function(tau,A)
      {
        return(f(tau,A)*J(tau))
      }
      
      LL<-0.8*fJ(U[j],A)
    }
  dd1<-matrix(results, nrow = 1000, ncol= 9)
  return(c(mean(dd1[,1]), mean(dd1[,2]),mean(dd1[,3]),mean(dd1[,4]),mean(dd1[,5]),mean(dd1[,6]),mean(dd1[,7]),mean(dd1[,8]),mean(dd1[,9])))
}
Ln(data)

## Bootstrap study 
## Bootstrap samples
B<-seq(50,600,100)

rmse1<-c()
rmse2<-c()
rmse3<-c()
rmse4<-c()
rmse5<-c()
rmse6<-c()
rmse7<-c()
rmse8<-c()
rmse9<-c()

bootstrap_sample <- function(X) {
  n <- nrow(X)
  indices <- sample(1:n, replace = TRUE)
  return(X[indices, ])
}

## number of bootstrap samples
for (b in B) {
  

Boot_samples<- lapply(1:b, function(ii)bootstrap_sample(data))

## Bootstrap L-estimates
 #Est_boot<-foreach(i = 1:b, .combine = c) %dopar%
   #{
   #  Est<- Ln(Boot_samples[[i]])
  #}

Est_boot<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]]) )

# Bootstrapping function

true_integral <- Ln(data)  # Estimate from the full data for comparison



## Bootstrap MSE
l1<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][1,]) )
l2<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][2,]) )
l3<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][3,]) )
l4<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][4,]) )
l5<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][5,]) )
l6<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][6,]) )
l7<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][7,]) )
l8<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][8,]) )
l9<-sapply(1:b, function(ii)Ln(Boot_samples[[ii]][9,]) )

## MSE 
rmse1<-append(rmse1,sqrt(mean((l1 - true_integral[1])^2)),after = length(rmse1))
rmse2<-append(rmse2,sqrt(mean((l2 - true_integral[2])^2)),after = length(rmse2))
rmse3<-append(rmse3,sqrt(mean((l3 - true_integral[3])^2)),after = length(rmse3))
rmse4<-append(rmse4,sqrt(mean((l4 - true_integral[4])^2)),after = length(rmse4))
rmse5<-append(rmse5,sqrt(mean((l5 - true_integral[5])^2)),after = length(rmse5))
rmse6<-append(rmse6,sqrt(mean((l6 - true_integral[6])^2)),after = length(rmse6))
rmse7<-append(rmse7,sqrt(mean((l7 - true_integral[7])^2)),after = length(rmse7))
rmse8<-append(rmse8,sqrt(mean((l8 - true_integral[8])^2)),after = length(rmse8))
rmse9<-append(rmse9,sqrt(mean((l9 - true_integral[9])^2)),after = length(rmse9))
}
rmse1
rmse2
rmse3
rmse4
rmse5
rmse6
rmse7
rmse8
rmse9
