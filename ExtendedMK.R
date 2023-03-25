

############################################
### Function for extended MK ############### 
#### Nikolakopoulos et al,  Extending the ##
#Mann-Kendall test to allow for ############
###  measurement  uncertainty ############## 
### 25-03-2023                ############## 
############################################ 

## Arguments: X: data, a numerical vector of
## observations
## d: Level of relevant difference
## The function returns a vector with elements names:

## S: The S statistic
## S.corr: The S corrected statistic used for inference
## Z: Standardized statistic
## VarS: Variance of S
## ta: tau a
## tb: tau b
## 2-sid p.value: 2-sdied p-value
## ties.p: Percentage of ties in the data, relavant for assessing
## normality assumption


MyMK20 <- function(X,d) {
  
  n <- length(X)
  SS <- Scor<- dim(1)
  i=1
  
  for (k in 1:(n-1)) {
    for (j in (k+1):n) {
      if (abs(X[k]-X[j])>d) {
        SS[i]<-sign(X[j]-X[k]) } else { SS[i]<-0}
      i=i+1            }
  }
  
  S <- sum(SS)
  if (S>0) {Scor<-S-1}
  if (S==0) {Scor<-0}
  if (S<0) {Scor<-S+1}
  
  
  Ui <- Vi <- dim(1) 
  for ( i in 1:n) {
    Ui[i] <- (sum(X[i] > X+d))  
    Vi[i] <- (sum(X[i] < X-d)) }
  
  VarS <- (1/3)*sum((Ui-Vi)^2) + (1/3)*sum(Ui)  
  
  maxS <- n*(n-1)/2
  Z <- Scor/sqrt(VarS)
  ta <- 2*S/(n*(n-1))
  tb <- 2*S/sqrt(2*sum(Ui)*n*(n-1))
  p.value <- 2*(pnorm(abs(Z),lower.tail = F))
  ties.p <- (maxS-sum(Ui))/maxS
  out <- c(S,Scor,Z,VarS,ta,tb,p.value,ties.p)
  names(out) <- c("S","S.corr","Z","VarS","ta","tb","2-sid p.value","ties.p")
  return(out)
  
}