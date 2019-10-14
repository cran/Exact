multinomialCode <-
function(data, alternative, npNumbers, np.interval, beta, method){
  
  #Function to calculate the test statistic for a given table:
  testStatistic <- function(method, i, j, k, alternative){
    
    if(method=="z-pooled"){
      TX <- (i/(i+k)-j/(j+(N-i-j-k)))/sqrt((i+j)/N*(1-(i+j)/N)*(1/(i+k)+1/(j+(N-i-j-k))))
    }
    
    if(method=="boschloo"){
      TX <- {}
      for(l in k){
        TX <- c(TX, fisher.2x2(matrix(c(i,j,l,N-i-j-l),nrow=2), alternative=alternative))
      }
    }
    
    if(method=="z-unpooled"){
      TX <- (i/(i+k)-j/(N-i-k))/sqrt(j/(N-i-k)*(1-j/(N-i-k))/(N-i-k)+i/(i+k)*(1-i/(i+k))/(i+k))
    }
    
    TX[is.na(TX)] <- 0
    TX <- signif(TX, 12) #Remove rounding errors
    return(TX)
  }
  
  N <- sum(data)
  
  #Observed test statistic:
  TXO <- testStatistic(method=method, data[1,1], data[1,2], data[2,1], alternative=alternative)
  
  #The p-value calculation function for Multinomial model:
  #Since function is symmetric, do not need to consider all values.
  MultiProb <- function(p1,p2){
    prob <- 0
    for(i in 0:N){
      for(j in 0:(N-i)){
        TXcrit <- testStatistic(method=method, i, j, 0:(N-i-j), alternative=alternative)
        
        if(alternative=="less" || method=="boschloo"){k <- which(!is.na(TXcrit)&TXcrit <= TXO)-1
        } else if(alternative=="greater"){k <- which(!is.na(TXcrit)&TXcrit >= TXO)-1
        } else if(alternative=="two.sided"){k <- which(!is.na(TXcrit)&abs(TXcrit) >= abs(TXO))-1}
        
        #Calculate probability even for vector of p2
        if(length(k) > 0){
          prob <- prob + colSums(exp(lgamma(N+1)-lgamma(i+1)-lgamma(j+1)-lgamma(k+1)-lgamma(N-i-j-k+1))*p1^(i+j)*
                                   matrix(rep(p2,each=length(k))^(i+rep(k,length(p2))),length(k),length(p2))*(1-p1)^(N-i-j)*
                                   matrix((1-rep(p2,each=length(k)))^(N-i-rep(k,length(p2))),length(k),length(p2)))
        }
      }
    }
    return(prob)
  }
  
  #Specify nuisance parameter range
  if (np.interval) {
    tempInt1 <- binom.CI(sum(data[1,]), N, conf.level=1-beta)
    int1 <- seq(max(c(0.00001, tempInt1[1])),min(c(0.99999, tempInt1[2])), length=npNumbers)
    tempInt2 <- binom.CI(sum(data[,1]), N, conf.level=1-beta)
    int2 <- seq(max(c(0.00001, tempInt2[1])),min(c(0.99999, tempInt2[2])), length=npNumbers)
  } else {
    int1 <- seq(0.00001, .99999, length=npNumbers)
    int2 <- seq(0.00001, .99999, length=npNumbers)
    beta <- 0
  }
  
  #Search for the maximum p-value (2 nuisance parameters):
  maxProb <- np1 <- np2 <- 0
  for(np1s in int1){
    if (np1s > max(int2)){int2temp <- int2
    } else if(min(int2) < min(int1)){int2temp <- c(seq(min(int2),min(int1),length=npNumbers),seq(np1s,max(int2), length=npNumbers))
    } else {int2temp <- seq(np1s,max(int2), length=npNumbers)}
    
    prob <- t(MultiProb(np1s, int2temp))
    if(max(prob) > maxProb){
      maxProb <- max(prob)
      pvalue <- maxProb + beta
      np1 <- np1s
      np2 <- int2temp[max.col(prob)]
    }
  }
  
  return(list(method=method, p.value=max(pvalue), test.statistic=TXO, np1=np1, np2=np2,
              np1.range=c(min(int1),max(int1)), np2.range=c(min(int2), max(int2))))
}
