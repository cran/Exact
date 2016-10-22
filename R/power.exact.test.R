power.exact.test <-
function(p1, p2, n1, n2, npNumbers=100, alpha=0.05,
                             alternative=c("two.sided", "less", "greater"), interval=FALSE, beta=0.001,
                             method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "csm approximate", "csm modified", "fisher"),
                             ref.pvalue=TRUE, simulation=FALSE, nsim = 100){
  
  stopifnot(is.logical(interval) && is.logical(ref.pvalue))
  if(alpha < 0 || alpha >= 0.5){stop("To improve code efficiency, alpha must be between 0 and 0.5")}
  if(interval && (beta < 0 || beta > 1)){stop("Beta must be between 0 and 1")}
  if(npNumbers < 1){stop("Total number of nuisance parameters considered must be at least 1")}
  if(nsim < 1){stop("Need at least one simulation")}
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  
  #Sometimes Z-pooled is called score and Z-unpooled is called wald statistic
  if(length(method)==1 && tolower(method)=="score"){method <- "z-pooled"}
  if(length(method)==1 && tolower(method)=="wald"){method <- "z-unpooled"}
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
                                         "csm", "csm approximate", "csm modified", "fisher"))
    
  if(p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1){stop("Probabilities must be between 0 and 1")}
  if(n1 <= 0 || n2 <= 0){stop("Fixed sample sizes must be greater than 0")}

  if(!simulation){
    prob <- matrix(0, n1+1, n2+1);
    #Consider all tables:
    for(i in 0:n1){
      for(j in 0:n2){
        tables <- matrix(c(i, n1-i, j, n2-j), 2, 2, byrow=TRUE)
        if((alternative=="greater" && tables[1,1]/n1 > tables[2,1]/n2) ||
             (alternative=="less" && tables[1,1]/n1 < tables[2,1]/n2) || 
             (alternative=="two.sided" && tables[1,1]/n1 != tables[2,1]/n2)){
          if(method=="fisher"){
            if(fisher.2x2(tables, alternative=alternative) < alpha){
              prob[i+1, j+1] <- dbinom(i, n1, p1)*dbinom(j, n2, p2)
            }
          } else {
            if(exact.test(tables, npNumbers=npNumbers, alternative=alternative, interval=interval,
                          beta=beta, method=method, to.plot=FALSE, ref.pvalue=ref.pvalue)$p.value < alpha){
              prob[i+1, j+1] <- dbinom(i, n1, p1)*dbinom(j, n2, p2)
            }
          }
        }
      }
      power <- sum(prob)
    }
  }
  if(simulation){
    #Randomly generate a table based on known proportions
    randA <- rbinom(nsim, size = n1, prob = p1)
    randC <- rbinom(nsim, size = n2, prob = p2)
    randTables <- cbind(randA, n1 - randA, randC, n2 - randC)
    p.value <- rep(1, nsim)
    for (i in 1:nsim){
      if((alternative=="greater" && randTables[i,1]/n1 > randTables[i,3]/n2) ||
           (alternative=="less" && randTables[i,1]/n1 < randTables[i,3]/n2) || 
           (alternative=="two.sided" && randTables[i,1]/n1 != randTables[i,3]/n2)){
        if(method=="fisher"){
          p.value[i] <- fisher.2x2(matrix(randTables[i,], 2, 2, byrow=TRUE), alternative=alternative)
        } else {
          p.value[i] <- exact.test(matrix(randTables[i,], 2, 2, byrow=TRUE), npNumbers=npNumbers, alternative=alternative,
                                   interval=interval, beta=beta, method=method, to.plot=FALSE, ref.pvalue=ref.pvalue)$p.value
        }
      }
    }
    power <- mean(p.value < alpha)
  }
  list(power=power, alternative=alternative, method=method)
}
