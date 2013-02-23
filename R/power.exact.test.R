power.exact.test <-
function (p1, p2, n1, n2, npNumbers=100, alpha = 0.05,alternative="two.sided",
interval=FALSE,beta=.001,method="Z-pooled",ref.pvalue=TRUE,
simulation=FALSE,nsim = 100){

if(p1<0|p1>1|p2<0|p2>1){stop("Probabilities must be between 0 and 1")}
if(n1 <=0 | n2 <=0){stop("fixed sample sizes must be greater than 0")}
if(alpha < 0 | alpha > 1){stop("Alpha must be between 0 and 1")}
if(nsim < 1){stop("Need at least one simulation")}
if(beta < 0 | beta > 1){stop("Beta must be between 0 and 1")};
if(npNumbers < 1){stop("Total number of nuisance parameters considered must be at least 1")};
if(!(tolower(alternative) %in% c("less","two.sided","greater"))){
stop("Set alternative to 'less', 'two.sided', or 'greater'")}
if(!(tolower(method) %in% c("z-pooled","pooled","score","z-unpooled","unpooled","boschloo","wald","fisher",
"santner and snell","santner","snell"))){
stop("Set method to 'Z-pooled', 'Z-unpooled', 'Boschloo', 'Fisher', or 'santner and snell'")}

if(!simulation){
    prob<-matrix(0,n1+1,n2+1);
    #Consider all tables:
    for( i in 0:n1){
    for( j in 0:n2){
    tables <- matrix(c(i,n1-i,j,n2-j),2,2,byrow=T)
    if(tolower(method)=="fisher"){
       if(fisher.2x2(tables,alternative=tolower(alternative))<alpha){
          prob[i+1,j+1]<-choose(n1,i)*p1^i*(1-p1)^(n1-i)*choose(n2,j)*p2^j*(1-p2)^(n2-j)}}
    if(tolower(method)!="fisher"){
    if(max(binomialCode(tables,npNumbers=npNumbers,alternative=alternative,interval=interval,beta=beta,
method=method,to.plot=F,ref.pvalue=ref.pvalue)$p.value) < alpha){
    prob[i+1,j+1]<-choose(n1,i)*p1^i*(1-p1)^(n1-i)*choose(n2,j)*p2^j*(1-p2)^(n2-j)}}}}
    power<-sum(prob)
}
if(simulation){
    #Randomly generate a table based on known proportions
    randA <- rbinom(nsim, size = n1, prob = p1)
    randC <- rbinom(nsim, size = n2, prob = p2)
    randTables <- cbind(randA, n1 - randA, randC, n2 - randC)
    p.value <- rep(0, nsim)
    for (i in 1:nsim){
    if(tolower(method)=="fisher"){
    p.value[i] <- fisher.2x2(matrix(randTables[i,],2,2,byrow=T),alternative=tolower(alternative))}
    if(tolower(method)!="fisher"){
    p.value[i] <- max(binomialCode(matrix(randTables[i,],2,2,byrow=T),npNumbers=npNumbers,alternative=alternative,
interval=interval,beta=beta,method=method,to.plot=F,ref.pvalue=ref.pvalue)$p.value)}}
    power <- mean(p.value < alpha)
}
list(power=power,alternative=alternative,method=method)
}
