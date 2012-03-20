
#Binomial Model:
binomialCode<-function(data,alternative="less",npNumbers=.01,beta=.001,interval=FALSE,method="Z-pooled",cond.row=TRUE){

#If conditioning on row, then transpose 2x2 table
if(cond.row){data<-t(data)}

Cs<-c(sum(data[,1]),sum(data[,2]))
N<-sum(Cs)

I<-matrix(0,Cs[1]+1,Cs[2]+1)
for( i in 1:(Cs[1]+1)){
	I[i,]<-rep(i,1,Cs[2]+1)-1}

J<-matrix(0,Cs[1]+1,Cs[2]+1)
for( j in 1:(Cs[2]+1)){
	J[,j]<-t(rep(j,1,Cs[1]+1)-1)}

#Calculate all the possible test statistics:
if(tolower(method)=="boschloo"){
TX<-matrix(0,Cs[1]+1,Cs[2]+1)
for( k in 1:(Cs[1]+1)){
for( l in 1:(Cs[2]+1)){
TX[k,l]<-fisher.test(matrix(c(I[k,l],J[k,l],Cs[1]-I[k,l],Cs[2]-J[k,l]),nrow=2),alternative=alternative)$p.value
}}}

if(tolower(method)=="z-pooled"|tolower(method)=="pooled"){
TX<-(I/Cs[1]-J/Cs[2])/sqrt(((I+J)/N)*(1-((I+J)/N))*sum(1/Cs))}

if(tolower(method)=="z-unpooled"|tolower(method)=="unpooled"|tolower(method)=="wald"){
TX<-(I/Cs[1]-J/Cs[2])/sqrt((J/Cs[2])*(1-J/Cs[2])/Cs[2]+(I/Cs[1])*(1-I/Cs[1])/Cs[1])}

TX[which(is.na(TX))]=0

#Observed Statistic:
TXO<-TX[data[1]+1,data[3]+1]

#Find tables that have a test statistic as or more extreme than the observed statistic:
if(tolower(method)!="boschloo"){
if(tolower(alternative)=="greater"){Tbls<-which(TX>=TXO,arr.ind=TRUE)}
if(tolower(alternative)=="less"){Tbls<-which(TX<=TXO,arr.ind=TRUE)}
if(tolower(alternative)=="two.sided"){Tbls<-which(abs(TX)>=abs(TXO),arr.ind=TRUE)}
} else {Tbls<-which(TX<=TXO,arr.ind=TRUE)}

#Specify nuisance parameter range
if(interval){
if(sum(data[1,])==0){
int<-seq(0.00001,((sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,]))))/
(N-sum(data[1,])+(sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,])))),by=npNumbers)
} else if(sum(data[1,])==N) {
int<-seq(sum(data[1,])/(sum(data[1,])+(N-sum(data[1,])+1)*
qf(1-beta/2,2*(N-sum(data[1,])+1),2*sum(data[1,]))),0.99999,by=npNumbers)
} else {
int<-seq(sum(data[1,])/(sum(data[1,])+(N-sum(data[1,])+1)*
qf(1-beta/2,2*(N-sum(data[1,])+1),2*sum(data[1,]))),
((sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,]))))/
(N-sum(data[1,])+(sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,])))),
by=npNumbers)}
} else {int<-seq(0.00001,.99999,by=npNumbers); beta<-0}

#Search for the maximum p-value:
prob<-apply(matrix(choose(Cs[1],Tbls[,1]-1)*choose(Cs[2],Tbls[,2]-1)*
rep(int,each=length(Tbls[,1]))^(Tbls[,1]+Tbls[,2]-2)*
(1-rep(int,each=length(Tbls[,1])))^(N-Tbls[,1]-Tbls[,2]+2),dim(Tbls)[1],length(int)),2,sum)

np<-int[which(prob==max(prob))]
pvalue<-max(prob)+beta

list<-list(method=method,p.value=pvalue,test.statistic=TXO,np=np,np.range=c(min(int),max(int)))
return(list)
}



#Multinomial Model:
multinomialCode<-function(data,alternative="less",npNumbers=.01,beta=.001,interval=FALSE,method="Z-pooled"){

#Function to calculate the test statistic for a given table:
testStatistic<-function(method,i,j,k,alternative){

if(tolower(method)=="z-pooled"|tolower(method)=="pooled"){
TX<-(i/(i+k)-j/(j+(N-i-j-k)))/sqrt((i+j)/N*(1-(i+j)/N)*(1/(i+k)+1/(j+(N-i-j-k))))}

if(tolower(method)=="boschloo"){
TX<-{}
for(l in k){
TX<-c(TX,fisher.test(matrix(c(i,j,l,N-i-j-l),nrow=2),alternative=alternative)$p.value)
}}

if(tolower(method)=="z-unpooled"|tolower(method)=="unpooled"|tolower(method)=="wald"){
TX<-(i/(i+k)-j/(N-i-k))/sqrt(j/(N-i-k)*(1-j/(N-i-k))/(N-i-k)+i/(i+k)*(1-i/(i+k))/(i+k))}

TX[which(is.na(TX))]=0
return(TX)}


N<-sum(data)

#Observed test statistic:
TXO<-testStatistic(method=method,data[1,1],data[1,2],data[2,1],alternative=alternative)

#The p-value calculation function for Multinomial model:
#Since function is symmetric, do not need to consider all values.
MultiProb<-function(p1,p2){
prob<-0
for( i in 0:N){
for( j in 0:(N-i)){
TXcrit<-testStatistic(method=method,i,j,0:(N-i-j),alternative=alternative)

if(alternative=="less" || tolower(method)=="boschloo"){k<-which(!is.na(TXcrit)&TXcrit <= TXO)-1
} else if(alternative=="greater"){k<-which(!is.na(TXcrit)&TXcrit >= TXO)-1
} else if(alternative=="two.sided"){k<-which(!is.na(TXcrit)&abs(TXcrit) >= abs(TXO))-1}

#Calculate probability even for vector of p2
if(length(k)>0){

prob<-prob+colSums(factorial(N)/(factorial(i)*factorial(j)*factorial(k)*factorial(N-i-j-k))*p1^(i+j)*
matrix(rep(p2,each=length(k))^(i+rep(k,length(p2))),length(k),length(p2))*(1-p1)^(N-i-j)*
matrix((1-rep(p2,each=length(k)))^(N-i-rep(k,length(p2))),length(k),length(p2)))
}}}
return(prob)
}

#Specify nuisance parameter range
if(interval){

if(sum(data[1,])==0){
int1<-seq(0.00001,((sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,]))))/
(N-sum(data[1,])+(sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,])))),
by=npNumbers)

} else if(sum(data[1,])==N) {

int1<-seq(sum(data[1,])/(sum(data[1,])+(N-sum(data[1,])+1)*
qf(1-beta/2,2*(N-sum(data[1,])+1),2*sum(data[1,]))),0.99999,by=npNumbers)

} else {

int1<-seq(sum(data[1,])/(sum(data[1,])+(N-sum(data[1,])+1)*
qf(1-beta/2,2*(N-sum(data[1,])+1),2*sum(data[1,]))),
((sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,]))))/
(N-sum(data[1,])+(sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,])))),
by=npNumbers)}


if(sum(data[,1])==0){

int2<-seq(0.00001,((sum(data[,1])+1)*qf(1-beta/2,2*(sum(data[,1])+1),2*(N-sum(data[,1]))))/
(N-sum(data[,1])+(sum(data[,1])+1)*qf(1-beta/2,2*(sum(data[,1])+1),2*(N-sum(data[,1])))),
by=npNumbers)

} else if(sum(data[,1])==N) {

int2<-seq(sum(data[,1])/(sum(data[,1])+(N-sum(data[,1])+1)*
qf(1-beta/2,2*(N-sum(data[,1])+1),2*sum(data[,1]))),0.99999,
by=npNumbers)

} else {

int2<-seq(sum(data[,1])/(sum(data[,1])+(N-sum(data[,1])+1)*
qf(1-beta/2,2*(N-sum(data[,1])+1),2*sum(data[,1]))),
((sum(data[,1])+1)*qf(1-beta/2,2*(sum(data[,1])+1),2*(N-sum(data[,1]))))/
(N-sum(data[,1])+(sum(data[,1])+1)*qf(1-beta/2,2*(sum(data[,1])+1),2*(N-sum(data[,1])))),
by=npNumbers)}

} else {

int1<-seq(0.00001,.99999,by=npNumbers)
int2<-seq(0.00001,.99999,by=npNumbers)
beta=0}


#Search for the maximum p-value (2 nuisance parameters):
maxProb=0; np1=0; np2=0
for(np1s in int1){

if (np1s>max(int2)){int2temp<-int2
} else if(min(int2)<min(int1)){int2temp<-c(seq(min(int2),min(int1),by=npNumbers),seq(np1s,max(int2),by=npNumbers))
} else {int2temp<-seq(np1s,max(int2),by=npNumbers)}

prob<-t(MultiProb(np1s,int2temp))
if(max(prob)>maxProb){
maxProb<-max(prob)
pvalue<-maxProb+beta
np1<-np1s
np2<-int2temp[max.col(prob)]}}

list<-list(method=method,p.value=pvalue,test.statistic=TXO,np1=np1,np2=np2,
np1.range=c(min(int1),max(int1)),np2.range=c(min(int2),max(int2)))
return(list)
}

#Combine Binomial and Multinomial Functions
exact.test<-function(data,alternative="less",npNumbers=.01,beta=.001,interval=FALSE,method="Z-pooled",model="Binomial"){

if(dim(data)[1]!=2 | dim(data)[2]!=2){stop("Input 2x2 table")};
if(length(which(data<0))>0){stop("Can't have negative entries")};
if(length(which(data==0))>2){stop("Can't have more than two 0 entries")};
if(beta < 0 | beta > 1){stop("Beta must be between 0 and 1")};
if(npNumbers <= 0 | npNumbers > 1){stop("Increase in nuisance parameters considered must be between 0 and 1")};
if(!(tolower(alternative) %in% c("less","two.sided","greater"))){
stop("Set alternative to 'less', 'two.sided', or 'greater'")}
if(!(tolower(method) %in% c("z-pooled","pooled","z-unpooled","pooled","boschloo","wald"))){
stop("Set method to 'Z-pooled', 'Z-unpooled', or 'Boschloo'")}
if(!(tolower(model) %in% c("binomial","multinomial"))){
stop("Set model to 'Binomial' or 'Multinomial'")}

if(tolower(model)=="binomial"){
results<-binomialCode(data,alternative=alternative,interval=interval,beta=beta,npNumbers=npNumbers,method=method)
list<-list(model=model,method=results$method,alternative=alternative,p.value=results$p.value,
	     test.statistic=results$test.statistic,np=results$np,np.range=results$np.range)}

if(tolower(model)=="multinomial"){
results<-multinomialCode(data,alternative=alternative,interval=interval,beta=beta,npNumbers=npNumbers,method=method)
list<-list(model=model,method=results$method,alternative=alternative,p.value=results$p.value,
	     test.statistic=results$test.statistic,np1=results$np1,np2=results$np2,
	     np1.range=results$np1.range,np2.range=results$np2.range)}

return(list)
}

#Power Calculations
power.exact.test<-function (p1, p2, n1, n2, npNumbers=.01, alpha = 0.05,alternative="two.sided",
					interval=FALSE,beta=.001,method="Z-pooled",simulation=FALSE,nsim = 100){

if(p1<0|p1>1|p2<0|p2>1){stop("Probabilities must be between 0 and 1")}
if(n1 <=0 | n2 <=0){stop("fixed sample sizes must be greater than 0")}
if(alpha < 0 | alpha > 1){stop("Alpha must be between 0 and 1")}
if(nsim < 1){stop("Need at least one simulation")}
if(beta < 0 | beta > 1){stop("Beta must be between 0 and 1")};
if(npNumbers <= 0 | npNumbers > 1){stop("Increase in nuisance parameters considered must be between 0 and 1")};
if(!(tolower(alternative) %in% c("less","two.sided","greater"))){
stop("Set alternative to 'less', 'two.sided', or 'greater'")}
if(!(tolower(method) %in% c("z-pooled","pooled","z-unpooled","pooled","boschloo","wald","fisher"))){
stop("Set method to 'Z-pooled', 'Z-unpooled', 'Boschloo', or 'Fisher'")}

if(!simulation){
    prob<-matrix(0,n1+1,n2+1);
    #Consider all tables:
    for( i in 0:n1){
    for( j in 0:n2){
    tables <- matrix(c(i,n1-i,j,n2-j),2,2,byrow=T)
    if(tolower(method)=="fisher"){
    if(fisher.test(tables,alternative=alternative)$p.value<alpha){
    prob[i+1,j+1]<-choose(n1,i)*p1^i*(1-p1)^(n1-i)*choose(n2,j)*p2^j*(1-p2)^(n2-j)}}
    if(tolower(method)!="fisher"){
    if(binomialCode(tables,npNumbers=npNumbers,alternative=alternative,
					interval=interval,beta=beta,method=method)$p.value < alpha){
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
    p.value[i] <- fisher.test(matrix(randTables[i,],2,2,byrow=T),alternative=alternative)$p.value}
    if(tolower(method)!="fisher"){
    p.value[i] <- binomialCode(matrix(randTables[i,],2,2,byrow=T),npNumbers=npNumbers,alternative=alternative,
					interval=interval,beta=beta,method=method)$p.value}}
    power <- mean(p.value < alpha)
}
list(power=power,alternative=alternative,method=method)
}



