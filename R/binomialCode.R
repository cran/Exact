binomialCode <-
function(data,alternative="two.sided",npNumbers=100,beta=0.001,interval=FALSE,method="Z-pooled",cond.row=TRUE,to.plot=TRUE,ref.pvalue=TRUE){

#If conditioning on row, then transpose 2x2 table
if(cond.row){data<-t(data)}

Ns<-.colSums(data,2,2)
N<-sum(Ns)

#Specify nuisance parameter range
if(interval){
   if(sum(data[1,])==0){
      int<-seq(0.00001,((sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,]))))/
      (N-sum(data[1,])+(sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,])))),length=npNumbers)
   } else if(sum(data[1,])==N) {
      int<-seq(sum(data[1,])/(sum(data[1,])+(N-sum(data[1,])+1)*
      qf(1-beta/2,2*(N-sum(data[1,])+1),2*sum(data[1,]))),0.99999,length=npNumbers)
   } else {
      int<-seq(sum(data[1,])/(sum(data[1,])+(N-sum(data[1,])+1)*
      qf(1-beta/2,2*(N-sum(data[1,])+1),2*sum(data[1,]))),
      ((sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,]))))/
      (N-sum(data[1,])+(sum(data[1,])+1)*qf(1-beta/2,2*(sum(data[1,])+1),2*(N-sum(data[1,])))),
      length=npNumbers)
   }
} else {int<-seq(0.00001,.99999,length=npNumbers); beta<-0}


x <- rep(0:Ns[1], each=(Ns[2]+1))
y <- rep.int(0:Ns[2], Ns[1]+1)
p1 <- x/Ns[1]
p2 <- y/Ns[2]

#Calculate the test statistics:
if(tolower(method)=="boschloo"){
   TX<-matrix(c(x,y,apply(matrix(c(x,Ns[1]-x,y,Ns[2]-y),(Ns[1]+1)*(Ns[2]+1),4),1,
                          FUN=function(tbls){fisher.2x2(matrix(tbls,2,2),alternative=tolower(alternative))})),
              (Ns[1]+1)*(Ns[2]+1),3)
}

if(tolower(method) %in% c("z-pooled","pooled","score")){
   TX<-matrix(c(x,y,(p1-p2)/sqrt(((x+y)/N)*(1-((x+y)/N))*sum(1/Ns))),(Ns[1]+1)*(Ns[2]+1),3)
}

if(tolower(method) %in% c("z-unpooled","unpooled","wald")){
   TX<-matrix(c(x,y,(p1-p2)/sqrt(p2*(1-p2)/Ns[2]+(p1)*(1-p1)/Ns[1])),(Ns[1]+1)*(Ns[2]+1),3)
}

if(tolower(method) %in% c("santner and snell","santner","snell")){
   TX<-matrix(c(x,y,(p1-p2)),(Ns[1]+1)*(Ns[2]+1),3)
}

#Observed Statistic:
TX[which(is.na(TX))]=0
TXO<-TX[(Ns[2]+1)*data[1,1]+(data[1,2]+1),3]

#Find tables that have a test statistic as or more extreme than the observed statistic:
if(tolower(method) != "boschloo"){
   if(tolower(alternative)=="greater"){Tbls<-TX[which(TX[,3]>=TXO),]
   } else if(tolower(alternative)=="less"){Tbls<-TX[which(TX[,3]<=TXO),]
   } else if(tolower(alternative)=="two.sided"){Tbls<-TX[which(abs(TX[,3])>=abs(TXO)),]}
} else {Tbls<-TX[which(TX[,3]<=TXO),]}


#Search for the maximum p-value:
#Instead of calculating the binomial probabilities several times for a more extreme cell,
#calculate the probability once and then combine with other more extreme cells:
if(!is.matrix(Tbls)){Tbls<-t(Tbls)}
xTbls<-Tbls[,1]
yTbls<-Tbls[,2]
nTbls<-length(xTbls) # Number of 'as or more extreme' tables

x.unique<-unique(xTbls)
y.unique<-unique(yTbls)
lxu<-length(x.unique)
lyu<-length(y.unique)

xnr<-max(xTbls)+1
A<-matrix(nrow=xnr, ncol=npNumbers)
cellnr<-rep.int(x.unique+1, npNumbers)+xnr*rep(seq(npNumbers)-1, each=lxu)
A[cellnr]<-dbinom(rep.int(x.unique, npNumbers), Ns[1], rep(int, each=lxu))

ynr<-max(yTbls)+1
B<-matrix(nrow=ynr, ncol=npNumbers)
cellnr<-rep.int(y.unique+1, npNumbers)+ynr*rep(seq(npNumbers)-1, each=lyu)
B[cellnr]<-dbinom(rep.int(y.unique, npNumbers), Ns[2], rep(int, each=lyu))

prob<-A[xTbls+1,]*B[yTbls+1,]
prob<-.colSums(prob, nTbls, npNumbers)

np<-int[which(prob==max(prob))]
pvalue<-max(prob)+beta


#Refine the p-value using the optimise function
if(ref.pvalue==TRUE){
   refPvalue<-rep(0,length(np))
   refNp<-rep(0,length(np))
   for(i in 1:length(np)){
      ref<-optimise(f=function(p){sum(dbinom(Tbls[,1],Ns[1],p)*dbinom(Tbls[,2],Ns[2],p))},
      interval=c(max(0.00001,np[i]-1/npNumbers),min(0.99999,np[i]+1/npNumbers)), maximum=TRUE)
      refPvalue[i]<-ref$objective+beta
      refNp[i]<-ref$maximum
   }
   np<-refNp[refPvalue==max(refPvalue)]
   pvalue<-refPvalue[refPvalue==max(refPvalue)]
}

#Plot p-value vs np
if(to.plot==TRUE){
   plot(int,prob+beta,xlim=c(floor(min(int)*10)/10,ceiling(max(int)*10)/10),
   ylim=c(0,max(pvalue)),xlab="np",ylab="P-value",main="P-value as a function of the nuisance parameter")
   points(np,pvalue,col="red",pch=21,bg="red")
}

list<-list(method=method,p.value=pvalue,test.statistic=TXO,np=np,np.range=c(min(int),max(int)))
return(list)

}
