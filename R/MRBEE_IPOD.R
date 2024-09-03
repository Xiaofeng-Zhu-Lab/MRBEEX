#' Mendelian randomization with bias-correction estimating equation: selecting horizontal pleiotropy via IPOD algorithm.
#'
#' Detailed description of the function goes here.
#'
#' @param by A vector (n x 1) of GWAS effect sizes for the outcome.
#' @param bX A matrix (n x p) of GWAS effect sizes for p exposures.
#' @param byse A vector (n x 1) of standard errors for the GWAS effect sizes of the outcome.
#' @param bXse A matrix (n x p) of standard errors for the GWAS effect sizes of the exposures.
#' @param LD A matrix representing the linkage disequilibrium (LD) among instrumental variables. This matrix should be a sparse matrix with lots of zero entries.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix including p exposures and the outcome. Outcome should be the last column.
#' @param cluster.index A vector indicating the cluster membership for each instrumental variable. This is used in standard error estimation.
#' @param Nmin Optional; the minimum sample size for the GWAS if not provided, defaults to the number of instrumental variables.
#' @param tauvec A vector of tuning parameters for penalizing horizontal pleiotropy in the IPOD algorithm.
#' @param max.iter The maximum number of iterations allowed for convergence of the causal effect estimates.
#' @param max.eps The tolerance level for convergence; iteration stops when changes are below this threshold.
#' @param ebic.gamma The penalty factor for extended Bayesian Information Criterion (eBIC) adjustments on pleiotropy.
#' @param reliability.thres A threshold on bias-correction term, defaults to 0.8.
#' @param rho The penalty multiplier used in the ADMM algorithm within the IPOD framework.
#' @param sampling.time The number of subsampling iterations used to estimate the standard error of the causal effect estimate. Defaults to 100. When set to 0, a sandwich formula is applied for the estimation.
#' @param sampling.frac The fraction of the data to be used in each subsampling iteration. Defaults to 0.5, meaning that 50\% of the data is used in each iteration.
#' @param sampling.iter The number of iteration of MRBEE_IPOD to be used in each subsampling iteration. Defaults to 5.
#' @param maxdiff The maximum allowed difference ratio between iterative causal estimates and initial estimations for stabilization.
#' @param theta.ini Initial estimates for the causal effects; defaults to FALSE, indicating automatic initialization.
#' @param gamma.ini Initial estimates for horizontal pleiotropy effects; also defaults to FALSE for automatic setup.
#' @return A list containing detailed results of the analysis, including estimated causal effects, pleiotropy effects, their respective standard errors, and Bayesian Information Criterion (BIC) scores, among other metrics.
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom varbvs varbvs
#' @importFrom Matrix Matrix solve chol
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

MRBEE_IPOD=function(by,bX,byse,bXse,LD="identity",Rxy,cluster.index=c(1:length(by)),Nmin=F,tauvec=seq(3,50,by=5),max.iter=100,max.eps=0.001,ebic.gamma=1,reliability.thres=0.8,rho=2,maxdiff=1.5,sampling.time=100,sampling.frac=0.5,sampling.iter=5,theta.ini=F,gamma.ini=F){
if(LD[1]=="identity"){
A=MRBEE_IPOD_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,Nmin=Nmin,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=rho,maxdiff=maxdiff,sampling.time=sampling.time,sampling.frac=sampling.frac,sampling.iter=sampling.iter,theta.ini=theta.ini,gamma.ini=gamma.ini)
}else{
if(is.vector(bX)==T){
A=MRBEE_IPOD_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Nmin=Nmin,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,maxdiff=maxdiff,theta.ini=theta.ini,gamma.ini=gamma.ini,reliability.thres=reliability.thres)
}else{
########################### Basic information #######################
cat("MRBEE IPOD:\n")
cat("preparing data -> ")
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
if(Nmin==F){
Nmin=m
}
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
bXinv=as.matrix(Theta%*%bX)
tilde.y=as.vector(TC%*%by)
tilde.X=as.matrix(TC%*%bX)
Bt=t(bXinv)
BtB=matrixMultiply(Bt,bX)
Thetarho=solve(LD+rho*diag(m))
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
if(theta.ini[1]==F){
if(length(tilde.y)<2000){
RC=as.matrix(TC%*%LD)
fit0=varbvs(X=RC,Z=tilde.X,y=tilde.y,verbose=F,maxiter=500)
gamma.ini=gamma.ini1=fit0$beta*(fit0$pip>0.8)
theta.ini=theta.ini1=fit0$beta.cov[-1]
}else{
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak")
gamma.ini=gamma.ini1=fit0$gamma
theta.ini=theta.ini1=fit0$theta
}
}else{
gamma.ini=gamma.ini/byse
theta.ini=theta.ini
}
############################## Tuning Parameter ######################
cat("selecting optimal tau in IPOD -> ")
w=length(tauvec)
Btheta=array(0,c(p,w))
Bgamma=array(0,c(m,w))
Bbic=c(1:w)
for(j in length(tauvec):1){
error=1
iter=1
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
Hinv=matrixInverse(BtB-Rxysum[1:p,1:p])
g=matrixVectorMultiply(Bt,as.vector(by-LD%*%gamma))-Rxysum[1:p,p+1]
theta=c(matrixVectorMultiply(Hinv,g))
if((norm(theta,"2")/norm(theta.ini,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini,"2")
}
########################### update gamma ############################
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)

iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta[,j]=theta
Bgamma[,j]=gamma1
df1=sum(gamma1!=0)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
rss=sum(res*as.vector(Theta%*%res))
Bbic[j]=log(rss)*Nmin+(log(Nmin)+ebic.gamma*log(m))*df1
}
######################## Inference #################################
cat("making inference\n")
jstar=which.min(Bbic)
theta=Btheta[,jstar]
gamma=Bgamma[,jstar]
error=1
iter=1
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma!=0)
indvalid=which(gamma==0)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
if(sampling.time==0){
bZ=cbind(bX,LD[,indgamma])
bZinv=Theta%*%bZ
bZ=as.matrix(bZ)
bZinv=as.matrix(bZinv)
H=matrixMultiply(t(bZ),as.matrix(Theta%*%bZ))
H[1:p,1:p]=H[1:p,1:p]-Rxysum[1:p,1:p]
H=matrixInverse(H)
Hat=matrixListProduct(list(bZ,H,t(bZinv)))
Hat=1-diag(Hat)
Hat[Hat<0.5]=0.5
S=LD*0
res=res/Hat
for(j in 1:max(cluster.index)){
a=which(cluster.index==j)
S[a,a]=outer(res[a],res[a])
}
GS=as.matrix(S%*%bZinv)
COV=matrixListProduct(list(H,t(bZinv),GS,H))*m/(m-ncol(bZ))
theta.cov=COV[1:p,1:p]
theta.se=sqrt(diag(theta.cov))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
}else{
ThetaList=matrix(0,sampling.time,p)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
cluster.sampling <- sample(1:max(cluster.index), round(sampling.frac * max(cluster.index)), replace = FALSE)
indj <- which(cluster.index %in% cluster.sampling)
indj <- sort(indj)
LDj <- Matrix(LD[indj, indj], sparse = TRUE)
Thetarhoj <- solve(LDj+rho*diag(length(indj)))
Thetaj <- solve(LDj)
Bt <- as.matrix(t(bX[indj, ]) %*% Thetaj)
BtB <- matrixMultiply(Bt, bX[indj, ])
thetaj=theta
gammaj=gamma1j=gamma
deltaj=0*gammaj
for(jiter in 1:sampling.iter){
indvalidj <- which(gamma1j==0)
indvalidj <- intersect(indvalidj, indj)
Rxysumj <- biasterm(RxyList = RxyList, indvalidj)
Hinv <- solve(BtB - Rxysumj[1:p, 1:p])
g <- matrixVectorMultiply(Bt, by[indj] - as.vector(LD[indj,indj]%*%gammaj[indj])) - Rxysumj[1:p, p + 1]
thetaj <- c(matrixVectorMultiply(Hinv, g))
if((norm(thetaj, "2") / norm(theta.ini, "2")) > maxdiff) {
thetaj <- thetaj / norm(thetaj, "2") * maxdiff * norm(theta.ini, "2")
}
gammaj[indj]=as.vector(Thetarhoj%*%(by[indj]-matrixVectorMultiply(bX[indj, ],thetaj)-deltaj[indj]+rho*gamma1j[indj]))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
ThetaList[j, ] <- thetaj
}
close(pb)
theta.se=colSD(ThetaList)
theta.cov=cov(ThetaList)
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
}
A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini*byse1
A$reliability.adjust=r
}
}
return(A)
}
