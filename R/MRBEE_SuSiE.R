#' Mendelian randomization with bias-correction estimating equation: selecting horizontal pleiotropy via IPOD and selecting exposures via SuSiE.
#'
#' Detailed description of the function goes here.
#'
#' @param by A vector (n x 1) of GWAS effect sizes for the outcome.
#' @param bX A matrix (n x p) of GWAS effect sizes for p exposures.
#' @param byse A vector (n x 1) of standard errors for the GWAS effect sizes of the outcome.
#' @param bXse A matrix (n x p) of standard errors for the GWAS effect sizes of the exposures.
#' @param LD A matrix representing the linkage disequilibrium (LD) among instrumental variables.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix including p exposures and the outcome. Outcome should be the last column.
#' @param Lvec The number of single effects used in SuSiE, defaults to c(1:min(5,nrow(bX))).
#' @param pip.thres The threshold of PIP in SuSiE, below which the causal effect estimate will be set to zero, defaults to 0.5.
#' @param cluster.index A vector indicating the cluster membership for each instrumental variable. This is used in standard error estimation.
#' @param Nmin Optional; the minimum sample size for the GWAS if not provided, defaults to the number of instrumental variables.
#' @param tauvec A vector of tuning parameters for penalizing horizontal pleiotropy in the IPOD algorithm.
#' @param max.iter The maximum number of iterations allowed for convergence of the causal effect estimates.
#' @param max.eps The tolerance level for convergence; iteration stops when changes are below this threshold.
#' @param susie.iter The maximum number of iterations allowed for convergence of SuSiE program, defaults to 50.
#' @param ebic.theta The penalty factor for extended Bayesian Information Criterion (eBIC) adjustments on causal effects
#' @param ebic.gamma The penalty factor for eBIC adjustments on pleiotropy.
#' @param empirical.variance.lower The lower boundary of empirical variance estimate of residual, defaults to 0.5.
#' @param empirical.variance.upper The upper boundary of empirical variance estimate of residual, defaults to 2.
#' @param reliability.thres A threshold on bias-correction term, defaults to 0.8.
#' @param rho The penalty multiplier used in the ADMM algorithm within the IPOD framework.
#' @param sampling.time The number of subsampling iterations used to estimate the standard error of the causal effect estimate. Defaults to 100. When set to 0, a sandwich formula is applied for the estimation.
#' @param sampling.frac The fraction of the data to be used in each subsampling iteration. Defaults to 0.5, meaning that 50\% of the data is used in each iteration.
#' @param sampling.iter The number of iteration of MRBEE_SuSiE to be used in each subsampling iteration. Defaults to 5.
#' @param maxdiff The maximum allowed difference ratio between iterative causal estimates and initial estimations for stabilization.
#' @param theta.ini Initial estimates for the causal effects; defaults to FALSE, indicating automatic initialization.
#' @param gamma.ini Initial estimates for horizontal pleiotropy effects; also defaults to FALSE for automatic setup.
#' @return A list containing detailed results of the analysis, including estimated causal effects, pleiotropy effects, their respective standard errors, and Bayesian Information Criterion (BIC) scores, among other metrics.
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom varbvs varbvs
#' @importFrom Matrix Matrix solve chol
#' @importFrom susieR susie_suff_stat coef.susie
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

MRBEE_SuSiE=function(by,bX,byse,bXse,LD="identity",Rxy,cluster.index=c(1:length(by)),Nmin=F,Lvec=c(1:min(10,nrow(bX))),pip.thres=0.5,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,susie.iter=100,ebic.theta=0,ebic.gamma=1,empirical.variance.lower=0.2,empirical.variance.upper=100,reliability.thres=0.8,rho=2,maxdiff=1.5,sampling.time=100,sampling.frac=0.5,sampling.iter=5,theta.ini=F,gamma.ini=F){
if(LD[1]=="identity"){
A=MRBEE_SuSiE_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,Nmin=Nmin,Lvec=Lvec,pip.thres=pip.thres,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,susie.iter=susie.iter,ebic.theta=ebic.theta,ebic.gamma=ebic.gamma,empirical.variance.lower=empirical.variance.lower,empirical.variance.upper=empirical.variance.upper,reliability.thres=reliability.thres,rho=rho,maxdiff=maxdiff,sampling.time=sampling.time,sampling.frac=sampling.frac,sampling.iter=sampling.iter,theta.ini=theta.ini,gamma.ini=gamma.ini)
}else{
if(is.vector(bX)==T){
A=MRBEE_IPOD_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Nmin=Nmin,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,maxdiff=maxdiff,theta.ini=theta.ini,gamma.ini=gamma.ini,reliability.thres=reliability.thres)
}else{
########################### Basic information #######################
cat("MRBEE SuSiE:\n")
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
BtB=(t(BtB)+BtB)/2
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
gamma.ini=gamma.ini1=gamma.ini/byse1
theta.ini=theta.ini1=theta.ini
}
############################## Tuning Parameter ######################
w=length(tauvec)
q=length(Lvec)
Btheta.ipod=array(0,c(p,w))
Bgamma.ipod=array(0,c(m,w))
Btheta.susie=array(0,c(p,q))
Bgamma.susie=array(0,c(m,q))
Bbic.ipod=c(1:w)
Bbic.susie=c(1:q)
cat("selecting optimal tau in IPOD and L in SuSiE -> ")
############################ Selecting Tau in Step 1 ###########################
for(j in length(tauvec):1){
error=1
iter=1
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
Hinv=matrixInverse(BtB-Rxysum[1:p,1:p])
g=matrixVectorMultiply(Bt,by-as.vector(LD%*%gamma1))-Rxysum[1:p,p+1]
theta=c(matrixVectorMultiply(Hinv,g))
if((norm(theta,"2")/norm(theta.ini,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta.ipod[,j]=theta
Bgamma.ipod[,j]=gamma1
df1=sum(gamma1!=0)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
rss=sum(res*as.vector(Theta%*%res))
Bbic.ipod[j]=log(rss)*Nmin+(log(Nmin)+ebic.gamma*log(m))*df1
}
jstar=which.min(Bbic.ipod)
theta.ini=Btheta.ipod[,jstar]
gamma.ini=Bgamma.ipod[,jstar]
indvalid=which(gamma.ini==0)
res=c(by-matrixVectorMultiply(bX,theta.ini)-as.vector(LD%*%gamma.ini))
rss=sum(res*(Theta%*%res))
empirical.variance.formula=min(rss/(m-sum(theta.ini!=0)-sum(gamma.ini!=0)),empirical.variance.upper)
empirical.variance.formula=max(empirical.variance.formula,empirical.variance.lower)
######################## selecting L in Step 1 ###############################
for(v in length(Lvec):1){
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
error=1
iter=1
empirical.variance=empirical.variance.formula
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-as.vector(LD%*%gamma)
XtX=BtB
Xty=matrixVectorMultiply(Bt,res.theta)
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[v],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=F,max_iter=susie.iter)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
delta=delta+rho*(gamma-gamma1)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma1))
rss=sum(res*(Theta%*%res))
empirical.variance=rss/(m-sum(gamma1!=0)-min(sum(theta!=0),Lvec[v]))
empirical.variance=max(empirical.variance,empirical.variance.lower)
empirical.variance=min(empirical.variance,empirical.variance.upper)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta.susie[,v]=theta
Bgamma.susie[,v]=gamma1
df1=sum(gamma1!=0)
df2=min(Lvec[v],sum(theta!=0))
Bbic.susie[v]=log(rss)*Nmin+df2*(log(Nmin)+log(p)*ebic.theta)+(log(Nmin)+ebic.gamma*log(m))*df1
}
vstar=which.min(Bbic.susie)
theta.ini=Btheta.susie[,vstar]
gamma.ini=Bgamma.susie[,vstar]
indvalid=which(gamma.ini==0)
res=c(by-matrixVectorMultiply(bX,theta.ini)-as.vector(LD%*%gamma.ini))
rss=sum(res*(Theta%*%res))
empirical.variance.formula=min(rss/(m-sum(theta.ini!=0)-sum(gamma.ini!=0)),empirical.variance.upper)
empirical.variance.formula=min(empirical.variance.formula,empirical.variance.upper)
######################## selecting L in Step 2 ###############################
for(j in length(tauvec):1){
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
error=1
iter=1
empirical.variance=empirical.variance.formula
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-as.vector(LD%*%gamma)
XtX=BtB
Xty=matrixVectorMultiply(Bt,res.theta)
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[vstar]+1,residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=F,max_iter=susie.iter)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma1))
rss=sum(res*(Theta%*%res))
empirical.variance=rss/(m-sum(gamma1!=0)-min(sum(theta!=0),Lvec[v]))
empirical.variance=max(empirical.variance,empirical.variance.lower)
empirical.variance=min(empirical.variance,empirical.variance.upper)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta.ipod[,j]=theta
Bgamma.ipod[,j]=gamma1
df1=sum(gamma1!=0)
df2=min(Lvec[vstar],sum(theta!=0))
Bbic.ipod[j]=log(rss)*Nmin+df2*(log(Nmin)+log(p)*ebic.theta)+(log(Nmin)+ebic.gamma*log(m))*df1
}
jstar=which.min(Bbic.ipod)
theta.ini=Btheta.ipod[,jstar]
gamma.ini=Bgamma.ipod[,jstar]
indvalid=which(gamma.ini==0)
res=c(by-matrixVectorMultiply(bX,theta.ini)-as.vector(LD%*%gamma.ini))
rss=sum(res*(Theta%*%res))
empirical.variance.formula=min(rss/(m-sum(theta.ini!=0)-sum(gamma.ini!=0)),empirical.variance.upper)
empirical.variance.formula=max(empirical.variance.formula,empirical.variance.lower)
######################## selecting L in Step 2 ###############################
for(v in length(Lvec):1){
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
error=1
iter=1
empirical.variance=empirical.variance.formula
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-as.vector(LD%*%gamma)
XtX=BtB
Xty=matrixVectorMultiply(Bt,res.theta)
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[v],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=F,max_iter=susie.iter)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
########################### update gamma ############################
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
delta=delta+rho*(gamma-gamma1)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma1))
rss=sum(res*(Theta%*%res))
empirical.variance=rss/(m-sum(gamma1!=0)-min(sum(theta!=0),Lvec[v]))
empirical.variance=max(empirical.variance,empirical.variance.lower)
empirical.variance=min(empirical.variance,empirical.variance.upper)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta.susie[,v]=theta
Bgamma.susie[,v]=gamma1
df1=sum(gamma1!=0)
df2=min(Lvec[v],sum(theta!=0))
Bbic.susie[v]=log(rss)*Nmin+df2*(log(Nmin)+log(p)*ebic.theta)+(log(Nmin)+ebic.gamma*log(m))*df1
}
######################## Final Estimate #################################
vstar=which.min(Bbic.susie)
error=1
iter=1
theta.ini=Btheta.susie[,vstar]
gamma.ini=Bgamma.susie[,vstar]
indvalid=which(gamma.ini==0)
res=c(by-matrixVectorMultiply(bX,theta.ini)-as.vector(LD%*%gamma.ini))
rss=sum(res*(Theta%*%res))
empirical.variance.formula=min(rss/(m-sum(theta.ini!=0)-sum(gamma.ini!=0)),empirical.variance.upper)
empirical.variance.formula=min(empirical.variance.formula,empirical.variance.upper)
gamma1=gamma
delta=gamma1
empirical.variance=empirical.variance.formula
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-as.vector(LD%*%gamma)
XtX=BtB
Xty=matrixVectorMultiply(Bt,res.theta)
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=F,max_iter=susie.iter)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
delta=delta+rho*(gamma-gamma1)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma1))
rss=sum(res*(Theta%*%res))
empirical.variance=rss/(m-sum(gamma1!=0)-min(sum(theta!=0),Lvec[vstar]))
empirical.variance=max(empirical.variance,empirical.variance.lower)
empirical.variance=min(empirical.variance,empirical.variance.upper)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
############################### inference #########################
cat("making inference\n")
gamma=gamma1
theta=theta
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma1!=0)
indvalid=which(gamma1==0)
res=by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma1)
if(sampling.time==0){
ThetaList=NULL
if(length(indtheta)==0){
indtheta=c(1:2)
}
pstar=length(indtheta)
bZ=cbind(bX[,indtheta],LD[,indgamma])
bZ=as.matrix(bZ)
bZinv=as.matrix(Theta%*%bZ)
H=matrixMultiply(t(bZ),as.matrix(Theta%*%bZ))
H[1:pstar,1:pstar]=H[1:pstar,1:pstar]-Rxysum[indtheta,indtheta]
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
theta.cov=diag(p)*0
theta.cov[indtheta,indtheta]=COV[1:pstar,1:pstar]
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
Thetaj <- solve(LDj)
Thetarhoj <- solve(LDj+rho*diag(length(indj)))
Btj <- as.matrix(t(bX[indj, ]) %*% Thetaj)
BtBj <- Btj%*%bX[indj, ]
BtBj=(t(BtBj)+BtBj)/2
thetaj=theta
gammaj=gamma1j=gamma
deltaj=gammaj*0
for(jiter in 1:sampling.iter){
indvalidj <- which(gamma1j==0)
indvalidj <- intersect(indvalidj, indj)
Rxysumj <- biasterm(RxyList = RxyList, indvalidj)
res.thetaj=by[indj]-as.vector(LD[indj,indj]%*%gammaj[indj])
XtXj=BtBj
Xtyj=matrixVectorMultiply(Btj,res.thetaj)
ytyj=sum(res.thetaj*(Thetaj%*%res.thetaj))
fit.thetaj=susie_suff_stat(XtX=BtBj,Xty=Xtyj,yty=ytyj,n=length(indvalidj),L=Lvec[vstar],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=F,max_iter=sampling.iter,s_init=fit.theta)
thetaj=coef.susie(fit.thetaj)[-1]*(fit.thetaj$pip>pip.thres)
indthetaj=which(thetaj!=0)
if(length(indthetaj)==1){
xtxj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
XtXj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
Xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=c(solve(XtXj)%*%Xtyj)
}
if((norm(thetaj, "2") / norm(theta.ini1, "2")) > maxdiff) {
thetaj <- thetaj / norm(thetaj, "2") * maxdiff * norm(theta.ini1, "2")
}
gammaj[indj]=as.vector(Thetarhoj%*%(by[indj]-matrixVectorMultiply(bX[indj, ],thetaj)-deltaj[indj]+rho*gamma1j[indj]))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
ThetaList[j, ] <- thetaj
}
close(pb)
theta.se=colSD(ThetaList)*sqrt((m-length(indtheta))/(m-length(indtheta)-length(indgamma)))
theta.cov=cov(ThetaList)*(m-length(indtheta))/(m-length(indtheta)-length(indgamma))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
}
A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$Bic.susie=Bbic.susie
A$Bic.ipod=Bbic.ipod
A$reliability.adjust=r
A$susie.theta=fit.theta
A$thetalist=ThetaList
}
}
return(A)
}
