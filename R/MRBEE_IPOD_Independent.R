MRBEE_IPOD_Independent=function(by,bX,byse,bXse,Rxy,Nmin=F,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,ebic.gamma=1,reliability.thres=0.5,rho=2,maxdiff=1.5,sampling.time=100,sampling.frac=0.5,sampling.iter=5,theta.ini=F,gamma.ini=F){
if(is.vector(bX)==T){
A=MRBEE_IPOD_UV_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,Nmin=Nmin,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,maxdiff=maxdiff,theta.ini=theta.ini,gamma.ini=gamma.ini,reliability.thres=reliability.thres)
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
BtB=matrixMultiply(t(bX),bX)
r=reliability.adj(bX,bXse,Theta="identity",thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak")
gamma.ini=gamma.ini1=fit0$gamma
theta.ini=theta.ini1=fit0$theta
}else{
gamma.ini=gamma.ini1=gamma.ini/byse1
theta.ini=theta.ini1=theta.ini
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
g=matrixVectorMultiply(t(bX),by-gamma)-Rxysum[1:p,p+1]
theta=c(matrixVectorMultiply(Hinv,g))
if((norm(theta,"2")/norm(theta.ini,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini,"2")
}
########################### update gamma ############################
gamma=(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1)/(1+rho)
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
res=c(by-matrixVectorMultiply(bX,theta)-gamma)
rss=sum(res^2)
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
res=c(by-matrixVectorMultiply(bX,theta)-gamma1)
if(sampling.time==0){
bZ=cbind(bX,diag(m)[,indgamma])
H=matrixMultiply(t(bZ),bZ)
H[1:p,1:p]=H[1:p,1:p]-Rxysum[1:p,1:p]
H=matrixInverse(H)
Hat=matrixListProduct(list(bZ,H,t(bZ)))
Hat=1-diag(Hat)
Hat[Hat<0.5]=0.5
E=-bZ*(res/Hat)
V=matrixMultiply(t(E),E)
COV=H%*%V%*%H
theta.cov=COV[1:p,1:p]
theta.se=sqrt(diag(theta.cov))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
}else{
ThetaList=matrix(0,sampling.time,p)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
indj=sample(m,sampling.frac*m,replace=F)
indj <- sort(indj)
BtB <- matrixMultiply(t(bX[indj, ]), bX[indj, ])
thetaj=theta
gammaj=gamma1j=gamma
deltaj=0*gammaj
for(jiter in 1:sampling.iter){
indvalidj <- which(gamma1j==0)
indvalidj <- intersect(indvalidj, indj)
Rxysumj <- biasterm(RxyList = RxyList, indvalidj)
Hinv <- solve(BtB - Rxysumj[1:p, 1:p])
g <- matrixVectorMultiply(t(bX[indj,]), by[indj] - gammaj[indj]) - Rxysumj[1:p, p + 1]
thetaj <- c(matrixVectorMultiply(Hinv, g))
if((norm(thetaj, "2") / norm(theta.ini, "2")) > maxdiff) {
thetaj <- thetaj / norm(thetaj, "2") * maxdiff * norm(theta.ini, "2")
}
gammaj[indj]=(by[indj]-matrixVectorMultiply(bX[indj, ],thetaj)-deltaj[indj]+rho*gamma1j[indj])/(1+rho)
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
A$gamma.ini=gamma.ini
A$reliability.adjust=r
}
return(A)
}
