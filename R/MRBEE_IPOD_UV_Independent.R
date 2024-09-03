MRBEE_IPOD_UV_Independent=function(by,bX,byse,bXse,Rxy,Nmin=F,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,ebic.gamma=1,maxdiff=1.5,sampling.time=100,sampling.frac=0.5,theta.ini=F,gamma.ini=F,reliability.thres=reliability.thres){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=length(bX)
if(Nmin==F){
Nmin=m
}
BtB=sum(bX^2)
r=reliability.adj.uv(bX,bXse,Theta="identity",thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
cluster.index=c(1:m)
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak")
gamma.ini=gamma.ini1=fit0$gamma
theta.ini=theta.ini1=fit0$theta
}
############################## Tuning Parameter ######################
w=length(tauvec)
Btheta=c(1:w)
Bgamma=matrix(0,m,w)
Bbic=tauvec
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
Hinv=1/(BtB-sum(bXse[indvalid]^2)*Rxy[1,1])
g=sum(bX*(by-gamma))-sum(bXse[indvalid])*Rxy[2,1]
theta=g*Hinv
########################### update gamma ############################
gamma=(by-bX*theta-delta+rho*gamma1)/(1+rho)
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)

iter=iter+1
if(iter>3){
error=abs(theta-theta1)
}
}
Btheta[j]=theta
Bgamma[,j]=gamma1
df1=sum(gamma1!=0)
res=by-bX*theta-gamma
rss=sum(res^2)
Bbic[j]=log(rss)*Nmin+(log(Nmin)+ebic.gamma*log(m))*df1
}
######################## Inference #################################
jstar=which.min(Bbic)
theta=Btheta[jstar]
gamma=Bgamma[,jstar]
error=1
iter=1
theta=theta
gamma=gamma
names(gamma)=rownames(bX)
indgamma=which(gamma!=0)
indvalid=which(gamma==0)
res=by-bX*theta-matrixVectorMultiply(LD,gamma)

ThetaList=c(1:sampling.time)
for(j in 1:sampling.time){
indj=sample(m,sampling.frac*m,replace=F)
indj=sort(indj)
gammaj=gamma1j=gamma
deltaj=delta*0
for(iterj in 1:5){
indvalidj=which(gamma1j==0)
indvalidj=intersect(indvalidj,indj)
BtB=sum(bX[indj]^2)
Hinv=1/(BtB-sum(bXse[indvalidj]^2)*Rxy[1,1])
g=sum(Bt*(by[indj]-gammaj[indj]))-sum(bXse[indvalidj])*Rxy[2,1]
thetaj=g*Hinv
resgammi=(by-bX*thetaj-deltaj+rho*gamma1j)/(1+rho)
gammaj[indj]=resgammai[indj]
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
ThetaList[j]=thetaj
}
theta.se=sd(ThetaList)*sqrt((m-length(theta))/(m-length(theta)-length(indgamma)))

A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini
A$reliability.adjust=r
return(A)
}
