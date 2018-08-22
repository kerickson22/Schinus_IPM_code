# Constructs a integral projection model using diameter and height as two continuous
# structuring state variables, where population transition probabilities are
# modeled on two size domains.

# This script is modified (by KDE and CCH) from the code provided by Ellner and Rees, 2006 that is available for download at: 
# https://www.jstor.org/stable/get_asset/10.1086/499438?supp_index=0&refreqid=excelsior%3Aa8885edea7fee9610a049fc3992a448a

#Here, we expand on Ellner and Rees' code by modeling the population dynamics on two continuous size domains 
# (See details in associated manuscript)

# Due to the large size of the objects created, this code can be slow to run. 
# To save disk space, whenever possible, large objects are saved and then removed from the workspace 

#clear everything
rm(list=ls(all=TRUE))

require(Matrix);
library(MASS)
library(mvtnorm)
library(msm)

#increase memory limit on Windows machines  (This will make everything run faster)
memory.limit(memory.limit() *2^30)

# Set matrix size: 
# m1: The number of categories of diameter that the seedling domain (D1) is divided into
# m2: The number of categories of height that the seedling domain (D1) is divided into
# m3: The number of categories of diameter that the large plant domain (D2) is divided into
# m4: The number of categories of height that the large plant domain (D2) is divided into

#The values used here were chosen by continuing to increase the dimensions of the matrix 
# until the calculated value of lambda stabilized 
m1=10
m2=m1+1
m3=100
m4=m3+1
tol=1.e-8; 




# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(200),...) {
	nx=length(x); ny=length(y); 
	x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
	y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}
#============================================================================================#
# (I) Parameters and demographic functions for computing the kernel 
#============================================================================================#

#p.vecs are defined externally in script 1demography_models.R
load('./Overall/p.vec_overall.RData')
load('./Eastern/p.vec_E.RData')
load('./Hybrid/p.vec_H.RData')
load('./Western/p.vec_W.RData')

#============================================================================# 
#  Define the kernels and iteration matrix:
#Domain 1: Seedling Domain
#	for diam=diameter in range [0, 1.6]
#	for height in range [0, 16]
#Domain 2: Larger Domain
#	 for diam=diameter in range [1.6,700]
#and for height=height in range [16, 800]
#============================================================================# 




#Survival-Growth of Seedlings (in D1)
pyx1=function(diamp,heightp,diam, height, params) { (1-((exp(params[12]+params[13]*diam + params[14]*height)/(1+exp(params[12]+params[13]*diam + params[14]*height))
)))*(exp(params[1]+params[2]*diam+params[3]*height)/(1+exp(params[1]+params[2]*diam+params[3]*height)))*dtnorm(diamp,mean=params[4] + params[5]*diam + params[6]*height,sd=params[7], lower=0, upper=1.6)*dtnorm(heightp,mean=params[8]+params[9]*diam + params[10]*height,sd=params[11], lower=0, upper=16)}


#Survival-Growth of Larger Plants (in D2)
#Note that D2 survival is multiplied by 0.997 to avoid 100% survival 

pyx2=function(diamp,heightp,diam, height, params) { 0.997*(exp(params[19]+params[20]*diam+params[21]*height)/(1+exp(params[19]+params[20]*diam+params[21]*height)))*dtnorm(diamp,mean=params[22] + params[23]*diam + params[24]*height,sd=params[25], lower=1.6, upper=700)*dtnorm(heightp,mean=params[26]+params[27]*diam + params[28]*height,sd=params[29], lower=16, upper=800)}


#Fecundity=P(fruiting)*# of Fruits Produced*P(survival of seeds)*Distribution of Seedling Diameters*Distribution of Seedling Heights
fyx=function(diamp,heightp,diam,height, params) {(exp(params[30]+params[31]*height)/(1+exp(params[30]+params[31]*height)))*(params[32]*diam*diam)*params[38]*params[37]*params[39]*dtnorm(diamp,mean=params[33],sd=params[34], lower=0, upper=1.6)*dtnorm(heightp,mean=params[35],sd=params[36], lower=0, upper=34)}

#Graduation = P(graduation)*Distribution of Graduates Diams*Distribution of Graduates Heights
gyx=function(diamp,heightp,diam,height, params) {((exp(params[1]+params[2]*diam+params[3]*height)/(1+exp(params[1]+params[2]*diam+params[3]*height))))*(exp(params[12]+params[13]*diam + params[14]*height)/(1+exp(params[12]+params[13]*diam + params[14]*height))
)*dtnorm(diamp,mean=params[15],sd=params[16], lower=1.6, upper=700)*dtnorm(heightp,mean=params[17],sd=params[18], lower=16, upper=800)}

#kyx=function(diamp,heightp,diam,height, params) {pyx(diamp,heightp, diam,height, params)+fyx(diamp,heightp,diam,height, params)}

# Compute meshpoints  

###Fix these to match domains 
h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
h3=(700-1.6)/m3;
y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
h4=(800-16)/m4
 y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2

# Compute the iteration matrix. With a bit of vectorizing it's not too slow,
# though you can probably do better if you need to. The shortcuts here have 
# been checked against the results from code that uses loops for everything. (comment from Ellner and Rees)

 
###Construct D1 (Seedling Domain):
plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in A 
Plop=outer(1:m1,1:m2,plop); 

D1_overall=matrix(0,m1*m2,m1*m2); 
D1_E=matrix(0,m1*m2,m1*m2); 
D1_H=matrix(0,m1*m2,m1*m2); 
D1_W=matrix(0,m1*m2,m1*m2); 

Kvals_D1_overall=array(0,c(m1,m2,m1,m2));  
Kvals_D1_E=array(0,c(m1,m2,m1,m2));  
Kvals_D1_H=array(0,c(m1,m2,m1,m2));  
Kvals_D1_W=array(0,c(m1,m2,m1,m2));  

#Overall
for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_overall)
				D1_overall[Plop[k,1:m2],Plop[i,j]]=kvals
				Kvals_D1_overall[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}
D1_overall=D1_overall*h1*h2 #multiply D1 by widths
save(D1_overall, file="./Overall/D1_overall.RData")
save(Kvals_D1_overall, file="./Overall/Kvals_D1_overall.RData")
rm(D1_overall, Kvals_D1_overall)
# Eastern
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_E)
      D1_E[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_E[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_E=D1_E*h1*h2 #multiply D1 by widths
save(D1_E, file="./Eastern/D1_E.RData")
save(Kvals_D1_E, file="./Eastern/Kvals_D1_E.RData")
rm(D1_E, Kvals_D1_E)
# Hybrid
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_H)
      D1_H[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_H[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_H=D1_H*h1*h2 #multiply D1 by widths
save(D1_H, file="./Hybrid/D1_H.RData")
save(Kvals_D1_H, file="./Hybrid/Kvals_D1_H.RData")
rm(D1_H, Kvals_D1_H)
# Western
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_W)
      D1_W[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_W[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_W=D1_W*h1*h2 #multiply D1 by widths
save(D1_W, file="./Western/D1_W.RData")
save(Kvals_D1_W, file="./Western/Kvals_D1_W.RData")
rm(D1_W, Kvals_D1_W)



		
###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 

D2_overall=matrix(0,m3*m4,m3*m4);
D2_E=matrix(0,m3*m4,m3*m4)
D2_H=matrix(0,m3*m4,m3*m4)
D2_W=matrix(0,m3*m4,m3*m4)

Kvals_D2_overall=array(0,c(m3,m4,m3,m4));  
Kvals_D2_E=array(0,c(m3,m4,m3,m4)); 
Kvals_D2_H=array(0,c(m3,m4,m3,m4)); 
Kvals_D2_W=array(0,c(m3,m4,m3,m4)); 

#overall
for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_overall)
				D2_overall[Plop[k,1:m4],Plop[i,j]]=kvals
				Kvals_D2_overall[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		
D2_overall=D2_overall*h3*h4 #Multiply D2 by widths
save(D2_overall, file="./Overall/D2_overall.RData")
save(Kvals_D2_overall, file="./Overall/Kvals_D2_overall")
rm(D2_overall, Kvals_D2_overall)

#Eastern
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_E)
      D2_E[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_E[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_E=D2_E*h3*h4 #Multiply D2 by widths
save(D2_E, file="./Eastern/D2_E.RData")
save(Kvals_D2_E, file="./Eastern/Kvals_D2_E.RData")
rm(D2_E, Kvals_D2_E)
#Hybrid
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_H)
      D2_H[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_H[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_H=D2_H*h3*h4 #Multiply D2 by widths
save(D2_H, file="./Hybrid/D2_H.RData")
save(Kvals_D2_H, file="./Hybrid/Kvals_D2_H.RData")
rm(D2_H, Kvals_D2_H)
#Western
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_W)
      D2_W[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_W[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_W=D2_W*h3*h4 #Multiply D2 by widths
save(D2_W, file="./Western/D2_W.RData")
save(Kvals_D2_W, file="./Western/Kvals_D2_W.RData")
rm(D2_W, Kvals_D2_W)




###Construct F (Fertility):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);
F_overall=matrix(0,m1*m2,m3*m4); 
F_E=matrix(0,m1*m2,m3*m4);
F_H=matrix(0,m1*m2,m3*m4);
F_W=matrix(0,m1*m2,m3*m4);

Kvals_F_overall=array(0, c(m1, m2, m3, m4))
Kvals_F_E=array(0, c(m1, m2, m3, m4))
Kvals_F_H=array(0, c(m1, m2, m3, m4))
Kvals_F_W=array(0, c(m1, m2, m3, m4))

#Overall
for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_overall)
			F_overall[Plop1[k, 1:m2], Plop2[i,j]]=kvals
			Kvals_F_overall[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}
F_overall=F_overall*h1*h2
save(F_overall, file="./Overall/F_overall.RData")
save(Kvals_F_overall, file="./Overall/Kvals_F_overall.RData")
rm(F_overall, Kvals_F_overall)

#Eastern
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_E)
      F_E[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_E[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_E=F_E*h1*h2
save(F_E, file="./Eastern/F_E.RData")
save(Kvals_F_E, file="./Eastern/Kvals_F_E.RData")
rm(F_E, Kvals_F_E)

#Hybrid
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_H)
      F_H[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_H[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_H=F_H*h1*h2
save(F_H, file="./Hybrid/F_H.RData")
save(Kvals_F_H, file="./Hybrid/Kvals_F_H.RData")
rm(F_H, Kvals_F_H)

#Western
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_W)
      F_W[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_W[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_W=F_W*h1*h2
save(F_W, file="./Western/F_W.RData")
save(Kvals_F_W, file="./Western/Kvals_F_W.RData")
rm(F_W, Kvals_F_W)



###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);
G_overall=matrix(0,m3*m4,m1*m2); 
G_E=matrix(0,m3*m4,m1*m2); 
G_H=matrix(0,m3*m4,m1*m2); 
G_W=matrix(0,m3*m4,m1*m2); 

Kvals_G_overall=array(0, c(m3, m4, m1, m2))
Kvals_G_E=array(0, c(m3, m4, m1, m2))
Kvals_G_H=array(0, c(m3, m4, m1, m2))
Kvals_G_W=array(0, c(m3, m4, m1, m2))

#Overall
for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_overall)
			G_overall[Plop1[k, 1:m4], Plop2[i,j]]=kvals
			Kvals_G_overall[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}
G_overall=G_overall*h3*h4
save(G_overall, file="./Overall/G_overall.RData")
save(Kvals_G_overall, file="./Overall/Kvals_G_overall.RData")
rm(G_overall, Kvals_G_overall)

#Eastern
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_E)
      G_E[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_E[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_E=G_E*h3*h4
save(G_E, file="./Eastern/G_E.RData")
save(Kvals_G_E, file="./Eastern/Kvals_G_E.RData")
rm(G_E, Kvals_G_E)

#Hybrid
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_H)
      G_H[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_H[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_H=G_H*h3*h4
save(G_H, file="./Hybrid/G_H.RData")
save(Kvals_G_H, file="./Hybrid/Kvals_G_H.RData")
rm(G_H, Kvals_G_H)

#Western
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_W)
      G_W[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_W[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_W=G_W*h3*h4
save(G_W, file="./Western/G_W.RData")
save(Kvals_G_W, file="./Western/Kvals_G_W.RData")
rm(G_W, Kvals_G_W)




#Putting the matrices together:

#Overall
rm(Plop, Plop1, Plop2, i, j, k, kvals ) #clear out workspace to maximize space available 
load("./Overall/D1_overall.RData")
load("./Overall/G_overall.RData")
left_side<-rbind(D1_overall, G_overall)
rm(D1_overall, G_overall)
load("./Overall/F_overall.RData")
load("./Overall/D2_overall.RData")
right_side<-rbind(F_overall, D2_overall)
rm(F_overall, D2_overall)
A_overall<-cbind(left_side, right_side)
save(A_overall, file="./Overall/A_overall.RData")

#Eastern
rm(A_overall, left_side, right_side) #clear out workspace to maximize space available 
load("./Eastern/D1_E.RData")
load("./Eastern/G_E.RData")
left_side<-rbind(D1_E, G_E)
rm(D1_E, G_E)
load("./Eastern/F_E.RData")
load("./Eastern/D2_E.RData")
right_side<-rbind(F_E, D2_E)
rm(F_E, D2_E)
A_E<-cbind(left_side, right_side)
save(A_E, file="./Eastern/A_E.RData")

#Hybrid
rm(A_E, left_side, right_side) #clear out workspace to maximize space available 
load("./Hybrid/D1_H.RData")
load("./Hybrid/G_H.RData")
left_side<-rbind(D1_H, G_H)
rm(D1_H, G_H)
load("./Hybrid/F_H.RData")
load("./Hybrid/D2_H.RData")
right_side<-rbind(F_H, D2_H)
rm(F_H, D2_H)
A_H<-cbind(left_side, right_side)
save(A_H, file="./Hybrid/A_H.RData")

#Western
rm(A_H, left_side, right_side) #clear out workspace to maximize space available 
load("./Western/D1_W.RData")
load("./Western/G_W.RData")
left_side<-rbind(D1_W, G_W)
rm(D1_W, G_W)
load("./Western/F_W.RData")
load("./Western/D2_W.RData")
right_side<-rbind(F_W, D2_W)
rm(F_W, D2_W)
A_W<-cbind(left_side, right_side)
save(A_W, file="./Western/A_W.RData")


rm(A_W, left_side, right_side)




#============================================================================# 
#  Find lambda, w by iteration 
#
#  Note: the Matrix package is used to iterate more quickly. The check  
#  for convergence requires extracting matrix entries via the @x slot
#  of a Matrix object. Matrix is S4-style -- see ?Matrix. 
#============================================================================# 

#Overall
load("./Overall/A_overall.RData")
A2=Matrix(A_overall); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	nt1=A2%*%nt;
	qmax=sum(abs((nt1-lam*nt)@x));  
	lam=sum(nt1@x); 
	nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
	cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
stable.dist=nt
lam.stable=lam;

# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_overall%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_overall<-lam.stable
stable.dist_overall<-stable.dist
save(lam.stable_overall, file="./Overall/lam.stable_overall.RData")
save(stable.dist_overall, file="./Overall/stable.dist_overall.RData")


##Alternate calculation of rv:

A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	nt1=A3%*%nt;
	qmax=sum(abs((nt1-lam*nt)@x));  
	lam=sum(nt1@x); 
	nt@x=(nt1@x)/lam;   
	cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
#stable.dist=nt
lam.stable.t=lam; 


# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_overall%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






#============================================================================# 
#  Find reproductive value function by iteration.  
#  Note that we use left-multiplication to do the transpose iteration,
#  which is equivalent to using the transpose kernel. 
#============================================================================# 
vt=Matrix(1,1,m1*m2+m3*m4); vt1=vt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	vt1=vt%*%A2;
	qmax=sum(abs((vt1-lam*vt)@x));  
	lam=sum(vt1@x); 
	vt@x=(vt1@x)/lam;   
	cat(lam,qmax,"\n");
} 
v=t(matrix(vt@x,1,m1*m2+m3*m4)); 
lam.stable.t=lam; 

v_overall<-v
save(v_overall, file="./Overall/v_overall.RData")




#Eastern
load("./Eastern/A_E.RData")
A2=Matrix(A_E); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  nt1=A2%*%nt;
  qmax=sum(abs((nt1-lam*nt)@x));  
  lam=sum(nt1@x); 
  nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
  cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
stable.dist=nt
lam.stable=lam;

# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_E%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_E<-lam.stable
stable.dist_E<-stable.dist
save(lam.stable_E, file="./Eastern/lam.stable_E.RData")
save(stable.dist_E, file="./Eastern/stable.dist_E.RData")


##Alternate calculation of rv:

A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  nt1=A3%*%nt;
  qmax=sum(abs((nt1-lam*nt)@x));  
  lam=sum(nt1@x); 
  nt@x=(nt1@x)/lam;   
  cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
#stable.dist=nt
lam.stable.t=lam; 


# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_E%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






#============================================================================# 
#  Find reproductive value function by iteration.  
#  Note that we use left-multiplication to do the transpose iteration,
#  which is equivalent to using the transpose kernel. 
#============================================================================# 
vt=Matrix(1,1,m1*m2+m3*m4); vt1=vt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  vt1=vt%*%A2;
  qmax=sum(abs((vt1-lam*vt)@x));  
  lam=sum(vt1@x); 
  vt@x=(vt1@x)/lam;   
  cat(lam,qmax,"\n");
} 
v=t(matrix(vt@x,1,m1*m2+m3*m4)); 
lam.stable.t=lam; 

v_E<-v
save(v_E, file="./Eastern/v_E.RData")


#Hybrid
load("./Hybrid/A_H.RData")
A2=Matrix(A_H); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  nt1=A2%*%nt;
  qmax=sum(abs((nt1-lam*nt)@x));  
  lam=sum(nt1@x); 
  nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
  cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
stable.dist=nt
lam.stable=lam;

# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_H%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_H<-lam.stable
stable.dist_H<-stable.dist
save(lam.stable_H, file="./Hybrid/lam.stable_H.RData")
save(stable.dist_H, file="./Hybrid/stable.dist_H.RData")


##Alternate calculation of rv:

A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  nt1=A3%*%nt;
  qmax=sum(abs((nt1-lam*nt)@x));  
  lam=sum(nt1@x); 
  nt@x=(nt1@x)/lam;   
  cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
#stable.dist=nt
lam.stable.t=lam; 


# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_H%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






#============================================================================# 
#  Find reproductive value function by iteration.  
#  Note that we use left-multiplication to do the transpose iteration,
#  which is equivalent to using the transpose kernel. 
#============================================================================# 
vt=Matrix(1,1,m1*m2+m3*m4); vt1=vt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  vt1=vt%*%A2;
  qmax=sum(abs((vt1-lam*vt)@x));  
  lam=sum(vt1@x); 
  vt@x=(vt1@x)/lam;   
  cat(lam,qmax,"\n");
} 
v=t(matrix(vt@x,1,m1*m2+m3*m4)); 
lam.stable.t=lam; 

v_H<-v
save(v_H, file="./Hybrid/v_H.RData")


#Western
load("./Western/A_W.RData")
A2=Matrix(A_overall); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  nt1=A2%*%nt;
  qmax=sum(abs((nt1-lam*nt)@x));  
  lam=sum(nt1@x); 
  nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
  cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
stable.dist=nt
lam.stable=lam;

# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_W%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_W<-lam.stable
stable.dist_W<-stable.dist
save(lam.stable_W, file="./Western/lam.stable_W.RData")
save(stable.dist_W, file="./Western/stable.dist_W.RData")


##Alternate calculation of rv:

A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  nt1=A3%*%nt;
  qmax=sum(abs((nt1-lam*nt)@x));  
  lam=sum(nt1@x); 
  nt@x=(nt1@x)/lam;   
  cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
#stable.dist=nt
lam.stable.t=lam; 


# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A_W%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






#============================================================================# 
#  Find reproductive value function by iteration.  
#  Note that we use left-multiplication to do the transpose iteration,
#  which is equivalent to using the transpose kernel. 
#============================================================================# 
vt=Matrix(1,1,m1*m2+m3*m4); vt1=vt; 

qmax=1000; lam=1; 
while(qmax>tol) {
  vt1=vt%*%A2;
  qmax=sum(abs((vt1-lam*vt)@x));  
  lam=sum(vt1@x); 
  vt@x=(vt1@x)/lam;   
  cat(lam,qmax,"\n");
} 
v=t(matrix(vt@x,1,m1*m2+m3*m4)); 
lam.stable.t=lam; 

v_W<-v
save(v_W, file="./Western/v_W.RData")



#####################################
####
#### ELASTICITY AND SENSITIVITY
####
#####################################




# Compute elasticity matrix 
#Ellner: 
#repro.val=matrix(v, m1, m2);
#stable.state=matrix(stable.dist, m1, m2);
#v.dot.w=sum(stable.state_smalls*repro.val_smalls)
#sens=outer(repro.val_smalls,stable.state_smalls)/v.dot.w
#elas=sens*Kvals/lam.stable;
# # Compute matrix of total(=integrated) elasticities for all transitions (x_i,q_j) -> anywhere 
#total.elas=h1*h2*apply(elas,c(3,4),sum);
# Checks
#cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
#cat("Integrated elasticity=",sum(h1*h2*h1*h2*elas_D1)," Should =1","\n"); 


####Method: Compute elasticity matrix first then separate:

#Overall
load("v_overall.RData")
load("stable.dist_overall.RData")
load("A_overall.RData")
load("lam.stable_overall.RData")

v.dot.w_overall<-sum(t(v_overall)%*%stable.dist_overall)
norm_v_overall<-v_overall/v.dot.w_overall
check_overall<-t(norm_v_overall)%*%stable.dist_overall #should be 1

rv_overall<-norm_v_overall

sens_overall<-norm_v_overall%*%t(stable.dist_overall)
elas_overall<-sens_overall*A_overall/lam.stable_overall

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_overall<-elas_overall[1:(m1*m2), 1:(m1*m2)]
elas_G_overall<-elas_overall[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_overall<-elas_overall[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_overall<-elas_overall[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_overall<-sens_overall[1:(m1*m2), 1:(m1*m2)]
sens_G_overall<-sens_overall[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_overall<-sens_overall[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_overall<-sens_overall[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_overall, file="elas_D1_overall.RData")
save(elas_G_overall, file="elas_G_overall.RData")
save(elas_F_overall, file="elas_F_overall.RData")
save(elas_D2_overall, file="elas_D2_overall.RData")
save(sens_D1_overall, file="sens_D1_overall.RData")
save(sens_G_overall, file="sens_G_overall.RData")
save(sens_F_overall, file="sens_F_overall.RData")
save(sens_D2_overall, file="sens_D2_overall.RData")
rm(list=ls())

#Eastern
load("v_E.RData")
load("stable.dist_E.RData")
load("A_E.RData")
load("lam.stable_E.RData")

v.dot.w_E<-sum(t(v_E)%*%stable.dist_E)
norm_v_E<-v_E/v.dot.w_E
check_E<-t(norm_v_E)%*%stable.dist_E #should be 1

rv_E<-norm_v_E

sens_E<-norm_v_E%*%t(stable.dist_E)
elas_E<-sens_E*A_E/lam.stable_E

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_E<-elas_E[1:(m1*m2), 1:(m1*m2)]
elas_G_E<-elas_E[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_E<-elas_E[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_E<-elas_E[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_E<-sens_E[1:(m1*m2), 1:(m1*m2)]
sens_G_E<-sens_E[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_E<-sens_E[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_E<-sens_E[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_E, file="elas_D1_E.RData")
save(elas_G_E, file="elas_G_E.RData")
save(elas_F_E, file="elas_F_E.RData")
save(elas_D2_E, file="elas_D2_E.RData")
save(sens_D1_E, file="sens_D1_E.RData")
save(sens_G_E, file="sens_G_E.RData")
save(sens_F_E, file="sens_F_E.RData")
save(sens_D2_E, file="sens_D2_E.RData")
rm(list=ls())

#Overall
load("v_H.RData")
load("stable.dist_H.RData")
load("A_H.RData")
load("lam.stable_H.RData")

v.dot.w_H<-sum(t(v_H)%*%stable.dist_H)
norm_v_H<-v_H/v.dot.w_H
check_H<-t(norm_v_H)%*%stable.dist_H #should be 1

rv_H<-norm_v_H

sens_H<-norm_v_H%*%t(stable.dist_H)
elas_H<-sens_H*A_H/lam.stable_H

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_H<-elas_H[1:(m1*m2), 1:(m1*m2)]
elas_G_H<-elas_H[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_H<-elas_H[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_H<-elas_H[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_H<-sens_H[1:(m1*m2), 1:(m1*m2)]
sens_G_H<-sens_H[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_H<-sens_H[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_H<-sens_H[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_H, file="elas_D1_H.RData")
save(elas_G_H, file="elas_G_H.RData")
save(elas_F_H, file="elas_F_H.RData")
save(elas_D2_H, file="elas_D2_H.RData")
save(sens_D1_H, file="sens_D1_H.RData")
save(sens_G_H, file="sens_G_H.RData")
save(sens_F_H, file="sens_F_H.RData")
save(sens_D2_H, file="sens_D2_H.RData")
rm(list=ls())

#Western
load("v_W.RData")
load("stable.dist_W.RData")
load("A_W.RData")
load("lam.stable_W.RData")

v.dot.w_W<-sum(t(v_W)%*%stable.dist_W)
norm_v_W<-v_W/v.dot.w_W
check_W<-t(norm_v_W)%*%stable.dist_W #should be 1

rv_W<-norm_v_W

sens_W<-norm_v_W%*%t(stable.dist_W)
elas_W<-sens_W*A_W/lam.stable_W

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_W<-elas_W[1:(m1*m2), 1:(m1*m2)]
elas_G_W<-elas_W[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_W<-elas_W[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_W<-elas_W[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_W<-sens_W[1:(m1*m2), 1:(m1*m2)]
sens_G_W<-sens_W[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_W<-sens_W[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_W<-sens_w[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_W, file="elas_D1_W.RData")
save(elas_G_W, file="elas_G_W.RData")
save(elas_F_W, file="elas_F_W.RData")
save(elas_D2_W, file="elas_D2_W.RData")
save(sens_D1_W, file="sens_D1_W.RData")
save(sens_G_W, file="sens_G_W.RData")
save(sens_F_W, file="sens_F_W.RData")
save(sens_D2_W, file="sens_D2_W.RData")
rm(list=ls())

########
###Now go back to unpack the elasticity matrices into components:
########

load("elas_D1_overall.RData")
load("elas_D1_E.RData")
load("elas_D1_H.RData")
load("elas_D1_W.RData")

plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in array
Plop=outer(1:m1,1:m2,plop); 

Kvals_elas_D1=array(0,c(m1,m2,m1,m2));  

for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals= elas_D1[Plop[k,1:m2],Plop[i,j]]
				Kvals_elas_D1[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_elas_D2=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= elas_D2[Plop[k,1:m4],Plop[i,j]]
				Kvals_elas_D2[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_elas_F=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=elas_F[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_elas_F[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_elas_G=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=elas_G[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_elas_G[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


##############
total.elas_D1<-apply(Kvals_elas_D1, c(3,4), sum);
total.elas_D2<-apply(Kvals_elas_D2, c(3,4), sum);
total.elas_F<-apply(Kvals_elas_F, c(3,4), sum);
total.elas_G<-apply(Kvals_elas_G, c(3,4), sum);

#Separate:
repro.val_smalls=matrix(v[1:(m1*m2)], m1,m2);
stable.state_smalls=matrix(stable.dist[1:(m1*m2)], m1, m2);

repro.val_larges=matrix(v[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges=matrix(stable.dist[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);


#============================================================================# 
#  Calculate stable state distribution and state-dependent total elasticity 
#============================================================================# 

ss_diam1<-apply(stable.state_smalls, 1, sum)
ss_diam2<-apply(stable.state_larges, 1, sum)
ss_diam<-c(ss_diam1, ss_diam2)
y_diam<-c(y1, y3)



ss_height1<-apply(stable.state_smalls, 2, sum)
ss_height2<-apply(stable.state_larges, 2, sum)
ss_height<-c(ss_height1, ss_height2)
y_height<-c(y2, y4)




rv_diam1<-apply(repro.val_smalls, 1, sum)
rv_diam2<-apply(repro.val_larges, 1, sum)
rv_diam<-c(rv_diam1, rv_diam2)
y_diam<-c(y1, y3)





rv_height1<-apply(repro.val_smalls, 2, sum)
rv_height2<-apply(repro.val_larges, 2, sum)
rv_height<-c(rv_height1, rv_height2)
y_height<-c(y2, y4)



total.elas_D1_overall<-total.elas_D1
total.elas_D2_overall<-total.elas_D2
total.elas_G_overall<-total.elas_G
total.elas_F_overall<-total.elas_F

save(total.elas_D1_all, file="total.elas_D1_all.RData")
save(total.elas_D2_all, file="total.elas_D2_all.RData")
save(total.elas_G_all, file="total.elas_G_all.RData")
save(total.elas_F_all, file="total.elas_F_all.RData")



setwd("/Users/ke2/Dropbox/matrix/outputs/sites/")


##Redo these with comparisons between overall and Eastern 
elas_D1_ref_E<-total.elas_D1_all-total.elas_D1_E
elas_D2_ref_E<-total.elas_D2_all-total.elas_D2_E
elas_G_ref_E<-total.elas_G_all-total.elas_G_E
elas_F_ref_E<-total.elas_F_all-total.elas_F_E

elas_D1_ref_H<-total.elas_D1_all-total.elas_D1_H
elas_D2_ref_H<-total.elas_D2_all-total.elas_D2_H
elas_G_ref_H<-total.elas_G_all-total.elas_G_H
elas_F_ref_H<-total.elas_F_all-total.elas_F_H

elas_D1_ref_W<-total.elas_D1_all-total.elas_D1_W
elas_D2_ref_W<-total.elas_D2_all-total.elas_D2_W
elas_G_ref_W<-total.elas_G_all-total.elas_G_W
elas_F_ref_W<-total.elas_F_all-total.elas_F_W


jpeg('marginal_elasticity_height_ref.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(colSums(elas_D1_ref_W), xlab="Height (cm)", ylim=c(-.01, .01), ylab="Change in Elasticity", main="D1",  type='l', col=col_W_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_D1_ref_H), col=col_H_dark, type='l', lwd=2)
lines(colSums(elas_D1_ref_E), col=col_E_dark, type='l', lwd=2)


plot(colSums(elas_F_ref_H), ylim=c(-.01, 0.01),  xlab="Height (cm)", ylab= "Change in Elasticity", main="F", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_F_ref_E), type='l', col=col_E_dark, lwd=2)
lines(colSums(elas_F_ref_W), type='l', col=col_W_dark, lwd=2)

plot(colSums(elas_G_ref_H), ylim=c(-0.01, 0.01), xlab="Height (cm)", ylab="Change in Elasticity", main="G", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_G_ref_E), type='l', col=col_E_dark, lwd=2)
lines(colSums(elas_G_ref_W), type='l', col=col_W_dark, lwd=2)

plot(colSums(elas_D2_ref_H), ylim=c(-0.01, 0.01), xlab="Height (cm)", ylab="Change in Elasticity", main="D2", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_D2_ref_E), type='l', col=col_E_dark, lwd=2)
lines(colSums(elas_D2_ref_W), type='l', col=col_W_dark, lwd=2)

dev.off()

###
jpeg('marginal_elasticity_diam_ref.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(rowSums(elas_D1_ref_W), xlab="Diameter (mm) ", ylim=c(-.01, .01), ylab="Change in Elasticity", main="D1",  type='l', col=col_W_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_D1_ref_H), col=col_H_dark, type='l', lwd=2)
lines(rowSums(elas_D1_ref_E), col=col_E_dark, type='l', lwd=2)


plot(rowSums(elas_F_ref_H), ylim=c(-.01, 0.01),  xlab="Diameter (mm)", ylab= "Change in Elasticity", main="F", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_F_ref_E), type='l', col=col_E_dark, lwd=2)
lines(rowSums(elas_F_ref_W), type='l', col=col_W_dark, lwd=2)

plot(rowSums(elas_G_ref_H), ylim=c(-0.01, 0.01), xlab="Diameter (mm)", ylab="Change in Elasticity", main="G", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_G_ref_E), type='l', col=col_E_dark, lwd=2)
lines(rowSums(elas_G_ref_W), type='l', col=col_W_dark, lwd=2)

plot(rowSums(elas_D2_ref_H), ylim=c(-0.01, 0.01), xlab="Diameter (mm)", ylab="Change in Elasticity", main="D2", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_D2_ref_E), type='l', col=col_E_dark, lwd=2)
lines(rowSums(elas_D2_ref_W), type='l', col=col_W_dark, lwd=2)

dev.off()



###


jpeg('marginal_elasticity_diam_all.jpeg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=36)
plot(rowSums(total.elas_D1_all), ylim=c(0, 0.07),xlab="Size", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
plot(rowSums(total.elas_F_all),  ylim=c(0, 0.07), xlab="Size", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
plot(rowSums(total.elas_G_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
plot(rowSums(total.elas_D2_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
dev.off()

jpeg('marginal_elasticity_diam_height_all.jpeg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=36)
plot(rowSums(total.elas_D1_all), ylim=c(0, 0.07),xlab="Size", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_D1_all), type='l', lwd=2, col="red")
plot(rowSums(total.elas_F_all),  ylim=c(0, 0.07), xlab="Size", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_F_all), type='l', lwd=2, col="red")
plot(rowSums(total.elas_G_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_G_all), type='l', lwd=2, col="red")
plot(rowSums(total.elas_D2_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_D2_all), type='l', lwd=2, col="red")
dev.off()







###

jpeg('marginal_elasticity_diam_lineage.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(rowSums(total.elas_D1_E), ylim=c(0, 0.07),xlab="Diameter (mm)", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_D1_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_D1_H), type='l', lwd=2, col=col_H_dark)

plot(rowSums(total.elas_F_E),  ylim=c(0, 0.07), xlab="Diameter (mm)", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_F_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_F_H), type='l', lwd=2, col=col_H_dark)

plot(rowSums(total.elas_G_E), ylim=c(0, 0.07),  xlab="Diameter (mm)", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_G_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_G_H), type='l', lwd=2, col=col_H_dark)

plot(rowSums(total.elas_D2_E), ylim=c(0, 0.07),  xlab="Diameter (mm)", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_D2_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_D2_H), type='l', lwd=2, col=col_H_dark)
dev.off()

###

jpeg('marginal_elasticity_height_lineage.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(colSums(total.elas_D1_E), ylim=c(0, 0.07),xlab="Height (cm)", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_D1_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_D1_H), type='l', lwd=2, col=col_H_dark)

plot(colSums(total.elas_F_E),  ylim=c(0, 0.07), xlab="Height (cm)", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_F_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_F_H), type='l', lwd=2, col=col_H_dark)

plot(colSums(total.elas_G_E), ylim=c(0, 0.07),  xlab="Height (cm)", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_G_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_G_H), type='l', lwd=2, col=col_H_dark)

plot(colSums(total.elas_D2_E), ylim=c(0, 0.07),  xlab="Height (cm)", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_D2_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_D2_H), type='l', lwd=2, col=col_H_dark)
dev.off()

