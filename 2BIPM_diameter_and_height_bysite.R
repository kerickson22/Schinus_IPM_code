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
load('./BC/p.vec_BC.RData')
load('./CC/p.vec_CC.RData')
load('./C/p.vec_C.RData')
load('./FP/p.vec_FP.RData')
load('./PG/p.vec_PG.RData')
load('./WT/p.vec_WT.RData')

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

D1_BC=matrix(0,m1*m2,m1*m2); 
D1_CC=matrix(0,m1*m2,m1*m2);
D1_C=matrix(0,m1*m2,m1*m2); 
D1_FP=matrix(0,m1*m2,m1*m2); 
D1_PG=matrix(0,m1*m2,m1*m2); 
D1_WT=matrix(0,m1*m2,m1*m2); 
 

Kvals_D1_BC=array(0,c(m1,m2,m1,m2));  
Kvals_D1_CC=array(0,c(m1,m2,m1,m2));  
Kvals_D1_C=array(0,c(m1,m2,m1,m2));  
Kvals_D1_FP=array(0,c(m1,m2,m1,m2));  
Kvals_D1_PG=array(0,c(m1,m2,m1,m2));
Kvals_D1_WT=array(0,c(m1,m2,m1,m2));  

#Big Cypress
for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_BC)
				D1_BC[Plop[k,1:m2],Plop[i,j]]=kvals
				Kvals_D1_BC[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}
D1_BC=D1_BC*h1*h2 #multiply D1 by widths
save(D1_BC, file="./BC/D1_BC.RData")
save(Kvals_D1_BC, file="./BC/Kvals_D1_BC.RData")
rm(D1_BC, Kvals_D1_BC)
#Cape Canaveral
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_CC)
      D1_CC[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_CC[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_CC=D1_CC*h1*h2 #multiply D1 by widths
save(D1_CC, file="./CC/D1_CC.RData")
save(Kvals_D1_CC, file="./CC/Kvals_D1_CC.RData")
rm(D1_CC, Kvals_D1_CC)
#Chekika
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_C)
      D1_C[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_C[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_C=D1_C*h1*h2 #multiply D1 by widths
save(D1_C, file="./C/D1_C.RData")
save(Kvals_D1_C, file="./C/Kvals_D1_C.RData")
rm(D1_C, Kvals_D1_C)
#Fort Pierce
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_FP)
      D1_FP[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_FP[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_FP=D1_FP*h1*h2 #multiply D1 by widths
save(D1_FP, file="./FP/D1_FP.RData")
save(Kvals_D1_FP, file="./FP/Kvals_D1_FP.RData")
rm(D1_FP, Kvals_D1_FP)
#Punta Gorda
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_PG)
      D1_PG[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_PG[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_PG=D1_PG*h1*h2 #multiply D1 by widths
save(D1_PG, file="./PG/D1_PG.RData")
save(Kvals_D1_PG, file="./PG/Kvals_D1_PG.RData")
rm(D1_PG, Kvals_D1_PG)
#Wild Turkey
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_WT)
      D1_WT[Plop[k,1:m2],Plop[i,j]]=kvals
      Kvals_D1_WT[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
D1_WT=D1_WT*h1*h2 #multiply D1 by widths
save(D1_WT, file="./WT/D1_WT.RData")
save(Kvals_D1_WT, file="./WT/Kvals_D1_WT.RData")
rm(D1_WT, Kvals_D1_WT)



		
###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 

D2_BC=matrix(0,m3*m4,m3*m4);
D2_CC=matrix(0,m3*m4,m3*m4)
D2_C=matrix(0,m3*m4,m3*m4)
D2_FP=matrix(0,m3*m4,m3*m4)
D2_PG=matrix(0,m3*m4,m3*m4)
D2_WT=matrix(0,m3*m4,m3*m4)

Kvals_D2_BC=array(0,c(m3,m4,m3,m4));  
Kvals_D2_CC=array(0,c(m3,m4,m3,m4)); 
Kvals_D2_C=array(0,c(m3,m4,m3,m4)); 
Kvals_D2_FP=array(0,c(m3,m4,m3,m4));
Kvals_D2_PG=array(0,c(m3,m4,m3,m4));
Kvals_D2_WT=array(0,c(m3,m4,m3,m4));

#Big Cypress
for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_BC)
				D2_BC[Plop[k,1:m4],Plop[i,j]]=kvals
				Kvals_D2_BC[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		
D2_BC=D2_BC*h3*h4 #Multiply D2 by widths
save(D2_BC, file="./BC/D2_BC.RData")
save(Kvals_D2_BC, file="./BC/Kvals_D2_BC")
rm(D2_BC, Kvals_D2_BC)

#Cape Canaveral
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_CC)
      D2_CC[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_CC[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_CC=D2_CC*h3*h4 #Multiply D2 by widths
save(D2_CC, file="./CC/D2_CC.RData")
save(Kvals_D2_CC, file="./CC/Kvals_D2_CC")
rm(D2_CC, Kvals_D2_CC)

#Chekika
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_C)
      D2_C[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_C[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_C=D2_C*h3*h4 #Multiply D2 by widths
save(D2_C, file="./C/D2_C.RData")
save(Kvals_D2_C, file="./C/Kvals_D2_C")
rm(D2_C, Kvals_D2_C)

#Fort Pierce
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_FP)
      D2_FP[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_FP[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_FP=D2_FP*h3*h4 #Multiply D2 by widths
save(D2_FP, file="./FP/D2_FP.RData")
save(Kvals_D2_FP, file="./FP/Kvals_D2_FP")
rm(D2_FP, Kvals_D2_FP)

#Punta Gorda
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_PG)
      D2_PG[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_PG[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_PG=D2_PG*h3*h4 #Multiply D2 by widths
save(D2_PG, file="./PG/D2_PG.RData")
save(Kvals_D2_PG, file="./PG/Kvals_D2_PG")
rm(D2_PG, Kvals_D2_PG)

#Wild Turkey
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_WT)
      D2_WT[Plop[k,1:m4],Plop[i,j]]=kvals
      Kvals_D2_WT[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		
D2_WT=D2_WT*h3*h4 #Multiply D2 by widths
save(D2_WT, file="./WT/D2_WT.RData")
save(Kvals_D2_WT, file="./WT/Kvals_D2_WT")
rm(D2_WT, Kvals_D2_WT)









###Construct F (Fertility):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

F_BC=matrix(0,m1*m2,m3*m4); 
F_CC=matrix(0,m1*m2,m3*m4);
F_C=matrix(0,m1*m2,m3*m4);
F_FP=matrix(0,m1*m2,m3*m4);
F_PG=matrix(0,m1*m2,m3*m4);
F_WT=matrix(0,m1*m2,m3*m4);

Kvals_F_BC=array(0, c(m1, m2, m3, m4))
Kvals_F_CC=array(0, c(m1, m2, m3, m4))
Kvals_F_C=array(0, c(m1, m2, m3, m4))
Kvals_F_FP=array(0, c(m1, m2, m3, m4))
Kvals_F_PG=array(0, c(m1, m2, m3, m4))
Kvals_F_WT=array(0, c(m1, m2, m3, m4))


#Big Cypress
for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_BC)
			F_BC[Plop1[k, 1:m2], Plop2[i,j]]=kvals
			Kvals_F_BC[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}
F_BC=F_BC*h1*h2
save(F_BC, file="./BC/F_BC.RData")
save(Kvals_F_BC, file="./BC/Kvals_F_BC.RData")
rm(F_BC, Kvals_F_BC)

#Cape Canaveral
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_CC)
      F_CC[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_CC[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_CC=F_CC*h1*h2
save(F_CC, file="./CC/F_CC.RData")
save(Kvals_F_CC, file="./CC/Kvals_F_CC.RData")
rm(F_CC, Kvals_F_CC)

#Chekika
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_C)
      F_C[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_C[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_C=F_C*h1*h2
save(F_C, file="./C/F_C.RData")
save(Kvals_F_C, file="./C/Kvals_F_C.RData")
rm(F_C, Kvals_F_C)

#Fort Pierce
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_FP)
      F_FP[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_FP[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_FP=F_FP*h1*h2
save(F_FP, file="./FP/F_FP.RData")
save(Kvals_F_FP, file="./FP/Kvals_F_FP.RData")
rm(F_FP, Kvals_F_FP)

#Punta Gorda
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_PG)
      F_PG[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_PG[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_PG=F_PG*h1*h2
save(F_PG, file="./PG/F_PG.RData")
save(Kvals_F_PG, file="./PG/Kvals_F_PG.RData")
rm(F_PG, Kvals_F_PG)

#Wild Turkey
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_WT)
      F_WT[Plop1[k, 1:m2], Plop2[i,j]]=kvals
      Kvals_F_WT[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
F_WT=F_WT*h1*h2
save(F_WT, file="./WT/F_WT.RData")
save(Kvals_F_WT, file="./WT/Kvals_F_WT.RData")
rm(F_WT, Kvals_F_WT)



###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

G_BC=matrix(0,m3*m4,m1*m2); 
G_CC=matrix(0,m3*m4,m1*m2); 
G_C=matrix(0,m3*m4,m1*m2); 
G_FP=matrix(0,m3*m4,m1*m2); 
G_PG=matrix(0,m3*m4,m1*m2);
G_WT=matrix(0,m3*m4,m1*m2);

Kvals_G_BC=array(0, c(m3, m4, m1, m2))
Kvals_G_CC=array(0, c(m3, m4, m1, m2))
Kvals_G_C=array(0, c(m3, m4, m1, m2))
Kvals_G_FP=array(0, c(m3, m4, m1, m2))
Kvals_G_PG=array(0, c(m3, m4, m1, m2))
Kvals_G_WT=array(0, c(m3, m4, m1, m2))

#Big Cypress
for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_BC)
			G_BC[Plop1[k, 1:m4], Plop2[i,j]]=kvals
			Kvals_G_BC[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}
G_BC=G_BC*h3*h4
save(G_BC, file="./BC/G_BC.RData")
save(Kvals_G_BC, file="./BC/Kvals_G_BC.RData")
rm(G_BC, Kvals_G_BC)

#Cape Canaveral
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_CC)
      G_CC[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_CC[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_CC=G_CC*h3*h4
save(G_CC, file="./CC/G_CC.RData")
save(Kvals_G_CC, file="./CC/Kvals_G_CC.RData")
rm(G_CC, Kvals_G_CC)

#Chekika
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_C)
      G_C[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_C[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_C=G_C*h3*h4
save(G_C, file="./C/G_C.RData")
save(Kvals_G_C, file="./C/Kvals_G_C.RData")
rm(G_C, Kvals_G_C)

#Fort Pierce
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_FP)
      G_FP[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_FP[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_FP=G_FP*h3*h4
save(G_FP, file="./FP/G_FP.RData")
save(Kvals_G_FP, file="./FP/Kvals_G_FP.RData")
rm(G_FP, Kvals_G_FP)

# Punta Gorda
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_PG)
      G_PG[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_PG[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_PG=G_PG*h3*h4
save(G_PG, file="./PG/G_PG.RData")
save(Kvals_G_PG, file="./PG/Kvals_G_PG.RData")
rm(G_PG, Kvals_G_PG)

#Wild Turkey
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_WT)
      G_WT[Plop1[k, 1:m4], Plop2[i,j]]=kvals
      Kvals_G_WT[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}
G_WT=G_WT*h3*h4
save(G_WT, file="./WT/G_WT.RData")
save(Kvals_G_WT, file="./WT/Kvals_G_WT.RData")
rm(G_WT, Kvals_G_WT)




#Putting the matrices together:

#Big Cypress
rm(Plop, Plop1, Plop2, i, j, k, kvals ) #clear out workspace to maximize space available 
load("./BC/D1_BC.RData")
load("./BC/G_BC.RData")
left_side<-rbind(D1_BC, G_BC)
rm(D1_BC, G_BC)
load("./BC/F_BC.RData")
load("./BC/D2_BC.RData")
right_side<-rbind(F_BC, D2_BC)
rm(F_BC, D2_BC)
A_BC<-cbind(left_side, right_side)
save(A_BC, file="./BC/A_BC.RData")

#Cape Canaveral
load("./CC/D1_CC.RData")
load("./CC/G_CC.RData")
left_side<-rbind(D1_CC, G_CC)
rm(D1_CC, G_CC)
load("./CC/F_CC.RData")
load("./CC/D2_CC.RData")
right_side<-rbind(F_CC, D2_CC)
rm(F_CC, D2_CC)
A_CC<-cbind(left_side, right_side)
save(A_CC, file="./CC/A_CC.RData")

#Chekika
load("./C/D1_C.RData")
load("./C/G_C.RData")
left_side<-rbind(D1_C, G_C)
rm(D1_C, G_C)
load("./C/F_C.RData")
load("./C/D2_C.RData")
right_side<-rbind(F_C, D2_C)
rm(F_C, D2_C)
A_C<-cbind(left_side, right_side)
save(A_C, file="./C/A_C.RData")

#Fort Pierce
load("./FP/D1_FP.RData")
load("./FP/G_FP.RData")
left_side<-rbind(D1_FP, G_FP)
rm(D1_FP, G_FP)
load("./FP/F_FP.RData")
load("./FP/D2_FP.RData")
right_side<-rbind(F_FP, D2_FP)
rm(F_FP, D2_FP)
A_FP<-cbind(left_side, right_side)
save(A_FP, file="./FP/A_FP.RData")

#Punta Gorda
load("./PG/D1_PG.RData")
load("./PG/G_PG.RData")
left_side<-rbind(D1_PG, G_PG)
rm(D1_PG, G_PG)
load("./PG/F_PG.RData")
load("./PG/D2_PG.RData")
right_side<-rbind(F_PG, D2_PG)
rm(F_PG, D2_PG)
A_PG<-cbind(left_side, right_side)
save(A_PG, file="./PG/A_PG.RData")

#Wild Turkey
load("./WT/D1_WT.RData")
load("./WT/G_WT.RData")
left_side<-rbind(D1_WT, G_WT)
rm(D1_WT, G_WT)
load("./WT/F_WT.RData")
load("./WT/D2_WT.RData")
right_side<-rbind(F_WT, D2_WT)
rm(F_WT, D2_WT)
A_WT<-cbind(left_side, right_side)
save(A_WT, file="./WT/A_WT.RData")




#============================================================================# 
#  Find lambda, w by iteration 
#
#  Note: the Matrix package is used to iterate more quickly. The check  
#  for convergence requires extracting matrix entries via the @x slot
#  of a Matrix object. Matrix is S4-style -- see ?Matrix. 
#============================================================================# 

#Big Cypress
load("./BC/A_BC.RData")
A2=Matrix(A_BC); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

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
qmax=sum(abs(lam*nt-A_BC%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_BC<-lam.stable
stable.dist_BC<-stable.dist
save(lam.stable_BC, file="./BC/lam.stable_BC.RData")
save(stable.dist_BC, file="./BC/stable.dist_BC.RData")


##Alternate calculation of rv:
# 
# A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
# 
# qmax=1000; lam=1; 
# while(qmax>tol) {
# 	nt1=A3%*%nt;
# 	qmax=sum(abs((nt1-lam*nt)@x));  
# 	lam=sum(nt1@x); 
# 	nt@x=(nt1@x)/lam;   
# 	cat(lam,qmax,"\n");
# } 
# nt=matrix(nt@x,m1*m2+m3*m4,1); 
# #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
# #stable.dist=nt
# lam.stable.t=lam; 
# 
# 
# # Check that the @bits worked as intended.   
# qmax=sum(abs(lam*nt-t(A_overall)%*%nt)); 
# cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






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

v_BC<-v
save(v_BC, file="./BC/v_BC.RData")

rm(A_overall, A2, A3, nt, nt1, stable.dist, stable.dist_overall, v, v_overall, vt, vt1,
   lam.stable, lam.stable_overall, lam.stable.t, qmax)


#Cape Canaveral
load("./CC/A_CC.RData")
A2=Matrix(A_CC); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

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
qmax=sum(abs(lam*nt-A_BC%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_CC<-lam.stable
stable.dist_CC<-stable.dist
save(lam.stable_CC, file="./CC/lam.stable_CC.RData")
save(stable.dist_CC, file="./CC/stable.dist_CC.RData")


##Alternate calculation of rv:
# 
# A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
# 
# qmax=1000; lam=1; 
# while(qmax>tol) {
# 	nt1=A3%*%nt;
# 	qmax=sum(abs((nt1-lam*nt)@x));  
# 	lam=sum(nt1@x); 
# 	nt@x=(nt1@x)/lam;   
# 	cat(lam,qmax,"\n");
# } 
# nt=matrix(nt@x,m1*m2+m3*m4,1); 
# #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
# #stable.dist=nt
# lam.stable.t=lam; 
# 
# 
# # Check that the @bits worked as intended.   
# qmax=sum(abs(lam*nt-t(A_overall)%*%nt)); 
# cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






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

v_CC<-v
save(v_CC, file="./CC/v_CC.RData")

#Chekika
load("./C/A_C.RData")
A2=Matrix(A_C); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

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
qmax=sum(abs(lam*nt-A_C%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_C<-lam.stable
stable.dist_C<-stable.dist
save(lam.stable_C, file="./C/lam.stable_C.RData")
save(stable.dist_C, file="./C/stable.dist_C.RData")


##Alternate calculation of rv:
# 
# A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
# 
# qmax=1000; lam=1; 
# while(qmax>tol) {
# 	nt1=A3%*%nt;
# 	qmax=sum(abs((nt1-lam*nt)@x));  
# 	lam=sum(nt1@x); 
# 	nt@x=(nt1@x)/lam;   
# 	cat(lam,qmax,"\n");
# } 
# nt=matrix(nt@x,m1*m2+m3*m4,1); 
# #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
# #stable.dist=nt
# lam.stable.t=lam; 
# 
# 
# # Check that the @bits worked as intended.   
# qmax=sum(abs(lam*nt-t(A_overall)%*%nt)); 
# cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






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

v_C<-v
save(v_C, file="./C/v_C.RData")


#Fort Pierce
load("./FP/A_FP.RData")
A2=Matrix(A_FP); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

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
qmax=sum(abs(lam*nt-A_FP%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_FP<-lam.stable
stable.dist_FP<-stable.dist
save(lam.stable_FP, file="./FP/lam.stable_FP.RData")
save(stable.dist_FP, file="./FP/stable.dist_FP.RData")


##Alternate calculation of rv:
# 
# A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
# 
# qmax=1000; lam=1; 
# while(qmax>tol) {
# 	nt1=A3%*%nt;
# 	qmax=sum(abs((nt1-lam*nt)@x));  
# 	lam=sum(nt1@x); 
# 	nt@x=(nt1@x)/lam;   
# 	cat(lam,qmax,"\n");
# } 
# nt=matrix(nt@x,m1*m2+m3*m4,1); 
# #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
# #stable.dist=nt
# lam.stable.t=lam; 
# 
# 
# # Check that the @bits worked as intended.   
# qmax=sum(abs(lam*nt-t(A_overall)%*%nt)); 
# cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






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

v_FP<-v
save(v_FP, file="./FP/v_FP.RData")



#Punta Gorda
load("./PG/A_PG.RData")
A2=Matrix(A_PG); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

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
qmax=sum(abs(lam*nt-A_PG%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_PG<-lam.stable
stable.dist_PG<-stable.dist
save(lam.stable_PG, file="./PG/lam.stable_PG.RData")
save(stable.dist_PG, file="./PG/stable.dist_PG.RData")


##Alternate calculation of rv:
# 
# A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
# 
# qmax=1000; lam=1; 
# while(qmax>tol) {
# 	nt1=A3%*%nt;
# 	qmax=sum(abs((nt1-lam*nt)@x));  
# 	lam=sum(nt1@x); 
# 	nt@x=(nt1@x)/lam;   
# 	cat(lam,qmax,"\n");
# } 
# nt=matrix(nt@x,m1*m2+m3*m4,1); 
# #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
# #stable.dist=nt
# lam.stable.t=lam; 
# 
# 
# # Check that the @bits worked as intended.   
# qmax=sum(abs(lam*nt-t(A_overall)%*%nt)); 
# cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






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

v_PG<-v
save(v_PG, file="./PG/v_PG.RData")

#Wild Turkey
load("./WT/A_WT.RData")
A2=Matrix(A_WT); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

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
qmax=sum(abs(lam*nt-A_WT%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 
lam.stable_WT<-lam.stable
stable.dist_WT<-stable.dist
save(lam.stable_WT, file="./WT/lam.stable_WT.RData")
save(stable.dist_WT, file="./WT/stable.dist_WT.RData")


##Alternate calculation of rv:
# 
# A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
# 
# qmax=1000; lam=1; 
# while(qmax>tol) {
# 	nt1=A3%*%nt;
# 	qmax=sum(abs((nt1-lam*nt)@x));  
# 	lam=sum(nt1@x); 
# 	nt@x=(nt1@x)/lam;   
# 	cat(lam,qmax,"\n");
# } 
# nt=matrix(nt@x,m1*m2+m3*m4,1); 
# #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
# #stable.dist=nt
# lam.stable.t=lam; 
# 
# 
# # Check that the @bits worked as intended.   
# qmax=sum(abs(lam*nt-t(A_overall)%*%nt)); 
# cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






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

v_WT<-v
save(v_WT, file="./WT/v_WT.RData")




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

#Big Cypress
load("./BC/v_BC.RData")
load("./BC/stable.dist_BC.RData")
load("./BC/A_BC.RData")
load("./BC/lam.stable_BC.RData")

v.dot.w_BC<-sum(t(v_BC)%*%stable.dist_BC)
norm_v_BC<-v_BC/v.dot.w_BC
check_BC<-t(norm_v_BC)%*%stable.dist_BC #should be 1

rv_BC<-norm_v_BC

sens_BC<-norm_v_BC%*%t(stable.dist_BC)
elas_BC<-sens_BC*A_BC/lam.stable_BC
save(sens_BC, file="./BC/sens_BC.RData")
save(elas_BC, file="./BC/elas_BC.RData")

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_BC<-elas_BC[1:(m1*m2), 1:(m1*m2)]
elas_G_BC<-elas_BC[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_BC<-elas_BC[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_BC<-elas_BC[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_BC<-sens_BC[1:(m1*m2), 1:(m1*m2)]
sens_G_BC<-sens_BC[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_BC<-sens_BC[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_BC<-sens_BC[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_BC, file="./BC/elas_D1_BC.RData")
save(elas_G_BC, file="./BC/elas_G_BC.RData")
save(elas_F_BC, file="./BC/elas_F_BC.RData")
save(elas_D2_BC, file="./BC/elas_D2_BC.RData")
save(sens_D1_BC, file="./BC/sens_D1_BC.RData")
save(sens_G_BC, file="./BC/sens_G_BC.RData")
save(sens_F_BC, file="./BC/sens_F_BC.RData")
save(sens_D2_BC, file="./BC/sens_D2_BC.RData")


#Cape Canaveral
load("./CC/v_CC.RData")
load("./CC/stable.dist_CC.RData")
load("./CC/A_CC.RData")
load("./CC/lam.stable_CC.RData")

v.dot.w_CC<-sum(t(v_CC)%*%stable.dist_CC)
norm_v_CC<-v_CC/v.dot.w_CC
check_CC<-t(norm_v_CC)%*%stable.dist_CC #should be 1

rv_CC<-norm_v_CC

sens_CC<-norm_v_CC%*%t(stable.dist_CC)
elas_CC<-sens_CC*A_CC/lam.stable_CC
save(sens_CC, file="./CC/sens_CC.RData")
save(elas_CC, file="./CC/elas_CC.RData")

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_CC<-elas_CC[1:(m1*m2), 1:(m1*m2)]
elas_G_CC<-elas_CC[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_CC<-elas_CC[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_CC<-elas_CC[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_CC<-sens_CC[1:(m1*m2), 1:(m1*m2)]
sens_G_CC<-sens_CC[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_CC<-sens_CC[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_CC<-sens_CC[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_CC, file="./CC/elas_D1_CC.RData")
save(elas_G_CC, file="./CC/elas_G_CC.RData")
save(elas_F_CC, file="./CC/elas_F_CC.RData")
save(elas_D2_CC, file="./CC/elas_D2_CC.RData")
save(sens_D1_CC, file="./CC/sens_D1_CC.RData")
save(sens_G_CC, file="./CC/sens_G_CC.RData")
save(sens_F_CC, file="./CC/sens_F_CC.RData")
save(sens_D2_CC, file="./CC/sens_D2_CC.RData")

#Chekika
load("./C/v_C.RData")
load("./C/stable.dist_C.RData")
load("./C/A_C.RData")
load("./C/lam.stable_C.RData")

v.dot.w_C<-sum(t(v_C)%*%stable.dist_C)
norm_v_C<-v_C/v.dot.w_C
check_C<-t(norm_v_C)%*%stable.dist_C #should be 1

rv_C<-norm_v_C

sens_C<-norm_v_C%*%t(stable.dist_C)
elas_C<-sens_C*A_C/lam.stable_C
save(sens_C, file="./C/sens_C.RData")
save(elas_C, file="./C/elas_C.RData")

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_C<-elas_C[1:(m1*m2), 1:(m1*m2)]
elas_G_C<-elas_C[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_C<-elas_C[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_C<-elas_C[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_C<-sens_C[1:(m1*m2), 1:(m1*m2)]
sens_G_C<-sens_C[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_C<-sens_C[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_C<-sens_C[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_C, file="./C/elas_D1_C.RData")
save(elas_G_C, file="./C/elas_G_C.RData")
save(elas_F_C, file="./C/elas_F_C.RData")
save(elas_D2_C, file="./C/elas_D2_C.RData")
save(sens_D1_C, file="./C/sens_D1_C.RData")
save(sens_G_C, file="./C/sens_G_C.RData")
save(sens_F_C, file="./C/sens_F_C.RData")
save(sens_D2_C, file="./C/sens_D2_C.RData")

#Fort Pierce
load("./FP/v_FP.RData")
load("./FP/stable.dist_FP.RData")
load("./FP/A_FP.RData")
load("./FP/lam.stable_FP.RData")

v.dot.w_FP<-sum(t(v_FP)%*%stable.dist_FP)
norm_v_FP<-v_FP/v.dot.w_FP
check_FP<-t(norm_v_FP)%*%stable.dist_FP #should be 1

rv_FP<-norm_v_FP

sens_FP<-norm_v_FP%*%t(stable.dist_FP)
elas_FP<-sens_FP*A_FP/lam.stable_FP
save(sens_FP, file="./FP/sens_FP.RData")
save(elas_FP, file="./FP/elas_FP.RData")

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_FP<-elas_FP[1:(m1*m2), 1:(m1*m2)]
elas_G_FP<-elas_FP[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_FP<-elas_FP[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_FP<-elas_FP[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_FP<-sens_FP[1:(m1*m2), 1:(m1*m2)]
sens_G_FP<-sens_FP[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_FP<-sens_FP[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_FP<-sens_FP[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_FP, file="./FP/elas_D1_FP.RData")
save(elas_G_FP, file="./FP/elas_G_FP.RData")
save(elas_F_FP, file="./FP/elas_F_FP.RData")
save(elas_D2_FP, file="./FP/elas_D2_FP.RData")
save(sens_D1_FP, file="./FP/sens_D1_FP.RData")
save(sens_G_FP, file="./FP/sens_G_FP.RData")
save(sens_F_FP, file="./FP/sens_F_FP.RData")
save(sens_D2_FP, file="./FP/sens_D2_FP.RData")


#Punta Gorda
load("./PG/v_PG.RData")
load("./PG/stable.dist_PG.RData")
load("./PG/A_PG.RData")
load("./PG/lam.stable_PG.RData")

v.dot.w_PG<-sum(t(v_PG)%*%stable.dist_PG)
norm_v_PG<-v_PG/v.dot.w_PG
check_PG<-t(norm_v_PG)%*%stable.dist_PG #should be 1

rv_PG<-norm_v_PG

sens_PG<-norm_v_PG%*%t(stable.dist_PG)
elas_PG<-sens_PG*A_PG/lam.stable_PG
save(sens_PG, file="./PG/sens_PG.RData")
save(elas_PG, file="./PG/elas_PG.RData")

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_PG<-elas_PG[1:(m1*m2), 1:(m1*m2)]
elas_G_PG<-elas_PG[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_PG<-elas_PG[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_PG<-elas_PG[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_PG<-sens_PG[1:(m1*m2), 1:(m1*m2)]
sens_G_PG<-sens_PG[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_PG<-sens_PG[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_PG<-sens_PG[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_PG, file="./PG/elas_D1_PG.RData")
save(elas_G_PG, file="./PG/elas_G_PG.RData")
save(elas_F_PG, file="./PG/elas_F_PG.RData")
save(elas_D2_PG, file="./PG/elas_D2_PG.RData")
save(sens_D1_PG, file="./PG/sens_D1_PG.RData")
save(sens_G_PG, file="./PG/sens_G_PG.RData")
save(sens_F_PG, file="./PG/sens_F_PG.RData")
save(sens_D2_PG, file="./PG/sens_D2_PG.RData")

#Wild Turkey
load("./WT/v_WT.RData")
load("./WT/stable.dist_WT.RData")
load("./WT/A_WT.RData")
load("./WT/lam.stable_WT.RData")

v.dot.w_WT<-sum(t(v_WT)%*%stable.dist_WT)
norm_v_WT<-v_WT/v.dot.w_WT
check_WT<-t(norm_v_WT)%*%stable.dist_WT #should be 1

rv_WT<-norm_v_WT

sens_WT<-norm_v_WT%*%t(stable.dist_WT)
elas_WT<-sens_WT*A_WT/lam.stable_WT
save(sens_WT, file="./WT/sens_WT.RData")
save(elas_WT, file="./WT/elas_WT.RData")

#plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1_WT<-elas_WT[1:(m1*m2), 1:(m1*m2)]
elas_G_WT<-elas_WT[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F_WT<-elas_WT[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2_WT<-elas_WT[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1_WT<-sens_WT[1:(m1*m2), 1:(m1*m2)]
sens_G_WT<-sens_WT[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F_WT<-sens_WT[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2_WT<-sens_WT[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]


save(elas_D1_WT, file="./WT/elas_D1_WT.RData")
save(elas_G_WT, file="./WT/elas_G_WT.RData")
save(elas_F_WT, file="./WT/elas_F_WT.RData")
save(elas_D2_WT, file="./WT/elas_D2_WT.RData")
save(sens_D1_WT, file="./WT/sens_D1_WT.RData")
save(sens_G_WT, file="./WT/sens_G_WT.RData")
save(sens_F_WT, file="./WT/sens_F_WT.RData")
save(sens_D2_WT, file="./WT/sens_D2_WT.RData")




########
###Now go back to unpack the elasticity matrices into components:
########

load("./BC/elas_D1_BC.RData")
load("./CC/elas_D1_CC.RData")
load("./C/elas_D1_C.RData")
load("./FP/elas_D1_FP.RData")
load("./PG/elas_D1_PG.RData")
load("./WT/elas_D1_WT.RData")


plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in array
Plop=outer(1:m1,1:m2,plop); 

Kvals_elas_D1_BC=array(0,c(m1,m2,m1,m2));  
Kvals_elas_D1_CC=array(0,c(m1,m2,m1,m2)); 
Kvals_elas_D1_C=array(0,c(m1,m2,m1,m2)); 
Kvals_elas_D1_FP=array(0,c(m1,m2,m1,m2)); 
Kvals_elas_D1_PG=array(0,c(m1,m2,m1,m2)); 
Kvals_elas_D1_WT=array(0,c(m1,m2,m1,m2));  

#Big Cypress
for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals= elas_D1_BC[Plop[k,1:m2],Plop[i,j]]
				Kvals_elas_D1_BC[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}
save(Kvals_elas_D1_BC, file="./BC/Kvals_elas_D1_BC.RData")

#Cape Canaveral
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals= elas_D1_CC[Plop[k,1:m2],Plop[i,j]]
      Kvals_elas_D1_CC[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
save(Kvals_elas_D1_CC, file="./CC/Kvals_elas_D1_CC.RData")

#Chekika
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals= elas_D1_C[Plop[k,1:m2],Plop[i,j]]
      Kvals_elas_D1_C[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
save(Kvals_elas_D1_C, file="./C/Kvals_elas_D1_C.RData")

#Fort Pierce
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals= elas_D1_FP[Plop[k,1:m2],Plop[i,j]]
      Kvals_elas_D1_FP[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
save(Kvals_elas_D1_FP, file="./FP/Kvals_elas_D1_FP.RData")

#Punta Gorda
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals= elas_D1_PG[Plop[k,1:m2],Plop[i,j]]
      Kvals_elas_D1_PG[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
save(Kvals_elas_D1_PG, file="./PG/Kvals_elas_D1_PG.RData")

#Wild Turkey
for(i in 1:m1){
  for(j in 1:m2){
    for(k in 1:m1){
      kvals= elas_D1_WT[Plop[k,1:m2],Plop[i,j]]
      Kvals_elas_D1_WT[k,1:m2,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}
save(Kvals_elas_D1_WT, file="./WT/Kvals_elas_D1_WT.RData")



rm(Kvals_elas_D1_BC, Kvals_elas_D1_CC, Kvals_elas_D1_C, Kvals_elas_D1_FP, Kvals_elaS_D1_PG,
   Kvals_elas_D1_WT, i, j, k, kvals, 
   elas_D1_BC, elas_D1_CC, elas_D1_C, elas_D1_FP, elas_D1_PG, elas_D1_WT, Plop)


###Construct D2 (Large Domain):

load("./BC/elas_D2_BC.RData")
load("./CC/elas_D2_CC.RData")
load("./C/elas_D2_C.RData")
load("./FP/elas_D2_FP.RData")
load("./PG/elas_D2_PG.RData")
load("./WT/elas_D2_WT.RData")

plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_elas_D2_BC=array(0,c(m3,m4,m3,m4));  
Kvals_elas_D2_CC=array(0,c(m3,m4,m3,m4));
Kvals_elas_D2_C=array(0,c(m3,m4,m3,m4));
Kvals_elas_D2_FP=array(0,c(m3,m4,m3,m4));
Kvals_elas_D2_PG=array(0,c(m3,m4,m3,m4));
Kvals_elas_D2_WT=array(0,c(m3,m4,m3,m4));

#Big Cypress
for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= elas_D2_BC[Plop[k,1:m4],Plop[i,j]]
				Kvals_elas_D2_BC[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		

save(Kvals_elas_D2_BC, file="./BC/Kvals_elas_D2_BC.RData")

#Cape Canaveral
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals= elas_D2_CC[Plop[k,1:m4],Plop[i,j]]
      Kvals_elas_D2_CC[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		

save(Kvals_elas_D2_CC, file="./CC/Kvals_elas_D2_CC.RData")

#Chekika
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals= elas_D2_C[Plop[k,1:m4],Plop[i,j]]
      Kvals_elas_D2_C[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		

save(Kvals_elas_D2_C, file="./C/Kvals_elas_D2_C.RData")

#Fort Pierce
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals= elas_D2_FP[Plop[k,1:m4],Plop[i,j]]
      Kvals_elas_D2_FP[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		

save(Kvals_elas_D2_FP, file="./FP/Kvals_elas_D2_FP.RData")

#Punta Gorda
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals= elas_D2_PG[Plop[k,1:m4],Plop[i,j]]
      Kvals_elas_D2_PG[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		

save(Kvals_elas_D2_PG, file="./PG/Kvals_elas_D2_PG.RData")

#Wild Turkey
for(i in 1:m3){
  for(j in 1:m4){
    for(k in 1:m3){
      kvals= elas_D2_WT[Plop[k,1:m4],Plop[i,j]]
      Kvals_elas_D2_WT[k,1:m4,i,j]=kvals
      
    }}
  cat(i,"\n"); 
}		

save(Kvals_elas_D2_WT, file="./WT/Kvals_elas_D2_WT.RData")


rm(elas_D2_BC, elas_D2_C, elas_D2_CC, elas_D2_FP, elas_D2_PG, elas_D2_WT, Plop)

###Construct F (Fecundity):
load("./BC/elas_F_BC.RData")
load("./CC/elas_F_CC.RData")
load("./C/elas_F_C.RData")
load("./FP/elas_F_FP.RData")
load("./PG/elas_F_PG.RData")
load("./WT/elas_F_WT.RData")

plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_elas_F_BC=array(0, c(m1, m2, m3, m4))
Kvals_elas_F_CC=array(0, c(m1, m2, m3, m4))
Kvals_elas_F_C=array(0, c(m1, m2, m3, m4))
Kvals_elas_F_FP=array(0, c(m1, m2, m3, m4))
Kvals_elas_F_PG=array(0, c(m1, m2, m3, m4))
Kvals_elas_F_WT=array(0, c(m1, m2, m3, m4))


#Big Cypress
for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=elas_F_BC[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_elas_F_BC[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}
save(Kvals_elas_F_BC, file="./BC/Kvals_elas_F_BC.RData")

#Cape Canaveral
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=elas_F_CC[Plop1[k, 1:m2], Plop2[i,j]]
      Kvals_elas_F_CC[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
save(Kvals_elas_F_CC, file="./CC/Kvals_elas_F_CC.RData")

#Chekika
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=elas_F_C[Plop1[k, 1:m2], Plop2[i,j]]
      Kvals_elas_F_C[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
save(Kvals_elas_F_C, file="./C/Kvals_elas_F_C.RData")

#Fort Pierce
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=elas_F_FP[Plop1[k, 1:m2], Plop2[i,j]]
      Kvals_elas_F_FP[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
save(Kvals_elas_F_FP, file="./FP/Kvals_elas_F_FP.RData")

#Punta Gorda
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=elas_F_PG[Plop1[k, 1:m2], Plop2[i,j]]
      Kvals_elas_F_PG[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
save(Kvals_elas_F_PG, file="./PG/Kvals_elas_F_PG.RData")

#Wild Turkey
for(i in 1:m3) {
  for (j in 1:m4) {
    for (k in 1:m1) {
      kvals=elas_F_WT[Plop1[k, 1:m2], Plop2[i,j]]
      Kvals_elas_F_WT[k, 1:m2, i, j]=kvals
    }}
  cat(i, "\n");
}
save(Kvals_elas_F_WT, file="./WT/Kvals_elas_F_WT.RData")






###Construct G (Graduation):

load("./BC/elas_G_BC.RData")
load("./CC/elas_G_CC.RData")
load("./C/elas_G_C.RData")
load("./FP/elas_G_FP.RData")
load("./PG/elas_G_PG.RData")
load("./WT/elas_G_WT.RData")

plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_elas_G_BC=array(0, c(m3, m4, m1, m2))
Kvals_elas_G_CC=array(0, c(m3, m4, m1, m2))
Kvals_elas_G_C=array(0, c(m3, m4, m1, m2))
Kvals_elas_G_FP=array(0, c(m3, m4, m1, m2))
Kvals_elas_G_PG=array(0, c(m3, m4, m1, m2))
Kvals_elas_G_WT=array(0, c(m3, m4, m1, m2))

#Big Cypress
for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=elas_G_BC[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_elas_G_BC[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}

save(Kvals_elas_G_BC, file="./BC/Kvals_elas_G_BC.RData")

#Cape Canaveral
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=elas_G_CC[Plop1[k, 1:m4], Plop2[i,j]]
      Kvals_elas_G_CC[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}

save(Kvals_elas_G_CC, file="./CC/Kvals_elas_G_CC.RData")
#Chekika
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=elas_G_C[Plop1[k, 1:m4], Plop2[i,j]]
      Kvals_elas_G_C[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}

save(Kvals_elas_G_C, file="./C/Kvals_elas_G_C.RData")

#Fort Pierce
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=elas_G_FP[Plop1[k, 1:m4], Plop2[i,j]]
      Kvals_elas_G_FP[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}

save(Kvals_elas_G_FP, file="./FP/Kvals_elas_G_FP.RData")

#Punta Gorda
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=elas_G_PG[Plop1[k, 1:m4], Plop2[i,j]]
      Kvals_elas_G_PG[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}

save(Kvals_elas_G_PG, file="./PG/Kvals_elas_G_PG.RData")

#Wild Turkey
for(i in 1:m1) {
  for (j in 1:m2) {
    for (k in 1:m3) {
      kvals=elas_G_WT[Plop1[k, 1:m4], Plop2[i,j]]
      Kvals_elas_G_WT[k, 1:m4, i, j]=kvals
    }}
  cat(i, "\n");
}

save(Kvals_elas_G_WT, file="./WT/Kvals_elas_G_WT.RData")




#####Create Total Elasticity Pieces 
load("./BC/Kvals_elas_D1_BC.RData")
load("./BC/Kvals_elas_D2_BC.RData")
load("./BC/Kvals_elas_F_BC.RData")
load("./BC/Kvals_elas_G_BC.RData")
##############
#Big Cypress
total.elas_D1_BC<-apply(Kvals_elas_D1_BC, c(3,4), sum);
total.elas_D2_BC<-apply(Kvals_elas_D2_BC, c(3,4), sum);
total.elas_F_BC<-apply(Kvals_elas_F_BC, c(3,4), sum);
total.elas_G_BC<-apply(Kvals_elas_G_BC, c(3,4), sum);
save(total.elas_D1_BC, total.elas_D2_BC, total.elas_F_BC, total.elas_G_BC, 
     file="./BC/total.elas.RData")
##############
#Cape Canaveral
load("./CC/Kvals_elas_D1_CC.RData")
load("./CC/Kvals_elas_D2_CC.RData")
load("./CC/Kvals_elas_F_CC.RData")
load("./CC/Kvals_elas_G_CC.RData")

total.elas_D1_CC<-apply(Kvals_elas_D1_CC, c(3,4), sum);
total.elas_D2_CC<-apply(Kvals_elas_D2_CC, c(3,4), sum);
total.elas_F_CC<-apply(Kvals_elas_F_CC, c(3,4), sum);
total.elas_G_CC<-apply(Kvals_elas_G_CC, c(3,4), sum);
save(total.elas_D1_CC, total.elas_D2_CC, total.elas_F_CC, total.elas_G_CC, 
     file="./CC/total.elas.RData")

##############
# Chekika
load("./C/Kvals_elas_D1_C.RData")
load("./C/Kvals_elas_D2_C.RData")
load("./C/Kvals_elas_F_C.RData")
load("./C/Kvals_elas_G_C.RData")

total.elas_D1_C<-apply(Kvals_elas_D1_C, c(3,4), sum);
total.elas_D2_C<-apply(Kvals_elas_D2_C, c(3,4), sum);
total.elas_F_C<-apply(Kvals_elas_F_C, c(3,4), sum);
total.elas_G_C<-apply(Kvals_elas_G_C, c(3,4), sum);
save(total.elas_D1_C, total.elas_D2_C, total.elas_F_C, total.elas_G_C, 
     file="./C/total.elas.RData")

##############
#Fort Pierce
load("./FP/Kvals_elas_D1_FP.RData")
load("./FP/Kvals_elas_D2_FP.RData")
load("./FP/Kvals_elas_F_FP.RData")
load("./FP/Kvals_elas_G_FP.RData")

total.elas_D1_FP<-apply(Kvals_elas_D1_FP, c(3,4), sum);
total.elas_D2_FP<-apply(Kvals_elas_D2_FP, c(3,4), sum);
total.elas_F_FP<-apply(Kvals_elas_F_FP, c(3,4), sum);
total.elas_G_FP<-apply(Kvals_elas_G_FP, c(3,4), sum);
save(total.elas_D1_FP, total.elas_D2_FP, total.elas_F_FP, total.elas_G_FP, 
     file="./FP/total.elas.RData")

##############
#Punta Gorda
load("./PG/Kvals_elas_D1_PG.RData")
load("./PG/Kvals_elas_D2_PG.RData")
load("./PG/Kvals_elas_F_PG.RData")
load("./PG/Kvals_elas_G_PG.RData")

total.elas_D1_PG<-apply(Kvals_elas_D1_PG, c(3,4), sum);
total.elas_D2_PG<-apply(Kvals_elas_D2_PG, c(3,4), sum);
total.elas_F_PG<-apply(Kvals_elas_F_PG, c(3,4), sum);
total.elas_G_PG<-apply(Kvals_elas_G_PG, c(3,4), sum);
save(total.elas_D1_PG, total.elas_D2_PG, total.elas_F_PG, total.elas_G_PG, 
     file="./PG/total.elas.RData")

##############
#Wild Turkey
load("./WT/Kvals_elas_D1_WT.RData")
load("./WT/Kvals_elas_D2_WT.RData")
load("./WT/Kvals_elas_F_WT.RData")
load("./WT/Kvals_elas_G_WT.RData")

total.elas_D1_WT<-apply(Kvals_elas_D1_WT, c(3,4), sum);
total.elas_D2_WT<-apply(Kvals_elas_D2_WT, c(3,4), sum);
total.elas_F_WT<-apply(Kvals_elas_F_WT, c(3,4), sum);
total.elas_G_WT<-apply(Kvals_elas_G_WT, c(3,4), sum);
save(total.elas_D1_WT, total.elas_D2_WT, total.elas_F_WT, total.elas_G_WT, 
     file="./WT/total.elas.RData")






#Separate:
load("./BC/v_BC.RData")
load("./CC/v_CC.RData")
load("./C/v_C.RData")
load("./FP/v_FP.RData")
load("./PG/v_PG.RData")
load("./WT/v_WT.RData")

load("./BC/stable.dist_BC.RData")
load("./CC/stable.dist_CC.RData")
load("./C/stable.dist_C.RData")
load("./FP/stable.dist_FP.RData")
load("./PG/stable.dist_PG.RData")
load("./WT/stable.dist_WT.RData")


repro.val_smalls_BC=matrix(v_BC[1:(m1*m2)], m1,m2);
stable.state_smalls_BC=matrix(stable.dist_BC[1:(m1*m2)], m1, m2);

repro.val_smalls_CC=matrix(v_CC[1:(m1*m2)], m1,m2);
stable.state_smalls_CC=matrix(stable.dist_CC[1:(m1*m2)], m1, m2);

repro.val_smalls_C=matrix(v_C[1:(m1*m2)], m1,m2);
stable.state_smalls_C=matrix(stable.dist_C[1:(m1*m2)], m1, m2);

repro.val_smalls_FP=matrix(v_FP[1:(m1*m2)], m1,m2);
stable.state_smalls_FP=matrix(stable.dist_FP[1:(m1*m2)], m1, m2);

repro.val_smalls_PG=matrix(v_PG[1:(m1*m2)], m1,m2);
stable.state_smalls_PG=matrix(stable.dist_PG[1:(m1*m2)], m1, m2);

repro.val_smalls_WT=matrix(v_WT[1:(m1*m2)], m1,m2);
stable.state_smalls_WT=matrix(stable.dist_WT[1:(m1*m2)], m1, m2);


repro.val_larges_BC=matrix(v_BC[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_BC=matrix(stable.dist_BC[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

repro.val_larges_CC=matrix(v_CC[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_CC=matrix(stable.dist_CC[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

repro.val_larges_C=matrix(v_C[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_C=matrix(stable.dist_C[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

repro.val_larges_FP=matrix(v_FP[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_FP=matrix(stable.dist_FP[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

repro.val_larges_PG=matrix(v_PG[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_PG=matrix(stable.dist_PG[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

repro.val_larges_WT=matrix(v_WT[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_WT=matrix(stable.dist_WT[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

#============================================================================# 
#  Calculate stable state distribution and state-dependent total elasticity 
#============================================================================# 

#Big Cypress
ss_diam1_BC<-apply(stable.state_smalls_BC, 1, sum)
ss_diam2_BC<-apply(stable.state_larges_BC, 1, sum)
ss_diam_BC<-c(ss_diam1_BC, ss_diam2_BC)
y_diam<-c(y1, y3)

save(y_diam, file="./BC/y_diam.RData")
save(ss_diam_BC, file="./BC/ss_diam_BC.Rdata")

#Cape Canaveral
ss_diam1_CC<-apply(stable.state_smalls_CC, 1, sum)
ss_diam2_CC<-apply(stable.state_larges_CC, 1, sum)
ss_diam_CC<-c(ss_diam1_CC, ss_diam2_CC)
y_diam<-c(y1, y3)

save(y_diam, file="./CC/y_diam.RData")
save(ss_diam_CC, file="./CC/ss_diam_CC.Rdata")

#Chekika
ss_diam1_C<-apply(stable.state_smalls_C, 1, sum)
ss_diam2_C<-apply(stable.state_larges_C, 1, sum)
ss_diam_C<-c(ss_diam1_C, ss_diam2_C)
y_diam<-c(y1, y3)

save(y_diam, file="./C/y_diam.RData")
save(ss_diam_C, file="./C/ss_diam_C.Rdata")

#Fort Pierce
ss_diam1_FP<-apply(stable.state_smalls_FP, 1, sum)
ss_diam2_FP<-apply(stable.state_larges_FP, 1, sum)
ss_diam_FP<-c(ss_diam1_FP, ss_diam2_FP)
y_diam<-c(y1, y3)

save(y_diam, file="./FP/y_diam.RData")
save(ss_diam_FP, file="./FP/ss_diam_FP.Rdata")

#Punta Gorda
ss_diam1_PG<-apply(stable.state_smalls_PG, 1, sum)
ss_diam2_PG<-apply(stable.state_larges_PG, 1, sum)
ss_diam_PG<-c(ss_diam1_PG, ss_diam2_PG)
y_diam<-c(y1, y3)

save(y_diam, file="./PG/y_diam.RData")
save(ss_diam_PG, file="./PG/ss_diam_PG.Rdata")

#Wild Turkey
ss_diam1_WT<-apply(stable.state_smalls_WT, 1, sum)
ss_diam2_WT<-apply(stable.state_larges_WT, 1, sum)
ss_diam_WT<-c(ss_diam1_WT, ss_diam2_WT)
y_diam<-c(y1, y3)

save(y_diam, file="./WT/y_diam.RData")
save(ss_diam_WT, file="./WT/ss_diam_WT.Rdata")

## Height
#Big Cypress
ss_height1_BC<-apply(stable.state_smalls_BC, 2, sum)
ss_height2_BC<-apply(stable.state_larges_BC, 2, sum)
ss_height_BC<-c(ss_height1_BC, ss_height2_BC)
y_height<-c(y2, y4)
save(y_height, file="./BC/y_height.RData")
save(ss_height_BC, file="./BC/ss_height_BC.RData")

#Cape Canaveral
ss_height1_CC<-apply(stable.state_smalls_CC, 2, sum)
ss_height2_CC<-apply(stable.state_larges_CC, 2, sum)
ss_height_CC<-c(ss_height1_CC, ss_height2_CC)
y_height<-c(y2, y4)
save(y_height, file="./CC/y_height.RData")
save(ss_height_CC, file="./CC/ss_height_CC.RData")

#Chekika
ss_height1_C<-apply(stable.state_smalls_C, 2, sum)
ss_height2_C<-apply(stable.state_larges_C, 2, sum)
ss_height_C<-c(ss_height1_C, ss_height2_C)
y_height<-c(y2, y4)
save(y_height, file="./C/y_height.RData")
save(ss_height_C, file="./C/ss_height_C.RData")

#Fort Pierce
ss_height1_FP<-apply(stable.state_smalls_FP, 2, sum)
ss_height2_FP<-apply(stable.state_larges_FP, 2, sum)
ss_height_FP<-c(ss_height1_FP, ss_height2_FP)
y_height<-c(y2, y4)
save(y_height, file="./FP/y_height.RData")
save(ss_height_FP, file="./FP/ss_height_FP.RData")

#Punta Gorda
ss_height1_PG<-apply(stable.state_smalls_PG, 2, sum)
ss_height2_PG<-apply(stable.state_larges_PG, 2, sum)
ss_height_PG<-c(ss_height1_PG, ss_height2_PG)
y_height<-c(y2, y4)
save(y_height, file="./PG/y_height.RData")
save(ss_height_PG, file="./PG/ss_height_PG.RData")

#Wild Turkey
ss_height1_WT<-apply(stable.state_smalls_WT, 2, sum)
ss_height2_WT<-apply(stable.state_larges_WT, 2, sum)
ss_height_WT<-c(ss_height1_WT, ss_height2_WT)
y_height<-c(y2, y4)
save(y_height, file="./WT/y_height.RData")
save(ss_height_WT, file="./WT/ss_height_WT.RData")







###Reproductive value
## Diameter 
#Big Cypress
rv_diam1_BC<-apply(repro.val_smalls_BC, 1, sum)
rv_diam2_BC<-apply(repro.val_larges_BC, 1, sum)
rv_diam_BC<-c(rv_diam1_BC, rv_diam2_BC)
save(rv_diam_BC, file="./BC/rv_diam_BC.RData")

#Cape Canaveral
rv_diam1_CC<-apply(repro.val_smalls_CC, 1, sum)
rv_diam2_CC<-apply(repro.val_larges_CC, 1, sum)
rv_diam_CC<-c(rv_diam1_CC, rv_diam2_CC)
save(rv_diam_CC, file="./CC/rv_diam_CC.RData")

#Chekika
rv_diam1_C<-apply(repro.val_smalls_C, 1, sum)
rv_diam2_C<-apply(repro.val_larges_C, 1, sum)
rv_diam_C<-c(rv_diam1_C, rv_diam2_C)
save(rv_diam_C, file="./C/rv_diam_C.RData")

#Fort Pierce
rv_diam1_FP<-apply(repro.val_smalls_FP, 1, sum)
rv_diam2_FP<-apply(repro.val_larges_FP, 1, sum)
rv_diam_FP<-c(rv_diam1_FP, rv_diam2_FP)
save(rv_diam_FP, file="./FP/rv_diam_FP.RData")

#Punta Gorda
rv_diam1_PG<-apply(repro.val_smalls_PG, 1, sum)
rv_diam2_PG<-apply(repro.val_larges_PG, 1, sum)
rv_diam_PG<-c(rv_diam1_PG, rv_diam2_PG)
save(rv_diam_PG, file="./PG/rv_diam_PG.RData")

#Wild Turkey
rv_diam1_WT<-apply(repro.val_smalls_WT, 1, sum)
rv_diam2_WT<-apply(repro.val_larges_WT, 1, sum)
rv_diam_WT<-c(rv_diam1_WT, rv_diam2_WT)
save(rv_diam_WT, file="./WT/rv_diam_WT.RData")

##Height
#Big Cypress
rv_height1_BC<-apply(repro.val_smalls_BC, 2, sum)
rv_height2_BC<-apply(repro.val_larges_BC, 2, sum)
rv_height_BC<-c(rv_height1_BC, rv_height2_BC)
save(rv_height_BC, file="./BC/rv_height_BC.RData")

#Cape Canaveral
rv_height1_CC<-apply(repro.val_smalls_CC, 2, sum)
rv_height2_CC<-apply(repro.val_larges_CC, 2, sum)
rv_height_CC<-c(rv_height1_CC, rv_height2_CC)
save(rv_height_CC, file="./CC/rv_height_CC.RData")

#Chekika
rv_height1_C<-apply(repro.val_smalls_C, 2, sum)
rv_height2_C<-apply(repro.val_larges_C, 2, sum)
rv_height_C<-c(rv_height1_C, rv_height2_C)
save(rv_height_C, file="./C/rv_height_C.RData")

#Fort Pierce
rv_height1_FP<-apply(repro.val_smalls_FP, 2, sum)
rv_height2_FP<-apply(repro.val_larges_FP, 2, sum)
rv_height_FP<-c(rv_height1_FP, rv_height2_FP)
save(rv_height_FP, file="./FP/rv_height_FP.RData")

#Punta Gorda
rv_height1_PG<-apply(repro.val_smalls_PG, 2, sum)
rv_height2_PG<-apply(repro.val_larges_PG, 2, sum)
rv_height_PG<-c(rv_height1_PG, rv_height2_PG)
save(rv_height_PG, file="./PG/rv_height_PG.RData")

#Wild Turkey
rv_height1_WT<-apply(repro.val_smalls_WT, 2, sum)
rv_height2_WT<-apply(repro.val_larges_WT, 2, sum)
rv_height_WT<-c(rv_height1_WT, rv_height2_WT)
save(rv_height_WT, file="./WT/rv_height_WT.RData")








