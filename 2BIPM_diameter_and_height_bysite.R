#***********************************
# IPM: Diameter and Height by Site *
#***********************************

#Constructs a integral projection model using diameter and height as
# two continuous structuring state variables, where population transition
# probabilities are modeled on two size domains.

# This script is modified (by KDE and CCH) from the code provided by
# Ellner and Rees, 2006 that is available for download at: 
# https://www.jstor.org/stable/get_asset/10.1086/499438?supp_index=0&refreqid=excelsior%3Aa8885edea7fee9610a049fc3992a448a

# Here, we expand on Ellner and Rees' code by modeling the population
# dynamics on two continuous size domains (See details in associated manuscript)

# Due to the large size of the objects created, this code can be slow
# to run. To save disk space, whenever possible, large objects are
# saved and then removed from the workspace 

#clear everything
rm(list=ls(all=TRUE))

require(Matrix);
library(MASS)
library(mvtnorm)
library(msm)

#increase memory limit on Windows machines  (This will make everything run faster)
memory.limit(memory.limit() *2^30)

# Set matrix size: 
# m1: # of categories seedling diameter divided into
# m2: # of categories seedling height divided into 
# m3a: # of categories small adult diameter divided into 
# m4a: # of categories small adult height divided into
# m3b: # of categories large adult diameter divided into
# m4b: # of categories large adult height divided into 

#The values used here were chosen by continuing to increase the dimensions of the matrix 
# until the calculated value of lambda stabilized 
m1=10
m2=m1+1
m3a=250
m4a=50

m3b = 120
m4b = 100
m3=m3a + m3b
m4=m4a + m4b
tol=1.e-8; 




# Function to do an image plot of a matrix in the usual orientation,
# A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(200),...) {
  nx=length(x); ny=length(y); 
  x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
  y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
  image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}
# (I) Parameters and demographic functions for computing the kernel #####

#p.vecs are defined externally in script 1demography_models.R
#load('./BC/p.vec_BC.RData')
#load('./CC/p.vec_CC.RData')
#load('./C/p.vec_C.RData')
#load('./FP/p.vec_FP.RData')
#load('./PG/p.vec_PG.RData')
#load('./WT/p.vec_WT.RData')
load("./Overall/p.vec_overall.RData")

#p.vec[1]:seedling survival intercept
#p.vec[2]:seedling survival diameter slope
#p.vec[3]:seedling survival height slope
#p.vec[4]:seedling diameter growth intercept
#p.vec[5]:seedling diameter growth diameter slope
#p.vec[6]:seedling diameter growth height intercept
#p.vec[7]:variance in seedling diameter growth 
#p.vec[8]:seedling height growth intercept
#p.vec[9]:seedling height growth diameter slope
#p.vec[10]:seedling height growth height slope 
#p.vec[11]:variance in seedling height growth 
#p.vec[12]: maturation intercept 
#p.vec[13]: maturation diameter slope 
#p.vec[14]: maturation height slope 
#p.vec[15]: mean graduates diameter
#p.vec[16]: sd graduates diameter 
#p.vec[17]: mean graduates height 
#p.vec[18]: sd graduates height 
#p.vec[19]: D2 survival intercept 
#p.vec[20]: D2 survival diameter slope 
#p.vec[21]: D2 survival height slope 
#p.vec[22]: D2 diameter growth intercept 
#p.vec[23]: D2 diameter growth diameter slope 
#p.vec[24]: D2 diameter growth height slope 
#p.vec[25]: variance in D2 diameter growth 
#p.vec[26]: D2 height growth intercept 
#p.vec[27]: D2 height growth diameter slope 
#p.vec[28]: D2 height growth height slope 
#p.vec[29]: variance in D2 height growth 
#p.vec[30]: probability of reproduction intercept 
#p.vec[31]: probability of reproduction height slope 
#p.vec[32]: fecundity diameter slope 
#p.vec[33]: mean diameter of new  recruits 
#p.vec[34]: sd diameter of new recruits
#p.vec[35]: mean height of new recruits
#p.vec[36]: sd height of new recruits 
#p.vec[37]: tau1 (pre-dispersal seed survival)
#p.vec[38]: delta (probability of dispersal)
#p.vec[39]: tau2 (post-dispersal seed survival)




# Define the kernels and iteration matrix:
#Domain 1: Seedling Domain
#	for diam=diameter in range [0, 1.6]
#	for height in range [0, 16]
#Domain 2: Larger Domain
#	 for diam=diameter in range [1.6,700]
#and for height=height in range [16, 800]
# Domain 2 is split into a 
#============================================================================# 





## Survival-Growth of Seedlings (in D1)
pyx1=function(diamp,heightp,diam, height, params) { (1-((exp(params[12]+params[13]*diam + params[14]*height)/(1+exp(params[12]+params[13]*diam + params[14]*height))
)))*(exp(params[1]+params[2]*diam+params[3]*height)/(1+exp(params[1]+params[2]*diam+params[3]*height)))*dtnorm(diamp,mean=params[4] + params[5]*diam + params[6]*height,sd=params[7], lower=0, upper=1.6)*dtnorm(heightp,mean=params[8]+params[9]*diam + params[10]*height,sd=params[11], lower=0, upper=16)}


# Survival-Growth of Larger Plants (in D2)
#Note that D2 survival is multiplied by 0.997 to avoid 100% survival 

pyx2=function(diamp,heightp,diam, height, params) { 0.997*(exp(params[19]+params[20]*diam+params[21]*height)/(1+exp(params[19]+params[20]*diam+params[21]*height)))*dtnorm(diamp,mean=params[22] + params[23]*diam + params[24]*height,sd=params[25], lower=1.6, upper=700)*dtnorm(heightp,mean=params[26]+params[27]*diam + params[28]*height,sd=params[29], lower=16, upper=800)}


# Fecundity=P(fruiting)*# of Fruits Produced*P(survival of seeds)*Distribution of Seedling Diameters*Distribution of Seedling Heights
fyx=function(diamp,heightp,diam,height, params) {(exp(params[30]+params[31]*height)/(1+exp(params[30]+params[31]*height)))*(params[32]*diam*diam)*params[38]*params[37]*params[39]*dtnorm(diamp,mean=params[33],sd=params[34], lower=0, upper=1.6)*dtnorm(heightp,mean=params[35],sd=params[36], lower=0, upper=34)}

#Maturation = P(maturation)*Distribution of Graduates Diams*Distribution of Graduates Heights
gyx=function(diamp,heightp,diam,height, params) {((exp(params[1]+params[2]*diam+params[3]*height)/(1+exp(params[1]+params[2]*diam+params[3]*height))))*(exp(params[12]+params[13]*diam + params[14]*height)/(1+exp(params[12]+params[13]*diam + params[14]*height))
)*dtnorm(diamp,mean=params[15],sd=params[16], lower=1.6, upper=700)*dtnorm(heightp,mean=params[17],sd=params[18], lower=16, upper=800)}

#kyx=function(diamp,heightp,diam,height, params) {pyx(diamp,heightp, diam,height, params)+fyx(diamp,heightp,diam,height, params)}

# Compute meshpoints  

###Fix these to match domains 
h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1

h3a = (150-1.6)/m3a; 

y3a=(h3a/2)*((0:(m3a-1))+(1:m3a))+1.6; #for diameter in D2a
h3b = (700-150)/m3b
y3b=(h3b/2)*((0:(m3b-1))+(1:m3b))+150; #for diameter in D2b
y3 <- c(y3a, y3b)

#h4=(800-16)/m4
h4a = (300-16)/m4a
y4a=(h4a/2)*((0:(m4a-1))+(1:m4a))+16;
h4b = (800-300)/m4b 
y4b=(h4b/2)*((0:(m4b-1))+(1:m4b))+300;
y4 = c(y4a, y4b)



# Compute the iteration matrix. With a bit of vectorizing it's not too slow,
# though you can probably do better if you need to. The shortcuts here have 
# been checked against the results from code that uses loops for everything. (comment from Ellner and Rees)


# Construct D1 (Seedling Domain): #####

build_D1 = function(p.vec) {
  plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in A 
  Plop=outer(1:m1,1:m2,plop) 
  
  D1=matrix(0,m1*m2,m1*m2)
  Kvals_D1=array(0,c(m1,m2,m1,m2))  
  
  for(i in 1:m1){
    for(j in 1:m2){
      for(k in 1:m1){
        kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec)
        D1[Plop[k,1:m2],Plop[i,j]]=kvals
        Kvals_D1[k,1:m2,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }
  
  D1=D1*h1*h2 #multiply D1 by widths
  return(list(D1 = D1, Kvals_D1 = Kvals_D1))
}

#Big Cypress
thing <- build_D1(p.vec_BC)
saveRDS(thing$D1, file = "./BC/D1_BC.rds")
saveRDS(thing$Kvals_D1, file ="./BC/Kvals_D1_BC.rds" )

#Cape Canaveral
thing <- build_D1(p.vec_CC)
saveRDS(thing$D1, file = "./CC/D1_CC.rds")
saveRDS(thing$Kvals_D1, file ="./CC/Kvals_D1_CC.rds" )

#Chekika
thing <- build_D1(p.vec_C)
saveRDS(thing$D1, file = "./C/D1_C.rds")
saveRDS(thing$Kvals_D1, file ="./C/Kvals_D1_C.rds" )

#Fort Pierce
thing <- build_D1(p.vec_FP)
saveRDS(thing$D1, file = "./FP/D1_FP.rds")
saveRDS(thing$Kvals_D1, file ="./FP/Kvals_D1_FP.rds" )

#Punta Gorda
thing <- build_D1(p.vec_PG)
saveRDS(thing$D1, file = "./PG/D1_PG.rds")
saveRDS(thing$Kvals_D1, file ="./PG/Kvals_D1_PG.rds" )

#Wild Turkey
thing <- build_D1(p.vec_WT)
saveRDS(thing$D1, file = "./WT/D1_WT.rds")
saveRDS(thing$Kvals_D1, file ="./WT/Kvals_D1_WT.rds" )

#Overall
thing <- build_D1(p.vec_overall)
saveRDS(thing$D1, file="./Overall/D1_overall.rds")
saveRDS(thing$Kvals_D1, file="./Overall/Kvals_D1_overall.rds")

# Construct D2AA (Large Domain): #####
build_D2AA = function(p.vec) {
  plop=function(i,j) {(j-1)*m3a+i} # for putting values in proper place in A 
  Plop=outer(1:m3a,1:m4a,plop); 
  
  D2AA=matrix(0,m3a*m4a,m3a*m4a);
  Kvals_D2AA=array(0,c(m3a,m4a,m3a,m4a));  
  
  for(i in 1:m3a){
    for(j in 1:m4a){
      for(k in 1:m3a){
        kvals=pyx2(y3a[k],y4a[1:m4a],y3a[i],y4a[j], p.vec)
        D2AA[Plop[k,1:m4a],Plop[i,j]]=kvals
        Kvals_D2AA[k,1:m4a,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }		
  D2AA=D2AA*h3a*h4a #Multiply D2 by widths
  return(list(D2AA = D2AA, Kvals_D2AA = Kvals_D2AA))
}

#Big Cypress
thing <- build_D2AA(p.vec_BC)
saveRDS(thing$D2AA, file = "./BC/D2AA_BC.rds")
saveRDS(thing$Kvals_D2AA, file ="./BC/Kvals_D2AA_BC.rds" )
#Cape Canaveral
thing <- build_D2AA(p.vec_CC)
saveRDS(thing$D2AA, file = "./CC/D2AA_CC.rds")
saveRDS(thing$Kvals_D2AA, file ="./CC/Kvals_D2AA_CC.rds" )
#Chekika
thing <- build_D2AA(p.vec_C)
saveRDS(thing$D2AA, file = "./C/D2AA_C.rds")
saveRDS(thing$Kvals_D2AA, file ="./C/Kvals_D2AA_C.rds" )
#Fort Pierce
thing <- build_D2AA(p.vec_FP)
saveRDS(thing$D2AA, file = "./FP/D2AA_FP.rds")
saveRDS(thing$Kvals_D2AA, file ="./FP/Kvals_D2AA_FP.rds" )
#Punta Gorda
thing <- build_D2AA(p.vec_PG)
saveRDS(thing$D2AA, file = "./PG/D2AA_PG.rds")
saveRDS(thing$Kvals_D2AA, file ="./PG/Kvals_D2AA_PG.rds" )
#Wild Turkey
thing <- build_D2AA(p.vec_WT)
saveRDS(thing$D2AA, file = "./WT/D2AA_WT.rds")
saveRDS(thing$Kvals_D2AA, file ="./WT/Kvals_D2AA_WT.rds" )

#Overall
thing <- build_D2AA(p.vec_overall)
saveRDS(thing$D2AA, file="./Overall/D2AA_overall.rds")
saveRDS(thing$Kvals_D2AA, file="./Overall/Kvals_D2AA_overall.rds")

# Construct D2BB (Large Domain): #####
build_D2BB = function(p.vec) {
  
  plop=function(i,j) {(j-1)*m3b+i} # for putting values in proper place in A 
  Plop=outer(1:m3b,1:m4b,plop); 
  
  D2BB=matrix(0,m3b*m4b,m3b*m4b);
  Kvals_D2BB=array(0,c(m3b,m4b,m3b,m4b));  
  
  for(i in 1:m3b){
    for(j in 1:m4b){
      for(k in 1:m3b){
        kvals=pyx2(y3b[k],y4b[1:m4b],y3b[i],y4b[j], p.vec)
        D2BB[Plop[k,1:m4b],Plop[i,j]]=kvals
        Kvals_D2BB[k,1:m4b,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }		
  D2BB=D2BB*h3b*h4b #Multiply D2 by widths
  return(list(D2BB = D2BB, Kvals_D2BB = Kvals_D2BB))
}

#Big Cypress
thing <- build_D2BB(p.vec_BC)
saveRDS(thing$D2BB, file = "./BC/D2BB_BC.rds")
saveRDS(thing$Kvals_D2BB, file ="./BC/Kvals_D2BB_BC.rds" )
#Cape Canaveral
thing <- build_D2BB(p.vec_CC)
saveRDS(thing$D2BB, file = "./CC/D2BB_CC.rds")
saveRDS(thing$Kvals_D2BB, file ="./CC/Kvals_D2BB_CC.rds" )
#Chekika
thing <- build_D2BB(p.vec_C)
saveRDS(thing$D2BB, file = "./C/D2BB_C.rds")
saveRDS(thing$Kvals_D2BB, file ="./C/Kvals_D2BB_C.rds" )
#Fort Pierce
thing <- build_D2BB(p.vec_FP)
saveRDS(thing$D2BB, file = "./FP/D2BB_FP.rds")
saveRDS(thing$Kvals_D2BB, file ="./FP/Kvals_D2BB_FP.rds" )
#Punta Gorda
thing <- build_D2BB(p.vec_PG)
saveRDS(thing$D2BB, file = "./PG/D2BB_PG.rds")
saveRDS(thing$Kvals_D2BB, file ="./PG/Kvals_D2BB_PG.rds" )
#Wild Turkey
thing <- build_D2BB(p.vec_WT)
saveRDS(thing$D2BB, file = "./WT/D2BB_WT.rds")
saveRDS(thing$Kvals_D2BB, file ="./WT/Kvals_D2BB_WT.rds" )

#Overall
thing <- build_D2BB(p.vec_overall)
saveRDS(thing$D2BB, file="./Overall/D2BB_overall.rds")
saveRDS(thing$Kvals_D2BB, file="./Overall/Kvals_D2BB_overall.rds")

# Construct D2BA (Large Domain): #####
build_D2BA = function(p.vec) {
  
  plop1=function(i, j) {(j-1)*m3a + i}
  plop2=function(i, j) {(j-1)*m3b + i}
  Plop1=outer(1:m3a,1:m4a,plop1); 
  Plop2=outer(1:m3b, 1:m4b, plop2);
  
  D2BA=matrix(0,m3a*m4a,m3b*m4b); 
  Kvals_D2BA=array(0, c(m3a, m4a, m3b, m4b))
  
  for(i in 1:m3b) {
    for (j in 1:m4b) {
      for (k in 1:m3a) {
        kvals=pyx2(y3a[k], y4a[1:m4a], y3b[i], y4b[j], p.vec)
        D2BA[Plop1[k, 1:m4a], Plop2[i,j]]=kvals
        Kvals_D2BA[k, 1:m4a, i, j]=kvals
      }}
    cat(i, "\n");
  }
  D2BA=D2BA*h3a*h4a
  return(list(D2BA = D2BA, Kvals_D2BA = Kvals_D2BA))
}

#Big Cypress
thing <- build_D2BA(p.vec_BC)
saveRDS(thing$D2BA, file = "./BC/D2BA_BC.rds")
saveRDS(thing$Kvals_D2BA, file ="./BC/Kvals_D2BA_BC.rds" )
#Cape Canaveral
thing <- build_D2BA(p.vec_CC)
saveRDS(thing$D2BA, file = "./CC/D2BA_CC.rds")
saveRDS(thing$Kvals_D2BA, file ="./CC/Kvals_D2BA_CC.rds" )
#Chekika
thing <- build_D2BA(p.vec_C)
saveRDS(thing$D2BA, file = "./C/D2BA_C.rds")
saveRDS(thing$Kvals_D2BA, file ="./C/Kvals_D2BA_C.rds" )
#Fort Pierce
thing <- build_D2BA(p.vec_FP)
saveRDS(thing$D2BA, file = "./FP/D2BA_FP.rds")
saveRDS(thing$Kvals_D2BA, file ="./FP/Kvals_D2BA_FP.rds" )
#Punta Gorda
thing <- build_D2BA(p.vec_PG)
saveRDS(thing$D2BA, file = "./PG/D2BA_PG.rds")
saveRDS(thing$Kvals_D2BA, file ="./PG/Kvals_D2BA_PG.rds" )
#Wild Turkey
thing <- build_D2BA(p.vec_WT)
saveRDS(thing$D2BA, file = "./WT/D2BA_WT.rds")
saveRDS(thing$Kvals_D2BA, file ="./WT/Kvals_D2BA_WT.rds" )

#Overall
thing <- build_D2BA(p.vec_overall)
saveRDS(thing$D2BA, file="./Overall/D2BA_overall.rds")
saveRDS(thing$Kvals_D2BA, file="./Overall/Kvals_D2BA_overall.rds")

# Construct D2AB (Large Domain): #####
build_D2AB = function(p.vec) {
  
  plop1=function(i, j) {(j-1)*m3b + i}
  plop2=function(i, j) {(j-1)*m3a + i}
  Plop1=outer(1:m3b,1:m4b,plop1); 
  Plop2=outer(1:m3a, 1:m4a, plop2);
  
  D2AB=matrix(0,m3b*m4b,m3a*m4a); 
  Kvals_D2AB=array(0, c(m3b, m4b, m3a, m4a))
  
  for(i in 1:m3a) {
    for (j in 1:m4a) {
      for (k in 1:m3b) {
        kvals=pyx2(y3b[k], y4b[1:m4b], y3a[i], y4a[j], p.vec)
        D2AB[Plop1[k, 1:m4b], Plop2[i,j]]=kvals
        Kvals_D2AB[k, 1:m4b, i, j]=kvals
      }}
    cat(i, "\n");
  }
  D2AB=D2AB*h3b*h4b
  return(list(D2AB = D2AB, Kvals_D2AB = Kvals_D2AB))
}

#Big Cypress
thing <- build_D2AB(p.vec_BC)
saveRDS(thing$D2AB, file = "./BC/D2AB_BC.rds")
saveRDS(thing$Kvals_D2AB, file ="./BC/Kvals_D2AB_BC.rds" )
#Cape Canaveral
thing <- build_D2AB(p.vec_CC)
saveRDS(thing$D2AB, file = "./CC/D2AB_CC.rds")
saveRDS(thing$Kvals_D2AB, file ="./CC/Kvals_D2AB_CC.rds" )
#Chekika
thing <- build_D2AB(p.vec_C)
saveRDS(thing$D2AB, file = "./C/D2AB_C.rds")
saveRDS(thing$Kvals_D2AB, file ="./C/Kvals_D2AB_C.rds" )
#Fort Pierce
thing <- build_D2AB(p.vec_FP)
saveRDS(thing$D2AB, file = "./FP/D2AB_FP.rds")
saveRDS(thing$Kvals_D2AB, file ="./FP/Kvals_D2AB_FP.rds" )
#Punta Gorda
thing <- build_D2AB(p.vec_PG)
saveRDS(thing$D2AB, file = "./PG/D2AB_PG.rds")
saveRDS(thing$Kvals_D2AB, file ="./PG/Kvals_D2AB_PG.rds" )
#Wild Turkey
thing <- build_D2AB(p.vec_WT)
saveRDS(thing$D2AB, file = "./WT/D2AB_WT.rds")
saveRDS(thing$Kvals_D2AB, file ="./WT/Kvals_D2AB_WT.rds" )

#Overall
thing <- build_D2AB(p.vec_overall)
saveRDS(thing$D2AB, file="./Overall/D2AB_overall.rds")
saveRDS(thing$Kvals_D2AB, file="./Overall/Kvals_D2AB_overall.rds")

# Put D2 together #####
#Big Cypress
D2AA_BC <- readRDS("./BC/D2AA_BC.rds")
D2BA_BC <- readRDS("./BC/D2BA_BC.rds")
D2AB_BC <- readRDS("./BC/D2AB_BC.rds")
D2BB_BC <- readRDS("./BC/D2BB_BC.rds")

D2_BC <- rbind(cbind(D2AA_BC, D2BA_BC), cbind(D2AB_BC, D2BB_BC))   
rm(D2AA_BC, D2BA_BC, D2AB_BC, D2BB_BC)
save(D2_BC, file="./BC/D2_BC.RData")

D2AA_CC <- readRDS("./CC/D2AA_CC.rds")
D2BA_CC <- readRDS("./CC/D2BA_CC.rds")
D2AB_CC <- readRDS("./CC/D2AB_CC.rds")
D2BB_CC <- readRDS("./CC/D2BB_CC.rds")

D2_CC <- rbind(cbind(D2AA_CC, D2BA_CC), cbind(D2AB_CC, D2BB_CC))   
rm(D2AA_CC, D2BA_CC, D2AB_CC, D2BB_CC)
save(D2_CC, file="./CC/D2_CC.RData")
rm(D2_CC)

#Chekika
D2AA_C <- readRDS("./C/D2AA_C.rds")
D2BA_C <- readRDS("./C/D2BA_C.rds")
D2AB_C <- readRDS("./C/D2AB_C.rds")
D2BB_C <- readRDS("./C/D2BB_C.rds")

D2_C <- rbind(cbind(D2AA_C, D2BA_C), cbind(D2AB_C, D2BB_C))   
rm(D2AA_C, D2BA_C, D2AB_C, D2BB_C)
save(D2_C, file="./C/D2_C.RData")
rm(D2_C)

#Fort Pierce
D2AA_FP <- readRDS("./FP/D2AA_FP.rds")
D2BA_FP <- readRDS("./FP/D2BA_FP.rds")
D2AB_FP <- readRDS("./FP/D2AB_FP.rds")
D2BB_FP <- readRDS("./FP/D2BB_FP.rds")

D2_FP <- rbind(cbind(D2AA_FP, D2BA_FP), cbind(D2AB_FP, D2BB_FP))   
rm(D2AA_FP, D2BA_FP, D2AB_FP, D2BB_FP)
save(D2_FP, file="./FP/D2_FP.RData")
rm(D2_FP)

#Punta Gorda
D2AA_PG <- readRDS("./PG/D2AA_PG.rds")
D2BA_PG <- readRDS("./PG/D2BA_PG.rds")
D2AB_PG <- readRDS("./PG/D2AB_PG.rds")
D2BB_PG <- readRDS("./PG/D2BB_PG.rds")

D2_PG <- rbind(cbind(D2AA_PG, D2BA_PG), cbind(D2AB_PG, D2BB_PG))   
rm(D2AA_PG, D2BA_PG, D2AB_PG, D2BB_PG)
save(D2_PG, file="./PG/D2_PG.RData")
rm(D2_PG)
#Wild Turkey
D2AA_WT <- readRDS("./WT/D2AA_WT.rds")
D2BA_WT <- readRDS("./WT/D2BA_WT.rds")
D2AB_WT <- readRDS("./WT/D2AB_WT.rds")
D2BB_WT <- readRDS("./WT/D2BB_WT.rds")

D2_WT <- rbind(cbind(D2AA_WT, D2BA_WT), cbind(D2AB_WT, D2BB_WT))   
rm(D2AA_WT, D2BA_WT, D2AB_WT, D2BB_WT)
save(D2_WT, file="./WT/D2_WT.RData")
rm(D2_WT)
#Overall
D2AA_overall <- readRDS("./Overall/D2AA_overall.rds")
D2BA_overall <- readRDS("./Overall/D2BA_overall.rds")
D2AB_overall <- readRDS("./Overall/D2AB_overall.rds")
D2BB_overall <- readRDS("./Overall/D2BB_overall.rds")

D2_overall <- rbind(cbind(D2AA_overall, D2BA_overall), cbind(D2AB_overall, D2BB_overall))   
rm(D2AA_overall, D2BA_overall, D2AB_overall, D2BB_overall)
save(D2_overall, file="./Overall/D2_overall.RData")
rm(D2_overall)


# Construct FA (Fertility): #####
build_FA = function(p.vec) {
  plop1=function(i, j) {(j-1)*m1 + i}
  plop2=function(i, j) {(j-1)*m3a + i}
  Plop1=outer(1:m1,1:m2,plop1); 
  Plop2=outer(1:m3a, 1:m4a, plop2);
  
  FA=matrix(0,m1*m2,m3a*m4a); 
  Kvals_FA=array(0, c(m1, m2, m3a, m4a))
  
  for(i in 1:m3a) {
    for (j in 1:m4a) {
      for (k in 1:m1) {
        kvals=fyx(y1[k], y2[1:m2], y3a[i], y4a[j], p.vec)
        FA[Plop1[k, 1:m2], Plop2[i,j]]=kvals
        Kvals_FA[k, 1:m2, i, j]=kvals
      }}
    cat(i, "\n");
  }
  FA=FA*h1*h2
  return(list(FA = FA, Kvals_FA = Kvals_FA))
}

# Construct FB (Fertility): #####
build_FB = function(p.vec) {
  plop1=function(i, j) {(j-1)*m1 + i}
  plop2=function(i, j) {(j-1)*m3b + i}
  Plop1=outer(1:m1,1:m2,plop1); 
  Plop2=outer(1:m3b, 1:m4b, plop2);
  
  FB=matrix(0,m1*m2,m3b*m4b); 
  Kvals_FB=array(0, c(m1, m2, m3b, m4b))
  
  for(i in 1:m3b) {
    for (j in 1:m4b) {
      for (k in 1:m1) {
        kvals=fyx(y1[k], y2[1:m2], y3b[i], y4b[j], p.vec)
        FB[Plop1[k, 1:m2], Plop2[i,j]]=kvals
        Kvals_FB[k, 1:m2, i, j]=kvals
      }}
    cat(i, "\n");
  }
  FB=FB*h1*h2
  return(list(FB = FB, Kvals_FB = Kvals_FB))
}

#Big Cypress
thing <- build_FA(p.vec_BC)
saveRDS(thing$FA, file = "./BC/FA_BC.rds")
saveRDS(thing$Kvals_FA, file ="./BC/Kvals_FA_BC.rds" )
FA_BC <- thing$FA
thing <- build_FB(p.vec_BC)
saveRDS(thing$FB, file = "./BC/FB_BC.rds")
saveRDS(thing$Kvals_FB, file ="./BC/Kvals_FB_BC.rds" )
FB_BC <- thing$FB
F_BC <- cbind(FA_BC, FB_BC)
save(F_BC, file="./BC/F_BC.RData")
rm (F_BC, FA_BC, FB_BC)



#Cape Canaveral
thing <- build_FA(p.vec_CC)
saveRDS(thing$FA, file = "./CC/FA_CC.rds")
saveRDS(thing$Kvals_FA, file ="./CC/Kvals_FA_CC.rds" )
FA_CC <- thing$FA
thing <- build_FB(p.vec_CC)
saveRDS(thing$FB, file = "./CC/FB_CC.rds")
saveRDS(thing$Kvals_FB, file ="./CC/Kvals_FB_CC.rds" )
FB_CC <- thing$FB
F_CC <- cbind(FA_CC, FB_CC)
save(F_CC, file="./CC/F_CC.RData")
rm (F_CC, FA_CC, FB_CC)

#Chekika
thing <- build_FA(p.vec_C)
saveRDS(thing$FA, file = "./C/FA_C.rds")
saveRDS(thing$Kvals_FA, file ="./C/Kvals_FA_C.rds" )
FA_C <- thing$FA
thing <- build_FB(p.vec_C)
saveRDS(thing$FB, file = "./C/FB_C.rds")
saveRDS(thing$Kvals_FB, file ="./C/Kvals_FB_C.rds" )
FB_C <- thing$FB
F_C <- cbind(FA_C, FB_C)
save(F_C, file="./C/F_C.RData")
rm (F_C, FA_C, FB_C)


#Fort Pierce
thing <- build_FA(p.vec_FP)
saveRDS(thing$FA, file = "./FP/FA_FP.rds")
saveRDS(thing$Kvals_FA, file ="./FP/Kvals_FA_FP.rds" )
FA_FP <- thing$FA
thing <- build_FB(p.vec_FP)
saveRDS(thing$FB, file = "./FP/FB_FP.rds")
saveRDS(thing$Kvals_FB, file ="./FP/Kvals_FB_FP.rds" )
FB_FP <- thing$FB
F_FP <- cbind(FA_FP, FB_FP)
save(F_FP, file="./FP/F_FP.RData")
rm (F_FP, FA_FP, FB_FP)

#Punta Gorda
thing <- build_FA(p.vec_PG)
saveRDS(thing$FA, file = "./PG/FA_PG.rds")
saveRDS(thing$Kvals_FA, file ="./PG/Kvals_FA_PG.rds" )
FA_PG <- thing$FA
thing <- build_FB(p.vec_PG)
saveRDS(thing$FB, file = "./PG/FB_PG.rds")
saveRDS(thing$Kvals_FB, file ="./PG/Kvals_FB_PG.rds" )
FB_PG <- thing$FB
F_PG <- cbind(FA_PG, FB_PG)
save(F_PG, file="./PG/F_PG.RData")
rm (F_PG, FA_PG, FB_PG)


#Wild Turkey
thing <- build_FA(p.vec_WT)
saveRDS(thing$FA, file = "./WT/FA_WT.rds")
saveRDS(thing$Kvals_FA, file ="./WT/Kvals_FA_WT.rds" )
FA_WT <- thing$FA
thing <- build_FB(p.vec_WT)
saveRDS(thing$FB, file = "./WT/FB_WT.rds")
saveRDS(thing$Kvals_FB, file ="./WT/Kvals_FB_WT.rds" )
FB_WT <- thing$FB
F_WT <- cbind(FA_WT, FB_WT)
save(F_WT, file="./WT/F_WT.RData")
rm (F_WT, FA_WT, FB_WT)


#Overall
thing <- build_FA(p.vec_overall)
saveRDS(thing$FA, file = "./Overall/FA_overall.rds")
saveRDS(thing$Kvals_FA, file ="./Overall/Kvals_FA_overall.rds" )
FA_overall <- thing$FA
thing <- build_FB(p.vec_overall)
saveRDS(thing$FB, file = "./Overall/FB_overall.rds")
saveRDS(thing$Kvals_FB, file ="./Overall/Kvals_FB_overall.rds" )
FB_overall <- thing$FB
F_overall <- cbind(FA_overall, FB_overall)
save(F_overall, file="./Overall/F_overall.RData")
rm (F_overall, FA_overall, FB_overall)

# Construct MA (Maturation): #####

build_GA = function(p.vec) {
  plop1=function(i, j) {(j-1)*m3a + i}
  plop2=function(i, j) {(j-1)*m1 + i}
  Plop1=outer(1:m3a,1:m4a,plop1); 
  Plop2=outer(1:m1, 1:m2, plop2);
  
  GA=matrix(0,m3a*m4a,m1*m2); 
  Kvals_GA=array(0, c(m3a, m4a, m1, m2))
  
  for(i in 1:m1) {
    for (j in 1:m2) {
      for (k in 1:m3a) {
        kvals=gyx(y3a[k], y4a[1:m4a], y1[i], y2[j], p.vec)
        GA[Plop1[k, 1:m4a], Plop2[i,j]]=kvals
        Kvals_GA[k, 1:m4a, i, j]=kvals
      }}
    cat(i, "\n");
  }
  GA=GA*h3a*h4a
  
  return(list(GA = GA, Kvals_GA = Kvals_GA))
}

build_GB = function(p.vec) {
  plop1=function(i, j) {(j-1)*m3b + i}
  plop2=function(i, j) {(j-1)*m1 + i}
  Plop1=outer(1:m3b,1:m4b,plop1); 
  Plop2=outer(1:m1, 1:m2, plop2);
  
  GB=matrix(0,m3b*m4b,m1*m2); 
  Kvals_GB=array(0, c(m3b, m4b, m1, m2))
  
  for(i in 1:m1) {
    for (j in 1:m2) {
      for (k in 1:m3b) {
        kvals=gyx(y3b[k], y4b[1:m4b], y1[i], y2[j], p.vec)
        GB[Plop1[k, 1:m4b], Plop2[i,j]]=kvals
        Kvals_GB[k, 1:m4b, i, j]=kvals
      }}
    cat(i, "\n");
  }
  GB=GB*h3b*h4b
  
  return(list(GB = GB, Kvals_GB = Kvals_GB))
}

#Big Cypress
thing <- build_GA(p.vec_BC)
saveRDS(thing$GA, file = "./BC/GA_BC.rds")
saveRDS(thing$Kvals_GA, file ="./BC/Kvals_GA_BC.rds" )
GA<-thing$GA

thing <- build_GB(p.vec_BC)
saveRDS(thing$GB, file="./BC/GB_BC.rds")
saveRDS(thing$Kvals_GB, file="./BC/Kvals_GB_BC.rds")
GB<-thing$GB


G_BC <- rbind(GA, GB)
save(G_BC, file="./BC/G_BC.RData")
rm(GA, GB, G_BC)

#Cape Canaveral
thing <- build_GA(p.vec_CC)
saveRDS(thing$GA, file = "./CC/GA_CC.rds")
saveRDS(thing$Kvals_GA, file ="./CC/Kvals_GA_CC.rds" )
GA_CC <- thing$GA

thing <- build_GB(p.vec_CC)
saveRDS(thing$GB, file="./CC/GB_CC.rds")
saveRDS(thing$Kvals_GB, file="./CC/Kvals_GB_CC.rds")
GB_CC <- thing$GB
G_CC <- rbind(GA_CC, GB_CC)
save(G_CC, file="./CC/G_CC.RData")
rm(GA_CC, GB_CC, G_CC)

#Chekika
thing <- build_GA(p.vec_C)
saveRDS(thing$GA, file = "./C/GA_C.rds")
saveRDS(thing$Kvals_GA, file ="./C/Kvals_GA_C.rds" )
GA_C <- thing$GA
thing <- build_GB(p.vec_C)
saveRDS(thing$GB, file="./C/GB_C.rds")
saveRDS(thing$Kvals_GB, file="./C/Kvals_GB_C.rds")
GB_C <- thing$GB
G_C <- rbind(GA_C, GB_C)
save(G_C, file="./C/G_C.RData")
rm(GA_C, GB_C, G_C)

#Fort Pierce
thing <- build_GA(p.vec_FP)
saveRDS(thing$GA, file = "./FP/GA_FP.rds")
saveRDS(thing$Kvals_GA, file ="./FP/Kvals_GA_FP.rds" )
GA_FP <- thing$GA

thing <- build_GB(p.vec_FP)
saveRDS(thing$GB, file="./FP/GB_FP.rds")
saveRDS(thing$Kvals_GB, file="./FP/Kvals_GB_FP.rds")
GB_FP <- thing$GB
G_FP <- rbind(GA_FP, GB_FP)
save(G_FP, file="./FP/G_FP.RData")
rm(GA_FP, GB_FP, G_FP)

#Punta Gorda
thing <- build_GA(p.vec_PG)
saveRDS(thing$GA, file = "./PG/GA_PG.rds")
saveRDS(thing$Kvals_GA, file ="./PG/Kvals_GA_PG.rds" )
GA_PG <- thing$GA

thing <- build_GB(p.vec_PG)
saveRDS(thing$GB, file="./PG/GB_PG.rds")
saveRDS(thing$Kvals_GB, file="./PG/Kvals_GB_PG.rds")
GB_PG <- thing$GB
G_PG <- rbind(GA_PG, GB_PG)
save(G_PG, file="./PG/G_PG.RData")
rm(GA_PG, GB_PG, G_PG)

#Wild Turkey
thing <- build_GA(p.vec_WT)
saveRDS(thing$GA, file = "./WT/GA_WT.rds")
saveRDS(thing$Kvals_GA, file ="./WT/Kvals_GA_WT.rds" )
GA_WT <- thing$GA

thing <- build_GB(p.vec_WT)
saveRDS(thing$GB, file="./WT/GB_WT.rds")
saveRDS(thing$Kvals_GB, file="./WT/Kvals_GB_WT.rds")
GB_WT <- thing$GB
G_WT <- rbind(GA_WT, GB_WT)
save(G_WT, file="./WT/G_WT.RData")
rm(GA_WT, GB_WT, G_WT)

#Overall
thing <- build_GA(p.vec_overall)
saveRDS(thing$GA, file = "./Overall/GA_overall.rds")
saveRDS(thing$Kvals_GA, file ="./Overall/Kvals_GA_overall.rds" )
GA_overall <- thing$GA

thing <- build_GB(p.vec_overall)
saveRDS(thing$GB, file="./Overall/GB_overall.rds")
saveRDS(thing$Kvals_GB, file="./Overall/Kvals_GB_overall.rds")
GB_overall <-thing$GB
G_overall <- rbind(GA_overall, GB_overall)
save(G_overall, file="./Overall/G_overall.RData")
rm(GA_overall, GB_overall, G_overall)

rm(thing)
gc()


# Assemble the matrices #####

#Big Cypress
D1_BC <- readRDS("./BC/D1_BC.rds")
load("./BC/G_BC.RData")
left_side<-rbind(D1_BC, G_BC)
rm(D1_BC, G_BC)
load("./BC/F_BC.RData")
load("./BC/D2_BC.RData")
right_side<-rbind(F_BC, D2_BC)
rm(F_BC, D2_BC)
gc()
A_BC<-cbind(left_side, right_side)
save(A_BC, file="./BC/A_BC.RData")
rm(A_BC)

#Cape Canaveral
D1_CC <-readRDS("./CC/D1_CC.rds")
load("./CC/G_CC.RData")
left_side<-rbind(D1_CC, G_CC)
rm(D1_CC, G_CC)
load("./CC/F_CC.RData")
load("./CC/D2_CC.RData")
right_side<-rbind(F_CC, D2_CC)
rm(F_CC, D2_CC)
A_CC<-cbind(left_side, right_side)
save(A_CC, file="./CC/A_CC.RData")
rm(A_CC)

#Chekika
D1_C <-readRDS("./C/D1_C.rds")
load("./C/G_C.RData")
left_side<-rbind(D1_C, G_C)
rm(D1_C, G_C)
load("./C/F_C.RData")
load("./C/D2_C.RData")
right_side<-rbind(F_C, D2_C)
rm(F_C, D2_C)
A_C<-cbind(left_side, right_side)
save(A_C, file="./C/A_C.RData")
rm(A_C)

#Fort Pierce
D1_FP <-readRDS("./FP/D1_FP.rds")
load("./FP/G_FP.RData")
left_side<-rbind(D1_FP, G_FP)
rm(D1_FP, G_FP)
load("./FP/F_FP.RData")
load("./FP/D2_FP.RData")
right_side<-rbind(F_FP, D2_FP)
rm(F_FP, D2_FP)
A_FP<-cbind(left_side, right_side)
save(A_FP, file="./FP/A_FP.RData")
rm(A_FP)

#Punta Gorda
D1_PG <-readRDS("./PG/D1_PG.rds")
load("./PG/G_PG.RData")
left_side<-rbind(D1_PG, G_PG)
rm(D1_PG, G_PG)
load("./PG/F_PG.RData")
load("./PG/D2_PG.RData")
right_side<-rbind(F_PG, D2_PG)
rm(F_PG, D2_PG)
A_PG<-cbind(left_side, right_side)
save(A_PG, file="./PG/A_PG.RData")
rm(A_PG)

#Wild Turkey
D1_WT <-readRDS("./WT/D1_WT.rds")
load("./WT/G_WT.RData")
left_side<-rbind(D1_WT, G_WT)
rm(D1_WT, G_WT)
load("./WT/F_WT.RData")
load("./WT/D2_WT.RData")
right_side<-rbind(F_WT, D2_WT)
rm(F_WT, D2_WT)
A_WT<-cbind(left_side, right_side)
save(A_WT, file="./WT/A_WT.RData")
rm(A_WT)

#Overall
D1_overall <-readRDS("./Overall/D1_overall.rds")
load("./overall/G_overall.RData")
left_side<-rbind(D1_overall, G_overall)
rm(D1_overall, G_overall)
load("./Overall/F_overall.RData")
load("./Overall/D2_overall.RData")
right_side<-rbind(F_overall, D2_overall)
rm(F_overall, D2_overall)
A_overall<-cbind(left_side, right_side)
save(A_overall, file="./Overall/A_overall.RData")
rm(A_overall)
rm(left_side, right_side)



#  Find lambda, w by iteration #####

#  Note: the Matrix package is used to iterate more quickly. The check  
#  for convergence requires extracting matrix entries via the @x slot
#  of a Matrix object. Matrix is S4-style -- see ?Matrix. 

find_lambda = function(A) {
  
  A2=Matrix(A); nt=Matrix(1,m1*m2+m3a*m4a + m3b*m4b,1); nt1=nt; 
  
  qmax=1000; lam=1; 
  while(qmax>tol) {
    nt1=A2%*%nt;
    qmax=sum(abs((nt1-lam*nt)@x));  
    lam=sum(nt1@x); 
    nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
    cat(lam,qmax,"\n");
  } 
  
  nt=matrix(nt@x,m1*m2+m3a*m4a + m3b*m4b,1); 
  #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
  stable.dist=nt
  lam.stable=lam;
  
  # Check that the @bits worked as intended.   
  qmax=sum(abs(lam*nt-A%*%nt)); 
  cat("Convergence: ",qmax," should be less than ",tol,"\n");
  
  #Find the reproductive value function by iteration
  vt=Matrix(1,1,m1*m2+m3a*m4a + m3b*m4b); vt1=vt; 
  
  qmax=1000; lam=1; 
  while(qmax>tol) {
    vt1=vt%*%A2;
    qmax=sum(abs((vt1-lam*vt)@x));  
    lam=sum(vt1@x); 
    vt@x=(vt1@x)/lam;   
    cat(lam,qmax,"\n");
  } 
  v=t(matrix(vt@x,1,m1*m2+m3a*m4a + m3b*m4b)); 
  lam.stable.t=lam; 
  
  return(list(lam.stable = lam.stable, stable.dist = stable.dist, v=v))
}

#Big Cypress
load("./BC/A_BC.RData")
thing <- find_lambda(A_BC)
saveRDS(thing$lam.stable, file = "./BC/lam.stable_BC.rds")
saveRDS(thing$stable.dist, file = "./BC/stable.dist_BC.rds")
saveRDS(thing$v, file = "./BC/v_BC.rds")
rm(A_BC)

#Cape Canaveral
load("./CC/A_CC.RData")
thing <- find_lambda(A_CC)
saveRDS(thing$lam.stable, file = "./CC/lam.stable_CC.rds")
saveRDS(thing$stable.dist, file = "./CC/stable.dist_CC.rds")
saveRDS(thing$v, file = "./CC/v_CC.rds")
rm(A_CC)

#Chekika
load("./C/A_C.RData")
thing <- find_lambda(A_C)
saveRDS(thing$lam.stable, file = "./C/lam.stable_C.rds")
saveRDS(thing$stable.dist, file = "./C/stable.dist_C.rds")
saveRDS(thing$v, file = "./C/v_C.rds")
rm(A_C)

#Fort Pierce
load("./FP/A_FP.RData")
thing <- find_lambda(A_FP)
saveRDS(thing$lam.stable, file = "./FP/lam.stable_FP.rds")
saveRDS(thing$stable.dist, file = "./FP/stable.dist_FP.rds")
saveRDS(thing$v, file = "./FP/v_FP.rds")
rm(A_FP)

#Punta Gorda
load("./PG/A_PG.RData")
thing <- find_lambda(A_PG)
saveRDS(thing$lam.stable, file = "./PG/lam.stable_PG.rds")
saveRDS(thing$stable.dist, file = "./PG/stable.dist_PG.rds")
saveRDS(thing$v, file = "./PG/v_PG.rds")
rm(A_PG)

#Wild Turkey
load("./WT/A_WT.RData")
thing <- find_lambda(A_WT)
saveRDS(thing$lam.stable, file = "./WT/lam.stable_WT.rds")
saveRDS(thing$stable.dist, file = "./WT/stable.dist_WT.rds")
saveRDS(thing$v, file = "./WT/v_WT.rds")
rm(A_WT)

#Overall
load("./Overall/A_overall.RData")
thing <- find_lambda(A_overall)
saveRDS(thing$lam.stable, file = "./Overall/lam.stable_overall.rds")
saveRDS(thing$stable.dist, file = "./Overall/stable.dist_overall.rds")
saveRDS(thing$v, file = "./Overall/v_overall.rds")
rm(A_overall)





# Elasticity and Sensitivity #####




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


#Method: Compute elasticity matrix first then separate:

elasticity = function(v, stable.dist, A, lam.stable) {
  v.dot.w<-sum(t(v)%*%stable.dist)
  norm_v<-v/v.dot.w
  #check<-t(norm_v)%*%stable.dist #should be 1
  
  rv<-norm_v
  
  sens<-norm_v%*%t(stable.dist)
  elas<-sens*A/lam.stable
  return(list(sens=sens, elas=elas))
} 

#Big Cypress
v_BC <- readRDS("./BC/v_BC.rds")
stable.dist_BC <- readRDS("./BC/stable.dist_BC.rds")
load("./BC/A_BC.RData")
lam.stable_BC <- readRDS("./BC/lam.stable_BC.rds")
thing <- elasticity(v_BC, stable.dist_BC, A_BC, lam.stable_BC)
rm(v_BC, stable.dist_BC, A_BC, lam.stable_BC)
saveRDS(thing$sens, file="./BC/sens_BC.rds")
saveRDS(thing$elas, file="./BC/elas_BC.rds")
rm(thing)

#Cape Canaveral
v_CC <- readRDS("./CC/v_CC.rds")
stable.dist_CC <- readRDS("./CC/stable.dist_CC.rds")
load("./CC/A_CC.RData")
lam.stable_CC <- readRDS("./CC/lam.stable_CC.rds")
thing <- elasticity(v_CC, stable.dist_CC, A_CC, lam.stable_CC)
rm(v_CC, stable.dist_CC, A_CC, lam.stable_CC)
saveRDS(thing$sens, file="./CC/sens_CC.rds")
saveRDS(thing$elas, file="./CC/elas_CC.rds")
rm(thing)
gc()

#Chekika
v_C <- readRDS("./C/v_C.rds")
stable.dist_C <- readRDS("./C/stable.dist_C.rds")
load("./C/A_C.RData")
lam.stable_C <- readRDS("./C/lam.stable_C.rds")
thing <- elasticity(v_C, stable.dist_C, A_C, lam.stable_C)
rm(v_C, stable.dist_C, A_C, lam.stable_C)
saveRDS(thing$sens, file="./C/sens_C.rds")
saveRDS(thing$elas, file="./C/elas_C.rds")
rm(thing)
gc()

#Fort Pierce
v_FP <- readRDS("./FP/v_FP.rds")
stable.dist_FP <- readRDS("./FP/stable.dist_FP.rds")
load("./FP/A_FP.RData")
lam.stable_FP <- readRDS("./FP/lam.stable_FP.rds")
thing <- elasticity(v_FP, stable.dist_FP, A_FP, lam.stable_FP)
rm(v_FP, stable.dist_FP, A_FP, lam.stable_FP)
saveRDS(thing$sens, file="./FP/sens_FP.rds")
saveRDS(thing$elas, file="./FP/elas_FP.rds")
rm(thing)
gc()

#Punta Gorda
v_PG <- readRDS("./PG/v_PG.rds")
stable.dist_PG <- readRDS("./PG/stable.dist_PG.rds")
load("./PG/A_PG.RData")
lam.stable_PG <- readRDS("./PG/lam.stable_PG.rds")
thing <- elasticity(v_PG, stable.dist_PG, A_PG, lam.stable_PG)
rm(v_PG, stable.dist_PG, A_PG, lam.stable_PG)
saveRDS(thing$sens, file="./PG/sens_PG.rds")
saveRDS(thing$elas, file="./PG/elas_PG.rds")
rm(thing)
gc()

#Wild Turkey
v_WT <- readRDS("./WT/v_WT.rds")
stable.dist_WT <- readRDS("./WT/stable.dist_WT.rds")
load("./WT/A_WT.RData")
lam.stable_WT <- readRDS("./WT/lam.stable_WT.rds")
thing <- elasticity(v_WT, stable.dist_WT, A_WT, lam.stable_WT)
rm(v_WT, stable.dist_WT, A_WT, lam.stable_WT)
saveRDS(thing$sens, file="./WT/sens_WT.rds")
saveRDS(thing$elas, file="./WT/elas_WT.rds")
rm(thing)
gc()

#Overall
v_overall <- readRDS("./Overall/v_overall.rds")
stable.dist_overall <- readRDS("./Overall/stable.dist_overall.rds")
load("./Overall/A_overall.RData")
lam.stable_overall <- readRDS("./Overall/lam.stable_overall.rds")
thing <- elasticity(v_overall, stable.dist_overall, A_overall, lam.stable_overall)
rm(v_overall, stable.dist_overall, A_overall, lam.stable_overall)
saveRDS(thing$sens, file="./Overall/sens_overall.rds")
saveRDS(thing$elas, file="./Overall/elas_overall.rds")
rm(thing)
gc()



decompose = function (mat) {
  #Break the elasticity matrix back into its component parts: 
  mat_D1<-mat[1:(m1*m2), 1:(m1*m2)]
  mat_G<-mat[((m1*m2)+1):(m1*m2+m3a*m4a + m3b*m4b), 1:(m1*m2)]
  mat_F <- mat[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3a*m4a + m3b*m4b)]
  mat_D2<-mat[((m1*m2)+1):(m1*m2+m3a*m4a + m3b*m4b), ((m1*m2)+1):(m1*m2 + m3a*m4a + m3b*m4b)]
  return(list(mat_D1 = mat_D1, mat_G = mat_G, mat_F = mat_F, mat_D2 = mat_D2))
}

#Big Cypress
elas_BC <- readRDS("./BC/elas_BC.rds")
thing <- decompose(elas_BC)
saveRDS(thing$mat_D1, file="./BC/elas_D1_BC.rds")
saveRDS(thing$mat_G, file="./BC/elas_G_BC.rds")
saveRDS(thing$mat_F, file="./BC/elas_F_BC.rds")
saveRDS(thing$mat_D2, file="./BC/elas_D2_BC.rds")
rm(elas_BC)
sens_BC <- readRDS("./BC/sens_BC.rds")
thing <- decompose(sens_BC)
saveRDS(thing$mat_D1, file="./BC/sens_D1_BC.rds")
saveRDS(thing$mat_G, file="./BC/sens_G_BC.rds")
saveRDS(thing$mat_F, file="./BC/sens_F_BC.rds")
saveRDS(thing$mat_D2, file="./BC/sens_D2_BC.rds")
rm(sens_BC)

#Cape Canaveral
elas_CC <- readRDS("./CC/elas_CC.rds")
thing <- decompose(elas_CC)
saveRDS(thing$mat_D1, file="./CC/elas_D1_CC.rds")
saveRDS(thing$mat_G, file="./CC/elas_G_CC.rds")
saveRDS(thing$mat_F, file="./CC/elas_F_CC.rds")
saveRDS(thing$mat_D2, file="./CC/elas_D2_CC.rds")
rm(elas_CC, thing)
sens_CC <- readRDS("./CC/sens_CC.rds")
thing <- decompose(sens_CC)
saveRDS(thing$mat_D1, file="./CC/sens_D1_CC.rds")
saveRDS(thing$mat_G, file="./CC/sens_G_CC.rds")
saveRDS(thing$mat_F, file="./CC/sens_F_CC.rds")
saveRDS(thing$mat_D2, file="./CC/sens_D2_CC.rds")
rm(sens_CC, thing)

#Chekika
elas_C <- readRDS("./C/elas_C.rds")
thing <- decompose(elas_C)
saveRDS(thing$mat_D1, file="./C/elas_D1_C.rds")
saveRDS(thing$mat_G, file="./C/elas_G_C.rds")
saveRDS(thing$mat_F, file="./C/elas_F_C.rds")
saveRDS(thing$mat_D2, file="./C/elas_D2_C.rds")
rm(elas_C, thing)
sens_C <- readRDS("./C/sens_C.rds")
thing <- decompose(sens_C)
saveRDS(thing$mat_D1, file="./C/sens_D1_C.rds")
saveRDS(thing$mat_G, file="./C/sens_G_C.rds")
saveRDS(thing$mat_F, file="./C/sens_F_C.rds")
saveRDS(thing$mat_D2, file="./C/sens_D2_C.rds")
rm(sens_C, thing)

#Fort Pierce
elas_FP <- readRDS("./FP/elas_FP.rds")
thing <- decompose(elas_FP)
saveRDS(thing$mat_D1, file="./FP/elas_D1_FP.rds")
saveRDS(thing$mat_G, file="./FP/elas_G_FP.rds")
saveRDS(thing$mat_F, file="./FP/elas_F_FP.rds")
saveRDS(thing$mat_D2, file="./FP/elas_D2_FP.rds")
rm(elas_FP, thing)
sens_FP <- readRDS("./FP/sens_FP.rds")
thing <- decompose(sens_FP)
saveRDS(thing$mat_D1, file="./FP/sens_D1_FP.rds")
saveRDS(thing$mat_G, file="./FP/sens_G_FP.rds")
saveRDS(thing$mat_F, file="./FP/sens_F_FP.rds")
saveRDS(thing$mat_D2, file="./FP/sens_D2_FP.rds")
rm(sens_FP, thing)

#Punta Gorda
elas_PG <- readRDS("./PG/elas_PG.rds")
thing <- decompose(elas_PG)
saveRDS(thing$mat_D1, file="./PG/elas_D1_PG.rds")
saveRDS(thing$mat_G, file="./PG/elas_G_PG.rds")
saveRDS(thing$mat_F, file="./PG/elas_F_PG.rds")
saveRDS(thing$mat_D2, file="./PG/elas_D2_PG.rds")
rm(elas_PG, thing)
sens_PG <- readRDS("./PG/sens_PG.rds")
thing <- decompose(sens_PG)
saveRDS(thing$mat_D1, file="./PG/sens_D1_PG.rds")
saveRDS(thing$mat_G, file="./PG/sens_G_PG.rds")
saveRDS(thing$mat_F, file="./PG/sens_F_PG.rds")
saveRDS(thing$mat_D2, file="./PG/sens_D2_PG.rds")
rm(sens_PG, thing)

#Wild Turkey
elas_WT <- readRDS("./WT/elas_WT.rds")
thing <- decompose(elas_WT)
saveRDS(thing$mat_D1, file="./WT/elas_D1_WT.rds")
saveRDS(thing$mat_G, file="./WT/elas_G_WT.rds")
saveRDS(thing$mat_F, file="./WT/elas_F_WT.rds")
saveRDS(thing$mat_D2, file="./WT/elas_D2_WT.rds")
rm(elas_WT, thing)
sens_WT <- readRDS("./WT/sens_WT.rds")
thing <- decompose(sens_WT)
saveRDS(thing$mat_D1, file="./WT/sens_D1_WT.rds")
saveRDS(thing$mat_G, file="./WT/sens_G_WT.rds")
saveRDS(thing$mat_F, file="./WT/sens_F_WT.rds")
saveRDS(thing$mat_D2, file="./WT/sens_D2_WT.rds")
rm(sens_WT, thing)

#Overall
elas_overall <- readRDS("./Overall/elas_overall.rds")
thing <- decompose(elas_overall)
saveRDS(thing$mat_D1, file="./Overall/elas_D1_overall.rds")
saveRDS(thing$mat_G, file="./Overall/elas_G_overall.rds")
saveRDS(thing$mat_F, file="./Overall/elas_F_overall.rds")
saveRDS(thing$mat_D2, file="./Overall/elas_D2_overall.rds")
rm(elas_overall, thing)
sens_overall<- readRDS("./Overall/sens_overall.rds")
thing <- decompose(sens_overall)
saveRDS(thing$mat_D1, file="./Overall/sens_D1_overall.rds")
saveRDS(thing$mat_G, file="./Overall/sens_G_overall.rds")
saveRDS(thing$mat_F, file="./Overall/sens_F_overall.rds")
saveRDS(thing$mat_D2, file="./Overall/sens_D2_overall.rds")
rm(sens_overall, thing)



# Convert elasticity components to Kvals #####

toKvals_D1 = function(elas_D1) {
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
  return(Kvals_elas_D1)
}

elas_D1_BC <- readRDS("./BC/elas_D1_BC.rds")
elas_D1_CC <- readRDS("./CC/elas_D1_CC.rds")
elas_D1_C  <- readRDS("./C/elas_D1_C.rds")
elas_D1_FP <- readRDS("./FP/elas_D1_FP.rds")
elas_D1_PG <- readRDS("./PG/elas_D1_PG.rds")
elas_D1_WT <- readRDS("./WT/elas_D1_WT.rds")
elas_D1_overall <- readRDS("./Overall/elas_D1_overall.rds")

Kvals_elas_D1_BC <- toKvals_D1(elas_D1_BC)
save(Kvals_elas_D1_BC, file="./BC/Kvals_elas_D1_BC.RData")
Kvals_elas_D1_CC <- toKvals_D1(elas_D1_CC)
save(Kvals_elas_D1_CC, file="./CC/Kvals_elas_D1_CC.RData")
Kvals_elas_D1_C <- toKvals_D1(elas_D1_C)
save(Kvals_elas_D1_C, file="./C/Kvals_elas_D1_C.RData")
Kvals_elas_D1_FP <- toKvals_D1(elas_D1_FP)
save(Kvals_elas_D1_FP, file="./FP/Kvals_elas_D1_FP.RData")
Kvals_elas_D1_PG <- toKvals_D1(elas_D1_PG)
save(Kvals_elas_D1_PG, file="./PG/Kvals_elas_D1_PG.RData")
Kvals_elas_D1_WT <- toKvals_D1(elas_D1_WT)
save(Kvals_elas_D1_WT, file="./WT/Kvals_elas_D1_WT.RData")
Kvals_elas_D1_overall <- toKvals_D1(elas_D1_overall)
save(Kvals_elas_D1_overall, file="./Overall/Kvals_elas_D1_overall.RData")




rm(Kvals_elas_D1_BC, Kvals_elas_D1_CC, Kvals_elas_D1_C, Kvals_elas_D1_FP,
   Kvals_elas_D1_PG, Kvals_elas_D1_WT, Kvals_elas_D1_overall)
rm(elas_D1_BC, elas_D1_CC, elas_D1_C, elas_D1_FP, elas_D1_PG, elas_D1_WT, elas_D1_overall)

###Construct D2 (Large Domain):


toKvals_D2AA = function(elas_D2AA) {
  plop=function(i,j) {(j-1)*m3a+i} # for putting values in proper place in A 
  Plop=outer(1:m3a,1:m4a,plop); 
  
  
  Kvals_elas_D2AA=array(0,c(m3a,m4a,m3a,m4a));  
  
  for(i in 1:m3a){
    for(j in 1:m4a){
      for(k in 1:m3a){
        kvals= elas_D2AA[Plop[k,1:m4a],Plop[i,j]]
        Kvals_elas_D2AA[k,1:m4a,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }		
  return(Kvals_elas_D2AA)
}

toKvals_D2BB = function(elas_D2BB) {
  plop=function(i,j) {(j-1)*m3b+i} # for putting values in proper place in A 
  Plop=outer(1:m3b,1:m4b,plop); 
  
  
  Kvals_elas_D2BB=array(0,c(m3b,m4b,m3b,m4b));  
  
  for(i in 1:m3b){
    for(j in 1:m4b){
      for(k in 1:m3b){
        kvals= elas_D2BB[Plop[k,1:m4b],Plop[i,j]]
        Kvals_elas_D2BB[k,1:m4b,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }		
  return(Kvals_elas_D2BB)
}

toKvals_D2BA = function(elas_D2BA) {
  plop1=function(i, j) {(j-1)*m3a + i}
  plop2=function(i, j) {(j-1)*m3b + i}
  Plop1=outer(1:m3a,1:m4a,plop1); 
  Plop2=outer(1:m3b, 1:m4b, plop2);

  
  Kvals_elas_D2BA=array(0, c(m3a, m4a, m3b, m4b))
  
  for(i in 1:m3b) {
    for (j in 1:m4b) {
      for (k in 1:m3a) {
        kvals=elas_D2BA[Plop1[k, 1:m4a], Plop2[i,j]]
        Kvals_elas_D2BA[k, 1:m4a, i, j]=kvals
      }}
    cat(i, "\n");
  }
  return(Kvals_elas_D2BA)
}

toKvals_D2AB = function(elas_D2AB) {
  plop1=function(i, j) {(j-1)*m3b + i}
  plop2=function(i, j) {(j-1)*m3a + i}
  Plop1=outer(1:m3b,1:m4b,plop1); 
  Plop2=outer(1:m3a, 1:m4a, plop2);
  
  Kvals_elas_D2AB=array(0, c(m3b, m4b, m3a, m4a))
  
  for(i in 1:m3a) {
    for (j in 1:m4a) {
      for (k in 1:m3b) {
        kvals=elas_D2AB[Plop1[k, 1:m4b], Plop2[i,j]]
        Kvals_elas_D2AB[k, 1:m4b, i, j]=kvals
      }}
    cat(i, "\n");
  }
  return(Kvals_elas_D2AB)
}




#START HERE: Double check this through 
decompose_D2 = function (mat) {
  #Break the elasticity matrix back into its component parts: 
  mat_D2AA<-mat[1:(m3a*m4a), 1:(m3a*m4a)]
  mat_D2BA<-mat[1:(m3a*m4a), ((m3a*m4a)+1):(m3a*m4a + m3b*m4b)]
  mat_D2AB<-mat[((m3a*m4a)+1):(m3a*m4a + m3b*m4b), 1:(m3a*m4a)]
  mat_D2BB<-mat[((m3a*m4a)+1):(m3a*m4a + m3b*m4b), ((m3a*m4a)+1):(m3a*m4a +m3b*m4b)]
  return(list(mat_D2AA = mat_D2AA,
              mat_D2BA = mat_D2BA,
              mat_D2AB = mat_D2AB,
              mat_D2BB = mat_D2BB))
}



elas_D2_BC <- readRDS("./BC/elas_D2_BC.rds")
thing <- decompose_D2(elas_D2_BC)
Kvals_elas_D2AA_BC <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_BC <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_BC <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_BC <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_BC, file="./BC/Kvals_elas_D2AA_BC.RData")
save(Kvals_elas_D2AB_BC, file="./BC/Kvals_elas_D2AB_BC.RData")
save(Kvals_elas_D2BA_BC, file="./BC/Kvals_elas_D2BA_BC.RData")
save(Kvals_elas_D2BB_BC, file="./BC/Kvals_elas_D2BB_BC.RData")
rm(Kvals_elas_D2AA_BC, Kvals_elas_D2AB_BC, Kvals_elas_D2BA_BC, 
   Kvals_elas_D2BB_BC)
rm(thing)
gc()

elas_D2_CC <- readRDS("./CC/elas_D2_CC.rds")
thing <- decompose_D2(elas_D2_CC)
Kvals_elas_D2AA_CC <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_CC <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_CC <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_CC <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_CC, file="./CC/Kvals_elas_D2AA_CC.RData")
save(Kvals_elas_D2AB_CC, file="./CC/Kvals_elas_D2AB_CC.RData")
save(Kvals_elas_D2BA_CC, file="./CC/Kvals_elas_D2BA_CC.RData")
save(Kvals_elas_D2BB_CC, file="./CC/Kvals_elas_D2BB_CC.RData")
rm(Kvals_elas_D2AA_CC, Kvals_elas_D2AB_CC, Kvals_elas_D2BA_CC, 
   Kvals_elas_D2BB_CC)
rm(elas_D2_CC, thing)
gc()

elas_D2_C  <- readRDS("./C/elas_D2_C.rds")
thing <- decompose_D2(elas_D2_C)
Kvals_elas_D2AA_C <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_C <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_C <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_C <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_C, file="./C/Kvals_elas_D2AA_C.RData")
save(Kvals_elas_D2AB_C, file="./C/Kvals_elas_D2AB_C.RData")
save(Kvals_elas_D2BA_C, file="./C/Kvals_elas_D2BA_C.RData")
save(Kvals_elas_D2BB_C, file="./C/Kvals_elas_D2BB_C.RData")
rm(Kvals_elas_D2AA_C, Kvals_elas_D2AB_C, Kvals_elas_D2BA_C, 
   Kvals_elas_D2BB_C)
rm(elas_D2_C, thing)
gc()

elas_D2_FP <- readRDS("./FP/elas_D2_FP.rds")
thing <- decompose_D2(elas_D2_FP)
Kvals_elas_D2AA_FP <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_FP <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_FP <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_FP <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_FP, file="./FP/Kvals_elas_D2AA_FP.RData")
save(Kvals_elas_D2AB_FP, file="./FP/Kvals_elas_D2AB_FP.RData")
save(Kvals_elas_D2BA_FP, file="./FP/Kvals_elas_D2BA_FP.RData")
save(Kvals_elas_D2BB_FP, file="./FP/Kvals_elas_D2BB_FP.RData")
rm(Kvals_elas_D2AA_FP, Kvals_elas_D2AB_FP, Kvals_elas_D2BA_FP, 
   Kvals_elas_D2BB_FP)
rm(elas_D2_FP, thing)
gc()


elas_D2_PG <- readRDS("./PG/elas_D2_PG.rds")
thing <- decompose_D2(elas_D2_PG)
Kvals_elas_D2AA_PG <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_PG <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_PG <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_PG <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_PG, file="./PG/Kvals_elas_D2AA_PG.RData")
save(Kvals_elas_D2AB_PG, file="./PG/Kvals_elas_D2AB_PG.RData")
save(Kvals_elas_D2BA_PG, file="./PG/Kvals_elas_D2BA_PG.RData")
save(Kvals_elas_D2BB_PG, file="./PG/Kvals_elas_D2BB_PG.RData")
rm(Kvals_elas_D2AA_PG, Kvals_elas_D2AB_PG, Kvals_elas_D2BA_PG, 
   Kvals_elas_D2BB_PG)
rm(elas_D2_PG, thing)
gc()

elas_D2_WT <- readRDS("./WT/elas_D2_WT.rds")
thing <- decompose_D2(elas_D2_WT)
Kvals_elas_D2AA_WT <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_WT <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_WT <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_WT <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_WT, file="./WT/Kvals_elas_D2AA_WT.RData")
save(Kvals_elas_D2AB_WT, file="./WT/Kvals_elas_D2AB_WT.RData")
save(Kvals_elas_D2BA_WT, file="./WT/Kvals_elas_D2BA_WT.RData")
save(Kvals_elas_D2BB_WT, file="./WT/Kvals_elas_D2BB_WT.RData")
rm(Kvals_elas_D2AA_WT, Kvals_elas_D2AB_WT, Kvals_elas_D2BA_WT, 
   Kvals_elas_D2BB_WT)
rm(elas_D2_WT, thing)
gc()

elas_D2_overall <- readRDS("./Overall/elas_D2_overall.rds")
thing <- decompose_D2(elas_D2_overall)
Kvals_elas_D2AA_overall <- toKvals_D2AA(thing$mat_D2AA)
Kvals_elas_D2AB_overall <- toKvals_D2AB(thing$mat_D2AB)
Kvals_elas_D2BA_overall <- toKvals_D2BA(thing$mat_D2BA)
Kvals_elas_D2BB_overall <- toKvals_D2BB(thing$mat_D2BB)

save(Kvals_elas_D2AA_overall, file="./Overall/Kvals_elas_D2AA_overall.RData")
save(Kvals_elas_D2AB_overall, file="./Overall/Kvals_elas_D2AB_overall.RData")
save(Kvals_elas_D2BA_overall, file="./Overall/Kvals_elas_D2BA_overall.RData")
save(Kvals_elas_D2BB_overall, file="./Overall/Kvals_elas_D2BB_overall.RData")
rm(Kvals_elas_D2AA_overall, Kvals_elas_D2AB_overall, Kvals_elas_D2BA_overall, 
   Kvals_elas_D2BB_overall)
rm(elas_D2_overall, thing)
gc()

###Construct F (Fecundity):

toKvals_FA = function(elas_FA) {
  plop1=function(i, j) {(j-1)*m1 + i}
  plop2=function(i, j) {(j-1)*m3a + i}
  Plop1=outer(1:m1,1:m2,plop1); 
  Plop2=outer(1:m3a, 1:m4a, plop2);
  
  Kvals_elas_FA=array(0, c(m1, m2, m3a, m4a))
  
  for(i in 1:m3a) {
    for (j in 1:m4a) {
      for (k in 1:m1) {
        kvals=elas_FA[Plop1[k, 1:m2], Plop2[i,j]]
        Kvals_elas_FA[k, 1:m2, i, j]=kvals
      }}
    cat(i, "\n");
  }
  return(Kvals_elas_FA)
}

toKvals_FB = function(elas_FB) {
  plop1=function(i, j) {(j-1)*m1 + i}
  plop2=function(i, j) {(j-1)*m3b + i}
  Plop1=outer(1:m1,1:m2,plop1); 
  Plop2=outer(1:m3b, 1:m4b, plop2);
  
  Kvals_elas_FB=array(0, c(m1, m2, m3b, m4b))
  
  for(i in 1:m3b) {
    for (j in 1:m4b) {
      for (k in 1:m1) {
        kvals=elas_FB[Plop1[k, 1:m2], Plop2[i,j]]
        Kvals_elas_FB[k, 1:m2, i, j]=kvals
      }}
    cat(i, "\n");
  }
  return(Kvals_elas_FB)
}

#START HERE: Double check this through 
decompose_F = function (mat) {
  #Break the elasticity matrix back into its component parts: 
  mat_FA<-mat[1:m1*m2, 1:(m3a*m4a)]
  mat_FB<-mat[1:m1*m2, ((m3a*m4a)+1):(m3a*m4a + m3b* 4b)]

  return(list(mat_D2AA = mat_D2AA,
              mat_D2BA = mat_D2BA,
              mat_D2AB = mat_D2AB,
              mat_D2BB = mat_D2BB))
}



elas_F_BC <- readRDS("./BC/elas_F_BC.rds")
elas_F_CC <- readRDS("./CC/elas_F_CC.rds")
elas_F_C  <- readRDS("./C/elas_F_C.rds")
elas_F_FP <- readRDS("./FP/elas_F_FP.rds")
elas_F_PG <- readRDS("./PG/elas_F_PG.rds")
elas_F_WT <- readRDS("./WT/elas_F_WT.rds")
elas_F_overall <- readRDS("./Overall/elas_F_overall.rds")


Kvals_elas_F_BC <- toKvals_F(elas_F_BC)
save(Kvals_elas_F_BC, file="./BC/Kvals_elas_F_BC.RData")

Kvals_elas_F_CC <- toKvals_F(elas_F_CC)
save(Kvals_elas_F_CC, file="./CC/Kvals_elas_F_CC.RData")
Kvals_elas_F_C <- toKvals_F(elas_F_C)
save(Kvals_elas_F_C, file="./C/Kvals_elas_F_C.RData")
Kvals_elas_F_FP <- toKvals_F(elas_F_FP)
save(Kvals_elas_F_FP, file="./FP/Kvals_elas_F_FP.RData")
Kvals_elas_F_PG <- toKvals_F(elas_F_PG)
save(Kvals_elas_F_PG, file="./PG/Kvals_elas_F_PG.RData")
Kvals_elas_F_WT <- toKvals_F(elas_F_WT)
save(Kvals_elas_F_WT, file="./WT/Kvals_elas_F_WT.RData")

Kvals_elas_F_overall <- toKvals_F(elas_F_overall)
save(Kvals_elas_F_overall, file="./Overall/Kvals_elas_F_overall.RData")

rm(elas_F_BC, elas_F_CC, elas_F_C, elas_F_FP, elas_F_PG, elas_F_WT, elas_F_overall)
rm(Kvals_elas_F_BC, Kvals_elas_F_CC, Kvals_elas_F_C, Kvals_elas_F_FP,
   Kvals_elas_F_PG, Kvals_elas_F_WT, Kvals_elas_F_overall)




###Construct G (Graduation):

toKvals_G = function(elas_G) {
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
  return(Kvals_elas_G)
}

elas_G_BC <- readRDS("./BC/elas_G_BC.rds")
elas_G_CC <- readRDS("./CC/elas_G_CC.rds")
elas_G_C  <- readRDS("./C/elas_G_C.rds")
elas_G_FP <- readRDS("./FP/elas_G_FP.rds")
elas_G_PG <- readRDS("./PG/elas_G_PG.rds")
elas_G_WT <- readRDS("./WT/elas_G_WT.rds")
elas_G_overall <- readRDS("./Overall/elas_G_overall.rds")

Kvals_elas_G_BC <- toKvals_G(elas_G_BC)
save(Kvals_elas_G_BC, file="./BC/Kvals_elas_G_BC.RData")
Kvals_elas_G_CC <- toKvals_G(elas_G_CC)
save(Kvals_elas_G_CC, file="./CC/Kvals_elas_G_CC.RData")
Kvals_elas_G_C <- toKvals_G(elas_G_C)
save(Kvals_elas_G_C, file="./C/Kvals_elas_G_C.RData")
Kvals_elas_G_FP <- toKvals_G(elas_G_FP)
save(Kvals_elas_G_FP, file="./FP/Kvals_elas_G_FP.RData")
Kvals_elas_G_PG <- toKvals_G(elas_G_PG)
save(Kvals_elas_G_PG, file="./PG/Kvals_elas_G_PG.RData")
Kvals_elas_G_WT <- toKvals_G(elas_G_WT)
save(Kvals_elas_G_WT, file="./WT/Kvals_elas_G_WT.RData")
Kvals_elas_G_overall <- toKvals_overall(elas_G_overall)
save(Kvals_elas_G_overall, file="./Overall/Kvals_elas_G_overall.RData")

rm(elas_G_BC, elas_G_CC, elas_G_C, elas_G_FP, elas_G_PG, elas_G_WT, elas_G_overall)
rm(Kvals_elas_G_BC, Kvals_elas_G_CC, Kvals_elas_G_C, Kvals_elas_G_FP,
   Kvals_elas_G_PG, Kvals_elas_G_WT, Kvals_elas_G_overall)



# Create Total Elasticity Pieces #####
load("./BC/Kvals_elas_D1_BC.RData")
load("./BC/Kvals_elas_D2_BC.RData")
load("./BC/Kvals_elas_F_BC.RData")
load("./BC/Kvals_elas_G_BC.RData")

# Big Cypress
total.elas_D1_BC<-apply(Kvals_elas_D1_BC, c(3,4), sum);
total.elas_D2_BC<-apply(Kvals_elas_D2_BC, c(3,4), sum);
total.elas_F_BC<-apply(Kvals_elas_F_BC, c(3,4), sum);
total.elas_G_BC<-apply(Kvals_elas_G_BC, c(3,4), sum);
save(total.elas_D1_BC, total.elas_D2_BC, total.elas_F_BC, total.elas_G_BC, 
     file="./BC/total.elas.RData")

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


#Overall
load("./Overall/Kvals_elas_D1_overall.RData")
load("./Overall/Kvals_elas_D2_overall.RData")
load("./Overall/Kvals_elas_F_overall.RData")
load("./Overall/Kvals_elas_G_overall.RData")

total.elas_D1_overall<-apply(Kvals_elas_D1_overall, c(3,4), sum);
total.elas_D2_overall<-apply(Kvals_elas_D2_overall, c(3,4), sum);
total.elas_F_overall<-apply(Kvals_elas_F_overall, c(3,4), sum);
total.elas_G_overall<-apply(Kvals_elas_G_overall, c(3,4), sum);
save(total.elas_D1_overall, total.elas_D2_overall, total.elas_F_overall, total.elas_G_overall, 
     file="./Overall/total.elas.RData")



# Separate rv and ssd #####
#Separate:
v_BC <- readRDS("./BC/v_BC.rds")
v_CC <- readRDS("./CC/v_CC.rds")
v_C  <- readRDS("./C/v_C.rds")
v_FP <- readRDS("./FP/v_FP.rds")
v_PG <- readRDS("./PG/v_PG.rds")
v_WT <- readRDS("./WT/v_WT.rds")
v_overall <- readRDS("./Overall/v_overall.rds")

stable.dist_BC <- readRDS("./BC/stable.dist_BC.rds")
stable.dist_CC <- readRDS("./CC/stable.dist_CC.rds")
stable.dist_C  <- readRDS("./C/stable.dist_C.rds")
stable.dist_FP <- readRDS("./FP/stable.dist_FP.rds")
stable.dist_PG <- readRDS("./PG/stable.dist_PG.rds")
stable.dist_WT <- readRDS("./WT/stable.dist_WT.rds")
stable.dist_overall <- readRDS("./Overall/stable.dist_overall.rds")


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

repro.val_smalls_overall=matrix(v_overall[1:(m1*m2)], m1,m2);
stable.state_smalls_overall=matrix(stable.dist_overall[1:(m1*m2)], m1, m2);


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

repro.val_larges_overall=matrix(v_overall[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges_overall=matrix(stable.dist_overall[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);

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

#Wild Turkey
ss_diam1_overall<-apply(stable.state_smalls_overall, 1, sum)
ss_diam2_overall<-apply(stable.state_larges_overall, 1, sum)
ss_diam_overall<-c(ss_diam1_overall, ss_diam2_overall)
y_diam<-c(y1, y3)

save(y_diam, file="./Overall/y_diam.RData")
save(ss_diam_overall, file="./Overall/ss_diam_overall.Rdata")


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

#Overall
ss_height1_overall<-apply(stable.state_smalls_overall, 2, sum)
ss_height2_overall<-apply(stable.state_larges_overall, 2, sum)
ss_height_overall<-c(ss_height1_overall, ss_height2_overall)
y_height<-c(y2, y4)
save(y_height, file="./Overall/y_height.RData")
save(ss_height_overall, file="./Overall/ss_height_overall.RData")






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

#Overall
rv_diam1_overall<-apply(repro.val_smalls_overall, 1, sum)
rv_diam2_overall<-apply(repro.val_larges_overall, 1, sum)
rv_diam_overall<-c(rv_diam1_overall, rv_diam2_overall)
save(rv_diam_overall, file="./Overall/rv_diam_Overall.RData")

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

#Overall
rv_height1_overall<-apply(repro.val_smalls_overall, 2, sum)
rv_height2_overall<-apply(repro.val_larges_overall, 2, sum)
rv_height_overall<-c(rv_height1_overall, rv_height2_overall)
save(rv_height_overall, file="./Overall/rv_height_overall.RData")









