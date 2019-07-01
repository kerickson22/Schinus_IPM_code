# Perform a Life Table Response Experiment analysis using the IPM constructed from 
# demography_sites.R. Makes use of code from Ellner and Rees. 
# Modified by Kelley D. Erickson and Carol C. Horvitz


rm(list=ls())

#LTRE
m1=10
m2=m1+1
m3=100
m4=m3+1
tol=1.e-8; 


h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
h3=(700-1.6)/m3;
y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
h4=(800-16)/m4
 y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2

#Calculate A_avg and associated quantities#####
 load("./BC/A_BC.RData")
 load("./CC/A_CC.RData")
 load("./C/A_C.RData")
 load("./FP/A_FP.RData")
 load("./PG/A_PG.RData")
 load("./WT/A_WT.RData")
 
 A_avg <- (A_BC + A_CC + A_C + A_FP + A_PG + A_WT)/6
 
 save(A_avg, file="A_avg.RData")
 rm(A_CC, A_C, A_FP, A_PG, A_WT)
 
#Now calculate the sensitivity of the average matrix: 

 find_lambda = function(A) {
   
   A2=Matrix(A); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
   
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
   qmax=sum(abs(lam*nt-A%*%nt)); 
   cat("Convergence: ",qmax," should be less than ",tol,"\n");
   
   #Find the reproductive value function by iteration
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
   
   return(list(lam.stable = lam.stable, stable.dist = stable.dist, v=v))
 }
 elasticity = function(v, stable.dist, A, lam.stable) {
   v.dot.w<-sum(t(v)%*%stable.dist)
   norm_v<-v/v.dot.w
   #check<-t(norm_v)%*%stable.dist #should be 1
   
   rv<-norm_v
   
   sens<-norm_v%*%t(stable.dist)
   elas<-sens*A/lam.stable
   return(list(sens=sens, elas=elas))
 } 
 
 thing_avg <-find_lambda(A_avg)
 thing_elas <- elasticity(thing_avg$v, thing_avg$stable.dist, A_avg, thing_avg$lam.stable)
 sens_avg <- thing_elas$elas
 save(sens_avg, file="./Average/sens_avg.RData")

# Calculate difference matrices ##### 
# Big Cypress 
load("./BC/A_BC.RData")
load("./Overall/A_overall.RData")
sens_overall <- readRDS("./Overall/sens_overall.rds")
D_BC_avg <-A_BC-A_avg #Difference matrix
D_BC_pool <- A_BC -A_overall 
rm(A_BC)
C_BC_avg<-D_BC_avg*sens_avg #Matrix of contributions
C_BC_pool <- D_BC_pool*sens_overall
rm(D_BC_avg, D_BC_pool)
save(C_BC_avg, file="./BC/C_BC_avg.RData")
save(C_BC_pool, file="./BC/C_BC_pool.RData")
rm(C_BC_avg, C_BC_pool)

# Cape Canaveral 
load("./CC/A_CC.RData")
D_CC_avg <-A_CC-A_avg #Difference matrix
D_CC_pool <- A_CC -A_overall 
rm(A_CC)
C_CC_avg<-D_CC_avg*sens_avg #Matrix of contributions
C_CC_pool <- D_CC_pool*sens_overall
rm(D_CC_avg, D_CC_pool)
save(C_CC_avg, file="./CC/C_CC_avg.RData")
save(C_CC_pool, file="./CC/C_CC_pool.RData")

#Chekika 
load("./C/A_C.RData")
D_C_avg <-A_C-A_avg #Difference matrix
D_C_pool <- A_C -A_overall 
rm(A_C)
C_C_avg<-D_C_avg*sens_avg #Matrix of contributions
C_C_pool <- D_C_pool*sens_overall
rm(D_C_avg, D_C_pool)
save(C_C_avg, file="./C/C_C_avg.RData")
save(C_C_pool, file="./C/C_C_pool.RData")

#Fort Pierce
load("./FP/A_FP.RData")
D_FP_avg <-A_FP-A_avg #Difference matrix
D_FP_pool <- A_FP -A_overall 
rm(A_FP)
C_FP_avg<-D_FP_avg*sens_avg #Matrix of contributions
C_FP_pool <- D_FP_pool*sens_overall
rm(D_FP_avg, D_FP_pool)
save(C_FP_avg, file="./FP/C_FP_avg.RData")
save(C_FP_pool, file="./FP/C_FP_pool.RData")

#Punta Gorda
load("./PG/A_PG.RData")
D_PG_avg <-A_PG-A_avg #Difference matrix
D_PG_pool <- A_PG -A_overall 
rm(A_PG)
C_PG_avg<-D_PG_avg*sens_avg #Matrix of contributions
C_PG_pool <- D_PG_pool*sens_overall
rm(D_PG_avg, D_PG_pool)
save(C_PG_avg, file="./PG/C_PG_avg.RData")
save(C_PG_pool, file="./PG/C_PG_pool.RData")

#Wild Turkey
#Punta Gorda
load("./WT/A_WT.RData")
D_WT_avg <-A_WT-A_avg #Difference matrix
D_WT_pool <- A_WT -A_overall 
rm(A_WT)
C_WT_avg<-D_WT_avg*sens_avg #Matrix of contributions
C_WT_pool <- D_WT_pool*sens_overall
rm(D_WT_avg, D_WT_pool)
save(C_WT_avg, file="./WT/C_WT_avg.RData")
save(C_WT_pool, file="./WT/C_WT_pool.RData")

# Decompose contribution matrices into components #####
decompose = function (mat) {
  #Break the elasticity matrix back into its component parts: 
  mat_D1<-mat[1:(m1*m2), 1:(m1*m2)]
  mat_G<-mat[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
  mat_F<-mat[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
  mat_D2<-mat[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]
  return(list(mat_D1 = mat_D1, mat_G = mat_G, mat_F = mat_F, mat_D2 = mat_D2))
}

load("./BC/C_BC_pool.RData")
load("./CC/C_CC_pool.RData")
load("./C/C_C_pool.RData")
load("./FP/C_FP_pool.RData")
load("./PG/C_PG_pool.RData")
load("./WT/C_WT_pool.RData")

load("./BC/C_BC_avg.RData")
load("./CC/C_CC_avg.RData")
load("./C/C_C_avg.RData")
load("./FP/C_FP_avg.RData")
load("./PG/C_PG_avg.RData")
load("./WT/C_WT_avg.RData")

#Big Cypress
thing <- decompose(C_BC_pool)
saveRDS(thing$mat_D1, file="./BC/C_D1_BC_pool.rds")
saveRDS(thing$mat_G, file="./BC/C_G_BC_pool.rds")
saveRDS(thing$mat_F, file="./BC/C_F_BC_pool.rds")
saveRDS(thing$mat_D2, file="./BC/C_D2_BC_pool.rds")

thing <- decompose(C_BC_avg)
saveRDS(thing$mat_D1, file="./BC/C_D1_BC_avg.rds")
saveRDS(thing$mat_G, file="./BC/C_G_BC_avg.rds")
saveRDS(thing$mat_F, file="./BC/C_F_BC_avg.rds")
saveRDS(thing$mat_D2, file="./BC/C_D2_BC_avg.rds")

#Cape Canaveral
thing <- decompose(C_CC_pool)
saveRDS(thing$mat_D1, file="./CC/C_D1_CC_pool.rds")
saveRDS(thing$mat_G, file="./CC/C_G_CC_pool.rds")
saveRDS(thing$mat_F, file="./CC/C_F_CC_pool.rds")
saveRDS(thing$mat_D2, file="./CC/C_D2_CC_pool.rds")

thing <- decompose(C_CC_avg)
saveRDS(thing$mat_D1, file="./CC/C_D1_CC_avg.rds")
saveRDS(thing$mat_G, file="./CC/C_G_CC_avg.rds")
saveRDS(thing$mat_F, file="./CC/C_F_CC_avg.rds")
saveRDS(thing$mat_D2, file="./CC/C_D2_CC_avg.rds")

#Chekika
thing <- decompose(C_C_pool)
saveRDS(thing$mat_D1, file="./C/C_D1_C_pool.rds")
saveRDS(thing$mat_G, file="./C/C_G_C_pool.rds")
saveRDS(thing$mat_F, file="./C/C_F_C_pool.rds")
saveRDS(thing$mat_D2, file="./C/C_D2_C_pool.rds")

thing <- decompose(C_C_avg)
saveRDS(thing$mat_D1, file="./C/C_D1_C_avg.rds")
saveRDS(thing$mat_G, file="./C/C_G_C_avg.rds")
saveRDS(thing$mat_F, file="./C/C_F_C_avg.rds")
saveRDS(thing$mat_D2, file="./C/C_D2_C_avg.rds")

#Fort Pierce
thing <- decompose(C_FP_pool)
saveRDS(thing$mat_D1, file="./FP/C_D1_FP_pool.rds")
saveRDS(thing$mat_G, file="./FP/C_G_FP_pool.rds")
saveRDS(thing$mat_F, file="./FP/C_F_FP_pool.rds")
saveRDS(thing$mat_D2, file="./FP/C_D2_FP_pool.rds")

thing <- decompose(C_FP_avg)
saveRDS(thing$mat_D1, file="./FP/C_D1_FP_avg.rds")
saveRDS(thing$mat_G, file="./FP/C_G_FP_avg.rds")
saveRDS(thing$mat_F, file="./FP/C_F_FP_avg.rds")
saveRDS(thing$mat_D2, file="./FP/C_D2_FP_avg.rds")

#Punta Gorda
thing <- decompose(C_PG_pool)
saveRDS(thing$mat_D1, file="./PG/C_D1_PG_pool.rds")
saveRDS(thing$mat_G, file="./PG/C_G_PG_pool.rds")
saveRDS(thing$mat_F, file="./PG/C_F_PG_pool.rds")
saveRDS(thing$mat_D2, file="./PG/C_D2_PG_pool.rds")

thing <- decompose(C_PG_avg)
saveRDS(thing$mat_D1, file="./PG/C_D1_PG_avg.rds")
saveRDS(thing$mat_G, file="./PG/C_G_PG_avg.rds")
saveRDS(thing$mat_F, file="./PG/C_F_PG_avg.rds")
saveRDS(thing$mat_D2, file="./PG/C_D2_PG_avg.rds")

#Wild Turkey
thing <- decompose(C_WT_pool)
saveRDS(thing$mat_D1, file="./WT/C_D1_WT_pool.rds")
saveRDS(thing$mat_G, file="./WT/C_G_WT_pool.rds")
saveRDS(thing$mat_F, file="./WT/C_F_WT_pool.rds")
saveRDS(thing$mat_D2, file="./WT/C_D2_WT_pool.rds")

thing <- decompose(C_WT_avg)
saveRDS(thing$mat_D1, file="./WT/C_D1_WT_avg.rds")
saveRDS(thing$mat_G, file="./WT/C_G_WT_avg.rds")
saveRDS(thing$mat_F, file="./WT/C_F_WT_avg.rds")
saveRDS(thing$mat_D2, file="./WT/C_D2_WT_avg.rds")



# Unpack decomposed matrix into pieces #####
toKvals_D1 = function(C_D1) {
  plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in array
  Plop=outer(1:m1,1:m2,plop); 
  
  Kvals_contrib_D1=array(0,c(m1,m2,m1,m2));  
  for(i in 1:m1){
    for(j in 1:m2){
      for(k in 1:m1){
        kvals= C_D1[Plop[k,1:m2],Plop[i,j]]
        Kvals_contrib_D1[k,1:m2,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }
  return(Kvals_contrib_D1)
}

#Big Cypress
C_D1_BC_pool <-readRDS("./BC/C_D1_BC_pool.rds")
C_D1_BC_avg  <-readRDS("./BC/C_D1_BC_avg.rds")
Kvals_contrib_D1_BC_pool <- toKvals_D1(C_D1_BC_pool)
Kvals_contrib_D1_BC_avg <- toKvals_D1(C_D1_BC_avg)
save(Kvals_contrib_D1_BC_pool, file="./BC/Kvals_contrib_D1_BC_pool.RData")
save(Kvals_contrib_D1_BC_avg, file="./BC/Kvals_contrib_D1_BC_avg.RData")

#Cape Canaveral
C_D1_CC_pool <-readRDS("./CC/C_D1_CC_pool.rds")
C_D1_CC_avg  <-readRDS("./CC/C_D1_CC_avg.rds")
Kvals_contrib_D1_CC_pool <- toKvals_D1(C_D1_CC_pool)
Kvals_contrib_D1_CC_avg <- toKvals_D1(C_D1_CC_avg)
save(Kvals_contrib_D1_CC_pool, file="./CC/Kvals_contrib_D1_CC_pool.RData")
save(Kvals_contrib_D1_CC_avg, file="./CC/Kvals_contrib_D1_CC_avg.RData")

#Chekika
C_D1_C_pool <-readRDS("./C/C_D1_C_pool.rds")
C_D1_C_avg  <-readRDS("./C/C_D1_C_avg.rds")
Kvals_contrib_D1_C_pool <- toKvals_D1(C_D1_C_pool)
Kvals_contrib_D1_C_avg <- toKvals_D1(C_D1_C_avg)
save(Kvals_contrib_D1_C_pool, file="./C/Kvals_contrib_D1_C_pool.RData")
save(Kvals_contrib_D1_C_avg, file="./C/Kvals_contrib_D1_C_avg.RData")

#Fort Pierce
C_D1_FP_pool <-readRDS("./FP/C_D1_FP_pool.rds")
C_D1_FP_avg  <-readRDS("./FP/C_D1_FP_avg.rds")
Kvals_contrib_D1_FP_pool <- toKvals_D1(C_D1_FP_pool)
Kvals_contrib_D1_FP_avg <- toKvals_D1(C_D1_FP_avg)
save(Kvals_contrib_D1_FP_pool, file="./FP/Kvals_contrib_D1_FP_pool.RData")
save(Kvals_contrib_D1_FP_avg, file="./FP/Kvals_contrib_D1_FP_avg.RData")

#Punta Gorda
C_D1_PG_pool <-readRDS("./PG/C_D1_PG_pool.rds")
C_D1_PG_avg  <-readRDS("./PG/C_D1_PG_avg.rds")
Kvals_contrib_D1_PG_pool <- toKvals_D1(C_D1_PG_pool)
Kvals_contrib_D1_PG_avg <- toKvals_D1(C_D1_PG_avg)
save(Kvals_contrib_D1_PG_pool, file="./PG/Kvals_contrib_D1_PG_pool.RData")
save(Kvals_contrib_D1_PG_avg, file="./PG/Kvals_contrib_D1_PG_avg.RData")

#Wild Turkey
C_D1_WT_pool <-readRDS("./WT/C_D1_WT_pool.rds")
C_D1_WT_avg  <-readRDS("./WT/C_D1_WT_avg.rds")
Kvals_contrib_D1_WT_pool <- toKvals_D1(C_D1_WT_pool)
Kvals_contrib_D1_WT_avg <- toKvals_D1(C_D1_WT_avg)
save(Kvals_contrib_D1_WT_pool, file="./WT/Kvals_contrib_D1_WT_pool.RData")
save(Kvals_contrib_D1_WT_avg, file="./WT/Kvals_contrib_D1_WT_avg.RData")




rm(list=ls())
m1=10
m2=m1+1
m3=100
m4=m3+1
tol=1.e-8; 


h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
h3=(700-1.6)/m3;
y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
h4=(800-16)/m4
y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2

toKvals_D2 = function(contrib_D2) {
  plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
  Plop=outer(1:m3,1:m4,plop); 
  
  
  Kvals_contrib_D2=array(0,c(m3,m4,m3,m4));  
  
  for(i in 1:m3){
    for(j in 1:m4){
      for(k in 1:m3){
        kvals= contrib_D2[Plop[k,1:m4],Plop[i,j]]
        Kvals_contrib_D2[k,1:m4,i,j]=kvals
        
      }}
    cat(i,"\n"); 
  }		
  return(Kvals_contrib_D2)
}


#Big Cypress
C_D2_BC_pool <-readRDS("./BC/C_D2_BC_pool.rds")
C_D2_BC_avg <- readRDS("./BC/C_D2_BC_avg.rds")
Kvals_contrib_D2_BC_pool <- toKvals_D2(C_D2_BC_pool)
Kvals_contrib_D2_BC_avg <- toKvals_D2(C_D2_BC_avg)
save(Kvals_contrib_D2_BC_pool, file="./BC/Kvals_contrib_D2_BC_pool.RData")
save(Kvals_contrib_D2_BC_avg, file="./BC/Kvals_contrib_D2_BC_avg.RData")

#Cape Canaveral
C_D2_CC_pool <-readRDS("./CC/C_D2_CC_pool.rds")
C_D2_CC_avg <- readRDS("./CC/C_D2_CC_avg.rds")
Kvals_contrib_D2_CC_pool <- toKvals_D2(C_D2_CC_pool)
Kvals_contrib_D2_CC_avg <- toKvals_D2(C_D2_CC_avg)
save(Kvals_contrib_D2_CC_pool, file="./CC/Kvals_contrib_D2_CC_pool.RData")
save(Kvals_contrib_D2_CC_avg, file="./CC/Kvals_contrib_D2_CC_avg.RData")

#Chekika
C_D2_C_pool <-readRDS("./C/C_D2_C_pool.rds")
C_D2_C_avg <- readRDS("./C/C_D2_C_avg.rds")
Kvals_contrib_D2_C_pool <- toKvals_D2(C_D2_C_pool)
Kvals_contrib_D2_C_avg <- toKvals_D2(C_D2_C_avg)
save(Kvals_contrib_D2_C_pool, file="./C/Kvals_contrib_D2_C_pool.RData")
save(Kvals_contrib_D2_C_avg, file="./C/Kvals_contrib_D2_C_avg.RData")

#Fort Pierce
C_D2_FP_pool <-readRDS("./FP/C_D2_FP_pool.rds")
C_D2_FP_avg <- readRDS("./FP/C_D2_FP_avg.rds")
Kvals_contrib_D2_FP_pool <- toKvals_D2(C_D2_FP_pool)
Kvals_contrib_D2_FP_avg <- toKvals_D2(C_D2_FP_avg)
save(Kvals_contrib_D2_FP_pool, file="./FP/Kvals_contrib_D2_FP_pool.RData")
save(Kvals_contrib_D2_FP_avg, file="./FP/Kvals_contrib_D2_FP_avg.RData")

#Punta Gorda
C_D2_PG_pool <-readRDS("./PG/C_D2_PG_pool.rds")
C_D2_PG_avg <- readRDS("./PG/C_D2_PG_avg.rds")
Kvals_contrib_D2_PG_pool <- toKvals_D2(C_D2_PG_pool)
Kvals_contrib_D2_PG_avg <- toKvals_D2(C_D2_PG_avg)
save(Kvals_contrib_D2_PG_pool, file="./PG/Kvals_contrib_D2_PG_pool.RData")
save(Kvals_contrib_D2_PG_avg, file="./PG/Kvals_contrib_D2_PG_avg.RData")

rm(C_D2_BC_avg, C_D2_BC_pool, C_D2_CC_avg, C_D2_CC_pool, C_D2_C_avg, 
   C_D2_C_pool, C_D2_FP_avg, C_D2_FP_pool, C_D2_PG_avg, C_D2_PG_pool)
#Wild Turkey
C_D2_WT_pool <-readRDS("./WT/C_D2_WT_pool.rds")
C_D2_WT_avg <- readRDS("./WT/C_D2_WT_avg.rds")
Kvals_contrib_D2_WT_pool <- toKvals_D2(C_D2_WT_pool)
Kvals_contrib_D2_WT_avg <- toKvals_D2(C_D2_WT_avg)
save(Kvals_contrib_D2_WT_pool, file="./WT/Kvals_contrib_D2_WT_pool.RData")
save(Kvals_contrib_D2_WT_avg, file="./WT/Kvals_contrib_D2_WT_avg.RData")


toKvals_F = function(elas_F) {
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
  return(Kvals_elas_F)
}

#Big Cypress
C_F_BC_pool <-readRDS("./BC/C_F_BC_pool.rds")
C_F_BC_avg  <-readRDS("./BC/C_F_BC_avg.rds")
Kvals_contrib_F_BC_pool <- toKvals_F(C_F_BC_pool)
Kvals_contrib_F_BC_avg  <- toKvals_F(C_F_BC_avg)
save(Kvals_contrib_F_BC_pool, file="./BC/Kvals_contrib_F_BC_pool.RData")
save(Kvals_contrib_F_BC_avg, file="./BC/Kvals_contrib_F_BC_avg.RData")

#Cape Canaveral
C_F_CC_pool <-readRDS("./CC/C_F_CC_pool.rds")
C_F_CC_avg  <-readRDS("./CC/C_F_CC_avg.rds")
Kvals_contrib_F_CC_pool <- toKvals_F(C_F_CC_pool)
Kvals_contrib_F_CC_avg  <- toKvals_F(C_F_CC_avg)
save(Kvals_contrib_F_CC_pool, file="./CC/Kvals_contrib_F_CC_pool.RData")
save(Kvals_contrib_F_CC_avg, file="./CC/Kvals_contrib_F_CC_avg.RData")

#Chekika
C_F_C_pool <-readRDS("./C/C_F_C_pool.rds")
C_F_C_avg  <-readRDS("./C/C_F_C_avg.rds")
Kvals_contrib_F_C_pool <- toKvals_F(C_F_C_pool)
Kvals_contrib_F_C_avg  <- toKvals_F(C_F_C_avg)
save(Kvals_contrib_F_C_pool, file="./C/Kvals_contrib_F_C_pool.RData")
save(Kvals_contrib_F_C_avg, file="./C/Kvals_contrib_F_C_avg.RData")

#Fort Pierce
C_F_FP_pool <-readRDS("./FP/C_F_FP_pool.rds")
C_F_FP_avg  <-readRDS("./FP/C_F_FP_avg.rds")
Kvals_contrib_F_FP_pool <- toKvals_F(C_F_FP_pool)
Kvals_contrib_F_FP_avg  <- toKvals_F(C_F_FP_avg)
save(Kvals_contrib_F_FP_pool, file="./FP/Kvals_contrib_F_FP_pool.RData")
save(Kvals_contrib_F_FP_avg, file="./FP/Kvals_contrib_F_FP_avg.RData")

#Punta Gorda
C_F_PG_pool <-readRDS("./PG/C_F_PG_pool.rds")
C_F_PG_avg  <-readRDS("./PG/C_F_PG_avg.rds")
Kvals_contrib_F_PG_pool <- toKvals_F(C_F_PG_pool)
Kvals_contrib_F_PG_avg  <- toKvals_F(C_F_PG_avg)
save(Kvals_contrib_F_PG_pool, file="./PG/Kvals_contrib_F_PG_pool.RData")
save(Kvals_contrib_F_PG_avg, file="./PG/Kvals_contrib_F_PG_avg.RData")

#Wild Turkey
C_F_WT_pool <-readRDS("./WT/C_F_WT_pool.rds")
C_F_WT_avg  <-readRDS("./WT/C_F_WT_avg.rds")
Kvals_contrib_F_WT_pool <- toKvals_F(C_F_WT_pool)
Kvals_contrib_F_WT_avg  <- toKvals_F(C_F_WT_avg)
save(Kvals_contrib_F_WT_pool, file="./WT/Kvals_contrib_F_WT_pool.RData")
save(Kvals_contrib_F_WT_avg, file="./WT/Kvals_contrib_F_WT_avg.RData")










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

#Big Cypress
C_G_BC_pool <-readRDS("./BC/C_G_BC_pool.rds")
C_G_BC_avg  <-readRDS("./BC/C_G_BC_avg.rds")
Kvals_contrib_G_BC_pool <- toKvals_G(C_G_BC_pool)
Kvals_contrib_G_BC_avg  <- toKvals_G(C_G_BC_avg)
save(Kvals_contrib_G_BC_pool, file="./BC/Kvals_contrib_G_BC_pool.RData")
save(Kvals_contrib_G_BC_avg, file="./BC/Kvals_contrib_G_BC_avg.RData")

#Cape Canaveral
C_G_CC_pool <-readRDS("./CC/C_G_CC_pool.rds")
C_G_CC_avg  <-readRDS("./CC/C_G_CC_avg.rds")
Kvals_contrib_G_CC_pool <- toKvals_G(C_G_CC_pool)
Kvals_contrib_G_CC_avg  <- toKvals_G(C_G_CC_avg)
save(Kvals_contrib_G_CC_pool, file="./CC/Kvals_contrib_G_CC_pool.RData")
save(Kvals_contrib_G_CC_avg, file="./CC/Kvals_contrib_G_CC_avg.RData")

#Chekika
C_G_C_pool <-readRDS("./C/C_G_C_pool.rds")
C_G_C_avg  <-readRDS("./C/C_G_C_avg.rds")
Kvals_contrib_G_C_pool <- toKvals_G(C_G_C_pool)
Kvals_contrib_G_C_avg  <- toKvals_G(C_G_C_avg)
save(Kvals_contrib_G_C_pool, file="./C/Kvals_contrib_G_C_pool.RData")
save(Kvals_contrib_G_C_avg, file="./C/Kvals_contrib_G_C_avg.RData")

#Fort Pierce
C_G_FP_pool <-readRDS("./FP/C_G_FP_pool.rds")
C_G_FP_avg  <-readRDS("./FP/C_G_FP_avg.rds")
Kvals_contrib_G_FP_pool <- toKvals_G(C_G_FP_pool)
Kvals_contrib_G_FP_avg  <- toKvals_G(C_G_FP_avg)
save(Kvals_contrib_G_FP_pool, file="./FP/Kvals_contrib_G_FP_pool.RData")
save(Kvals_contrib_G_FP_avg, file="./FP/Kvals_contrib_G_FP_avg.RData")

#Punta Gorda
C_G_PG_pool <-readRDS("./PG/C_G_PG_pool.rds")
C_G_PG_avg  <-readRDS("./PG/C_G_PG_avg.rds")
Kvals_contrib_G_PG_pool <- toKvals_G(C_G_PG_pool)
Kvals_contrib_G_PG_avg  <- toKvals_G(C_G_PG_avg)
save(Kvals_contrib_G_PG_pool, file="./PG/Kvals_contrib_G_PG_pool.RData")
save(Kvals_contrib_G_PG_avg, file="./PG/Kvals_contrib_G_PG_avg.RData")

#Wild Turkey
C_G_WT_pool <-readRDS("./WT/C_G_WT_pool.rds")
C_G_WT_avg  <-readRDS("./WT/C_G_WT_avg.rds")
Kvals_contrib_G_WT_pool <- toKvals_G(C_G_WT_pool)
Kvals_contrib_G_WT_avg  <- toKvals_G(C_G_WT_avg)
save(Kvals_contrib_G_WT_pool, file="./WT/Kvals_contrib_G_WT_pool.RData")
save(Kvals_contrib_G_WT_avg, file="./WT/Kvals_contrib_G_WT_avg.RData")










rm(C_D2_WT_avg, C_D2_WT_pool, C_F_BC_avg, C_F_BC_pool, C_F_CC_avg, 
   C_F_CC_pool, C_F_C_avg, C_F_C_pool, C_F_FP_avg, C_F_FP_pool,
   C_F_PG_avg, C_F_PG_pool, C_F_WT_avg, C_F_WT_pool, C_G_BC_avg,
   C_G_BC_pool, C_G_CC_avg, C_G_CC_pool, C_G_C_avg, C_G_C_pool,
   C_G_FP_avg, C_G_FP_pool, C_G_PG_avg, C_G_PG_pool, C_G_WT_avg, 
   C_G_WT_pool)

# Convert to totals: #####
#Big Cypress

load("./BC/Kvals_contrib_D1_BC_pool.RData")
load("./BC/Kvals_contrib_D2_BC_pool.RData")
load("./BC/Kvals_contrib_F_BC_pool.RData")
load("./BC/Kvals_contrib_G_BC_pool.RData")
total.contrib_D1_BC_pool<-apply(Kvals_contrib_D1_BC_pool, c(3,4), sum);
total.contrib_D2_BC_pool<-apply(Kvals_contrib_D2_BC_pool, c(3,4), sum);
total.contrib_F_BC_pool<-apply(Kvals_contrib_F_BC_pool, c(3,4), sum);
total.contrib_G_BC_pool<-apply(Kvals_contrib_G_BC_pool, c(3,4), sum);
save(total.contrib_D1_BC_pool, total.contrib_D2_BC_pool, total.contrib_F_BC_pool,
     total.contrib_G_BC_pool, 
     file="./BC/total.contrib_BC_pool.RData")
#rm(Kvals_contrib_D1_BC_pool, Kvals_contrib_D2_BC_pool, Kvals_contrib_F_BC_pool,
 #  Kvals_contrib_G_BC_pool)


load("./BC/Kvals_contrib_D1_BC_avg.RData")
load("./BC/Kvals_contrib_D2_BC_avg.RData")
load("./BC/Kvals_contrib_F_BC_avg.RData")
load("./BC/Kvals_contrib_G_BC_avg.RData")
total.contrib_D1_BC_avg<-apply(Kvals_contrib_D1_BC_avg, c(3,4), sum);
total.contrib_D2_BC_avg<-apply(Kvals_contrib_D2_BC_avg, c(3,4), sum);
total.contrib_F_BC_avg<-apply(Kvals_contrib_F_BC_avg, c(3,4), sum);
total.contrib_G_BC_avg<-apply(Kvals_contrib_G_BC_avg, c(3,4), sum);
save(total.contrib_D1_BC_avg, total.contrib_D2_BC_avg, total.contrib_F_BC_avg,
     total.contrib_G_BC_avg, 
     file="./BC/total.contrib_BC_avg.RData")
#rm(Kvals_contrib_D1_BC_avg, Kvals_contrib_D2_BC_avg, Kvals_contrib_F_BC_avg,
   #Kvals_contrib_G_BC_avg)


#Cape Canaveral
load("./CC/Kvals_contrib_D1_CC_pool.RData")
load("./CC/Kvals_contrib_D2_CC_pool.RData")
load("./CC/Kvals_contrib_F_CC_pool.RData")
load("./CC/Kvals_contrib_G_CC_pool.RData")
total.contrib_D1_CC_pool<-apply(Kvals_contrib_D1_CC_pool, c(3,4), sum);
total.contrib_D2_CC_pool<-apply(Kvals_contrib_D2_CC_pool, c(3,4), sum);
total.contrib_F_CC_pool<-apply(Kvals_contrib_F_CC_pool, c(3,4), sum);
total.contrib_G_CC_pool<-apply(Kvals_contrib_G_CC_pool, c(3,4), sum);
save(total.contrib_D1_CC_pool, total.contrib_D2_CC_pool, total.contrib_F_CC_pool,
     total.contrib_G_CC_pool, 
     file="./CC/total.contrib_CC_pool.RData")
#rm(Kvals_contrib_D1_CC_pool, Kvals_contrib_D2_CC_pool, Kvals_contrib_F_CC_pool,
  # Kvals_contrib_G_CC_pool)

load("./CC/Kvals_contrib_D1_CC_avg.RData")
load("./CC/Kvals_contrib_D2_CC_avg.RData")
load("./CC/Kvals_contrib_F_CC_avg.RData")
load("./CC/Kvals_contrib_G_CC_avg.RData")
total.contrib_D1_CC_avg<-apply(Kvals_contrib_D1_CC_avg, c(3,4), sum);
total.contrib_D2_CC_avg<-apply(Kvals_contrib_D2_CC_avg, c(3,4), sum);
total.contrib_F_CC_avg<-apply(Kvals_contrib_F_CC_avg, c(3,4), sum);
total.contrib_G_CC_avg<-apply(Kvals_contrib_G_CC_avg, c(3,4), sum);
save(total.contrib_D1_CC_avg, total.contrib_D2_CC_avg, total.contrib_F_CC_avg,
     total.contrib_G_CC_avg, 
     file="./CC/total.contrib_CC_avg.RData")
#rm(Kvals_contrib_D1_CC_avg, Kvals_contrib_D2_CC_avg, Kvals_contrib_F_CC_avg,
 #  Kvals_contrib_G_CC_avg)

#Chekika
load("./C/Kvals_contrib_D1_C_pool.RData")
load("./C/Kvals_contrib_D2_C_pool.RData")
load("./C/Kvals_contrib_F_C_pool.RData")
load("./C/Kvals_contrib_G_C_pool.RData")
total.contrib_D1_C_pool<-apply(Kvals_contrib_D1_C_pool, c(3,4), sum);
total.contrib_D2_C_pool<-apply(Kvals_contrib_D2_C_pool, c(3,4), sum);
total.contrib_F_C_pool<-apply(Kvals_contrib_F_C_pool, c(3,4), sum);
total.contrib_G_C_pool<-apply(Kvals_contrib_G_C_pool, c(3,4), sum);
save(total.contrib_D1_C_pool, total.contrib_D2_C_pool, total.contrib_F_C_pool,
     total.contrib_G_C_pool, 
     file="./C/total.contrib_C_pool.RData")
#rm(Kvals_contrib_D1_C_pool, Kvals_contrib_D2_C_pool, Kvals_contrib_F_C_pool,
 #  Kvals_contrib_G_C_pool)

load("./C/Kvals_contrib_D1_C_avg.RData")
load("./C/Kvals_contrib_D2_C_avg.RData")
load("./C/Kvals_contrib_F_C_avg.RData")
load("./C/Kvals_contrib_G_C_avg.RData")
total.contrib_D1_C_avg<-apply(Kvals_contrib_D1_C_avg, c(3,4), sum);
total.contrib_D2_C_avg<-apply(Kvals_contrib_D2_C_avg, c(3,4), sum);
total.contrib_F_C_avg<-apply(Kvals_contrib_F_C_avg, c(3,4), sum);
total.contrib_G_C_avg<-apply(Kvals_contrib_G_C_avg, c(3,4), sum);
save(total.contrib_D1_C_avg, total.contrib_D2_C_avg, total.contrib_F_C_avg,
     total.contrib_G_C_avg, 
     file="./C/total.contrib_C_avg.RData")
#rm(Kvals_contrib_D1_C_avg, Kvals_contrib_D2_C_avg, Kvals_contrib_F_C_avg,
  # Kvals_contrib_G_C_avg)

#Fort Pierce
load("./FP/Kvals_contrib_D1_FP_pool.RData")
load("./FP/Kvals_contrib_D2_FP_pool.RData")
load("./FP/Kvals_contrib_F_FP_pool.RData")
load("./FP/Kvals_contrib_G_FP_pool.RData")
total.contrib_D1_FP_pool<-apply(Kvals_contrib_D1_FP_pool, c(3,4), sum);
total.contrib_D2_FP_pool<-apply(Kvals_contrib_D2_FP_pool, c(3,4), sum);
total.contrib_F_FP_pool<-apply(Kvals_contrib_F_FP_pool, c(3,4), sum);
total.contrib_G_FP_pool<-apply(Kvals_contrib_G_FP_pool, c(3,4), sum);
save(total.contrib_D1_FP_pool, total.contrib_D2_FP_pool, total.contrib_F_FP_pool,
     total.contrib_G_FP_pool, 
     file="./FP/total.contrib_FP_pool.RData")
#rm(Kvals_contrib_D1_FP_pool, Kvals_contrib_D2_FP_pool, Kvals_contrib_F_FP_pool,
 #  Kvals_contrib_G_FP_pool)

load("./FP/Kvals_contrib_D1_FP_avg.RData")
load("./FP/Kvals_contrib_D2_FP_avg.RData")
load("./FP/Kvals_contrib_F_FP_avg.RData")
load("./FP/Kvals_contrib_G_FP_avg.RData")
total.contrib_D1_FP_avg<-apply(Kvals_contrib_D1_FP_avg, c(3,4), sum);
total.contrib_D2_FP_avg<-apply(Kvals_contrib_D2_FP_avg, c(3,4), sum);
total.contrib_F_FP_avg<-apply(Kvals_contrib_F_FP_avg, c(3,4), sum);
total.contrib_G_FP_avg<-apply(Kvals_contrib_G_FP_avg, c(3,4), sum);
save(total.contrib_D1_FP_avg, total.contrib_D2_FP_avg, total.contrib_F_FP_avg,
     total.contrib_G_FP_avg, 
     file="./FP/total.contrib_FP_avg.RData")
#rm(Kvals_contrib_D1_FP_avg, Kvals_contrib_D2_FP_avg, Kvals_contrib_F_FP_avg,
  # Kvals_contrib_G_FP_avg)

#Punta Gorda
load("./PG/Kvals_contrib_D1_PG_pool.RData")
load("./PG/Kvals_contrib_D2_PG_pool.RData")
load("./PG/Kvals_contrib_F_PG_pool.RData")
load("./PG/Kvals_contrib_G_PG_pool.RData")
total.contrib_D1_PG_pool<-apply(Kvals_contrib_D1_PG_pool, c(3,4), sum);
total.contrib_D2_PG_pool<-apply(Kvals_contrib_D2_PG_pool, c(3,4), sum);
total.contrib_F_PG_pool<-apply(Kvals_contrib_F_PG_pool, c(3,4), sum);
total.contrib_G_PG_pool<-apply(Kvals_contrib_G_PG_pool, c(3,4), sum);
save(total.contrib_D1_PG_pool, total.contrib_D2_PG_pool, total.contrib_F_PG_pool,
     total.contrib_G_PG_pool, 
     file="./PG/total.contrib_PG_pool.RData")
#rm(Kvals_contrib_D1_PG_pool, Kvals_contrib_D2_PG_pool, Kvals_contrib_F_PG_pool,
  # Kvals_contrib_G_PG_pool)

load("./PG/Kvals_contrib_D1_PG_avg.RData")
load("./PG/Kvals_contrib_D2_PG_avg.RData")
load("./PG/Kvals_contrib_F_PG_avg.RData")
load("./PG/Kvals_contrib_G_PG_avg.RData")
total.contrib_D1_PG_avg<-apply(Kvals_contrib_D1_PG_avg, c(3,4), sum);
total.contrib_D2_PG_avg<-apply(Kvals_contrib_D2_PG_avg, c(3,4), sum);
total.contrib_F_PG_avg<-apply(Kvals_contrib_F_PG_avg, c(3,4), sum);
total.contrib_G_PG_avg<-apply(Kvals_contrib_G_PG_avg, c(3,4), sum);
save(total.contrib_D1_PG_avg, total.contrib_D2_PG_avg, total.contrib_F_PG_avg,
     total.contrib_G_PG_avg, 
     file="./PG/total.contrib_PG_avg.RData")
#rm(Kvals_contrib_D1_PG_avg, Kvals_contrib_D2_PG_avg, Kvals_contrib_F_PG_avg,
 #  Kvals_contrib_G_PG_avg)

#Wild Turkey
load("./WT/Kvals_contrib_D1_WT_pool.RData")
load("./WT/Kvals_contrib_D2_WT_pool.RData")
load("./WT/Kvals_contrib_F_WT_pool.RData")
load("./WT/Kvals_contrib_G_WT_pool.RData")
total.contrib_D1_WT_pool<-apply(Kvals_contrib_D1_WT_pool, c(3,4), sum);
total.contrib_D2_WT_pool<-apply(Kvals_contrib_D2_WT_pool, c(3,4), sum);
total.contrib_F_WT_pool<-apply(Kvals_contrib_F_WT_pool, c(3,4), sum);
total.contrib_G_WT_pool<-apply(Kvals_contrib_G_WT_pool, c(3,4), sum);
save(total.contrib_D1_WT_pool, total.contrib_D2_WT_pool, total.contrib_F_WT_pool,
     total.contrib_G_WT_pool, 
     file="./WT/total.contrib_WT_pool.RData")
#rm(Kvals_contrib_D1_WT_pool, Kvals_contrib_D2_WT_pool, Kvals_contrib_F_WT_pool,
  # Kvals_contrib_G_WT_pool)

load("./WT/Kvals_contrib_D1_WT_avg.RData")
load("./WT/Kvals_contrib_D2_WT_avg.RData")
load("./WT/Kvals_contrib_F_WT_avg.RData")
load("./WT/Kvals_contrib_G_WT_avg.RData")
total.contrib_D1_WT_avg<-apply(Kvals_contrib_D1_WT_avg, c(3,4), sum);
total.contrib_D2_WT_avg<-apply(Kvals_contrib_D2_WT_avg, c(3,4), sum);
total.contrib_F_WT_avg<-apply(Kvals_contrib_F_WT_avg, c(3,4), sum);
total.contrib_G_WT_avg<-apply(Kvals_contrib_G_WT_avg, c(3,4), sum);
save(total.contrib_D1_WT_avg, total.contrib_D2_WT_avg, total.contrib_F_WT_avg,
     total.contrib_G_WT_avg, 
     file="./WT/total.contrib_WT_avg.RData")
#rm(Kvals_contrib_D1_WT_avg, Kvals_contrib_D2_WT_avg, Kvals_contrib_F_WT_avg,
  # Kvals_contrib_G_WT_avg)


#Collapse contributions by site: #####


 BC_pool<-c( 
 sum(rowSums(total.contrib_D1_BC_pool)), 
 sum(rowSums(total.contrib_G_BC_pool)), 
 sum(rowSums(total.contrib_F_BC_pool)), 
 sum(rowSums(total.contrib_D2_BC_pool)) )
 save(BC_pool, file="./BC/BC_pool.RData")
 
 BC_avg<-c( 
   sum(rowSums(total.contrib_D1_BC_avg)), 
   sum(rowSums(total.contrib_G_BC_avg)), 
   sum(rowSums(total.contrib_F_BC_avg)), 
   sum(rowSums(total.contrib_D2_BC_avg)) )
 save(BC_avg, file="./BC/BC_avg.RData")
 
 
 
 CC_pool<-c( 
   sum(rowSums(total.contrib_D1_CC_pool)), 
   sum(rowSums(total.contrib_G_CC_pool)), 
   sum(rowSums(total.contrib_F_CC_pool)), 
   sum(rowSums(total.contrib_D2_CC_pool)) )
 save(CC_pool, file="./CC/CC_pool.RData")
 
 CC_avg<-c( 
   sum(rowSums(total.contrib_D1_CC_avg)), 
   sum(rowSums(total.contrib_G_CC_avg)), 
   sum(rowSums(total.contrib_F_CC_avg)), 
   sum(rowSums(total.contrib_D2_CC_avg)) )
 save(CC_avg, file="./CC/CC_avg.RData")
 
 C_pool<-c( 
   sum(rowSums(total.contrib_D1_C_pool)), 
   sum(rowSums(total.contrib_G_C_pool)), 
   sum(rowSums(total.contrib_F_C_pool)), 
   sum(rowSums(total.contrib_D2_C_pool)) )
 save(C_pool, file="./C/C_pool.RData")
 
 C_avg<-c( 
   sum(rowSums(total.contrib_D1_C_avg)), 
   sum(rowSums(total.contrib_G_C_avg)), 
   sum(rowSums(total.contrib_F_C_avg)), 
   sum(rowSums(total.contrib_D2_C_avg)) )
 save(C_avg, file="./C/C_avg.RData")
 
 FP_pool<-c( 
   sum(rowSums(total.contrib_D1_FP_pool)), 
   sum(rowSums(total.contrib_G_FP_pool)), 
   sum(rowSums(total.contrib_F_FP_pool)), 
   sum(rowSums(total.contrib_D2_FP_pool)) )
 save(FP_pool, file="./FP/FP_pool.RData")
 
 FP_avg<-c( 
   sum(rowSums(total.contrib_D1_FP_avg)), 
   sum(rowSums(total.contrib_G_FP_avg)), 
   sum(rowSums(total.contrib_F_FP_avg)), 
   sum(rowSums(total.contrib_D2_FP_avg)) )
 save(FP_avg, file="./FP/FP_avg.RData")
 
 PG_pool<-c( 
   sum(rowSums(total.contrib_D1_PG_pool)), 
   sum(rowSums(total.contrib_G_PG_pool)), 
   sum(rowSums(total.contrib_F_PG_pool)), 
   sum(rowSums(total.contrib_D2_PG_pool)) )
 save(PG_pool, file="./PG/PG_pool.RData")
 
 PG_avg<-c( 
   sum(rowSums(total.contrib_D1_PG_avg)), 
   sum(rowSums(total.contrib_G_PG_avg)), 
   sum(rowSums(total.contrib_F_PG_avg)), 
   sum(rowSums(total.contrib_D2_PG_avg)) )
 save(PG_avg, file="./PG/PG_avg.RData")
 
 WT_pool<-c( 
   sum(rowSums(total.contrib_D1_WT_pool)), 
   sum(rowSums(total.contrib_G_WT_pool)), 
   sum(rowSums(total.contrib_F_WT_pool)), 
   sum(rowSums(total.contrib_D2_WT_pool)) )
 save(WT_pool, file="./WT/WT_pool.RData")
 
 WT_avg<-c( 
   sum(rowSums(total.contrib_D1_WT_avg)), 
   sum(rowSums(total.contrib_G_WT_avg)), 
   sum(rowSums(total.contrib_F_WT_avg)), 
   sum(rowSums(total.contrib_D2_WT_avg)) )
 save(WT_avg, file="./WT/WT_avg.RData")
 
# Variability in lambda #####
 #With six sites to compare, this part of the
 #code takes a long time to run (over 48 hours at least)
# Calculate matrices of coefficients of variation #####

	
D1_BC<-readRDS("./BC/D1_BC.rds")
D1_PG <- readRDS("./PG/D1_PG.rds")
D1_C<-readRDS("./C/D1_C.rds")
D1_FP<-readRDS("./FP/D1_FP.rds")
D1_PG<-readRDS("./PG/D1_PG.rds")
D1_WT<-readRDS("./WT/D1_WT.rds")
	
CV<-function(x) {
	return(sd(x)/mean(x))
	
}



	
cv_D1<-matrix(nrow=(m1*m2), ncol=(m1*m2))

for(i in 1:(m1*m2)) {
	for(j in 1:(m1*m2)) {
		cv_D1[i,j]<-CV(c(D1_BC[i,j], D1_PG[i,j], D1_C[i,j], D1_FP[i,j], 
		                 D1_PG[i,j], D1_WT[i,j]))
		if(is.na(cv_D1[i,j])) {cv_D1[i,j]<-0}
	}
}


D1_avg=(D1_BC + D1_PG + D1_C + D1_FP + D1_PG + D1_WT)/6
rm(D1_BC, D1_PG, D1_C, D1_FP, D1_PG, D1_WT)

#G
G_BC<-readRDS("./BC/G_BC.rds")
G_PG<-readRDS("./PG/G_PG.rds")
G_C <- readRDS("./C/G_C.rds")
G_FP <- readRDS("./FP/G_FP.rds")
G_PG <- readRDS("./PG/G_PG.rds")
G_WT <- readRDS("./WT/G_WT.rds")

cv_G<-matrix(nrow=(m3*m4), ncol=(m1*m2))


for(i in 1:(m3*m4)) {
	for(j in 1:(m1*m2)) {
		cv_G[i,j]<-CV(c(G_BC[i,j], G_PG[i,j], G_C[i,j], G_FP[i,j],
		                G_PG[i,j], G_WT[i,j]))
		if(is.na(cv_G[i,j])) {cv_G[i,j]<-0}
	}
}

G_avg<-(G_BC+G_PG+G_C + G_FP + G_PG + G_WT)/6

rm(G_BC, G_PG, G_C, G_FP, G_PG, G_WT)

#F
F_BC <- readRDS("./BC/F_BC.rds")
F_PG <- readRDS("./PG/F_PG.rds")
F_C  <- readRDS("./C/F_C.rds")
F_FP <- readRDS("./FP/F_FP.rds")
F_PG <- readRDS("./PG/F_PG.rds")
F_WT <- readRDS("./WT/F_WT.rds")

cv_F<-matrix(nrow=(m1*m2), ncol=(m3*m4))
for(i in 1:(m1*m2)) {
	for(j in 1:(m3*m4)) {
		cv_F[i,j]<-CV(c(F_BC[i,j], F_PG[i,j], F_C[i,j], F_FP[i,j],
		              F_PG[i,j], F_WT[i,j]))
		if(is.na(cv_F[i,j])) {cv_F[i,j]<-0}
	}
}

F_avg<-(F_BC + F_PG + F_C + F_FP + F_PG + F_WT)/6

rm(F_BC, F_PG, F_C, F_FP, F_PG, F_WT)


save(cv_D1, file="cv_D1.RData")
save(cv_F, file="cv_F.RData")
save(cv_G, file="cv_G.RData")
rm(cv_D1, cv_F, cv_G)
#D2
D2_BC <- readRDS("./BC/D2_BC.rds")
D2_PG <- readRDS("./PG/D2_PG.rds")
D2_C  <- readRDS("./C/D2_C.rds")
D2_FP <- readRDS("./FP/D2_FP.rds")
D2_PG <- readRDS("./PG/D2_PG.rds")
D2_WT <- readRDS("./WT/D2_WT.rds")

#Lines 568: 577 are run on a separate computer because they take
# a long time to run (X hours)
cv_D2<-matrix(nrow=(m3*m4), ncol=(m3*m4))
for(i in 1:(m3*m4)) {
	for(j in 1:(m3*m4)) {
		cv_D2[i,j]<-CV(c(D2_BC[i,j], D2_PG[i,j], D2_C[i,j], D2_FP,
		                 D2_PG, D2_WT))
		if(is.na(cv_D2[i,j])) {cv_D2[i,j]<-0}
	}
  cat("i=", i, "\n")
}
save(cv_D2, file="cv_D2.RData")

D2_avg<-(D2_E+D2_H+D2_W)/3
rm(D2_E)
rm(D2_W)
rm(D2_H)


A_cv<-cbind(rbind(cv_D1, cv_G), rbind(cv_F, cv_D2))
rm(cv_D1, cv_G, cv_F, cv_D2)

save(A_cv, file="./Overall/A_cv.RData")

load("./Overall/elas_overall.RData")
thing<-A_cv*elas_overall

save(thing, file="./Overall/thing.RData")

#Break thing down into components ##### 

thing_D1<-thing[1:(m1*m2), 1:(m1*m2)]
thing_G<-thing[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
thing_F<-thing[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
thing_D2<-thing[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

#Unpack matrix of contributions into pieces #####

	
	#D1#
	plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in 	array
	Plop=outer(1:m1,1:m2,plop); 

	Kvals_thing_D1=array(0,c(m1,m2,m1,m2));  

	for(i in 1:m1){
		for(j in 1:m2){
		for(k in 1:m1){
				kvals= thing_D1[Plop[k,1:m2],Plop[i,j]]
				Kvals_thing_D1[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_thing_D2=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= thing_D2[Plop[k,1:m4],Plop[i,j]]
				Kvals_thing_D2[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_thing_F=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=thing_F[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_thing_F[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_thing_G=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=thing_G[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_thing_G[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


# Sum up total contributions #####
total.thing_D1<-apply(Kvals_thing_D1, c(3,4), sum);
total.thing_D2<-apply(Kvals_thing_D2, c(3,4), sum);
total.thing_F<-apply(Kvals_thing_F, c(3,4), sum);
total.thing_G<-apply(Kvals_thing_G, c(3,4), sum);
	
save(total.thing_D1, total.thing_D2, total.thing_F, total.thing_G, 
     file="./Overall/total.thing.RData")








