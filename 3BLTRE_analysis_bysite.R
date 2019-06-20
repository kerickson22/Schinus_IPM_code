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
D_BC<-A_BC-A_avg #Difference matrix
rm(A_BC)
C_BC<-D_BC*sens_avg #Matrix of contributions
rm(D_BC)
save(C_BC, file="./BC/C_BC.RData")


#Cape Canaveral
load("./CC/A_CC.RData")
D_CC<-A_CC-A_avg #Difference matrix
rm(A_CC)
C_CC<-D_CC*sens_avg #Matrix of contributions
rm(D_CC)
save(C_CC, file="./CC/C_CC.RData")

#Chekika
load("./C/A_C.RData")
D_C<-A_C-A_avg #Difference matrix
rm(A_C)
C_C<-D_C*sens_avg #Matrix of contributions
rm(D_C)
save(C_C, file="./C/C_C.RData")

#Fort Pierce
load("./FP/A_FP.RData")
D_FP<-A_FP-A_avg #Difference matrix
rm(A_FP)
C_FP<-D_FP*sens_avg #Matrix of contributions
rm(D_FP)
save(C_FP, file="./FP/C_FP.RData")

#Punta Gorda
load("./PG/A_PG.RData")
D_PG<-A_PG-A_avg #Difference matrix
rm(A_PG)
C_PG<-D_PG*sens_avg #Matrix of contributions
rm(D_PG)
save(C_PG, file="./PG/C_PG.RData")

#Wild Turkey
load("./WT/A_WT.RData")
D_WT<-A_WT-A_avg #Difference matrix
rm(A_WT)
C_WT<-D_WT*sens_avg #Matrix of contributions
rm(D_WT)
save(C_WT, file="./WT/C_WT.RData")

# Decompose contribution matrices into components #####
decompose = function (mat) {
  #Break the elasticity matrix back into its component parts: 
  mat_D1<-mat[1:(m1*m2), 1:(m1*m2)]
  mat_G<-mat[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
  mat_F<-mat[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
  mat_D2<-mat[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]
  return(list(mat_D1 = mat_D1, mat_G = mat_G, mat_F = mat_F, mat_D2 = mat_D2))
}

load("./BC/C_BC.RData")
load("./CC/C_CC.RData")
load("./C/C_C.RData")
load("./FP/C_FP.RData")
load("./PG/C_PG.RData")
load("./WT/C_WT.RData")

#Big Cypress
thing <- decompose(C_BC)
saveRDS(thing$mat_D1, file="./BC/C_D1_BC.rds")
saveRDS(thing$mat_G, file="./BC/C_G_BC.rds")
saveRDS(thing$mat_F, file="./BC/C_F_BC.rds")
saveRDS(thing$mat_D2, file="./BC/C_D2_BC.rds")
#Cape Canaveral
thing <- decompose(C_CC)
saveRDS(thing$mat_D1, file="./CC/C_D1_CC.rds")
saveRDS(thing$mat_G, file="./CC/C_G_CC.rds")
saveRDS(thing$mat_F, file="./CC/C_F_CC.rds")
saveRDS(thing$mat_D2, file="./CC/C_D2_CC.rds")
#Chekika
thing <- decompose(C_C)
saveRDS(thing$mat_D1, file="./C/C_D1_C.rds")
saveRDS(thing$mat_G, file="./C/C_G_C.rds")
saveRDS(thing$mat_F, file="./C/C_F_C.rds")
saveRDS(thing$mat_D2, file="./C/C_D2_C.rds")
#Fort Pierce
thing <- decompose(C_FP)
saveRDS(thing$mat_D1, file="./FP/C_D1_FP.rds")
saveRDS(thing$mat_G, file="./FP/C_G_FP.rds")
saveRDS(thing$mat_F, file="./FP/C_F_FP.rds")
saveRDS(thing$mat_D2, file="./FP/C_D2_FP.rds")
#Punta Gorda
thing <- decompose(C_PG)
saveRDS(thing$mat_D1, file="./PG/C_D1_PG.rds")
saveRDS(thing$mat_G, file="./PG/C_G_PG.rds")
saveRDS(thing$mat_F, file="./PG/C_F_PG.rds")
saveRDS(thing$mat_D2, file="./PG/C_D2_PG.rds")
#Wild Turkey
thing <- decompose(C_WT)
saveRDS(thing$mat_D1, file="./WT/C_D1_WT.rds")
saveRDS(thing$mat_G, file="./WT/C_G_WT.rds")
saveRDS(thing$mat_F, file="./WT/C_F_WT.rds")
saveRDS(thing$mat_D2, file="./WT/C_D2_WT.rds")

rm(C_BC, C_C, C_CC, C_FP, C_PG, C_WT)

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
C_D1_BC <-readRDS("./BC/C_D1_BC.rds")
Kvals_contrib_D1_BC <- toKvals_D1(C_D1_BC)
save(Kvals_contrib_D1_BC, file="./BC/Kvals_contrib_D1_BC.RData")
#Cape Canaveral
C_D1_CC <-readRDS("./CC/C_D1_CC.rds")
Kvals_contrib_D1_CC <- toKvals_D1(C_D1_CC)
save(Kvals_contrib_D1_CC, file="./CC/Kvals_contrib_D1_CC.RData")
#Chekika
C_D1_C <- readRDS("./C/C_D1_C.rds")
Kvals_contrib_D1_C <- toKvals_D1(C_D1_C)
save(Kvals_contrib_D1_C, file="./C/Kvals_contrib_D1_C.RData")
#Fort Pierce
C_D1_FP <- readRDS("./FP/C_D1_FP.rds")
Kvals_contrib_D1_FP <- toKvals_D1(C_D1_FP)
save(Kvals_contrib_D1_FP, file="./FP/Kvals_contrib_D1_FP.RData")
#Punta Gorda
C_D1_PG <- readRDS("./PG/C_D1_PG.rds")
Kvals_contrib_D1_PG <- toKvals_D1(C_D1_PG)
save(Kvals_contrib_D1_PG, file="./PG/Kvals_contrib_D1_PG.RData")
#Wild Turkey
C_D1_WT <- readRDS("./WT/C_D1_WT.rds")
Kvals_contrib_D1_WT <- toKvals_D1(C_D1_WT)
save(Kvals_contrib_D1_WT, file="./WT/Kvals_contrib_D1_WT.RData")





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
C_D2_BC <-readRDS("./BC/C_D2_BC.rds")
Kvals_contrib_D2_BC <- toKvals_D2(C_D2_BC)
save(Kvals_contrib_D2_BC, file="./BC/Kvals_contrib_D2_BC.RData")
#Cape Canaveral
C_D2_CC <-readRDS("./CC/C_D2_CC.rds")
Kvals_contrib_D2_CC <- toKvals_D2(C_D2_CC)
save(Kvals_contrib_D2_CC, file="./CC/Kvals_contrib_D2_CC.RData")
#Chekika
C_D2_C <- readRDS("./C/C_D2_C.rds")
Kvals_contrib_D2_C <- toKvals_D2(C_D2_C)
save(Kvals_contrib_D2_C, file="./C/Kvals_contrib_D2_C.RData")
#Fort Pierce
C_D2_FP <- readRDS("./FP/C_D2_FP.rds")
Kvals_contrib_D2_FP <- toKvals_D2(C_D2_FP)
save(Kvals_contrib_D2_FP, file="./FP/Kvals_contrib_D2_FP.RData")
#Punta Gorda
C_D2_PG <- readRDS("./PG/C_D2_PG.rds")
Kvals_contrib_D2_PG <- toKvals_D2(C_D2_PG)
save(Kvals_contrib_D2_PG, file="./PG/Kvals_contrib_D2_PG.RData")
#Wild Turkey
C_D2_WT <- readRDS("./WT/C_D2_WT.rds")
Kvals_contrib_D2_WT <- toKvals_D2(C_D2_WT)
save(Kvals_contrib_D2_WT, file="./WT/Kvals_contrib_D2_WT.RData")

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
C_F_BC <-readRDS("./BC/C_F_BC.rds")
Kvals_contrib_F_BC <- toKvals_F(C_F_BC)
save(Kvals_contrib_F_BC, file="./BC/Kvals_contrib_F_BC.RData")
#Cape Canaveral
C_F_CC <-readRDS("./CC/C_F_CC.rds")
Kvals_contrib_F_CC <- toKvals_F(C_F_CC)
save(Kvals_contrib_F_CC, file="./CC/Kvals_contrib_F_CC.RData")
#Chekika
C_F_C <- readRDS("./C/C_F_C.rds")
Kvals_contrib_F_C <- toKvals_F(C_F_C)
save(Kvals_contrib_F_C, file="./C/Kvals_contrib_F_C.RData")
#Fort Pierce
C_F_FP <- readRDS("./FP/C_F_FP.rds")
Kvals_contrib_F_FP <- toKvals_F(C_F_FP)
save(Kvals_contrib_F_FP, file="./FP/Kvals_contrib_F_FP.RData")
#Punta Gorda
C_F_PG <- readRDS("./PG/C_F_PG.rds")
Kvals_contrib_F_PG <- toKvals_F(C_F_PG)
save(Kvals_contrib_F_PG, file="./PG/Kvals_contrib_F_PG.RData")
#Wild Turkey
C_F_WT <- readRDS("./WT/C_F_WT.rds")
Kvals_contrib_F_WT <- toKvals_F(C_F_WT)
save(Kvals_contrib_F_WT, file="./WT/Kvals_contrib_F_WT.RData")


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
C_G_BC <-readRDS("./BC/C_G_BC.rds")
Kvals_contrib_G_BC <- toKvals_G(C_G_BC)
save(Kvals_contrib_G_BC, file="./BC/Kvals_contrib_G_BC.RData")
#Cape Canaveral
C_G_CC <-readRDS("./CC/C_G_CC.rds")
Kvals_contrib_G_CC <- toKvals_G(C_G_CC)
save(Kvals_contrib_G_CC, file="./CC/Kvals_contrib_G_CC.RData")
#Chekika
C_G_C <- readRDS("./C/C_G_C.rds")
Kvals_contrib_G_C <- toKvals_G(C_G_C)
save(Kvals_contrib_G_C, file="./C/Kvals_contrib_G_C.RData")
#Fort Pierce
C_G_FP <- readRDS("./FP/C_G_FP.rds")
Kvals_contrib_G_FP <- toKvals_G(C_G_FP)
save(Kvals_contrib_G_FP, file="./FP/Kvals_contrib_G_FP.RData")
#Punta Gorda
C_G_PG <- readRDS("./PG/C_G_PG.rds")
Kvals_contrib_G_PG <- toKvals_G(C_G_PG)
save(Kvals_contrib_G_PG, file="./PG/Kvals_contrib_G_PG.RData")
#Wild Turkey
C_G_WT <- readRDS("./WT/C_G_WT.rds")
Kvals_contrib_G_WT <- toKvals_G(C_G_WT)
save(Kvals_contrib_G_WT, file="./WT/Kvals_contrib_G_WT.RData")


# Convert to totals: #####
#Big Cypress
total.contrib_D1_BC<-apply(Kvals_contrib_D1_BC, c(3,4), sum);
total.contrib_D2_BC<-apply(Kvals_contrib_D2_BC, c(3,4), sum);
total.contrib_F_BC<-apply(Kvals_contrib_F_BC, c(3,4), sum);
total.contrib_G_BC<-apply(Kvals_contrib_G_BC, c(3,4), sum);
save(total.contrib_D1_BC, total.contrib_D2_BC, total.contrib_F_BC, total.contrib_G_BC, 
     file="./BC/total.contrib_BC.RData")

#Cape Canaveral
total.contrib_D1_CC<-apply(Kvals_contrib_D1_CC, c(3,4), sum);
total.contrib_D2_CC<-apply(Kvals_contrib_D2_CC, c(3,4), sum);
total.contrib_F_CC<-apply(Kvals_contrib_F_CC, c(3,4), sum);
total.contrib_G_CC<-apply(Kvals_contrib_G_CC, c(3,4), sum);
save(total.contrib_D1_CC, total.contrib_D2_CC, total.contrib_F_CC, total.contrib_G_CC, 
     file="./CC/total.contrib_CC.RData")

#Chekika
total.contrib_D1_C<-apply(Kvals_contrib_D1_C, c(3,4), sum);
total.contrib_D2_C<-apply(Kvals_contrib_D2_C, c(3,4), sum);
total.contrib_F_C<-apply(Kvals_contrib_F_C, c(3,4), sum);
total.contrib_G_C<-apply(Kvals_contrib_G_C, c(3,4), sum);
save(total.contrib_D1_C, total.contrib_D2_C, total.contrib_F_C, total.contrib_G_C, 
     file="./C/total.contrib_C.RData")

#Fort Pierce
total.contrib_D1_FP<-apply(Kvals_contrib_D1_FP, c(3,4), sum);
total.contrib_D2_FP<-apply(Kvals_contrib_D2_FP, c(3,4), sum);
total.contrib_F_FP<-apply(Kvals_contrib_F_FP, c(3,4), sum);
total.contrib_G_FP<-apply(Kvals_contrib_G_FP, c(3,4), sum);
save(total.contrib_D1_FP, total.contrib_D2_FP, total.contrib_F_FP, total.contrib_G_FP, 
     file="./FP/total.contrib_FP.RData")

#Punta Gorda
total.contrib_D1_PG<-apply(Kvals_contrib_D1_PG, c(3,4), sum);
total.contrib_D2_PG<-apply(Kvals_contrib_D2_PG, c(3,4), sum);
total.contrib_F_PG<-apply(Kvals_contrib_F_PG, c(3,4), sum);
total.contrib_G_PG<-apply(Kvals_contrib_G_PG, c(3,4), sum);
save(total.contrib_D1_PG, total.contrib_D2_PG, total.contrib_F_PG, total.contrib_G_PG, 
     file="./PG/total.contrib_PG.RData")

#Wild Turkey
total.contrib_D1_WT<-apply(Kvals_contrib_D1_WT, c(3,4), sum);
total.contrib_D2_WT<-apply(Kvals_contrib_D2_WT, c(3,4), sum);
total.contrib_F_WT<-apply(Kvals_contrib_F_WT, c(3,4), sum);
total.contrib_G_WT<-apply(Kvals_contrib_G_WT, c(3,4), sum);
save(total.contrib_D1_WT, total.contrib_D2_WT, total.contrib_F_WT, total.contrib_G_WT, 
     file="./WT/total.contrib_WT.RData")


#Collapse contributions by site: #####


 BC<-c( 
 sum(rowSums(total.contrib_D1_BC)), 
 sum(rowSums(total.contrib_G_BC)), 
 sum(rowSums(total.contrib_F_BC)), 
 sum(rowSums(total.contrib_D2_BC)) )
 save(BC, file="./BC/BC.RData")
 
 
 CC<-c( 
   sum(rowSums(total.contrib_D1_CC)), 
   sum(rowSums(total.contrib_G_CC)), 
   sum(rowSums(total.contrib_F_CC)), 
   sum(rowSums(total.contrib_D2_CC)) )
 save(CC, file="./CC/CC.RData")
 
 C<-c( 
   sum(rowSums(total.contrib_D1_C)), 
   sum(rowSums(total.contrib_G_C)), 
   sum(rowSums(total.contrib_F_C)), 
   sum(rowSums(total.contrib_D2_C)) )
 save(C, file="./C/C.RData")
 
 FP<-c( 
   sum(rowSums(total.contrib_D1_FP)), 
   sum(rowSums(total.contrib_G_FP)), 
   sum(rowSums(total.contrib_F_FP)), 
   sum(rowSums(total.contrib_D2_FP)) )
 save(FP, file="./FP/FP.RData")
 
 PG<-c( 
   sum(rowSums(total.contrib_D1_PG)), 
   sum(rowSums(total.contrib_G_PG)), 
   sum(rowSums(total.contrib_F_PG)), 
   sum(rowSums(total.contrib_D2_PG)) )
 save(PG, file="./PG/PG.RData")
 
 
 WT<-c( 
   sum(rowSums(total.contrib_D1_WT)), 
   sum(rowSums(total.contrib_G_WT)), 
   sum(rowSums(total.contrib_F_WT)), 
   sum(rowSums(total.contrib_D2_WT)) )
 save(WT, file="./WT/WT.RData")
# Calculate matrices of coefficients of variation #####

	
D1_BC<-readRDS("./BC/D1_BC.rds")
D1_CC <- readRDS("./CC/D1_CC.rds")
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
		cv_D1[i,j]<-CV(c(D1_BC[i,j], D1_CC[i,j], D1_C[i,j], D1_FP[i,j], 
		                 D1_PG[i,j], D1_WT[i,j]))
		if(is.na(cv_D1[i,j])) {cv_D1[i,j]<-0}
	}
}


D1_avg=(D1_BC + D1_CC + D1_C + D1_FP + D1_PG + D1_WT)/6
rm(D1_BC, D1_CC, D1_C, D1_FP, D1_PG, D1_WT)

#G
G_BC<-readRDS("./BC/G_BC.rds")
G_CC<-readRDS("./CC/G_CC.rds")
G_C <- readRDS("./C/G_C.rds")
G_FP <- readRDS("./FP/G_FP.rds")
G_PG <- readRDS("./PG/G_PG.rds")
G_WT <- readRDS("./WT/G_WT.rds")

cv_G<-matrix(nrow=(m3*m4), ncol=(m1*m2))


for(i in 1:(m3*m4)) {
	for(j in 1:(m1*m2)) {
		cv_G[i,j]<-CV(c(G_BC[i,j], G_CC[i,j], G_C[i,j], G_FP[i,j],
		                G_PG[i,j], G_WT[i,j]))
		if(is.na(cv_G[i,j])) {cv_G[i,j]<-0}
	}
}

G_avg<-(G_BC+G_CC+G_C + G_FP + G_PG + G_WT)/6

rm(G_BC, G_CC, G_C, G_FP, G_PG, G_WT)

#F
F_BC <- readRDS("./BC/F_BC.rds")
F_CC <- readRDS("./CC/F_CC.rds")
F_C  <- readRDS("./C/F_C.rds")
F_FP <- readRDS("./FP/F_FP.rds")
F_PG <- readRDS("./PG/F_PG.rds")
F_WT <- readRDS("./WT/F_WT.rds")

cv_F<-matrix(nrow=(m1*m2), ncol=(m3*m4))
for(i in 1:(m1*m2)) {
	for(j in 1:(m3*m4)) {
		cv_F[i,j]<-CV(c(F_BC[i,j], F_CC[i,j], F_C[i,j], F_FP[i,j],
		              F_PG[i,j], F_WT[i,j]))
		if(is.na(cv_F[i,j])) {cv_F[i,j]<-0}
	}
}

F_avg<-(F_BC + F_CC + F_C + F_FP + F_PG + F_WT)/6

rm(F_BC, F_CC, F_C, F_FP, F_PG, F_WT)

#D2
D2_BC <- readRDS("./BC/D2_BC.rds")
D2_CC <- readRDS("./CC/D2_CC.rds")
D2_C  <- readRDS("./C/D2_C.rds")
D2_FP <- readRDS("./FP/D2_FP.rds")
D2_PG <- readRDS("./PG/D2_PG.rds")
D2_WT <- readRDS("./WT/D2_WT.rds")


cv_D2<-matrix(nrow=(m3*m4), ncol=(m3*m4))
for(i in 1:(m3*m4)) {
	for(j in 1:(m3*m4)) {
		cv_D2[i,j]<-CV(c(D2_BC[i,j], D2_CC[i,j], D2_C[i,j], D2_FP,
		                 D2_PG, D2_WT))
		if(is.na(cv_D2[i,j])) {cv_D2[i,j]<-0}
	}
}

save(cv_D1, file="cv_D1.RData")
save(cv_D2, file="cv_D2.RData")
save(cv_F_, file="cv_F.RData")
save(cv_G, file="cv_G.RData")

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








