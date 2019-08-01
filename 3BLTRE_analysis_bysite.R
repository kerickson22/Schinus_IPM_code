# Perform a Life Table Response Experiment analysis using the IPM constructed from 
# demography_sites.R. Makes use of code from Ellner and Rees. 
# Modified by Kelley D. Erickson and Carol C. Horvitz


rm(list=ls())

#LTRE
m1=10
m2=m1+1
m3a=250
m4a=50

m3b = 120
m4b = 100
m3=m3a + m3b
m4=m4a + m4b
tol=1.e-8; 




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


#Calculate A_avg and associated quantities#####
#Computational limit: can only load two matrices at a time 
 load("./BC/A_BC.RData")
 load("./CC/A_CC.RData")
 A_temp <- A_BC + A_CC
 rm(A_BC, A_CC)
 gc()
 load("./C/A_C.RData")
 A_temp <- A_temp + A_C
 rm(A_C)
 gc()

 load("./FP/A_FP.RData")
 A_temp <- A_temp + A_FP
 rm(A_FP)
 gc()
 
 load("./PG/A_PG.RData")
 A_temp <- A_temp + A_PG
 rm(A_PG)
 gc()
 
 load("./WT/A_WT.RData")
 A_temp <- A_temp + A_WT
 rm(A_WT)
 gc()
 
 A_avg <- (A_temp)/6
 
 save(A_avg, file="./Average/temp/A_avg.RData")

 
#Now calculate the sensitivity of the average matrix: 

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
 thing_elas <- elasticity(thing_avg$v, thing_avg$stable.dist, A_avg, 
                          thing_avg$lam.stable)
 sens_avg <- thing_elas$sens #fix this! (it was elas)
 save(sens_avg, file="./Average/temp/sens_avg.RData")

# Calculate difference matrices ##### 
# Big Cypress 
#(a) Overall
load("./BC/A_BC.RData")
load("./Overall/A_overall.RData")
D_BC_pool <- A_BC -A_overall 
rm(A_BC, A_overall)
load("./Overall/sens_pool.RData")
C_BC_pool <- D_BC_pool*sens_pool
rm(D_BC_pool, sens_pool)
save(C_BC_pool, file="./BC/C_BC_pool.RData")
rm(C_BC_pool)

#(b) average matrix 
load("./BC/A_BC.RData")
load("./Average/temp/A_avg.RData")
D_BC_avg <-A_BC-A_avg #Difference matrix
rm(A_BC, A_avg)
load("./Average/temp/sens_avg.RData")
C_BC_avg<-D_BC_avg*sens_avg #Matrix of contributions
rm(D_BC_avg, sens_avg)
save(C_BC_avg, file="./BC/C_BC_avg.RData")
rm(C_BC_avg)



# Cape Canaveral 
#(a) Overall
load("./CC/A_CC.RData")
load("./Overall/A_overall.RData")
D_CC_pool <- A_CC -A_overall 
rm(A_CC, A_overall)
load("./Overall/sens_pool.RData")
C_CC_pool <- D_CC_pool*sens_pool
rm(D_CC_pool, sens_pool)
save(C_CC_pool, file="./CC/C_CC_pool.RData")
rm(C_CC_pool)

#(b) average matrix 
load("./CC/A_CC.RData")
load("./Average/temp/A_avg.RData")
D_CC_avg <-A_CC-A_avg #Difference matrix
rm(A_CC, A_avg)
load("./Average/temp/sens_avg.RData")
C_CC_avg<-D_CC_avg*sens_avg #Matrix of contributions
rm(D_CC_avg, sens_avg)
save(C_CC_avg, file="./CC/C_CC_avg.RData")
rm(C_CC_avg)


#Chekika 
#(a) Overall
load("./C/A_C.RData")
load("./Overall/A_overall.RData")
D_C_pool <- A_C -A_overall 
rm(A_C, A_overall)
load("./Overall/sens_pool.RData")
C_C_pool <- D_C_pool*sens_pool
rm(D_C_pool, sens_pool)
save(C_C_pool, file="./C/C_C_pool.RData")
rm(C_C_pool)

#(b) average matrix 
load("./C/A_C.RData")
load("./Average/temp/A_avg.RData")
D_C_avg <-A_C-A_avg #Difference matrix
rm(A_C, A_avg)
load("./Average/temp/sens_avg.RData")
C_C_avg<-D_C_avg*sens_avg #Matrix of contributions
rm(D_C_avg, sens_avg)
save(C_C_avg, file="./C/C_C_avg.RData")
rm(C_C_avg)

#Fort Pierce
#(a) Overall
load("./FP/A_FP.RData")
load("./Overall/A_overall.RData")
D_FP_pool <- A_FP -A_overall 
rm(A_FP, A_overall)
load("./Overall/sens_pool.RData")
C_FP_pool <- D_FP_pool*sens_pool
rm(D_FP_pool, sens_pool)
save(C_FP_pool, file="./FP/C_FP_pool.RData")
rm(C_FP_pool)

#(b) average matrix 
load("./FP/A_FP.RData")
load("./Average/temp/A_avg.RData")
D_FP_avg <-A_FP-A_avg #Difference matrix
rm(A_FP, A_avg)
load("./Average/temp/sens_avg.RData")
C_FP_avg<-D_FP_avg*sens_avg #Matrix of contributions
rm(D_FP_avg, sens_avg)
save(C_FP_avg, file="./FP/C_FP_avg.RData")
rm(C_FP_avg)

#Punta Gorda
#(a) Overall
load("./PG/A_PG.RData")
load("./Overall/A_overall.RData")
D_PG_pool <- A_PG -A_overall 
rm(A_PG, A_overall)
load("./Overall/sens_pool.RData")
C_PG_pool <- D_PG_pool*sens_pool
rm(D_PG_pool, sens_pool)
save(C_PG_pool, file="./PG/C_PG_pool.RData")
rm(C_PG_pool)

#(b) average matrix 
load("./PG/A_PG.RData")
load("./Average/temp/A_avg.RData")
D_PG_avg <-A_PG-A_avg #Difference matrix
rm(A_PG, A_avg)
load("./Average/temp/sens_avg.RData")
C_PG_avg<-D_PG_avg*sens_avg #Matrix of contributions
rm(D_PG_avg, sens_avg)
save(C_PG_avg, file="./PG/C_PG_avg.RData")
rm(C_PG_avg)

#Wild Turkey
#(a) Overall
load("./WT/A_WT.RData")
load("./Overall/A_overall.RData")
D_WT_pool <- A_WT -A_overall 
rm(A_WT, A_overall)
load("./Overall/sens_pool.RData")
C_WT_pool <- D_WT_pool*sens_pool
rm(D_WT_pool, sens_pool)
save(C_WT_pool, file="./WT/C_WT_pool.RData")
rm(C_WT_pool)

#(b) average matrix 
load("./WT/A_WT.RData")
load("./Average/temp/A_avg.RData")
D_WT_avg <-A_WT-A_avg #Difference matrix
rm(A_WT, A_avg)
load("./Average/temp/sens_avg.RData")
C_WT_avg<-D_WT_avg*sens_avg #Matrix of contributions
rm(D_WT_avg, sens_avg)
save(C_WT_avg, file="./WT/C_WT_avg.RData")
rm(C_WT_avg)
# Decompose contribution matrices into components #####
decompose = function (mat) {
  #Break the elasticity matrix back into its component parts: 
  mat_D1 <- mat[1:(m1*m2), 1:(m1*m2)]
  mat_FA <- mat[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3a*m4a)]
  mat_FB <- mat[1:(m1*m2), 
                ((m1*m2 + m3a*m4a)+1):(m1*m2 +m3a*m4a + m3b*m4b)]
  mat_GA <- mat[((m1*m2)+1): (m1*m2 + m3a*m4a), 1:(m1*m2)]
  mat_GB <- mat[((m1*m2 + m3a*m4a)+1):(m1*m2 + m3a*m4a + m3b*m4b),
                1:(m1*m2)]
  mat_D2AA <- mat[((m1*m2)+1):(m1*m2 + m3a*m4a), 
                  ((m1*m2)+1):(m1*m2 + m3a*m4a)]
  mat_D2BA <- mat[((m1*m2)+1):(m1*m2 + m3a*m4a), 
                  ((m1*m2 + m3a*m4a)+1):(m1*m2 + m3a*m4a + m3b*m4b)]
  mat_D2AB <- mat[((m1*m2 + m3a*m4a)+1):(m1*m2 + m3a*m4a+m3b*m4b),
                  ((m1*m2)+1):(m1*m2 + m3a*m4a)]
  mat_D2BB <- mat[((m1*m2 + m3a*m4a)+1):(m1*m2 + m3a*m4a + m3b*m4b),
                    ((m1*m2 + m3a*m4a)+1):(m1*m2 + m3a*m4a + m3b*m4b)]
  return(list(mat_D1 = mat_D1, mat_GA = mat_GA, mat_GB = mat_GB,
              mat_FA = mat_FA, mat_FB = mat_FB,  mat_D2AA = mat_D2AA,
              mat_D2AB = mat_D2AB, mat_D2BA = mat_D2BA, 
              mat_D2BB = mat_D2BB))
}


#Big Cypress
load("./BC/C_BC_avg.RData")
thing <- decompose(C_BC_avg)
saveRDS(thing$mat_D1, file="./BC/C_D1_BC_avg.rds")
saveRDS(thing$mat_GA, file="./BC/C_GA_BC_avg.rds")
saveRDS(thing$mat_GB, file="./BC/C_GB_BC_avg.rds")
saveRDS(thing$mat_FA, file="./BC/C_FA_BC_avg.rds")
saveRDS(thing$mat_FB, file="./BC/C_FB_BC_avg.rds")
saveRDS(thing$mat_D2AA, file="./BC/C_D2AA_BC_avg.rds")
saveRDS(thing$mat_D2AB, file="./BC/C_D2AB_BC_avg.rds")
saveRDS(thing$mat_D2BA, file="./BC/C_D2BA_BC_avg.rds")
saveRDS(thing$mat_D2BB, file="./BC/C_D2BB_BC_avg.rds")
rm(thing, C_BC_avg)

#(b) Pool
load("./BC/C_BC_pool.RData")
thing <- decompose(C_BC_pool)
saveRDS(thing$mat_D1, file="./BC/C_D1_BC_pool.rds")
saveRDS(thing$mat_GA, file="./BC/C_GA_BC_pool.rds")
saveRDS(thing$mat_GB, file="./BC/C_GB_BC_pool.rds")
saveRDS(thing$mat_FA, file="./BC/C_FA_BC_pool.rds")
saveRDS(thing$mat_FB, file="./BC/C_FB_BC_pool.rds")
saveRDS(thing$mat_D2AA, file="./BC/C_D2AA_BC_pool.rds")
saveRDS(thing$mat_D2AB, file="./BC/C_D2AB_BC_pool.rds")
saveRDS(thing$mat_D2BA, file="./BC/C_D2BA_BC_pool.rds")
saveRDS(thing$mat_D2BB, file="./BC/C_D2BB_BC_pool.rds")
rm(thing, C_BC_pool)


#Cape Canaveral
load("./CC/C_CC_avg.RData")
thing <- decompose(C_CC_avg)
saveRDS(thing$mat_D1, file="./CC/C_D1_CC_avg.rds")
saveRDS(thing$mat_GA, file="./CC/C_GA_CC_avg.rds")
saveRDS(thing$mat_GB, file="./CC/C_GB_CC_avg.rds")
saveRDS(thing$mat_FA, file="./CC/C_FA_CC_avg.rds")
saveRDS(thing$mat_FB, file="./CC/C_FB_CC_avg.rds")
saveRDS(thing$mat_D2AA, file="./CC/C_D2AA_CC_avg.rds")
saveRDS(thing$mat_D2AB, file="./CC/C_D2AB_CC_avg.rds")
saveRDS(thing$mat_D2BA, file="./CC/C_D2BA_CC_avg.rds")
saveRDS(thing$mat_D2BB, file="./CC/C_D2BB_CC_avg.rds")
rm(thing, C_CC_avg)

#(b) Pool
load("./CC/C_CC_pool.RData")
thing <- decompose(C_CC_pool)
saveRDS(thing$mat_D1, file="./CC/C_D1_CC_pool.rds")
saveRDS(thing$mat_GA, file="./CC/C_GA_CC_pool.rds")
saveRDS(thing$mat_GB, file="./CC/C_GB_CC_pool.rds")
saveRDS(thing$mat_FA, file="./CC/C_FA_CC_pool.rds")
saveRDS(thing$mat_FB, file="./CC/C_FB_CC_pool.rds")
saveRDS(thing$mat_D2AA, file="./CC/C_D2AA_CC_pool.rds")
saveRDS(thing$mat_D2AB, file="./CC/C_D2AB_CC_pool.rds")
saveRDS(thing$mat_D2BA, file="./CC/C_D2BA_CC_pool.rds")
saveRDS(thing$mat_D2BB, file="./CC/C_D2BB_CC_pool.rds")
rm(thing, C_CC_pool)

#Chekika
load("./C/C_C_avg.RData")
thing <- decompose(C_C_avg)
saveRDS(thing$mat_D1, file="./C/C_D1_C_avg.rds")
saveRDS(thing$mat_GA, file="./C/C_GA_C_avg.rds")
saveRDS(thing$mat_GB, file="./C/C_GB_C_avg.rds")
saveRDS(thing$mat_FA, file="./C/C_FA_C_avg.rds")
saveRDS(thing$mat_FB, file="./C/C_FB_C_avg.rds")
saveRDS(thing$mat_D2AA, file="./C/C_D2AA_C_avg.rds")
saveRDS(thing$mat_D2AB, file="./C/C_D2AB_C_avg.rds")
saveRDS(thing$mat_D2BA, file="./C/C_D2BA_C_avg.rds")
saveRDS(thing$mat_D2BB, file="./C/C_D2BB_C_avg.rds")
rm(thing, C_C_avg)

#(b) Pool
load("./C/C_C_pool.RData")
thing <- decompose(C_C_pool)
saveRDS(thing$mat_D1, file="./C/C_D1_C_pool.rds")
saveRDS(thing$mat_GA, file="./C/C_GA_C_pool.rds")
saveRDS(thing$mat_GB, file="./C/C_GB_C_pool.rds")
saveRDS(thing$mat_FA, file="./C/C_FA_C_pool.rds")
saveRDS(thing$mat_FB, file="./C/C_FB_C_pool.rds")
saveRDS(thing$mat_D2AA, file="./C/C_D2AA_C_pool.rds")
saveRDS(thing$mat_D2AB, file="./C/C_D2AB_C_pool.rds")
saveRDS(thing$mat_D2BA, file="./C/C_D2BA_C_pool.rds")
saveRDS(thing$mat_D2BB, file="./C/C_D2BB_C_pool.rds")
rm(thing, C_C_pool)

#Fort Pierce
load("./FP/C_FP_avg.RData")
thing <- decompose(C_FP_avg)
saveRDS(thing$mat_D1, file="./FP/C_D1_FP_avg.rds")
saveRDS(thing$mat_GA, file="./FP/C_GA_FP_avg.rds")
saveRDS(thing$mat_GB, file="./FP/C_GB_FP_avg.rds")
saveRDS(thing$mat_FA, file="./FP/C_FA_FP_avg.rds")
saveRDS(thing$mat_FB, file="./FP/C_FB_FP_avg.rds")
saveRDS(thing$mat_D2AA, file="./FP/C_D2AA_FP_avg.rds")
saveRDS(thing$mat_D2AB, file="./FP/C_D2AB_FP_avg.rds")
saveRDS(thing$mat_D2BA, file="./FP/C_D2BA_FP_avg.rds")
saveRDS(thing$mat_D2BB, file="./FP/C_D2BB_FP_avg.rds")
rm(thing, C_FP_avg)

#(b) Pool
load("./FP/C_FP_pool.RData")
thing <- decompose(C_FP_pool)
saveRDS(thing$mat_D1, file="./FP/C_D1_FP_pool.rds")
saveRDS(thing$mat_GA, file="./FP/C_GA_FP_pool.rds")
saveRDS(thing$mat_GB, file="./FP/C_GB_FP_pool.rds")
saveRDS(thing$mat_FA, file="./FP/C_FA_FP_pool.rds")
saveRDS(thing$mat_FB, file="./FP/C_FB_FP_pool.rds")
saveRDS(thing$mat_D2AA, file="./FP/C_D2AA_FP_pool.rds")
saveRDS(thing$mat_D2AB, file="./FP/C_D2AB_FP_pool.rds")
saveRDS(thing$mat_D2BA, file="./FP/C_D2BA_FP_pool.rds")
saveRDS(thing$mat_D2BB, file="./FP/C_D2BB_FP_pool.rds")
rm(thing, C_FP_pool)

#Punta Gorda
load("./PG/C_PG_avg.RData")
thing <- decompose(C_PG_avg)
saveRDS(thing$mat_D1, file="./PG/C_D1_PG_avg.rds")
saveRDS(thing$mat_GA, file="./PG/C_GA_PG_avg.rds")
saveRDS(thing$mat_GB, file="./PG/C_GB_PG_avg.rds")
saveRDS(thing$mat_FA, file="./PG/C_FA_PG_avg.rds")
saveRDS(thing$mat_FB, file="./PG/C_FB_PG_avg.rds")
saveRDS(thing$mat_D2AA, file="./PG/C_D2AA_PG_avg.rds")
saveRDS(thing$mat_D2AB, file="./PG/C_D2AB_PG_avg.rds")
saveRDS(thing$mat_D2BA, file="./PG/C_D2BA_PG_avg.rds")
saveRDS(thing$mat_D2BB, file="./PG/C_D2BB_PG_avg.rds")
rm(thing, C_PG_avg)

#(b) Pool
load("./PG/C_PG_pool.RData")
thing <- decompose(C_PG_pool)
saveRDS(thing$mat_D1, file="./PG/C_D1_PG_pool.rds")
saveRDS(thing$mat_GA, file="./PG/C_GA_PG_pool.rds")
saveRDS(thing$mat_GB, file="./PG/C_GB_PG_pool.rds")
saveRDS(thing$mat_FA, file="./PG/C_FA_PG_pool.rds")
saveRDS(thing$mat_FB, file="./PG/C_FB_PG_pool.rds")
saveRDS(thing$mat_D2AA, file="./PG/C_D2AA_PG_pool.rds")
saveRDS(thing$mat_D2AB, file="./PG/C_D2AB_PG_pool.rds")
saveRDS(thing$mat_D2BA, file="./PG/C_D2BA_PG_pool.rds")
saveRDS(thing$mat_D2BB, file="./PG/C_D2BB_PG_pool.rds")
rm(thing, C_PG_pool)

#Wild Turkey
load("./WT/C_WT_avg.RData")
thing <- decompose(C_WT_avg)
saveRDS(thing$mat_D1, file="./WT/C_D1_WT_avg.rds")
saveRDS(thing$mat_GA, file="./WT/C_GA_WT_avg.rds")
saveRDS(thing$mat_GB, file="./WT/C_GB_WT_avg.rds")
saveRDS(thing$mat_FA, file="./WT/C_FA_WT_avg.rds")
saveRDS(thing$mat_FB, file="./WT/C_FB_WT_avg.rds")
saveRDS(thing$mat_D2AA, file="./WT/C_D2AA_WT_avg.rds")
saveRDS(thing$mat_D2AB, file="./WT/C_D2AB_WT_avg.rds")
saveRDS(thing$mat_D2BA, file="./WT/C_D2BA_WT_avg.rds")
saveRDS(thing$mat_D2BB, file="./WT/C_D2BB_WT_avg.rds")
rm(thing, C_WT_avg)

#(b) Pool
load("./WT/C_WT_pool.RData")
thing <- decompose(C_WT_pool)
saveRDS(thing$mat_D1, file="./WT/C_D1_WT_pool.rds")
saveRDS(thing$mat_GA, file="./WT/C_GA_WT_pool.rds")
saveRDS(thing$mat_GB, file="./WT/C_GB_WT_pool.rds")
saveRDS(thing$mat_FA, file="./WT/C_FA_WT_pool.rds")
saveRDS(thing$mat_FB, file="./WT/C_FB_WT_pool.rds")
saveRDS(thing$mat_D2AA, file="./WT/C_D2AA_WT_pool.rds")
saveRDS(thing$mat_D2AB, file="./WT/C_D2AB_WT_pool.rds")
saveRDS(thing$mat_D2BA, file="./WT/C_D2BA_WT_pool.rds")
saveRDS(thing$mat_D2BB, file="./WT/C_D2BB_WT_pool.rds")
rm(thing, C_WT_pool)

# Unpack decomposed matrix into pieces #####

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


#Big Cypress
C_D1_BC_avg  <-readRDS("./BC/C_D1_BC_avg.rds")
Kvals_contrib_D1_BC_avg <- toKvals_D1(C_D1_BC_avg)
save(Kvals_contrib_D1_BC_avg, file="./BC/Kvals_contrib_D1_BC_avg.RData")

C_D1_BC_pool  <-readRDS("./BC/C_D1_BC_pool.rds")
Kvals_contrib_D1_BC_pool <- toKvals_D1(C_D1_BC_pool)
save(Kvals_contrib_D1_BC_pool, file="./BC/Kvals_contrib_D1_BC_pool.RData")

#Cape Canaveral
C_D1_CC_avg  <-readRDS("./CC/C_D1_CC_avg.rds")
Kvals_contrib_D1_CC_avg <- toKvals_D1(C_D1_CC_avg)
save(Kvals_contrib_D1_CC_avg, file="./CC/Kvals_contrib_D1_CC_avg.RData")

C_D1_CC_pool  <-readRDS("./CC/C_D1_CC_pool.rds")
Kvals_contrib_D1_CC_pool <- toKvals_D1(C_D1_CC_pool)
save(Kvals_contrib_D1_CC_pool, file="./CC/Kvals_contrib_D1_CC_pool.RData")

#Chekika
C_D1_C_avg  <-readRDS("./C/C_D1_C_avg.rds")
Kvals_contrib_D1_C_avg <- toKvals_D1(C_D1_C_avg)
save(Kvals_contrib_D1_C_avg, file="./C/Kvals_contrib_D1_C_avg.RData")

C_D1_C_pool  <-readRDS("./C/C_D1_C_pool.rds")
Kvals_contrib_D1_C_pool <- toKvals_D1(C_D1_C_pool)
save(Kvals_contrib_D1_C_pool, file="./C/Kvals_contrib_D1_C_pool.RData")

#Fort Pierce
C_D1_FP_avg  <-readRDS("./FP/C_D1_FP_avg.rds")
Kvals_contrib_D1_FP_avg <- toKvals_D1(C_D1_FP_avg)
save(Kvals_contrib_D1_FP_avg, file="./FP/Kvals_contrib_D1_FP_avg.RData")

C_D1_FP_pool  <-readRDS("./FP/C_D1_FP_pool.rds")
Kvals_contrib_D1_FP_pool <- toKvals_D1(C_D1_FP_pool)
save(Kvals_contrib_D1_FP_pool, file="./FP/Kvals_contrib_D1_FP_pool.RData")

#Punta Gorda
C_D1_PG_avg  <-readRDS("./PG/C_D1_PG_avg.rds")
Kvals_contrib_D1_PG_avg <- toKvals_D1(C_D1_PG_avg)
save(Kvals_contrib_D1_PG_avg, file="./PG/Kvals_contrib_D1_PG_avg.RData")

C_D1_PG_pool  <-readRDS("./PG/C_D1_PG_pool.rds")
Kvals_contrib_D1_PG_pool <- toKvals_D1(C_D1_PG_pool)
save(Kvals_contrib_D1_PG_pool, file="./PG/Kvals_contrib_D1_PG_pool.RData")

#Wild Turkey
C_D1_WT_avg  <-readRDS("./WT/C_D1_WT_avg.rds")
Kvals_contrib_D1_WT_avg <- toKvals_D1(C_D1_WT_avg)
save(Kvals_contrib_D1_WT_avg, file="./WT/Kvals_contrib_D1_WT_avg.RData")

C_D1_WT_pool  <-readRDS("./WT/C_D1_WT_pool.rds")
Kvals_contrib_D1_WT_pool <- toKvals_D1(C_D1_WT_pool)
save(Kvals_contrib_D1_WT_pool, file="./WT/Kvals_contrib_D1_WT_pool.RData")

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


#Big Cypress
C_D2AA_BC_avg <- readRDS("./BC/C_D2AA_BC_avg.rds")
Kvals_contrib_D2AA_BC_avg <- toKvals_D2AA(C_D2AA_BC_avg)
save(Kvals_contrib_D2AA_BC_avg, file="./BC/Kvals_contrib_D2AA_BC_avg.RData")
rm(Kvals_contrib_D2AA_BC_avg, C_D2AA_BC_avg)

C_D2AB_BC_avg <- readRDS("./BC/C_D2AB_BC_avg.rds")
Kvals_contrib_D2AB_BC_avg <- toKvals_D2AB(C_D2AB_BC_avg)
save(Kvals_contrib_D2AB_BC_avg, file="./BC/Kvals_contrib_D2AB_BC_avg.RData")
rm(Kvals_contrib_D2AB_BC_avg, C_D2AB_BC_avg)

C_D2BA_BC_avg <- readRDS("./BC/C_D2BA_BC_avg.rds")
Kvals_contrib_D2BA_BC_avg <- toKvals_D2BA(C_D2BA_BC_avg)
save(Kvals_contrib_D2BA_BC_avg, file="./BC/Kvals_contrib_D2BA_BC_avg.RData")
rm(Kvals_contrib_D2BA_BC_avg, C_D2BA_BC_avg)

C_D2BB_BC_avg <- readRDS("./BC/C_D2BB_BC_avg.rds")
Kvals_contrib_D2BB_BC_avg <- toKvals_D2BB(C_D2BB_BC_avg)
save(Kvals_contrib_D2BB_BC_avg, file="./BC/Kvals_contrib_D2BB_BC_avg.RData")
rm(Kvals_contrib_D2BB_BC_avg, C_D2BB_BC_avg)

#(b) Pool
C_D2AA_BC_pool <- readRDS("./BC/C_D2AA_BC_pool.rds")
Kvals_contrib_D2AA_BC_pool <- toKvals_D2AA(C_D2AA_BC_pool)
save(Kvals_contrib_D2AA_BC_pool, file="./BC/Kvals_contrib_D2AA_BC_pool.RData")
rm(Kvals_contrib_D2AA_BC_pool, C_D2AA_BC_pool)

C_D2AB_BC_pool <- readRDS("./BC/C_D2AB_BC_pool.rds")
Kvals_contrib_D2AB_BC_pool <- toKvals_D2AB(C_D2AB_BC_pool)
save(Kvals_contrib_D2AB_BC_pool, file="./BC/Kvals_contrib_D2AB_BC_pool.RData")
rm(Kvals_contrib_D2AB_BC_pool, C_D2AB_BC_pool)

C_D2BA_BC_pool <- readRDS("./BC/C_D2BA_BC_pool.rds")
Kvals_contrib_D2BA_BC_pool <- toKvals_D2BA(C_D2BA_BC_pool)
save(Kvals_contrib_D2BA_BC_pool, file="./BC/Kvals_contrib_D2BA_BC_pool.RData")
rm(Kvals_contrib_D2BA_BC_pool, C_D2BA_BC_pool)

C_D2BB_BC_pool <- readRDS("./BC/C_D2BB_BC_pool.rds")
Kvals_contrib_D2BB_BC_pool <- toKvals_D2BB(C_D2BB_BC_pool)
save(Kvals_contrib_D2BB_BC_pool, file="./BC/Kvals_contrib_D2BB_BC_pool.RData")
rm(Kvals_contrib_D2BB_BC_pool, C_D2BB_BC_pool)


#Cape Canaveral
C_D2AA_CC_avg <- readRDS("./CC/C_D2AA_CC_avg.rds")
Kvals_contrib_D2AA_CC_avg <- toKvals_D2AA(C_D2AA_CC_avg)
save(Kvals_contrib_D2AA_CC_avg, file="./CC/Kvals_contrib_D2AA_CC_avg.RData")
rm(Kvals_contrib_D2AA_CC_avg, C_D2AA_CC_avg)

C_D2AB_CC_avg <- readRDS("./CC/C_D2AB_CC_avg.rds")
Kvals_contrib_D2AB_CC_avg <- toKvals_D2AB(C_D2AB_CC_avg)
save(Kvals_contrib_D2AB_CC_avg, file="./CC/Kvals_contrib_D2AB_CC_avg.RData")
rm(Kvals_contrib_D2AB_CC_avg, C_D2AB_CC_avg)

C_D2BA_CC_avg <- readRDS("./CC/C_D2BA_CC_avg.rds")
Kvals_contrib_D2BA_CC_avg <- toKvals_D2BA(C_D2BA_CC_avg)
save(Kvals_contrib_D2BA_CC_avg, file="./CC/Kvals_contrib_D2BA_CC_avg.RData")
rm(Kvals_contrib_D2BA_CC_avg, C_D2BA_CC_avg)

C_D2BB_CC_avg <- readRDS("./CC/C_D2BB_CC_avg.rds")
Kvals_contrib_D2BB_CC_avg <- toKvals_D2BB(C_D2BB_CC_avg)
save(Kvals_contrib_D2BB_CC_avg, file="./CC/Kvals_contrib_D2BB_CC_avg.RData")
rm(Kvals_contrib_D2BB_CC_avg, C_D2BB_CC_avg)

C_D2AA_BC_avg <- readRDS("./BC/C_D2AA_BC_avg.rds")
Kvals_contrib_D2AA_BC_avg <- toKvals_D2AA(C_D2AA_BC_avg)
save(Kvals_contrib_D2AA_BC_avg, file="./BC/Kvals_contrib_D2AA_BC_avg.RData")
rm(Kvals_contrib_D2AA_BC_avg, C_D2AA_BC_avg)

C_D2AB_BC_avg <- readRDS("./BC/C_D2AB_BC_avg.rds")
Kvals_contrib_D2AB_BC_avg <- toKvals_D2AB(C_D2AB_BC_avg)
save(Kvals_contrib_D2AB_BC_avg, file="./BC/Kvals_contrib_D2AB_BC_avg.RData")
rm(Kvals_contrib_D2AB_BC_avg, C_D2AB_BC_avg)

C_D2BA_BC_avg <- readRDS("./BC/C_D2BA_BC_avg.rds")
Kvals_contrib_D2BA_BC_avg <- toKvals_D2BA(C_D2BA_BC_avg)
save(Kvals_contrib_D2BA_BC_avg, file="./BC/Kvals_contrib_D2BA_BC_avg.RData")
rm(Kvals_contrib_D2BA_BC_avg, C_D2BA_BC_avg)

C_D2BB_BC_avg <- readRDS("./BC/C_D2BB_BC_avg.rds")
Kvals_contrib_D2BB_BC_avg <- toKvals_D2BB(C_D2BB_BC_avg)
save(Kvals_contrib_D2BB_BC_avg, file="./BC/Kvals_contrib_D2BB_BC_avg.RData")
rm(Kvals_contrib_D2BB_BC_avg, C_D2BB_BC_avg)

#(b) Pool
C_D2AA_CC_pool <- readRDS("./CC/C_D2AA_CC_pool.rds")
Kvals_contrib_D2AA_CC_pool <- toKvals_D2AA(C_D2AA_CC_pool)
save(Kvals_contrib_D2AA_CC_pool, file="./CC/Kvals_contrib_D2AA_CC_pool.RData")
rm(Kvals_contrib_D2AA_CC_pool, C_D2AA_CC_pool)

C_D2AB_CC_pool <- readRDS("./CC/C_D2AB_CC_pool.rds")
Kvals_contrib_D2AB_CC_pool <- toKvals_D2AB(C_D2AB_CC_pool)
save(Kvals_contrib_D2AB_CC_pool, file="./CC/Kvals_contrib_D2AB_CC_pool.RData")
rm(Kvals_contrib_D2AB_CC_pool, C_D2AB_CC_pool)

C_D2BA_CC_pool <- readRDS("./CC/C_D2BA_CC_pool.rds")
Kvals_contrib_D2BA_CC_pool <- toKvals_D2BA(C_D2BA_CC_pool)
save(Kvals_contrib_D2BA_CC_pool, file="./CC/Kvals_contrib_D2BA_CC_pool.RData")
rm(Kvals_contrib_D2BA_CC_pool, C_D2BA_CC_pool)

C_D2BB_CC_pool <- readRDS("./CC/C_D2BB_CC_pool.rds")
Kvals_contrib_D2BB_CC_pool <- toKvals_D2BB(C_D2BB_CC_pool)
save(Kvals_contrib_D2BB_CC_pool, file="./CC/Kvals_contrib_D2BB_CC_pool.RData")
rm(Kvals_contrib_D2BB_CC_pool, C_D2BB_CC_pool)

#Chekika
C_D2AA_C_avg <- readRDS("./C/C_D2AA_C_avg.rds")
Kvals_contrib_D2AA_C_avg <- toKvals_D2AA(C_D2AA_C_avg)
save(Kvals_contrib_D2AA_C_avg, file="./C/Kvals_contrib_D2AA_C_avg.RData")
rm(Kvals_contrib_D2AA_C_avg, C_D2AA_C_avg)

C_D2AB_C_avg <- readRDS("./C/C_D2AB_C_avg.rds")
Kvals_contrib_D2AB_C_avg <- toKvals_D2AB(C_D2AB_C_avg)
save(Kvals_contrib_D2AB_C_avg, file="./C/Kvals_contrib_D2AB_C_avg.RData")
rm(Kvals_contrib_D2AB_C_avg, C_D2AB_C_avg)

C_D2BA_C_avg <- readRDS("./C/C_D2BA_C_avg.rds")
Kvals_contrib_D2BA_C_avg <- toKvals_D2BA(C_D2BA_C_avg)
save(Kvals_contrib_D2BA_C_avg, file="./C/Kvals_contrib_D2BA_C_avg.RData")
rm(Kvals_contrib_D2BA_C_avg, C_D2BA_C_avg)

C_D2BB_C_avg <- readRDS("./C/C_D2BB_C_avg.rds")
Kvals_contrib_D2BB_C_avg <- toKvals_D2BB(C_D2BB_C_avg)
save(Kvals_contrib_D2BB_C_avg, file="./C/Kvals_contrib_D2BB_C_avg.RData")
rm(Kvals_contrib_D2BB_C_avg, C_D2BB_C_avg)

#(b) Pool
C_D2AA_C_pool <- readRDS("./C/C_D2AA_C_pool.rds")
Kvals_contrib_D2AA_C_pool <- toKvals_D2AA(C_D2AA_C_pool)
save(Kvals_contrib_D2AA_C_pool, file="./C/Kvals_contrib_D2AA_C_pool.RData")
rm(Kvals_contrib_D2AA_C_pool, C_D2AA_C_pool)

C_D2AB_C_pool <- readRDS("./C/C_D2AB_C_pool.rds")
Kvals_contrib_D2AB_C_pool <- toKvals_D2AB(C_D2AB_C_pool)
save(Kvals_contrib_D2AB_C_pool, file="./C/Kvals_contrib_D2AB_C_pool.RData")
rm(Kvals_contrib_D2AB_C_pool, C_D2AB_C_pool)

C_D2BA_C_pool <- readRDS("./C/C_D2BA_C_pool.rds")
Kvals_contrib_D2BA_C_pool <- toKvals_D2BA(C_D2BA_C_pool)
save(Kvals_contrib_D2BA_C_pool, file="./C/Kvals_contrib_D2BA_C_pool.RData")
rm(Kvals_contrib_D2BA_C_pool, C_D2BA_C_pool)

C_D2BB_C_pool <- readRDS("./C/C_D2BB_C_pool.rds")
Kvals_contrib_D2BB_C_pool <- toKvals_D2BB(C_D2BB_C_pool)
save(Kvals_contrib_D2BB_C_pool, file="./C/Kvals_contrib_D2BB_C_pool.RData")
rm(Kvals_contrib_D2BB_C_pool, C_D2BB_C_pool)

#Fort Pierce
C_D2AA_FP_avg <- readRDS("./FP/C_D2AA_FP_avg.rds")
Kvals_contrib_D2AA_FP_avg <- toKvals_D2AA(C_D2AA_FP_avg)
save(Kvals_contrib_D2AA_FP_avg, file="./FP/Kvals_contrib_D2AA_FP_avg.RData")
rm(Kvals_contrib_D2AA_FP_avg, C_D2AA_FP_avg)

C_D2AB_FP_avg <- readRDS("./FP/C_D2AB_FP_avg.rds")
Kvals_contrib_D2AB_FP_avg <- toKvals_D2AB(C_D2AB_FP_avg)
save(Kvals_contrib_D2AB_FP_avg, file="./FP/Kvals_contrib_D2AB_FP_avg.RData")
rm(Kvals_contrib_D2AB_FP_avg, C_D2AB_FP_avg)

C_D2BA_FP_avg <- readRDS("./FP/C_D2BA_FP_avg.rds")
Kvals_contrib_D2BA_FP_avg <- toKvals_D2BA(C_D2BA_FP_avg)
save(Kvals_contrib_D2BA_FP_avg, file="./FP/Kvals_contrib_D2BA_FP_avg.RData")
rm(Kvals_contrib_D2BA_FP_avg, C_D2BA_FP_avg)

C_D2BB_FP_avg <- readRDS("./FP/C_D2BB_FP_avg.rds")
Kvals_contrib_D2BB_FP_avg <- toKvals_D2BB(C_D2BB_FP_avg)
save(Kvals_contrib_D2BB_FP_avg, file="./FP/Kvals_contrib_D2BB_FP_avg.RData")
rm(Kvals_contrib_D2BB_FP_avg, C_D2BB_FP_avg)

#(b) Pool
C_D2AA_FP_pool <- readRDS("./FP/C_D2AA_FP_pool.rds")
Kvals_contrib_D2AA_FP_pool <- toKvals_D2AA(C_D2AA_FP_pool)
save(Kvals_contrib_D2AA_FP_pool, file="./FP/Kvals_contrib_D2AA_FP_pool.RData")
rm(Kvals_contrib_D2AA_FP_pool, C_D2AA_FP_pool)

C_D2AB_FP_pool <- readRDS("./FP/C_D2AB_FP_pool.rds")
Kvals_contrib_D2AB_FP_pool <- toKvals_D2AB(C_D2AB_FP_pool)
save(Kvals_contrib_D2AB_FP_pool, file="./FP/Kvals_contrib_D2AB_FP_pool.RData")
rm(Kvals_contrib_D2AB_FP_pool, C_D2AB_FP_pool)

C_D2BA_FP_pool <- readRDS("./FP/C_D2BA_FP_pool.rds")
Kvals_contrib_D2BA_FP_pool <- toKvals_D2BA(C_D2BA_FP_pool)
save(Kvals_contrib_D2BA_FP_pool, file="./FP/Kvals_contrib_D2BA_FP_pool.RData")
rm(Kvals_contrib_D2BA_FP_pool, C_D2BA_FP_pool)

C_D2BB_FP_pool <- readRDS("./FP/C_D2BB_FP_pool.rds")
Kvals_contrib_D2BB_FP_pool <- toKvals_D2BB(C_D2BB_FP_pool)
save(Kvals_contrib_D2BB_FP_pool, file="./FP/Kvals_contrib_D2BB_FP_pool.RData")
rm(Kvals_contrib_D2BB_FP_pool, C_D2BB_FP_pool)

#Punta Gorda
C_D2AA_PG_avg <- readRDS("./PG/C_D2AA_PG_avg.rds")
Kvals_contrib_D2AA_PG_avg <- toKvals_D2AA(C_D2AA_PG_avg)
save(Kvals_contrib_D2AA_PG_avg, file="./PG/Kvals_contrib_D2AA_PG_avg.RData")
rm(Kvals_contrib_D2AA_PG_avg, C_D2AA_PG_avg)

C_D2AB_PG_avg <- readRDS("./PG/C_D2AB_PG_avg.rds")
Kvals_contrib_D2AB_PG_avg <- toKvals_D2AB(C_D2AB_PG_avg)
save(Kvals_contrib_D2AB_PG_avg, file="./PG/Kvals_contrib_D2AB_PG_avg.RData")
rm(Kvals_contrib_D2AB_PG_avg, C_D2AB_PG_avg)

C_D2BA_PG_avg <- readRDS("./PG/C_D2BA_PG_avg.rds")
Kvals_contrib_D2BA_PG_avg <- toKvals_D2BA(C_D2BA_PG_avg)
save(Kvals_contrib_D2BA_PG_avg, file="./PG/Kvals_contrib_D2BA_PG_avg.RData")
rm(Kvals_contrib_D2BA_PG_avg, C_D2BA_PG_avg)

C_D2BB_PG_avg <- readRDS("./PG/C_D2BB_PG_avg.rds")
Kvals_contrib_D2BB_PG_avg <- toKvals_D2BB(C_D2BB_PG_avg)
save(Kvals_contrib_D2BB_PG_avg, file="./PG/Kvals_contrib_D2BB_PG_avg.RData")
rm(Kvals_contrib_D2BB_PG_avg, C_D2BB_PG_avg)

#(b) Pool
C_D2AA_PG_pool <- readRDS("./PG/C_D2AA_PG_pool.rds")
Kvals_contrib_D2AA_PG_pool <- toKvals_D2AA(C_D2AA_PG_pool)
save(Kvals_contrib_D2AA_PG_pool, file="./PG/Kvals_contrib_D2AA_PG_pool.RData")
rm(Kvals_contrib_D2AA_PG_pool, C_D2AA_PG_pool)

C_D2AB_PG_pool <- readRDS("./PG/C_D2AB_PG_pool.rds")
Kvals_contrib_D2AB_PG_pool <- toKvals_D2AB(C_D2AB_PG_pool)
save(Kvals_contrib_D2AB_PG_pool, file="./PG/Kvals_contrib_D2AB_PG_pool.RData")
rm(Kvals_contrib_D2AB_PG_pool, C_D2AB_PG_pool)

C_D2BA_PG_pool <- readRDS("./PG/C_D2BA_PG_pool.rds")
Kvals_contrib_D2BA_PG_pool <- toKvals_D2BA(C_D2BA_PG_pool)
save(Kvals_contrib_D2BA_PG_pool, file="./PG/Kvals_contrib_D2BA_PG_pool.RData")
rm(Kvals_contrib_D2BA_PG_pool, C_D2BA_PG_pool)

C_D2BB_PG_pool <- readRDS("./PG/C_D2BB_PG_pool.rds")
Kvals_contrib_D2BB_PG_pool <- toKvals_D2BB(C_D2BB_PG_pool)
save(Kvals_contrib_D2BB_PG_pool, file="./PG/Kvals_contrib_D2BB_PG_pool.RData")
rm(Kvals_contrib_D2BB_PG_pool, C_D2BB_PG_pool)

#Wild Turkey
C_D2AA_WT_avg <- readRDS("./WT/C_D2AA_WT_avg.rds")
Kvals_contrib_D2AA_WT_avg <- toKvals_D2AA(C_D2AA_WT_avg)
save(Kvals_contrib_D2AA_WT_avg, file="./WT/Kvals_contrib_D2AA_WT_avg.RData")
rm(Kvals_contrib_D2AA_WT_avg, C_D2AA_WT_avg)

C_D2AB_WT_avg <- readRDS("./WT/C_D2AB_WT_avg.rds")
Kvals_contrib_D2AB_WT_avg <- toKvals_D2AB(C_D2AB_WT_avg)
save(Kvals_contrib_D2AB_WT_avg, file="./WT/Kvals_contrib_D2AB_WT_avg.RData")
rm(Kvals_contrib_D2AB_WT_avg, C_D2AB_WT_avg)

C_D2BA_WT_avg <- readRDS("./WT/C_D2BA_WT_avg.rds")
Kvals_contrib_D2BA_WT_avg <- toKvals_D2BA(C_D2BA_WT_avg)
save(Kvals_contrib_D2BA_WT_avg, file="./WT/Kvals_contrib_D2BA_WT_avg.RData")
rm(Kvals_contrib_D2BA_WT_avg, C_D2BA_WT_avg)

C_D2BB_WT_avg <- readRDS("./WT/C_D2BB_WT_avg.rds")
Kvals_contrib_D2BB_WT_avg <- toKvals_D2BB(C_D2BB_WT_avg)
save(Kvals_contrib_D2BB_WT_avg, file="./WT/Kvals_contrib_D2BB_WT_avg.RData")
rm(Kvals_contrib_D2BB_WT_avg, C_D2BB_WT_avg)

#(b) Pool
C_D2AA_WT_pool <- readRDS("./WT/C_D2AA_WT_pool.rds")
Kvals_contrib_D2AA_WT_pool <- toKvals_D2AA(C_D2AA_WT_pool)
save(Kvals_contrib_D2AA_WT_pool, file="./WT/Kvals_contrib_D2AA_WT_pool.RData")
rm(Kvals_contrib_D2AA_WT_pool, C_D2AA_WT_pool)

C_D2AB_WT_pool <- readRDS("./WT/C_D2AB_WT_pool.rds")
Kvals_contrib_D2AB_WT_pool <- toKvals_D2AB(C_D2AB_WT_pool)
save(Kvals_contrib_D2AB_WT_pool, file="./WT/Kvals_contrib_D2AB_WT_pool.RData")
rm(Kvals_contrib_D2AB_WT_pool, C_D2AB_WT_pool)

C_D2BA_WT_pool <- readRDS("./WT/C_D2BA_WT_pool.rds")
Kvals_contrib_D2BA_WT_pool <- toKvals_D2BA(C_D2BA_WT_pool)
save(Kvals_contrib_D2BA_WT_pool, file="./WT/Kvals_contrib_D2BA_WT_pool.RData")
rm(Kvals_contrib_D2BA_WT_pool, C_D2BA_WT_pool)

C_D2BB_WT_pool <- readRDS("./WT/C_D2BB_WT_pool.rds")
Kvals_contrib_D2BB_WT_pool <- toKvals_D2BB(C_D2BB_WT_pool)
save(Kvals_contrib_D2BB_WT_pool, file="./WT/Kvals_contrib_D2BB_WT_pool.RData")
rm(Kvals_contrib_D2BB_WT_pool, C_D2BB_WT_pool)

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
#Big Cypress

C_FA_BC_avg  <-readRDS("./BC/C_FA_BC_avg.rds")
Kvals_contrib_FA_BC_avg  <- toKvals_FA(C_FA_BC_avg)
save(Kvals_contrib_FA_BC_avg, file="./BC/Kvals_contrib_FA_BC_avg.RData")
rm(C_FA_BC_avg, Kvals_contrib_FA_BC_avg)

C_FB_BC_avg  <-readRDS("./BC/C_FB_BC_avg.rds")
Kvals_contrib_FB_BC_avg  <- toKvals_FB(C_FB_BC_avg)
save(Kvals_contrib_FB_BC_avg, file="./BC/Kvals_contrib_FB_BC_avg.RData")
rm(C_FB_BC_avg, Kvals_contrib_FB_BC_avg)

#(b) Pool
C_FA_BC_pool  <-readRDS("./BC/C_FA_BC_pool.rds")
Kvals_contrib_FA_BC_pool  <- toKvals_FA(C_FA_BC_pool)
save(Kvals_contrib_FA_BC_pool, file="./BC/Kvals_contrib_FA_BC_pool.RData")
rm(C_FA_BC_pool, Kvals_contrib_FA_BC_pool)

C_FB_BC_pool  <-readRDS("./BC/C_FB_BC_pool.rds")
Kvals_contrib_FB_BC_pool  <- toKvals_FB(C_FB_BC_pool)
save(Kvals_contrib_FB_BC_pool, file="./BC/Kvals_contrib_FB_BC_pool.RData")
rm(C_FB_BC_pool, Kvals_contrib_FB_BC_pool)

#Cape Canaveral

C_FA_CC_avg  <-readRDS("./CC/C_FA_CC_avg.rds")
Kvals_contrib_FA_CC_avg  <- toKvals_FA(C_FA_CC_avg)
save(Kvals_contrib_FA_CC_avg, file="./CC/Kvals_contrib_FA_CC_avg.RData")
rm(C_FA_CC_avg, Kvals_contrib_FA_CC_avg)

C_FB_CC_avg  <-readRDS("./CC/C_FB_CC_avg.rds")
Kvals_contrib_FB_CC_avg  <- toKvals_FB(C_FB_CC_avg)
save(Kvals_contrib_FB_CC_avg, file="./CC/Kvals_contrib_FB_CC_avg.RData")
rm(C_FB_CC_avg, Kvals_contrib_FB_CC_avg)

#(b) Pool
C_FA_CC_pool  <-readRDS("./CC/C_FA_CC_pool.rds")
Kvals_contrib_FA_CC_pool  <- toKvals_FA(C_FA_CC_pool)
save(Kvals_contrib_FA_CC_pool, file="./CC/Kvals_contrib_FA_CC_pool.RData")
rm(C_FA_CC_pool, Kvals_contrib_FA_CC_pool)

C_FB_CC_pool  <-readRDS("./CC/C_FB_CC_pool.rds")
Kvals_contrib_FB_CC_pool  <- toKvals_FB(C_FB_CC_pool)
save(Kvals_contrib_FB_CC_pool, file="./CC/Kvals_contrib_FB_CC_pool.RData")
rm(C_FB_CC_pool, Kvals_contrib_FB_CC_pool)

#Chekika

C_FA_C_avg  <-readRDS("./C/C_FA_C_avg.rds")
Kvals_contrib_FA_C_avg  <- toKvals_FA(C_FA_C_avg)
save(Kvals_contrib_FA_C_avg, file="./C/Kvals_contrib_FA_C_avg.RData")
rm(C_FA_C_avg, Kvals_contrib_FA_C_avg)

C_FB_C_avg  <-readRDS("./C/C_FB_C_avg.rds")
Kvals_contrib_FB_C_avg  <- toKvals_FB(C_FB_C_avg)
save(Kvals_contrib_FB_C_avg, file="./C/Kvals_contrib_FB_C_avg.RData")
rm(C_FB_C_avg, Kvals_contrib_FB_C_avg)

#(b) Pool
C_FA_C_pool  <-readRDS("./C/C_FA_C_pool.rds")
Kvals_contrib_FA_C_pool  <- toKvals_FA(C_FA_C_pool)
save(Kvals_contrib_FA_C_pool, file="./C/Kvals_contrib_FA_C_pool.RData")
rm(C_FA_C_pool, Kvals_contrib_FA_C_pool)

C_FB_C_pool  <-readRDS("./C/C_FB_C_pool.rds")
Kvals_contrib_FB_C_pool  <- toKvals_FB(C_FB_C_pool)
save(Kvals_contrib_FB_C_pool, file="./C/Kvals_contrib_FB_C_pool.RData")
rm(C_FB_C_pool, Kvals_contrib_FB_C_pool)

#Fort Pierce

C_FA_FP_avg  <-readRDS("./FP/C_FA_FP_avg.rds")
Kvals_contrib_FA_FP_avg  <- toKvals_FA(C_FA_FP_avg)
save(Kvals_contrib_FA_FP_avg, file="./FP/Kvals_contrib_FA_FP_avg.RData")
rm(C_FA_FP_avg, Kvals_contrib_FA_FP_avg)

C_FB_FP_avg  <-readRDS("./FP/C_FB_FP_avg.rds")
Kvals_contrib_FB_FP_avg  <- toKvals_FB(C_FB_FP_avg)
save(Kvals_contrib_FB_FP_avg, file="./FP/Kvals_contrib_FB_FP_avg.RData")
rm(C_FB_FP_avg, Kvals_contrib_FB_FP_avg)

#(b) Pool
C_FA_FP_pool  <-readRDS("./FP/C_FA_FP_pool.rds")
Kvals_contrib_FA_FP_pool  <- toKvals_FA(C_FA_FP_pool)
save(Kvals_contrib_FA_FP_pool, file="./FP/Kvals_contrib_FA_FP_pool.RData")
rm(C_FA_FP_pool, Kvals_contrib_FA_FP_pool)

C_FB_FP_pool  <-readRDS("./FP/C_FB_FP_pool.rds")
Kvals_contrib_FB_FP_pool  <- toKvals_FB(C_FB_FP_pool)
save(Kvals_contrib_FB_FP_pool, file="./FP/Kvals_contrib_FB_FP_pool.RData")
rm(C_FB_FP_pool, Kvals_contrib_FB_FP_pool)

#Punta Gorda

C_FA_PG_avg  <-readRDS("./PG/C_FA_PG_avg.rds")
Kvals_contrib_FA_PG_avg  <- toKvals_FA(C_FA_PG_avg)
save(Kvals_contrib_FA_PG_avg, file="./PG/Kvals_contrib_FA_PG_avg.RData")
rm(C_FA_PG_avg, Kvals_contrib_FA_PG_avg)

C_FB_PG_avg  <-readRDS("./PG/C_FB_PG_avg.rds")
Kvals_contrib_FB_PG_avg  <- toKvals_FB(C_FB_PG_avg)
save(Kvals_contrib_FB_PG_avg, file="./PG/Kvals_contrib_FB_PG_avg.RData")
rm(C_FB_PG_avg, Kvals_contrib_FB_PG_avg)

#(b) Pool
C_FA_PG_pool  <-readRDS("./PG/C_FA_PG_pool.rds")
Kvals_contrib_FA_PG_pool  <- toKvals_FA(C_FA_PG_pool)
save(Kvals_contrib_FA_PG_pool, file="./PG/Kvals_contrib_FA_PG_pool.RData")
rm(C_FA_PG_pool, Kvals_contrib_FA_PG_pool)

C_FB_PG_pool  <-readRDS("./PG/C_FB_PG_pool.rds")
Kvals_contrib_FB_PG_pool  <- toKvals_FB(C_FB_PG_pool)
save(Kvals_contrib_FB_PG_pool, file="./PG/Kvals_contrib_FB_PG_pool.RData")
rm(C_FB_PG_pool, Kvals_contrib_FB_PG_pool)

#Wild Turkey

C_FA_WT_avg  <-readRDS("./WT/C_FA_WT_avg.rds")
Kvals_contrib_FA_WT_avg  <- toKvals_FA(C_FA_WT_avg)
save(Kvals_contrib_FA_WT_avg, file="./WT/Kvals_contrib_FA_WT_avg.RData")
rm(C_FA_WT_avg, Kvals_contrib_FA_WT_avg)

C_FB_WT_avg  <-readRDS("./WT/C_FB_WT_avg.rds")
Kvals_contrib_FB_WT_avg  <- toKvals_FB(C_FB_WT_avg)
save(Kvals_contrib_FB_WT_avg, file="./WT/Kvals_contrib_FB_WT_avg.RData")
rm(C_FB_WT_avg, Kvals_contrib_FB_WT_avg)

#(b) Pool
C_FA_WT_pool  <-readRDS("./WT/C_FA_WT_pool.rds")
Kvals_contrib_FA_WT_pool  <- toKvals_FA(C_FA_WT_pool)
save(Kvals_contrib_FA_WT_pool, file="./WT/Kvals_contrib_FA_WT_pool.RData")
rm(C_FA_WT_pool, Kvals_contrib_FA_WT_pool)

C_FB_WT_pool  <-readRDS("./WT/C_FB_WT_pool.rds")
Kvals_contrib_FB_WT_pool  <- toKvals_FB(C_FB_WT_pool)
save(Kvals_contrib_FB_WT_pool, file="./WT/Kvals_contrib_FB_WT_pool.RData")
rm(C_FB_WT_pool, Kvals_contrib_FB_WT_pool)






toKvals_GA = function(elas_GA) {
  plop1=function(i, j) {(j-1)*m3a + i}
  plop2=function(i, j) {(j-1)*m1 + i}
  Plop1=outer(1:m3a,1:m4a,plop1); 
  Plop2=outer(1:m1, 1:m2, plop2);
  
  Kvals_elas_GA=array(0, c(m3a, m4a, m1, m2))
  
  for(i in 1:m1) {
    for (j in 1:m2) {
      for (k in 1:m3a) {
        kvals=elas_GA[Plop1[k, 1:m4a], Plop2[i,j]]
        Kvals_elas_GA[k, 1:m4a, i, j]=kvals
      }}
    cat(i, "\n");
  }
  return(Kvals_elas_GA)
}

toKvals_GB = function(elas_GB) {
  plop1=function(i, j) {(j-1)*m3b + i}
  plop2=function(i, j) {(j-1)*m1 + i}
  Plop1=outer(1:m3b,1:m4b,plop1); 
  Plop2=outer(1:m1, 1:m2, plop2);
  
  Kvals_elas_GB=array(0, c(m3b, m4b, m1, m2))
  
  for(i in 1:m1) {
    for (j in 1:m2) {
      for (k in 1:m3b) {
        kvals=elas_GB[Plop1[k, 1:m4b], Plop2[i,j]]
        Kvals_elas_GB[k, 1:m4b, i, j]=kvals
      }}
    cat(i, "\n");
  }
  return(Kvals_elas_GB)
}

#Big Cypress
C_GA_BC_avg  <-readRDS("./BC/C_GA_BC_avg.rds")
Kvals_contrib_GA_BC_avg  <- toKvals_GA(C_GA_BC_avg)
save(Kvals_contrib_GA_BC_avg, file="./BC/Kvals_contrib_GA_BC_avg.RData")
rm(C_GA_BC_avg, Kvals_contrib_GA_BC_avg)

C_GB_BC_avg <- readRDS("./BC/C_GB_BC_avg.rds")
Kvals_contrib_GB_BC_avg <- toKvals_GB(C_GB_BC_avg)
save(Kvals_contrib_GB_BC_avg, file="./BC/Kvals_contrib_GB_BC_avg.RData")
rm(C_GB_BC_avg, Kvals_contrib_GB_BC_avg)

#(b) Pool
C_GA_BC_pool  <-readRDS("./BC/C_GA_BC_pool.rds")
Kvals_contrib_GA_BC_pool  <- toKvals_GA(C_GA_BC_pool)
save(Kvals_contrib_GA_BC_pool, file="./BC/Kvals_contrib_GA_BC_pool.RData")
rm(C_GA_BC_pool, Kvals_contrib_GA_BC_pool)

C_GB_BC_pool <- readRDS("./BC/C_GB_BC_pool.rds")
Kvals_contrib_GB_BC_pool <- toKvals_GB(C_GB_BC_pool)
save(Kvals_contrib_GB_BC_pool, file="./BC/Kvals_contrib_GB_BC_pool.RData")
rm(C_GB_BC_pool, Kvals_contrib_GB_BC_pool)


#Cape Canaveral
C_GA_CC_avg  <-readRDS("./CC/C_GA_CC_avg.rds")
Kvals_contrib_GA_CC_avg  <- toKvals_GA(C_GA_CC_avg)
save(Kvals_contrib_GA_CC_avg, file="./CC/Kvals_contrib_GA_CC_avg.RData")
rm(C_GA_CC_avg, Kvals_contrib_GA_CC_avg)

C_GB_CC_avg <- readRDS("./CC/C_GB_CC_avg.rds")
Kvals_contrib_GB_CC_avg <- toKvals_GB(C_GB_CC_avg)
save(Kvals_contrib_GB_CC_avg, file="./CC/Kvals_contrib_GB_CC_avg.RData")
rm(C_GB_CC_avg, Kvals_contrib_GB_CC_avg)

#(b) Pool
C_GA_CC_pool  <-readRDS("./CC/C_GA_CC_pool.rds")
Kvals_contrib_GA_CC_pool  <- toKvals_GA(C_GA_CC_pool)
save(Kvals_contrib_GA_CC_pool, file="./CC/Kvals_contrib_GA_CC_pool.RData")
rm(C_GA_CC_pool, Kvals_contrib_GA_CC_pool)

C_GB_CC_pool <- readRDS("./CC/C_GB_CC_pool.rds")
Kvals_contrib_GB_CC_pool <- toKvals_GB(C_GB_CC_pool)
save(Kvals_contrib_GB_CC_pool, file="./CC/Kvals_contrib_GB_CC_pool.RData")
rm(C_GB_CC_pool, Kvals_contrib_GB_CC_pool)

#Chekika
C_GA_C_avg  <-readRDS("./C/C_GA_C_avg.rds")
Kvals_contrib_GA_C_avg  <- toKvals_GA(C_GA_C_avg)
save(Kvals_contrib_GA_C_avg, file="./C/Kvals_contrib_GA_C_avg.RData")
rm(C_GA_C_avg, Kvals_contrib_GA_C_avg)

C_GB_C_avg <- readRDS("./C/C_GB_C_avg.rds")
Kvals_contrib_GB_C_avg <- toKvals_GB(C_GB_C_avg)
save(Kvals_contrib_GB_C_avg, file="./C/Kvals_contrib_GB_C_avg.RData")
rm(C_GB_C_avg, Kvals_contrib_GB_C_avg)

#(b) Pool
C_GA_C_pool  <-readRDS("./C/C_GA_C_pool.rds")
Kvals_contrib_GA_C_pool  <- toKvals_GA(C_GA_C_pool)
save(Kvals_contrib_GA_C_pool, file="./C/Kvals_contrib_GA_C_pool.RData")
rm(C_GA_C_pool, Kvals_contrib_GA_C_pool)

C_GB_C_pool <- readRDS("./C/C_GB_C_pool.rds")
Kvals_contrib_GB_C_pool <- toKvals_GB(C_GB_C_pool)
save(Kvals_contrib_GB_C_pool, file="./C/Kvals_contrib_GB_C_pool.RData")
rm(C_GB_C_pool, Kvals_contrib_GB_C_pool)

#Fort Pierce
C_GA_FP_avg  <-readRDS("./FP/C_GA_FP_avg.rds")
Kvals_contrib_GA_FP_avg  <- toKvals_GA(C_GA_FP_avg)
save(Kvals_contrib_GA_FP_avg, file="./FP/Kvals_contrib_GA_FP_avg.RData")
rm(C_GA_FP_avg, Kvals_contrib_GA_FP_avg)

C_GB_FP_avg <- readRDS("./FP/C_GB_FP_avg.rds")
Kvals_contrib_GB_FP_avg <- toKvals_GB(C_GB_FP_avg)
save(Kvals_contrib_GB_FP_avg, file="./FP/Kvals_contrib_GB_FP_avg.RData")
rm(C_GB_FP_avg, Kvals_contrib_GB_FP_avg)

#(b) Pool
C_GA_FP_pool  <-readRDS("./FP/C_GA_FP_pool.rds")
Kvals_contrib_GA_FP_pool  <- toKvals_GA(C_GA_FP_pool)
save(Kvals_contrib_GA_FP_pool, file="./FP/Kvals_contrib_GA_FP_pool.RData")
rm(C_GA_FP_pool, Kvals_contrib_GA_FP_pool)

C_GB_FP_pool <- readRDS("./FP/C_GB_FP_pool.rds")
Kvals_contrib_GB_FP_pool <- toKvals_GB(C_GB_FP_pool)
save(Kvals_contrib_GB_FP_pool, file="./FP/Kvals_contrib_GB_FP_pool.RData")
rm(C_GB_FP_pool, Kvals_contrib_GB_FP_pool)

#Punta Gorda
C_GA_PG_avg  <-readRDS("./PG/C_GA_PG_avg.rds")
Kvals_contrib_GA_PG_avg  <- toKvals_GA(C_GA_PG_avg)
save(Kvals_contrib_GA_PG_avg, file="./PG/Kvals_contrib_GA_PG_avg.RData")
rm(C_GA_PG_avg, Kvals_contrib_GA_PG_avg)

C_GB_PG_avg <- readRDS("./PG/C_GB_PG_avg.rds")
Kvals_contrib_GB_PG_avg <- toKvals_GB(C_GB_PG_avg)
save(Kvals_contrib_GB_PG_avg, file="./PG/Kvals_contrib_GB_PG_avg.RData")
rm(C_GB_PG_avg, Kvals_contrib_GB_PG_avg)
#(b) Pool
C_GA_PG_pool  <-readRDS("./PG/C_GA_PG_pool.rds")
Kvals_contrib_GA_PG_pool  <- toKvals_GA(C_GA_PG_pool)
save(Kvals_contrib_GA_PG_pool, file="./PG/Kvals_contrib_GA_PG_pool.RData")
rm(C_GA_PG_pool, Kvals_contrib_GA_PG_pool)

C_GB_PG_pool <- readRDS("./PG/C_GB_PG_pool.rds")
Kvals_contrib_GB_PG_pool <- toKvals_GB(C_GB_PG_pool)
save(Kvals_contrib_GB_PG_pool, file="./PG/Kvals_contrib_GB_PG_pool.RData")
rm(C_GB_PG_pool, Kvals_contrib_GB_PG_pool)

#Wild Turkey
C_GA_WT_avg  <-readRDS("./WT/C_GA_WT_avg.rds")
Kvals_contrib_GA_WT_avg  <- toKvals_GA(C_GA_WT_avg)
save(Kvals_contrib_GA_WT_avg, file="./WT/Kvals_contrib_GA_WT_avg.RData")
rm(C_GA_WT_avg, Kvals_contrib_GA_WT_avg)

C_GB_WT_avg <- readRDS("./WT/C_GB_WT_avg.rds")
Kvals_contrib_GB_WT_avg <- toKvals_GB(C_GB_WT_avg)
save(Kvals_contrib_GB_WT_avg, file="./WT/Kvals_contrib_GB_WT_avg.RData")
rm(C_GB_WT_avg, Kvals_contrib_GB_WT_avg)

#(b) Pool
C_GA_WT_pool  <-readRDS("./WT/C_GA_WT_pool.rds")
Kvals_contrib_GA_WT_pool  <- toKvals_GA(C_GA_WT_pool)
save(Kvals_contrib_GA_WT_pool, file="./WT/Kvals_contrib_GA_WT_pool.RData")
rm(C_GA_WT_pool, Kvals_contrib_GA_WT_pool)

C_GB_WT_pool <- readRDS("./WT/C_GB_WT_pool.rds")
Kvals_contrib_GB_WT_pool <- toKvals_GB(C_GB_WT_pool)
save(Kvals_contrib_GB_WT_pool, file="./WT/Kvals_contrib_GB_WT_pool.RData")
rm(C_GB_WT_pool, Kvals_contrib_GB_WT_pool)

# Convert to totals: #####
#Big Cypress

load("./BC/Kvals_contrib_D1_BC_avg.RData")
load("./BC/Kvals_contrib_D2AA_BC_avg.RData")
load("./BC/Kvals_contrib_D2AB_BC_avg.RData")
load("./BC/Kvals_contrib_D2BA_BC_avg.RData")
load("./BC/Kvals_contrib_D2BB_BC_avg.RData")
load("./BC/Kvals_contrib_FA_BC_avg.RData")
load("./BC/Kvals_contrib_FB_BC_avg.RData")
load("./BC/Kvals_contrib_GA_BC_avg.RData")
load("./BC/Kvals_contrib_GB_BC_avg.RData")

total.contrib_D1_BC_avg<-apply(Kvals_contrib_D1_BC_avg, c(3,4), sum);

total.contrib_D2AA_BC_avg<-apply(Kvals_contrib_D2AA_BC_avg, c(3,4), sum);
total.contrib_D2AB_BC_avg<-apply(Kvals_contrib_D2AB_BC_avg, c(3,4),sum)
total.contrib_D2BA_BC_avg<-apply(Kvals_contrib_D2BA_BC_avg, c(3,4),sum)
total.contrib_D2BB_BC_avg<-apply(Kvals_contrib_D2BB_BC_avg, c(3,4), sum)
total.contrib_FA_BC_avg<-apply(Kvals_contrib_FA_BC_avg, c(3,4), sum);
total.contrib_FB_BC_avg<-apply(Kvals_contrib_FB_BC_avg, c(3,4), sum)
total.contrib_GA_BC_avg<-apply(Kvals_contrib_GA_BC_avg, c(3,4), sum);
total.contrib_GB_BC_avg<-apply(Kvals_contrib_GB_BC_avg, c(3,4), sum)

save(total.contrib_D1_BC_avg, total.contrib_D2AA_BC_avg,
     total.contrib_D2AB_BC_avg, total.contrib_D2BA_BC_avg,
     total.contrib_D2BB_BC_avg, total.contrib_FA_BC_avg, 
     total.contrib_FB_BC_avg,total.contrib_GA_BC_avg, 
     total.contrib_GB_BC_avg,
     file="./BC/total.contrib_BC_avg.RData")

BC_avg<-c( 
  sum(rowSums(total.contrib_D1_BC_avg)), 
  sum(rowSums(total.contrib_GA_BC_avg)),
  sum(rowSums(total.contrib_GB_BC_avg)),
  sum(rowSums(total.contrib_FA_BC_avg)),
  sum(rowSums(total.contrib_FB_BC_avg)),
  sum(rowSums(total.contrib_D2AA_BC_avg)),
  sum(rowSums(total.contrib_D2AB_BC_avg)),
  sum(rowSums(total.contrib_D2BA_BC_avg)),
  sum(rowSums(total.contrib_D2BB_BC_avg))
)
save(BC_avg, file="./BC/BC_avg.RData")
rm(total.contrib_D1_BC_avg, total.contrib_D2AA_BC_avg,
   total.contrib_D2AB_BC_avg, total.contrib_D2BA_BC_avg,
   total.contrib_D2BB_BC_avg, total.contrib_FA_BC_avg, 
   total.contrib_FB_BC_avg,total.contrib_GA_BC_avg, 
   total.contrib_GB_BC_avg)
rm(Kvals_contrib_D1_BC_avg, Kvals_contrib_D2AA_BC_avg, 
   Kvals_contrib_D2AB_BC_avg, Kvals_contrib_D2BA_BC_avg,
   Kvals_contrib_D2BB_BC_avg, Kvals_contrib_FA_BC_avg,
   Kvals_contrib_FB_BC_avg, Kvals_contrib_GA_BC_avg,
   Kvals_contrib_GB_BC_avg)

#(b) Pool
load("./BC/Kvals_contrib_D1_BC_pool.RData")
load("./BC/Kvals_contrib_D2AA_BC_pool.RData")
load("./BC/Kvals_contrib_D2AB_BC_pool.RData")
load("./BC/Kvals_contrib_D2BA_BC_pool.RData")
load("./BC/Kvals_contrib_D2BB_BC_pool.RData")
load("./BC/Kvals_contrib_FA_BC_pool.RData")
load("./BC/Kvals_contrib_FB_BC_pool.RData")
load("./BC/Kvals_contrib_GA_BC_pool.RData")
load("./BC/Kvals_contrib_GB_BC_pool.RData")

total.contrib_D1_BC_pool<-apply(Kvals_contrib_D1_BC_pool, c(3,4), sum);
total.contrib_D2AA_BC_pool<-apply(Kvals_contrib_D2AA_BC_pool, c(3,4), sum);
total.contrib_D2AB_BC_pool<-apply(Kvals_contrib_D2AB_BC_pool, c(3,4),sum)
total.contrib_D2BA_BC_pool<-apply(Kvals_contrib_D2BA_BC_pool, c(3,4),sum)
total.contrib_D2BB_BC_pool<-apply(Kvals_contrib_D2BB_BC_pool, c(3,4), sum)
total.contrib_FA_BC_pool<-apply(Kvals_contrib_FA_BC_pool, c(3,4), sum);
total.contrib_FB_BC_pool<-apply(Kvals_contrib_FB_BC_pool, c(3,4), sum)
total.contrib_GA_BC_pool<-apply(Kvals_contrib_GA_BC_pool, c(3,4), sum);
total.contrib_GB_BC_pool<-apply(Kvals_contrib_GB_BC_pool, c(3,4), sum)

save(total.contrib_D1_BC_pool, total.contrib_D2AA_BC_pool,
     total.contrib_D2AB_BC_pool, total.contrib_D2BA_BC_pool,
     total.contrib_D2BB_BC_pool, total.contrib_FA_BC_pool, 
     total.contrib_FB_BC_pool,total.contrib_GA_BC_pool, 
     total.contrib_GB_BC_pool,
     file="./BC/total.contrib_BC_pool.RData")

BC_pool<-c( 
  sum(rowSums(total.contrib_D1_BC_pool)), 
  sum(rowSums(total.contrib_GA_BC_pool)),
  sum(rowSums(total.contrib_GB_BC_pool)),
  sum(rowSums(total.contrib_FA_BC_pool)),
  sum(rowSums(total.contrib_FB_BC_pool)),
  sum(rowSums(total.contrib_D2AA_BC_pool)),
  sum(rowSums(total.contrib_D2AB_BC_pool)),
  sum(rowSums(total.contrib_D2BA_BC_pool)),
  sum(rowSums(total.contrib_D2BB_BC_pool))
)
save(BC_pool, file="./BC/BC_pool.RData")
rm(total.contrib_D1_BC_pool, total.contrib_D2AA_BC_pool,
   total.contrib_D2AB_BC_pool, total.contrib_D2BA_BC_pool,
   total.contrib_D2BB_BC_pool, total.contrib_FA_BC_pool, 
   total.contrib_FB_BC_pool,total.contrib_GA_BC_pool, 
   total.contrib_GB_BC_pool)
rm(Kvals_contrib_D1_BC_pool, Kvals_contrib_D2AA_BC_pool, 
   Kvals_contrib_D2AB_BC_pool, Kvals_contrib_D2BA_BC_pool,
   Kvals_contrib_D2BB_BC_pool, Kvals_contrib_FA_BC_pool,
   Kvals_contrib_FB_BC_pool, Kvals_contrib_GA_BC_pool,
   Kvals_contrib_GB_BC_pool)


#Cape Canaveral

load("./CC/Kvals_contrib_D1_CC_avg.RData")
load("./CC/Kvals_contrib_D2AA_CC_avg.RData")
load("./CC/Kvals_contrib_D2AB_CC_avg.RData")
load("./CC/Kvals_contrib_D2BA_CC_avg.RData")
load("./CC/Kvals_contrib_D2BB_CC_avg.RData")
load("./CC/Kvals_contrib_FA_CC_avg.RData")
load("./CC/Kvals_contrib_FB_CC_avg.RData")
load("./CC/Kvals_contrib_GA_CC_avg.RData")
load("./CC/Kvals_contrib_GB_CC_avg.RData")

total.contrib_D1_CC_avg<-apply(Kvals_contrib_D1_CC_avg, c(3,4), sum);
total.contrib_D2AA_CC_avg<-apply(Kvals_contrib_D2AA_CC_avg, c(3,4), sum);
total.contrib_D2AB_CC_avg<-apply(Kvals_contrib_D2AB_CC_avg, c(3,4),sum)
total.contrib_D2BA_CC_avg<-apply(Kvals_contrib_D2BA_CC_avg, c(3,4),sum)
total.contrib_D2BB_CC_avg<-apply(Kvals_contrib_D2BB_CC_avg, c(3,4), sum)
total.contrib_FA_CC_avg<-apply(Kvals_contrib_FA_CC_avg, c(3,4), sum);
total.contrib_FB_CC_avg<-apply(Kvals_contrib_FB_CC_avg, c(3,4), sum)
total.contrib_GA_CC_avg<-apply(Kvals_contrib_GA_CC_avg, c(3,4), sum);
total.contrib_GB_CC_avg<-apply(Kvals_contrib_GB_CC_avg, c(3,4), sum)

save(total.contrib_D1_CC_avg, total.contrib_D2AA_CC_avg,
     total.contrib_D2AB_CC_avg, total.contrib_D2BA_CC_avg,
     total.contrib_D2BB_CC_avg, total.contrib_FA_CC_avg, 
     total.contrib_FB_CC_avg,total.contrib_GA_CC_avg, 
     total.contrib_GB_CC_avg,
     file="./CC/total.contrib_CC_avg.RData")

CC_avg<-c( 
  sum(rowSums(total.contrib_D1_CC_avg)), 
  sum(rowSums(total.contrib_GA_CC_avg)),
  sum(rowSums(total.contrib_GB_CC_avg)),
  sum(rowSums(total.contrib_FA_CC_avg)),
  sum(rowSums(total.contrib_FB_CC_avg)),
  sum(rowSums(total.contrib_D2AA_CC_avg)),
  sum(rowSums(total.contrib_D2AB_CC_avg)),
  sum(rowSums(total.contrib_D2BA_CC_avg)),
  sum(rowSums(total.contrib_D2BB_CC_avg))
)
save(CC_avg, file="./CC/CC_avg.RData")
rm(total.contrib_D1_CC_avg, total.contrib_D2AA_CC_avg,
   total.contrib_D2AB_CC_avg, total.contrib_D2BA_CC_avg,
   total.contrib_D2BB_CC_avg, total.contrib_FA_CC_avg, 
   total.contrib_FB_CC_avg,total.contrib_GA_CC_avg, 
   total.contrib_GB_CC_avg)
rm(Kvals_contrib_D1_CC_avg, Kvals_contrib_D2AA_CC_avg, 
   Kvals_contrib_D2AB_CC_avg, Kvals_contrib_D2BA_CC_avg,
   Kvals_contrib_D2BB_CC_avg, Kvals_contrib_FA_CC_avg,
   Kvals_contrib_FB_CC_avg, Kvals_contrib_GA_CC_avg,
   Kvals_contrib_GB_CC_avg)

#(b) Pool
load("./CC/Kvals_contrib_D1_CC_pool.RData")
load("./CC/Kvals_contrib_D2AA_CC_pool.RData")
load("./CC/Kvals_contrib_D2AB_CC_pool.RData")
load("./CC/Kvals_contrib_D2BA_CC_pool.RData")
load("./CC/Kvals_contrib_D2BB_CC_pool.RData")
load("./CC/Kvals_contrib_FA_CC_pool.RData")
load("./CC/Kvals_contrib_FB_CC_pool.RData")
load("./CC/Kvals_contrib_GA_CC_pool.RData")
load("./CC/Kvals_contrib_GB_CC_pool.RData")

total.contrib_D1_CC_pool<-apply(Kvals_contrib_D1_CC_pool, c(3,4), sum);
total.contrib_D2AA_CC_pool<-apply(Kvals_contrib_D2AA_CC_pool, c(3,4), sum);
total.contrib_D2AB_CC_pool<-apply(Kvals_contrib_D2AB_CC_pool, c(3,4),sum)
total.contrib_D2BA_CC_pool<-apply(Kvals_contrib_D2BA_CC_pool, c(3,4),sum)
total.contrib_D2BB_CC_pool<-apply(Kvals_contrib_D2BB_CC_pool, c(3,4), sum)
total.contrib_FA_CC_pool<-apply(Kvals_contrib_FA_CC_pool, c(3,4), sum);
total.contrib_FB_CC_pool<-apply(Kvals_contrib_FB_CC_pool, c(3,4), sum)
total.contrib_GA_CC_pool<-apply(Kvals_contrib_GA_CC_pool, c(3,4), sum);
total.contrib_GB_CC_pool<-apply(Kvals_contrib_GB_CC_pool, c(3,4), sum)

save(total.contrib_D1_CC_pool, total.contrib_D2AA_CC_pool,
     total.contrib_D2AB_CC_pool, total.contrib_D2BA_CC_pool,
     total.contrib_D2BB_CC_pool, total.contrib_FA_CC_pool, 
     total.contrib_FB_CC_pool,total.contrib_GA_CC_pool, 
     total.contrib_GB_CC_pool,
     file="./CC/total.contrib_CC_pool.RData")

CC_pool<-c( 
  sum(rowSums(total.contrib_D1_CC_pool)), 
  sum(rowSums(total.contrib_GA_CC_pool)),
  sum(rowSums(total.contrib_GB_CC_pool)),
  sum(rowSums(total.contrib_FA_CC_pool)),
  sum(rowSums(total.contrib_FB_CC_pool)),
  sum(rowSums(total.contrib_D2AA_CC_pool)),
  sum(rowSums(total.contrib_D2AB_CC_pool)),
  sum(rowSums(total.contrib_D2BA_CC_pool)),
  sum(rowSums(total.contrib_D2BB_CC_pool))
)
save(CC_pool, file="./CC/CC_pool.RData")
rm(total.contrib_D1_CC_pool, total.contrib_D2AA_CC_pool,
   total.contrib_D2AB_CC_pool, total.contrib_D2BA_CC_pool,
   total.contrib_D2BB_CC_pool, total.contrib_FA_CC_pool, 
   total.contrib_FB_CC_pool,total.contrib_GA_CC_pool, 
   total.contrib_GB_CC_pool)
rm(Kvals_contrib_D1_CC_pool, Kvals_contrib_D2AA_CC_pool, 
   Kvals_contrib_D2AB_CC_pool, Kvals_contrib_D2BA_CC_pool,
   Kvals_contrib_D2BB_CC_pool, Kvals_contrib_FA_CC_pool,
   Kvals_contrib_FB_CC_pool, Kvals_contrib_GA_CC_pool,
   Kvals_contrib_GB_CC_pool)


#Chekika
load("./C/Kvals_contrib_D1_C_avg.RData")
load("./C/Kvals_contrib_D2AA_C_avg.RData")
load("./C/Kvals_contrib_D2AB_C_avg.RData")
load("./C/Kvals_contrib_D2BA_C_avg.RData")
load("./C/Kvals_contrib_D2BB_C_avg.RData")
load("./C/Kvals_contrib_FA_C_avg.RData")
load("./C/Kvals_contrib_FB_C_avg.RData")
load("./C/Kvals_contrib_GA_C_avg.RData")
load("./C/Kvals_contrib_GB_C_avg.RData")

total.contrib_D1_C_avg<-apply(Kvals_contrib_D1_C_avg, c(3,4), sum);

total.contrib_D2AA_C_avg<-apply(Kvals_contrib_D2AA_C_avg, c(3,4), sum);
total.contrib_D2AB_C_avg<-apply(Kvals_contrib_D2AB_C_avg, c(3,4),sum)
total.contrib_D2BA_C_avg<-apply(Kvals_contrib_D2BA_C_avg, c(3,4),sum)
total.contrib_D2BB_C_avg<-apply(Kvals_contrib_D2BB_C_avg, c(3,4), sum)
total.contrib_FA_C_avg<-apply(Kvals_contrib_FA_C_avg, c(3,4), sum);
total.contrib_FB_C_avg<-apply(Kvals_contrib_FB_C_avg, c(3,4), sum)
total.contrib_GA_C_avg<-apply(Kvals_contrib_GA_C_avg, c(3,4), sum);
total.contrib_GB_C_avg<-apply(Kvals_contrib_GB_C_avg, c(3,4), sum)

save(total.contrib_D1_C_avg, total.contrib_D2AA_C_avg,
     total.contrib_D2AB_C_avg, total.contrib_D2BA_C_avg,
     total.contrib_D2BB_C_avg, total.contrib_FA_C_avg, 
     total.contrib_FB_C_avg,total.contrib_GA_C_avg, 
     total.contrib_GB_C_avg,
     file="./C/total.contrib_C_avg.RData")

C_avg<-c( 
  sum(rowSums(total.contrib_D1_C_avg)), 
  sum(rowSums(total.contrib_GA_C_avg)),
  sum(rowSums(total.contrib_GB_C_avg)),
  sum(rowSums(total.contrib_FA_C_avg)),
  sum(rowSums(total.contrib_FB_C_avg)),
  sum(rowSums(total.contrib_D2AA_C_avg)),
  sum(rowSums(total.contrib_D2AB_C_avg)),
  sum(rowSums(total.contrib_D2BA_C_avg)),
  sum(rowSums(total.contrib_D2BB_C_avg))
)
save(C_avg, file="./C/C_avg.RData")
rm(total.contrib_D1_C_avg, total.contrib_D2AA_C_avg,
   total.contrib_D2AB_C_avg, total.contrib_D2BA_C_avg,
   total.contrib_D2BB_C_avg, total.contrib_FA_C_avg, 
   total.contrib_FB_C_avg,total.contrib_GA_C_avg, 
   total.contrib_GB_C_avg)
rm(Kvals_contrib_D1_C_avg, Kvals_contrib_D2AA_C_avg, 
   Kvals_contrib_D2AB_C_avg, Kvals_contrib_D2BA_C_avg,
   Kvals_contrib_D2BB_C_avg, Kvals_contrib_FA_C_avg,
   Kvals_contrib_FB_C_avg, Kvals_contrib_GA_C_avg,
   Kvals_contrib_GB_C_avg)

#(b) Pool
load("./C/Kvals_contrib_D1_C_pool.RData")
load("./C/Kvals_contrib_D2AA_C_pool.RData")
load("./C/Kvals_contrib_D2AB_C_pool.RData")
load("./C/Kvals_contrib_D2BA_C_pool.RData")
load("./C/Kvals_contrib_D2BB_C_pool.RData")
load("./C/Kvals_contrib_FA_C_pool.RData")
load("./C/Kvals_contrib_FB_C_pool.RData")
load("./C/Kvals_contrib_GA_C_pool.RData")
load("./C/Kvals_contrib_GB_C_pool.RData")

total.contrib_D1_C_pool<-apply(Kvals_contrib_D1_C_pool, c(3,4), sum);
total.contrib_D2AA_C_pool<-apply(Kvals_contrib_D2AA_C_pool, c(3,4), sum);
total.contrib_D2AB_C_pool<-apply(Kvals_contrib_D2AB_C_pool, c(3,4),sum)
total.contrib_D2BA_C_pool<-apply(Kvals_contrib_D2BA_C_pool, c(3,4),sum)
total.contrib_D2BB_C_pool<-apply(Kvals_contrib_D2BB_C_pool, c(3,4), sum)
total.contrib_FA_C_pool<-apply(Kvals_contrib_FA_C_pool, c(3,4), sum);
total.contrib_FB_C_pool<-apply(Kvals_contrib_FB_C_pool, c(3,4), sum)
total.contrib_GA_C_pool<-apply(Kvals_contrib_GA_C_pool, c(3,4), sum);
total.contrib_GB_C_pool<-apply(Kvals_contrib_GB_C_pool, c(3,4), sum)

save(total.contrib_D1_C_pool, total.contrib_D2AA_C_pool,
     total.contrib_D2AB_C_pool, total.contrib_D2BA_C_pool,
     total.contrib_D2BB_C_pool, total.contrib_FA_C_pool, 
     total.contrib_FB_C_pool,total.contrib_GA_C_pool, 
     total.contrib_GB_C_pool,
     file="./C/total.contrib_C_pool.RData")

C_pool<-c( 
  sum(rowSums(total.contrib_D1_C_pool)), 
  sum(rowSums(total.contrib_GA_C_pool)),
  sum(rowSums(total.contrib_GB_C_pool)),
  sum(rowSums(total.contrib_FA_C_pool)),
  sum(rowSums(total.contrib_FB_C_pool)),
  sum(rowSums(total.contrib_D2AA_C_pool)),
  sum(rowSums(total.contrib_D2AB_C_pool)),
  sum(rowSums(total.contrib_D2BA_C_pool)),
  sum(rowSums(total.contrib_D2BB_C_pool))
)
save(C_pool, file="./C/C_pool.RData")
rm(total.contrib_D1_C_pool, total.contrib_D2AA_C_pool,
   total.contrib_D2AB_C_pool, total.contrib_D2BA_C_pool,
   total.contrib_D2BB_C_pool, total.contrib_FA_C_pool, 
   total.contrib_FB_C_pool,total.contrib_GA_C_pool, 
   total.contrib_GB_C_pool)
rm(Kvals_contrib_D1_C_pool, Kvals_contrib_D2AA_C_pool, 
   Kvals_contrib_D2AB_C_pool, Kvals_contrib_D2BA_C_pool,
   Kvals_contrib_D2BB_C_pool, Kvals_contrib_FA_C_pool,
   Kvals_contrib_FB_C_pool, Kvals_contrib_GA_C_pool,
   Kvals_contrib_GB_C_pool)

#Fort Pierce
load("./FP/Kvals_contrib_D1_FP_avg.RData")
load("./FP/Kvals_contrib_D2AA_FP_avg.RData")
load("./FP/Kvals_contrib_D2AB_FP_avg.RData")
load("./FP/Kvals_contrib_D2BA_FP_avg.RData")
load("./FP/Kvals_contrib_D2BB_FP_avg.RData")
load("./FP/Kvals_contrib_FA_FP_avg.RData")
load("./FP/Kvals_contrib_FB_FP_avg.RData")
load("./FP/Kvals_contrib_GA_FP_avg.RData")
load("./FP/Kvals_contrib_GB_FP_avg.RData")

total.contrib_D1_FP_avg<-apply(Kvals_contrib_D1_FP_avg, c(3,4), sum);

total.contrib_D2AA_FP_avg<-apply(Kvals_contrib_D2AA_FP_avg, c(3,4), sum);
total.contrib_D2AB_FP_avg<-apply(Kvals_contrib_D2AB_FP_avg, c(3,4),sum)
total.contrib_D2BA_FP_avg<-apply(Kvals_contrib_D2BA_FP_avg, c(3,4),sum)
total.contrib_D2BB_FP_avg<-apply(Kvals_contrib_D2BB_FP_avg, c(3,4), sum)
total.contrib_FA_FP_avg<-apply(Kvals_contrib_FA_FP_avg, c(3,4), sum);
total.contrib_FB_FP_avg<-apply(Kvals_contrib_FB_FP_avg, c(3,4), sum)
total.contrib_GA_FP_avg<-apply(Kvals_contrib_GA_FP_avg, c(3,4), sum);
total.contrib_GB_FP_avg<-apply(Kvals_contrib_GB_FP_avg, c(3,4), sum)

save(total.contrib_D1_FP_avg, total.contrib_D2AA_FP_avg,
     total.contrib_D2AB_FP_avg, total.contrib_D2BA_FP_avg,
     total.contrib_D2BB_FP_avg, total.contrib_FA_FP_avg, 
     total.contrib_FB_FP_avg,total.contrib_GA_FP_avg, 
     total.contrib_GB_FP_avg,
     file="./FP/total.contrib_FP_avg.RData")

FP_avg<-c( 
  sum(rowSums(total.contrib_D1_FP_avg)), 
  sum(rowSums(total.contrib_GA_FP_avg)),
  sum(rowSums(total.contrib_GB_FP_avg)),
  sum(rowSums(total.contrib_FA_FP_avg)),
  sum(rowSums(total.contrib_FB_FP_avg)),
  sum(rowSums(total.contrib_D2AA_FP_avg)),
  sum(rowSums(total.contrib_D2AB_FP_avg)),
  sum(rowSums(total.contrib_D2BA_FP_avg)),
  sum(rowSums(total.contrib_D2BB_FP_avg))
)
save(FP_avg, file="./FP/FP_avg.RData")
rm(total.contrib_D1_FP_avg, total.contrib_D2AA_FP_avg,
   total.contrib_D2AB_FP_avg, total.contrib_D2BA_FP_avg,
   total.contrib_D2BB_FP_avg, total.contrib_FA_FP_avg, 
   total.contrib_FB_FP_avg,total.contrib_GA_FP_avg, 
   total.contrib_GB_FP_avg)
rm(Kvals_contrib_D1_FP_avg, Kvals_contrib_D2AA_FP_avg, 
   Kvals_contrib_D2AB_FP_avg, Kvals_contrib_D2BA_FP_avg,
   Kvals_contrib_D2BB_FP_avg, Kvals_contrib_FA_FP_avg,
   Kvals_contrib_FB_FP_avg, Kvals_contrib_GA_FP_avg,
   Kvals_contrib_GB_FP_avg)

#(b) Pool
load("./FP/Kvals_contrib_D1_FP_pool.RData")
load("./FP/Kvals_contrib_D2AA_FP_pool.RData")
load("./FP/Kvals_contrib_D2AB_FP_pool.RData")
load("./FP/Kvals_contrib_D2BA_FP_pool.RData")
load("./FP/Kvals_contrib_D2BB_FP_pool.RData")
load("./FP/Kvals_contrib_FA_FP_pool.RData")
load("./FP/Kvals_contrib_FB_FP_pool.RData")
load("./FP/Kvals_contrib_GA_FP_pool.RData")
load("./FP/Kvals_contrib_GB_FP_pool.RData")

total.contrib_D1_FP_pool<-apply(Kvals_contrib_D1_FP_pool, c(3,4), sum);
total.contrib_D2AA_FP_pool<-apply(Kvals_contrib_D2AA_FP_pool, c(3,4), sum);
total.contrib_D2AB_FP_pool<-apply(Kvals_contrib_D2AB_FP_pool, c(3,4),sum)
total.contrib_D2BA_FP_pool<-apply(Kvals_contrib_D2BA_FP_pool, c(3,4),sum)
total.contrib_D2BB_FP_pool<-apply(Kvals_contrib_D2BB_FP_pool, c(3,4), sum)
total.contrib_FA_FP_pool<-apply(Kvals_contrib_FA_FP_pool, c(3,4), sum);
total.contrib_FB_FP_pool<-apply(Kvals_contrib_FB_FP_pool, c(3,4), sum)
total.contrib_GA_FP_pool<-apply(Kvals_contrib_GA_FP_pool, c(3,4), sum);
total.contrib_GB_FP_pool<-apply(Kvals_contrib_GB_FP_pool, c(3,4), sum)

save(total.contrib_D1_FP_pool, total.contrib_D2AA_FP_pool,
     total.contrib_D2AB_FP_pool, total.contrib_D2BA_FP_pool,
     total.contrib_D2BB_FP_pool, total.contrib_FA_FP_pool, 
     total.contrib_FB_FP_pool,total.contrib_GA_FP_pool, 
     total.contrib_GB_FP_pool,
     file="./FP/total.contrib_FP_pool.RData")

FP_pool<-c( 
  sum(rowSums(total.contrib_D1_FP_pool)), 
  sum(rowSums(total.contrib_GA_FP_pool)),
  sum(rowSums(total.contrib_GB_FP_pool)),
  sum(rowSums(total.contrib_FA_FP_pool)),
  sum(rowSums(total.contrib_FB_FP_pool)),
  sum(rowSums(total.contrib_D2AA_FP_pool)),
  sum(rowSums(total.contrib_D2AB_FP_pool)),
  sum(rowSums(total.contrib_D2BA_FP_pool)),
  sum(rowSums(total.contrib_D2BB_FP_pool))
)
save(FP_pool, file="./FP/FP_pool.RData")
rm(total.contrib_D1_FP_pool, total.contrib_D2AA_FP_pool,
   total.contrib_D2AB_FP_pool, total.contrib_D2BA_FP_pool,
   total.contrib_D2BB_FP_pool, total.contrib_FA_FP_pool, 
   total.contrib_FB_FP_pool,total.contrib_GA_FP_pool, 
   total.contrib_GB_FP_pool)
rm(Kvals_contrib_D1_FP_pool, Kvals_contrib_D2AA_FP_pool, 
   Kvals_contrib_D2AB_FP_pool, Kvals_contrib_D2BA_FP_pool,
   Kvals_contrib_D2BB_FP_pool, Kvals_contrib_FA_FP_pool,
   Kvals_contrib_FB_FP_pool, Kvals_contrib_GA_FP_pool,
   Kvals_contrib_GB_FP_pool)



#Punta Gorda
load("./PG/Kvals_contrib_D1_PG_avg.RData")
load("./PG/Kvals_contrib_D2AA_PG_avg.RData")
load("./PG/Kvals_contrib_D2AB_PG_avg.RData")
load("./PG/Kvals_contrib_D2BA_PG_avg.RData")
load("./PG/Kvals_contrib_D2BB_PG_avg.RData")
load("./PG/Kvals_contrib_FA_PG_avg.RData")
load("./PG/Kvals_contrib_FB_PG_avg.RData")
load("./PG/Kvals_contrib_GA_PG_avg.RData")
load("./PG/Kvals_contrib_GB_PG_avg.RData")

total.contrib_D1_PG_avg<-apply(Kvals_contrib_D1_PG_avg, c(3,4), sum);

total.contrib_D2AA_PG_avg<-apply(Kvals_contrib_D2AA_PG_avg, c(3,4), sum);
total.contrib_D2AB_PG_avg<-apply(Kvals_contrib_D2AB_PG_avg, c(3,4),sum)
total.contrib_D2BA_PG_avg<-apply(Kvals_contrib_D2BA_PG_avg, c(3,4),sum)
total.contrib_D2BB_PG_avg<-apply(Kvals_contrib_D2BB_PG_avg, c(3,4), sum)
total.contrib_FA_PG_avg<-apply(Kvals_contrib_FA_PG_avg, c(3,4), sum);
total.contrib_FB_PG_avg<-apply(Kvals_contrib_FB_PG_avg, c(3,4), sum)
total.contrib_GA_PG_avg<-apply(Kvals_contrib_GA_PG_avg, c(3,4), sum);
total.contrib_GB_PG_avg<-apply(Kvals_contrib_GB_PG_avg, c(3,4), sum)

save(total.contrib_D1_PG_avg, total.contrib_D2AA_PG_avg,
     total.contrib_D2AB_PG_avg, total.contrib_D2BA_PG_avg,
     total.contrib_D2BB_PG_avg, total.contrib_FA_PG_avg, 
     total.contrib_FB_PG_avg,total.contrib_GA_PG_avg, 
     total.contrib_GB_PG_avg,
     file="./PG/total.contrib_PG_avg.RData")

PG_avg<-c( 
  sum(rowSums(total.contrib_D1_PG_avg)), 
  sum(rowSums(total.contrib_GA_PG_avg)),
  sum(rowSums(total.contrib_GB_PG_avg)),
  sum(rowSums(total.contrib_FA_PG_avg)),
  sum(rowSums(total.contrib_FB_PG_avg)),
  sum(rowSums(total.contrib_D2AA_PG_avg)),
  sum(rowSums(total.contrib_D2AB_PG_avg)),
  sum(rowSums(total.contrib_D2BA_PG_avg)),
  sum(rowSums(total.contrib_D2BB_PG_avg))
)
save(PG_avg, file="./PG/PG_avg.RData")
rm(total.contrib_D1_PG_avg, total.contrib_D2AA_PG_avg,
   total.contrib_D2AB_PG_avg, total.contrib_D2BA_PG_avg,
   total.contrib_D2BB_PG_avg, total.contrib_FA_PG_avg, 
   total.contrib_FB_PG_avg,total.contrib_GA_PG_avg, 
   total.contrib_GB_PG_avg)
rm(Kvals_contrib_D1_PG_avg, Kvals_contrib_D2AA_PG_avg, 
   Kvals_contrib_D2AB_PG_avg, Kvals_contrib_D2BA_PG_avg,
   Kvals_contrib_D2BB_PG_avg, Kvals_contrib_FA_PG_avg,
   Kvals_contrib_FB_PG_avg, Kvals_contrib_GA_PG_avg,
   Kvals_contrib_GB_PG_avg)

#(b) Pool
load("./PG/Kvals_contrib_D1_PG_pool.RData")
load("./PG/Kvals_contrib_D2AA_PG_pool.RData")
load("./PG/Kvals_contrib_D2AB_PG_pool.RData")
load("./PG/Kvals_contrib_D2BA_PG_pool.RData")
load("./PG/Kvals_contrib_D2BB_PG_pool.RData")
load("./PG/Kvals_contrib_FA_PG_pool.RData")
load("./PG/Kvals_contrib_FB_PG_pool.RData")
load("./PG/Kvals_contrib_GA_PG_pool.RData")
load("./PG/Kvals_contrib_GB_PG_pool.RData")

total.contrib_D1_PG_pool<-apply(Kvals_contrib_D1_PG_pool, c(3,4), sum);
total.contrib_D2AA_PG_pool<-apply(Kvals_contrib_D2AA_PG_pool, c(3,4), sum);
total.contrib_D2AB_PG_pool<-apply(Kvals_contrib_D2AB_PG_pool, c(3,4),sum)
total.contrib_D2BA_PG_pool<-apply(Kvals_contrib_D2BA_PG_pool, c(3,4),sum)
total.contrib_D2BB_PG_pool<-apply(Kvals_contrib_D2BB_PG_pool, c(3,4), sum)
total.contrib_FA_PG_pool<-apply(Kvals_contrib_FA_PG_pool, c(3,4), sum);
total.contrib_FB_PG_pool<-apply(Kvals_contrib_FB_PG_pool, c(3,4), sum)
total.contrib_GA_PG_pool<-apply(Kvals_contrib_GA_PG_pool, c(3,4), sum);
total.contrib_GB_PG_pool<-apply(Kvals_contrib_GB_PG_pool, c(3,4), sum)

save(total.contrib_D1_PG_pool, total.contrib_D2AA_PG_pool,
     total.contrib_D2AB_PG_pool, total.contrib_D2BA_PG_pool,
     total.contrib_D2BB_PG_pool, total.contrib_FA_PG_pool, 
     total.contrib_FB_PG_pool,total.contrib_GA_PG_pool, 
     total.contrib_GB_PG_pool,
     file="./PG/total.contrib_PG_pool.RData")

PG_pool<-c( 
  sum(rowSums(total.contrib_D1_PG_pool)), 
  sum(rowSums(total.contrib_GA_PG_pool)),
  sum(rowSums(total.contrib_GB_PG_pool)),
  sum(rowSums(total.contrib_FA_PG_pool)),
  sum(rowSums(total.contrib_FB_PG_pool)),
  sum(rowSums(total.contrib_D2AA_PG_pool)),
  sum(rowSums(total.contrib_D2AB_PG_pool)),
  sum(rowSums(total.contrib_D2BA_PG_pool)),
  sum(rowSums(total.contrib_D2BB_PG_pool))
)
save(PG_pool, file="./PG/PG_pool.RData")
rm(total.contrib_D1_PG_pool, total.contrib_D2AA_PG_pool,
   total.contrib_D2AB_PG_pool, total.contrib_D2BA_PG_pool,
   total.contrib_D2BB_PG_pool, total.contrib_FA_PG_pool, 
   total.contrib_FB_PG_pool,total.contrib_GA_PG_pool, 
   total.contrib_GB_PG_pool)
rm(Kvals_contrib_D1_PG_pool, Kvals_contrib_D2AA_PG_pool, 
   Kvals_contrib_D2AB_PG_pool, Kvals_contrib_D2BA_PG_pool,
   Kvals_contrib_D2BB_PG_pool, Kvals_contrib_FA_PG_pool,
   Kvals_contrib_FB_PG_pool, Kvals_contrib_GA_PG_pool,
   Kvals_contrib_GB_PG_pool)

#Wild Turkey
load("./WT/Kvals_contrib_D1_WT_avg.RData")
load("./WT/Kvals_contrib_D2AA_WT_avg.RData")
load("./WT/Kvals_contrib_D2AB_WT_avg.RData")
load("./WT/Kvals_contrib_D2BA_WT_avg.RData")
load("./WT/Kvals_contrib_D2BB_WT_avg.RData")
load("./WT/Kvals_contrib_FA_WT_avg.RData")
load("./WT/Kvals_contrib_FB_WT_avg.RData")
load("./WT/Kvals_contrib_GA_WT_avg.RData")
load("./WT/Kvals_contrib_GB_WT_avg.RData")

total.contrib_D1_WT_avg<-apply(Kvals_contrib_D1_WT_avg, c(3,4), sum);

total.contrib_D2AA_WT_avg<-apply(Kvals_contrib_D2AA_WT_avg, c(3,4), sum);
total.contrib_D2AB_WT_avg<-apply(Kvals_contrib_D2AB_WT_avg, c(3,4),sum)
total.contrib_D2BA_WT_avg<-apply(Kvals_contrib_D2BA_WT_avg, c(3,4),sum)
total.contrib_D2BB_WT_avg<-apply(Kvals_contrib_D2BB_WT_avg, c(3,4), sum)
total.contrib_FA_WT_avg<-apply(Kvals_contrib_FA_WT_avg, c(3,4), sum);
total.contrib_FB_WT_avg<-apply(Kvals_contrib_FB_WT_avg, c(3,4), sum)
total.contrib_GA_WT_avg<-apply(Kvals_contrib_GA_WT_avg, c(3,4), sum);
total.contrib_GB_WT_avg<-apply(Kvals_contrib_GB_WT_avg, c(3,4), sum)

save(total.contrib_D1_WT_avg, total.contrib_D2AA_WT_avg,
     total.contrib_D2AB_WT_avg, total.contrib_D2BA_WT_avg,
     total.contrib_D2BB_WT_avg, total.contrib_FA_WT_avg, 
     total.contrib_FB_WT_avg,total.contrib_GA_WT_avg, 
     total.contrib_GB_WT_avg,
     file="./WT/total.contrib_WT_avg.RData")

WT_avg<-c( 
  sum(rowSums(total.contrib_D1_WT_avg)), 
  sum(rowSums(total.contrib_GA_WT_avg)),
  sum(rowSums(total.contrib_GB_WT_avg)),
  sum(rowSums(total.contrib_FA_WT_avg)),
  sum(rowSums(total.contrib_FB_WT_avg)),
  sum(rowSums(total.contrib_D2AA_WT_avg)),
  sum(rowSums(total.contrib_D2AB_WT_avg)),
  sum(rowSums(total.contrib_D2BA_WT_avg)),
  sum(rowSums(total.contrib_D2BB_WT_avg))
)
save(WT_avg, file="./WT/WT_avg.RData")
rm(total.contrib_D1_WT_avg, total.contrib_D2AA_WT_avg,
   total.contrib_D2AB_WT_avg, total.contrib_D2BA_WT_avg,
   total.contrib_D2BB_WT_avg, total.contrib_FA_WT_avg, 
   total.contrib_FB_WT_avg,total.contrib_GA_WT_avg, 
   total.contrib_GB_WT_avg)
rm(Kvals_contrib_D1_WT_avg, Kvals_contrib_D2AA_WT_avg, 
   Kvals_contrib_D2AB_WT_avg, Kvals_contrib_D2BA_WT_avg,
   Kvals_contrib_D2BB_WT_avg, Kvals_contrib_FA_WT_avg,
   Kvals_contrib_FB_WT_avg, Kvals_contrib_GA_WT_avg,
   Kvals_contrib_GB_WT_avg)

#(b) Pool
load("./WT/Kvals_contrib_D1_WT_pool.RData")
load("./WT/Kvals_contrib_D2AA_WT_pool.RData")
load("./WT/Kvals_contrib_D2AB_WT_pool.RData")
load("./WT/Kvals_contrib_D2BA_WT_pool.RData")
load("./WT/Kvals_contrib_D2BB_WT_pool.RData")
load("./WT/Kvals_contrib_FA_WT_pool.RData")
load("./WT/Kvals_contrib_FB_WT_pool.RData")
load("./WT/Kvals_contrib_GA_WT_pool.RData")
load("./WT/Kvals_contrib_GB_WT_pool.RData")

total.contrib_D1_WT_pool<-apply(Kvals_contrib_D1_WT_pool, c(3,4), sum);
total.contrib_D2AA_WT_pool<-apply(Kvals_contrib_D2AA_WT_pool, c(3,4), sum);
total.contrib_D2AB_WT_pool<-apply(Kvals_contrib_D2AB_WT_pool, c(3,4),sum)
total.contrib_D2BA_WT_pool<-apply(Kvals_contrib_D2BA_WT_pool, c(3,4),sum)
total.contrib_D2BB_WT_pool<-apply(Kvals_contrib_D2BB_WT_pool, c(3,4), sum)
total.contrib_FA_WT_pool<-apply(Kvals_contrib_FA_WT_pool, c(3,4), sum);
total.contrib_FB_WT_pool<-apply(Kvals_contrib_FB_WT_pool, c(3,4), sum)
total.contrib_GA_WT_pool<-apply(Kvals_contrib_GA_WT_pool, c(3,4), sum);
total.contrib_GB_WT_pool<-apply(Kvals_contrib_GB_WT_pool, c(3,4), sum)

save(total.contrib_D1_WT_pool, total.contrib_D2AA_WT_pool,
     total.contrib_D2AB_WT_pool, total.contrib_D2BA_WT_pool,
     total.contrib_D2BB_WT_pool, total.contrib_FA_WT_pool, 
     total.contrib_FB_WT_pool,total.contrib_GA_WT_pool, 
     total.contrib_GB_WT_pool,
     file="./WT/total.contrib_WT_pool.RData")

WT_pool<-c( 
  sum(rowSums(total.contrib_D1_WT_pool)), 
  sum(rowSums(total.contrib_GA_WT_pool)),
  sum(rowSums(total.contrib_GB_WT_pool)),
  sum(rowSums(total.contrib_FA_WT_pool)),
  sum(rowSums(total.contrib_FB_WT_pool)),
  sum(rowSums(total.contrib_D2AA_WT_pool)),
  sum(rowSums(total.contrib_D2AB_WT_pool)),
  sum(rowSums(total.contrib_D2BA_WT_pool)),
  sum(rowSums(total.contrib_D2BB_WT_pool))
)
save(WT_pool, file="./WT/WT_pool.RData")
rm(total.contrib_D1_WT_pool, total.contrib_D2AA_WT_pool,
   total.contrib_D2AB_WT_pool, total.contrib_D2BA_WT_pool,
   total.contrib_D2BB_WT_pool, total.contrib_FA_WT_pool, 
   total.contrib_FB_WT_pool,total.contrib_GA_WT_pool, 
   total.contrib_GB_WT_pool)
rm(Kvals_contrib_D1_WT_pool, Kvals_contrib_D2AA_WT_pool, 
   Kvals_contrib_D2AB_WT_pool, Kvals_contrib_D2BA_WT_pool,
   Kvals_contrib_D2BB_WT_pool, Kvals_contrib_FA_WT_pool,
   Kvals_contrib_FB_WT_pool, Kvals_contrib_GA_WT_pool,
   Kvals_contrib_GB_WT_pool)














 
# # Variability in lambda #####
#  #With six sites to compare, this part of the
#  #code takes a long time to run (over 48 hours at least)
# # Calculate matrices of coefficients of variation #####
# 
# 	
# D1_BC<-readRDS("./BC/D1_BC.rds")
# D1_PG <- readRDS("./PG/D1_PG.rds")
# D1_C<-readRDS("./C/D1_C.rds")
# D1_FP<-readRDS("./FP/D1_FP.rds")
# D1_PG<-readRDS("./PG/D1_PG.rds")
# D1_WT<-readRDS("./WT/D1_WT.rds")
# 	
# CV<-function(x) {
# 	return(sd(x)/mean(x))
# 	
# }
# 
# 
# 
# 	
# cv_D1<-matrix(nrow=(m1*m2), ncol=(m1*m2))
# 
# for(i in 1:(m1*m2)) {
# 	for(j in 1:(m1*m2)) {
# 		cv_D1[i,j]<-CV(c(D1_BC[i,j], D1_PG[i,j], D1_C[i,j], D1_FP[i,j], 
# 		                 D1_PG[i,j], D1_WT[i,j]))
# 		if(is.na(cv_D1[i,j])) {cv_D1[i,j]<-0}
# 	}
# }
# 
# 
# D1_avg=(D1_BC + D1_PG + D1_C + D1_FP + D1_PG + D1_WT)/6
# rm(D1_BC, D1_PG, D1_C, D1_FP, D1_PG, D1_WT)
# 
# #G
# G_BC<-readRDS("./BC/G_BC.rds")
# G_PG<-readRDS("./PG/G_PG.rds")
# G_C <- readRDS("./C/G_C.rds")
# G_FP <- readRDS("./FP/G_FP.rds")
# G_PG <- readRDS("./PG/G_PG.rds")
# G_WT <- readRDS("./WT/G_WT.rds")
# 
# cv_G<-matrix(nrow=(m3*m4), ncol=(m1*m2))
# 
# 
# for(i in 1:(m3*m4)) {
# 	for(j in 1:(m1*m2)) {
# 		cv_G[i,j]<-CV(c(G_BC[i,j], G_PG[i,j], G_C[i,j], G_FP[i,j],
# 		                G_PG[i,j], G_WT[i,j]))
# 		if(is.na(cv_G[i,j])) {cv_G[i,j]<-0}
# 	}
# }
# 
# G_avg<-(G_BC+G_PG+G_C + G_FP + G_PG + G_WT)/6
# 
# rm(G_BC, G_PG, G_C, G_FP, G_PG, G_WT)
# 
# #F
# F_BC <- readRDS("./BC/F_BC.rds")
# F_PG <- readRDS("./PG/F_PG.rds")
# F_C  <- readRDS("./C/F_C.rds")
# F_FP <- readRDS("./FP/F_FP.rds")
# F_PG <- readRDS("./PG/F_PG.rds")
# F_WT <- readRDS("./WT/F_WT.rds")
# 
# cv_F<-matrix(nrow=(m1*m2), ncol=(m3*m4))
# for(i in 1:(m1*m2)) {
# 	for(j in 1:(m3*m4)) {
# 		cv_F[i,j]<-CV(c(F_BC[i,j], F_PG[i,j], F_C[i,j], F_FP[i,j],
# 		              F_PG[i,j], F_WT[i,j]))
# 		if(is.na(cv_F[i,j])) {cv_F[i,j]<-0}
# 	}
# }
# 
# F_avg<-(F_BC + F_PG + F_C + F_FP + F_PG + F_WT)/6
# 
# rm(F_BC, F_PG, F_C, F_FP, F_PG, F_WT)
# 
# 
# save(cv_D1, file="cv_D1.RData")
# save(cv_F, file="cv_F.RData")
# save(cv_G, file="cv_G.RData")
# rm(cv_D1, cv_F, cv_G)
# #D2
# D2_BC <- readRDS("./BC/D2_BC.rds")
# D2_PG <- readRDS("./PG/D2_PG.rds")
# D2_C  <- readRDS("./C/D2_C.rds")
# D2_FP <- readRDS("./FP/D2_FP.rds")
# D2_PG <- readRDS("./PG/D2_PG.rds")
# D2_WT <- readRDS("./WT/D2_WT.rds")
# 
# #Lines 568: 577 are run on a separate computer because they take
# # a long time to run (X hours)
# cv_D2<-matrix(nrow=(m3*m4), ncol=(m3*m4))
# for(i in 1:(m3*m4)) {
# 	for(j in 1:(m3*m4)) {
# 		cv_D2[i,j]<-CV(c(D2_BC[i,j], D2_PG[i,j], D2_C[i,j], D2_FP,
# 		                 D2_PG, D2_WT))
# 		if(is.na(cv_D2[i,j])) {cv_D2[i,j]<-0}
# 	}
#   cat("i=", i, "\n")
# }
# save(cv_D2, file="cv_D2.RData")
# 
# D2_avg<-(D2_E+D2_H+D2_W)/3
# rm(D2_E)
# rm(D2_W)
# rm(D2_H)
# 
# 
# A_cv<-cbind(rbind(cv_D1, cv_G), rbind(cv_F, cv_D2))
# rm(cv_D1, cv_G, cv_F, cv_D2)
# 
# save(A_cv, file="./Overall/A_cv.RData")
# 
# load("./Overall/elas_overall.RData")
# thing<-A_cv*elas_overall
# 
# save(thing, file="./Overall/thing.RData")
# 
# #Break thing down into components ##### 
# 
# thing_D1<-thing[1:(m1*m2), 1:(m1*m2)]
# thing_G<-thing[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
# thing_F<-thing[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
# thing_D2<-thing[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]
# 
# #Unpack matrix of contributions into pieces #####
# 
# 	
# 	#D1#
# 	plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in 	array
# 	Plop=outer(1:m1,1:m2,plop); 
# 
# 	Kvals_thing_D1=array(0,c(m1,m2,m1,m2));  
# 
# 	for(i in 1:m1){
# 		for(j in 1:m2){
# 		for(k in 1:m1){
# 				kvals= thing_D1[Plop[k,1:m2],Plop[i,j]]
# 				Kvals_thing_D1[k,1:m2,i,j]=kvals
# 			
# 	}}
# 	cat(i,"\n"); 
# }
# 
# 
# ###Construct D2 (Large Domain):
# plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
# Plop=outer(1:m3,1:m4,plop); 
# 
# 
# Kvals_thing_D2=array(0,c(m3,m4,m3,m4));  
# 
# 
# 
# for(i in 1:m3){
# 	for(j in 1:m4){
# 		for(k in 1:m3){
# 				kvals= thing_D2[Plop[k,1:m4],Plop[i,j]]
# 				Kvals_thing_D2[k,1:m4,i,j]=kvals
# 			
# 	}}
# 	cat(i,"\n"); 
# }		
# 
# 
# ###Construct F (Fecundity):
# plop1=function(i, j) {(j-1)*m1 + i}
# plop2=function(i, j) {(j-1)*m3 + i}
# Plop1=outer(1:m1,1:m2,plop1); 
# Plop2=outer(1:m3, 1:m4, plop2);
# 
# Kvals_thing_F=array(0, c(m1, m2, m3, m4))
# 
# for(i in 1:m3) {
# 	for (j in 1:m4) {
# 		for (k in 1:m1) {
# 			kvals=thing_F[Plop1[k, 1:m2], Plop2[i,j]]
# 			Kvals_thing_F[k, 1:m2, i, j]=kvals
# 		}}
# 		cat(i, "\n");
# }
# 
# 
# ###Construct G (Graduation):
# plop1=function(i, j) {(j-1)*m3 + i}
# plop2=function(i, j) {(j-1)*m1 + i}
# Plop1=outer(1:m3,1:m4,plop1); 
# Plop2=outer(1:m1, 1:m2, plop2);
# 
# Kvals_thing_G=array(0, c(m3, m4, m1, m2))
# 
# for(i in 1:m1) {
# 	for (j in 1:m2) {
# 		for (k in 1:m3) {
# 			kvals=thing_G[Plop1[k, 1:m4], Plop2[i,j]]
# 			Kvals_thing_G[k, 1:m4, i, j]=kvals
# 		}}
# 		cat(i, "\n");
# }
# 
# 
# # Sum up total contributions #####
# total.thing_D1<-apply(Kvals_thing_D1, c(3,4), sum);
# total.thing_D2<-apply(Kvals_thing_D2, c(3,4), sum);
# total.thing_F<-apply(Kvals_thing_F, c(3,4), sum);
# total.thing_G<-apply(Kvals_thing_G, c(3,4), sum);
# 	
# save(total.thing_D1, total.thing_D2, total.thing_F, total.thing_G, 
#      file="./Overall/total.thing.RData")
# 
# 
# 
# 
# 
# 
# 
# 
