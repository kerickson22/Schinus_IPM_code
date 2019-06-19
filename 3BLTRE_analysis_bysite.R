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

#Big Cypress
thing <- decompose(C_BC)
saveRDS(thing$D1, file="./BC/C_D1_BC.rds")
saveRDS(thing$G, file="./BC/C_G_BC.rds")
saveRDS(thing$F, file="./BC/C_F_BC.rds")
saveRDS(thing$D2, file="./BC/C_D2_BC.rds")
#Cape Canaveral
thing <- decompose(C_CC)
saveRDS(thing$D1, file="./CC/C_D1_CC.rds")
saveRDS(thing$G, file="./CC/C_G_CC.rds")
saveRDS(thing$F, file="./CC/C_F_CC.rds")
saveRDS(thing$D2, file="./CC/C_D2_CC.rds")
#Chekika
thing <- decompose(C_C)
saveRDS(thing$D1, file="./C/C_D1_C.rds")
saveRDS(thing$G, file="./C/C_G_C.rds")
saveRDS(thing$F, file="./C/C_F_C.rds")
saveRDS(thing$D2, file="./C/C_D2_C.rds")
#Fort Pierce
thing <- decompose(C_FP)
saveRDS(thing$D1, file="./FP/C_D1_FP.rds")
saveRDS(thing$G, file="./FP/C_G_FP.rds")
saveRDS(thing$F, file="./FP/C_F_FP.rds")
saveRDS(thing$D2, file="./FP/C_D2_FP.rds")
#Punta Gorda
thing <- decompose(C_PG)
saveRDS(thing$D1, file="./PG/C_D1_PG.rds")
saveRDS(thing$G, file="./PG/C_G_PG.rds")
saveRDS(thing$F, file="./PG/C_F_PG.rds")
saveRDS(thing$D2, file="./PG/C_D2_PG.rds")
#Wild Turkey
thing <- decompose(C_WT)
saveRDS(thing$D1, file="./WT/C_D1_WT.rds")
saveRDS(thing$G, file="./WT/C_G_WT.rds")
saveRDS(thing$F, file="./WT/C_F_WT.rds")
saveRDS(thing$D2, file="./WT/C_D2_WT.rds")


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

toKvals_D2 = function(elas_D2) {
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
  return(Kvals_elas_D2)
}

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

# Convert to totals: #####

total.contrib_D1_BC<-apply(Kvals_contrib_D1_BC, c(3,4), sum);
total.contrib_D2_BC<-apply(Kvals_contrib_D2_BC, c(3,4), sum);
total.contrib_F_BC<-apply(Kvals_contrib_F_BC, c(3,4), sum);
total.contrib_G_BC<-apply(Kvals_contrib_G_BC, c(3,4), sum);
save(total.contrib_D1_BC, total.contrib_D2_BC, total.contrib_F_BC, total.contrib_G_BC, 
     file="./BC/total.contrib_BC.RData")




#Collapse contributions by site: #####	
load("./Western/total.contrib_W.RData")
 W<-c( 
 sum(rowSums(total.contrib_D1_W)), 
 sum(rowSums(total.contrib_G_W)), 
 sum(rowSums(total.contrib_F_W)), 
 sum(rowSums(total.contrib_D2_W)) )
 save(W, file="./Western/W.RData")
 
load("./Hybrid/total.contrib_H.RData")
 H<-c(
 sum(rowSums(total.contrib_D1_H)), 
 sum(rowSums(total.contrib_G_H)), 
 sum(rowSums(total.contrib_F_H)), 
 sum(rowSums(total.contrib_D2_H)) )
 save(H, file="./Hybrid/H.RData")
 
 load("./Eastern/total.contrib_E.RData")
 E<-c( 
 sum(rowSums(total.contrib_D1_E)), 
 sum(rowSums(total.contrib_G_E)), 
 sum(rowSums(total.contrib_F_E)), 
 sum(rowSums(total.contrib_D2_E)) )
 save(E, file="./Eastern/E.RData")


# Calculate matrices of coefficients of variation #####

	
load("./Eastern/D1_E.RData")
load("./Hybrid/D1_H.RData")
load("./Western/D1_W.RData")
	
CV<-function(x) {
	return(sd(x)/mean(x))
	
}



	
cv_D1<-matrix(nrow=(m1*m2), ncol=(m1*m2))

for(i in 1:(m1*m2)) {
	for(j in 1:(m1*m2)) {
		cv_D1[i,j]<-CV(c(D1_E[i,j], D1_W[i,j], D1_H[i,j]))
		if(is.na(cv_D1[i,j])) {cv_D1[i,j]<-0}
	}
}


D1_avg=(D1_E+D1_H+D1_W)/3
rm(D1_E, D1_H, D1_W)

#G
load("./Eastern/G_E.RData")
load("./Hybrid/G_H.RData")
load("./Western/G_W.RData")

cv_G<-matrix(nrow=(m3*m4), ncol=(m1*m2))


for(i in 1:(m3*m4)) {
	for(j in 1:(m1*m2)) {
		cv_G[i,j]<-CV(c(G_E[i,j], G_W[i,j], G_H[i,j]))
		if(is.na(cv_G[i,j])) {cv_G[i,j]<-0}
	}
}

G_avg<-(G_E+G_H+G_W)/3

rm(G_E)
rm(G_H)
rm(G_W)

#F
load("./Eastern/F_E.RData")
load("./Hybrid/F_H.RData")
load("./Western/F_W.RData")
cv_F<-matrix(nrow=(m1*m2), ncol=(m3*m4))
for(i in 1:(m1*m2)) {
	for(j in 1:(m3*m4)) {
		cv_F[i,j]<-CV(c(F_E[i,j], F_W[i,j], F_H[i,j]))
		if(is.na(cv_F[i,j])) {cv_F[i,j]<-0}
	}
}

F_avg<-(F_E + F_W + F_H)/3

rm(F_E)
rm(F_W)
rm(F_H)

#D2
load("./Eastern/D2_E.RData")
load("./Hybrid/D2_H.RData")
load("./Western/D2_W.RData")


cv_D2<-matrix(nrow=(m3*m4), ncol=(m3*m4))
for(i in 1:(m3*m4)) {
	for(j in 1:(m3*m4)) {
		cv_D2[i,j]<-CV(c(D2_E[i,j], D2_W[i,j], D2_H[i,j]))
		if(is.na(cv_D2[i,j])) {cv_D2[i,j]<-0}
	}
}


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








