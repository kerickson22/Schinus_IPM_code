library(msm)


results <- data.frame("m1" = rep(10, 6),
                      "m2" = rep(11, 6),
                       "m3a" = rep(NA, 6),
                       "m4a" = rep(50, 6),
                       "m3b" = rep(50, 6),
                       "m4b" = rep(50, 6), 
                      "lambda" = rep(NA, 6))
results$m3a <- c(50, 100, 150, 200, 250, 300)
results2 <- data.frame("m1" = rep(10, 6),
                      "m2" = rep(11, 6),
                      "m3a" = rep(50, 6),
                      "m4a" = rep(NA, 6),
                      "m3b" = rep(50, 6),
                      "m4b" = rep(50, 6), 
                      "lambda" = rep(NA, 6))

results2$m4a <- c(50, 100, 150, 200, 250, 300)

results3<- data.frame("m1" = rep(10, 6),
                       "m2" = rep(11, 6),
                       "m3a" = rep(50, 6),
                       "m4a" = rep(50, 6),
                       "m3b" = rep(NA, 6),
                       "m4b" = rep(50, 6), 
                       "lambda" = rep(NA, 6))

results3$m3b <- c(50, 100, 150, 200, 250, 300)

results4<- data.frame("m1" = rep(10, 6),
                      "m2" = rep(11, 6),
                      "m3a" = rep(50, 6),
                      "m4a" = rep(50, 6),
                      "m3b" = rep(50, 6),
                      "m4b" = rep(NA, 6), 
                      "lambda" = rep(NA, 6))

results4$m4b <- c(50, 100, 150, 200, 250, 300)

temp <- c(10, 11, 75, 75, 75, 75, NA)


results<-rbind(results, results2, results3, results4)
results <- rbind(results, temp)


temp <- c(10, 11, 50, 50, 200, 100, NA)

results <- rbind(results, temp)

results5<- data.frame("m1" = rep(NA, 12),
                      "m2" = rep(NA, 12),
                      "m3a" = rep(50, 12),
                      "m4a" = rep(50, 12),
                      "m3b" = rep(200, 12),
                      "m4b" = rep(100, 12), 
                      "lambda" = rep(NA, 12))

results5$m1 <- c(5, 15, 20, 25, 30, 35, 10, 10, 10, 10, 10, 10)
results5$m2 <- c(11, 11, 11, 11, 11, 11, 5, 15, 20, 25, 30, 35)

results <- rbind(results, results5)
load("./Overall/p.vec_overall.RData")
p.vec <- p.vec_overall


for(ii in 27:nrow(results)) {

# Set matrix size (to show up errors) and convergence tolerance. 
# Make m1 and m2 smaller to make this run faster.  
m1=results$m1[ii]
m2=results$m2[ii]

m3a <- results$m3a[ii]
m3b <- results$m3b[ii]
m3 <- m3a + m3b
m4a <- results$m4a[ii]
m4b <- results$m4b[ii]
m4=m4a + m4b 
tol=1.e-8; 

# Compute meshpoints  

###Fix these to match domains 
h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1

h3a = (150-1.6)/m3a; 
#h3=(700-1.6)/m3;
#h3=3
y3a=(h3a/2)*((0:(m3a-1))+(1:m3a))+1.6; #for diameter in D2a
h3b = (700-150)/m3b
y3b=(h3b/2)*((0:(m3b-1))+(1:m3b))+150; #for diameter in D2b
y3 <- c(y3a, y3b)

#h4=(800-16)/m4
h4a = (300-16)/m4a
y4a=(h4a/2)*((0:(m4a-1))+(1:m4a))+150;

h4b = (800-300)/m4b 
y4b=(h4b/2)*((0:(m4b-1))+(1:m4b))+300;

y4 = c(y4a, y4b)







# Define the kernels and iteration matrix:
#Domain 1: Seedling Domain
#	for diam=diameter in range [0, 1.6]
#	for height in range [0, 16]
#Domain 2: Larger Domain
#2a: diam= diameter in range [1.6, 150]
#     height in range [16, 300]
#2b for diam=diameter in range [150,700]
#and for height=height in range [300, 800]
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



# Compute the iteration matrix. With a bit of vectorizing it's not too slow,
# though you can probably do better if you need to. The shortcuts here have 
# been checked against the results from code that uses loops for everything. (comment from Ellner and Rees)


# Construct D1 (Seedling Domain): #####

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
  
# Construct D2AA (Large Domain): #####

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

# Construct D2BB #####
  
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
  
# Construct D2BA #####

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
# Construct D2AB 
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

D2 <- rbind(cbind(D2AA, D2BA), cbind(D2AB, D2BB))   
rm(D2AA, D2BA, D2AB, D2BB)
# Construct FA (Fertility): #####

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
#Construct FB #####
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
F <- cbind(FA, FB)
rm (FA, FB)

# Construct MA (Maturation): #####

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
  
# Construct MB #####
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
G <- rbind(GA, GB)
rm(GA, GB)

# Assemble the matrices #####


left_side<-rbind(D1, G)
rm(D1, G)
right_side<-rbind(F, D2)
rm(F, D2)
A<-cbind(left_side, right_side)
rm(left_side, right_side)
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

results$lambda[ii] = lam.stable
rm(A)
rm(A2)
}

write.csv(results, file="meshpoint_results.csv")
