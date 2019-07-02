# Sensitivity analysis: how sensitive is lambda to changes in underlying parameters? 

#Load the parameter vectors calculated from the dataset where the one problematic individual 
# at Big Cypress whose future diameter was recorded incorrectly. 


load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/BC/p.vec_BC.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/CC/p.vec_CC.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/C/p.vec_C.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/FP/p.vec_FP.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/PG/p.vec_PG.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/WT/p.vec_WT.RData")


sim_grad_diams_mu <- c(1.7, 2, 2.3, 2.6, 2.9)

sim_grad_diams_sd <- c(0.25, 0.375, 0.5)
lambda <- rep(NA, 6)
p.vecs <- list(p.vec_BC, p.vec_CC, p.vec_C, p.vec_FP, p.vec_PG, p.vec_WT)
getLambdas <- function(p.vecs) {
  
  # Part II: Building the IPM #####
  m1=10
  m2=m1+1
  m3=100
  m4=m3+1
  tol=1.e-8; 
  
  # Define the kernels and iteration matrix:
  #Domain 1: Seedling Domain
  #	for diam=diameter in range [0, 1.6]
  #	for height in range [0, 16]
  #Domain 2: Larger Domain
  #	 for diam=diameter in range [1.6,700]
  #and for height=height in range [16, 800]
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
  h3=(700-1.6)/m3;
  y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
  h4=(800-16)/m4
  y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2
  
  # Compute the iteration matrix. With a bit of vectorizing it's not too slow,
  # though you can probably do better if you need to. The shortcuts here have 
  # been checked against the results from code that uses loops for everything. (comment from Ellner and Rees)

  
  for (j in 1:6) {
    
    p.vec <- p.vecs[[j]]
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
        
      }
      
      D1=D1*h1*h2 #multiply D1 by widths
      return(list(D1 = D1, Kvals_D1 = Kvals_D1))
    }
    
    thing <- build_D1(p.vec)
    D1 <- thing$D1
    
    # Construct D2 (Large Domain): #####
    build_D2 = function(p.vec) {
      plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
      Plop=outer(1:m3,1:m4,plop); 
      
      D2=matrix(0,m3*m4,m3*m4);
      Kvals_D2=array(0,c(m3,m4,m3,m4));  
      
      for(i in 1:m3){
        for(j in 1:m4){
          for(k in 1:m3){
            kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec)
            D2[Plop[k,1:m4],Plop[i,j]]=kvals
            Kvals_D2[k,1:m4,i,j]=kvals
            
          }}
        
      }		
      D2=D2*h3*h4 #Multiply D2 by widths
      return(list(D2 = D2, Kvals_D2 = Kvals_D2))
    }
    
    thing <- build_D2(p.vec)
    D2 <-thing$D2
    
    # Construct F (Fertility): #####
    build_F = function(p.vec) {
      plop1=function(i, j) {(j-1)*m1 + i}
      plop2=function(i, j) {(j-1)*m3 + i}
      Plop1=outer(1:m1,1:m2,plop1); 
      Plop2=outer(1:m3, 1:m4, plop2);
      
      F=matrix(0,m1*m2,m3*m4); 
      Kvals_F=array(0, c(m1, m2, m3, m4))
      
      for(i in 1:m3) {
        for (j in 1:m4) {
          for (k in 1:m1) {
            kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec)
            F[Plop1[k, 1:m2], Plop2[i,j]]=kvals
            Kvals_F[k, 1:m2, i, j]=kvals
          }}
      }
      F=F*h1*h2
      return(list(F = F, Kvals_F = Kvals_F))
    }
    
    
    thing <- build_F(p.vec)
    F<- thing$F
    
    # Construct M (Maturation): #####
    build_G = function(p.vec) {
      plop1=function(i, j) {(j-1)*m3 + i}
      plop2=function(i, j) {(j-1)*m1 + i}
      Plop1=outer(1:m3,1:m4,plop1); 
      Plop2=outer(1:m1, 1:m2, plop2);
      
      G=matrix(0,m3*m4,m1*m2); 
      Kvals_G=array(0, c(m3, m4, m1, m2))
      
      for(i in 1:m1) {
        for (j in 1:m2) {
          for (k in 1:m3) {
            kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec)
            G[Plop1[k, 1:m4], Plop2[i,j]]=kvals
            Kvals_G[k, 1:m4, i, j]=kvals
          }}
      }
      G=G*h3*h4
      return(list(G = G, Kvals_G = Kvals_G))
    }
    
    thing <- build_G(p.vec)
    G <- thing$G
    
    rm(thing)
    gc()
    # Assemble the matrix #####
    
    A <- cbind(rbind(D1, G), rbind(F, D2))
    
    
    #  Find lambda, w by iteration #####
    
    #  Note: the Matrix package is used to iterate more quickly. The check  
    #  for convergence requires extracting matrix entries via the @x slot
    #  of a Matrix object. Matrix is S4-style -- see ?Matrix. 
    
    find_lambda = function(A) {
      
      A2=Matrix(A); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 
      
      qmax=1000; lam=1; 
      while(qmax>tol) {
        nt1=A2%*%nt;
        qmax=sum(abs((nt1-lam*nt)@x));  
        lam=sum(nt1@x); 
        nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
        
      } 
      nt=matrix(nt@x,m1*m2+m3*m4,1); 
      #stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
      stable.dist=nt
      lam.stable=lam;
      
      # Check that the @bits worked as intended.   
      qmax=sum(abs(lam*nt-A%*%nt)); 
      
      
      #Find the reproductive value function by iteration
      vt=Matrix(1,1,m1*m2+m3*m4); vt1=vt; 
      
      qmax=1000; lam=1; 
      while(qmax>tol) {
        vt1=vt%*%A2;
        qmax=sum(abs((vt1-lam*vt)@x));  
        lam=sum(vt1@x); 
        vt@x=(vt1@x)/lam;   
      } 
      v=t(matrix(vt@x,1,m1*m2+m3*m4)); 
      lam.stable.t=lam; 
      
      return(list(lam.stable = lam.stable, stable.dist = stable.dist, v=v))
    }
    
    thing <- find_lambda(A)
    lambda[j] <- thing$lam.stable
    cat("Site: ", j,"\n");
  }
  return(lambda)
}

#First calculate the "real" values of lambda
thing <- getLambdas(p.vecs)
results_mu <- data.frame( mu_diam = rep(NA, 6),
                          lambda_BC = rep(NA, 6),
                          lambda_CC = rep(NA, 6),
                          lambda_C  = rep(NA, 6),
                          lambda_FP = rep(NA, 6),
                          lambda_PG = rep(NA, 6),
                          lambda_WT = rep(NA, 6))
results_mu[1,] <- c(p.vec_BC[15], thing)

#Now go through and see what happens as the parameter value 
# for mean graduate diameter is changed


for(i in 1:length(sim_grad_diams_mu)) {
  p.vec_BC_new <- p.vec_BC
  p.vec_CC_new <- p.vec_CC
  p.vec_C_new  <- p.vec_C
  p.vec_FP_new <- p.vec_FP
  p.vec_PG_new <- p.vec_PG
  p.vec_WT_new <- p.vec_WT
  p.vecs <- list(p.vec_BC_new, p.vec_CC_new, p.vec_C_new, 
                 p.vec_FP_new, p.vec_PG_new, p.vec_WT_new)
  
  for (k in 1:6){
    p.vecs[[k]][15] <- sim_grad_diams_mu[i]
  }
  thing <- getLambdas(p.vecs)
  results_mu[i+1,] <- c(sim_grad_diams_mu[i], thing)
  cat("Param: ", i,"\n"); 
}



#And finally for the standard deviation of graduate diameter: 
#Now go through and see what happens as the parameter value 
# for mean graduate diameter is changed

results_sd <- data.frame( sd_diam = rep(NA, 4),
                          lambda_BC = rep(NA, 4),
                          lambda_CC = rep(NA, 4),
                          lambda_C  = rep(NA, 4),
                          lambda_FP = rep(NA, 4),
                          lambda_PG = rep(NA, 4),
                          lambda_WT = rep(NA, 4))

results_sd[1,1] <- p.vec_BC[16]
results_sd[, 2:7] <- results_mu[, 2:7]

for(i in 1:length(sim_grad_diams_sd)) {
  p.vec_BC_new <- p.vec_BC
  p.vec_CC_new <- p.vec_CC
  p.vec_C_new  <- p.vec_C
  p.vec_FP_new <- p.vec_FP
  p.vec_PG_new <- p.vec_PG
  p.vec_WT_new <- p.vec_WT
  p.vecs <- list(p.vec_BC_new, p.vec_CC_new, p.vec_C_new, 
                 p.vec_FP_new, p.vec_PG_new, p.vec_WT_new)
  
  for (k in 1:6){
    p.vecs[[k]][16] <- sim_grad_diams_sd[i]
  }
  thing <- getLambdas(p.vecs)
  results_sd[i+1,] <- c(sim_grad_diams_sd[i], thing)
  cat("Param: ", i,"\n"); 
}

#Plot the relationship between parameter choice for mean graduate diameter and lambda
par(mfrow=c(2,3))
plot(sim_grad_diams_mu, results_mu[2:6,2], main = "BC", 
     xlab="Mean diameter of graduates", ylab="Lambda")
points(results_mu[1,1], results_mu[1,2], col="red", pch=19)

plot(sim_grad_diams_mu, results_mu[2:6,3], main = "CC",
     xlab="Mean diameter of graduates", ylab="Lambda")
points(results_mu[1,1], results_mu[1,3], col="red", pch=19)

plot(sim_grad_diams_mu, results_mu[2:6,4], main = "C",
     xlab="Mean diameter of graduates", ylab="Lambda")
points(results_mu[1,1], results_mu[1,4], col="red", pch=19)

plot(sim_grad_diams_mu, results_mu[2:6,5], main = "FP",
     xlab="Mean diameter of graduates", ylab="Lambda")
points(results_mu[1,1], results_mu[1,5], col="red", pch=19)

plot(sim_grad_diams_mu, results_mu[2:6,6], main = "PG",
     xlab="Mean diameter of graduates", ylab="Lambda")
points(results_mu[1,1], results_mu[1,6], col="red", pch=19)

plot(sim_grad_diams_mu, results_mu[2:6,7], main = "WT",
     xlab="Mean diameter of graduates", ylab="Lambda")
points(results_mu[1,1], results_mu[1,7], col="red", pch=19)



#Plot the relationship between parameter choice for standard deviation of graduate diameter and lambda
par(mfrow=c(2,3))
plot(sim_grad_diams_sd, results_sd[2:4,2], main = "BC", 
     xlab="sd of graduate diameter", ylab="Lambda")
points(results_sd[1,1], results_sd[1,2], col="red", pch=19)

plot(sim_grad_diams_sd, results_sd[2:4,3], main = "CC",
     xlab="sd of graduate diameter", ylab="Lambda")
points(results_sd[1,1], results_sd[1,3], col="red", pch=19)

plot(sim_grad_diams_sd, results_sd[2:4,4], main = "C",
     xlab="sd of graduate diameter", ylab="Lambda")
points(results_sd[1,1], results_sd[1,4], col="red", pch=19)

plot(sim_grad_diams_sd, results_sd[2:4,5], main = "FP",
     xlab="sd of graduate diameter", ylab="Lambda")
points(results_sd[1,1], results_sd[1,5], col="red", pch=19)

plot(sim_grad_diams_sd, results_sd[2:4,6], main = "PG",
     xlab="sd of graduate diameter", ylab="Lambda")
points(results_sd[1,1], results_sd[1,6], col="red", pch=19)

plot(sim_grad_diams_sd, results_sd[2:4,7], main = "WT",
     xlab="sd of graduate diameter", ylab="Lambda")
points(results_sd[1,1], results_sd[1,7], col="red", pch=19)


sim_grad_diams_mu2 <- c(3, 6, 9, 12, 15,18)

for(i in 1:length(sim_grad_diams_mu)) {
  p.vec_BC_new <- p.vec_BC
  p.vec_CC_new <- p.vec_CC
  p.vec_C_new  <- p.vec_C
  p.vec_FP_new <- p.vec_FP
  p.vec_PG_new <- p.vec_PG
  p.vec_WT_new <- p.vec_WT
  p.vecs <- list(p.vec_BC_new, p.vec_CC_new, p.vec_C_new, 
                 p.vec_FP_new, p.vec_PG_new, p.vec_WT_new)
  
  for (k in 1:6){
    p.vecs[[k]][15] <- sim_grad_diams_mu[i]
  }
  thing <- getLambdas(p.vecs)
  results_mu[i+1,] <- c(sim_grad_diams_mu2[i], thing)
  cat("Param: ", i,"\n"); 
}

