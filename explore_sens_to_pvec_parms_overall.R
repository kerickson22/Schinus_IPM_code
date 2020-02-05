# Sensitivity analysis: how sensitive is lambda to changes in underlying parameters? 

#Load the parameter vectors calculated from the dataset where the one problematic individual 
# at Big Cypress whose future diameter was recorded incorrectly is removed . 
require(Matrix);
library(MASS)
library(mvtnorm)
library(msm)


load("./BC/p.vec_BC.RData")
load("./CC/p.vec_CC.RData")
load("./C/p.vec_C.RData")
load("./FP/p.vec_FP.RData")
load("./PG/p.vec_PG.RData")
load("./WT/p.vec_WT.RData")
#results <- read.csv("parms_to_explore.csv", head=T)
#results2 <- read.csv("parms_to_explore2.csv", head=T)
#results <- rbind(results, results2)
#names(results) <- c("pos", "value")

percentages <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5)

newparm <- c(0.2, 0.02, 0.001, 0.003, 0.0002)
results <- data.frame(pos  = rep(37, length(newparm)),
                      newparm = newparm,
                      lambda_BC = rep(NA, length(newparm)),
                      lambda_CC = rep(NA, length(newparm)),
                      lambda_C = rep(NA, length(newparm)),
                      lambda_FP = rep(NA, length(newparm)),
                      lambda_PG = rep(NA, length(newparm)),
                      lambda_WT = rep(NA, length(newparm)))





getLambda <- function(p.vec) {
  
  # Part II: Building the IPM #####
  m1=10
  m2=m1+1
  m3a=250
  m4a=50
  
  m3b = 120
  m4b = 100
  m3=m3a + m3b
  m4=m4a + m4b
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
        
      }
      
      D1=D1*h1*h2 #multiply D1 by widths
      return(list(D1 = D1, Kvals_D1 = Kvals_D1))
    }
    
    thing <- build_D1(p.vec)
    D1 <- thing$D1
    
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
    
    thing <- build_D2AA(p.vec)
    D2AA <- thing$D2AA
    
    # D2BB #####
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
    
    thing <- build_D2BB(p.vec)
    D2BB <- thing$D2BB
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
    
    
    thing <- build_D2BA(p.vec)
    D2BA <- thing$D2BA
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
    
    
    thing <- build_D2AB(p.vec)
    D2AB <- thing$D2AB
    D2 <- rbind(cbind(D2AA, D2BA), cbind(D2AB, D2BB)) 
    rm(D2AA, D2BA, D2AB, D2BB)
    # Construct F #####
    
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
    
    thing <- build_FA(p.vec)
    FA <- thing$FA
    thing <- build_FB(p.vec)
    FB <-thing$FB
    
    F<- cbind(FA, FB)
    rm(FA, FB)
    # Construct M #####
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
    
    
    thing <- build_GA(p.vec)
    GA<-thing$GA
    
    thing <- build_GB(p.vec)
    
    GB<-thing$GB
    
    
    G <- rbind(GA, GB)
    rm(GA, GB)
    rm(thing)
    gc()
    # Assemble the matrix #####
    
    A <- cbind(rbind(D1, G), rbind(F, D2))
    rm(D1, G, F, D2)
    
    
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
      
      lam.stable=lam;
      
      return(list(lam.stable = lam.stable))
    }
    thing <- find_lambda(A)
    return(list(lam.stable = thing$lam.stable))
}    
    



time_start <- proc.time()
for(jj in 2:6) {
for(ii in 1:5) {
  cat("Param set: ", ii, "of ", nrow(results),"Site:", jj, ": \n");
  p.vec_new_BC <- p.vec_BC
 p.vec_new_CC <- p.vec_CC
   p.vec_new_C  <- p.vec_C 
   p.vec_new_FP <- p.vec_FP
   p.vec_new_PG <- p.vec_PG
   p.vec_new_WT <- p.vec_WT

  
  if(jj==1){p.vec_new_BC[results$pos[ii]] <- results$newparm[ii]
  thing <- getLambda(p.vec_new_BC)
  results[ii, jj+2] <- thing$lam.stable}
   
  if(jj==2){p.vec_new_CC[results$pos[ii]] <- results$newparm[ii]
  thing <- getLambda(p.vec_new_CC)
  results[ii, jj+2] <- thing$lam.stable}
   
  if(jj==3){ p.vec_new_C[results$pos[ii]] <- results$newparm[ii]
  thing <- getLambda(p.vec_new_C)
  results[ii, jj+2] <- thing$lam.stable}
   
  if(jj==4){ p.vec_new_FP[results$pos[ii]] <- results$newparm[ii]
  thing <- getLambda(p.vec_new_FP)
  results[ii, jj+2] <- thing$lam.stable}
   
  if(jj==5){p.vec_new_PG[results$pos[ii]] <- results$newparm[ii]
  thing <- getLambda(p.vec_new_PG)
  results[ii, jj+2] <- thing$lam.stable}
   
  if(jj==6){p.vec_new_WT[results$pos[ii]] <- results$newparm[ii]
  thing <- getLambda(p.vec_new_WT)
  results[ii, jj+2] <- thing$lam.stable}
}
}
time_end <- proc.time() - time_start 

results$delta_BC <- results$lambda_BC - lam.stable_BC
results$delta_CC <- results$lambda_CC - lam.stable_CC
results$delta_C <- results$lambda_C - lam.stable_C
results$delta_FP <- results$lambda_FP - lam.stable_FP
results$delta_PG <- results$lambda_PG - lam.stable_PG
results$delta_WT <- results$lambda_WT - lam.stable_WT
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/tau.eps", width=width.cm/2.54, 
           height=width.cm/(2.54), pointsize=pointsize,  encoding = "TeXtext.enc")
plot(log10(results$newparm), results$lambda_BC, ylim=c(0.98, 1.75),
     col=cols[1], pch=19, xlab = expression(paste("log(", tau[1], ")")),
     ylab = expression(paste(lambda)))
points(log10(0.002), lam.stable_BC, col=cols[1], pch=4)
points(log10(results$newparm), results$lambda_CC, col=cols[2], pch=19)
points(log10(0.002), lam.stable_CC, col=cols[2], pch=4)
points(log10(results$newparm), results$lambda_C, col=cols[3],  pch=19)
points(log10(0.002), lam.stable_C, col=cols[3], pch=4)

points(log10(results$newparm), results$lambda_FP, col=cols[4],  pch=19)
points(log10(0.002), lam.stable_FP, col=cols[4], pch=4)
points(log10(results$newparm), results$lambda_PG, col=cols[5],  pch=19)
points(log10(0.002), lam.stable_PG, col=cols[5], pch=4)
points(log10(results$newparm), results$lambda_WT, col=cols[6],  pch=19)
points(log10(0.002), lam.stable_WT, col=cols[6], pch=4)
legend("topleft", col=cols, pch=19, 
       c("BC", "CC", "C", "FP", "PG", "WT"))
dev.off()




lam.stable_overall <- readRDS("./Overall/lam.stable_overall.rds")
lam.stable_BC <- readRDS("./BC/lam.stable_BC.rds")
lam.stable_CC <- readRDS("./CC/lam.stable_CC.rds")
lam.stable_C  <- readRDS("./C/lam.stable_C.rds")
lam.stable_FP <- readRDS("./C/lam.stable_C.rds")
lam.stable_PG <- readRDS("./FP/lam.stable_FP.rds")
lam.stable_WT <- readRDS("./WT/lam.stable_WT.rds")

load("./BC/p.vec_BC.RData")
load("./CC/p.vec_CC.RData")
load("./C/p.vec_C.RData")
load("./FP/p.vec_FP.RData")
load("./PG/p.vec_PG.RData")
load("./WT/p.vec_WT.RData")

results$delta <- results$lambda - lam.stable_overall

delt_BC <- lam.stable_BC - lam.stable_overall
delt_CC <- lam.stable_CC - lam.stable_overall
delt_C <- lam.stable_C - lam.stable_overall
delt_FP <- lam.stable_FP - lam.stable_overall
delt_PG <- lam.stable_PG - lam.stable_overall
delt_WT <- lam.stable_WT - lam.stable_overall

# -Seedling Survival  #####
#logit(Ss) = B0 + B1*diam + B2*height
par(mfrow=c(2,2))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=F)
text(3,8, "Seedling Survival", cex=1)
expression('title'[2])
text(4, 7, expression("Logit (S"[s]*") = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
     cex=0.8)

plot(results$new_par[1:10], results$delta[1:10],
     xlim=c(-4.1, -1),
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[1], col="red")
abline(h =0)
abline(v = p.vec_BC[1], lty=3, col=cols[1])
text(p.vec_BC[1], 1e-06, "BC")
abline(v = p.vec_CC[1], lty=3, col=cols[2])
text(p.vec_CC[1], 1e-06, "CC")
abline(v = p.vec_C[1],  lty=3, col=cols[3])
text(p.vec_C[1], 1e-06, "C")
abline(v = p.vec_FP[1], lty=3, col=cols[4])
text(p.vec_FP[1], 1e-06, "FP")
abline(v = p.vec_PG[1], lty=3, col=cols[5])
text(p.vec_PG[1], 8.5e-07, "PG")
abline(v = p.vec_WT[1], lty=3, col=cols[6])
text(p.vec_WT[1], 7e-07, "WT")

plot(results$new_par[11:20], results$delta[11:20],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[11], col="red")
abline(h =0)
abline(v = p.vec_BC[2], lty=3, col=cols[1])
text(p.vec_BC[2], 3e-7, "BC")
text(p.vec_CC[2], 2.5e-7, "CC")
text(p.vec_C[2], 2e-7, "C")
text(p.vec_FP[2], 1.5e-7, "FP")
text(p.vec_PG[2], 1e-7, "PG")
text(p.vec_WT[2], 6.5e-8, "WT")

plot(results$new_par[21:30], results$delta[21:30],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[21], col="red")
abline(h =0)

abline(v = p.vec_BC[3], lty=3, col=cols[1])
text(p.vec_BC[3], 4e-7, "BC")
text(p.vec_CC[3], 3.5e-7, "CC")
text(p.vec_C[3], 3e-7, "C")
text(p.vec_FP[3], 2.5e-7, "FP")
text(p.vec_PG[3], 2e-7, "PG")
text(p.vec_WT[3], 1.5e-7, "WT")
dev.copy(pdf, 'seedling_survival_sens.pdf')
dev.off()

#-- Seedling Growth #####

par(mfrow=c(4, 2))


#Seedling Survival logit(Ss)  = B0 + B1*diam + B2*height
x11(width=11)
par(mfrow=c(2,4))

#Diameter
plot(results$new_par[31:40], results$delta[31:40],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[31], col="red")
abline(h =0)
abline(v = p.vec_BC[4], lty=3, col=cols[1])
text(p.vec_BC[4], 1.5e-07, "BC")
abline(v = p.vec_CC[4], lty=3, col=cols[2])
text(p.vec_CC[4], 1.5e-07, "CC")
abline(v = p.vec_C[4],  lty=3, col=cols[3])
text(p.vec_C[4], 1e-07, "C")
abline(v = p.vec_FP[4], lty=3, col=cols[4])
text(p.vec_FP[4], 1e-07, "FP")
abline(v = p.vec_PG[4], lty=3, col=cols[5])
text(p.vec_PG[4], 1.25e-07, "PG")
abline(v = p.vec_WT[4], lty=3, col=cols[6])
text(p.vec_WT[4], 0.5e-07, "WT")

plot(results$new_par[41:50], results$delta[41:50],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
mtext(side=3, expression(mu[diam]*" = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
      cex=0.8)
abline(v = results$par_original[41], col="red")
abline(h =0)
abline(v = p.vec_BC[5], lty=3, col=cols[1])
text(p.vec_BC[5],1.5e-7, "BC")
text(p.vec_CC[5], 1.25e-7, "CC")
text(p.vec_C[5], 10e-8, "C")
text(p.vec_FP[5], 8e-8, "FP")
text(p.vec_PG[5], 6.5e-8, "PG")
text(p.vec_WT[5], 5e-8, "WT")

plot(results$new_par[51:60], results$delta[51:60],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[51], col="red")
abline(h =0)

abline(v = p.vec_BC[6], lty=3, col=cols[1])
text(p.vec_BC[6], 4e-8, "BC")
text(p.vec_CC[6], 3.5e-8, "CC")
text(p.vec_C[6], 3e-8, "C")
text(p.vec_FP[6], 2.5e-8, "FP")
text(p.vec_PG[6], 2e-8, "PG")
text(p.vec_WT[6], 1.5e-8, "WT")

plot(results$new_par[61:70], results$delta[61:70],
     xlab = expression(sigma[diam]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[61], col="red")
abline(h =0)

abline(v = p.vec_BC[7], lty=3, col=cols[1])
text(p.vec_BC[7], -1.5e-9, "BC")
text(p.vec_CC[7], -2e-9, "CC")
text(p.vec_C[7], -2.5e-9, "C")
text(p.vec_FP[7], -3e-9, "FP")
text(p.vec_PG[7], -3.5e-9, "PG")
text(p.vec_WT[7], -4e-9, "WT")

#Height
plot(results$new_par[71:80], results$delta[71:80],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[71], col="red")
abline(h =0)
abline(v = p.vec_BC[8], lty=3, col=cols[1])
text(p.vec_BC[8], 2e-07, "BC")
abline(v = p.vec_CC[8], lty=3, col=cols[2])
text(p.vec_CC[8], 2e-07, "CC")
abline(v = p.vec_C[8],  lty=3, col=cols[3])
text(p.vec_C[8], 1.5e-07, "C")
abline(v = p.vec_FP[8], lty=3, col=cols[4])
text(p.vec_FP[8], 1.25e-07, "FP")
abline(v = p.vec_PG[8], lty=3, col=cols[5])
text(p.vec_PG[8], 1e-07, "PG")
abline(v = p.vec_WT[8], lty=3, col=cols[6])
text(p.vec_WT[8], 2e-07, "WT")

plot(results$new_par[81:90], results$delta[81:90],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
mtext(side=3, expression(mu[height]*" = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
      cex=0.8)
abline(v = results$par_original[81], col="red")
abline(h =0)
abline(v = p.vec_BC[9], lty=3, col=cols[1])
text(p.vec_BC[9],4e-8, "BC")
text(p.vec_CC[9], 3.5e-8, "CC")
text(p.vec_C[9], 3e-8, "C")
text(p.vec_FP[9], 2.5e-8, "FP")
text(p.vec_PG[9], 2e-8, "PG")
text(p.vec_WT[9], 1.5e-8, "WT")

plot(results$new_par[91:100], results$delta[91:100],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[91], col="red")
abline(h =0)

abline(v = p.vec_BC[10], lty=3, col=cols[1])
text(p.vec_BC[10], 9e-8, "BC")
text(p.vec_CC[10], 7.5e-8, "CC")
text(p.vec_C[10], 6e-8, "C")
text(p.vec_FP[10], 4.5e-8, "FP")
text(p.vec_PG[10], 3e-8, "PG")
text(p.vec_WT[10], 1.5e-8, "WT")

plot(results$new_par[101:110], results$delta[101:110],
     xlab = expression(sigma[height]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[101], col="red")
abline(h =0)

abline(v = p.vec_BC[11], lty=3, col=cols[1])
text(p.vec_BC[11], -1e-9, "BC")
text(p.vec_CC[11], -2e-9, "CC")
text(p.vec_C[11], -3e-9, "C")
text(p.vec_FP[11], -4e-9, "FP")
text(p.vec_PG[11], -5e-9, "PG")
text(p.vec_WT[11], -6e-9, "WT")

dev.copy(pdf, 'seedling_growth_sens.pdf')
dev.off()

## Maturation ##### 
x11(width=11)
par(mfrow=c(2,4))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=F)
text(3,8, "Maturation", cex=1)
expression('title'[2])
text(4, 7, expression("Logit (m) = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
     cex=0.8)

plot(results$new_par[111:120], results$delta[111:120],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[111], col="red")
abline(h =0)
abline(v = p.vec_BC[12], lty=3, col=cols[1])
text(p.vec_BC[12], 1e-06, "BC")
abline(v = p.vec_CC[12], lty=3, col=cols[2])
text(p.vec_CC[12], 9e-07, "CC")
abline(v = p.vec_C[12],  lty=3, col=cols[3])
text(p.vec_C[12], 1e-06, "C")
abline(v = p.vec_FP[12], lty=3, col=cols[4])
text(p.vec_FP[12], 7.5e-07, "FP")
abline(v = p.vec_PG[12], lty=3, col=cols[5])
text(p.vec_PG[12], 6e-07, "PG")
abline(v = p.vec_WT[12], lty=3, col=cols[6])
text(p.vec_WT[12], 1e-06, "WT")

plot(results$new_par[121:130], results$delta[121:130],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[121], col="red")
abline(h =0)
abline(v = p.vec_BC[13], lty=3, col=cols[1])
text(p.vec_BC[13], 4e-7, "BC")
text(p.vec_CC[13], 3.5e-7, "CC")
text(p.vec_C[13], 3e-7, "C")
text(p.vec_FP[13], 2.5e-7, "FP")
text(p.vec_PG[13], 2e-7, "PG")
text(p.vec_WT[13], 1.5e-7, "WT")

plot(results$new_par[131:140], results$delta[131:140],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[131], col="red")
abline(h =0)

abline(v = p.vec_BC[14], lty=3, col=cols[1])
text(p.vec_BC[14], 3e-7, "BC")
text(p.vec_CC[14], 2.5e-7, "CC")
text(p.vec_C[14], 2e-7, "C")
text(p.vec_FP[14], 1.5e-7, "FP")
text(p.vec_PG[14], 1e-7, "PG")
text(p.vec_WT[14], .5e-7, "WT")

plot(results$new_par[144:150], results$delta[144:150],
     xlab = expression(mu[diam]), 
     ylab=expression(Delta*lambda), main = "Distribution of Diameter")
abline(v = results$par_original[144], col="red")
abline(h =0)
abline(v = p.vec_BC[15], lty=3, col=cols[1])
text(p.vec_BC[15], 0.025, "BC")
text(p.vec_CC[15], 0.0225, "CC")
text(p.vec_C[15], 0.02, "C")
text(p.vec_FP[15], 0.0175, "FP")
text(p.vec_PG[15], 0.015, "PG")
text(p.vec_PG[15], 0.0125, "WT")

plot(results$new_par[151:160], results$delta[151:160],
     xlab = expression(sd[diam]^2), 
     ylab=expression(Delta*lambda), main = "")
abline(v = results$par_original[151], col="red")
abline(h =0)
abline(v = p.vec_BC[16], lty=3, col=cols[1])
text(p.vec_BC[16], 0.012, "BC")
text(p.vec_CC[16], 0.011, "CC")
text(p.vec_C[16], 0.01, "C")
text(p.vec_FP[16], 0.009, "FP")
text(p.vec_PG[16], 0.008, "PG")
text(p.vec_PG[16], 0.007, "WT")

plot(results$new_par[163:170], results$delta[163:170],
     xlab = expression(mu[height]), 
     ylab=expression(Delta*lambda), main = "Distribution of Height")
abline(v = results$par_original[163], col="red")
abline(h =0)
abline(v = p.vec_BC[17], lty=3, col=cols[1])
text(p.vec_BC[17], 2e-9, "BC")
text(p.vec_CC[17], 1.5e-9, "CC")
text(p.vec_C[17], 1e-9, "C")
text(p.vec_FP[17], 5e-10, "FP")
text(p.vec_PG[17], 0, "PG")
text(p.vec_PG[17], -5e-10, "WT")

plot(results$new_par[171:180], results$delta[171:180],
     xlab = expression(sd[height]^2), 
     ylab=expression(Delta*lambda), main = "")
abline(v = results$par_original[171], col="red")
abline(h =0)
abline(v = p.vec_BC[18], lty=3, col=cols[1])
text(p.vec_BC[18], -0.5e-8, "BC")
text(p.vec_CC[18], -1e-8, "CC")
text(p.vec_C[18], -1.5e-8, "C")
text(p.vec_FP[18], -2e-8, "FP")
text(p.vec_PG[18], -2.5e-8, "PG")
text(p.vec_PG[18], -3e-8, "WT")

dev.copy(pdf, 'maturation_sens.pdf')
dev.off()

# -Adult Survival  #####
#logit(Ss) = B0 + B1*diam + B2*height
x11()
par(mfrow=c(2,2))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=F)
text(3,8, "Adult Survival", cex=1)
expression('title'[2])
text(4, 7, expression("Logit (S"[A]*") = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
     cex=0.8)

plot(results$new_par[181:190], results$delta[181:190],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda), xlim=c(0,2.5 ))
abline(v = results$par_original[19], col="red")
abline(h =0)
abline(v = p.vec_BC[19], lty=3, col=cols[1])
text(p.vec_BC[19], 4.5e-08, "BC")
abline(v = p.vec_CC[19], lty=3, col=cols[2])
text(p.vec_CC[19], 3e-08, "CC")
abline(v = p.vec_C[19],  lty=3, col=cols[3])
text(p.vec_C[19], 4.5e-08, "C")
abline(v = p.vec_FP[19], lty=3, col=cols[4])
text(p.vec_FP[19], 4.5e-08, "FP")
abline(v = p.vec_PG[19], lty=3, col=cols[5])
text(p.vec_PG[19], 1.5e-8, "PG")
abline(v = p.vec_WT[19], lty=3, col=cols[6])
text(p.vec_WT[19], 4.5e-8, "WT")

plot(results$new_par[191:200], results$delta[191:200],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[191], col="red")
abline(h =0)
abline(v = p.vec_BC[20], lty=3, col=cols[1])
text(p.vec_BC[20], 2e-8, "BC")
text(p.vec_CC[20], 1.25e-8, "CC")
text(p.vec_C[20], 5e-09, "C")
text(p.vec_FP[20], -2.5e-9, "FP")
text(p.vec_PG[20], -1e-8, "PG")
text(p.vec_WT[20], -1.75e-8, "WT")

plot(results$new_par[201:210], results$delta[201:210],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[201], col="red")
abline(h =0)

abline(v = p.vec_BC[21], lty=3, col=cols[1])
text(p.vec_BC[21], 2e-8, "BC")
text(p.vec_CC[21], 1.25e-8, "CC")
text(p.vec_C[21], 5e-09, "C")
text(p.vec_FP[21], -2.5e-9, "FP")
text(p.vec_PG[21], -1e-8, "PG")
text(p.vec_WT[21], -1.75e-8, "WT")
dev.copy(pdf, 'adult_survival_sens.pdf')
dev.off()

#-- Adult Growth #####


#Adult Survival logit(Ss)  = B0 + B1*diam + B2*height
x11(width=11)
par(mfrow=c(2,4))

#Diameter
plot(results$new_par[211:220], results$delta[211:220],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda),
     xlim=c(-0.4, 9.1))
abline(v = results$par_original[211], col="red")
abline(h =0)
abline(v = p.vec_BC[22], lty=3, col=cols[1])
text(p.vec_BC[22], -1e-5, "BC")
abline(v = p.vec_CC[22], lty=3, col=cols[2])
text(p.vec_CC[22], -1e-5, "CC")
abline(v = p.vec_C[22],  lty=3, col=cols[3])
text(p.vec_C[22], -1e-5, "C")
abline(v = p.vec_FP[22], lty=3, col=cols[4])
text(p.vec_FP[22], -1e-5, "FP")
abline(v = p.vec_PG[22], lty=3, col=cols[5])
text(p.vec_PG[22], -2e-5, "PG")
abline(v = p.vec_WT[22], lty=3, col=cols[6])
text(p.vec_WT[22], -2e-5, "WT")

plot(results$new_par[221:230], results$delta[221:230],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
mtext(side=3, expression(mu[diam]*" = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
      cex=0.8)
abline(v = results$par_original[221], col="red")
abline(h =0)
abline(v = p.vec_BC[23], lty=3, col=cols[1])
text(p.vec_BC[23],-0.03, "BC")
text(p.vec_CC[23], -0.045, "CC")
text(p.vec_C[23], -0.06, "C")
text(p.vec_FP[23], -0.075, "FP")
text(p.vec_PG[23], -0.09, "PG")
text(p.vec_WT[23], -0.105, "WT")

plot(results$new_par[231:240], results$delta[231:240],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[231], col="red")
abline(h =0)

abline(v = p.vec_BC[24], lty=3, col=cols[1])
text(p.vec_BC[24], 2e-4, "BC")
text(p.vec_CC[24], 1e-4, "CC")
text(p.vec_C[24], 0, "C")
text(p.vec_FP[24], -1e-4, "FP")
text(p.vec_PG[24], -2e-4, "PG")
text(p.vec_WT[24], -3e-4, "WT")

plot(results$new_par[241:250], results$delta[241:250],
     xlab = expression(sigma[diam]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[241], col="red")
abline(h =0)

abline(v = p.vec_BC[25], lty=3, col=cols[1])
text(p.vec_BC[25], -0.01, "BC")
text(p.vec_CC[25], -0.015, "CC")
text(p.vec_C[25], -0.02, "C")
text(p.vec_FP[25], -0.025, "FP")
text(p.vec_PG[25], -0.03, "PG")
text(p.vec_WT[25], -0.035, "WT")

#Height
plot(results$new_par[251:260], results$delta[251:260],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda),
     xlim=c(12, 65))
abline(v = results$par_original[251], col="red")
abline(h =0)
abline(v = p.vec_BC[26], lty=3, col=cols[1])
text(p.vec_BC[26], 2e-4, "BC")
abline(v = p.vec_CC[26], lty=3, col=cols[2])
text(p.vec_CC[26], 2e-04, "CC")
abline(v = p.vec_C[26],  lty=3, col=cols[3])
text(p.vec_C[26], 2e-4, "C")
abline(v = p.vec_FP[26], lty=3, col=cols[4])
text(p.vec_FP[26], 2e-4, "FP")
abline(v = p.vec_PG[26], lty=3, col=cols[5])
text(p.vec_PG[26], 2e-04, "PG")
abline(v = p.vec_WT[26], lty=3, col=cols[6])
text(p.vec_WT[26], 2e-04, "WT")

plot(results$new_par[261:270], results$delta[261:270],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
mtext(side=3, expression(mu[height]*" = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
      cex=0.8)
abline(v = results$par_original[261], col="red")
abline(h =0)
abline(v = p.vec_BC[27], lty=3, col=cols[1])
text(p.vec_BC[27],-5e-4, "BC")
text(p.vec_CC[27], -7e-4, "CC")
text(p.vec_C[27], -10e-4, "C")
text(p.vec_FP[27], -13e-4, "FP")
text(p.vec_PG[27], -16e-4, "PG")
text(p.vec_WT[27], -19e-4, "WT")

plot(results$new_par[271:280], results$delta[271:280],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[271], col="red")
abline(h =0)

abline(v = p.vec_BC[28], lty=3, col=cols[1])
text(p.vec_BC[28], -0.005, "BC")
text(p.vec_CC[28], -0.006, "CC")
text(p.vec_C[28], -0.007, "C")
text(p.vec_FP[28], -0.008, "FP")
text(p.vec_PG[28], -0.009, "PG")
text(p.vec_WT[28], -0.010, "WT")

plot(results$new_par[281:290], results$delta[281:290],
     xlab = expression(sigma[height]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[281], col="red")
abline(h =0)

abline(v = p.vec_BC[29], lty=3, col=cols[1])
text(p.vec_BC[29], -0.005, "BC")
text(p.vec_CC[29], -0.006, "CC")
text(p.vec_C[29], -0.007, "C")
text(p.vec_FP[29], -0.008, "FP")
text(p.vec_PG[29], -0.009, "PG")
text(p.vec_WT[29], -0.010, "WT")

dev.copy(pdf, 'adult_growth_sens.pdf')
dev.off()

##- Reproduction #####

x11(width=11)
par(mfrow=c(2, 5))
plot(results$new_par[291:300], results$delta[291:300],
     xlab = expression(B[0]), 
     ylab=expression(Delta*lambda), 
     xlim=c(-4.5, 1.6),
     main=expression("logit("*p[f]*") = B"[0]*"+ B"[1]*" * height"))
abline(v = results$par_original[291], col="red")
abline(h =0)

abline(v = p.vec_BC[30], lty=3, col=cols[1])
text(p.vec_BC[30], -1.5e-8, "BC")
abline(v = p.vec_CC[30], lty=3, col=cols[2])
text(p.vec_CC[30], -2e-8, "CC")
abline(v = p.vec_C[30], lty=3, col=cols[3])
text(p.vec_C[30], -1.5e-8, "C")
abline(v=p.vec_FP[30], lty=3, col=cols[4])
text(p.vec_FP[30], -2.5e-8, "FP")
abline(v=p.vec_PG[30], lty=3, col=cols[5])
text(p.vec_PG[30], -1.5e-8, "PG")
abline(v=p.vec_WT[30], lty=3, col=cols[6])
text(p.vec_WT[30], -1.5e-8, "WT")

plot(results$new_par[301:310], results$delta[301:310],
     xlab = expression(B[1]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[301], col="red")
abline(h =0)

abline(v = p.vec_BC[31], lty=3, col=cols[1])
text(p.vec_BC[31], -1.5e-8, "BC")
text(p.vec_CC[31], -2.5e-8, "CC")
text(p.vec_C[31], -3.5e-8, "C")
text(p.vec_FP[31], -4.5e-8, "FP")
text(p.vec_PG[31], -5.5e-8, "PG")
text(p.vec_WT[31], -6.5e-8, "WT")

plot(results$new_par[311:320], results$delta[311:320],
     xlab = expression(B[0]), 
     ylab=expression(Delta*lambda),
     xlim=c(1.5, 7),
     main=expression("logit(F) = B"[0]*x^2))
abline(v = results$par_original[311], col="red")
abline(h =0)

abline(v = p.vec_BC[32], lty=3, col=cols[1])
text(p.vec_BC[32], -.5e-7, "BC")
abline(v = p.vec_CC[32], lty=3, col=cols[2])
text(p.vec_CC[32], -.5e-7, "CC")
abline(v = p.vec_C[32], lty=3, col=cols[3])
text(p.vec_C[32], -.5e-7, "C")
abline(v=p.vec_FP[32], lty=3, col=cols[4])
text(p.vec_FP[32], -.5e-7, "FP")
abline(v=p.vec_PG[32], lty=3, col=cols[5])
text(p.vec_PG[32], -.5e-7, "PG")
abline(v=p.vec_WT[32], lty=3, col=cols[6])
text(p.vec_WT[32], -.5e-7, "WT")

plot(results$new_par[321:330], results$delta[321:330],
     xlab = expression(mu[diam]), 
     ylab=expression(Delta*lambda),
     main="Seedling diameter distribution")
abline(v = results$par_original[321], col="red")
abline(h =0)

abline(v = p.vec_BC[33], lty=3, col=cols[1])
text(p.vec_BC[33], 4e-7, "BC")
text(p.vec_CC[33], 3.5e-7, "CC")
text(p.vec_C[33], 3e-7, "C")
text(p.vec_FP[33], 2.5e-7, "FP")
text(p.vec_PG[33], 2e-7, "PG")
text(p.vec_WT[33], 1.5e-7, "WT")

plot(results$new_par[331:340], results$delta[331:340],
     xlab = expression(sigma[diam]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[331], col="red")
abline(h =0)

abline(v = p.vec_BC[34], lty=3, col=cols[1])
text(p.vec_BC[34], 4e-8, "BC")
text(p.vec_CC[34], 3.5e-8, "CC")
text(p.vec_C[34], 3e-8, "C")
text(p.vec_FP[34], 2.5e-8, "FP")
text(p.vec_PG[34], 2e-8, "PG")
text(p.vec_WT[34], 1.5e-8, "WT")

plot(results$new_par[341:350], results$delta[341:350],
     xlab = expression(mu[height]), 
     ylab=expression(Delta*lambda),
     main="Seedling height distribution")
abline(v = results$par_original[341], col="red")
abline(h =0)

abline(v = p.vec_BC[35], lty=3, col=cols[1])
text(p.vec_BC[35], 3e-7, "BC")
text(p.vec_CC[35], 2.5e-7, "CC")
text(p.vec_C[35], 2e-7, "C")
text(p.vec_FP[35], 1.5e-7, "FP")
text(p.vec_PG[35], 1e-7, "PG")
text(p.vec_WT[35], .5e-7, "WT")

plot(results$new_par[351:360], results$delta[351:360],
     xlab = expression(sigma[height]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[351], col="red")
abline(h =0)

abline(v = p.vec_BC[36], lty=3, col=cols[1])
text(p.vec_BC[36], 4e-8, "BC")
text(p.vec_CC[36], 3.5e-8, "CC")
text(p.vec_C[36], 3e-8, "C")
text(p.vec_FP[36], 2.5e-8, "FP")
text(p.vec_PG[36], 2e-8, "PG")
text(p.vec_WT[36], 1.5e-8, "WT")


plot(results$new_par[361:370], results$delta[361:370],
     xlab = expression(tau[1]), 
     ylab=expression(Delta*lambda),
     main="Pre-dispersal seed survival")
abline(v = results$par_original[361], col="red")
abline(h =0)

abline(v = p.vec_BC[37], lty=3, col=cols[1])
text(p.vec_BC[37], 2e-7, "BC")
text(p.vec_CC[37], 1.8e-7, "CC")
text(p.vec_C[37], 1.6e-7, "C")
text(p.vec_FP[37], 1.4e-7, "FP")
text(p.vec_PG[37], 1.2e-7, "PG")
text(p.vec_WT[37], 1.e-7, "WT")

plot(results$new_par[371:380], results$delta[371:380],
     xlab = expression(delta), 
     ylab=expression(Delta*lambda),
     main="P(dispersal)")
abline(v = results$par_original[371], col="red")
abline(h =0)

abline(v = p.vec_BC[38], lty=3, col=cols[1])
text(p.vec_BC[38], 2e-7, "BC")
text(p.vec_CC[38], 1.8e-7, "CC")
text(p.vec_C[38], 1.6e-7, "C")
text(p.vec_FP[38], 1.4e-7, "FP")
text(p.vec_PG[38], 1.2e-7, "PG")
text(p.vec_WT[38], 1.e-7, "WT")

plot(results$new_par[381:390], results$delta[381:390],
     xlab = expression(tau[2]), 
     ylab=expression(Delta*lambda),
     main="Post-dispersal seed survival")
abline(v = results$par_original[381], col="red")
abline(h =0)

abline(v = p.vec_BC[39], lty=3, col=cols[1])
text(p.vec_BC[39], 2e-7, "BC")
text(p.vec_CC[39], 1.8e-7, "CC")
text(p.vec_C[39], 1.6e-7, "C")
text(p.vec_FP[39], 1.4e-7, "FP")
text(p.vec_PG[39], 1.2e-7, "PG")
text(p.vec_WT[39], 1.e-7, "WT")
dev.copy(pdf, 'reproduction_sens.pdf')
dev.off()


#Some additional params 
thing1 <- c("min", NA, 22, -0.4, NA, NA)
thing2<-c("max", NA, 22, 9.1, NA, NA)
thing3 <- c("min", NA, 26, 12, NA, NA)
thing4 <- c("max", NA, 26, 64, NA, NA)
#results <- rbind(results, thing1, thing2, thing3, thing4)


#Extremes of p.vec[22]
for(i in 391:nrow(results)) {
  cat("Param: ", i, "of ", nrow(results), ": \n");
  p.vec_new <- p.vec_overall
  p.vec_new[as.numeric(results$pos[i])] <- as.numeric(results$new_par[i])
  results$lambda[i] <- getLambda(p.vec_new)
}

results$delta <- as.numeric(results$lambda) - lam.stable_overall

results$new_par <- as.numeric(results$new_par)
## REDO ADULT GROWTH #####
x11(width=11)
par(mfrow=c(2,4))

#Diameter
plot(results$new_par[211:220], results$delta[211:220],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda),
     xlim=c(-1, 9.5), ylim=c(-0.002, 0.0003))
points(results$new_par[391], results$delta[391])
points(results$new_par[392], results$delta[392])
abline(v = results$par_original[211], col="red")
abline(h =0)
abline(v = p.vec_BC[22], lty=3, col=cols[1])
#text(p.vec_BC[22], -1e-5, "BC")
abline(v = p.vec_CC[22], lty=3, col=cols[2])
#text(p.vec_CC[22], -1e-5, "CC")
abline(v = p.vec_C[22],  lty=3, col=cols[3])
#text(p.vec_C[22], -1e-5, "C")
abline(v = p.vec_FP[22], lty=3, col=cols[4])
#text(p.vec_FP[22], -1e-5, "FP")
abline(v = p.vec_PG[22], lty=3, col=cols[5])
#text(p.vec_PG[22], -2e-5, "PG")
abline(v = p.vec_WT[22], lty=3, col=cols[6])
#text(p.vec_WT[22], -2e-5, "WT")

plot(results$new_par[221:230], results$delta[221:230],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
mtext(side=3, expression(mu[diam]*" = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
      cex=0.8)
abline(v = results$par_original[221], col="red")
abline(h =0)
abline(v = p.vec_BC[23], lty=3, col=cols[1])
text(p.vec_BC[23],-0.03, "BC")
text(p.vec_CC[23], -0.045, "CC")
text(p.vec_C[23], -0.06, "C")
text(p.vec_FP[23], -0.075, "FP")
text(p.vec_PG[23], -0.09, "PG")
text(p.vec_WT[23], -0.105, "WT")

plot(results$new_par[231:240], results$delta[231:240],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[231], col="red")
abline(h =0)

abline(v = p.vec_BC[24], lty=3, col=cols[1])
text(p.vec_BC[24], 2e-4, "BC")
text(p.vec_CC[24], 1e-4, "CC")
text(p.vec_C[24], 0, "C")
text(p.vec_FP[24], -1e-4, "FP")
text(p.vec_PG[24], -2e-4, "PG")
text(p.vec_WT[24], -3e-4, "WT")

plot(results$new_par[241:250], results$delta[241:250],
     xlab = expression(sigma[diam]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[241], col="red")
abline(h =0)

abline(v = p.vec_BC[25], lty=3, col=cols[1])
text(p.vec_BC[25], -0.01, "BC")
text(p.vec_CC[25], -0.015, "CC")
text(p.vec_C[25], -0.02, "C")
text(p.vec_FP[25], -0.025, "FP")
text(p.vec_PG[25], -0.03, "PG")
text(p.vec_WT[25], -0.035, "WT")

#Height
plot(results$new_par[251:260], results$delta[251:260],
     xlab = expression("B"[0]), 
     ylab=expression(Delta*lambda),
     xlim=c(12, 65),
     ylim=c(-0.00075, 0.0004))
points(results$new_par[393], results$delta[393])
points(results$new_par[394], results$delta[394])
abline(v = results$par_original[251], col="red")
abline(h =0)
abline(v = p.vec_BC[26], lty=3, col=cols[1])
text(p.vec_BC[26], 2e-4, "BC")
abline(v = p.vec_CC[26], lty=3, col=cols[2])
text(p.vec_CC[26], 2e-04, "CC")
abline(v = p.vec_C[26],  lty=3, col=cols[3])
text(p.vec_C[26], 2e-4, "C")
abline(v = p.vec_FP[26], lty=3, col=cols[4])
text(p.vec_FP[26], 2e-4, "FP")
abline(v = p.vec_PG[26], lty=3, col=cols[5])
text(p.vec_PG[26], 2e-04, "PG")
abline(v = p.vec_WT[26], lty=3, col=cols[6])
text(p.vec_WT[26], 2e-04, "WT")

plot(results$new_par[261:270], results$delta[261:270],
     xlab = expression("B"[1]), 
     ylab=expression(Delta*lambda))
mtext(side=3, expression(mu[height]*" = B"[0]*" + B"[1]*"*diam + B"[2]*" * height"),
      cex=0.8)
abline(v = results$par_original[261], col="red")
abline(h =0)
abline(v = p.vec_BC[27], lty=3, col=cols[1])
text(p.vec_BC[27],-5e-4, "BC")
text(p.vec_CC[27], -7e-4, "CC")
text(p.vec_C[27], -10e-4, "C")
text(p.vec_FP[27], -13e-4, "FP")
text(p.vec_PG[27], -16e-4, "PG")
text(p.vec_WT[27], -19e-4, "WT")

plot(results$new_par[271:280], results$delta[271:280],
     xlab = expression("B"[2]), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[271], col="red")
abline(h =0)

abline(v = p.vec_BC[28], lty=3, col=cols[1])
text(p.vec_BC[28], -0.005, "BC")
text(p.vec_CC[28], -0.006, "CC")
text(p.vec_C[28], -0.007, "C")
text(p.vec_FP[28], -0.008, "FP")
text(p.vec_PG[28], -0.009, "PG")
text(p.vec_WT[28], -0.010, "WT")

plot(results$new_par[281:290], results$delta[281:290],
     xlab = expression(sigma[height]^2), 
     ylab=expression(Delta*lambda))
abline(v = results$par_original[281], col="red")
abline(h =0)

abline(v = p.vec_BC[29], lty=3, col=cols[1])
text(p.vec_BC[29], -0.005, "BC")
text(p.vec_CC[29], -0.006, "CC")
text(p.vec_C[29], -0.007, "C")
text(p.vec_FP[29], -0.008, "FP")
text(p.vec_PG[29], -0.009, "PG")
text(p.vec_WT[29], -0.010, "WT")

dev.copy(pdf, 'adult_growth_sens.pdf')
dev.off()

save(results, file="results_sens_pvec.RData")
