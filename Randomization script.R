library(Matrix)
library(MASS)
library(mvtnorm)
library(msm)
library(boot)

n.samp <- 200 #for now, to test
demog<-read.csv("demography_15_clean.csv", head=T)
load("Stems.RData")
Stems <- subset(Stems, Stems$Site != "Holiday Park")
lam.boot_random <- array(NA, dim=c(n.samp, 6))

start_time <- proc.time()
for (boots in 1:n.samp) {
  temp <- demog
  #reshuffle the labels for site
  temp$Site <- sample(temp$Site)
  temp2 <-Stems 
  temp2$Site <- sample(temp2$Site)
  
  
  ### Part I: Running the statistical models to determine parameters #####
  
  p.vec_BC<-rep(0, 39)
  p.vec_C <- rep(0, 39)
  p.vec_FP <- rep(0, 39)
  p.vec_WT <- rep(0, 39)
  p.vec_PG <- rep(0, 39)
  p.vec_CC <- rep(0, 39)
  
  
  #Gather together:
  Death_Cause<-c(temp$Death_Cause_2, temp$Death_Cause_3, temp$Death_Cause_4, temp$Death_Cause_5, temp$Death_Cause_6, temp$Death_Cause_7)
  Status<-c(temp$Status_1, temp$Status_2, temp$Status_3, temp$Status_4, temp$Status_5, temp$Status_6, temp$Status_7)
  Diams<-cbind(temp$Diam_1, temp$Diam_2, temp$Diam_3, temp$Diam_4, temp$Diam_5, temp$Diam_6, temp$Diam_7) 
  row.names(Diams)<-temp$Plant_ID
  Heights<-cbind(temp$Height_1, temp$Height_2, temp$Height_3, temp$Height_4, temp$Height_5, temp$Height_6, temp$Height_7)
  row.names(Heights)<-temp$Plant_ID
  
  #Stack data:
  Diams_pooled<-rbind( cbind(Diams[,1], Diams[,2]), cbind(Diams[,2], Diams[,3]), cbind(Diams[, 3], Diams[,4]), cbind(Diams[,4], Diams[,5]), cbind(Diams[,5], Diams[,6]), cbind(Diams[, 6], Diams[, 7]))
  Diams_pooled<-data.frame(Diams_pooled[,1], Diams_pooled[,2])
  names(Diams_pooled)<-c("t", "t_plus_one")
  
  Sites_pooled<-rep(temp$Site, 6)
  ID_pooled<-rep(temp$Plant_ID, 6)
  Location<-rep(temp$Location, 6)
  
  Heights_pooled<-rbind( cbind(Heights[,1], Heights[,2]), cbind(Heights[,2], Heights[,3]), cbind(Heights[, 3], Heights[,4]), cbind(Heights[,4], Heights[,5]), cbind(Heights[,5], Heights[,6]), cbind(Heights[, 6], Heights[, 7]))
  Heights_pooled<-data.frame(Heights_pooled[,1], Heights_pooled[,2])
  names(Heights_pooled)<-c("t", "t_plus_one")
  
  ###Use status matrix to get survival matrix
  statmat<-as.matrix(temp[, c(6, 14, 22, 30, 38, 46, 54)])
  survmat<-statmat
  survmat<- ifelse(survmat=="prestudy"|survmat=="missing"|
                     survmat=="still_dead", NA, survmat)
  survmat<- ifelse(survmat =="tagged"|survmat=="alive", 1, 0)
  
  survmat_pooled<-rbind( cbind(survmat[,1], survmat[,2]), cbind(survmat[,2], survmat[,3]), cbind(survmat[, 3], survmat[,4]), cbind(survmat[,4], survmat[,5]), cbind(survmat[,5], survmat[,6]), cbind(survmat[, 6], survmat[, 7]))
  names(survmat_pooled)<-c("surv_t", "surv_t_plus_one")
  
  survmat_pooled_bysize<-cbind(Diams_pooled[, 1], survmat_pooled[,2])
  names(survmat_pooled_bysize)<-c("size_t", "surv_t_plus_one")
  
  #Use reproduction status matrix to get reproduction matrix
  repstatmat<-as.matrix(temp[, c(8, 16, 24, 32, 40, 48, 56)])
  repmat<-repstatmat
  repmat<- ifelse(repmat =="female", 1, 0)
  repmat_pooled<-rbind( cbind(repmat[,1], repmat[,2]), cbind(repmat[,2], repmat[,3]), cbind(repmat[, 3], repmat[,4]), cbind(repmat[,4], repmat[,5]), cbind(repmat[,5], repmat[,6]), cbind(repmat[, 6], 
                                                                                                                                                                                          repmat[, 7]))
  
  pooled<-cbind(ID_pooled, Sites_pooled, Location, Death_Cause, Diams_pooled, Heights_pooled, repmat_pooled, survmat_pooled)
  names(pooled)<-c("ID", "Site","Location", "Death_Cause", "Diameter_t", "Diameter_tplus1", "Height_t", "Height_tplus1", "Rep_t", "Rep_tplus1", "Surv_t", "Surv_tplus1")
  
  #####Cleaning up the data
  
  #Remove individual deaths that were sprayed or fire or grinding
  restrict<-subset(pooled, is.na(Death_Cause))
  
  #Remove outliers (very large diameters)
  restrict2<-subset(restrict, Diameter_t<800)
  #Remove problematic individual
  restrict2 <- subset(restrict2, restrict2$ID != 5)
  
  #Divide dataset into size domains (seedlings and larges)
  seedlings<-subset(restrict2, restrict2$Diameter_t<1.6)
  seedlings<-subset(seedlings, seedlings$Height_t<16)
  larges<-subset(restrict2, restrict2$Diameter_t>1.6)
  larges<-subset(larges, larges$Height_t>16)
  
  # **Seedling Survival #####
  s1<-glm(seedlings$Surv_tplus1~ seedlings$Diameter_t + seedlings$Height_t+seedlings$Site, family=binomial)
  
  p.vec_BC[1]<-s1$coeff[1]
  p.vec_BC[2]<-s1$coeff[2]
  p.vec_BC[3]<-s1$coeff[3]
  
  p.vec_CC[1]<-s1$coeff[1]+s1$coeff[4]
  p.vec_CC[2]<-s1$coeff[2]
  p.vec_CC[3]<-s1$coeff[3]
  
  p.vec_C[1]<-s1$coeff[1] + s1$coeff[5]
  p.vec_C[2]<-s1$coeff[2]
  p.vec_C[3]<-s1$coeff[3]
  
  p.vec_FP[1] <- s1$coeff[1] + s1$coeff[6]
  p.vec_FP[2] <- s1$coeff[2]
  p.vec_FP[3] <- s1$coeff[3]
  
  p.vec_PG[1] <- s1$coeff[1] + s1$coeff[7]
  p.vec_PG[2] <- s1$coeff[2]
  p.vec_PG[3] <- s1$coeff[3]
  
  p.vec_WT[1] <- s1$coeff[1] + s1$coeff[8]
  p.vec_WT[2] <- s1$coeff[2]
  p.vec_WT[3] <- s1$coeff[3]
  
  # ** Seedling Growth #####
  # (Create seedlings2 to only use individuals that REMAIN as seedlings: 
  seedlings2<-subset(seedlings, seedlings$Diameter_tplus1<1.6)
  seedlings2<-subset(seedlings2, seedlings2$Height_tplus1<16)
  
  #To model variance in D1 growth for diameter and height, use separate linear models and store the variance: 
  #First, model variance in growth for future diameter
  g1_diam<-lm(seedlings2$Diameter_tplus1 ~ seedlings2$Diameter_t+seedlings2$Height_t + seedlings2$Site)
  
  
  #Then, model variance in growth for future height
  g1_height<-lm(seedlings2$Height_tplus1 ~ seedlings2$Diameter_t + seedlings2$Height_t + seedlings2$Site)
  
  
  #By Site: 
  g1<-manova(cbind(seedlings2$Diameter_tplus1, seedlings2$Height_tplus1) ~ seedlings2$Diameter_t+seedlings2$Height_t+seedlings2$Site)
  
  
  #Capture parameters associated with predicting future diameter 
  p.vec_BC[4]<-g1$coeff[1,1]
  p.vec_BC[5]<-g1$coeff[2,1]
  p.vec_BC[6]<-g1$coeff[3,1]
  p.vec_BC[7]<-summary(g1_diam)$sigma
  
  p.vec_CC[4]<-g1$coeff[1,1] + g1$coeff[4, 1]
  p.vec_CC[5]<-g1$coeff[2,1]
  p.vec_CC[6]<-g1$coeff[3,1]
  p.vec_CC[7]<-summary(g1_diam)$sigma
  
  p.vec_C[4]<-g1$coeff[1,1] + g1$coeff[5, 1]
  p.vec_C[5]<-g1$coeff[2,1]
  p.vec_C[6]<-g1$coeff[3,1]
  p.vec_C[7]<-summary(g1_diam)$sigma
  
  p.vec_FP[4]<-g1$coeff[1,1] + g1$coeff[6, 1]
  p.vec_FP[5]<-g1$coeff[2,1]
  p.vec_FP[6]<-g1$coeff[3,1]
  p.vec_FP[7]<-summary(g1_diam)$sigma
  
  p.vec_PG[4]<-g1$coeff[1,1] + g1$coeff[7, 1]
  p.vec_PG[5]<-g1$coeff[2,1]
  p.vec_PG[6]<-g1$coeff[3,1]
  p.vec_PG[7]<-summary(g1_diam)$sigma
  
  p.vec_WT[4]<-g1$coeff[1,1] + g1$coeff[8, 1]
  p.vec_WT[5]<-g1$coeff[2,1]
  p.vec_WT[6]<-g1$coeff[3,1]
  p.vec_WT[7]<-summary(g1_diam)$sigma
  
  
  #Grab parameters associated with predicting future height
  p.vec_BC[8]<-g1$coeff[1,2]
  p.vec_BC[9]<-g1$coeff[2,2]
  p.vec_BC[10]<-g1$coeff[3,2]
  p.vec_BC[11]<-summary(g1_height)$sigma
  
  p.vec_CC[8]<-g1$coeff[1,2] + g1$coeff[4,2]
  p.vec_CC[9]<-g1$coeff[2,2]
  p.vec_CC[10]<-g1$coeff[3,2]
  p.vec_CC[11]<-summary(g1_height)$sigma
  
  p.vec_C[8]<-g1$coeff[1,2] + g1$coeff[5,2]
  p.vec_C[9]<-g1$coeff[2,2]
  p.vec_C[10]<-g1$coeff[3,2]
  p.vec_C[11]<-summary(g1_height)$sigma
  
  p.vec_FP[8]<-g1$coeff[1,2] + g1$coeff[6,2]
  p.vec_FP[9]<-g1$coeff[2,2]
  p.vec_FP[10]<-g1$coeff[3,2]
  p.vec_FP[11]<-summary(g1_height)$sigma
  
  p.vec_PG[8]<-g1$coeff[1,2] + g1$coeff[7,2]
  p.vec_PG[9]<-g1$coeff[2,2]
  p.vec_PG[10]<-g1$coeff[3,2]
  p.vec_PG[11]<-summary(g1_height)$sigma
  
  p.vec_WT[8]<-g1$coeff[1,2] + g1$coeff[8,2]
  p.vec_WT[9]<-g1$coeff[2,2]
  p.vec_WT[10]<-g1$coeff[3,2]
  p.vec_WT[11]<-summary(g1_height)$sigma
  
  # ** Maturation #####
  grad_status<-rep(NA, length(seedlings$ID))
  for (i in 1:length(seedlings$ID)) {
    if (is.na(seedlings$Diameter_tplus1[i])) {
      grad_status[i]<-NA
    }
    else if (is.na(seedlings$Height_tplus1[i])) {
      grad_status[i]<-NA
    }
    else if (seedlings$Diameter_tplus1[i]>1.6& seedlings$Height_tplus1[i]>16) { grad_status[i]<-1}
    else {grad_status[i]<-0}
  }
  
  seedlings<-cbind(seedlings, grad_status)
  
  
  #Model distribution of sizes of graduates 
  
  graduates<-subset(seedlings, seedlings$grad_status==1)
  
  
  #Get distribution of new recruits sizes (these are not modeled by biotype) 
  mu_grad_diam<-mean(graduates$Diameter_tplus1)
  sd_grad_diam<-sd(graduates$Diameter_tplus1)
  
  mu_grad_height<-mean(graduates$Height_tplus1)
  sd_grad_height<-sd(graduates$Height_tplus1)
  
  
  p.vec_BC[15]<-mu_grad_diam
  p.vec_BC[16]<-sd_grad_diam
  p.vec_BC[17]<-mu_grad_height
  p.vec_BC[18]<-sd_grad_height
  
  p.vec_CC[15]<-mu_grad_diam
  p.vec_CC[16]<-sd_grad_diam
  p.vec_CC[17]<-mu_grad_height
  p.vec_CC[18]<-sd_grad_height
  
  p.vec_C[15]<-mu_grad_diam
  p.vec_C[16]<-sd_grad_diam
  p.vec_C[17]<-mu_grad_height
  p.vec_C[18]<-sd_grad_height
  
  p.vec_FP[15]<-mu_grad_diam
  p.vec_FP[16]<-sd_grad_diam
  p.vec_FP[17]<-mu_grad_height
  p.vec_FP[18]<-sd_grad_height
  
  p.vec_PG[15]<-mu_grad_diam
  p.vec_PG[16]<-sd_grad_diam
  p.vec_PG[17]<-mu_grad_height
  p.vec_PG[18]<-sd_grad_height
  
  p.vec_WT[15]<-mu_grad_diam
  p.vec_WT[16]<-sd_grad_diam
  p.vec_WT[17]<-mu_grad_height
  p.vec_WT[18]<-sd_grad_height
  
  
  #(b) Model graduation by site
  m <-glm(seedlings$grad_status~ seedlings$Diameter_t + seedlings$Height_t + seedlings$Site, family=binomial)
  
  p.vec_BC[12]<-m$coeff[1]
  p.vec_BC[13]<-m$coeff[2]
  p.vec_BC[14]<-m$coeff[3]
  
  p.vec_CC[12]<-m$coeff[1] + m$coeff[4]
  p.vec_CC[13]<-m$coeff[2]
  p.vec_CC[14]<-m$coeff[3]
  
  p.vec_C[12]<-m$coeff[1] + m$coeff[5]
  p.vec_C[13]<-m$coeff[2]
  p.vec_C[14]<-m$coeff[3]
  
  p.vec_FP[12]<-m$coeff[1] + m$coeff[6]
  p.vec_FP[13]<-m$coeff[2]
  p.vec_FP[14]<-m$coeff[3]
  
  p.vec_PG[12]<-m$coeff[1] + m$coeff[7]
  p.vec_PG[13]<-m$coeff[2]
  p.vec_PG[14]<-m$coeff[3]
  
  p.vec_WT[12]<-m$coeff[1] + m$coeff[8]
  p.vec_WT[13]<-m$coeff[2]
  p.vec_WT[14]<-m$coeff[3]
  
  # **D2 Survival #####
  s2<-glm(larges$Surv_tplus1~larges$Diameter_t + larges$Height_t +larges$Site, family=binomial)
  
  p.vec_BC[19]<-s2$coeff[1]
  p.vec_BC[20]<-s2$coeff[2]
  p.vec_BC[21]<-s2$coeff[3]
  
  p.vec_CC[19]<-s2$coeff[1] + s2$coeff[4]
  p.vec_CC[20]<-s2$coeff[2]
  p.vec_CC[21]<-s2$coeff[3]
  
  p.vec_C[19]<-s2$coeff[1] + s2$coeff[5]
  p.vec_C[20]<-s2$coeff[2]
  p.vec_C[21]<-s2$coeff[3]
  
  p.vec_FP[19]<-s2$coeff[1] + s2$coeff[6]
  p.vec_FP[20]<-s2$coeff[2]
  p.vec_FP[21]<-s2$coeff[3]
  
  p.vec_PG[19]<-s2$coeff[1] + s2$coeff[7]
  p.vec_PG[20]<-s2$coeff[2]
  p.vec_PG[21]<-s2$coeff[3]
  
  p.vec_WT[19]<-s2$coeff[1] + s2$coeff[8]
  p.vec_WT[20]<-s2$coeff[2]
  p.vec_WT[21]<-s2$coeff[3]
  
  # **D2 Growth #####
  g2_diam<-lm(larges$Diameter_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Site)
  
  #Model variance in growth in height: 
  g2_height<-lm(larges$Height_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Site)
  
  g2<-manova(cbind(larges$Diameter_tplus1, larges$Height_tplus1) ~ larges$Diameter_t+larges$Height_t+larges$Site)
  
  #Store parameters associated with growth of D2 individuals in diameter
  p.vec_BC[22]<-g2$coeff[1,1]
  p.vec_BC[23]<-g2$coeff[2,1]
  p.vec_BC[24]<-g2$coeff[3,1]
  
  p.vec_CC[22]<-g2$coeff[1,1] + g2$coeff[4,1]
  p.vec_CC[23]<-g2$coeff[2,1]
  p.vec_CC[24]<-g2$coeff[3,1]
  
  p.vec_C[22]<-g2$coeff[1,1] + g2$coeff[5,1]
  p.vec_C[23]<-g2$coeff[2,1]
  p.vec_C[24]<-g2$coeff[3,1]
  
  p.vec_FP[22]<-g2$coeff[1,1] + g2$coeff[6,1]
  p.vec_FP[23]<-g2$coeff[2,1]
  p.vec_FP[24]<-g2$coeff[3,1]
  
  p.vec_PG[22]<-g2$coeff[1,1] + g2$coeff[7,1]
  p.vec_PG[23]<-g2$coeff[2,1]
  p.vec_PG[24]<-g2$coeff[3,1]
  
  p.vec_WT[22]<-g2$coeff[1,1] + g2$coeff[8,1]
  p.vec_WT[23]<-g2$coeff[2,1]
  p.vec_WT[24]<-g2$coeff[3,1]
  
  p.vec_BC[25]<-summary(g2_diam)$sigma
  p.vec_CC[25]<-summary(g2_diam)$sigma
  p.vec_C[25]<-summary(g2_diam)$sigma
  p.vec_FP[25]<-summary(g2_diam)$sigma
  p.vec_PG[25]<-summary(g2_diam)$sigma
  p.vec_WT[25]<-summary(g2_diam)$sigma
  
  #Store parameters associated with modeling growth of D2 individuals in height
  p.vec_BC[26]<-g2$coeff[1,2]
  p.vec_BC[27]<-g2$coeff[2,2]
  p.vec_BC[28]<-g2$coeff[3,2]
  
  p.vec_CC[26]<-g2$coeff[1,2] + g2$coeff[4,2]
  p.vec_CC[27]<-g2$coeff[2,2]
  p.vec_CC[28]<-g2$coeff[3,2]
  
  p.vec_C[26]<-g2$coeff[1,2] + g2$coeff[5,2]
  p.vec_C[27]<-g2$coeff[2,2]
  p.vec_C[28]<-g2$coeff[3,2]
  
  p.vec_FP[26]<-g2$coeff[1,2] + g2$coeff[6,2]
  p.vec_FP[27]<-g2$coeff[2,2]
  p.vec_FP[28]<-g2$coeff[3,2]
  
  p.vec_PG[26]<-g2$coeff[1,2] + g2$coeff[7,2]
  p.vec_PG[27]<-g2$coeff[2,2]
  p.vec_PG[28]<-g2$coeff[3,2]
  
  p.vec_WT[26]<-g2$coeff[1,2] + g2$coeff[8,2]
  p.vec_WT[27]<-g2$coeff[2,2]
  p.vec_WT[28]<-g2$coeff[3,2]
  
  p.vec_BC[29]<-summary(g2_height)$sigma
  p.vec_CC[29]<-summary(g2_height)$sigma
  p.vec_C[29]<-summary(g2_height)$sigma
  p.vec_FP[29]<-summary(g2_height)$sigma
  p.vec_PG[29]<-summary(g2_height)$sigma
  p.vec_WT[29]<-summary(g2_height)$sigma
  
  # **Reproduction #####
  #(b) Model probability of being reproductive     (by site)
  p_f1<-glm(larges$Rep_tplus1 ~ larges$Diameter_t+larges$Height_t + larges$Site, family=binomial)
  
  #Diameter term is not significant, so remove it: 
  p_f<-glm(larges$Rep_tplus1 ~ larges$Height_t +larges$Site, family=binomial)
  
  
  p.vec_BC[30]<-p_f$coeff[1]
  p.vec_BC[31]<-p_f$coeff[2]
  
  p.vec_CC[30]<-p_f$coeff[1] + p_f$coeff[3]
  p.vec_CC[31]<-p_f$coeff[2]
  
  p.vec_C[30]<-p_f$coeff[1] + p_f$coeff[4]
  p.vec_C[31]<-p_f$coeff[2]
  
  p.vec_FP[30]<-p_f$coeff[1] + p_f$coeff[5]
  p.vec_FP[31]<-p_f$coeff[2]
  
  p.vec_PG[30]<-p_f$coeff[1] + p_f$coeff[6]
  p.vec_PG[31]<-p_f$coeff[2]
  
  p.vec_WT[30]<-p_f$coeff[1] + p_f$coeff[7]
  p.vec_WT[31]<-p_f$coeff[2]
  
  # **Fecundity #####
  
  
  #(b) Model fecundity (by biotype)
  #These values come from rerunning the allometry analysis using site instead of biotype. They are fixed because the data comes from a different data set. 
  
  
  
  x<-temp2$diam_base*10 #convert from cm to mm
  x <- x*x
  f<-lm(temp2$Seed_No ~ 0 + x:temp2$Site)
  summary(f)
  
  f_overall <- lm(temp2$Seed_No ~ 0 + x)
  summary(f_overall)
  
  
  #(b) Model fecundity (by biotype)
  p.vec_BC[32]<-f$coeff[1]
  p.vec_CC[32]<-f$coeff[2]
  p.vec_C[32]<-f$coeff[3]
  p.vec_FP[32] <- f$coeff[4]
  p.vec_PG[32] <- f$coeff[5]
  p.vec_WT[32] <- f$coeff[6]
  
  
  #Determine new recruits size
  recruit4_diam<-temp$Diam_4[temp$Status_4=="tagged"&
                               temp$Location=="SeedlingPlot"]
  recruit5_diam<-temp$Diam_5[temp$Status_5=="tagged"&
                               temp$Location=="SeedlingPlot"]
  recruit6_diam<-temp$Diam_6[temp$Status_6=="tagged"&
                               temp$Location=="SeedlingPlot"]
  recruits_diam<-c(recruit4_diam, recruit5_diam, recruit6_diam)
  
  recruit4_height<-temp$Height_4[temp$Status_4=="tagged"&
                                   temp$Location=="SeedlingPlot"]
  recruit5_height<-temp$Height_5[temp$Status_5=="tagged"&
                                   temp$Location=="SeedlingPlot"]
  recruit6_height<-temp$Height_6[temp$Status_6=="tagged"&
                                   temp$Location=="SeedlingPlot"]
  recruits_height<-c(recruit4_height, recruit5_height, recruit6_height)
  
  
  
  
  
  
  #####Summary statistics for just the 476 individuals with a diameter<1.6mm
  mu_diam<-mean(recruits_diam[recruits_diam<1.6], na.rm=T)
  sd_diam<-sd(recruits_diam[recruits_diam<1.6], na.rm=T)
  var_diam<-var(recruits_diam[recruits_diam<1.6], na.rm=T) 
  
  ###Summary statistics for just the 481 individuals with a height<16cm
  mu_height<-mean(recruits_height[recruits_height<16], na.rm=T) 
  sd_height<-sd(recruits_height[recruits_height<16],na.rm=T)
  var_height<-var(recruits_height[recruits_height<16],na.rm=T) 
  
  #Parameters associated with distribution of recruits diameter
  p.vec_BC[33]<-mu_diam
  p.vec_BC[34]<-sd_diam
  
  p.vec_CC[33]<-mu_diam
  p.vec_CC[34]<-sd_diam
  
  p.vec_C[33]<-mu_diam
  p.vec_C[34]<-sd_diam
  
  p.vec_FP[33]<-mu_diam
  p.vec_FP[34]<-sd_diam
  
  p.vec_PG[33]<-mu_diam
  p.vec_PG[34]<-sd_diam
  
  p.vec_WT[33]<-mu_diam
  p.vec_WT[34]<-sd_diam
  
  #Parameters associated with distribution of recruits heights
  
  p.vec_BC[35]<-mu_height
  p.vec_BC[36]<-sd_height
  
  p.vec_CC[35]<-mu_height
  p.vec_CC[36]<-sd_height
  
  p.vec_C[35]<-mu_height
  p.vec_C[36]<-sd_height
  
  p.vec_FP[35]<-mu_height
  p.vec_FP[36]<-sd_height
  
  p.vec_PG[35]<-mu_height
  p.vec_PG[36]<-sd_height
  
  p.vec_WT[35]<-mu_height
  p.vec_WT[36]<-sd_height
  
  
  ##### Some additional (fixed parameters)
  # TAU_1: Pre-dispersal seed survival 
  
  p.vec_BC[37]<-0.002
  p.vec_CC[37]<-0.002
  p.vec_C[37]<-0.002
  p.vec_FP[37]<-0.002
  p.vec_PG[37]<-0.002
  p.vec_WT[37]<-0.002
  
  
  #Delta  
  
  p.vec_BC[38]<-0.19
  p.vec_CC[38]<-0.19
  p.vec_C[38]<-0.19
  p.vec_FP[38]<-0.19
  p.vec_PG[38]<-0.19
  p.vec_WT[38]<-0.19
  
  #TAU_2: Post-dispersal seed survival
  
  #p.vec_overall[39]<- 0.5072 #Selected the highest survival probability (could have selected lowest) from the Geiger et al. study
  p.vec_BC[39]<-0.5072
  p.vec_CC[39]<-0.5072
  p.vec_C[39]<-0.5072
  p.vec_FP[39]<-0.5072
  p.vec_PG[39]<-0.5072
  p.vec_WT[39]<-0.5072
  
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
  pyx1=function(diamp,heightp,diam, height, params) { (1-((exp(params[12]+params[13]*diam + params[14]*height)/(1+exp(params[12]+params[13]*diam + params[14]*height)))))*
      (exp(params[1]+params[2]*diam+params[3]*height)/
         (1+exp(params[1]+params[2]*diam+params[3]*height)))*
      dtnorm(diamp,mean=params[4] + params[5]*diam + params[6]*height,sd=params[7], lower=0, upper=1.6)*dtnorm(heightp,mean=params[8]+params[9]*diam + params[10]*height,sd=params[11], lower=0, upper=16)}
  
  
  # Survival-Growth of Larger Plants (in D2)
  #Note that D2 survival is multiplied by 0.997 to avoid 100% survival 
  
  pyx2=function(diamp,heightp,diam, height, params) { 0.997*(exp(params[19]+params[20]*diam+params[21]*height)/
                                                               (1+exp(params[19]+params[20]*diam+params[21]*height)))*
      dtnorm(diamp,mean=params[22] + params[23]*diam + params[24]*height,sd=params[25], lower=1.6, upper=700)*dtnorm(heightp,mean=params[26]+params[27]*diam + params[28]*height,sd=params[29], lower=16, upper=800)}
  
  
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
  p.vecs <- list(p.vec_BC, p.vec_CC, p.vec_C, p.vec_FP, p.vec_PG, p.vec_WT)
  
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
    rm(D2, G, F, D2)
    
    
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
    
    
    thing <- find_lambda(A)
    lam.boot_random[boots,j] <- thing$lam.stable
    cat("Bootstrap:", boots,"Site:", j,"\n");
  }
}
end_time <- proc.time() - start_time
saveRDS(lam.boot_random, file = "lam.boot_random.rds")
#Time elapsed:
# user    system   elapsed 
# 734381.92  27859.35 769453.59