
#clear everything
rm(list=ls(all=TRUE))

require(Matrix);
library(MASS)
library(mvtnorm)
library(msm)

#increase memory limit on Windows machines  (This will make everything run faster)
memory.limit(memory.limit() *2^30)

m1=10
m2=m1+1
m3=100
m4=m3+1
tol=1.e-8; 

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

matrix.image=function(x,y,A,col=topo.colors(200),...) {
  nx=length(x); ny=length(y); 
  x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
  y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
  image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}

boot.lam <- function(dataset, dataset_allometry, sample.index) {

boot_data <- dataset[sample.index, ] 

# Subset the data into two size domains (seedlings and adults) #####
boot_data_seedlings <- subset(boot_data, boot_data$diameter_t < 1.6 & 
                                boot_data$height_t <16)
#Subset the seedling data correctly for modeling D1 growth: it should only include
# individuals who remain in the seedling domain at time t + 1
boot_data_seedlings_g <- subset(boot_data_seedlings, 
                                boot_data_seedlings$diameter_tplus1 < 1.6 & 
                                boot_data$height_tplus1 < 16)

boot_data_larges <- subset(boot_data, boot_data$diameter_t > 1.6 & 
                             boot_data$height_t >16)




# Initialize vectors to hold parameters #####
p.vec_BC<-rep(0, 39)
p.vec_C <- rep(0, 39)
p.vec_FP <- rep(0, 39)
p.vec_WT <- rep(0, 39)
p.vec_PG <- rep(0, 39)
p.vec_CC <- rep(0, 39)

#fit the functions and fill out the parameters: 
# Seedling Survival #####
s1<-glm(boot_data_seedlings$Surv_tplus1~ boot_data_seedlings$Diameter_t +
          boot_data_seedlings$Height_t+boot_data_seedlings$Site, family=binomial)

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

# Seedling Growth #####

g1_diam<-lm(boot_data_seedlings_g$Diameter_tplus1 ~ 
              boot_data_seedlings_g$Diameter_t+
              boot_data_seedlings_g$Height_t +
              boot_data_seedlings_g$Site)


#Then, model variance in growth for future height
g1_height<-lm(boot_data_seedlings_g$Height_tplus1 ~ 
                boot_data_seedlings_g$Diameter_t + 
                boot_data_seedlings_g$Height_t + 
                boot_data_seedlings_g$Site)


g1<-manova(cbind(boot_data_seedlings_g$Diameter_tplus1, 
                 boot_data_seedlings_g$Height_tplus1) ~ 
             boot_data_seedlings_g$Diameter_t+
             boot_data_seedlings_g$Height_t+
             boot_data_seedlings_g$Site)


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

# Maturation #####
graduates<-subset(boot_data_seedlings,
                  boot_data_seedlings$grad_status==1)


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
m <-glm(boot_data_seedlings$grad_status~
          boot_data_seedlings$Diameter_t + 
          boot_data_seedlings$Height_t + 
          boot_data_seedlings$Site, family=binomial)

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

# Adult survival #####

s2<-glm(boot_data_larges$Surv_tplus1~
          boot_data_larges$Diameter_t +
          boot_data_larges$Height_t +
          boot_data_larges$Site, family=binomial)

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

# Adult growth #####
#Model variance in growth in diameter: 

g2_diam<-lm(boot_data_larges$Diameter_tplus1 ~ 
          boot_data_larges$Diameter_t + 
            boot_data_larges$Height_t + boot_data_larges$Site)



#Model variance in growth in height: 

g2_height<-lm(boot_data_larges$Height_tplus1 ~ 
                boot_data_larges$Diameter_t + 
                boot_data_larges$Height_t +
                boot_data_larges$Site)




#(d) Modeling growth in D2 by site

g2<-manova(cbind(boot_data_larges$Diameter_tplus1,
                 boot_data_larges$Height_tplus1) ~ 
             boot_data_larges$Diameter_t+
             boot_data_larges$Height_t+
             boot_data_larges$Site)


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




#Model variance in growth of D2 individuals in diameter


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

# Reproduction #####
p_f<-glm(boot_data_larges$Rep_tplus1 ~ 
           boot_data_larges$Height_t +
           boot_data_larges$Site, family=binomial)

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


# Fecundity #####

x<-dataset_allometry$diam_base*10 #convert from cm to mm
x <- x*x
f<-lm(dataset_allometry$Seed_No ~ 0 + x:dataset_allometry$Site)

#(b) Model fecundity (by biotype)
p.vec_BC[32]<-f$coeff[1]
p.vec_CC[32]<-f$coeff[2]
p.vec_C[32]<-f$coeff[3]
p.vec_FP[32] <- f$coeff[4]
p.vec_PG[32] <- f$coeff[5]
p.vec_WT[32] <- f$coeff[6]

# Distribution of new recruits size #####

recruit4_diam<-demog$Diam_4[demog$Status_4=="tagged"&demog$Location=="SeedlingPlot"]
recruit5_diam<-demog$Diam_5[demog$Status_5=="tagged"&demog$Location=="SeedlingPlot"]
recruit6_diam<-demog$Diam_6[demog$Status_6=="tagged"&demog$Location=="SeedlingPlot"]
recruits_diam<-c(recruit4_diam, recruit5_diam, recruit6_diam)

recruit4_height<-demog$Height_4[demog$Status_4=="tagged"&demog$Location=="SeedlingPlot"]
recruit5_height<-demog$Height_5[demog$Status_5=="tagged"&demog$Location=="SeedlingPlot"]
recruit6_height<-demog$Height_6[demog$Status_6=="tagged"&demog$Location=="SeedlingPlot"]
recruits_height<-c(recruit4_height, recruit5_height, recruit6_height)

mu_diam<-mean(recruits_diam[recruits_diam<1.6], na.rm=T)
sd_diam<-sd(recruits_diam[recruits_diam<1.6], na.rm=T)
var_diam<-var(recruits_diam[recruits_diam<1.6], na.rm=T) 

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

# Additional fixed parameters ##### 
# TAU_1: Pre-dispersal seed survival 
#Rethinking this parameter value: Isn't this already encapsulated in tau_2? 

p.vec_BC[37]<-0.002
p.vec_CC[37]<-0.002
p.vec_C[37]<-0.002
p.vec_FP[37]<-0.002
p.vec_PG[37]<-0.002
p.vec_WT[37]<-0.002


p.vec_BC[38]<-0.19
p.vec_CC[38]<-0.19
p.vec_C[38]<-0.19
p.vec_FP[38]<-0.19
p.vec_PG[38]<-0.19
p.vec_WT[38]<-0.19

p.vec_BC[39]<-0.5072
p.vec_CC[39]<-0.5072
p.vec_C[39]<-0.5072
p.vec_FP[39]<-0.5072
p.vec_PG[39]<-0.5072
p.vec_WT[39]<-0.5072


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

}
