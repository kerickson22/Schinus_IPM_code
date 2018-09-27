#####Code for performing demographic analysis for overall model as well as by biotype
#####This code initializes a vector (called p.vec) that will be used in constructed the integral projection models
### By Kelley D. Erickson and Carol C. Horvitz


#Load demographic data
demog<-read.csv("demography_15_clean.csv", head=T)

##This is a csv file where each row is a an individual 
# The column names are : 
# Obs : Observation number (each plant should have a different observation number; these were assigned post data 
#      entry when we discovered that some of the numbers for individuals were reused between sites
# Site: Which site (Big Cypress, Cape Canaveral, Chekika, Forty Pierce, Punta Gorda, Wild Turkey)
# Plant ID: The tag number of an individual (each individual in a site has a different number, but numbers were
#      reused across sites)
# Location: SeedlingPlot (individual located in a seedling plot) or NA (individual not in seedling plot)
# Date_i: The date the census took place for censuses i = 1:7
# Status_i: One of four categories for censuses i=1:7: 
#     alive: individual was alive when censused
#     firstdead: the individual was recorded as dead at this census and was alive at the last census 
#     missing: the individual was not located at this census
#     prestudy: the individual had not been tagged yet 
#     stilldead: the individual was already marked dead at a previous census
# Death_Cause_i: One of three categories for censuses i = 1:7
#     NA: for individuals that died, there was no reason to suspect their deaths were due to management actions
#     sprayed: individual death attributed to being sprayed with herbicide (these individuals should be excluded)
#     fire: individual death attributed to fire (these individuals should be excluded from analysis)
# Rep_stat_i: One of 4 categories for censuses  i = 1:7
#     NA: for individuals whose reproductive status is unknown at census i 
#     female: for individuals who have produced fruit at any previous census up to an including census i 
#     male: for individuals who have not produced fruit at any previous census (may have observed male flowers)
#     pre-reproductive: for individuals who have not produced fruit at any previous census 
# Diam_i: Diameter (in mm) of individual at census i 
# Height_i: Height (in cm) of individual at census i 
# Max_i: Maximum canopy width in cm of individual at census i 
# Min_i: Minimum canopy width in cm of individual at census i 
# Multistem: For each individual takes on one of two values: 
#     NA: individuals that are not multistem
#     FLAG: Individuals that are multistem 



##### Initialize pvec's to store parameters that will go into the IPMs:
#Initialize p.vec's to store output from model (these will feed into the IPM)
p.vec_E<-rep(0, 39);
p.vec_H<-rep(0, 39);
p.vec_W<-rep(0, 39);
p.vec_overall<-rep(0, 39)

#### Restructure the data in a few different ways: 

Death_Cause<-c(demog$Death_Cause_2, demog$Death_Cause_3, demog$Death_Cause_4, demog$Death_Cause_5, demog$Death_Cause_6, demog$Death_Cause_7)

Status<-c(demog$Status_1, demog$Status_2, demog$Status_3, demog$Status_4, demog$Status_5, demog$Status_6, demog$Status_7)

Diams<-cbind(demog$Diam_1, demog$Diam_2, demog$Diam_3, demog$Diam_4, demog$Diam_5, demog$Diam_6, demog$Diam_7) 

row.names(Diams)<-demog$Plant_ID

Heights<-cbind(demog$Height_1, demog$Height_2, demog$Height_3, demog$Height_4, demog$Height_5, demog$Height_6, demog$Height_7)

row.names(Heights)<-demog$Plant_ID

Mins<-cbind(demog$Min_1, demog$Min_2, demog$Min_3, demog$Min_4, demog$Min_5, demog$Min_6, demog$Min_7)

row.names(Mins)<-demog$Plant_ID

Maxs<-cbind(demog$Max_1, demog$Max_2, demog$Max_3, demog$Max_4, demog$Max_5, demog$Max_6, demog$Max_7)

row.names(Maxs)<-demog$Plant_ID

Diams_pooled<-rbind( cbind(Diams[,1], Diams[,2]), cbind(Diams[,2], Diams[,3]), cbind(Diams[, 3], Diams[,4]), cbind(Diams[,4], Diams[,5]), cbind(Diams[,5], Diams[,6]), cbind(Diams[, 6], Diams[, 7]))

Diams_pooled<-data.frame(Diams_pooled[,1], Diams_pooled[,2])

names(Diams_pooled)<-c("t", "t_plus_one")

Sites_pooled<-rep(demog$Site, 6)

ID_pooled<-rep(demog$Plant_ID, 6)

Location<-rep(demog$Location, 6)

Genetic_type_pooled<-rep(0, length(Sites_pooled))

#code which sites are which biotype
for (i in 1:length(Genetic_type_pooled)) { 
	if (Sites_pooled[i]=="Big Cypress") Genetic_type_pooled[i]<-"hybrid"
	else if (Sites_pooled[i]=="Cape Canaveral") Genetic_type_pooled[i]<-"hybrid"
	else if (Sites_pooled[i]=="Chekika") Genetic_type_pooled[i]<-"eastern"
	else if (Sites_pooled[i]=="Fort Pierce") Genetic_type_pooled[i]<-"eastern"
	else if (Sites_pooled[i]=="Punta Gorda") Genetic_type_pooled[i]<-"western"
	else Genetic_type_pooled[i]<-"western"
	
}



Heights_pooled<-rbind( cbind(Heights[,1], Heights[,2]), cbind(Heights[,2], Heights[,3]), cbind(Heights[, 3], Heights[,4]), cbind(Heights[,4], Heights[,5]), cbind(Heights[,5], Heights[,6]), cbind(Heights[, 6], Heights[, 7]))

Heights_pooled<-data.frame(Heights_pooled[,1], Heights_pooled[,2])
names(Heights_pooled)<-c("t", "t_plus_one")



#Finish for maxs and mins 

Maxs_pooled<-rbind(cbind(Maxs[,1], Maxs[,2]), cbind(Maxs[,2], Maxs[,3]), 
	cbind(Maxs[,3], Maxs[,4]), cbind(Maxs[,4], Maxs[,5]), 
	cbind(Maxs[,5], Maxs[,6]), cbind(Maxs[,6], Maxs[,7]))
	
Maxs_pooled<-data.frame(Maxs_pooled[,1], Maxs_pooled[,2])
names(Maxs_pooled)<-c("t", "t_plus_one")



Mins_pooled<-rbind( cbind(Mins[,1], Mins[,2]), cbind(Mins[,2], Mins[,3]), 
	cbind(Mins[,3], Mins[,4]), cbind(Mins[,4], Mins[,5]), 
	cbind(Mins[, 5], Mins[,6]), cbind(Mins[,6], Mins[,7]))
	
Mins_pooled<-data.frame(Mins_pooled[,1], Mins_pooled[,2])
names(Mins_pooled)<-c("t", "t_plus_one")




###Use status matrix to get survival matrix
statmat<-as.matrix(demog[, c(6, 14, 22, 30, 38, 46, 54)])
survmat<-statmat
survmat<- ifelse(survmat=="prestudy"|survmat=="missing"|survmat=="still_dead", NA, survmat)
survmat<- ifelse(survmat =="tagged"|survmat=="alive", 1, 0)

survmat_pooled<-rbind( cbind(survmat[,1], survmat[,2]), cbind(survmat[,2], survmat[,3]), cbind(survmat[, 3], survmat[,4]), cbind(survmat[,4], survmat[,5]), cbind(survmat[,5], survmat[,6]), cbind(survmat[, 6], survmat[, 7]))
names(survmat_pooled)<-c("surv_t", "surv_t_plus_one")


survmat_pooled_bysize<-cbind(Diams_pooled[, 1], survmat_pooled[,2])
names(survmat_pooled_bysize)<-c("size_t", "surv_t_plus_one")



#Use reproduction status matrix to get reproduction matrix
repstatmat<-as.matrix(demog[, c(8, 16, 24, 32, 40, 48, 56)])
repmat<-repstatmat
repmat<- ifelse(repmat =="female", 1, 0)
repmat_pooled<-rbind( cbind(repmat[,1], repmat[,2]), cbind(repmat[,2], repmat[,3]), cbind(repmat[, 3], repmat[,4]), cbind(repmat[,4], repmat[,5]), cbind(repmat[,5], repmat[,6]), cbind(repmat[, 6], repmat[, 7]))


pooled<-cbind(ID_pooled, Sites_pooled, Location, Death_Cause, Genetic_type_pooled, Mins_pooled, Maxs_pooled, Diams_pooled, Heights_pooled, repmat_pooled, survmat_pooled)
names(pooled)<-c("ID", "Site","Location", "Death_Cause", "Genetic_type", "Min_t", "Min_tplus1", "Max_t", "Max_tplus1", "Diameter_t", "Diameter_tplus1", "Height_t", "Height_tplus1", "Rep_t", "Rep_tplus1", "Surv_t", "Surv_tplus1")


##### Figuring out which state variables to use: 
predictors<-cbind(pooled$Min_t, pooled$Diameter_t, pooled$Height_t, pooled$Max_t)

# What is the correlation amongst the predictor variables? 
correlation<-cor(predictors[, 1:4], use="pairwise.complete.obs")
colnames(correlation)<-c("Min", "Diameter", "Height", "Max")
rownames(correlation)<-c("Min", "Diameter", "Height", "Max")

correlation<-round(correlation, 2)
correlations<-write.csv(correlation, "correlations.csv")
#Over the entire dataset, there is the lowest correlation between diameter and height



#####Cleaning up the data

#Remove individual deaths that were sprayed or fire or grinding
restrict<-subset(pooled, is.na(Death_Cause))

#Remove outliers (very large diameters)
restrict2<-subset(restrict, Diameter_t<800)

#Range of diameter: [0, 800]
0.9*max(restrict2$Diameter_t) #700
#Range of height:  [0, 800]
0.9*max(restrict2$Height_t, na.rm=T) #800



#Divide dataset into size domains (seedlings and larges)
seedlings<-subset(restrict2, restrict2$Diameter_t<1.6)
seedlings<-subset(seedlings, seedlings$Height_t<16)
larges<-subset(restrict2, restrict2$Diameter_t>1.6)
larges<-subset(larges, larges$Height_t>16)



x1seq_larges<-seq(1.6, 800, length.out=50)
x2seq_larges<-seq(16, 800, length.out=50)

x1seq_seedlings<-seq(0, 1.6, length.out=50)
x2seq_seedlings<-seq(0, 16, length.out=50)


##### SURVIVAL of Seedlings (D1)

#(a) Overall model 
surv_mod_seedlings<-glm(seedlings$Surv_tplus1~ seedlings$Diameter_t + seedlings$Height_t, family=binomial)
summary(surv_mod_seedlings)
#All terms are significant

p.vec_overall[1]<-surv_mod_seedlings$coeff[1]
p.vec_overall[2]<-surv_mod_seedlings$coeff[2]
p.vec_overall[3]<-surv_mod_seedlings$coeff[3]

#(b) By Biotype 

surv_mod_seedlings_site<-glm(seedlings$Surv_tplus1~ seedlings$Diameter_t + seedlings$Height_t+seedlings$Site, family=binomial)
summary(surv_mod_seedlings_site)
#Because I have low sample size at some sites, I decided to group by biotype instead of focusing on differences between sites.

surv_mod_seedlings_lineage<-glm(seedlings$Surv_tplus1 ~ seedlings$Diameter_t + seedlings$Height_t + seedlings$Genetic_type, family=binomial)
summary(surv_mod_seedlings_lineage)
#Hybrid is different from E+W, but Eastern and Western are not different from eachother

#Create a new dataframe where Eastern and Western biotypes are lumped together 
lineage2<-rep(0, length(seedlings$Genetic_type))
for (i in 1:length(seedlings$Genetic_type)) {
  if (seedlings$Genetic_type[i]=="western" || seedlings$Genetic_type[i]=="eastern") {lineage2[i]<-"E-W"}
  else{lineage2[i]<-"Hybrid"}
}

seedlings2<-cbind(seedlings, lineage2)

surv_mod_seedlings_lineage2<-glm(seedlings2$Surv_tplus1 ~ seedlings2$Diameter_t + seedlings2$Height_t + seedlings2$lineage2, family=binomial)
summary(surv_mod_seedlings_lineage2)



p.vec_E[1]<-surv_mod_seedlings_lineage2$coeff[1]
p.vec_E[2]<-surv_mod_seedlings_lineage2$coeff[2]
p.vec_E[3]<-surv_mod_seedlings_lineage2$coeff[3]

p.vec_H[1]<-surv_mod_seedlings_lineage2$coeff[1]+surv_mod_seedlings_lineage2$coeff[4]
p.vec_H[2]<-surv_mod_seedlings_lineage2$coeff[2]
p.vec_H[3]<-surv_mod_seedlings_lineage2$coeff[3]

p.vec_W[1]<-surv_mod_seedlings_lineage2$coeff[1]
p.vec_W[2]<-surv_mod_seedlings_lineage2$coeff[2]
p.vec_W[3]<-surv_mod_seedlings_lineage2$coeff[3]



##### Seedling(D1) GROWTH 

#Size at t+1 is a function of height_t AND diam_t

#(a) Overall model 
growth_mod_seedlings<-manova(cbind(seedlings$Diameter_tplus1, seedlings$Height_tplus1) ~ seedlings$Diameter_t+seedlings$Height_t)
summary(growth_mod_seedlings)

#Grab parameters associated with future diameter
p.vec_overall[4]<-growth_mod_seedlings$coeff[1,1]
p.vec_overall[5]<-growth_mod_seedlings$coeff[2,1]
p.vec_overall[6]<-growth_mod_seedlings$coeff[3,1]

#To model variance in D1 growth for diameter and height, use separate linear models and store the variance: 
#First, model variance in growth for future diameter
growth_mod_diam_seedlings<-lm(seedlings$Diameter_tplus1 ~ seedlings$Diameter_t+seedlings$Height_t)
summary(growth_mod_diam_seedlings)
p.vec_overall[7]<-summary(growth_mod_diam_seedlings)$sigma

#Grab parameters associated with future height
p.vec_overall[8]<-growth_mod_seedlings$coeff[1,2]
p.vec_overall[9]<-growth_mod_seedlings$coeff[2,2]
p.vec_overall[10]<-growth_mod_seedlings$coeff[3,2]


#Then, model variance in growth for future height
growth_mod_height_seedlings<-lm(seedlings$Height_tplus1 ~ seedlings$Diameter_t + seedlings$Height_t)
summary(growth_mod_height_seedlings)
p.vec_overall[11]<-summary(growth_mod_height_seedlings)$sigma


#(b) By Biotype: 
growth_mod_seedlings_lineage<-manova(cbind(seedlings$Diameter_tplus1, seedlings$Height_tplus1) ~ seedlings$Diameter_t+seedlings$Height_t+seedlings$Genetic_type)
summary(growth_mod_seedlings_lineage)
#Lineage does not appear to have an effect on seedling growth, so use overall model for all biotypes

#Capture parameters associated with predicting future diameter 
p.vec_E[4]<-growth_mod_seedlings$coeff[1,1]
p.vec_E[5]<-growth_mod_seedlings$coeff[2,1]
p.vec_E[6]<-growth_mod_seedlings$coeff[3,1]
p.vec_E[7]<-summary(growth_mod_diam_seedlings)$sigma

p.vec_H[4]<-growth_mod_seedlings$coeff[1,1]
p.vec_H[5]<-growth_mod_seedlings$coeff[2,1]
p.vec_H[6]<-growth_mod_seedlings$coeff[3,1]
p.vec_H[7]<-summary(growth_mod_diam_seedlings)$sigma

p.vec_W[4]<-growth_mod_seedlings$coeff[1,1]
p.vec_W[5]<-growth_mod_seedlings$coeff[2,1]
p.vec_W[6]<-growth_mod_seedlings$coeff[3,1]
p.vec_W[7]<-summary(growth_mod_diam_seedlings)$sigma

#Grab parameters associated with predicting future height
p.vec_E[8]<-growth_mod_seedlings$coeff[1,2]
p.vec_E[9]<-growth_mod_seedlings$coeff[2,2]
p.vec_E[10]<-growth_mod_seedlings$coeff[3,2]
p.vec_E[11]<-summary(growth_mod_height_seedlings)$sigma

p.vec_H[8]<-growth_mod_seedlings$coeff[1,2]
p.vec_H[9]<-growth_mod_seedlings$coeff[2,2]
p.vec_H[10]<-growth_mod_seedlings$coeff[3,2]
p.vec_H[11]<-summary(growth_mod_height_seedlings)$sigma

p.vec_W[8]<-growth_mod_seedlings$coeff[1,2]
p.vec_W[9]<-growth_mod_seedlings$coeff[2,2]
p.vec_W[10]<-growth_mod_seedlings$coeff[3,2]
p.vec_W[11]<-summary(growth_mod_height_seedlings)$sigma

##### "GRADUATION" from the seedling(D1) domain into the larger plant domain (D2)



#Create a new variable that tracks whether individuals transition out of D1 into D2
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

#Probability of leaving the seedling domain: 


#(a) Overall model

sdlng_grad_mod<-glm(seedlings$grad_status~ seedlings$Diameter_t + seedlings$Height_t, family=binomial)
summary(sdlng_grad_mod)

p.vec_overall[12]<-sdlng_grad_mod$coeff[1]
p.vec_overall[13]<-sdlng_grad_mod$coeff[2]
p.vec_overall[14]<-sdlng_grad_mod$coeff[3]


#Model distribution of sizes of graduates 

graduates<-subset(seedlings, seedlings$grad_status==1)


#Get distribution of new recruits sizes (these are not modeled by biotype) 
summary(graduates$Diameter_tplus1)
hist(graduates$Diameter_tplus1)
mu_grad_diam<-mean(graduates$Diameter_tplus1)
sd_grad_diam<-sd(graduates$Diameter_tplus1)

summary(graduates$Height_tplus1)
hist(graduates$Height_tplus1)
mu_grad_height<-mean(graduates$Height_tplus1)
sd_grad_height<-sd(graduates$Height_tplus1)


p.vec_overall[15]<-mu_grad_diam
p.vec_overall[16]<-sd_grad_diam

p.vec_overall[17]<-mu_grad_height
p.vec_overall[18]<-sd_grad_height

#(b) Model graduation by biotype
sdlng_grad_mod_lineage<-glm(seedlings$grad_status~ seedlings$Diameter_t + seedlings$Height_t + seedlings$Genetic_type, family=binomial)
summary(sdlng_grad_mod_lineage)
#The effect of hybrids is only moderatately significantly diff from E, so collapse E and H

lineage4<-rep(0, length(seedlings$Genetic_type))
for (i in 1:length(seedlings$Genetic_type)) {
  if (seedlings$Genetic_type[i]=="eastern" || seedlings$Genetic_type[i]=="hybrid") {lineage4[i]<-"E-H"}
  else{lineage4[i]<-"Western"}
}

seedlings3<-cbind(seedlings, lineage4)

grad_mod_lineage2<-glm(seedlings3$grad_status ~ seedlings3$Diameter_t + seedlings3$Height_t + seedlings3$lineage4, family=binomial)
summary(grad_mod_lineage2)


p.vec_E[12]<-grad_mod_lineage2$coeff[1]
p.vec_E[13]<-grad_mod_lineage2$coeff[2]
p.vec_E[14]<-grad_mod_lineage2$coeff[3]

p.vec_H[12]<-grad_mod_lineage2$coeff[1]
p.vec_H[13]<-grad_mod_lineage2$coeff[2]
p.vec_H[14]<-grad_mod_lineage2$coeff[3]

p.vec_W[12]<-grad_mod_lineage2$coeff[1]+grad_mod_lineage2$coeff[4]
p.vec_W[13]<-grad_mod_lineage2$coeff[2]
p.vec_W[14]<-grad_mod_lineage2$coeff[3]

#The distribution of graduates sizes does not differ by biotype: 
p.vec_E[15]<-mu_grad_diam
p.vec_E[16]<-sd_grad_diam

p.vec_H[15]<-mu_grad_diam
p.vec_H[16]<-sd_grad_diam

p.vec_W[15]<-mu_grad_diam
p.vec_W[16]<-sd_grad_diam

p.vec_E[17]<-mu_grad_height
p.vec_E[18]<-sd_grad_height

p.vec_H[17]<-mu_grad_height
p.vec_H[18]<-sd_grad_height

p.vec_W[17]<-mu_grad_height
p.vec_W[18]<-sd_grad_height

##### SURVIVAL of larger plants (in the D2 domain)
#(a) Overall model 

surv_mod_larges<-glm(larges$Surv_tplus1~larges$Diameter_t + larges$Height_t, family=binomial)
summary(surv_mod_larges)

p.vec_overall[19]<-surv_mod_larges$coeff[1]
p.vec_overall[20]<-surv_mod_larges$coeff[2]
p.vec_overall[21]<-surv_mod_larges$coeff[3]

#(b) Model D2-survival (by biotype)

surv_mod_larges_lineage<-glm(larges$Surv_tplus1~larges$Diameter_t + larges$Height_t +larges$Genetic_type, family=binomial)
summary(surv_mod_larges_lineage)
#Hybrid is different from E+W, but Eastern and Western are not different from eachother

lineage3<-rep(0, length(larges$Genetic_type))
for (i in 1:length(larges$Genetic_type)) {
  if (larges$Genetic_type[i]=="western" || larges$Genetic_type[i]=="eastern") {lineage3[i]<-"E-W"}
  else{lineage3[i]<-"Hybrid"}
}

larges2<-cbind(larges, lineage3)

surv_mod_larges_lineage2<-glm(larges2$Surv_tplus1 ~ larges2$Diameter_t + larges2$Height_t+larges2$lineage3, family=binomial)
summary(surv_mod_larges_lineage2)

p.vec_E[19]<-surv_mod_larges_lineage2$coeff[1]
p.vec_E[20]<-surv_mod_larges_lineage2$coeff[2]
p.vec_E[21]<-surv_mod_larges_lineage2$coeff[3]

p.vec_H[19]<-surv_mod_larges_lineage2$coeff[1] + surv_mod_larges_lineage2$coeff[4]
p.vec_H[20]<-surv_mod_larges_lineage2$coeff[2]
p.vec_H[21]<-surv_mod_larges_lineage2$coeff[3]

p.vec_W[19]<-surv_mod_larges_lineage2$coeff[1]
p.vec_W[20]<-surv_mod_larges_lineage2$coeff[2]
p.vec_W[21]<-surv_mod_larges_lineage2$coeff[3]

##### GROWTH of larger plants (in the D2 domain)

#(a) Overall model 
growth_mod_larges<-manova(cbind(larges$Diameter_tplus1, larges$Height_tplus1) ~ larges$Diameter_t+larges$Height_t)
summary(growth_mod_larges)

#Store parameters associated with modeling future growth in diameter
p.vec_overall[22]<-growth_mod_larges$coeff[1,1]
p.vec_overall[23]<-growth_mod_larges$coeff[2,1]
p.vec_overall[24]<-growth_mod_larges$coeff[3,1]

#Model variance in growth in diameter: 

growth_mod_diam_larges<-lm(larges$Diameter_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Genetic_type)
summary(growth_mod_diam_larges)

p.vec_overall[25]<-summary(growth_mod_diam_larges)$sigma

#Store parameters associated with modeling future growth in height
p.vec_overall[26]<-growth_mod_larges$coeff[1,2]
p.vec_overall[27]<-growth_mod_larges$coeff[2,2]
p.vec_overall[28]<-growth_mod_larges$coeff[3,2]

#Model variance in growth in height: 

growth_mod_height_larges<-lm(larges$Height_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Genetic_type)
summary(growth_mod_height_larges)
p.vec_overall[29]<-summary(growth_mod_height_larges)$sigma


#(d) Modeling growth in D2 by biotype

growth_mod_larges_lineage<-manova(cbind(larges$Diameter_tplus1, larges$Height_tplus1) ~ larges$Diameter_t+larges$Height_t+larges$Genetic_type)
summary(growth_mod_larges_lineage)
#BUT lineage DOES have an effect on the growth of the larger individuals in D2

#Store parameters associated with growth of D2 individuals in diameter
p.vec_E[22]<-growth_mod_larges_lineage$coeff[1,1]
p.vec_E[23]<-growth_mod_larges_lineage$coeff[2,1]
p.vec_E[24]<-growth_mod_larges_lineage$coeff[3,1]

p.vec_H[22]<-growth_mod_larges_lineage$coeff[1,1] + growth_mod_larges_lineage$coeff[4,1]
p.vec_H[23]<-growth_mod_larges_lineage$coeff[2,1]
p.vec_H[24]<-growth_mod_larges_lineage$coeff[3,1]

p.vec_W[22]<-growth_mod_larges_lineage$coeff[1,1] + growth_mod_larges_lineage$coeff[5, 1]
p.vec_W[23]<-growth_mod_larges_lineage$coeff[2,1]
p.vec_W[24]<-growth_mod_larges_lineage$coeff[3,1]


#Model variance in growth of D2 individuals in diameter
growth_mod_diam_larges_lineage<-lm(larges$Diameter_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Genetic_type)
summary(growth_mod_diam_larges_lineage)

p.vec_E[25]<-summary(growth_mod_diam_larges_lineage)$sigma
p.vec_H[25]<-summary(growth_mod_diam_larges_lineage)$sigma
p.vec_W[25]<-summary(growth_mod_diam_larges_lineage)$sigma

#Store parameters associated with modeling growth of D2 individuals in height
p.vec_E[26]<-growth_mod_larges_lineage$coeff[1,2]
p.vec_E[27]<-growth_mod_larges_lineage$coeff[2,2]
p.vec_E[28]<-growth_mod_larges_lineage$coeff[3,2]

p.vec_H[26]<-growth_mod_larges_lineage$coeff[1,2] + growth_mod_larges_lineage$coeff[4,2]
p.vec_H[27]<-growth_mod_larges_lineage$coeff[2,2]
p.vec_H[28]<-growth_mod_larges_lineage$coeff[3,2]

p.vec_W[26]<-growth_mod_larges_lineage$coeff[1,2] + growth_mod_larges_lineage$coeff[5,2]
p.vec_W[27]<-growth_mod_larges_lineage$coeff[2,2]
p.vec_W[28]<-growth_mod_larges_lineage$coeff[3,2]

#Model variance in growth of D2 individuals in height
growth_mod_height_larges_lineage<-lm(larges$Height_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Genetic_type)
summary(growth_mod_height_larges_lineage)

p.vec_E[29]<-summary(growth_mod_height_larges_lineage)$sigma
p.vec_H[29]<-summary(growth_mod_height_larges_lineage)$sigma
p.vec_W[29]<-summary(growth_mod_height_larges_lineage)$sigma

##### Reproduction

#(a) Model probability of being reproductive (overall)
mod_repro1<-glm(larges$Rep_tplus1 ~ larges$Diameter_t+larges$Height_t, family=binomial)
summary(mod_repro)
#Diameter term is not significant, so remove it: 
mod_repro<-glm(larges$Rep_tplus1 ~ larges$Height_t, family=binomial)

p.vec_overall[30]<-mod_repro$coeff[1]
p.vec_overall[31]<-mod_repro$coeff[2]

#(b) Model probability of being reproductive     (by biotype)
mod_repro_lineage1<-glm(larges$Rep_tplus1 ~ larges$Diameter_t+larges$Height_t + larges$Genetic_type, family=binomial)
summary(mod_repro_lineage1)
#Diameter term is not significant, so remove it: 
mod_repro_lineage<-glm(larges$Rep_tplus1 ~ larges$Height_t +larges$Genetic_type, family=binomial)
#Hybrid and Eastern are not different from eachother so collapse them: 

lineage4<-rep(0, length(larges$Genetic_type))
for (i in 1:length(larges$Genetic_type)) {
  if (larges$Genetic_type[i]=="hybrid" || larges$Genetic_type[i]=="eastern") {lineage3[i]<-"E-H"}
  else{lineage4[i]<-"Western"}
}

larges3<-cbind(larges, lineage4)

mod_repro_lineage2<-glm(larges3$Rep_tplus1~ larges3$Height_t + larges3$lineage4, family=binomial)
summary(mod_repro_lineage2)

#Eastern and Hybrid have same probabilty of reproducing, but W differs

p.vec_E[30]<-mod_repro_lineage2$coeff[1]
p.vec_E[31]<-mod_repro_lineage2$coeff[2]

p.vec_H[30]<-mod_repro_lineage2$coeff[1]
p.vec_H[31]<-mod_repro_lineage2$coeff[2]

p.vec_W[30]<-mod_repro_lineage2$coeff[1]+mod_repro_lineage2$coeff[3]
p.vec_W[31]<-mod_repro_lineage2$coeff[2]


##### Fecundity (given that an individual is reproductive )
#These values are taken from Erickson et al. 2017 and rescaled (in that paper diameter was 
# measured in cm instead of mm)
#(a) Overall 

p.vec_overall[32]<-4.89

#(b) Model fecundity (by biotype)
p.vec_E[32]<-4
p.vec_H[32]<-4.25
p.vec_W[32]<-5.52


##### Size distribution of new recruits 
#Figure out distribution of newly tagged individuals in seedling plots 
#Construct a dataframe just holding recruits (individuals that go from 'untagged' to 'tagged' 
#AND are located within a seedlingplot)

recruit4_diam<-demog$Diam_4[demog$Status_4=="tagged"&demog$Location=="SeedlingPlot"]
recruit5_diam<-demog$Diam_5[demog$Status_5=="tagged"&demog$Location=="SeedlingPlot"]
recruit6_diam<-demog$Diam_6[demog$Status_6=="tagged"&demog$Location=="SeedlingPlot"]
recruits_diam<-c(recruit4_diam, recruit5_diam, recruit6_diam)

recruit4_height<-demog$Height_4[demog$Status_4=="tagged"&demog$Location=="SeedlingPlot"]
recruit5_height<-demog$Height_5[demog$Status_5=="tagged"&demog$Location=="SeedlingPlot"]
recruit6_height<-demog$Height_6[demog$Status_6=="tagged"&demog$Location=="SeedlingPlot"]
recruits_height<-c(recruit4_height, recruit5_height, recruit6_height)


#Not confident that ALL of these are actually new recruits (some are way too tall), so: 


#How does R calculate the whiskers?:
#https://www.r-bloggers.com/about-boxplot/

#min (max(recruits_height, na.rm=T), 10+1.5*(10-6))

#upper whisker=min(max(x), Q_3 + 1.5 * IQR)
#lower whisker = max(min(x), Q_1 – 1.5 * IQR) 
#IQR=Q3-Q1

#From ?boxplot help: 

#stat     # a matrix, each column contains the extreme of the lower whisker, the lower hinge, the median, the upper hinge and the extreme of the upper whisker for one group/plot. If all the inputs have the same class attribute, so will this component.

#The upper whisker extends up to 1.6mm


#Solution: Restrict to only the part of the data that falls below the upper hinge. 
#For diameter, this is all individuals <1.6mm
#For height, this is all individuals <16cm




#####Summary statistics for just the 476 individuals with a diameter<1.6mm
boxplot(recruits_diam[recruits_diam<1.6])
hist(recruits_diam[recruits_diam<1.6])
summary(recruits_diam[recruits_diam<1.6])
mu_diam<-mean(recruits_diam[recruits_diam<1.6], na.rm=T)
sd_diam<-sd(recruits_diam[recruits_diam<1.6], na.rm=T)
var_diam<-var(recruits_diam[recruits_diam<1.6], na.rm=T) 
length(recruits_diam[recruits_diam<1.6])




###Summary statistics for just the 481 individuals with a height<16cm
summary(recruits_height[recruits_height<16])
boxplot(recruits_height[recruits_height<16])
thing3<-hist(recruits_height[recruits_height<16])
mu_height<-mean(recruits_height[recruits_height<16], na.rm=T) 
sd_height<-sd(recruits_height[recruits_height<16],na.rm=T)
var_height<-var(recruits_height[recruits_height<16],na.rm=T) 
length(recruits_height[recruits_height<16])


recruits<-cbind(recruits_diam, recruits_height)
recruits<-as.data.frame(recruits)
names(recruits)<-c("diam", "height")
recruits_1<-subset(recruits, recruits$diam<1.6)
recruits_2<-subset(recruits_1, recruits_1$height<16)


#Parameters associated with distribution of recruits diameter
p.vec_overall[33]<-mu_diam
p.vec_overall[34]<-sd_diam

p.vec_E[33]<-mu_diam
p.vec_E[34]<-sd_diam

p.vec_H[33]<-mu_diam
p.vec_H[34]<-sd_diam

p.vec_W[33]<-mu_diam
p.vec_W[34]<-sd_diam


#Parameters associated with distribution of recruits heights
p.vec_overall[35]<-mu_height
p.vec_overall[36]<-sd_height

p.vec_E[35]<-mu_height
p.vec_E[36]<-sd_height

p.vec_H[35]<-mu_height
p.vec_H[36]<-sd_height

p.vec_W[35]<-mu_height
p.vec_W[36]<-sd_height


##### Some additional (fixed parameters)
# TAU_1: Pre-dispersal seed survival 
#Rethinking this parameter value: Isn't this already encapsulated in tau_2? 
#p.vec_overall[37]<-0.002
#p.vec_E[37]<-0.002
#p.vec_H[37]<-0.002
#p.vec_W[37]<-0.002

#At some point could write this out of the code, but for now will eliminate it by setting it to 1
p.vec_overall[37]<- 1
p.vec_E[37]<-1
p.vec_H[37]<-1
p.vec_W[37]<-1

#DELTA: Probability of dispersal 
# p.vec_overall[38]<-0.005
# p.vec_E[38]<-0.005
# p.vec_H[38]<-0.005
# p.vec_W[38]<-0.005

#Revised based on calculations for appendix 
p.vec_overall[38]<-0.19
p.vec_E[38]<-0.19
p.vec_H[38]<-0.19
p.vec_W[38]<-0.19

#TAU_2: Post-dispersal seed survival
# Estimated from Geiger et al. 2011 

# Geiger et al. (2011) report 8mo survival rate from seeds sourced from the six sites grown in a 
# common garden experiment in Davie: 
# "A greater proportion of hybrid seedlings survived than did western seedlings
# (63.6% vs. 53.2%; see fig. 5)." (Geiger et al. 2011) 



month<-seq(0, 12, by=1)
hybrid_seed_surv<-rep(100, 12)
eastwest_seed_surv<-rep(100, 12)


#Assumption: assume mortality is exponential: 
#For hybrids:
#	e^(-8*lambda)=0.636
#	-8*lambda = ln(63.6)
#	lambda=(ln(63.6))/-8
lambda_hybrid<-(log(0.636))/(-8)
lambda_eastwest<-(log(0.532))/(-8)


y_hybrid<-exp(-1*lambda_hybrid*month)
y_eastwest<-exp(-1*lambda_eastwest*month)
plot(month, 100*y_hybrid, type='b', ylim=c(0,100), col="black", xlab="Months", ylab="% Alive")
points(month, 100*y_eastwest, type='b', col="red")
legend("bottomleft", legend=c("Hybrid", "Eastern and Western"), col=c("black", "red"), lty=c(1,1))

#By month 12, 50.72% of hybrid seeds have survived and germinated
#By month 12, 38.8% of Eastern and Western seeds have survived 

seed_surv<-.5072

tau_2_hybrid<-0.5072
tau_2_E_W<-0.388

p.vec_overall[39]<- 0.5072 #Selected the highest survival probability (could have selected lowest)
p.vec_E[39]<-0.388
p.vec_H[39]<-0.5072
p.vec_W[39]<-0.388


#Save p.vec's for later use: 
save(p.vec_overall, file="./Overall/p.vec_overall.RData")
save(p.vec_E, file="./Eastern/p.vec_E.RData")
save(p.vec_H, file="./Hybrid/p.vec_H.RData")
save(p.vec_W, file="./Western/p.vec_W.RData")

save.image(file="demography_models.RData")
