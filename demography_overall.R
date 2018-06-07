#####Code for performing demographic analysis 
### By Kelley D. Erickson and Carol C. Horvitz

#Load packages 
library(MASS)
library(mvtnorm)

##Set working directory to load data: 
#MAC
setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 
#WINDOWS:
setwd("C:/Users/KE2/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 


load("demography_overall.RData")


demog<-read.csv("demography_15_clean.csv", head=T)


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










#Max and Min are highly correlated with eachother


###Use status matrix to get survival matrix
statmat<-as.matrix(demog[, c(6, 14, 22, 30, 38, 46, 54)])
survmat<-statmat
survmat<- ifelse(survmat=="prestudy"|survmat=="missing"|survmat=="still_dead", NA, survmat)
survmat<- ifelse(survmat =="tagged"|survmat=="alive", 1, 0)

survmat_pooled<-rbind( cbind(survmat[,1], survmat[,2]), cbind(survmat[,2], survmat[,3]), cbind(survmat[, 3], survmat[,4]), cbind(survmat[,4], survmat[,5]), cbind(survmat[,5], survmat[,6]), cbind(survmat[, 6], survmat[, 7]))

names(survmat_pooled)<-c("surv_t", "surv_t_plus_one")


survmat_pooled_bysize<-cbind(Diams_pooled[, 1], survmat_pooled[,2])
names(survmat_pooled_bysize)<-c("size_t", "surv_t_plus_one")



###Use reproduction status matrix to get reproduction matrix
repstatmat<-as.matrix(demog[, c(8, 16, 24, 32, 40, 48, 56)])
repmat<-repstatmat
repmat<- ifelse(repmat =="female", 1, 0)
repmat_pooled<-rbind( cbind(repmat[,1], repmat[,2]), cbind(repmat[,2], repmat[,3]), cbind(repmat[, 3], repmat[,4]), cbind(repmat[,4], repmat[,5]), cbind(repmat[,5], repmat[,6]), cbind(repmat[, 6], repmat[, 7]))



pooled<-cbind(ID_pooled, Sites_pooled, Location, Death_Cause, Genetic_type_pooled, Mins_pooled, Maxs_pooled, Diams_pooled, Heights_pooled, repmat_pooled, survmat_pooled)



names(pooled)<-c("ID", "Site","Location", "Death_Cause", "Genetic_type", "Min_t", "Min_tplus1", "Max_t", "Max_tplus1", "Diameter_t", "Diameter_tplus1", "Height_t", "Height_tplus1", "Rep_t", "Rep_tplus1", "Surv_t", "Surv_tplus1")

predictors<-cbind(pooled$Min_t, pooled$Diameter_t, pooled$Height_t, pooled$Max_t)

correlation<-cor(predictors[, 1:4], use="pairwise.complete.obs")
colnames(correlation)<-c("Min", "Diameter", "Height", "Max")
rownames(correlation)<-c("Min", "Diameter", "Height", "Max")

correlation
#Over the entire dataset, there is the lowest correlation between diameter and height
#(In lab group we were looking at correlation matrix for census 1 only )
write.csv(correlation, file="correlation_matrix_size.csv")


#Boxplot of sizes of individuals that survive to t plus 1 and those that do not
boxplot(pooled$Diameter_t[pooled$Surv_tplus1==0], pooled$Diameter_t[pooled$Surv_tplus1==1], names=c("Don't survive", "Survive to t+1"), ylab="Diameter (cm) at time t")


#Boxplots of sizes of individuals who reproduce and those that do not 
boxplot(pooled$Diameter_t[pooled$Rep_t==0], pooled$Diameter_t[pooled$Rep_t==1], names=c("Not a repro. female at t", "repro. female at t"), ylab="Diameter (mm) at time t")


Diam_cat_t<-rep(0, 9072 )
Height_cat_t<-rep(0, 9072)
Min_cat_t<-rep(0, 9072)

for (i in 1:9072) { 
	if (is.na(pooled$Diameter_t[i])) Diam_cat_t[i]<-NA 
	else if (pooled$Diameter_t[i]<100) Diam_cat_t[i]<-"(0, 100mm)" 
	else if (pooled$Diameter_t[i]<200) Diam_cat_t[i]<-"[100, 200mm)" 
	else if (pooled$Diameter_t[i]<300) Diam_cat_t[i]<-"[200, 300mm)" 
	else if (pooled$Diameter_t[i]<400) Diam_cat_t[i]<-"[300, 400mm)"
	else if (pooled$Diameter_t[i]<500) Diam_cat_t[i]<-"[400, 500mm)"
	else if (pooled$Diameter_t[i]<600) Diam_cat_t[i]<-"[500, 600mm)"
	else if (pooled$Diameter_t[i]<700) Diam_cat_t[i]<-"[600, 700mm)" 
	else if (pooled$Diameter_t[i]<800) Diam_cat_t[i]<-"[700, 800mm)"
	else Diam_cat_t[i]<-"<=800mm"
}

for (i in 1:9072) { 
	if (is.na(pooled$Height_t[i])) Height_cat_t[i]<-NA 
	else if (pooled$Height_t[i]<100) Height_cat_t[i]<-"(0, 100cm)" 
	else if (pooled$Height_t[i]<200) Height_cat_t[i]<-"[100, 200cm)" 
	else if (pooled$Height_t[i]<300) Height_cat_t[i]<-"[200, 300cm)" 
	else if (pooled$Height_t[i]<400) Height_cat_t[i]<-"[300, 400cm)"
	else if (pooled$Height_t[i]<500) Height_cat_t[i]<-"[400, 500cm)"
	else if (pooled$Height_t[i]<600) Height_cat_t[i]<-"[500, 600cm)"
	else if (pooled$Height_t[i]<700) Height_cat_t[i]<-"[600, 700cm)" 
	else if (pooled$Height_t[i]<800) Height_cat_t[i]<-"[700, 800cm)"
	else Height_cat_t[i]<-"<=800cm"
}

for (i in 1:9072) { 
	if (is.na(pooled$Min_t[i])) Min_cat_t[i]<-NA 
	else if (pooled$Min_t[i]<100) Min_cat_t[i]<-"(0, 100cm)" 
	else if (pooled$Min_t[i]<200) Min_cat_t[i]<-"[100, 200cm)" 
	else if (pooled$Min_t[i]<300) Min_cat_t[i]<-"[200, 300cm)" 
	else if (pooled$Min_t[i]<400) Min_cat_t[i]<-"[300, 400cm)"
	else if (pooled$Min_t[i]<500) Min_cat_t[i]<-"[400, 500cm)"
	else if (pooled$Min_t[i]<600) Min_cat_t[i]<-"[500, 600cm)"
	else if (pooled$Min_t[i]<700) Min_cat_t[i]<-"[600, 700cm)" 
	else if (pooled$Min_t[i]<800) Min_cat_t[i]<-"[700, 800cm)"
	else Min_cat_t[i]<-"<=800cm"
}

pooled2<-cbind(pooled, Diam_cat_t, Height_cat_t, Min_cat_t)








#Remove individual deaths that were sprayed or fire or grinding
restrict<-subset(pooled2, is.na(Death_Cause))

#Remove outliers (very large diameters)
restrict2<-subset(restrict, Diameter_t<800)

diam_table<-table(restrict2$Surv_tplus1, restrict2$Diam_cat_t)
surv_by_diam<-diam_table[2,]/(diam_table[1,]+diam_table[2,])


#In these tables below I have already removed deaths by "unnatural" causes as well as outliers (very large diameters >800mm)
height_table<-table(restrict2$Surv_tplus1, restrict2$Height_cat_t)
surv_by_height<-height_table[2,]/(height_table[1,] + height_table[2,])

min_table<-table(restrict2$Surv_tplus1, restrict2$Min_cat_t)
surv_by_min<-min_table[2,]/(min_table[1,] + min_table[2,])






write.csv(height_table, file="height_table.csv")
write.csv(min_table, file="min_table.csv")
write.csv(diam_table, file="diam_table.csv")

write.csv(surv_by_diam, file="surv_by_diam.csv")
write.csv(surv_by_height, file="surv_by_height.csv")
write.csv(surv_by_min, file="surv_by_min.csv")

#Divide dataset into size domains (seedlings and larges)
seedlings<-subset(restrict2, restrict2$Diameter_t<1.6)
seedlings<-subset(seedlings, seedlings$Height_t<16)
larges<-subset(restrict2, restrict2$Diameter_t>1.6)
larges<-subset(larges, larges$Height_t>16)

x1seq_larges<-seq(1.6, 800, length.out=50)
x2seq_larges<-seq(16, 800, length.out=50)

x1seq_seedlings<-seq(0, 1.6, length.out=50)
x2seq_seedlings<-seq(0, 16, length.out=50)



#Create a function for plotting logistic equations (TODO: Do I actually use this function?)
plot_logistic=function(model, xlab, ylab, x1seq, x2seq){
b0<-model$coeff[1]
b1<-model$coeff[2]
b2<-model$coeff[3]


z1<-outer(x1seq, x2seq, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]

x11()
persp(x1seq, x2seq, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.55, nticks=4, xlab=xlab, ylab=ylab, border=NA, col = color[facetcol])
}




#Predicting survival as a function of diameter and growth (removed outliers with large diameters, as well as "unnatural" death, ie fire and spraying)
# I do not have enough data for an interaction term 

#Full model without interaction:

#Diameter and Height

#Model survival of seedlings in the D1 domain(not separated by biotype)
surv_mod_seedlings<-glm(seedlings$Surv_tplus1~ seedlings$Diameter_t + seedlings$Height_t, family=binomial)

summary(surv_mod_seedlings)
#All terms are significant


b0<-surv_mod_seedlings$coeff[1]
b1<-surv_mod_seedlings$coeff[2]
b2<-surv_mod_seedlings$coeff[3]


z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]



setwd("/Users/curculion/Dropbox/UM_Dissertation_LaTeX-master/untitled folder/Figures/New_July")
png(file="seedling_surv_overall.png", width=10, height=10, units="in", res=300)
par(ps=24)
persp(x1seq_seedlings, x2seq_seedlings, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.55, nticks=4, xlab="\n Diameter (mm)", ylab="\n \n Height (cm)", col = color[facetcol])
dev.off()


#Model D2 survival (not by biotype)
surv_mod_larges<-glm(larges$Surv_tplus1~larges$Diameter_t + larges$Height_t, family=binomial)
summary(surv_mod_larges)

b0<-surv_mod_larges$coeff[1]
b1<-surv_mod_larges$coeff[2]
b2<-surv_mod_larges$coeff[3]


z1<-outer(x1seq_larges, x2seq_larges, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]


png(file="D2_surv_overall.png", width=10, height=10, units="in", res=300)
par(ps=24)
persp(x1seq_larges, x2seq_larges, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.55, nticks=4, xlab="\n Diameter (mm)", ylab="\n \n Height (cm)", col = color[facetcol])
dev.off()


####GROWTH #####

#Height at t+1 is a function of height_t AND diam_t


#Model growth on the D1 domain: 
growth_mod_seedlings<-manova(cbind(seedlings$Diameter_tplus1, seedlings$Height_tplus1) ~ seedlings$Diameter_t+seedlings$Height_t)

summary(growth_mod_seedlings)


#Plot diameter tplus1
b0<-growth_mod_seedlings$coeff[1,1]
b1<-growth_mod_seedlings$coeff[2,1]
b2<-growth_mod_seedlings$coeff[3,1]


z_diam<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) (b0+b1*a + b2*b))

#Plot height tplus1
b0<-growth_mod_seedlings$coeff[1,2]
b1<-growth_mod_seedlings$coeff[2,2]
b2<-growth_mod_seedlings$coeff[3,2]

z_height<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) (b0+b1*a + b2*b))


nrz<-nrow(z_diam)
ncz<-ncol(z_diam)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color_diam<-jet.colors(nbcol)
zfacet_diam<-z_diam[-1, -1] + z_diam[-1, -ncz] + z_diam[-nrz, -1] + z_diam[-nrz, -ncz]
facetcol_diam<-cut(zfacet_diam, nbcol)


nrz<-nrow(z_height)
ncz<-ncol(z_height)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color_height<-jet.colors(nbcol)
zfacet_height<-z_height[-1, -1] + z_height[-1, -ncz] + z_height[-nrz, -1] + z_height[-nrz, -ncz]
facetcol_height<-cut(zfacet_height, nbcol)


png(file="growth_seedlings_overall.png", width=8, height=16, units="in", res=300)
par(mfrow=c(2, 1), ps=24)
persp(x1seq_seedlings, x2seq_seedlings, z_diam, theta=-30, xlab="\n Diameter at t (mm)", ylab="\n \n Height at t (cm)",ticktype="detailed", col = color_diam[facetcol_diam], zlab="\n \n Diameter at t + 1", main="(A) Diameter at t + 1", nticks=4) 

persp(x1seq_seedlings, x2seq_seedlings, z_height, theta=-30, xlab="\n Diameter at t (mm)", ylab="\n \n Height at t (cm)", ticktype="detailed",  col= color_height[facetcol_height], zlab="\n \n Height at t + 1", main="(B) Height at t + 1", nticks=4)  
dev.off()


#Model growth of individuals in the D2 domain (overall)

growth_mod_larges<-manova(cbind(larges$Diameter_tplus1, larges$Height_tplus1) ~ larges$Diameter_t+larges$Height_t)

summary(growth_mod_larges)

#Plot diameter tplus1
b0<-growth_mod_larges$coeff[1,1]
b1<-growth_mod_larges$coeff[2,1]
b2<-growth_mod_larges$coeff[3,1]


z_diam<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))

#Plot height tplus1
b0<-growth_mod_larges$coeff[1,2]
b1<-growth_mod_larges$coeff[2,2]
b2<-growth_mod_larges$coeff[3,2]

z_height<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))


nrz<-nrow(z_diam)
ncz<-ncol(z_diam)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color_diam<-jet.colors(nbcol)
zfacet_diam<-z_diam[-1, -1] + z_diam[-1, -ncz] + z_diam[-nrz, -1] + z_diam[-nrz, -ncz]
facetcol_diam<-cut(zfacet_diam, nbcol)


nrz<-nrow(z_height)
ncz<-ncol(z_height)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color_height<-jet.colors(nbcol)
zfacet_height<-z_height[-1, -1] + z_height[-1, -ncz] + z_height[-nrz, -1] + z_height[-nrz, -ncz]
facetcol_height<-cut(zfacet_height, nbcol)



png(file="growth_D2_overall.png", width=8, height=16, units="in", res=300)
par(mfrow=c(2, 1), ps=24)
persp(x1seq_larges, x2seq_larges, z_diam, theta=-30, xlab="\n Diameter at t (mm)", ylab="\n \n Height at t (cm)",ticktype="detailed", col = color_diam[facetcol_diam], zlab="\n \n Diameter at t + 1", main="(A) Diameter at t + 1", nticks=4) 

persp(x1seq_larges, x2seq_larges, z_height, theta=-30, xlab="\n Diameter at t (mm)", ylab="\n \n Height at t (cm)", ticktype="detailed",  col= color_height[facetcol_height], zlab="\n \n Height at t + 1", main=" (B) Height at t + 1", nticks=4)  
dev.off()


#Capture variance in growth for modeling size distribution in D2 at time t+1:
growth_mod_diam_seedlings<-lm(seedlings$Diameter_tplus1 ~ seedlings$Diameter_t+seedlings$Height_t)
summary(growth_mod_diam_seedlings)

growth_mod_diam_seedlings2<-lm(seedlings$Diameter_tplus1~ 0+seedlings$Diameter_t)
summary(growth_mod_diam_seedlings2)

sigma_growth_diam_seedlings<-summary(growth_mod_diam_seedlings)$sigma
sigma_growth_diam_seedlings2<-summary(growth_mod_diam_seedlings2)$sigma

growth_mod_diam_larges<-lm(larges$Diameter_tplus1 ~ larges$Diameter_t+larges$Height_t)
summary(growth_mod_diam_larges)

growth_mod_diam_larges2<-lm(larges$Diameter_tplus1~ larges$Diameter_t)
summary(growth_mod_diam_larges2)

sigma_growth_diam_larges<-summary(growth_mod_diam_larges)$sigma
sigma_growth_diam_larges2<-summary(growth_mod_diam_larges2)$sigma


growth_mod_height_seedlings<-lm(seedlings$Height_tplus1 ~ seedlings$Diameter_t + seedlings$Height_t)
summary(growth_mod_height_seedlings)

growth_mod_height_larges<-lm(larges$Height_tplus1 ~ larges$Diameter_t + larges$Height_t)
summary(growth_mod_height_larges)


sigma_growth_height_seedlings<-summary(growth_mod_height_seedlings)$sigma
sigma_growth_height_larges<-summary(growth_mod_height_larges)$sigma




#########Reproduction########
setwd("/Users/curculion/Dropbox/UM_Dissertation_LaTeX-master/untitled folder/Figures/New_July")

#(i) Probability of being a reproductive female at time t: 


mod_repro<-glm(larges$Rep_tplus1 ~ larges$Diameter_t+larges$Height_t, family=binomial)

summary(mod_repro)
#Diameter term is not significant, so remove it: 
mod_repro<-glm(larges$Rep_tplus1 ~ larges$Height_t, family=binomial)

setwd("/Users/curculion/Dropbox/UM_Dissertation_LaTeX-master/untitled folder/Figures/New_July")


b0<-mod_repro$coefficients[1]
b1<-mod_repro$coefficients[2]
y<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))

png(file="reproduction_overall.png", width=10, height=10, units="in", res=300)
par(ps=24)
par(mar=c(5, 5, 4, 2))
plot(x2seq_larges, y, xlab="Height at time t (cm)", ylab="\n \n Probability of Reproducing", lwd=3, type='l', lty=5)
dev.off()




#(ii) Given that an individual is a reproductive female of size[diameter=x, height=y] and genetic type z[known from site], how many fruits do they produce? 
#Depends on diameter and genetic type only (from biomass allocation chapter)

setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Allometry")
load("biomass_170315.RData")

setwd("/Users/curculion/Dropbox/UM_Dissertation_LaTeX-master/untitled folder/Figures/New_July")

LHS$diam_base<-LHS$diam_base*10 #convert from cm to mm
fruit_mod<-lm(LHS$Fruit_No~LHS$diam_base + LHS$Height)
summary(fruit_mod)
#Only diameter:
fruit_mod2<-lm(LHS$Fruit_No~0+ LHS$diam_base)


b1<-fruit_mod2$coeff[1]
y<-x1seq_larges*x1seq_larges*b1


png(file="fecundity_overall.png", width=10, height=10, units="in", res=300)
par(ps=24)
par(mar=c(5, 5, 4, 2))
plot(x1seq_larges, y, xlab="Diameter at time t (mm)", ylab="\n \n Number of Offspring", lwd=3, type='l')
dev.off()


#(iii) Recruitment of seedlings from fruits..... 
#(a) P(survival) from seed at time t to seedling at time t + 1: 

#*Do have literature based report of 8mo survival rate from seeds sourced from the six sites grown in a common garden experiment in Davie: 
# "A greater proportion of hybrid seedlings survived than did western seedlings (63.6% vs. 53.2%; see fig. 5)." (Geiger et al. 2011) 


##Assumption: Assume mortality is linear
##Monthly mortality: 
#	Hybrids: (100-63.6)/8 = 4.55
#	Eastern/Western: (100-53.2)/8 = 5.85

month<-seq(0, 12, by=1)
hybrid_seed_surv<-rep(100, 12)
eastwest_seed_surv<-rep(100, 12)

for(i in 1:12) {
	hybrid_seed_surv[i+1]<-hybrid_seed_surv[i]-4.55
	eastwest_seed_surv[i+1]<-eastwest_seed_surv[i]-5.85
	}



par(mfrow=c(1, 2))
plot(month, hybrid_seed_surv, type='b', col="black", xlab="Months", ylim=c(0,100), ylab="% Alive")
points(month, eastwest_seed_surv, type='b', col="red")
legend("bottomleft", legend=c("Hybrid", "Eastern and Western"), col=c("black", "red"), lty=c(1,1))

#By month 12, 45.4% of hybrid seeds have survived and germinated
#By month 12, 29.8% of Eastern and Western seeds have survived and germinated 



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


###Notes from meeting with Carol on 3/21: In a logistic regression, the residuals are calculated as -2 * log deviance (the log of a large number is small!) Because predicted values are either 0s or 1s, we wind up with the two lines (as I see). There is no cause for alarm in my model checking plots. (There are probably more deep things going on, but will not worry about them too much right now). 
# http://stats.stackexchange.com/questions/1432/what-do-the-residuals-in-a-logistic-regression-mean
# https://www.r-bloggers.com/residuals-from-a-logistic-regression/ 



#Estimating size distribution of seedlings at time t +1: 
#Figure out distribution of newly tagged individuals in seedling plots 


recruit4_diam<-demog$Diam_4[demog$Status_4=="tagged"&demog$Location=="SeedlingPlot"]
recruit5_diam<-demog$Diam_5[demog$Status_5=="tagged"&demog$Location=="SeedlingPlot"]
recruit6_diam<-demog$Diam_6[demog$Status_6=="tagged"&demog$Location=="SeedlingPlot"]
recruits_diam<-c(recruit4_diam, recruit5_diam, recruit6_diam)

recruit4_height<-demog$Height_4[demog$Status_4=="tagged"&demog$Location=="SeedlingPlot"]
recruit5_height<-demog$Height_5[demog$Status_5=="tagged"&demog$Location=="SeedlingPlot"]
recruit6_height<-demog$Height_6[demog$Status_6=="tagged"&demog$Location=="SeedlingPlot"]
recruits_height<-c(recruit4_height, recruit5_height, recruit6_height)


#Not confident that ALL of these are actually new recruits (some are way too tall)

#thing<-boxplot(recruits_diam, ylab="Diameter (mm)")
#thing2<-boxplot(recruits_height, ylab="Height(cm)", notch=T)
#thing2$stats
#boxplot.stats(recruits_height)$stats[5]
#quantile(recruits_height, na.rm=T)

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


#Graduation into large domain: 
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

sdlng_grad_mod<-glm(seedlings$grad_status~ seedlings$Diameter_t + seedlings$Height_t, family=binomial)

summary(sdlng_grad_mod)

b0<-sdlng_grad_mod $coeff[1]
b1<-sdlng_grad_mod $coeff[2]
b2<-sdlng_grad_mod $coeff[3]

z_sdling_grad<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

nrz<-nrow(z_sdling_grad)
ncz<-ncol(z_sdling_grad)
nbcol<-100
jet.colors<-colorRampPalette(c("blue", "green"))
color<-jet.colors(nbcol)
zfacet<-z_sdling_grad[-1, -1] + z_sdling_grad[-1, -ncz] + z_sdling_grad[-nrz, -1] + z_sdling_grad[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


png(file="graduation_overall.png", width=10, height=10, units="in", res=300)
par(ps=24)
persp(x1seq_seedlings, x2seq_seedlings, z_sdling_grad, ticktype="detailed", theta=-30, zlim=c(0, 1.02), xlab="\n Diameter (mm)", ylab="\n Height (cm)", zlab="\n P(graduating)", col = color[facetcol])
dev.off()



graduates<-subset(seedlings, seedlings$grad_status==1)


#Get distribution of new recruits sizes 
summary(graduates$Diameter_tplus1)
hist(graduates$Diameter_tplus1)
mu_grad_diam<-mean(graduates$Diameter_tplus1)
sd_grad_diam<-sd(graduates$Diameter_tplus1)

summary(graduates$Height_tplus1)
hist(graduates$Height_tplus1)
mu_grad_height<-mean(graduates$Height_tplus1)
sd_grad_height<-sd(graduates$Height_tplus1)

save.image("demography_overall.RData")