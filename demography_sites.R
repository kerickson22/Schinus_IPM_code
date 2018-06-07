#Code for running the demographic analyses incorporating biotype differences 
#Code written by Kelley D. Erickson and Carol C. Horvitz



library(MASS)
library(mvtnorm)
library(RColorBrewer)
rm(list=ls(all=TRUE))
#MAC
setwd("/Users/kelley/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 

#WINDOWS:
setwd("C:/Users/KE2/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 

colors<-brewer.pal("Paired", n=10)
display.brewer.pal("Paired", n=10)

col_E_light<-colors[1]
col_E_dark<-colors[2]

col_W_light<-colors[5]
col_W_dark<-colors[6]

col_H_light<-colors[9]
col_H_dark<-colors[10]




data<-read.csv("demography_15_clean.csv", head=T)


Death_Cause<-c(data$Death_Cause_2, data$Death_Cause_3, data$Death_Cause_4, data$Death_Cause_5, data$Death_Cause_6, data$Death_Cause_7)


Status<-c(data$Status_1, data$Status_2, data$Status_3, data$Status_4, data$Status_5, data$Status_6, data$Status_7)
#Status2<-cbind.data.frame(data$Status_1, data$Status_2, data$Status_3, data$Status_4, data$Status_5, data$Status_6, data$Status_7)
# row.names(Status)<-data$Plant_ID

#Why is cbind converting to numeric? #Have to use cbind.data.frame 
#Since there were some duplicated tag#s, went into Excel and recoded (A lot of the tag numbers at Cocoa were repeated at other sites 

Diams<-cbind(data$Diam_1, data$Diam_2, data$Diam_3, data$Diam_4, data$Diam_5, data$Diam_6, data$Diam_7) 
row.names(Diams)<-data$Plant_ID

Heights<-cbind(data$Height_1, data$Height_2, data$Height_3, data$Height_4, data$Height_5, data$Height_6, data$Height_7)
row.names(Heights)<-data$Plant_ID

Mins<-cbind(data$Min_1, data$Min_2, data$Min_3, data$Min_4, data$Min_5, data$Min_6, data$Min_7)
row.names(Mins)<-data$Plant_ID

Maxs<-cbind(data$Max_1, data$Max_2, data$Max_3, data$Max_4, data$Max_5, data$Max_6, data$Max_7)
row.names(Maxs)<-data$Plant_ID


Diams_pooled<-rbind( cbind(Diams[,1], Diams[,2]), cbind(Diams[,2], Diams[,3]), cbind(Diams[, 3], Diams[,4]), cbind(Diams[,4], Diams[,5]), cbind(Diams[,5], Diams[,6]), cbind(Diams[, 6], Diams[, 7]))


Diams_pooled<-data.frame(Diams_pooled[,1], Diams_pooled[,2])
names(Diams_pooled)<-c("t", "t_plus_one")


Sites_pooled<-rep(data$Site, 6)
ID_pooled<-rep(data$Plant_ID, 6)
Location<-rep(data$Location, 6)

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
statmat<-as.matrix(data[, c(6, 14, 22, 30, 38, 46, 54)])
survmat<-statmat
survmat<- ifelse(survmat=="prestudy"|survmat=="missing"|survmat=="still_dead", NA, survmat)
survmat<- ifelse(survmat =="tagged"|survmat=="alive", 1, 0)

survmat_pooled<-rbind( cbind(survmat[,1], survmat[,2]), cbind(survmat[,2], survmat[,3]), cbind(survmat[, 3], survmat[,4]), cbind(survmat[,4], survmat[,5]), cbind(survmat[,5], survmat[,6]), cbind(survmat[, 6], survmat[, 7]))

names(survmat_pooled)<-c("surv_t", "surv_t_plus_one")


survmat_pooled_bysize<-cbind(Diams_pooled[, 1], survmat_pooled[,2])
names(survmat_pooled_bysize)<-c("size_t", "surv_t_plus_one")



###Use reproduction status matrix to get reproduction matrix
repstatmat<-as.matrix(data[, c(8, 16, 24, 32, 40, 48, 56)])
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
par(ps=24)
persp(x1seq, x2seq, z1, ticktype="detailed", theta=-40, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.55, nticks=4, xlab=xlab, ylab=ylab, border=NA, col = color[facetcol])
}


p.vec_E<-rep(0, 39);
p.vec_H<-rep(0, 39);
p.vec_W<-rep(0, 39);


#Predicting survival as a function of diameter and growth (removed outliers with large diameters, as well as "unnatural" death, ie fire and spraying)
# I do not have enough data for an interaction term 

#Full model without interaction:

#Diameter and Height
surv_mod_seedlings_site<-glm(seedlings$Surv_tplus1~ seedlings$Diameter_t + seedlings$Height_t+seedlings$Site, family=binomial)
summary(surv_mod_seedlings_site)

surv_mod_seedlings_lineage<-glm(seedlings$Surv_tplus1 ~ seedlings$Diameter_t + seedlings$Height_t + seedlings$Genetic_type, family=binomial)
summary(surv_mod_seedlings_lineage)
#Hybrid is different from E+W, but Eastern and Western are not different from eachother

lineage2<-rep(0, length(seedlings$Genetic_type))
for (i in 1:length(seedlings$Genetic_type)) {
	if (seedlings$Genetic_type[i]=="western" || seedlings$Genetic_type[i]=="eastern") {lineage2[i]<-"E-W"}
	else{lineage2[i]<-"Hybrid"}
}

seedlings2<-cbind(seedlings, lineage2)

surv_mod_seedlings_lineage2<-glm(seedlings2$Surv_tplus1 ~ seedlings2$Diameter_t + seedlings2$Height_t + seedlings2$lineage2, family=binomial)
summary(surv_mod_seedlings_lineage2)


surv_mod_seedlings_lineage2$coeff

p.vec_E[20]<-surv_mod_seedlings_lineage2$coeff[1]
p.vec_E[21]<-surv_mod_seedlings_lineage2$coeff[2]
p.vec_E[22]<-surv_mod_seedlings_lineage2$coeff[3]

p.vec_H[20]<-surv_mod_seedlings_lineage2$coeff[1]+surv_mod_seedlings_lineage2$coeff[4]
p.vec_H[21]<-surv_mod_seedlings_lineage2$coeff[2]
p.vec_H[22]<-surv_mod_seedlings_lineage2$coeff[3]

p.vec_W[20]<-surv_mod_seedlings_lineage2$coeff[1]
p.vec_W[21]<-surv_mod_seedlings_lineage2$coeff[2]
p.vec_W[22]<-surv_mod_seedlings_lineage2$coeff[3]

###Plotting Eastern seedling survival
x1seq_seedling<-seq(0, 1.6, length.out=50)
x2seq_seedling<-seq(0, 16, length.out=50)


png(file="D1_survival_biotype.png", width=10, height=20, units="in", res=300)
par(ps=24, mfrow=c(2, 1))
b0<-p.vec_E[20]
b1<-p.vec_E[21]
b2<-p.vec_E[22]


z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c(col_E_light, col_E_dark))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]


persp(x1seq_seedling, x2seq_seedling, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", border=col_E_dark, col = color[facetcol],  lwd=3, xaxs="i", main="Eastern/Western", cex.main=1)




###Plotting Hybrid seedling survival
x1seq_seedling<-seq(0, 1.6, length.out=50)
x2seq_seedling<-seq(0, 16, length.out=50)

b0<-p.vec_H[20]
b1<-p.vec_H[21]
b2<-p.vec_H[22]


z1<-outer(x1seq_seedling, x2seq_seedling, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c(col_H_light, col_H_dark))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]

persp(x1seq_seedling, x2seq_seedling, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", border=col_H_dark, col = color[facetcol],  lwd=3, xaxs="i", main="Hybrid", cex.main=1)

dev.off()

#####Larges Survival by Lineage

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

p.vec_E[1]<-surv_mod_larges_lineage2$coeff[1]
p.vec_E[2]<-surv_mod_larges_lineage2$coeff[2]
p.vec_E[3]<-surv_mod_larges_lineage2$coeff[3]

p.vec_H[1]<-surv_mod_larges_lineage2$coeff[1] + surv_mod_larges_lineage2$coeff[4]
p.vec_H[2]<-surv_mod_larges_lineage2$coeff[2]
p.vec_H[3]<-surv_mod_larges_lineage2$coeff[3]

p.vec_W[1]<-surv_mod_larges_lineage2$coeff[1]
p.vec_W[2]<-surv_mod_larges_lineage2$coeff[2]
p.vec_W[3]<-surv_mod_larges_lineage2$coeff[3]

###
x1seq_larges<-seq(1.6, 700,  length.out=50)
x2seq_larges<-seq(16, 800, length.out=50)

png(file="D2_survival_biotype.png", width=10, height=20, units="in", res=300)
par(ps=24, mfrow=c(2, 1))
b0<-p.vec_E[1]
b1<-p.vec_E[2]
b2<-p.vec_E[3]


z1<-outer(x1seq_larges, x2seq_larges, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c(col_E_light, col_E_dark))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]


persp(x1seq_larges, x2seq_larges, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", border=col_E_dark, col = color[facetcol],  lwd=3, xaxs="i", main="Eastern/Western", cex.main=1)




###Plotting Hybrid seedling survival


b0<-p.vec_H[1]
b1<-p.vec_H[2]
b2<-p.vec_H[3]


z1<-outer(x1seq_larges, x2seq_larges, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c(col_H_light, col_H_dark))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]

persp(x1seq_larges, x2seq_larges, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", border=col_H_dark, col = color[facetcol],  lwd=3, xaxs="i", main="Hybrid", cex.main=1)

dev.off()





####GROWTH #####

#Height at t+1 is a function of height_t AND diam_t



growth_mod_seedlings_lineage<-manova(cbind(seedlings$Diameter_tplus1, seedlings$Height_tplus1) ~ seedlings$Diameter_t+seedlings$Height_t+seedlings$Genetic_type)
summary(growth_mod_seedlings_lineage)
#Lineage does not appear to have an effect on seedling growth

#Try grouping E+W again:

growth_mod_seedlings_lineage2<-manova(cbind(seedlings$Diameter_tplus1, seedlings$Height_tplus1) ~ seedlings$Diameter_t + seedlings$Height_t + seedlings2$lineage2)
summary(growth_mod_seedlings_lineage2)
#Lineage does not have an effect on seedling growth, so remove it: 


growth_mod_seedlings<-manova(cbind(seedlings$Diameter_tplus1, seedlings$Height_tplus1) ~ seedlings$Diameter_t + seedlings$Height_t)

summary(growth_mod_seedlings)

p.vec_E[23]<-growth_mod_seedlings$coeff[1,1]
p.vec_E[24]<-growth_mod_seedlings$coeff[2,1]
p.vec_E[25]<-growth_mod_seedlings$coeff[3,1]

p.vec_H[23]<-growth_mod_seedlings$coeff[1,1]
p.vec_H[24]<-growth_mod_seedlings$coeff[2,1]
p.vec_H[25]<-growth_mod_seedlings$coeff[3,1]

p.vec_W[23]<-growth_mod_seedlings$coeff[1,1]
p.vec_W[24]<-growth_mod_seedlings$coeff[2,1]
p.vec_W[25]<-growth_mod_seedlings$coeff[3,1]











growth_mod_diam_seedlings<-lm(seedlings$Diameter_tplus1 ~ seedlings$Diameter_t+seedlings$Height_t)
summary(growth_mod_diam_seedlings)
sigma_growth_diam_seedlings<-summary(growth_mod_diam_seedlings)$sigma

p.vec_E[26]<-sigma_growth_diam_seedlings
p.vec_H[26]<-sigma_growth_diam_seedlings
p.vec_W[26]<-sigma_growth_diam_seedlings


p.vec_E[27]<-growth_mod_seedlings$coeff[1,2]
p.vec_E[28]<-growth_mod_seedlings$coeff[2,2]
p.vec_E[29]<-growth_mod_seedlings$coeff[3,2]

p.vec_H[27]<-growth_mod_seedlings$coeff[1,2]
p.vec_H[28]<-growth_mod_seedlings$coeff[2,2]
p.vec_H[29]<-growth_mod_seedlings$coeff[3,2]

p.vec_W[27]<-growth_mod_seedlings$coeff[1,2]
p.vec_W[28]<-growth_mod_seedlings$coeff[2,2]
p.vec_W[29]<-growth_mod_seedlings$coeff[3,2]



growth_mod_height_seedlings<-lm(seedlings$Height_tplus1 ~ seedlings$Diameter_t + seedlings$Height_t)
summary(growth_mod_height_seedlings)

sigma_growth_height_seedlings<-summary(growth_mod_height_seedlings)$sigma

p.vec_E[30]<-sigma_growth_height_seedlings
p.vec_H[30]<-sigma_growth_height_seedlings
p.vec_W[30]<-sigma_growth_height_seedlings

### Plot Growth Mod Seedlings 

#Plot diameter tplus1
b0<-p.vec_E[23]
b1<-p.vec_E[24]
b2<-p.vec_E[25]


z_diam<-outer(x1seq_seedling, x2seq_seedling, function(a,b) (b0+b1*a + b2*b))

#Plot height tplus1
b0<-p.vec_E[27]
b1<-p.vec_E[28]
b2<-p.vec_E[29]

z_height<-outer(x1seq_seedling, x2seq_seedling, function(a,b) (b0+b1*a + b2*b))


nrz<-nrow(z_diam)
ncz<-ncol(z_diam)
nbcol<-100
jet.colors<-colorRampPalette(c("grey", "black"))
color_diam<-jet.colors(nbcol)
zfacet_diam<-z_diam[-1, -1] + z_diam[-1, -ncz] + z_diam[-nrz, -1] + z_diam[-nrz, -ncz]
facetcol_diam<-cut(zfacet_diam, nbcol)


nrz<-nrow(z_height)
ncz<-ncol(z_height)
nbcol<-100
jet.colors<-colorRampPalette(c("grey", "black"))
color_height<-jet.colors(nbcol)
zfacet_height<-z_height[-1, -1] + z_height[-1, -ncz] + z_height[-nrz, -1] + z_height[-nrz, -ncz]
facetcol_height<-cut(zfacet_height, nbcol)


x11(width=10, height=10)
par(ps=42)
par(mar=c(5, 9, 4, 2))
persp(x1seq_seedling, x2seq_seedling, z_diam, ticktype="detailed", theta=-30,  zlab="\n \n Diameter \n at t plus 1 (mm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_diam[facetcol_diam], cex.lab=1, cex.axis=0.8, lwd=1, xpd=T)
dev.copy(png, "growth_seedlings_diam.png", width=10, height=10, units="in", res=300)
dev.off()

x11(width=10, height=10)
par(ps=42)
par(mar=c(5, 9, 4, 2))
persp(x1seq_seedling, x2seq_seedling, z_height, ticktype="detailed", theta=-30,  zlab="\n \n Height \n at t plus 1 (cm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_height[facetcol_height], cex.lab=1, cex.axis=0.8, lwd=1, xpd=T)
dev.copy(png, "growth_seedlings_height.png", width=10, height=10, units="in", res=300)
dev.off()


###
















growth_mod_larges_lineage<-manova(cbind(larges$Diameter_tplus1, larges$Height_tplus1) ~ larges$Diameter_t+larges$Height_t+larges$Genetic_type)
summary(growth_mod_larges_lineage)
#BUT lineage DOES have an effect on the growth of the larger individuals in D2

p.vec_E[4]<-growth_mod_larges_lineage$coeff[1,1]
p.vec_E[5]<-growth_mod_larges_lineage$coeff[2,1]
p.vec_E[6]<-growth_mod_larges_lineage$coeff[3,1]

p.vec_H[4]<-growth_mod_larges_lineage$coeff[1,1] + growth_mod_larges_lineage$coeff[4,1]
p.vec_H[5]<-growth_mod_larges_lineage$coeff[2,1]
p.vec_H[6]<-growth_mod_larges_lineage$coeff[3,1]

p.vec_W[4]<-growth_mod_larges_lineage$coeff[1,1] + growth_mod_larges_lineage$coeff[5, 1]
p.vec_W[5]<-growth_mod_larges_lineage$coeff[2,1]
p.vec_W[6]<-growth_mod_larges_lineage$coeff[3,1]


growth_mod_diam_larges_lineage<-lm(larges$Diameter_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Genetic_type)
summary(growth_mod_diam_larges_lineage)
sigma_growth_diameter_larges_lineage<-summary(growth_mod_diam_larges_lineage)$sigma

p.vec_E[7]<-sigma_growth_diameter_larges_lineage
p.vec_H[7]<-sigma_growth_diameter_larges_lineage
p.vec_W[7]<-sigma_growth_diameter_larges_lineage


p.vec_E[8]<-growth_mod_larges_lineage$coeff[1,2]
p.vec_E[9]<-growth_mod_larges_lineage$coeff[2,2]
p.vec_E[10]<-growth_mod_larges_lineage$coeff[3,2]

p.vec_H[8]<-growth_mod_larges_lineage$coeff[1,2] + growth_mod_larges_lineage$coeff[4,2]
p.vec_H[9]<-growth_mod_larges_lineage$coeff[2,2]
p.vec_H[10]<-growth_mod_larges_lineage$coeff[3,2]

p.vec_W[8]<-growth_mod_larges_lineage$coeff[1,2] + growth_mod_larges_lineage$coeff[5,2]
p.vec_W[9]<-growth_mod_larges_lineage$coeff[2,2]
p.vec_W[10]<-growth_mod_larges_lineage$coeff[3,2]


growth_mod_height_larges_lineage<-lm(larges$Height_tplus1 ~ larges$Diameter_t + larges$Height_t + larges$Genetic_type)
summary(growth_mod_height_larges_lineage)
sigma_growth_height_larges_lineage<-summary(growth_mod_height_larges_lineage)$sigma


p.vec_E[11]<-sigma_growth_height_larges_lineage
p.vec_H[11]<-sigma_growth_height_larges_lineage
p.vec_W[11]<-sigma_growth_height_larges_lineage






###

### Plot Growth Mod Eastern Larges

png(file="growth_mod_biotype.png", width=20, height=30, units="in", res=300)
par(mfrow=c(3, 2), ps=42, mar=c(5.1, 11.1, 4.1, 2.1))
#Plot diameter tplus1
b0<-p.vec_E[4]
b1<-p.vec_E[5]
b2<-p.vec_E[6]


z_diam<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))

#Plot height tplus1
b0<-p.vec_E[8]
b1<-p.vec_E[9]
b2<-p.vec_E[10]

z_height<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))


nrz<-nrow(z_diam)
ncz<-ncol(z_diam)
nbcol<-100
jet.colors<-colorRampPalette(c(col_E_light, col_E_dark))
color_diam<-jet.colors(nbcol)
zfacet_diam<-z_diam[-1, -1] + z_diam[-1, -ncz] + z_diam[-nrz, -1] + z_diam[-nrz, -ncz]
facetcol_diam<-cut(zfacet_diam, nbcol)


nrz<-nrow(z_height)
ncz<-ncol(z_height)
nbcol<-100
jet.colors<-colorRampPalette(c(col_E_light, col_E_dark))
color_height<-jet.colors(nbcol)
zfacet_height<-z_height[-1, -1] + z_height[-1, -ncz] + z_height[-nrz, -1] + z_height[-nrz, -ncz]
facetcol_height<-cut(zfacet_height, nbcol)


persp(x1seq_larges, x2seq_larges, z_diam, ticktype="detailed", theta=-30,  zlab="\n \n Diameter \n at t plus 1 (mm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_diam[facetcol_diam], lwd=1, xpd=T, main="(A)", cex.main=1)

persp(x1seq_larges, x2seq_larges, z_height, ticktype="detailed", theta=-30,  zlab="\n \n \n Height \n at t plus 1 (cm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_height[facetcol_height], main="(B)", cex.main=1,  lwd=1, xpd=T)

###

### Plot Growth Mod Hybrid Larges

#Plot diameter tplus1
b0<-p.vec_H[4]
b1<-p.vec_H[5]
b2<-p.vec_H[6]


z_diam<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))

#Plot height tplus1
b0<-p.vec_H[8]
b1<-p.vec_H[9]
b2<-p.vec_H[10]

z_height<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))


nrz<-nrow(z_diam)
ncz<-ncol(z_diam)
nbcol<-100
jet.colors<-colorRampPalette(c(col_H_light, col_H_dark))
color_diam<-jet.colors(nbcol)
zfacet_diam<-z_diam[-1, -1] + z_diam[-1, -ncz] + z_diam[-nrz, -1] + z_diam[-nrz, -ncz]
facetcol_diam<-cut(zfacet_diam, nbcol)


nrz<-nrow(z_height)
ncz<-ncol(z_height)
nbcol<-100
jet.colors<-colorRampPalette(c(col_H_light, col_H_dark))
color_height<-jet.colors(nbcol)
zfacet_height<-z_height[-1, -1] + z_height[-1, -ncz] + z_height[-nrz, -1] + z_height[-nrz, -ncz]
facetcol_height<-cut(zfacet_height, nbcol)



persp(x1seq_larges, x2seq_larges, z_diam, ticktype="detailed", theta=-30,  zlab="\n \n Diameter \n at t plus 1 (mm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_diam[facetcol_diam], main="(C)", cex.main=1,  lwd=1, xpd=T)

persp(x1seq_larges, x2seq_larges, z_height, ticktype="detailed", theta=-30,  zlab="\n \n \n Height \n at t plus 1 (cm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_height[facetcol_height], main="(D)", cex.main=1,  lwd=1, xpd=T)

###

### Plot Growth Mod Western Larges

#Plot diameter tplus1
b0<-p.vec_W[4]
b1<-p.vec_W[5]
b2<-p.vec_W[6]


z_diam<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))

#Plot height tplus1
b0<-p.vec_W[8]
b1<-p.vec_W[9]
b2<-p.vec_W[10]

z_height<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0+b1*a + b2*b))


nrz<-nrow(z_diam)
ncz<-ncol(z_diam)
nbcol<-100
jet.colors<-colorRampPalette(c(col_W_light, col_W_dark))
color_diam<-jet.colors(nbcol)
zfacet_diam<-z_diam[-1, -1] + z_diam[-1, -ncz] + z_diam[-nrz, -1] + z_diam[-nrz, -ncz]
facetcol_diam<-cut(zfacet_diam, nbcol)


nrz<-nrow(z_height)
ncz<-ncol(z_height)
nbcol<-100
jet.colors<-colorRampPalette(c(col_W_light, col_W_dark))
color_height<-jet.colors(nbcol)
zfacet_height<-z_height[-1, -1] + z_height[-1, -ncz] + z_height[-nrz, -1] + z_height[-nrz, -ncz]
facetcol_height<-cut(zfacet_height, nbcol)

persp(x1seq_larges, x2seq_larges, z_diam, ticktype="detailed", theta=-30,  zlab="\n \n Diameter \n at t plus 1 (mm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_diam[facetcol_diam],main="(E)", cex.main=1,  lwd=1, xpd=T)


persp(x1seq_larges, x2seq_larges, z_height, ticktype="detailed", theta=-30,  zlab="\n \n \n Height \n at t plus 1 (cm) ", shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", border="black", col = color_height[facetcol_height], main="(F)", cex.main=1,  lwd=1, xpd=T)

dev.off()
###


#########Reproduction########


(i) Probability of being a reproductive female at time t: 


mod_repro_lineage<-glm(larges$Rep_tplus1 ~ larges$Diameter_t+larges$Height_t + larges$Genetic_type, family=binomial)

summary(mod_repro_lineage)
#Diameter term is not significant, so remove it: 
mod_repro_lineage<-glm(larges$Rep_tplus1 ~ larges$Height_t +larges$Genetic_type, family=binomial)
#Hybrid and Eastern are not different from eachother

lineage4<-rep(0, length(larges$Genetic_type))
for (i in 1:length(larges$Genetic_type)) {
	if (larges$Genetic_type[i]=="hybrid" || larges$Genetic_type[i]=="eastern") {lineage3[i]<-"E-H"}
	else{lineage4[i]<-"Western"}
}

larges3<-cbind(larges, lineage4)

mod_repro_lineage2<-glm(larges3$Rep_tplus1~ larges3$Height_t + larges3$lineage4, family=binomial)
summary(mod_repro_lineage2)
#Eastern and Hybrid have same probabilty of reproducing, but W differs

p.vec_E[12]<-mod_repro_lineage2$coeff[1]
p.vec_E[13]<-mod_repro_lineage2$coeff[2]

p.vec_H[12]<-mod_repro_lineage2$coeff[1]
p.vec_H[13]<-mod_repro_lineage2$coeff[2]

p.vec_W[12]<-mod_repro_lineage2$coeff[1]+mod_repro_lineage2$coeff[3]
p.vec_W[13]<-mod_repro_lineage2$coeff[2]



###Plotting Reproduction probability by site



b0<-p.vec_E[12]
b1<-p.vec_E[13]
y_E<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))

b0<-p.vec_H[12]
b1<-p.vec_H[13]
y_H<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))

b0<-p.vec_W[12]
b1<-p.vec_W[13]
y_W<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))




x11(width=10, height=10)

png(file="repro_by_site.png", width=10, height=10, units="in", res=300)
par(ps=24)
par(mar=c(5, 5, 4, 2))
plot(x2seq_larges, y_E, xlab="Height at time t (cm)", ylab="\n \n Probability of Reproducing", col=col_E_dark, lwd=3, type='l')
lines(x2seq_larges, y_H, col=col_H_dark, lwd=3, type='l')
lines(x2seq_larges, y_W, col=col_W_dark, lwd=3, type='l')
legend(600, 0.3, c("Eastern", "Hybrid", "Western"), col=c(col_E_dark, col_H_dark, col_W_dark), lty=c(1, 1, 1), y.intersp=2, lwd=c(3, 3, 3))
dev.off()






(ii) Given that an individual is a reproductive female of size[diameter=x, height=y] and genetic type z[known from site], how many fruits do they produce? 
#Depends on diameter and genetic type only (from biomass allocation chapter)

setwd("/Users/kelley/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Allometry")
load("biomass_170315.RData")
setwd("/Users/ke2/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data")

diam_base2<-Stems$diam_base*LHS$diam_base
mod<-lm(Stems$Seed_No ~ 0+diam_base2)

diam_base<-Stems$diam_base*10 #convert to cm
diam_base2<-diam_base*diam_base #square
mod<-lm(Stems$Seed_No ~ 0 + diam_base2)

mod_fruit_production<-lm(Stems$Seed_No ~ 0 + diam_base2:Stems$Genetic_Type)
summary(mod_fruit_production)

p.vec_E[14]<-mod_fruit_production$coeff[1]
p.vec_H[14]<-mod_fruit_production$coeff[2]
p.vec_W[14]<-mod_fruit_production$coeff[3]


y_E<-x1seq_larges*x1seq_larges*p.vec_E[14]
y_W<-x1seq_larges*x1seq_larges*p.vec_W[14]
y_H<-x1seq_larges*x1seq_larges*p.vec_H[14]

setwd("/Users/kelley/Dropbox/UM_Dissertation_LaTeX-master/untitled folder/Figures/biotype")

png(file="fecundity_biotype.png", width=10, height=10, units="in", res=300)
par(ps=24)
par(mar=c(5, 5, 4, 2))
plot(x1seq_larges, y_E, xlab="Diameter at time t (mm)", ylab="\n \n Number of Offspring", col=col_E_dark, lwd=3, type='l')
lines(x1seq_larges, y_H, col=col_H_dark, lwd=3, type='l')
lines(x1seq_larges, y_W, col=col_W_dark, lwd=3, type='l')
legend(550, 500000, c("Eastern", "Hybrid", "Western"), col=c(col_E_dark, col_H_dark, col_W_dark), lty=c(1, 1, 1), y.intersp=2, lwd=c(3, 3, 3))
dev.off()



png(file="repro_by_site.png", width=10, height=10, units="in", res=300)
par(ps=24)
par(mar=c(5, 5, 4, 2))
plot(x2seq_larges, y_E, xlab="Height at time t (cm)", ylab="\n \n Probability of Reproducing", col=col_E_dark, lwd=3, type='l')
lines(x2seq_larges, y_H, col=col_H_dark, lwd=3, type='l')
lines(x2seq_larges, y_W, col=col_W_dark, lwd=3, type='l')
legend(600, 0.3, c("Eastern", "Hybrid", "Western"), col=c(col_E_dark, col_H_dark, col_W_dark), lty=c(1, 1, 1), y.intersp=2, lwd=c(3, 3, 3))
dev.off()







(iii) Recruitment of seedlings from fruits..... 
(a) P(survival) from seed at time t to seedling at time t + 1: 

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


x11()
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

tau_2_hybrid<-0.5072
tau_2_E-W<-0.388

p.vec_E[19]<-0.388
p.vec_H[19]<-0.5072
p.vec_W[19]<-0.388


#Estimating size distribution of seedlings at time t +1: 
#Figure out distribution of newly tagged individuals in seedling plots 

#Need to reload demography data if you have loaded in the biomass data 

recruit4_diam<-data$Diam_4[data$Status_4=="tagged"&data$Location=="SeedlingPlot"]
recruit5_diam<-data$Diam_5[data$Status_5=="tagged"&data$Location=="SeedlingPlot"]
recruit6_diam<-data$Diam_6[data$Status_6=="tagged"&data$Location=="SeedlingPlot"]
recruits_diam<-c(recruit4_diam, recruit5_diam, recruit6_diam)

recruit4_height<-data$Height_4[data$Status_4=="tagged"&data$Location=="SeedlingPlot"]
recruit5_height<-data$Height_5[data$Status_5=="tagged"&data$Location=="SeedlingPlot"]
recruit6_height<-data$Height_6[data$Status_6=="tagged"&data$Location=="SeedlingPlot"]
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
boxplot(recruits_diam[recruits_diam<1.6], na.rm=T)
hist(recruits_diam[recruits_diam<1.6])
summary(recruits_diam[recruits_diam<1.6])
mu_diam<-mean(recruits_diam[recruits_diam<1.6], na.rm=T)
sd_diam<-sd(recruits_diam[recruits_diam<1.6], na.rm=T)
var_diam<-var(recruits_diam[recruits_diam<1.6], na.rm=T) 
length(recruits_diam[recruits_diam<1.6])

p.vec_E[15]<-mu_diam
p.vec_E[16]<-sd_diam

p.vec_H[15]<-mu_diam
p.vec_H[16]<-sd_diam

p.vec_W[15]<-mu_diam
p.vec_W[16]<-sd_diam


###Summary statistics for just the 481 individuals with a height<16cm
summary(recruits_height[recruits_height<16])
boxplot(recruits_height[recruits_height<16])
thing3<-hist(recruits_height[recruits_height<16])
mu_height<-mean(recruits_height[recruits_height<16], na.rm=T) 
sd_height<-sd(recruits_height[recruits_height<16],na.rm=T)
var_height<-var(recruits_height[recruits_height<16],na.rm=T) 
length(recruits_height[recruits_height<16])


p.vec_E[17]<-mu_height
p.vec_E[18]<-sd_height

p.vec_H[17]<-mu_height
p.vec_H[18]<-sd_height

p.vec_W[17]<-mu_height
p.vec_W[18]<-sd_height




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

sdlng_grad_mod_lineage<-glm(seedlings$grad_status~ seedlings$Diameter_t + seedlings$Height_t + seedlings$Genetic_type, family=binomial)
summary(sdlng_grad_mod_lineage)
#The effect of hybrids is only moderatately significantly diff from E, so collapse E and H




lineage4<-rep(0, length(seedlings$Genetic_type))
for (i in 1:length(seedlings$Genetic_type)) {
	if (seedlings$Genetic_type[i]=="eastern" || seedlings$Genetic_type[i]=="hybrid") {lineage4[i]<-"E-H"}
	else{lineage4[i]<-"Western"}
}

seedlings3<-cbind(seedlings, lineage4)

grad_mod_lineage<-glm(seedlings3$grad_status ~ seedlings3$Diameter_t + seedlings3$Height_t + seedlings3$lineage4, family=binomial)
summary(grad_mod_lineage)



p.vec_E[31]<-grad_mod_lineage$coeff[1]
p.vec_E[32]<-grad_mod_lineage$coeff[2]
p.vec_E[33]<-grad_mod_lineage$coeff[3]

p.vec_H[31]<-grad_mod_lineage$coeff[1]
p.vec_H[32]<-grad_mod_lineage$coeff[2]
p.vec_H[33]<-grad_mod_lineage$coeff[3]

p.vec_W[31]<-grad_mod_lineage$coeff[1]+grad_mod_lineage$coeff[4]
p.vec_W[32]<-grad_mod_lineage$coeff[2]
p.vec_W[33]<-grad_mod_lineage$coeff[3]




#Plot Graduation by Lineage 

png(file="graduation_biotype.png", width=10, height=20, units="in", res=300)
par(mfrow=c(2,1), ps=42, mar=c(5.1, 8.1, 4.1, 2.1))

#Eastern and Hybrid
b0<-p.vec_E[31]
b1<-p.vec_E[32]
b2<-p.vec_E[33]



z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c(col_E_light, col_E_dark))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]

persp(x1seq_seedlings, x2seq_seedlings, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Graduation)", shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", border="black", col = color[facetcol], main="(A) Eastern and Hybrid", cex.main=1, lwd=1, xaxs="i")


b0<-p.vec_W[31]
b1<-p.vec_W[32]
b2<-p.vec_W[33]



z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


nrz<-nrow(z1)
ncz<-ncol(z1)
nbcol<-100
jet.colors<-colorRampPalette(c(col_W_light, col_W_dark))
color<-jet.colors(nbcol)
zfacet<-z1[-1, -1] + z1[-1, -ncz] + z1[-nrz, -1] + z1[-nrz, -ncz]
facetcol<-cut(zfacet, nbcol)


col = color[facetcol]

persp(x1seq_seedlings, x2seq_seedlings, z1, ticktype="detailed", theta=-30, zlim=c(0, 1.25), zlab="\n P(Graduation)", shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", border="black", col = color[facetcol], main=" (B) Western", cex.main=1, lwd=1, xaxs="i")

dev.off()











graduates<-subset(seedlings, seedlings$grad_status==1)

summary(graduates$Diameter_tplus1)
hist(graduates$Diameter_tplus1)
mu_grad_diam<-mean(graduates$Diameter_tplus1)
sd_grad_diam<-sd(graduates$Diameter_tplus1)


p.vec_E[34]<-mu_grad_diam
p.vec_E[35]<-sd_grad_diam

p.vec_H[34]<-mu_grad_diam
p.vec_H[35]<-sd_grad_diam

p.vec_W[34]<-mu_grad_diam
p.vec_W[35]<-sd_grad_diam

summary(graduates$Height_tplus1)
hist(graduates$Height_tplus1)
mu_grad_height<-mean(graduates$Height_tplus1)
sd_grad_height<-sd(graduates$Height_tplus1)

p.vec_E[36]<-mu_grad_height
p.vec_E[37]<-sd_grad_height

p.vec_H[36]<-mu_grad_height
p.vec_H[37]<-sd_grad_height

p.vec_W[36]<-mu_grad_height
p.vec_W[37]<-sd_grad_height


p.vec_E[38]<-0.005
p.vec_E[39]<-0.002

p.vec_H[38]<-0.005
p.vec_H[39]<-0.002

p.vec_W[38]<-0.005
p.vec_W[39]<-0.002

save(p.vec_E, file="p.vec_E.RData")
save(p.vec_H, file="p.vec_H.RData")
save(p.vec_W, file="p.vec_W.RData")

save.image(file="170502_sites.RData")