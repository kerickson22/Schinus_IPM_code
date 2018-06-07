#Code for running the demographic analyses incorporating biotype differences 
#Code written by Kelley D. Erickson and Carol C. Horvitz

library(RColorBrewer)
rm(list=ls(all=TRUE))
#MAC
setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 

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

setwd("/Users/curculion/Documents/GitHub/Schinus_IPM_code")
load("demography_overall.RData")
setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 





pooled2<-cbind(pooled, Diam_cat_t, Height_cat_t, Min_cat_t)




#Initialize p.vec's to store output from model (these will feed into the IPM)
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


#(i) Probability of being a reproductive female at time t: 


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






png(file="repro_by_site.png", width=10, height=10, units="in", res=300)
par(ps=24)
par(mar=c(5, 5, 4, 2))
plot(x2seq_larges, y_E, xlab="Height at time t (cm)", ylab="\n \n Probability of Reproducing", col=col_E_dark, lwd=3, type='l')
lines(x2seq_larges, y_H, col=col_H_dark, lwd=3, type='l')
lines(x2seq_larges, y_W, col=col_W_dark, lwd=3, type='l')
legend(600, 0.3, c("Eastern", "Hybrid", "Western"), col=c(col_E_dark, col_H_dark, col_W_dark), lty=c(1, 1, 1), y.intersp=2, lwd=c(3, 3, 3))
dev.off()






#(ii) Given that an individual is a reproductive female of size[diameter=x, height=y] and genetic type z[known from site], how many fruits do they produce? 
#Depends on diameter and genetic type only (from biomass allocation chapter)

mod_fruit_production<-lm(LHS$Fruit_No ~ 0 + LHS$Height:LHS$Genetic_Type)
summary(mod_fruit_production)

p.vec_E[14]<-mod_fruit_production$coeff[1]
p.vec_H[14]<-mod_fruit_production$coeff[2]
p.vec_W[14]<-mod_fruit_production$coeff[3]


y_E<-x1seq_larges*x1seq_larges*p.vec_E[14]
y_W<-x1seq_larges*x1seq_larges*p.vec_W[14]
y_H<-x1seq_larges*x1seq_larges*p.vec_H[14]

setwd("/Users/curculion/Dropbox/UM_Dissertation_LaTeX-master/untitled folder/Figures/biotype")

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







seed_surv<-.5072

tau_2_hybrid<-0.5072
tau_2_E-W<-0.388

p.vec_E[19]<-0.388
p.vec_H[19]<-0.5072
p.vec_W[19]<-0.388



p.vec_E[15]<-mu_diam
p.vec_E[16]<-sd_diam

p.vec_H[15]<-mu_diam
p.vec_H[16]<-sd_diam

p.vec_W[15]<-mu_diam
p.vec_W[16]<-sd_diam



p.vec_E[17]<-mu_height
p.vec_E[18]<-sd_height

p.vec_H[17]<-mu_height
p.vec_H[18]<-sd_height

p.vec_W[17]<-mu_height
p.vec_W[18]<-sd_height

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


p.vec_E[34]<-mu_grad_diam
p.vec_E[35]<-sd_grad_diam

p.vec_H[34]<-mu_grad_diam
p.vec_H[35]<-sd_grad_diam

p.vec_W[34]<-mu_grad_diam
p.vec_W[35]<-sd_grad_diam

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

save.image(file="demography_sites.RData")