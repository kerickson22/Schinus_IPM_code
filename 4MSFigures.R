# Code to produce figures in Erickson et al. integral projection manuscript. 
# Uses code provided by Ellner and Rees 
# Modifications to Ellner and Rees scripts added by Kelley D. Erickson


#Default width requested by journal: 
width.cm <- 12.9
width.cm_onepanel<-8.4
width.cm_onepanel_small<-3.9
height.cm <- 6.45
pointsize <- 8
x11(width = width.cm/2.54, height = height.cm/2.54, 
    pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.3,0.3,0.1,0.1))


m1=10
m2=m1+1
m3=100
m4=m3+1

h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
h3=(700-1.6)/m3;
y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
h4=(800-16)/m4
y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2



#setwd("/Users/kelley/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data/MSFigures")
#setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data/MSFigures")

#load('MSFigures_1.RData') 

###Map of Study Sites

# #sites <- read.csv('sites.csv', as.is=TRUE)
# 
# 
# sitesSpatial <- SpatialPointsDataFrame(
# 	coords= cbind(sites$Longitude, sites$Latitude),
# 	data=sites,
# 	proj4string=CRS('+proj=longlat +datum=WGS84 +no defs +ellps=WGS84 +towgs84=0,0,0'))
# 
# 
# states <- rgdal::readOGR(dsn='./USA_adm_shp', 
# layer='USA_adm2')
# 
# florida<-subset(florida, NAME_1=='Florida')
# 
# 
# 
# 
# 
# xs<-c(-87.63723, -78.5)
# ys<-c(31.00211, 24.52042)
# 
# 
# 
# png(file="map.png", width=13.5, height=10, units='in', res=50)
# par(ps=24, oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
# plot(xs, ys, pch=26)
# 
# sp::plot(florida, col='gray80',add=TRUE, xpd=NA)
# points(sitesSpatial, pch=c(15, 15, 16, 16, 17, 17), cex=4, bg=c('purple', 'purple', 'blue', 'blue', 'red', 'red'))
# 
# #points(-87.63723, 31.00211, pch=15, col="red", cex=4)
# #points(-78.5, 24.52042, pch=15, col="red", cex=4)
# legend(-87, 28, pch=c(15, 16, 17), pt.cex=4, y.intersp=3, x.intersp=2.5, legend=c("Hybrid", "Eastern", "Western"), bty='n')
# #text(sitesSpatial, labels=sites$Site)
# rect(-87.2, 26.5, -85.5, 28)
# 
# #xleft, ybottom, xright, ytop
# 
# 
# text(-82, 26, "Big Cypress")
# text(-82.8, 26.3, "Wild Turkey")
# text(-83.2, 26.9, "Punta Gorda")
# text(-79.6, 25.6, "Chekika")
# text(-79.5, 27.4, "Fort Pierce")
# text(-79.5, 28.4, "Cape Canaveral")
# dev.off()


x1seq_seedlings<-seq(0, 1.6, length.out=50)

x2seq_seedlings<-seq(0, 16, length.out=50)
x1seq_larges <- seq(0, 800, length.out = 50 )
x2seq_larges <- seq (0, 800, length.out = 50)


###FIGURE 1: D1 Survival by Biotype
load("./Overall/p.vec_overall.RData")
load("./Eastern/p.vec_E.RData")
load("./Hybrid/p.vec_H.RData")
load("./Western/p.vec_W.RData")

#png(file="D1_survival.png", width=20, height=10, units="in", res=300)
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D1_survival.eps", width=width.cm/2.54, 
           height=height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
 #x11(width = width.cm/2.54, height = height.cm/2.54, 
  #   pointsize = pointsize)
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(1, 2))
b0<-p.vec_E[20]
b1<-p.vec_E[21]
b2<-p.vec_E[22]


z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

persp(x1seq_seedling, x2seq_seedling, z1, ticktype="detailed",
      theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)",
      shade=0.1, nticks=4, xlab="\n Diameter (mm)",
      ylab="\n Height (cm)",  lwd=1, xaxs="i",
      main="(a) Eastern and Western", cex.main=1)

###Plotting Hybrid seedling survival
x1seq_seedling<-seq(0, 1.6, length.out=50)
x2seq_seedling<-seq(0, 16, length.out=50)

b0<-p.vec_H[20]
b1<-p.vec_H[21]
b2<-p.vec_H[22]


z1<-outer(x1seq_seedling, x2seq_seedling, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

persp(x1seq_seedling, x2seq_seedling, z1, ticktype="detailed",
      theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", 
      shade=0.1, nticks=4, xlab="\n Diameter (mm)", 
      ylab="\n Height (cm)",  lwd=1, xaxs="i", 
      main="(b) Hybrid", cex.main=1)

#This creates an EPS figure! That LaTeX can read

dev.off()





### FIGURE 2 D1 GROWTH 
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D1_growth.eps", width=width.cm/2.54, 
           height=height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(1, 2))
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



persp(x1seq_seedling, x2seq_seedling, z_diam, ticktype="detailed",
      theta=-30,  zlab="\n \n Diameter \n at t plus 1 (mm) ", 
      shade=0.1, nticks=4, xlab="\n Diameter at t (mm)",
      ylab="\n Height at t (cm)", 
      lwd=1)


persp(x1seq_seedling, x2seq_seedling, z_height, ticktype="detailed",
      theta=-30,  zlab="\n \n Height \n at t plus 1 (cm) ",
      shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", 
      ylab="\n Height at t (cm)", lwd=1)

dev.off()

###Figure 3: D2 Survival
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D2_survival.eps", width=width.cm/2.54, 
           height=height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(1, 2))
b0<-p.vec_E[1]
b1<-p.vec_E[2]
b2<-p.vec_E[3]

z1<-outer(x1seq_larges, x2seq_larges, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

persp(x1seq_larges, x2seq_larges, z1, ticktype="detailed",
      theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", 
      shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)",
      lwd=1, main="(a) Eastern and Western", cex.main=1) #xaxs="i"
#Plotting Hybrid seedling survival
b0<-p.vec_H[1]
b1<-p.vec_H[2]
b2<-p.vec_H[3]

z1<-outer(x1seq_larges, x2seq_larges, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

persp(x1seq_larges, x2seq_larges, z1, ticktype="detailed", 
      theta=-30, zlim=c(0, 1.25), zlab="\n P(Survival)", 
      shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)",  
      lwd=1, main="(b) Hybrid", cex.main=1) #xaxs="i"

dev.off()


###FIGURE 4: GROWTH OF D2 INDIVIDUALS

### Plot Growth Mod Eastern Larges

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D2_growth.eps", width=width.cm/2.54, 
           height=3*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(3, 2))
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

persp(x1seq_larges, x2seq_larges, z_diam, ticktype="detailed", 
      theta=-30,  zlab="\n \n Diameter \n at t plus 1 (mm) ", 
      shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)",
      lwd=1, cex.lab=1.5)
mtext(" (a) Eastern" , side=2, line = 4)

persp(x1seq_larges, x2seq_larges, z_height, ticktype="detailed",
      theta=-30,  zlab="\n \n \n Height \n at t plus 1 (cm) ", 
      shade=0.1, nticks=4, xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)",
      lwd=1,  cex.lab=1.5)

#Plot Growth Mod Hybrid Larges

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

persp(x1seq_larges, x2seq_larges, z_diam, ticktype="detailed", theta=-30,
      zlab="\n \n Diameter \n at t plus 1 (mm) ", shade=0.1, nticks=4,
      xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)", 
      cex.lab=1.5,  lwd=1)
mtext(" (b) Hybrid" , side=2, line = 4)

persp(x1seq_larges, x2seq_larges, z_height, ticktype="detailed", theta=-30,  
      zlab="\n \n \n Height \n at t plus 1 (cm) ", shade=0.1, nticks=4,
      xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)",
      cex.lab=1.5,  lwd=1)

#Plot Growth Mod Western Larges

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

persp(x1seq_larges, x2seq_larges, z_diam, ticktype="detailed", theta=-30,
      zlab="\n \n Diameter \n at t plus 1 (mm) ", shade=0.1, nticks=4,
      xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)",
      cex.lab=1.5,  lwd=1)
mtext(" (c) Western" , side=2, line = 4)

persp(x1seq_larges, x2seq_larges, z_height, ticktype="detailed", theta=-30, 
      zlab="\n \n \n Height \n at t plus 1 (cm) ", shade=0.1, nticks=4, 
      xlab="\n Diameter at t (mm)", ylab="\n Height at t (cm)",
      cex.lab=1.5,  lwd=1)
dev.off()


###FIGURE 5: GRADUATION
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/maturation.eps", width=width.cm/2.54, 
           height=height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(1, 2))
#Eastern and Hybrid
b0<-p.vec_E[31]
b1<-p.vec_E[32]
b2<-p.vec_E[33]



z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

persp(x1seq_seedlings, x2seq_seedlings, z1, ticktype="detailed",
      theta=-30, zlim=c(0, 1.25), zlab="\n P(Maturation)", 
      shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)", 
      main="(a) Eastern and Hybrid", lwd=1)

b0<-p.vec_W[31]
b1<-p.vec_W[32]
b2<-p.vec_W[33]

z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

persp(x1seq_seedlings, x2seq_seedlings, z1, ticktype="detailed",
      theta=-30, zlim=c(0, 1.25), zlab="\n P(Maturation)",
      shade=0.1, nticks=4, xlab="\n Diameter (mm)", ylab="\n Height (cm)",
      main=" (b) Western", lwd=1)

dev.off()

###FIGURE 6: REPRODUCTION
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/reproduction.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
 # x11(width = width.cm_onepanel/2.54, height = width.cm_onepanel/2.54, 
 #     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))
b0<-p.vec_E[12]
b1<-p.vec_E[13]
y_E<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))

b0<-p.vec_H[12]
b1<-p.vec_H[13]
y_H<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))

b0<-p.vec_W[12]
b1<-p.vec_W[13]
y_W<-((exp(b0+b1*x2seq_larges)/(1+exp(b0+b1*x2seq_larges))))

plot(x2seq_larges, y_E, xlab="Height at time t (cm)", 
     ylab="\n \n Probability of Reproducing", lwd=1, type='l', lty=1)
lines(x2seq_larges, y_W,  lwd=1, type='l', lty=3)
legend(5, 0.98, c("Eastern and Hybrid", "Western"), 
       lty=c(1, 3), lwd=c(1, 1), seg.len=1.5)
dev.off()

###FIGURE 7: FECUNDITY
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/fecundity.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))

y_overall<- 4.89*x1seq_larges*x1seq_larges
y_E<- 4.00*x1seq_larges*x1seq_larges
y_H<- 4.25*x1seq_larges*x1seq_larges
y_W<-5.52*x1seq_larges*x1seq_larges


plot(x1seq_larges, y_E, xlab="Diameter at time t (mm)", 
     ylab="\n \n Number of Offspring",  lwd=1, type='l', lty=5)
lines(x1seq_larges, y_H, lwd=1, type='l', lty=1)
lines(x1seq_larges, y_W,  lwd=1, type='l', lty=3)
#lines(x1seq_larges, y_overall, lwd=3, type='l', lty=2, col="red")
legend(5, 2500000, c("Eastern", "Hybrid", "Western"),   lty=c(5, 1, 3, 2), seg.len=3)
dev.off()




save.image(file="2018_06_27.RData")







#rm(list=ls())
load("./Eastern/total.contrib_D1_E.RData")
load("./Eastern/total.contrib_D2_E.RData")
load("./Eastern/total.contrib_F_E.RData")
load("./Eastern/total.contrib_G_E.RData")

load("y1.RData")
load("y2.RData")
load("y3.RData")
load("y4.RData")

load("total.contrib_D1_H.RData")
load("total.contrib_D2_H.RData")
load("total.contrib_F_H.RData")
load("total.contrib_G_H.RData")

load("total.contrib_D1_W.RData")
load("total.contrib_D2_W.RData")
load("total.contrib_F_W.RData")
load("total.contrib_G_W.RData")

load("lambda_diff_E.RData")
load("lambda_diff_H.RData")
load("lambda_diff_W.RData")

load("./Eastern/E.RData")
load("./Hybrid/H.RData")
load("./Western/W.RData")

load("lambda_all.RData")
load("lambda_E.RData")
load("lambda_H.RData")
load("lambda_W.RData")

#thing<-A_cv*elas_all
load("total.thing_D1.RData")
load("total.thing_D2.RData")
load("total.thing_F.RData")
load("total.thing_G.RData")

#Figure 8: Population Growth Rates
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/lambdas.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))
barplot(c(lambda_all, lambda_E, lambda_H, lambda_W), 
        col=c("grey", "grey", "grey", "grey"),
        ylab=expression(paste(lambda)), 
        names.arg=c("Overall", "Eastern", "Hybrid", "Western"),
        ylim=c(1, 1.1), xpd=F)
dev.off()



#Figure 9: Barplot of contribution of kernel components to differences in lambda
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/barplot_contributions.eps", width=width.cm_onepanel/2.54, 
           height=3*height.cm/2.54, pointsize=pointsize)
#png(file="./Figures/barplot_contributions.png")
 #x11(width = width.cm_onepanel/2.54, height = 3*width.cm_onepanel/2.54, 
  #   pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(0, 0, 0), # Distance of axis tickmark labels (second value)
    tcl = 0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.2,.5,0.2,0.1),
    oma=c(0.1, 0, 0.1, 0),
    mfrow=c(3, 1))

barplot(W, col=c("grey", "grey", "grey"), ylim=c(-0.03, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), axes=T, cex.names=1.6, cex.axis=1.6)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["W"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(H, col=c("grey", "grey", "grey"), ylim=c(-0.03, 0.06), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=1.6, cex.axis=1.6)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["H"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(E, col=c("grey", "grey", "grey"), ylim=c(-0.03, 0.06), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=1.6, cex.axis=1.6)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["E"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)
dev.off()

#Figure A1: SSD and RV 

# #m1=10
# m2=m1+1
# m3=100
# m4=m3+1

# h1=1.6/m1; 
# y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
# h2=16/m2
# y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
# h3=(700-1.6)/m3;
# y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
# h4=(800-16)/m4
 # y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2



# setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Eastern")
# load("v.RData")
# v_E<-v
# rm(v)
# load("stable.dist.RData")
# stable.dist_E<-stable.dist
# rm(stable.dist)

# setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Hybrid")
# load("v.RData")
# v_H<-v
# rm(v)
# load("stable.dist.RData")
# stable.dist_H<-stable.dist
# rm(stable.dist)


# setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Western")
# load("v.RData")
# v_W<-v
# rm(v)
# load("stable.dist.RData")
# stable.dist_W<-stable.dist
# rm(stable.dist)






# repro.val_smalls_E=matrix(v_E[1:(m1*m2)], m1,m2);
# repro.val_smalls_H=matrix(v_H[1:(m1*m2)], m1,m2);
# repro.val_smalls_W=matrix(v_W[1:(m1*m2)], m1,m2);




# stable.state_smalls_E=matrix(stable.dist_E[1:(m1*m2)], m1, m2);
# stable.state_smalls_H=matrix(stable.dist_H[1:(m1*m2)], m1, m2);
# stable.state_smalls_W=matrix(stable.dist_W[1:(m1*m2)], m1, m2);

# repro.val_larges_E=matrix(v_E[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
# repro.val_larges_H=matrix(v_H[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
# repro.val_larges_W=matrix(v_W[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);

# stable.state_larges_E=matrix(stable.dist_E[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);
# stable.state_larges_H=matrix(stable.dist_H[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);
# stable.state_larges_W=matrix(stable.dist_W[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);


# ss_diam1_E<-apply(stable.state_smalls_E, 1, sum)
# ss_diam1_H<-apply(stable.state_smalls_H, 1, sum)
# ss_diam1_W<-apply(stable.state_smalls_W, 1, sum)

# ss_diam2_E<-apply(stable.state_larges_E, 1, sum)
# ss_diam2_H<-apply(stable.state_larges_H, 1, sum)
# ss_diam2_W<-apply(stable.state_larges_W, 1, sum)


# ss_diam_E<-c(ss_diam1_E, ss_diam2_E)
# ss_diam_H<-c(ss_diam1_H, ss_diam2_H)
# ss_diam_W<-c(ss_diam1_W, ss_diam2_W)



# y_diam<-c(y1, y3)

# ss_height1_E<-apply(stable.state_smalls_E, 2, sum)
# ss_height1_H<-apply(stable.state_smalls_H, 2, sum)
# ss_height1_W<-apply(stable.state_smalls_W, 2, sum)



# ss_height2_E<-apply(stable.state_larges_E, 2, sum)
# ss_height2_H<-apply(stable.state_larges_H, 2, sum)
# ss_height2_W<-apply(stable.state_larges_W, 2, sum)

# ss_height_E<-c(ss_height1_E, ss_height2_E)
# ss_height_H<-c(ss_height1_H, ss_height2_H)
# ss_height_W<-c(ss_height1_W, ss_height2_W)


# y_height<-c(y2, y4)


# #CONVERT BELOW FOR ALL BIOTYPES

# rv_diam1_E<-apply(repro.val_smalls_E, 1, sum)
# rv_diam1_H<-apply(repro.val_smalls_H, 1, sum)
# rv_diam1_W<-apply(repro.val_smalls_W, 1, sum)


# rv_diam2_E<-apply(repro.val_larges_E, 1, sum)
# rv_diam2_H<-apply(repro.val_larges_H, 1, sum)
# rv_diam2_W<-apply(repro.val_larges_W, 1, sum)

# rv_diam_E<-c(rv_diam1_E, rv_diam2_E)
# rv_diam_H<-c(rv_diam1_H, rv_diam2_H)
# rv_diam_W<-c(rv_diam1_W, rv_diam2_W)

# y_diam<-c(y1, y3)


# rv_height1_E<-apply(repro.val_smalls_E, 2, sum)
# rv_height1_H<-apply(repro.val_smalls_H, 2, sum)
# rv_height1_W<-apply(repro.val_smalls_W, 2, sum)


# rv_height2_E<-apply(repro.val_larges_E, 2, sum)
# rv_height2_H<-apply(repro.val_larges_H, 2, sum)
# rv_height2_W<-apply(repro.val_larges_W, 2, sum)

# rv_height_E<-c(rv_height1_E, rv_height2_E)
# rv_height_H<-c(rv_height1_H, rv_height2_H)
# rv_height_W<-c(rv_height1_W, rv_height2_W)


# y_height<-c(y2, y4)

# setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data/MSFigures")

# save(y_diam,y_height, ss_diam_E, ss_diam_H, ss_diam_W, y_height, ss_height_E, ss_height_H, ss_height_W, rv_diam_E, rv_diam_H, rv_diam_W, rv_height_E, rv_height_H, rv_height_W, file="ssd_rv_figure.RData")

load("ssd_rv_figure.RData")


setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/ssd_rv.eps", width=width.cm/2.54, 
           height=2*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
#png(file="./Figures/ssd_rv.png")
  #x11(width = width.cm/2.54, height = 2*height.cm/2.54, 
   #   pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.4,0.3,0.2,0.1),
    mfrow=c(2, 2))

load("./Overall/y_diam.RData")
load("./Overall/y_height.RData")
load("./Eastern/ss_diam_E.RData")
load("./Hybrid/ss_diam_H.RData")
load("./Western/ss_diam_W.RData")
load("./Eastern/ss_height_E.RData")
load("./Hybrid/ss_height_H.RData")
load("./Western/ss_height_W.RData")
load("./Eastern/rv_diam_E.RData")
load("./Eastern/rv_height_E.RData")
load("./Hybrid/rv_diam_H.RData")
load("./Hybrid/rv_height_H.RData")
load("./Western/rv_diam_W.RData")
load("./Western/rv_height_W.RData")

plot(y_diam, ss_diam_E,xlab="Diameter (mm)",ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=5, lwd=1, ylim=c(0, 0.2)); 
lines(y_diam, ss_diam_W, lwd=1, lty=3)
lines(y_diam, ss_diam_H, lwd=1, lty=1)

title(main="(a) Stable diameter distribution", cex.main=1);


plot(y_height,ss_height_E,xlab="Height (cm)",ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=5, lwd=1, ylim=c(0, 0.2)); 
lines(y_height, ss_height_W, lwd=1, lty=3)
lines(y_height, ss_height_H, lwd=1, lty=1)
title(main="(b) Stable height distribution", cex.main=1); 


plot(y_diam, rv_diam_H, xlab="Diameter(mm)", ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=1, lwd=1, ylim=c(0, 0.02))
lines(y_diam, rv_diam_E, lwd=1, lty=5)
lines(y_diam, rv_diam_W, lwd=1, lty=3)
lines(y_diam, rv_diam_H, lwd=1, lty=1)
title(main="(c) Reproductive value: by diameter", cex.main=1);

 
plot(y_height, rv_height_H, xlab="Height (cm)", ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=1, lwd=1, ylim=c(0, 0.02))
lines(y_height, rv_height_E, lwd=1, lty=5)
lines(y_height, rv_height_W, lwd=1, lty=3)
title(main="(d) Reproductive value: by height", cex.main=1); 

legend(350, 0.005, lty=c(3, 1, 5), lwd=1,
       c("Western", "Hybrid", "Eastern"), cex=1,  seg.len=3)

dev.off()




#Figure A2: Marginal elasticity by diameter
#setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Eastern")
# load("total.elas_D1_E.RData")
# load("total.elas_F_E.RData")
# load("total.elas_G_E.RData")
# load("total.elas_D2_E.RData")

load("./Eastern/total.elas.RData")

# #setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Hybrid")
# load("total.elas_D1_H.RData")
# load("total.elas_F_H.RData")
# load("total.elas_G_H.RData")
# load("total.elas_D2_H.RData")

load("./Hybrid/total.elas.RData")

# #setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Western")
# load("total.elas_D1_W.RData")
# load("total.elas_F_W.RData")
# load("total.elas_G_W.RData")
# load("total.elas_D2_W.RData")

load("./Western/total.elas.RData")

#setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data/MSFigures")

#### Figure Marginal elasticities 
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/marginal_elasticity.eps", width=width.cm/2.54, 
           height=1.5*width.cm/(2.54), pointsize=pointsize,  encoding = "TeXtext.enc")
 #png(file="./Figures/marginal_elasticity.png")
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, .5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.3,0.4,0.35,0.1),
    mfrow=c(4, 2))

x1_seedlings <- seq(0, 1.6, length.out=10)
x1_adults <- seq(1.6, 800, length.out=100)
x2_seedlings <- seq(0, 16, length.out=11)
x2_adults <- seq(16, 800, length.out=101)


plot(x1_seedlings, rowSums(total.elas_D1_E), ylim=c(0, 0.03), xlim=c(0, 1.6), xlab="", 
     ylab="", main="",  type='l', lwd=1, lty=5, cex.axis=1.5)
text(0.2, 0.049, "Seedling", cex=1.5)
lines(x1_seedlings, rowSums(total.elas_D1_H), type='l', lwd=1, lty=1)
lines(x1_seedlings, rowSums(total.elas_D1_W), type='l', lwd=1, lty=3)
mtext(side=3, "(a) Marginal elasticity over diameter", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.elas_F_E), ylim=c(0, 0.03), xlim=c(1.6, 800), xlab="", 
     ylab= "", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
text(65, 0.049, "Fertility", cex=1.5)
lines(x1_adults, rowSums(total.elas_F_H), type='l', lwd=1, lty=1)
lines(x1_adults, rowSums(total.elas_F_W), type='l', lwd=1, lty=3)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_seedlings, rowSums(total.elas_G_E), ylim= c(0, 0.03), xlim=c(0, 1.6), xlab="",
     ylab="", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
text(0.2, 0.049, "Maturation", cex=1.5)
lines(x1_seedlings, rowSums(total.elas_G_H), type='l', lty=1)
lines(x1_seedlings, rowSums(total.elas_G_W), type='l', lty=3)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.elas_D2_E), ylim= c(0, 0.03), xlim=c(1.6, 800), xlab=" ",
     ylab="", main="",  type='l',  lwd=1, lty=5, cex.axis=1.5)
text(65, 0.049, "Adult", cex=1.5)
lines(x1_adults, rowSums(total.elas_D2_H), type='l', lwd=1, lty=1)
lines(x1_adults, rowSums(total.elas_D2_W), type='l', lwd=1, lty=3)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x2_seedlings, colSums(total.elas_D1_E), ylim=c(0, 0.03), xlim=c(0, 16), xlab=" ", 
     ylab="", main="",  type='l', lwd=1, lty=5, cex.axis=1.5)
text(2, 0.049, "Seedling", cex=1.5)
lines(x2_seedlings, colSums(total.elas_D1_H), type='l', lwd=1, lty=1)
lines(x2_seedlings, colSums(total.elas_D1_W), type='l', lwd=1, lty=3)
mtext(side=3, "(b) Marginal elasticity over height", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.elas_F_E), ylim=c(0, 0.03), xlim=c(16,800),  xlab=" ",
     ylab= "", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
text(65, 0.049, "Fertility", cex=1.5)
lines(x2_adults, colSums(total.elas_F_H), type='l', lwd=1, lty=1)
lines(x2_adults, colSums(total.elas_F_W), type='l', lwd=1, lty=3)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)


plot(x2_seedlings, colSums(total.elas_G_E), ylim= c(0, 0.03), xlim=c(0, 16), xlab="",
     ylab="", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
lines(x2_seedlings, colSums(total.elas_G_H), type='l', lwd=1, lty=1)
lines(x2_seedlings, colSums(total.elas_G_W), type='l', lwd=1, lty=3)
text(2, 0.049, "Maturation", cex=1.5)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.elas_D2_E), ylim= c(0, 0.03), xlim=c(16, 800),  xlab=" ", 
     ylab="", main="",  type='l',  lwd=1, lty=5, cex.axis=1.5)
lines(x2_adults, colSums(total.elas_D2_H), type='l', lwd=1, lty=1)
lines(x2_adults, colSums(total.elas_D2_W), type='l', lwd=1, lty=3)
text(65, 0.049, "Adult", cex=1.5)
legend(250, 0.03, lty=c(3, 1, 5), lwd=1, c("Western", "Hybrid", "Eastern"),
       seg.len=2)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)
dev.off() #x11(width = width.cm/2.54, height = 1.5*width.cm/(2.54), 
  #   pointsize = pointsize)




###FIGURE B1: Contribution of matrix elements to differences in lambda by diameter
load("./Eastern/total.contrib_E.RData")
load("./Hybrid/total.contrib_H.RData")
load("./Western/total.contrib_W.RData")

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/marginal_contributions.eps", width=width.cm/2.54, 
           height=1.5*width.cm/(2.54), pointsize=pointsize,  encoding = "TeXtext.enc")
  #x11(width = width.cm/2.54, height = 1.5*width.cm/(2.54), 
   #  pointsize = pointsize)
#png(file="./Figures/marginal_contributions.png")
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, .5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.3,0.4,0.35,0.1),
    mfrow=c(4, 2))


plot(x1_seedlings, rowSums(total.contrib_D1_E), ylim=c(-0.003, 0.01),
     xlim=c(0, 1.6), xlab="", ylab="", main="",  type='l', lwd=1, 
     lty=5, cex.axis=1.5)
text(0.2, 0.01, "Seedling", cex=1.5)
lines(x1_seedlings, rowSums(total.contrib_D1_H), type='l', lwd=1, lty=1)
lines(x1_seedlings, rowSums(total.contrib_D1_W), type='l', lwd=1, lty=3)
mtext(side=3, "(a) Contributions by diameter", line=1, cex=1.2, adj=0)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.contrib_F_E), ylim=c(-0.003, 0.01), xlim=c(1.6, 800), xlab="", 
     ylab= "", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
text(65, 0.01, "Fertility", cex=1.5)
lines(x1_adults, rowSums(total.contrib_F_H), type='l', lwd=1, lty=1)
lines(x1_adults, rowSums(total.contrib_F_W), type='l', lwd=1, lty=3)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x1_seedlings, rowSums(total.contrib_G_E), ylim=c(-0.003, 0.01), xlim=c(0, 1.6), xlab="",
     ylab="", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
text(0.2, 0.01, "Maturation", cex=1.5)
lines(x1_seedlings, rowSums(total.contrib_G_H), type='l', lty=1)
lines(x1_seedlings, rowSums(total.contrib_G_W), type='l', lty=3)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.contrib_D2_E), ylim=c(-0.003, 0.01), xlim=c(1.6, 800), xlab=" ",
     ylab="", main="",  type='l',  lwd=1, lty=5, cex.axis=1.5)
text(65, 0.01, "Adult", cex=1.5)
lines(x1_adults, rowSums(total.contrib_D2_H), type='l', lwd=1, lty=1)
lines(x1_adults, rowSums(total.contrib_D2_W), type='l', lwd=1, lty=3)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x2_seedlings, colSums(total.contrib_D1_E), ylim=c(-0.003, 0.01), xlim=c(0, 16), xlab=" ", 
     ylab="", main="",  type='l', lwd=1, lty=5, cex.axis=1.5)
text(2, 0.01, "Seedling", cex=1.5)
lines(x2_seedlings, colSums(total.contrib_D1_H), type='l', lwd=1, lty=1)
lines(x2_seedlings, colSums(total.contrib_D1_W), type='l', lwd=1, lty=3)
mtext(side=3, "(b) Contributions by height", line=1, cex=1.2, adj=0)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.contrib_F_E), ylim=c(-0.003, 0.01), xlim=c(16,800),  xlab=" ",
     ylab= "", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
text(65, 0.01, "Fertility", cex=1.5)
lines(x2_adults, colSums(total.contrib_F_H), type='l', lwd=1, lty=1)
lines(x2_adults, colSums(total.contrib_F_W), type='l', lwd=1, lty=3)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)


plot(x2_seedlings, colSums(total.contrib_G_E), ylim=c(-0.003, 0.01), xlim=c(0, 16), xlab="",
     ylab="", main="", type='l', lwd=1, lty=5, cex.axis=1.5)
lines(x2_seedlings, colSums(total.contrib_G_H), type='l', lwd=1, lty=1)
lines(x2_seedlings, colSums(total.contrib_G_W), type='l', lwd=1, lty=3)
text(2, 0.01, "Maturation", cex=1.5)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.contrib_D2_E), ylim=c(-0.003, 0.01), xlim=c(16, 800),  xlab=" ", 
     ylab="", main="",  type='l',  lwd=1, lty=5, cex.axis=1.5)
lines(x2_adults, colSums(total.contrib_D2_H), type='l', lwd=1, lty=1)
lines(x2_adults, colSums(total.contrib_D2_W), type='l', lwd=1, lty=3)
text(65, 0.01, "Adult", cex=1.5)
legend(250, 0.01, lty=c(3, 1, 5), lwd=1, c("Western", "Hybrid", "Eastern"),
       seg.len=2)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)
dev.off()



#Figure B3: Contributions of matrix elements to variation in lambda by diameter
load("./Overall/total.thing.RData")

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/variability.eps", width=width.cm/2.54, 
           height=1.5*width.cm/(2.54), pointsize=pointsize,  encoding = "TeXtext.enc")
  #x11(width = width.cm/2.54, height = 1.5*width.cm/(2.54), 
   # pointsize = pointsize)
#png(file="./Figures/variability.png")
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, .5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.3,0.4,0.35,0.1),
    mfrow=c(4, 2))
plot(x1_seedlings, rowSums(total.thing_D1), ylim=c(0, 0.006),
     xlim=c(0, 1.6), xlab="", ylab="", main="",  type='l', lwd=1, 
     lty=1, cex.axis=1.5)
text(0.2, 0.0055, "Seedling", cex=1.5)
mtext(side=3, "(a) Contributions by diameter", line=1, cex=1.2, adj=0)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.thing_F), ylim=c(0, 0.006), xlim=c(1.6, 800), xlab="", 
     ylab= "", main="", type='l', lwd=1, lty=1, cex.axis=1.5)
text(65, 0.0055, "Fertility", cex=1.5)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x1_seedlings, rowSums(total.thing_G), ylim= c(0, 0.006), xlim=c(0, 1.6), xlab="",
     ylab="", main="", type='l', lwd=1, lty=1, cex.axis=1.5)
text(0.2, 0.0055, "Maturation", cex=1.5)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.thing_D2), ylim= c(0, 0.006), xlim=c(1.6, 800), xlab=" ",
     ylab="", main="",  type='l',  lwd=1, lty=1, cex.axis=1.5)
text(65, 0.0055, "Adult", cex=1.5)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x2_seedlings, colSums(total.thing_D1), ylim=c(0, 0.006), xlim=c(0, 16), xlab=" ", 
     ylab="", main="",  type='l', lwd=1, lty=1, cex.axis=1.5)
text(2, 0.0055, "Seedling", cex=1.5)
mtext(side=3, "(b) Contributions by height", line=1, cex=1.2, adj=0)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.thing_F), ylim=c(0, 0.006), xlim=c(16,800),  xlab=" ",
     ylab= "", main="", type='l', lwd=1, lty=1, cex.axis=1.5)
text(65, 0.0055, "Fertility", cex=1.5)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)


plot(x2_seedlings, colSums(total.thing_G), ylim=c(0, 0.006), xlim=c(0, 16), xlab="",
     ylab="", main="", type='l', lwd=1, lty=1, cex.axis=1.5)
text(2, 0.0055, "Maturation", cex=1.5)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.thing_D2), ylim=c(0, 0.006), xlim=c(16, 800),  xlab=" ", 
     ylab="", main="",  type='l',  lwd=1, lty=1, cex.axis=1.5)
text(65, 0.0055, "Adult", cex=1.5)
mtext(side=2, "CV * elas", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)
dev.off()


##### Original versions of 
