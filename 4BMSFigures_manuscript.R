# Code to produce figures in Erickson et al. integral projection manuscript. 
# Uses code provided by Ellner and Rees 
# Modifications to Ellner and Rees scripts added by Kelley D. Erickson

library(RColorBrewer)
library(viridis)
library(ggplot2)
library(arm)
library(gridExtra)
library(grid)
setwd("~/Documents/GitHub/Schinus_IPM_code")

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

cols <- viridis(6)
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


seedlings2<-subset(seedlings, seedlings$Diameter_tplus1<1.6)
seedlings2<-subset(seedlings2, seedlings2$Height_tplus1<16)

display.brewer.all(colorblindFriendly=TRUE)
my.palette <- brewer.pal(11, "RdYlBu")
my.palette2 <- brewer.pal(9, "YlGnBu")
#setwd("/Users/kelley/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data/MSFigures")
#setwd("/Users/curculion/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data/MSFigures")

#load('MSFigures_1.RData') 

###Map of Study Sites
#FIG 1: Map of Study Sites #####
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

#FIG 2: Constructed using tikz picture in document #####

x1seq_seedlings<-seq(0, 1.6, length.out=50)

x2seq_seedlings<-seq(0, 16, length.out=50)
x1seq_larges <- seq(1.6, 800, length.out = 50 )
x2seq_larges <- seq (16, 800, length.out = 50)


# FIG 3: D1 Survival by Biotype #####
load("./BC/p.vec_BC.RData")
load("./CC/p.vec_CC.RData")
load("./C/p.vec_C.RData")
load("./FP/p.vec_FP.RData")
load("./PG/p.vec_PG.RData")
load("./WT/p.vec_WT.RData")

df <- expand.grid(x1seq_seedlings, x2seq_seedlings)
#Punta Gorda, Wild Turkey and Big Cypress
b0 <- p.vec_BC[1]
b1 <- p.vec_BC[2]
b2 <- p.vec_BC[3]
df$BC_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_BCs<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= BC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(a) PG, WT, BC") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Big Cypress" |
                              seedlings$Site == "Punta Gorda" |
                              seedlings$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t, shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


#Chekika
b0 <- p.vec_C[1]
b1 <- p.vec_C[2]
b2 <- p.vec_C[3]
df$C_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_Cs <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= C_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(b) C ") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t, shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Fort Pierce
b0 <- p.vec_FP[1]
b1 <- p.vec_FP[2]
b2 <- p.vec_FP[3]
df$FP_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_FPs <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= FP_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(c) FP") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Fort Pierce")),
             aes(x=Diameter_t, y=Height_t, shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Cape Canaveral
b0 <- p.vec_CC[1]
b1 <- p.vec_CC[2]
b2 <- p.vec_CC[3]
df$CC_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_CCs<-ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= CC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(d) CC") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Cape Canaveral")),
             aes(x=Diameter_t, y=Height_t, shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

g <- ggplotGrob(plot_BCs + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot_BCs2 <- plot_BCs + theme(legend.position ="none")
plot_CCs2 <- plot_CCs + theme(legend.position ="none")
plot_Cs2  <- plot_Cs  + theme(legend.position = "none")
plot_FPs2 <- plot_FPs + theme(legend.position = "none")

postscript("./Figures/D1_survival2.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 5, 5),
             c(1, 1, 2, 2, 5, 5),
             c(3, 3, 4, 4, 5, 5),
             c(3, 3, 4, 4, 5, 5))
grid.arrange(plot_BCs2, plot_Cs2,plot_FPs2, plot_CCs2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()

###FIG 4: MATURATION #####

#PG, BC, FP, CC
b0 <- p.vec_BC[12]
b1 <- p.vec_BC[13]
b2 <- p.vec_BC[14]
df$BC_m <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

#Big Cypress, Cape Canaveral, Fort Pierce, Punta Gorda
plot_BCm<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= BC_m)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(a) PG, BC, FP, CC") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$grad_status) &
                           (seedlings$Site == "Big Cypress" |
                              seedlings$Site == "Cape Canaveral" |
                              seedlings$Site == "Fort Pierce" |
                              seedlings$Site == "Punta Gorda")),
             aes(x=Diameter_t, y=Height_t, shape = factor(grad_status))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n maturation \n at time \n t + 1") +
  labs(shape = "Measured \n maturation \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


#WT
b0 <- p.vec_WT[12]
b1 <- p.vec_WT[13]
b2 <- p.vec_WT[14]
df$WT_m <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

#Wild Turkey
plot_WTm<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= WT_m)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(b) WT") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$grad_status) &
                           (seedlings$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t, shape = factor(grad_status))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n maturation \n at time \n t + 1") +
  labs(shape = "Measured \n maturation \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


#Chekika
b0 <- p.vec_C[12]
b1 <- p.vec_C[13]
b2 <- p.vec_C[14]
df$C_m <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

#Chekika
plot_Cm<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= C_m)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(c) C") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$grad_status) &
                           (seedlings$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t, shape = factor(grad_status))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n maturation \n at time \n t + 1") +
  labs(shape = "Measured \n maturation \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


g <- ggplotGrob(plot_BCm + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot_BCm2 <- plot_BCm + theme(legend.position ="none")
plot_Cm2<- plot_Cm + theme(legend.position ="none")
plot_WTm2  <- plot_WTm  + theme(legend.position = "none")




postscript("./Figures/Maturation.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2),
             c(1, 1, 2, 2),
             c(3, 3, 4, 4),
             c(3, 3, 4, 4))
grid.arrange(plot_BCm2, plot_WTm2,plot_Cm2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1))
dev.off()

###FIG 5: REPRODUCTION #####
#Punta Gorda
b0<-p.vec_PG[30]
b1<-p.vec_PG[31]
pf_PG <- invlogit(b0 + b1*x2seq_larges)

#Wild Turkey
b0<-p.vec_WT[30]
b1<-p.vec_WT[31]
pf_WT <- invlogit(b0 + b1*x2seq_larges)

#Big Cypress, Fort Pierce, Cape Canaveral
b0<-p.vec_BC[30]
b1<-p.vec_BC[31]
pf_BC <- invlogit(b0 + b1*x2seq_larges)

#Chekika
b0<-p.vec_C[30]
b1<-p.vec_C[31]
pf_C <- invlogit(b0 + b1*x2seq_larges)


setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/reproduction2.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm_onepanel/2.54, height = width.cm_onepanel/2.54, 
#     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))



plot(x2seq_larges, pf_PG, xlab="Height at time t (cm)", ylim=c(0,1),
     ylab="\n \n Probability of Reproducing", lwd=1, type='l', lty=1, col=cols[1])
points(larges$Height_t[larges$Site == "Big Cypress" | 
                         larges$Site == "Cape Canaveral" |
                         larges$Site == "Fort Pierce"], 
       larges$Rep_tplus1[larges$Site == "Big Cypress" | 
                           larges$Site == "Cape Canaveral" |
                           larges$Site == "Fort Pierce"],
       col=cols[1])
lines(x2seq_larges, pf_WT, lwd=2, type='l', lty=2, col=cols[2])
points(larges$Height_t[larges$Site == "Wild Turkey"], 
       larges$Rep_tplus1[larges$Site == "Wild Turkey"],
       col=cols[2])
lines(x2seq_larges, pf_BC, lwd=2, type='l', lty=3, col=cols[3])
points(larges$Height_t[larges$Site == "Big Cypress" |
                         larges$Site == "Fort Pierce" |
                         larges$Site == "Cape Canaveral"],
       col = cols[3])
lines(x2seq_larges, pf_C, lwd=2, type='l', lty=4, col=cols[4])
points(larges$Height_t[larges$Site == "Chekika"], 
       larges$Rep_tplus1[larges$Site == "Chekika"],
       col=cols[4])
legend(5, 0.98, c("PG", "WT", "BC-FP-CC", "C"), 
       lty=c(1, 2, 3, 4), lwd=2, seg.len=2, col=cols[c(1,2, 3, 4)])
dev.off()

###FIG 6: FECUNDITY #####

y_BC <- p.vec_BC[32]*x1seq_larges*x1seq_larges
y_CC <- p.vec_CC[32]*x1seq_larges*x1seq_larges
y_C  <- p.vec_C[32]*x1seq_larges*x1seq_larges
y_FP <- p.vec_FP[32]*x1seq_larges*x1seq_larges
y_PG <- p.vec_PG[32]*x1seq_larges*x1seq_larges
y_WT <- p.vec_WT[32]*x1seq_larges*x1seq_larges

x<-Stems$diam_base*10 #convert from cm to mm
#x <- x*x

Stems <- cbind(Stems, x)

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/fecundity2.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))

plot(x1seq_larges, y_PG, xlab="Diameter at time t (mm)",
     ylab="\n \n Number of Offspring",  lwd=2, type='l', lty=1, col=cols[1], ylim=c(0, 4000000))
points(Stems$x[Stems$Site == "Punta Gorda"], 
       Stems$Seed_No[Stems$Site == "Punta Gorda"], col=cols[1])
lines(x1seq_larges, y_WT, lwd=2, type='l', lty=2, col=cols[2])
points(Stems$x[Stems$Site == "Wild Turkey"],
       Stems$Seed_No[Stems$Site == "Wild Turkey"], col=cols[2])
lines(x1seq_larges, y_BC, lwd=2, type='l', lty=3, col=cols[3])
points(Stems$x[Stems$Site == "Big Cypress"],
       Stems$Seed_No[Stems$Site == "Big Cypress"], col=cols[3])
lines(x1seq_larges, y_C, lwd=2, type='l', lty=4, col=cols[4])
points(Stems$x[Stems$Site == "Chekika"],
       Stems$Seed_No[Stems$Site == "Chekika"], col=cols[4])
lines(x1seq_larges, y_FP, lwd=2, type='l', lty=5, col=cols[5])
points(Stems$x[Stems$Site == "Fort Pierce"],
       Stems$Seed_No[Stems$Site == "Fort Pierce"], col=cols[5])

lines(x1seq_larges, y_CC, lwd=2, type='l', lty=2, col=cols[6])
points(Stems$x[Stems$Site == "Cape Canaveral"],
       Stems$Seed_No[Stems$Site == "Cape Canaveral"], col=cols[6])
legend(5, 2500000, lwd=2,  c("PG", "WT", "BC", "C", "FP", "CC"),
       lty=1:6, seg.len=3, col=cols)
dev.off()

#FIG 7: Population Growth Rates #####
lam.stable_BC <- readRDS("./BC/lam.stable_BC.rds")
lam.stable_CC <- readRDS("./CC/lam.stable_CC.rds")
lam.stable_C  <- readRDS("./C/lam.stable_C.rds")
lam.stable_FP <- readRDS("./FP/lam.stable_FP.rds")
lam.stable_PG <- readRDS("./PG/lam.stable_PG.rds")
lam.stable_WT <- readRDS("./WT/lam.stable_WT.rds")

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
barCenters <- plot(c(lam.stable_PG, lam.stable_WT, lam.stable_BC,
                     lam.stable_C, lam.stable_FP, lam.stable_CC), 
                   col=cols,
                   cex=2,
                   ylab=expression(paste(lambda)),
                   axes= F,
                   xlab= "Site",
                   pch=19,
                   xpd=F, 
                   ylim =c(0.98, 1.26))
axis(side=1, labels = c("PG", "WT", "BC", "C",
                        "FP", "CC"), at = 1:6 )
axis(side=2)
abline(h=1, lty=2)
segments(3, 1.062, 3, 1.206) #Big Cypress
arrows(3, 1.062, 3, 1.206, angle = 90, 
       code = 3, length = 0.05)

segments(6, 1.055, 6, 1.142) #Cape Canaveral
arrows(6, 1.055, 6, 1.142, angle = 90, 
       code = 3, length = 0.05)

segments(4, 0.986, 4, 1.05) #Chekika
arrows(4, 0.986, 4, 1.05, angle = 90, 
       code = 3, length = 0.05)

segments(5, 1.079, 5, 1.221) #Fort Pierce
arrows(5, 1.079, 5, 1.221, angle = 90, code = 3, length = 0.05)

segments(1, 0.995, 1, 1.140) #Punta Gorda
arrows(1, 0.995, 1, 1.140, angle=90, code = 3, length = 0.05)

segments(2, 1.062, 2, 1.251) #Wild Turkey
arrows(2, 1.062, 2, 1.251, angle = 90, code =3, length = 0.05)
dev.off()

#FIG 8: Marginal elasticity by diameter ######

load("./BC/total.elas.RData")
load("./CC/total.elas.RData")
load("./C/total.elas.RData")
load("./FP/total.elas.RData")
load("./PG/total.elas.RData")
load("./WT/total.elas.RData")


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
x1_adults <- seq(1.6, 800, length.out=370)
x2_seedlings <- seq(0, 16, length.out=11)
x2_adults <- seq(16, 800, length.out=150)

#Marginal elasticity seedlings by diameter
plot(x1_seedlings, rowSums(total.elas_D1_PG), ylim=c(0, 0.03),
     xlim=c(0, 1.6), xlab="", ylab="", main="",  type='l',
     lwd=1, lty=1, cex.axis=1.5, col=cols[1])
text(0.2, 0.029, "Seedling", cex=1.5)
lines(x1_seedlings, rowSums(total.elas_D1_WT), type='l', lwd=1,
      lty=2, col=cols[2])
lines(x1_seedlings, rowSums(total.elas_D1_BC), type='l', lwd=1, 
      lty=3, col=cols[3])
lines(x1_seedlings, rowSums(total.elas_D1_C), type='l', lwd=1, 
      lty=4, col=cols[4])
lines(x1_seedlings, rowSums(total.elas_D1_FP), type='l', lwd=1,
      lty=5, col=cols[5])
lines(x1_seedlings, rowSums(total.elas_D1_CC), type='l', lwd=1, 
      lty=6, col=cols[6])
mtext(side=3, "(a)", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

#Marginal elasticity of seedlings by height
plot(x2_seedlings, colSums(total.elas_D1_PG), ylim=c(0, 0.03),
     xlim=c(0, 16), xlab=" ", ylab="", main="",  type='l', lwd=1,
     lty=1, cex.axis=1.5, col=cols[1])
text(2, 0.029, "Seedling", cex=1.5)
lines(x2_seedlings, colSums(total.elas_D1_WT), type='l', lwd=1,
      lty=2, col=cols[2])
lines(x2_seedlings, colSums(total.elas_D1_BC), type='l', lwd=1,
      lty=3, col=cols[3])
lines(x2_seedlings, colSums(total.elas_D1_C), type='l', lwd=1,
      lty=4, col=cols[4])
lines(x2_seedlings, colSums(total.elas_D1_FP), type='l', lwd=1,
      lty=5, col=cols[5])
lines(x2_seedlings, colSums(total.elas_D1_CC), type='l', lwd=1, 
      lty=6, col=cols[6])
mtext(side=3, "(b)", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

#marginal elasticity of maturation by diameter
plot(x1_seedlings, rowSums(total.elas_GA_PG,total.elas_GB_PG), 
     ylim= c(0, 0.03), xlim=c(0, 1.6), xlab="",
     ylab="", main="", type='l', lwd=1, lty=1, cex.axis=1.5, 
     col=cols[1])
text(0.2, 0.029, "Maturation", cex=1.5)
lines(x1_seedlings, rowSums(total.elas_GA_WT, total.elas_GB_WT),
      type='l', lty=2, col=cols[2])
lines(x1_seedlings, rowSums(total.elas_GA_BC, total.elas_GB_BC), 
      type='l', lty=3, col=cols[3])
lines(x1_seedlings, rowSums(total.elas_GA_C, total.elas_GB_C), 
      type='l', lty=4, col=cols[4])
lines(x1_seedlings, rowSums(total.elas_GA_FP, total.elas_GB_FP),
      type='l', lty=5, col=cols[5])
lines(x1_seedlings, rowSums(total.elas_GA_CC, total.elas_GB_CC),
      type='l', lty=6, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=3, "(c)", line=1, cex=1.2, adj=0)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

#marginal elasticity maturation by height
plot(x2_seedlings, colSums(total.elas_GA_PG, total.elas_GB_PG),
     ylim= c(0, 0.03), xlim=c(0, 16), xlab="",
     ylab="", main="", type='l', lwd=1, lty=1, cex.axis=1.5, 
     col=cols[1])
lines(x2_seedlings, colSums(total.elas_GA_WT, total.elas_GB_WT),
      type='l', lwd=1, lty=2, col=cols[2])
lines(x2_seedlings, colSums(total.elas_GA_BC, total.elas_GB_BC), 
      type='l', lwd=1, lty=3, col=cols[3])
lines(x2_seedlings, colSums(total.elas_GA_C, total.elas_GB_C),
      type='l', lwd=1, lty=4, col=cols[4])
lines(x2_seedlings, colSums(total.elas_GA_FP, total.elas_GB_FP),
      type='l', lwd=1, lty=5, col=cols[5])
lines(x2_seedlings, colSums(total.elas_GA_CC, total.elas_GB_CC),
      type='l', lwd=1, lty=6, col=cols[6])
text(2, 0.029, "Maturation", cex=1.5)
mtext(side=3, "(d)", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

#marginal elasticity of fertility by diameter
plot(x1_adults, c(rowSums(total.elas_FA_PG), 
                  rowSums(total.elas_FB_PG)),
     ylim=c(0, 0.03), xlim=c(1.6, 800), xlab="", 
     ylab= "", main="", type='l', lwd=1, lty=1, cex.axis=1.5, 
     col=cols[1])
text(65, 0.029, "Fertility", cex=1.5)
lines(x1_adults, c(rowSums(total.elas_FA_WT)
                   , rowSums(total.elas_FB_WT)), type='l', lwd=1,
      lty=2, col=cols[2])
lines(x1_adults, c(rowSums(total.elas_FA_BC),
                   rowSums(total.elas_FB_BC)), type='l', lwd=1,
      lty=3, col=cols[3])
lines(x1_adults, c(rowSums(total.elas_FA_C), 
                   rowSums(total.elas_FB_C)), type='l', lwd=1,
      lty=4, col=cols[4])
lines(x1_adults, c(rowSums(total.elas_FA_FP),
                   rowSums(total.elas_FB_FP)), type='l', lwd=1,
      lty=5, col=cols[5])
lines(x1_adults, c(rowSums(total.elas_FA_CC), 
                   rowSums(total.elas_FB_CC)), type='l', lwd=1, 
      lty=6, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=3, "(e)", line=1, cex=1.2, adj=0)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

#marginal elasticity fertility by height
plot(x2_adults, c(colSums(total.elas_FA_PG), 
                  colSums(total.elas_FB_PG)), ylim=c(0, 0.03),
     xlim=c(16,800),  xlab=" ", ylab= "", main="", type='l', lwd=1,
     lty=1, cex.axis=1.5, col=cols[1])
text(65, 0.029, "Fertility", cex=1.5)
lines(x2_adults, c(colSums(total.elas_FA_WT), 
                   colSums(total.elas_FB_WT)), type='l', lwd=1, 
      lty=2, col=cols[2])
lines(x2_adults, c(colSums(total.elas_FA_BC), 
                   colSums(total.elas_FB_BC)), type='l', lwd=1,
      lty=3, col=cols[3])
lines(x2_adults, c(colSums(total.elas_FA_C), 
                   colSums(total.elas_FB_C)), type='l', lwd=1,
      lty=4, col=cols[4])
lines(x2_adults, c(colSums(total.elas_FA_FP),
                   colSums(total.elas_FB_FP)), type='l', lwd=1, 
      lty=5, col=cols[5])
lines(x2_adults, c(colSums(total.elas_FA_FP), 
                   colSums(total.elas_FB_FP)), type='l', lwd=1,
      lty=6, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=3, "(f)", line=1, cex=1.2, adj=0)
mtext(side=1, "Height (cm)", line=2,cex=1)

#marginal elasticity adults by diameter
plot(x1_adults, c(rowSums(total.elas_D2AA_PG, total.elas_D2AB_PG), 
                  rowSums(total.elas_D2BA_PG, total.elas_D2BB_PG)),
     ylim= c(0, 0.03), xlim=c(1.6, 800), xlab=" ",
     ylab="", main="",  type='l',  lwd=1, lty=1, cex.axis=1.5,
     col=cols[1])
text(65, 0.029, "Adult", cex=1.5)
lines(x1_adults, c(rowSums(total.elas_D2AA_WT, total.elas_D2AB_WT), 
                   rowSums(total.elas_D2BA_WT, total.elas_D2BB_WT)),
      type='l', lwd=1, lty=2, col=cols[2])
lines(x1_adults, c(rowSums(total.elas_D2AA_BC, total.elas_D2AB_BC), 
                   rowSums(total.elas_D2BA_BC, total.elas_D2BB_BC)),
      type='l', lwd=1, lty=3, col=cols[3])
lines(x1_adults, c(rowSums(total.elas_D2AA_C, total.elas_D2AB_C), 
                   rowSums(total.elas_D2BA_C, total.elas_D2BB_C)),
      type='l', lwd=1, lty=4, col=cols[4])
lines(x1_adults, c(rowSums(total.elas_D2AA_FP, total.elas_D2AB_FP), 
                   rowSums(total.elas_D2BA_FP, total.elas_D2BB_FP)),
      type='l', lwd=1, lty=5, col=cols[5])
lines(x1_adults, c(rowSums(total.elas_D2AA_CC, total.elas_D2AB_CC), 
                   rowSums(total.elas_D2BA_CC, total.elas_D2BB_CC)),
      type='l', lwd=1, lty=6, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=3, "(g)", line=1, cex=1.2, adj=0)
mtext(side=1, "Diameter (mm)", line=2,cex=1)
#marginal elasticity adults by height
plot(x2_adults, c(colSums(total.elas_D2AA_PG, total.elas_D2AB_PG), 
                  colSums(total.elas_D2BA_PG, total.elas_D2BB_PG)),
     ylim= c(0, 0.03), xlim=c(16, 800),  xlab=" ", ylab="", main="", 
     type='l',  lwd=1, lty=1, cex.axis=1.5, col=cols[1])
lines(x2_adults, c(colSums(total.elas_D2AA_WT, total.elas_D2AB_WT), 
                   colSums(total.elas_D2BA_WT, total.elas_D2BB_WT)), 
      type='l', lwd=1, lty=2, col=cols[2])
lines(x2_adults, c(colSums(total.elas_D2AA_BC, total.elas_D2AB_BC), 
                   colSums(total.elas_D2BA_BC, total.elas_D2BB_BC)), 
      type='l', lwd=1, lty=3, col=cols[3])
lines(x2_adults, c(colSums(total.elas_D2AA_C, total.elas_D2AB_C), 
                   colSums(total.elas_D2BA_C, total.elas_D2BB_C)),
      type='l', lwd=1, lty=4, col=cols[4])
lines(x2_adults, c(colSums(total.elas_D2AA_FP, total.elas_D2AB_FP), 
                   colSums(total.elas_D2BA_FP, total.elas_D2BB_FP)),
      type='l', lwd=1, lty=5, col=cols[5])
lines(x2_adults, c(colSums(total.elas_D2AA_CC, total.elas_D2AB_CC), 
                   colSums(total.elas_D2BA_CC, total.elas_D2BB_CC)),
      type='l', lwd=1, lty=6, col=cols[6])
text(75, 0.029, "Adult", cex=1.5)
legend(250, 0.03, lty=c(1, 2, 3, 4, 5, 6), col=cols, lwd=1,
       c("PG", "WT", "BC", "C", "FP", "CC"),
       seg.len=2)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=3, "(h)", line=1, cex=1.2, adj=0)
mtext(side=1, "Height (cm)", line=2,cex=1)
dev.off() 

#FIG 9: BARPLOT CONTRIBUTIONS #####

#Barplot of the contributions of eahc kernel component to the 
# differences in analytical population growth rate 
load("./BC/BC_pool.RData")
load("./CC/CC_pool.RData")
load("./C/C_pool.RData")
load("./FP/FP_pool.RData")
load("./PG/PG_pool.RData")
load("./WT/WT_pool.RData")

#(b) Pooled figure 
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/barplot_contribution.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize)
par(mar = c(3, 1, 1, 1), # Margins
    mgp = c(0, 0, 0), # Distance of axis tickmark labels (second value)
    tcl = 0.3, # Length of axis tickmarks
    xpd=F,
    #mai=c(0.2,.5,0.2,0.1),
    oma=c(1, 5, 1, 0),
    mfrow=c(1, 6))

BC <- c(BC_pool[1], sum(BC_pool[2:3]), sum(BC_pool[4:5]), sum(BC_pool[6:9]))
CC <- c(CC_pool[1], sum(CC_pool[2:3]), sum(CC_pool[4:5]), sum(CC_pool[6:9]))
C  <- c(C_pool[1],  sum(C_pool[2:3]),  sum(C_pool[4:5]),  sum(C_pool[6:9]))
FP <- c(FP_pool[1], sum(FP_pool[2:3]), sum(FP_pool[4:5]), sum(FP_pool[6:9]))
PG <- c(PG_pool[1], sum(PG_pool[2:3]), sum(PG_pool[4:5]), sum(PG_pool[6:9]))
WT <- c(WT_pool[1], sum(WT_pool[2:3]), sum(WT_pool[4:5]), sum(WT_pool[6:9]))


#1) PG
barplot(PG, col=cols[1], ylim=c(-0.05, 0.09),
        axes=F, cex.names=0.6, cex.axis=1, horiz=F)
axis(side=2)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["site"],"-",lambda["pool"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
mtext(side=1, "PG")
text(x=xx[1], y = PG[1], label = "D1", pos=3, cex=1)
text(x=xx[2], y = PG[2], label = "M",  pos=3, cex=1)
text(x=xx[3], y = PG[3], label = "F",  pos=1, cex=1)
text(x=xx[4], y = PG[4], label = "D2", pos=1, cex = 1)
abline(h=0)
abline(v=6.5, lty=2, xpd=T)
mtext(side=3, "(a)", adj=0)

#2) WT
xx <- barplot(WT, col=cols[2], ylim=c(-0.05, 0.09), 
              axes=F, cex.names=0.6, cex.axis=1)
text(x=xx[1], y = WT[1], label = "D1", pos=1, cex = 1)
text(x=xx[2], y = WT[2], label = "M",  pos=3, cex = 1)
text(x=xx[3], y = WT[3], label = "F",  pos=1, cex = 1)
text(x=xx[4], y = WT[4], label = "D2", pos=3, cex = 1)
mtext(side=1, "WT")
abline(h=0)
abline(v=6.5, lty=2, xpd=T)
mtext(side=3, "(b)", adj=0)

#3) BC
barplot(BC, col=cols[3], ylim=c(-0.05, 0.09),
        axes=F, cex.names=0.6, cex.axis=1, horiz = F)
mtext(side=1, "BC")
text(x=xx[1], y = BC[1], label = "D1", pos=1, cex=1)
text(x=xx[2], y = BC[2], label = "M",  pos=3, cex=1)
text(x=xx[3], y = BC[3], label = "F",  pos=3, cex=1)
text(x=xx[4], y = BC[4], label = "D2", pos=1, cex = 1)
abline(h=0)
abline(v=6.5, lty=2, xpd=T)
mtext(side=3, "(c)", adj=0)

#)4 C
barplot(C, col= cols[4], ylim=c(-0.05, 0.09), 
        axes=F, cex.names=0.6, cex.axis=1, horiz=F)
mtext(side=1, "C")
text(x=xx[1], y = C[1], label = "D1", pos=1, cex=1)
text(x=xx[2], y = C[2], label = "M",  pos=1, cex=1)
text(x=xx[3], y = C[3], label = "F",  pos=1, cex=1)
text(x=xx[4], y = C[4], label = "D2", pos=1, cex = 1)
abline(h=0)
abline(v=6.5, lty=2, xpd=T)
mtext(side=3, "(d)", adj=0)

#5) FP
xx<- barplot(FP, col=cols[5], ylim=c(-0.05, 0.09), 
             axes=F, cex.names=0.6, cex.axis=1, horiz = F)
text(x=xx[1], y = FP[1], label = "D1", pos=3, cex=1)
text(x=xx[2], y = FP[2], label = "M",  pos=3, cex=1)
text(x=xx[3], y = FP[3], label = "F",  pos=1, cex=1)
text(x=xx[4], y = FP[4], label = "D2", pos=1, cex = 1)
mtext(side=1, "FP")
abline(h=0)
abline(v=6.5, lty=2, xpd=T)
mtext(side=3, "(e)", adj=0)

#6)CC
xx<-barplot(CC, col=cols[6], ylim=c(-0.05, 0.09), 
            axes=F, cex.names=0.6, cex.axis=1, horiz = F)
mtext(side=1, "CC")
text(x=xx[1], y = CC[1], label = "D1", pos=3, cex=1)
text(x=xx[2], y = CC[2], label = "M",  pos=3, cex=1)
text(x=xx[3], y = CC[3], label = "F",  pos=1, cex=1)
text(x=xx[4], y = CC[4], label = "D2", pos=1, cex = 1)
abline(h=0)
mtext(side=3, "(f)", adj=0)
dev.off()

# FIG A.1,A.2 D1 GROWTH  #####
##Punta Gorda
#Grab parameters associated with future diameter
b0<-p.vec_PG[4]
b1<-p.vec_PG[5]
b2<-p.vec_PG[6]

#Grab parameters associated with future height
b3<-p.vec_PG[8]
b4<-p.vec_PG[9]
b5<-p.vec_PG[10]

df$diam_PG <- b0 + b1*df$Var1 + b2*df$Var2
df$height_PG <- b3 + b4*df$Var1 + b5*df$Var2

#Punta Gorda
plot1_PG<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_PG)) + 
  scale_fill_viridis_c(limits = c(0, 1.7)) +
  ggtitle("(a) Punta Gorda") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Punta Gorda"),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(0.3, 1.5)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_PG<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_PG)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(a) Punta Gorda") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Punta Gorda"), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(3, 15)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Wild Turkey
#Grab parameters associated with future diameter
b0<-p.vec_WT[4]
b1<-p.vec_WT[5]
b2<-p.vec_WT[6]

#Grab parameters associated with future height
b3<-p.vec_WT[8]
b4<-p.vec_WT[9]
b5<-p.vec_WT[10]

df$diam_WT <- b0 + b1*df$Var1 + b2*df$Var2
df$height_WT <- b3 + b4*df$Var1 + b5*df$Var2

#Wild Turkey
plot1_WT<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_WT)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(b) Wild Turkey") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Wild Turkey"), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1),
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(0.3, 1.5)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_WT<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_WT)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(b) Wild Turkey") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Wild Turkey"), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(3, 15)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Big Cypress
#Grab parameters associated with future diameter
b0<-p.vec_BC[4]
b1<-p.vec_BC[5]
b2<- p.vec_BC[6]

#Grab parameters associated with future height
b3<-p.vec_BC[8]
b4<-p.vec_BC[9]
b5<-p.vec_BC[10]

df$diam_BC <- b0 + b1*df$Var1 + b2*df$Var2
df$height_BC <- b3 + b4*df$Var1 + b5*df$Var2

#Big Cypress
plot1_BC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_BC)) +
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(c) Big Cypress") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Big Cypress"), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1),
             shape=1, show.legend=F)+
  scale_color_viridis(limits=c(0,1.6)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_BC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_BC)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(c) Big Cypress") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Big Cypress"), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(3, 15)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Chekika
#Grab parameters associated with future diameter
b0<-p.vec_C[4]
b1<-p.vec_C[5]
b2<-p.vec_C[6]
#Grab parameters associated with future height
b3<-p.vec_C[8]
b4<-p.vec_C[9]
b5<-p.vec_C[10]

df$diam_C <- b0 + b1*df$Var1 + b2*df$Var2
df$height_C <- b3 + b4*df$Var1 + b5*df$Var2

#Chekika
plot1_C<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_C)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(d) Chekika") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Chekika"),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend=F)+
  scale_color_viridis(limits=c(0.3, 1.5)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))


plot2_C<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_C)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(d) Chekika") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Chekika"), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(3, 15)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Fort Pierce
#Grab parameters associated with future diameter
b0<-p.vec_FP[4]
b1<-p.vec_FP[5]
b2<-p.vec_FP[6]
#Grab parameters associated with future height
b3<-p.vec_FP[8]
b4<-p.vec_FP[9]
b5<-p.vec_FP[10]
df$diam_FP <- b0 + b1*df$Var1 + b2*df$Var2
df$height_FP <- b3 + b4*df$Var1 + b5*df$Var2

#Fort Pierce
plot1_FP<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_FP)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(e) Fort Pierce") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Fort Pierce"),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(0.3, 1.5)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_FP<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_FP)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(e) Fort Pierce") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Fort Pierce"), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(3, 15)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Cape Canaveral
#Grab parameters associated with future diameter
b0 <- p.vec_CC[4]
b1 <- p.vec_CC[5]
b2 <- p.vec_CC[6]

#Grab parameters associated with future height
b3<- p.vec_CC[8]
b4<- p.vec_CC[9]
b5<- p.vec_CC[10]

df$diam_CC <- b0 + b1*df$Var1 + b2*df$Var2
df$height_CC <- b3 + b4*df$Var1 + b5*df$Var2

#Cape Canaveral
plot1_CC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_CC)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(f) Cape Canaveral") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Cape Canaveral"), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend = F)+
  scale_color_viridis(limits=c(0.3, 1.5)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_CC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_CC)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(f) Cape Canaveral") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Cape Canaveral"), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F, shape=1 ) +
  scale_color_viridis(limits=c(3, 15)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

g <- ggplotGrob(plot1_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot1_BC2 <- plot1_BC + theme(legend.position ="none")
plot1_CC2 <- plot1_CC + theme(legend.position ="none")
plot1_C2  <- plot1_C  + theme(legend.position = "none")
plot1_FP2 <- plot1_FP + theme(legend.position = "none")
plot1_PG2 <- plot1_PG + theme(legend.position = "none")
plot1_WT2 <- plot1_WT + theme(legend.position = "none")

postscript("./Figures/D1_growth_diameter2.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot1_PG2, plot1_WT2,plot1_BC2, plot1_C2,
             plot1_FP2, plot1_CC2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()


g <- ggplotGrob(plot2_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot2_BC2 <- plot2_BC + theme(legend.position ="none")
plot2_CC2 <- plot2_CC + theme(legend.position ="none")
plot2_C2  <- plot2_C  + theme(legend.position = "none")
plot2_FP2 <- plot2_FP + theme(legend.position = "none")
plot2_PG2 <- plot2_PG + theme(legend.position = "none")
plot2_WT2 <- plot2_WT + theme(legend.position = "none")


postscript("./Figures/D1_growth_height2.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot2_PG2, plot2_WT2,plot2_BC2, plot2_C2,
             plot2_FP2, plot2_CC2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()

### FIG A.3,A.4: D2 Survival ##### 
df2 <- expand.grid(x1seq_larges, x2seq_larges)
df3 <- expand.grid(x1seq_larges, x2seq_larges[1])
df4 <- expand.grid(x1seq_larges[1], x2seq_larges)
df5 <- expand.grid(x1seq_larges[10], x2seq_larges)
df7 <- expand.grid(x1seq_larges[20], x2seq_larges)
df8 <- expand.grid(x1seq_larges[30], x2seq_larges)
df9 <- expand.grid(x1seq_larges[40], x2seq_larges)
df10 <- expand.grid(x1seq_larges[50], x2seq_larges)
#PG, BC, CC
b0 <- p.vec_BC[19]
b1 <- p.vec_BC[20]
b2 <- p.vec_BC[21]
df2$BC_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

plot_BCs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= BC_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
<<<<<<< HEAD
  ggtitle("(a) PG, BC, CC") + 
=======
  ggtitle("(a) BC, CC, PG") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
 guides(fill = guide_colorbar(title.position = "right", order =2))

plot_BCs2_meas <- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
geom_point(data=subset(larges,
                !is.na(larges$Surv_tplus1) & (larges$Site == "Big Cypress" |
                            larges$Site == "Cape Canaveral" |
                            larges$Site == "Punta Gorda")),
           aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
           shape=1, size=0.8) +
  ggtitle("(a) ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Chekika
plot_Cs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= C_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(b) C") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot_Cs2_meas <- ggplot()  + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
              !is.na(larges$Surv_tplus1) & (larges$Site == "Chekika")),
              aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
               shape=1, size=0.8) +
  ggtitle(" (b) ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Fort Pierce
plot_FPs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= FP_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(c) FP") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


plot_FPs2_meas <- ggplot()  + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
                  !is.na(larges$Surv_tplus1) &
                    (larges$Site == "Fort Pierce")),
                      aes(x=Diameter_t, y=Height_t, 
                      col = factor(Surv_tplus1)), 
                      shape=1, size=0.8) +
  ggtitle("(c)") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Wild Turkey
plot_WTs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= WT_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(d) WT") + 
>>>>>>> 3e0678db558215c0c2ffe614aa1e814956e9ac49
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))
<<<<<<< HEAD
=======

plot_WTs2_meas <- ggplot()  + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
        !is.na(larges$Surv_tplus1) & (larges$Site == "Wild Turkey")),
        aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
        shape=1, size=0.8) +
  ggtitle("(d) ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1))


g1 <- ggplotGrob(plot_BCs2_pred + theme(legend.position = "right"))$grobs
legend1 <- g1[[which(sapply(g1, function(x) x$name) == "guide-box")]]

g2 <- ggplotGrob(plot_BCs2_meas + theme(legend.position = "right"))$grobs
legend2<- g2[[which(sapply(g2, function(x) x$name) == "guide-box")]]

plot_BCs2_meas <- plot_BCs2_meas + theme(legend.position ="none")
plot_Cs2_meas <- plot_Cs2_meas + theme(legend.position ="none")
plot_FPs2_meas  <- plot_FPs2_meas  + theme(legend.position = "none")
plot_WTs2_meas <- plot_WTs2_meas+ theme(legend.position = "none")

plot_BCs2_pred <- plot_BCs2_pred + theme(legend.position ="none")
plot_Cs2_pred <- plot_Cs2_pred + theme(legend.position ="none")
plot_FPs2_pred   <- plot_FPs2_pred  + theme(legend.position = "none")
plot_WTs2_pred <- plot_WTs2_pred+ theme(legend.position = "none")



postscript("./Figures/D2_survival2.eps", width=width.cm/2.54, 
           height=(4*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 9, NA),
             c(1, 1, 2, 2, 9, NA),
             c(3, 3, 4, 4, 9, NA),
             c(3, 3, 4, 4, 9, 9),
             c(5, 5, 6, 6, 10, 10),
             c(5, 5, 6, 6, 10, NA),
             c(7, 7, 8, 8, 10, NA),
             c(7, 7, 8, 8, 10, NA))
grid.arrange(plot_BCs2_pred, plot_BCs2_meas,
             plot_Cs2_pred, plot_Cs2_meas,
             plot_FPs2_pred, plot_FPs2_meas,
             plot_WTs2_pred, plot_WTs2_meas,
             legend1, legend2, 
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()




#*-- New Separate #####


postscript("./Figures/D2_survival_predicted.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 5, 5),
             c(1, 1, 2, 2, 5, 5),
             c(3, 3, 4, 4, 5, 5),
             c(3, 3, 4, 4, 5, 5))
        
grid.arrange(plot_BCs2_pred, 
             plot_Cs2_pred, 
             plot_FPs2_pred, 
             plot_WTs2_pred, 
             legend1, 
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, 0.7))
dev.off()

postscript("./Figures/D2_survival_measured.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 5, 5),
             c(1, 1, 2, 2, 5, 5),
             c(3, 3, 4, 4, 5, 5),
             c(3, 3, 4, 4, 5, 5))

grid.arrange(plot_BCs2_meas, 
             plot_Cs2_meas, 
             plot_FPs2_meas, 
             plot_WTs2_meas, 
             legend2, 
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, 0.7))
dev.off()
#####

 ggplot() + xlim(1.6, 100) + ylim(15, 18) + 
  geom_tile(data = df3, aes(x=Var1, y = Var2, fill= overall_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
###FIG 4:D2 GROWTH OF D2 INDIVIDUALS #####

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
#*--- ggplot #####


##Big Cypress
#Grab parameters associated with future diameter
b0<-p.vec_BC[22]
b1<-p.vec_BC[23]
b2<- p.vec_BC[24]

#Grab parameters associated with future height
b3<-p.vec_BC[26]
b4<-p.vec_BC[27]
b5<-p.vec_BC[28]

df2$diam2_BC <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_BC <- b3 + b4*df2$Var1 + b5*df2$Var2

##Cape Canaveral
#Grab parameters associated with future diameter
b0 <- p.vec_CC[22]
b1 <- p.vec_CC[23]
b2 <- p.vec_CC[24]

#Grab parameters associated with future height
b3<- p.vec_CC[26]
b4<- p.vec_CC[27]
b5<- p.vec_CC[28]

df2$diam2_CC <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_CC <- b3 + b4*df2$Var1 + b5*df2$Var2
##Chekika
#Grab parameters associated with future diameter
b0<-p.vec_C[22]
b1<-p.vec_C[23]
b2<-p.vec_C[24]
#Grab parameters associated with future height
b3<-p.vec_C[26]
b4<-p.vec_C[27]
b5<-p.vec_C[28]

df2$diam2_C <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_C <- b3 + b4*df2$Var1 + b5*df2$Var2

##Fort Pierce
#Grab parameters associated with future diameter
b0<-p.vec_FP[22]
b1<-p.vec_FP[23]
b2<-p.vec_FP[24]
#Grab parameters associated with future height
b3<-p.vec_FP[26]
b4<-p.vec_FP[27]
b5<-p.vec_FP[28]
df2$diam2_FP <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_FP <- b3 + b4*df2$Var1 + b5*df2$Var2

##Punta Gorda
#Grab parameters associated with future diameter
b0<-p.vec_PG[22]
b1<-p.vec_PG[23]
b2<-p.vec_PG[24]

#Grab parameters associated with future height
b3<-p.vec_PG[26]
b4<-p.vec_PG[27]
b5<-p.vec_PG[28]

df2$diam2_PG <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_PG <- b3 + b4*df2$Var1 + b5*df2$Var2

##Wild Turkey
#Grab parameters associated with future diameter
b0<-p.vec_WT[22]
b1<-p.vec_WT[23]
b2<-p.vec_WT[24]

#Grab parameters associated with future height
b3<-p.vec_WT[26]
b4<-p.vec_WT[27]
b5<-p.vec_WT[28]

df2$diam_WT <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height_WT <- b3 + b4*df2$Var1 + b5*df2$Var2

#overall
b0<-p.vec_overall[26]
b1<-p.vec_overall[27]
b2<-p.vec_overall[28]

df3$diam_WT <- b0 + b1*df3$Var1 + b2*df3$Var2

df4$height <-b0 + b1*df4$Var1 + b2*df4$Var2 
df5$height <-b0 + b1*df5$Var1 + b2*df5$Var2 
df6$height <-b0 + b1*df6$Var1 + b2*df6$Var2 
df7$height <-b0 + b1*df7$Var1 + b2*df7$Var2 
df8$height <-b0 + b1*df8$Var1 + b2*df8$Var2 
df9$height<- b0 + b1*df9$Var1 + b2*df9$Var2
df10$height <-b0 + b1*df10$Var1 + b2*df10$Var2


plot(df4$Var2, df4$height, ylim=c(0,900), col=cols[1])
points(df5$Var2, df5$height, col=cols[2])
points(df6$Var2, df6$height, col=cols[3])
points(df7$Var2, df7$height, col=cols[4])
points(df8$Var2, df8$height, col=cols[5])
points(df9$Var2, df9$height, col=cols[6])
points(df10$Var2, df10$height, col=cols[7])



#x11(width = width.cm/2.54, height = height.cm/2.54, 
#pointsize = pointsize)


#Big Cypress
plot1_BC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_BC)) +
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(a) Big Cypress") + 
  geom_point(data=subset(larges, larges$Site == "Big Cypress"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(limits=c(1.6, 800), range= c(0,8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_BC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_BC)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(a) Big Cypress") + 
  geom_point(data=subset(larges,larges$Site == "Big Cypress"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(limits=c(16, 800), 
                        range = c(0, 8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Cape Canaveral
plot1_CC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_CC)) +
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(b) Cape Canaveral") + 
  geom_point(data=subset(larges, larges$Site == "Cape Canaveral"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(limits=c(1.6, 800), range= c(0,8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_CC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_CC)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(b) Cape Canaveral") + 
  geom_point(data=subset(larges,larges$Site == "Cape Canaveral"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(limits=c(16, 800), 
                        range = c(0, 8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Chekika
plot1_C<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_C)) +
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(c) Chekika") + 
  geom_point(data=subset(larges, larges$Site == "Chekika"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(limits=c(1.6, 800), range= c(0,8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_C<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_C)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(c) Chekika") + 
  geom_point(data=subset(larges,larges$Site == "Chekika"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(limits=c(16, 800), 
                        range = c(0, 8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Fort Pierce
plot1_FP<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_FP)) +
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(d) Fort Pierce") + 
  geom_point(data=subset(larges, larges$Site == "Fort Pierce"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(limits=c(1.6, 800), range= c(0,8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_FP<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_FP)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(d) Fort Pierce") + 
  geom_point(data=subset(larges,larges$Site == "Fort Pierce"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(limits=c(16, 800), 
                        range = c(0, 8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Punta Gorda
plot1_PG<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_PG)) +
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(e) Punta Gorda") + 
  geom_point(data=subset(larges, larges$Site == "Punta Gorda"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(limits=c(1.6, 800), range= c(0,8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_PG<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_PG)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(e) Punta Gorda") + 
  geom_point(data=subset(larges,larges$Site == "Punta Gorda"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(limits=c(16, 800), 
                        range = c(0, 8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))


#Wild Turkey
plot1_WT<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam_WT)) +
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(f) Wild Turkey") + 
  geom_point(data=subset(larges, larges$Site == "Wild Turkey"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(limits=c(1.6, 800), range= c(0,8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_WT<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height_WT)) + 
  scale_fill_viridis_c(limits = c(16, 1000)) +
  ggtitle("(f) Wild Turkey") + 
  geom_point(data=subset(larges,larges$Site == "Wild Turkey"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(limits=c(16, 800), 
                        range = c(0, 8)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))


g <- ggplotGrob(plot1_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot1_BC2 <- plot1_BC + theme(legend.position ="none")
plot1_CC2 <- plot1_CC + theme(legend.position ="none")
plot1_C2  <- plot1_C  + theme(legend.position = "none")
plot1_FP2 <- plot1_FP + theme(legend.position = "none")
plot1_PG2 <- plot1_PG + theme(legend.position = "none")
plot1_WT2 <- plot1_WT + theme(legend.position = "none")

postscript("./Figures/D2_growth_diameter.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot1_BC2, plot1_CC2,plot1_C2, plot1_FP2,
             plot1_PG2, plot1_WT2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()


g <- ggplotGrob(plot2_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot2_BC2 <- plot2_BC + theme(legend.position ="none")
plot2_CC2 <- plot2_CC + theme(legend.position ="none")
plot2_C2  <- plot2_C  + theme(legend.position = "none")
plot2_FP2 <- plot2_FP + theme(legend.position = "none")
plot2_PG2 <- plot2_PG + theme(legend.position = "none")
plot2_WT2 <- plot2_WT + theme(legend.position = "none")



postscript("./Figures/D2_growth_height.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot2_BC2, plot2_CC2,plot2_C2, plot2_FP2,
             plot2_PG2, plot2_WT2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()

#*-- Version 2: Color #####


#Big Cypress
plot1_BC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_BC)) +
  scale_fill_viridis_c(limits = c(1.6, 860)) +
  ggtitle("(a) Big Cypress") + 
  geom_point(data=subset(larges, larges$Site == "Big Cypress" &
                           !is.na(larges$Diameter_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend=F)+
  scale_color_viridis(limits=c(1, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_BC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_BC)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(a) Big Cypress") + 
  geom_point(data=subset(larges,larges$Site == "Big Cypress"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

#Cape Canaveral
plot1_CC<- ggplot() + xlim(1.6, 800) + ylim(1.6, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_CC)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(b) Cape Canaveral") + 
  geom_point(data=subset(larges, larges$Site == "Cape Canaveral"&
                           !is.na(larges$Diameter_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend = F)+
  scale_color_viridis(limits=c(1, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_CC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_CC)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(b) Cape Canaveral") + 
  geom_point(data=subset(larges,larges$Site == "Cape Canaveral"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F, shape=1 ) +
  scale_color_viridis(limits=c(6.5, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

#Chekika
plot1_C<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_C)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(c) Chekika") + 
  geom_point(data=subset(larges, larges$Site == "Chekika"&
                           !is.na(larges$Diameter_tplus1)),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend=F)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))


plot2_C<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_C)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(c) Chekika") + 
  geom_point(data=subset(larges,larges$Site == "Chekika"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

#Fort Pierce
plot1_FP<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_FP)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(d) Fort Pierce") + 
  geom_point(data=subset(larges, larges$Site == "Fort Pierce"&
                           !is.na(larges$Diameter_tplus1)),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_FP<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_FP)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(d) Fort Pierce") + 
  geom_point(data=subset(larges,larges$Site == "Fort Pierce"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

#Punta Gorda
plot1_PG<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_PG)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(e) Punta Gorda") + 
  geom_point(data=subset(larges, larges$Site == "Punta Gorda"&
                           !is.na(larges$Diameter_tplus1)),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_PG<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_PG)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(e) Punta Gorda") + 
  geom_point(data=subset(larges,larges$Site == "Punta Gorda"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

#Wild Turkey
plot1_WT<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam_WT)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(f) Wild Turkey") + 
  geom_point(data=subset(larges, larges$Site == "Wild Turkey"&
                           !is.na(larges$Diameter_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1),
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_WT<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height_WT)) + 
  scale_fill_viridis_c(limits = c(16, 971)) +
  ggtitle("(f) Wild Turkey") + 
  geom_point(data=subset(larges,larges$Site == "Wild Turkey"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))



g <- ggplotGrob(plot1_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot1_BC2 <- plot1_BC + theme(legend.position ="none")
plot1_CC2 <- plot1_CC + theme(legend.position ="none")
plot1_C2  <- plot1_C  + theme(legend.position = "none")
plot1_FP2 <- plot1_FP + theme(legend.position = "none")
plot1_PG2 <- plot1_PG + theme(legend.position = "none")
plot1_WT2 <- plot1_WT + theme(legend.position = "none")

postscript("./Figures/D2_growth_diameter2.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot1_BC2, plot1_CC2,plot1_C2, plot1_FP2,
             plot1_PG2, plot1_WT2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()


g <- ggplotGrob(plot2_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot2_BC2 <- plot2_BC + theme(legend.position ="none")
plot2_CC2 <- plot2_CC + theme(legend.position ="none")
plot2_C2  <- plot2_C  + theme(legend.position = "none")
plot2_FP2 <- plot2_FP + theme(legend.position = "none")
plot2_PG2 <- plot2_PG + theme(legend.position = "none")
plot2_WT2 <- plot2_WT + theme(legend.position = "none")



postscript("./Figures/D2_growth_height2.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot2_BC2, plot2_CC2,plot2_C2, plot2_FP2,
             plot2_PG2, plot2_WT2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()

#*-- Version 3: Persp

#Grab parameters associated with future diameter
b0<-p.vec_BC[22]
b1<-p.vec_BC[23]
b2<-p.vec_BC[24]
z1_BC<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0 + b1*a + b2*b))

b0<-p.vec_CC[22]
b1<-p.vec_CC[23]
b2<-p.vec_CC[24]
z1_CC<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0 + b1*a + b2*b))

b0<-p.vec_C[22]
b1<-p.vec_C[23]
b2<-p.vec_C[24]
z1_C<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0 + b1*a + b2*b))

b0<-p.vec_FP[22]
b1<-p.vec_FP[23]
b2<-p.vec_FP[24]
z1_FP<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0 + b1*a + b2*b))

b0<-p.vec_PG[22]
b1<-p.vec_PG[23]
b2<-p.vec_PG[24]
z1_PG<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0 + b1*a + b2*b))

b0<-p.vec_WT[22]
b1<-p.vec_WT[23]
b2<-p.vec_WT[24]
z1_WT<-outer(x1seq_larges, x2seq_larges, function(a,b) (b0 + b1*a + b2*b))








setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D2_growth_diameter3.eps", width=width.cm/2.54, 
           height=3*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
#x11(width = width.cm/2.54, height = 3*height.cm/2.54, 
#  pointsize = pointsize)
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(3, 2))

#A) Big Cypress
pmat<-persp(x1seq_larges, x2seq_larges, z1_BC, ticktype="detailed",
            theta=-30, zlim=c(1.6, 1200), zlab="\n Diameter at t+1 (mm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(a) Big Cypress", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Big Cypress"],
                    larges$Height_t[larges$Site == "Big Cypress"],
                    larges$Diameter_tplus1[larges$Site == "Big Cypress"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#B) Cape Canaveral
pmat<-persp(x1seq_larges, x2seq_larges, z1_CC, ticktype="detailed",
            theta=-30, zlim=c(1.6, 1200), zlab="\n Diameter at t+1 (mm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(b) Cape Canaveral", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Cape Canaveral"],
                    larges$Height_t[larges$Site == "Cape Canaveral"],
                    larges$Diameter_tplus1[larges$Site == "Cape Canaveral"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#C) Chekika
pmat<-persp(x1seq_larges, x2seq_larges, z1_C, ticktype="detailed",
            theta=-30, zlim=c(1.6, 1200), zlab="\n Diameter at t+1 (mm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(c) Chekika", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Chekika"],
                    larges$Height_t[larges$Site == "Chekika"],
                    larges$Diameter_tplus1[larges$Site == "Chekika"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#D) Fort Pierce
pmat<-persp(x1seq_larges, x2seq_larges, z1_FP, ticktype="detailed",
            theta=-30, zlim=c(1.6, 1200), zlab="\n Diameter at t+1 (mm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(d) Fort Pierce", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Fort Pierce"],
                    larges$Height_t[larges$Site == "Fort Pierce"],
                    larges$Diameter_tplus1[larges$Site == "Fort Pierce"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#E) Punta Gorda
pmat<-persp(x1seq_larges, x2seq_larges, z1_PG, ticktype="detailed",
            theta=-30, zlim=c(1.6, 1200), zlab="\n Diameter at t+1 (mm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(e) Punta Gorda", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Punta Gorda"],
                    larges$Height_t[larges$Site == "Punta Gorda"],
                    larges$Diameter_tplus1[larges$Site == "Punta Gorda"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#F) Wild Turkey
pmat<-persp(x1seq_larges, x2seq_larges, z1_WT, ticktype="detailed",
            theta=-30, zlim=c(1.6, 1200), zlab="\n Diameter at t+1 (mm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(f) Wild Turkey", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Wild Turkey"],
                    larges$Height_t[larges$Site == "Wild Turkey"],
                    larges$Diameter_tplus1[larges$Site == "Wild Turkey"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)
dev.off()




#Now for height
#Grab parameters associated with future height
b3<-p.vec_BC[26]
b4<-p.vec_BC[27]
b5<-p.vec_BC[28]
z2_BC<- outer(x1seq_larges, x2seq_larges, function(a,b) (b3+ b4*a + b5*b))

b3<-p.vec_CC[26]
b4<-p.vec_CC[27]
b5<-p.vec_CC[28]
z2_CC<- outer(x1seq_larges, x2seq_larges, function(a,b) (b3+ b4*a + b5*b))

b3<-p.vec_C[26]
b4<-p.vec_C[27]
b5<-p.vec_C[27]
z2_C<- outer(x1seq_larges, x2seq_larges, function(a,b) (b3+ b4*a + b5*b))

b3<-p.vec_FP[26]
b4<-p.vec_FP[27]
b5<-p.vec_FP[28]
z2_FP<- outer(x1seq_larges, x2seq_larges, function(a,b) (b3+ b4*a + b5*b))

b3<-p.vec_PG[26]
b4<-p.vec_PG[27]
b5<-p.vec_PG[28]
z2_PG<- outer(x1seq_larges, x2seq_larges, function(a,b) (b3+ b4*a + b5*b))

b3<-p.vec_WT[26]
b4<-p.vec_WT[27]
b5<-p.vec_WT[28]
z2_WT<- outer(x1seq_larges, x2seq_larges, function(a,b) (b3+ b4*a + b5*b))




setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D2_growth_height3.eps", width=width.cm/2.54, 
           height=3*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
#x11(width = width.cm/2.54, height = 3*height.cm/2.54, 
#  pointsize = pointsize)
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.1,0.5,0.1,0.1),
    mfrow=c(3, 2))

#A) Big Cypress
pmat<-persp(x1seq_larges, x2seq_larges, z2_BC, ticktype="detailed",
            theta=-30, zlim=c(16, 950), zlab="\n Height at t+1 (cm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(a) Big Cypress", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Big Cypress"],
                    larges$Height_t[larges$Site == "Big Cypress"],
                    larges$Height_tplus1[larges$Site == "Big Cypress"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#B) Cape Canaveral
pmat<-persp(x1seq_larges, x2seq_larges, z2_CC, ticktype="detailed",
            theta=-30, zlim=c(16, 950), zlab="\n Height at t+1 (cm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(b) Cape Canaveral", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Cape Canaveral"],
                    larges$Height_t[larges$Site == "Cape Canaveral"],
                    larges$Height_tplus1[larges$Site == "Cape Canaveral"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#C) Chekika
pmat<-persp(x1seq_larges, x2seq_larges, z2_C, ticktype="detailed",
            theta=-30, zlim=c(16, 950), zlab="\n Height at t+1 (cm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(c) Chekika", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Chekika"],
                    larges$Height_t[larges$Site == "Chekika"],
                    larges$Height_tplus1[larges$Site == "Chekika"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#D) Fort Pierce
pmat<-persp(x1seq_larges, x2seq_larges, z2_FP, ticktype="detailed",
            theta=-30, zlim=c(16, 950), zlab="\n Height at t+1 (cm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(d) Fort Pierce", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Fort Pierce"],
                    larges$Height_t[larges$Site == "Fort Pierce"],
                    larges$Height_tplus1[larges$Site == "Fort Pierce"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#E) Punta Gorda
pmat<-persp(x1seq_larges, x2seq_larges, z2_PG, ticktype="detailed",
            theta=-30, zlim=c(16, 950), zlab="\n Height at t+1 (cm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(e) Punta Gorda", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Punta Gorda"],
                    larges$Height_t[larges$Site == "Punta Gorda"],
                    larges$Height_tplus1[larges$Site == "Punta Gorda"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)

#F) Wild Turkey
pmat<-persp(x1seq_larges, x2seq_larges, z2_WT, ticktype="detailed",
            theta=-30, zlim=c(16, 950), zlab="\n Height at t+1 (cm)", 
            shade=0.1, nticks=4, xlab="\n Diameter at t(mm)", 
            ylab="\n Height at t(cm)",  lwd=1, xaxs="i", 
            main="", cex.main=1)
mtext(side=3, "(f) Wild Turkey", line=-3)

# from 3D to 2D coordinates
mypoints <- trans3d(larges$Diameter_t[larges$Site == "Wild Turkey"],
                    larges$Height_t[larges$Site == "Wild Turkey"],
                    larges$Height_tplus1[larges$Site == "Wild Turkey"],
                    pmat=pmat)
# plot in 2D space with pointsize related to distance
points(mypoints, col=4)
dev.off()






###FIGURE 5: GRADUATION #####
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

#*--- in ggplot #####
#BC, CC, FP and PG
b0 <- p.vec_BC[12]
b1 <- p.vec_BC[13]
b2 <- p.vec_BC[14]
df$BC_m <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

#Chekika
b0 <- p.vec_C[12]
b1 <- p.vec_C[13]
b2 <- p.vec_C[14]
df$C_m <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

#WT
b0 <- p.vec_WT[12]
b1 <- p.vec_WT[13]
b2 <- p.vec_WT[14]
df$WT_m <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

#Big Cypress, Cape Canaveral, Fort Pierce, Punta Gorda
plot_BCm<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= BC_m)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(a) BC, CC, FP, PG") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$grad_status) &
                           (seedlings$Site == "Big Cypress" |
                              seedlings$Site == "Cape Canaveral" |
                              seedlings$Site == "Fort Pierce" |
                              seedlings$Site == "Punta Gorda")),
             aes(x=Diameter_t, y=Height_t, shape = factor(grad_status))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n maturation \n at time \n t + 1") +
  labs(shape = "Measured \n maturation \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Chekika
plot_Cm<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= C_m)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(b) C") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$grad_status) &
                           (seedlings$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t, shape = factor(grad_status))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n maturation \n at time \n t + 1") +
  labs(shape = "Measured \n maturation \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Wild Turkey
plot_WTm<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= WT_m)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(c) WT") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$grad_status) &
                           (seedlings$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t, shape = factor(grad_status))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n maturation \n at time \n t + 1") +
  labs(shape = "Measured \n maturation \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

g <- ggplotGrob(plot_BCm + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot_BCm2 <- plot_BCm + theme(legend.position ="none")
plot_Cm2<- plot_Cm + theme(legend.position ="none")
plot_WTm2  <- plot_WTm  + theme(legend.position = "none")




postscript("./Figures/Maturation.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2),
             c(1, 1, 2, 2),
             c(3, 3, 4, 4),
             c(3, 3, 4, 4))
grid.arrange(plot_BCm2, plot_Cm2,plot_WTm2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1))
dev.off()




###FIGURE 6: REPRODUCTION #####



#BC, CC, and FP
b0<-p.vec_BC[30]
b1<-p.vec_BC[31]
pf_BC <- invlogit(b0 + b1*x2seq_larges)

#Chekika
b0<-p.vec_C[30]
b1<-p.vec_C[31]
pf_C <- invlogit(b0 + b1*x2seq_larges)

#Punta Gorda
b0<-p.vec_PG[30]
b1<-p.vec_PG[31]
pf_PG <- invlogit(b0 + b1*x2seq_larges)

#Wild Turkey
b0<-p.vec_WT[30]
b1<-p.vec_WT[31]
pf_WT <- invlogit(b0 + b1*x2seq_larges)

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

# TODO: Add datapoints to this plot
plot(x2seq_larges, pf_BC, xlab="Height at time t (cm)", 
     ylab="\n \n Probability of Reproducing", lwd=1, type='l', lty=1, col=cols[1])
lines(x2seq_larges, pf_C, lwd=2, type='l', lty=3, col=cols[3])
lines(x2seq_larges, pf_PG, lwd=2, type='l', lty=5, col=cols[5])
lines(x2seq_larges, pf_WT, lwd=2, type='l', lty=6, col=cols[6])
legend(5, 0.98, c("BC-CC-FP", "C", "PG", "WT"), 
       lty=c(1, 3,  5, 6), lwd=2, seg.len=1.5, col=cols[c(1,3,5,6)])
dev.off()


#*---With datapoints #####
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/reproduction2.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm_onepanel/2.54, height = width.cm_onepanel/2.54, 
#     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))


# TODO: Add datapoints to this plot
plot(x2seq_larges, pf_BC, xlab="Height at time t (cm)", ylim=c(0,1),
     ylab="\n \n Probability of Reproducing", lwd=1, type='l', lty=1, col=cols[1])
points(larges$Height_t[larges$Site == "Big Cypress" | 
                      larges$Site == "Cape Canaveral" |
                      larges$Site == "Fort Pierce"], 
       larges$Rep_tplus1[larges$Site == "Big Cypress" | 
       larges$Site == "Cape Canaveral" |
       larges$Site == "Fort Pierce"],
       col=cols[1])
lines(x2seq_larges, pf_C, lwd=2, type='l', lty=3, col=cols[3])
points(larges$Height_t[larges$Site == "Chekika"], 
       larges$Rep_tplus1[larges$Site == "Chekika"],
       col=cols[3])
lines(x2seq_larges, pf_PG, lwd=2, type='l', lty=5, col=cols[5])
points(larges$Height_t[larges$Site == "Punta Gorda"], 
       larges$Rep_tplus1[larges$Site == "Punta Gorda"],
       col=cols[5])
lines(x2seq_larges, pf_WT, lwd=2, type='l', lty=6, col=cols[6])
points(larges$Height_t[larges$Site == "Wild Turkey"], 
       larges$Rep_tplus1[larges$Site == "Wild Turkey"],
       col=cols[6])
legend(5, 0.98, c("BC-CC-FP", "C", "PG", "WT"), 
       lty=c(1, 3,  5, 6), lwd=2, seg.len=2, col=cols[c(1,3,5,6)])
dev.off()


###FIGURE 7: FECUNDITY #####
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

y_BC <- p.vec_BC[32]*x1seq_larges*x1seq_larges
y_CC <- p.vec_CC[32]*x1seq_larges*x1seq_larges
y_C  <- p.vec_C[32]*x1seq_larges*x1seq_larges
y_FP <- p.vec_FP[32]*x1seq_larges*x1seq_larges
y_PG <- p.vec_PG[32]*x1seq_larges*x1seq_larges
y_WT <- p.vec_WT[32]*x1seq_larges*x1seq_larges

x<-Stems$diam_base*10 #convert from cm to mm
x <- x*x

plot(x1seq_larges, y_BC, xlab="Diameter at time t (mm)", 
     ylab="\n \n Number of Offspring",  lwd=2, type='l', lty=1, col=cols[1])
lines(x1seq_larges, y_CC, lwd=2, type='l', lty=2, col=cols[2])
lines(x1seq_larges, y_C,  lwd=2, type='l', lty=3, col=cols[3])
lines(x1seq_larges, y_FP, lwd=2, type='l', lty=4, col=cols[4])
lines(x1seq_larges, y_PG, lwd=2, type='l', lty=5, col=cols[5])
lines(x1seq_larges, y_WT, lwd=2, type='l', lty=6, col=cols[6])
#lines(x1seq_larges, y_overall, lwd=3, type='l', lty=2, col="red")
legend(5, 2500000, lwd=2,  c("BC", "CC", "C","FP", "PG","WT" ),   lty=1:6, seg.len=3, col=cols)
dev.off()


#*---With datapoints #####
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/fecundity2.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))

x<-Stems$diam_base*10 #convert from cm to mm


Stems <- cbind(Stems, x)

plot(x1seq_larges, y_BC, xlab="Diameter at time t (mm)",
     ylab="\n \n Number of Offspring",  lwd=2, type='l', lty=1, col=cols[1])
points(Stems$x[Stems$Site == "Big Cypress"], 
       Stems$Seed_No[Stems$Site == "Big Cypress"], col=cols[1])
lines(x1seq_larges, y_CC, lwd=2, type='l', lty=2, col=cols[2])
points(Stems$x[Stems$Site == "Cape Canaveral"],
       Stems$Seed_No[Stems$Site == "Cape Canaveral"], col=cols[2])

lines(x1seq_larges, y_C,  lwd=2, type='l', lty=3, col=cols[3])
points(Stems$x[Stems$Site == "Chekika"],
       Stems$Seed_No[Stems$Site == "Chekika"], col=cols[3])

lines(x1seq_larges, y_FP, lwd=2, type='l', lty=4, col=cols[4])
points(Stems$x[Stems$Site == "Fort Pierce"],
       Stems$Seed_No[Stems$Site == "Fort Pierce"], col=cols[4])

lines(x1seq_larges, y_PG, lwd=2, type='l', lty=5, col=cols[5])
points(Stems$x[Stems$Site == "Punta Gorda"],
       Stems$Seed_No[Stems$Site == "Punta Gorda"], col=cols[5])

lines(x1seq_larges, y_WT, lwd=2, type='l', lty=6, col=cols[6])
points(Stems$x[Stems$Site == "Wild Turkey"],
       Stems$Seed_No[Stems$Site == "Wild Turkey"], col=cols[6])

#lines(x1seq_larges, y_overall, lwd=3, type='l', lty=2, col="red")
legend(5, 2500000, lwd=2,  c("BC", "CC", "C","FP", "PG","WT" ),   lty=1:6, seg.len=3, col=cols)
dev.off()







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

load("./BC/BC.RData")
load("./CC/CC.RData")
load("./C/C.RData")
load("./FP/FP.RData")
load("./PG/PG.RData")
load("./WT/WT.RData")


load("lambda_all.RData")
load("lambda_E.RData")
load("lambda_H.RData")
load("lambda_W.RData")

#thing<-A_cv*elas_all
load("total.thing_D1.RData")
load("total.thing_D2.RData")
load("total.thing_F.RData")
load("total.thing_G.RData")


lam.stable_BC <- readRDS("./BC/lam.stable_BC.rds")
lam.stable_CC <- readRDS("./CC/lam.stable_CC.rds")
lam.stable_C  <- readRDS("./C/lam.stable_C.rds")
lam.stable_FP <- readRDS("./FP/lam.stable_FP.rds")
lam.stable_PG <- readRDS("./PG/lam.stable_PG.rds")
lam.stable_WT <- readRDS("./WT/lam.stable_WT.rds")
#Figure 8: Population Growth Rates #####
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
barCenters <- plot(c(lam.stable_BC, lam.stable_CC, lam.stable_C, lam.stable_FP,
          lam.stable_PG, lam.stable_WT), 
        col=cols,
        cex=2,
        ylab=expression(paste(lambda)),
        axes= F,
        xlab= "Site",
        pch=19,
         xpd=F, 
        ylim =c(0.98, 1.26))
axis(side=1, labels = c("BC", "CC", "C", "FP",
                        "PG", "WT"), at = 1:6 )
axis(side=2)

segments(1, 1.062, 1, 1.206) #Big Cypress
arrows(1, 1.062, 1, 1.206, angle = 90, 
       code = 3, length = 0.05)
segments(2, 1.055, 2, 1.142) #Cape Canaveral
arrows(2, 1.055, 2, 1.142, angle = 90, 
       code = 3, length = 0.05)
segments(3, 0.986, 3, 1.05) #Chekika
arrows(3, 0.986, 3, 1.05, angle = 90, 
       code = 3, length = 0.05)
segments(4, 1.079, 4, 1.221) #Fort Pierce
arrows(4, 1.079, 4, 1.221, angle = 90, code = 3, length = 0.05)
segments(5, 0.995, 5, 1.140) #Punta Gorda
arrows(5, 0.995, 5, 1.140, angle=90, code = 3, length = 0.05)
segments(6, 1.062, 6, 1.251) #Wild Turkey
arrows(6, 1.062, 6, 1.251, angle = 90, code =3, length = 0.05)

dev.off()

load("./BC/BC_pool.RData")
load("./CC/CC_pool.RData")
load("./C/C_pool.RData")
load("./FP/FP_pool.RData")
load("./PG/PG_pool.RData")
load("./WT/WT_pool.RData")
load("./BC/BC_avg.RData")

#Figure 9: Barplot of contribution of kernel components to differences in lambda #####
#(a) Average 
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/barplot_contributions_average.eps", width=width.cm_onepanel/2.54, 
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
    mfrow=c(3, 2))

BC <- c(BC_avg[1], sum(BC_avg[2:3]), sum(BC_avg[4:5]), sum(BC_avg[6:9]))
CC <- c(CC_avg[1], sum(CC_avg[2:3]), sum(CC_avg[4:5]), sum(CC_avg[6:9]))
C  <- c(C_avg[1],  sum(C_avg[2:3]),  sum(C_avg[4:5]),  sum(C_avg[6:9]))
FP <- c(FP_avg[1], sum(FP_avg[2:3]), sum(FP_avg[4:5]), sum(FP_avg[6:9]))
PG <- c(PG_avg[1], sum(PG_avg[2:3]), sum(PG_avg[4:5]), sum(PG_avg[6:9]))
WT <- c(WT_avg[1], sum(WT_avg[2:3]), sum(WT_avg[4:5]), sum(WT_avg[6:9]))


barplot(BC, col=cols[1], ylim=c(-0.05, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), 
        axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["BC"],"-",bar(lambda))))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(CC, col=cols[2], ylim=c(-0.05, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), 
        axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["CC"],"-",bar(lambda))))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(C, col=cols[3], ylim=c(-0.05, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), 
        axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["C"],"-",bar(lambda))))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(FP, col=cols[4], ylim=c(-0.05, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), 
        axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["FP"],"-",bar(lambda))))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(BC, col=cols[5], ylim=c(-0.05, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), 
        axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["PG"],"-",bar(lambda))))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(WT, col=cols[6], ylim=c(-0.05, 0.06),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), 
        axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["WT"],"-",bar(lambda))))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

dev.off()

#(b) Pooled figure 
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/barplot_contributions_pool_dashed_letters.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize)
#png(file="./Figures/barplot_contributions.png")
#x11(width = width.cm_onepanel/2.54, height = width.cm_onepanel/2.54, 
#   pointsize = pointsize)

par(mar = c(3, 1, 1, 1), # Margins
    mgp = c(0, 0, 0), # Distance of axis tickmark labels (second value)
    tcl = 0.3, # Length of axis tickmarks
    xpd=F,
    #mai=c(0.2,.5,0.2,0.1),
    oma=c(1, 5, 1, 0),
    mfrow=c(1, 6))
>>>>>>> 3e0678db558215c0c2ffe614aa1e814956e9ac49

plot_BCs2_meas <- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
                         !is.na(larges$Surv_tplus1) & (larges$Site == "Big Cypress" |
                                                         larges$Site == "Cape Canaveral" |
                                                         larges$Site == "Punta Gorda")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
             shape=1, size=0.8) +
  ggtitle("(a) PG, BC, CC ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Wild Turkey
b0 <- p.vec_WT[19]
b1 <- p.vec_WT[20]
b2 <- p.vec_WT[21]
df2$WT_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

plot_WTs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= WT_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(b) WT") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot_WTs2_meas <- ggplot()  + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
                         !is.na(larges$Surv_tplus1) & (larges$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
             shape=1, size=0.8) +
  ggtitle("(b) WT ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1))

#Chekika
b0 <- p.vec_C[19]
b1 <- p.vec_C[20]
b2 <- p.vec_C[21]
df2$C_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

plot_Cs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= C_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(c) C") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot_Cs2_meas <- ggplot()  + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
                         !is.na(larges$Surv_tplus1) & (larges$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
             shape=1, size=0.8) +
  ggtitle("(c) C ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 


#Fort Pierce
b0 <- p.vec_FP[19]
b1 <- p.vec_FP[20]
b2 <- p.vec_FP[21]
df2$FP_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

plot_FPs2_pred<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= FP_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(d) FP") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


plot_FPs2_meas <- ggplot()  + xlim(1.6, 800) + ylim(16, 800) + 
  geom_point(data=subset(larges,
                  !is.na(larges$Surv_tplus1) &
                    (larges$Site == "Fort Pierce")),
                      aes(x=Diameter_t, y=Height_t, 
                      col = factor(Surv_tplus1)), 
                      shape=1, size=0.8) +
  ggtitle("(d) FP ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 


g1 <- ggplotGrob(plot_BCs2_pred + theme(legend.position = "right"))$grobs
legend1 <- g1[[which(sapply(g1, function(x) x$name) == "guide-box")]]

g2 <- ggplotGrob(plot_BCs2_meas + theme(legend.position = "right"))$grobs
legend2<- g2[[which(sapply(g2, function(x) x$name) == "guide-box")]]

plot_BCs2_meas <- plot_BCs2_meas + theme(legend.position ="none")
plot_Cs2_meas <- plot_Cs2_meas + theme(legend.position ="none")
plot_FPs2_meas  <- plot_FPs2_meas  + theme(legend.position = "none")
plot_WTs2_meas <- plot_WTs2_meas+ theme(legend.position = "none")

plot_BCs2_pred <- plot_BCs2_pred + theme(legend.position ="none")
plot_Cs2_pred <- plot_Cs2_pred + theme(legend.position ="none")
plot_FPs2_pred   <- plot_FPs2_pred  + theme(legend.position = "none")
plot_WTs2_pred <- plot_WTs2_pred+ theme(legend.position = "none")


postscript("./Figures/D2_survival_predicted.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 5, 5),
             c(1, 1, 2, 2, 5, 5),
             c(3, 3, 4, 4, 5, 5),
             c(3, 3, 4, 4, 5, 5))
        
grid.arrange(plot_BCs2_pred, 
             plot_WTs2_pred, 
             plot_Cs2_pred, 
             plot_FPs2_pred, 
             legend1, 
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, 0.7))
dev.off()

postscript("./Figures/D2_survival_measured.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 5, 5),
             c(1, 1, 2, 2, 5, 5),
             c(3, 3, 4, 4, 5, 5),
             c(3, 3, 4, 4, 5, 5))

grid.arrange(plot_BCs2_meas, 
             plot_WTs2_meas, 
             plot_Cs2_meas, 
             plot_FPs2_meas, 
             legend2, 
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, 0.7))
dev.off()

###FIG A.5, A.6: D2 GROWTH OF D2 INDIVIDUALS #####
##Punta Gorda
#Grab parameters associated with future diameter
b0<-p.vec_PG[22]
b1<-p.vec_PG[23]
b2<-p.vec_PG[24]

#Grab parameters associated with future height
b3<-p.vec_PG[26]
b4<-p.vec_PG[27]
b5<-p.vec_PG[28]

df2$diam2_PG <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_PG <- b3 + b4*df2$Var1 + b5*df2$Var2

#Punta Gorda
plot1_PG<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_PG)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(a) Punta Gorda") + 
  geom_point(data=subset(larges, larges$Site == "Punta Gorda"&
                           !is.na(larges$Diameter_tplus1)),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_PG<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_PG)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(a) Punta Gorda") + 
  geom_point(data=subset(larges,larges$Site == "Punta Gorda"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Wild Turkey
#Grab parameters associated with future diameter
b0<-p.vec_WT[22]
b1<-p.vec_WT[23]
b2<-p.vec_WT[24]

#Grab parameters associated with future height
b3<-p.vec_WT[26]
b4<-p.vec_WT[27]
b5<-p.vec_WT[28]

df2$diam_WT <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height_WT <- b3 + b4*df2$Var1 + b5*df2$Var2

#Wild Turkey
plot1_WT<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam_WT)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(b) Wild Turkey") + 
  geom_point(data=subset(larges, larges$Site == "Wild Turkey"&
                           !is.na(larges$Diameter_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1),
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_WT<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height_WT)) + 
  scale_fill_viridis_c(limits = c(16, 971)) +
  ggtitle("(b) Wild Turkey") + 
  geom_point(data=subset(larges,larges$Site == "Wild Turkey"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Big Cypress
#Grab parameters associated with future diameter
b0<-p.vec_BC[22]
b1<-p.vec_BC[23]
b2<- p.vec_BC[24]

#Grab parameters associated with future height
b3<-p.vec_BC[26]
b4<-p.vec_BC[27]
b5<-p.vec_BC[28]

df2$diam2_BC <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_BC <- b3 + b4*df2$Var1 + b5*df2$Var2

#Big Cypress
plot1_BC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_BC)) +
  scale_fill_viridis_c(limits = c(1.6, 860)) +
  ggtitle("(c) Big Cypress") + 
  geom_point(data=subset(larges, larges$Site == "Big Cypress" &
                           !is.na(larges$Diameter_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend=F)+
  scale_color_viridis(limits=c(1, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_BC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_BC)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(c) Big Cypress") + 
  geom_point(data=subset(larges,larges$Site == "Big Cypress"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Chekika
#Grab parameters associated with future diameter
b0<-p.vec_C[22]
b1<-p.vec_C[23]
b2<-p.vec_C[24]
#Grab parameters associated with future height
b3<-p.vec_C[26]
b4<-p.vec_C[27]
b5<-p.vec_C[28]

df2$diam2_C <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_C <- b3 + b4*df2$Var1 + b5*df2$Var2

#Chekika
plot1_C<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_C)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(d) Chekika") + 
  geom_point(data=subset(larges, larges$Site == "Chekika"&
                           !is.na(larges$Diameter_tplus1)),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend=F)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))


plot2_C<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_C)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(d) Chekika") + 
  geom_point(data=subset(larges,larges$Site == "Chekika"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Fort Pierce
#Grab parameters associated with future diameter
b0<-p.vec_FP[22]
b1<-p.vec_FP[23]
b2<-p.vec_FP[24]
#Grab parameters associated with future height
b3<-p.vec_FP[26]
b4<-p.vec_FP[27]
b5<-p.vec_FP[28]
df2$diam2_FP <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_FP <- b3 + b4*df2$Var1 + b5*df2$Var2

#Fort Pierce
plot1_FP<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_FP)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(e) Fort Pierce") + 
  geom_point(data=subset(larges, larges$Site == "Fort Pierce"&
                           !is.na(larges$Diameter_tplus1)),
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             show.legend=F,shape=1)+
  scale_color_viridis(limits=c(1.6, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_FP<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_FP)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(e) Fort Pierce") + 
  geom_point(data=subset(larges,larges$Site == "Fort Pierce"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1),
             show.legend=F,shape=1 ) +
  scale_color_viridis(limits=c(16, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

##Cape Canaveral
#Grab parameters associated with future diameter
b0 <- p.vec_CC[22]
b1 <- p.vec_CC[23]
b2 <- p.vec_CC[24]

#Grab parameters associated with future height
b3<- p.vec_CC[26]
b4<- p.vec_CC[27]
b5<- p.vec_CC[28]

df2$diam2_CC <- b0 + b1*df2$Var1 + b2*df2$Var2
df2$height2_CC <- b3 + b4*df2$Var1 + b5*df2$Var2

#Cape Canaveral
plot1_CC<- ggplot() + xlim(1.6, 800) + ylim(1.6, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= diam2_CC)) + 
  scale_fill_viridis_c(limits = c(1.6, 950)) +
  ggtitle("(f) Cape Canaveral") + 
  geom_point(data=subset(larges, larges$Site == "Cape Canaveral"&
                           !is.na(larges$Diameter_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Diameter_tplus1), 
             shape=1, show.legend = F)+
  scale_color_viridis(limits=c(1, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order =2))

plot2_CC<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= height2_CC)) + 
  scale_fill_viridis_c(limits = c(16, 950)) +
  ggtitle("(f) Cape Canaveral") + 
  geom_point(data=subset(larges,larges$Site == "Cape Canaveral"&
                           !is.na(larges$Height_tplus1)), 
             aes(x=Diameter_t, y=Height_t, colour = Height_tplus1), 
             show.legend=F, shape=1 ) +
  scale_color_viridis(limits=c(6.5, 800)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "top", order = 2))

g <- ggplotGrob(plot1_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot1_BC2 <- plot1_BC + theme(legend.position ="none")
plot1_CC2 <- plot1_CC + theme(legend.position ="none")
plot1_C2  <- plot1_C  + theme(legend.position = "none")
plot1_FP2 <- plot1_FP + theme(legend.position = "none")
plot1_PG2 <- plot1_PG + theme(legend.position = "none")
plot1_WT2 <- plot1_WT + theme(legend.position = "none")

postscript("./Figures/D2_growth_diameter2.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot1_PG2, plot1_WT2,plot1_BC2, plot1_C2,
             plot1_FP2, plot1_CC2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()


g <- ggplotGrob(plot2_BC + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot2_BC2 <- plot2_BC + theme(legend.position ="none")
plot2_CC2 <- plot2_CC + theme(legend.position ="none")
plot2_C2  <- plot2_C  + theme(legend.position = "none")
plot2_FP2 <- plot2_FP + theme(legend.position = "none")
plot2_PG2 <- plot2_PG + theme(legend.position = "none")
plot2_WT2 <- plot2_WT + theme(legend.position = "none")



postscript("./Figures/D2_growth_height2.eps", width=width.cm/2.54, 
           height=(3*height.cm)/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 7, 7),
             c(1, 1, 2, 2, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(3, 3, 4, 4, 7, 7),
             c(5, 5, 6, 6, 7, 7),
             c(5, 5, 6, 6, 7, 7))
grid.arrange(plot2_PG2, plot2_WT2,plot2_BC2, plot2_C2,
             plot2_FP2, plot2_CC2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()








#FIG A.7 SSD ##### 
y_diam <- c(y1, y3)
y_height <- c(y2, y4)


load("./BC/ss_diam_BC.RData")
load("./CC/ss_diam_CC.RData")
load("./C/ss_diam_C.RData")
load("./FP/ss_diam_FP.RData")
load("./PG/ss_diam_PG.RData")
load("./WT/ss_diam_WT.RData")

load("./BC/ss_height_BC.RData")
load("./CC/ss_height_CC.RData")
load("./C/ss_height_C.RData")
load("./FP/ss_height_FP.RData")
load("./PG/ss_height_PG.RData")
load("./WT/ss_height_WT.RData")

postscript("./Figures/ss.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
ss_diamBC <-c( sum(ss_diam_BC[1:2]), #0-0.32
              sum(ss_diam_BC[3:4]),  #0.32 - 0.64
              sum(ss_diam_BC[5:6]),  #0.64 - 0.96
              sum(ss_diam_BC[7:8]),  #0.96 - 1.28
              sum(ss_diam_BC[9:10]), #1.28- 1.6 (y1[10] + h1/2)
              sum(ss_diam_BC[11:380])) #D2 1.6 - 700
ss_diamCC <-c( sum(ss_diam_CC[1:2]), 
               sum(ss_diam_CC[3:4]), 
               sum(ss_diam_CC[5:6]), 
               sum(ss_diam_CC[7:8]), 
               sum(ss_diam_CC[9:10]), 
               sum(ss_diam_CC[11:380]))
ss_diamC <-c( sum(ss_diam_C[1:2]), 
               sum(ss_diam_C[3:4]), 
               sum(ss_diam_C[5:6]), 
               sum(ss_diam_C[7:8]), 
               sum(ss_diam_C[9:10]), 
               sum(ss_diam_C[11:380]))

ss_diamFP <-c( sum(ss_diam_FP[1:2]), 
               sum(ss_diam_FP[3:4]), 
               sum(ss_diam_FP[5:6]), 
               sum(ss_diam_FP[7:8]), 
               sum(ss_diam_FP[9:10]), 
               sum(ss_diam_FP[11:380]))

ss_diamPG <-c( sum(ss_diam_PG[1:2]), 
               sum(ss_diam_PG[3:4]), 
               sum(ss_diam_PG[5:6]), 
               sum(ss_diam_PG[7:8]), 
               sum(ss_diam_PG[9:10]), 
               sum(ss_diam_PG[11:380]))

ss_diamWT <-c( sum(ss_diam_WT[1:2]), 
               sum(ss_diam_WT[3:4]), 
               sum(ss_diam_WT[5:6]), 
               sum(ss_diam_WT[7:8]), 
               sum(ss_diam_WT[9:10]), 
               sum(ss_diam_WT[11:380]))


ss_heightBC  <- c(sum(ss_height_BC[1:2]), #(0, 2.91) #y2[2] + h2/2
                 sum(ss_height_BC[3:4]), #(2.91, 5.82)
                 sum(ss_height_BC[5:6]), #(5.82, 8.73)
                 sum(ss_height_BC[7:8]), #(8.73, 11.64)
                 sum(ss_height_BC[9:11]),#(11.64, 16)
                 sum(ss_height_BC[12:161])) #(16-800) y4b[100] + h4b/2

ss_heightCC  <- c(sum(ss_height_CC[1:2]),
                   sum(ss_height_CC[3:4]),
                   sum(ss_height_CC[5:6]),
                   sum(ss_height_CC[7:8]),
                   sum(ss_height_CC[9:11]),
                   sum(ss_height_CC[12:161]))

ss_heightC  <- c(sum(ss_height_C[1:2]),
                   sum(ss_height_C[3:4]),
                   sum(ss_height_C[5:6]),
                   sum(ss_height_C[7:8]),
                   sum(ss_height_C[9:11]),
                   sum(ss_height_C[12:161]))

ss_heightFP  <- c(sum(ss_height_FP[1:2]),
                   sum(ss_height_FP[3:4]),
                   sum(ss_height_FP[5:6]),
                   sum(ss_height_FP[7:8]),
                   sum(ss_height_FP[9:11]),
                   sum(ss_height_FP[12:161]))

ss_heightPG  <- c(sum(ss_height_PG[1:2]),
                   sum(ss_height_PG[3:4]),
                   sum(ss_height_PG[5:6]),
                   sum(ss_height_PG[7:8]),
                   sum(ss_height_PG[9:11]),
                   sum(ss_height_PG[12:161]))

ss_heightWT  <- c(sum(ss_height_WT[1:2]),
                   sum(ss_height_WT[3:4]),
                   sum(ss_height_WT[5:6]),
                   sum(ss_height_WT[7:8]),
                   sum(ss_height_WT[9:11]),
                   sum(ss_height_WT[12:161]))
  
site = c(rep("BC",6), rep("C",6), rep("CC", 6), rep("FP", 6), 
         rep("PG", 6), rep("WT", 6))        
size = rep(c("0 - 0.32", "0.32 - 0.64", "0.64 - 0.96", 
             "0.96-1.28", "1.28 - 1.6", "1.6 - 700"), 6)
size2 = rep(c("0 - 2.91", "2.91 - 5.82", "5.82 - 8.73", 
              "8.73 - 11.64", "11.64 - 16", "16 - 800"), 6)

value = c(ss_diamPG, ss_diamWT, ss_diamBC, ss_diamC, ss_diamFP,
          ss_diamCC)
value2 = c(ss_heightPG, ss_heightWT, ss_heightBC, ss_heightC, 
           ss_heightFP,ss_heightCC)

data = data.frame(site, size, value, size2, value2)

levels(data$size2) <- c("0 - 2.91", "2.91 - 5.82", "5.82 - 8.73", 
                        "8.73 - 11.64", "11.64 - 16", "16 - 800")
levels(data$site) <- c("PG", "WT", "BC", "C", "FP", "CC")

# Stacked
p1 <- ggplot(data, aes(fill=size, y=value, x=site)) + 
  geom_bar( stat="identity") +
  scale_fill_brewer(palette = "YlGnBu") +
  labs(fill = "Diameter (mm) ") + 
  labs(x = "Site") + 
  labs(y = "Stable Stage Distribution") +
  ggtitle("(a)") + theme(plot.title = element_text(hjust=0))



p2 <- ggplot(data, aes(fill=size2, y=value2, x=site)) + 
  geom_bar( stat="identity") +
  scale_fill_brewer(palette = "YlGnBu") +
  labs(fill = "Height (cm) ") + 
  labs(x = "Site") + 
  labs(y = "Stable Stage Distribution") +
  ggtitle("(b)")

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/ssd.eps", width=width.cm/2.54, 
           height=2*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
grid.arrange(p1, p2, ncol=1)
dev.off()


#FIG A.8 RV #####
load("./BC/rv_diam_BC.RData")
load("./BC/rv_height_BC.RData")
load("./CC/rv_diam_CC.RData")
load("./CC/rv_height_CC.RData")
load("./C/rv_diam_C.RData")
load("./C/rv_height_C.RData")
load("./FP/rv_diam_FP.RData")
load("./FP/rv_height_FP.RData")
load("./PG/rv_diam_PG.RData")
load("./PG/rv_height_PG.RData")
load("./WT/rv_diam_WT.RData")
load("./WT/rv_height_WT.RData")




rv_diamBC <- c(sum(rv_diam_BC[1:10]),  #0-1.6
              sum(rv_diam_BC[11:245]),  #1.6 - 141.10
              sum(rv_diam_BC[246:288]), # 141.10 - 278.33 
              sum(rv_diam_BC[289: 319]),  # 278.33 - 420.42
              sum(rv_diam_BC[320:350]),  # 420.42 - 562.5
              sum(rv_diam_BC[351:380]))  # 562.5 - 700

rv_diamCC <- c(sum(rv_diam_CC[1:10]),  #0-1.6
               sum(rv_diam_CC[11:245]),  #1.6 - 141.10
               sum(rv_diam_CC[246:288]), # 141.10 - 278.33 
               sum(rv_diam_CC[289: 319]),  # 278.33 - 420.42
               sum(rv_diam_CC[320:350]),  # 420.42 - 562.5
               sum(rv_diam_CC[351:380]))  # 562.5 - 700     

rv_diamC <- c(sum(rv_diam_C[1:10]),  #0-1.6
               sum(rv_diam_C[11:245]),  #1.6 - 141.10
               sum(rv_diam_C[246:288]), # 141.10 - 278.33 
               sum(rv_diam_C[289: 319]),  # 278.33 - 420.42
               sum(rv_diam_C[320:350]),  # 420.42 - 562.5
               sum(rv_diam_C[351:380]))  # 562.5 - 700

rv_diamFP <- c(sum(rv_diam_FP[1:10]),  #0-1.6
               sum(rv_diam_FP[11:245]),  #1.6 - 141.10
               sum(rv_diam_FP[246:288]), # 141.10 - 278.33 
               sum(rv_diam_FP[289: 319]),  # 278.33 - 420.42
               sum(rv_diam_FP[320:350]),  # 420.42 - 562.5
               sum(rv_diam_FP[351:380]))  # 562.5 - 700

rv_diamPG <- c(sum(rv_diam_PG[1:10]),  #0-1.6
               sum(rv_diam_PG[11:245]),  #1.6 - 141.10
               sum(rv_diam_PG[246:288]), # 141.10 - 278.33 
               sum(rv_diam_PG[289: 319]),  # 278.33 - 420.42
               sum(rv_diam_PG[320:350]),  # 420.42 - 562.5
               sum(rv_diam_PG[351:380]))  # 562.5 - 700

rv_diamWT <- c(sum(rv_diam_WT[1:10]),  #0-1.6
               sum(rv_diam_WT[11:245]),  #1.6 - 141.10
               sum(rv_diam_WT[246:288]), # 141.10 - 278.33 
               sum(rv_diam_WT[289: 319]),  # 278.33 - 420.42
               sum(rv_diam_WT[320:350]),  # 420.42 - 562.5
               sum(rv_diam_WT[351:380]))  # 562.5 - 700

rv_heightBC <- c(sum(rv_height_BC[1:11]),  #0-16
               sum(rv_height_BC[12:39]),  #16 - 175
               sum(rv_height_BC[40:67]), # 175 - 330 
               sum(rv_height_BC[68:98]),  # 330 - 485
               sum(rv_height_BC[99:130]),  # 485 - 645
               sum(rv_height_BC[131:161]))  #645-800

rv_heightCC <- c(sum(rv_height_CC[1:11]),  #0-16
                 sum(rv_height_CC[12:39]),  #16 - 175
                 sum(rv_height_CC[40:67]), # 175 - 330 
                 sum(rv_height_CC[68:98]),  # 330 - 485
                 sum(rv_height_CC[99:130]),  # 485 - 645
                 sum(rv_height_CC[131:161]))  #645-800

rv_heightC <- c(sum(rv_height_C[1:11]),  #0-16
                 sum(rv_height_C[12:39]),  #16 - 175
                 sum(rv_height_C[40:67]), # 175 - 330 
                 sum(rv_height_C[68:98]),  # 330 - 485
                 sum(rv_height_C[99:130]),  # 485 - 645
                 sum(rv_height_C[131:161]))  #645-800

rv_heightFP <- c(sum(rv_height_FP[1:11]),  #0-16
                 sum(rv_height_FP[12:39]),  #16 - 175
                 sum(rv_height_FP[40:67]), # 175 - 330 
                 sum(rv_height_FP[68:98]),  # 330 - 485
                 sum(rv_height_FP[99:130]),  # 485 - 645
                 sum(rv_height_FP[131:161]))  #645-800

rv_heightPG <- c(sum(rv_height_PG[1:11]),  #0-16
                 sum(rv_height_PG[12:39]),  #16 - 175
                 sum(rv_height_PG[40:67]), # 175 - 330 
                 sum(rv_height_PG[68:98]),  # 330 - 485
                 sum(rv_height_PG[99:130]),  # 485 - 645
                 sum(rv_height_PG[131:161]))  #645-800

rv_heightWT <- c(sum(rv_height_WT[1:11]),  #0-16
                 sum(rv_height_WT[12:39]),  #16 - 175
                 sum(rv_height_WT[40:67]), # 175 - 330 
                 sum(rv_height_WT[68:98]),  # 330 - 485
                 sum(rv_height_WT[99:130]),  # 485 - 645
                 sum(rv_height_WT[131:161]))  #645-800

site = c(rep("BC",6), rep("C",6), rep("CC", 6), rep("FP", 6), 
         rep("PG", 6),rep("WT", 6))        
size = rep(c("0 - 1.6", "1.6-141", "141-278", 
             "278 - 420", "420 - 563", "563 - 700"), 6)
size2 = rep(c("0 - 16", "16 - 175", "175 - 330", 
              "330 - 485 ", "485 - 645", "645 - 800"), 6)


value = c(rv_diamPG, rv_diamWT, rv_diamBC, rv_diamC, rv_diamFP,
          rv_diamCC)
value2 = c(rv_heightPG, rv_heightWT, rv_heightBC, rv_heightC, 
           rv_heightFP,rv_heightCC)

data = data.frame(site, size, value, size2, value2)
levels(data$site) <- c("PG", "WT", "BC", "C", "FP", "CC")




# Stacked
p1 <- ggplot(data, aes(fill=size, y=value, x=site)) + 
  geom_bar( stat="identity") +
  scale_fill_brewer(palette = "RdYlBu", direction= -1) +
  labs(fill = "Diameter (mm) ") + 
  labs(x = "Site") + 
  labs(y = "Reproductive Value") +
  ggtitle("(a)") + theme(plot.title = element_text(hjust=0))

p2 <- ggplot(data, aes(fill=size2, y=value2, x=site)) + 
  geom_bar( stat="identity") +
  scale_fill_brewer(palette = "RdYlBu", direction= -1) +
  labs(fill = "Height (cm) ") + 
  labs(x = "Site") + 
  labs(y = "Reproductive Value") +
  ggtitle("(b)") + theme(plot.title = element_text(hjust=0))


setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/rv.eps", width=width.cm/2.54, 
           height=2*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
grid.arrange(p1, p2, ncol=1)
dev.off()


## FIG A.9 MARGINAL ELASTICITY #####
#Big Cypress
D1_elas_diamBC<- sum(rowSums(total.elas_D1_BC))
G_elas_diamBC <- sum(rowSums(total.elas_GA_BC),
                     rowSums(total.elas_GB_BC))
F_elas_diamBC <- sum(rowSums(total.elas_FA_BC), 
                      rowSums(total.elas_FB_BC))
D2_elas_diamBC <- sum(rowSums(total.elas_D2AA_BC),
                       rowSums(total.elas_D2AB_BC),
                       rowSums(total.elas_D2BA_BC),
                       rowSums(total.elas_D2BB_BC))
elas_diamBC <- c(D1_elas_diamBC, G_elas_diamBC, F_elas_diamBC,
                 D2_elas_diamBC)

D1_elas_heightBC<- sum(colSums(total.elas_D1_BC))
G_elas_heightBC <- sum(colSums(total.elas_GA_BC),
                       colSums(total.elas_GB_BC))
F_elas_heightBC <- sum(colSums(total.elas_FA_BC), 
                       colSums(total.elas_FB_BC))
D2_elas_heightBC <- sum(colSums(total.elas_D2AA_BC),
                        colSums(total.elas_D2AB_BC),
                        colSums(total.elas_D2BA_BC),
                        colSums(total.elas_D2BB_BC))
elas_heightBC <- c(D1_elas_heightBC, G_elas_heightBC, F_elas_heightBC,
                 D2_elas_heightBC)

#Cape Canaveral
D1_elas_diamCC<- sum(rowSums(total.elas_D1_CC))
G_elas_diamCC <- sum(rowSums(total.elas_GA_CC),
                     rowSums(total.elas_GB_CC))
F_elas_diamCC <- sum(rowSums(total.elas_FA_CC), 
                     rowSums(total.elas_FB_CC))
D2_elas_diamCC <- sum(rowSums(total.elas_D2AA_CC),
                      rowSums(total.elas_D2AB_CC),
                      rowSums(total.elas_D2BA_CC),
                      rowSums(total.elas_D2BB_CC))
elas_diamCC <- c(D1_elas_diamCC, G_elas_diamCC, F_elas_diamCC,
                 D2_elas_diamCC)

D1_elas_heightCC<- sum(colSums(total.elas_D1_CC))
G_elas_heightCC <- sum(colSums(total.elas_GA_CC),
                       colSums(total.elas_GB_CC))
F_elas_heightCC <- sum(colSums(total.elas_FA_CC), 
                       colSums(total.elas_FB_CC))
D2_elas_heightCC <- sum(colSums(total.elas_D2AA_CC),
                        colSums(total.elas_D2AB_CC),
                        colSums(total.elas_D2BA_CC),
                        colSums(total.elas_D2BB_CC))
elas_heightCC <- c(D1_elas_heightCC, G_elas_heightCC, F_elas_heightCC,
                   D2_elas_heightCC)

#Chekika
D1_elas_diamC<- sum(rowSums(total.elas_D1_C))
G_elas_diamC <- sum(rowSums(total.elas_GA_C),
                     rowSums(total.elas_GB_C))
F_elas_diamC <- sum(rowSums(total.elas_FA_C), 
                     rowSums(total.elas_FB_C))
D2_elas_diamC <- sum(rowSums(total.elas_D2AA_C),
                      rowSums(total.elas_D2AB_C),
                      rowSums(total.elas_D2BA_C),
                      rowSums(total.elas_D2BB_C))
elas_diamC <- c(D1_elas_diamC, G_elas_diamC, F_elas_diamC,
                 D2_elas_diamC)

D1_elas_heightC<- sum(colSums(total.elas_D1_C))
G_elas_heightC <- sum(colSums(total.elas_GA_C),
                       colSums(total.elas_GB_C))
F_elas_heightC <- sum(colSums(total.elas_FA_C), 
                       colSums(total.elas_FB_C))
D2_elas_heightC <- sum(colSums(total.elas_D2AA_C),
                        colSums(total.elas_D2AB_C),
                        colSums(total.elas_D2BA_C),
                        colSums(total.elas_D2BB_C))
elas_heightC <- c(D1_elas_heightC, G_elas_heightC, F_elas_heightC,
                   D2_elas_heightC)

#Fort Pierce
D1_elas_diamFP<- sum(rowSums(total.elas_D1_FP))
G_elas_diamFP <- sum(rowSums(total.elas_GA_FP),
                     rowSums(total.elas_GB_FP))
F_elas_diamFP <- sum(rowSums(total.elas_FA_FP), 
                     rowSums(total.elas_FB_FP))
D2_elas_diamFP <- sum(rowSums(total.elas_D2AA_FP),
                      rowSums(total.elas_D2AB_FP),
                      rowSums(total.elas_D2BA_FP),
                      rowSums(total.elas_D2BB_FP))
elas_diamFP <- c(D1_elas_diamFP, G_elas_diamFP, F_elas_diamFP,
                 D2_elas_diamFP)

D1_elas_heightFP<- sum(colSums(total.elas_D1_FP))
G_elas_heightFP <- sum(colSums(total.elas_GA_FP),
                       colSums(total.elas_GB_FP))
F_elas_heightFP <- sum(colSums(total.elas_FA_FP), 
                       colSums(total.elas_FB_FP))
D2_elas_heightFP <- sum(colSums(total.elas_D2AA_FP),
                        colSums(total.elas_D2AB_FP),
                        colSums(total.elas_D2BA_FP),
                        colSums(total.elas_D2BB_FP))
elas_heightFP <- c(D1_elas_heightFP, G_elas_heightFP, F_elas_heightFP,
                   D2_elas_heightFP)

#Punta Gorda
D1_elas_diamPG<- sum(rowSums(total.elas_D1_PG))
G_elas_diamPG <- sum(rowSums(total.elas_GA_PG),
                     rowSums(total.elas_GB_PG))
F_elas_diamPG <- sum(rowSums(total.elas_FA_PG), 
                     rowSums(total.elas_FB_PG))
D2_elas_diamPG <- sum(rowSums(total.elas_D2AA_PG),
                      rowSums(total.elas_D2AB_PG),
                      rowSums(total.elas_D2BA_PG),
                      rowSums(total.elas_D2BB_PG))
elas_diamPG <- c(D1_elas_diamPG, G_elas_diamPG, F_elas_diamPG,
                 D2_elas_diamPG)

D1_elas_heightPG<- sum(colSums(total.elas_D1_PG))
G_elas_heightPG <- sum(colSums(total.elas_GA_PG),
                       colSums(total.elas_GB_PG))
F_elas_heightPG <- sum(colSums(total.elas_FA_PG), 
                       colSums(total.elas_FB_PG))
D2_elas_heightPG <- sum(colSums(total.elas_D2AA_PG),
                        colSums(total.elas_D2AB_PG),
                        colSums(total.elas_D2BA_PG),
                        colSums(total.elas_D2BB_PG))
elas_heightPG <- c(D1_elas_heightPG, G_elas_heightPG, F_elas_heightPG,
                   D2_elas_heightPG)

#Wild Turkey
D1_elas_diamWT<- sum(rowSums(total.elas_D1_WT))
G_elas_diamWT <- sum(rowSums(total.elas_GA_WT),
                     rowSums(total.elas_GB_WT))
F_elas_diamWT <- sum(rowSums(total.elas_FA_WT), 
                     rowSums(total.elas_FB_WT))
D2_elas_diamWT <- sum(rowSums(total.elas_D2AA_WT),
                      rowSums(total.elas_D2AB_WT),
                      rowSums(total.elas_D2BA_WT),
                      rowSums(total.elas_D2BB_WT))
elas_diamWT <- c(D1_elas_diamWT, G_elas_diamWT, F_elas_diamWT,
                 D2_elas_diamWT)

D1_elas_heightWT<- sum(colSums(total.elas_D1_WT))
G_elas_heightWT <- sum(colSums(total.elas_GA_WT),
                       colSums(total.elas_GB_WT))
F_elas_heightWT <- sum(colSums(total.elas_FA_WT), 
                       colSums(total.elas_FB_WT))
D2_elas_heightWT <- sum(colSums(total.elas_D2AA_WT),
                        colSums(total.elas_D2AB_WT),
                        colSums(total.elas_D2BA_WT),
                        colSums(total.elas_D2BB_WT))
elas_heightWT <- c(D1_elas_heightWT, G_elas_heightWT, F_elas_heightWT,
                   D2_elas_heightWT)


site = c(rep("BC",4), rep("C",4), rep("CC", 4), rep("FP", 4), 
         rep("PG", 4), rep("WT", 4))  

domain = rep(c("Seedling", "Maturation", "Fertility", "Adult"), 6)
levels(domain) <- c("Seedling", "Maturation", "Fertility", "Adult")

value = c(elas_diamPG, elas_diamWT, elas_diamBC, elas_diamC,
          elas_diamFP,elas_diamCC)
value2 = c(elas_heightPG, elas_heightWT, elas_heightBC,
           elas_heightC,elas_heightFP, elas_heightCC)

data = data.frame(site, domain, value, value2)

levels(data$site) <- c("PG", "WT", "BC", "C", "FP", "CC")


# Stacked

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/total_elasticity.eps", width=width.cm/2.54, 
           height=2*height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
ggplot(data, aes(fill=domain, y=value, x=site)) + 
  geom_bar( stat="identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = "") + 
  labs(x = "Site") + 
  labs(y = "Elasticity") 
dev.off()

##### FIG E.1 Randomization #####
lam.stable_BC <- readRDS("./BC/lam.stable_BC.rds")
lam.stable_CC <- readRDS("./CC/lam.stable_CC.rds")
lam.stable_C  <- readRDS("./C/lam.stable_C.rds")
lam.stable_FP <- readRDS("./FP/lam.stable_FP.rds")
lam.stable_PG <- readRDS("./PG/lam.stable_PG.rds")
lam.stable_WT <- readRDS("./WT/lam.stable_WT.rds")

PG_WT_obs <- abs(lam.stable_PG - lam.stable_WT)
PG_BC_obs <- abs(lam.stable_PG - lam.stable_BC)
PG_C_obs  <- abs(lam.stable_PG - lam.stable_C)
PG_FP_obs <- abs(lam.stable_PG - lam.stable_FP)
PG_CC_obs <- abs(lam.stable_PG - lam.stable_CC)

WT_BC_obs <- abs(lam.stable_WT - lam.stable_BC)
WT_C_obs  <- abs(lam.stable_WT - lam.stable_C)
WT_FP_obs <- abs(lam.stable_WT - lam.stable_FP)
WT_CC_obs <- abs(lam.stable_WT - lam.stable_CC)

BC_C_obs  <- abs(lam.stable_BC - lam.stable_C)
BC_FP_obs <- abs(lam.stable_BC - lam.stable_FP)
BC_CC_obs <- abs(lam.stable_BC - lam.stable_CC)

C_FP_obs  <- abs(lam.stable_C  - lam.stable_FP)
C_CC_obs  <- abs(lam.stable_C  - lam.stable_CC)

FP_CC_obs <- abs(lam.stable_FP - lam.stable_CC)
lam.boot <- readRDS("lam.boot_271.rds")
#lam.boot
colnames(lam.boot) <- c("BC", "CC", "C", "FP", "PG", "WT")
differences <- data.frame(PG_WT = rep(NA, nrow(lam.boot)),
                          PG_BC  = rep(NA, nrow(lam.boot)),
                          PG_C = rep(NA, nrow(lam.boot)),
                          PG_FP = rep(NA, nrow(lam.boot)),
                          PG_CC = rep(NA, nrow(lam.boot)),
                          WT_BC  = rep(NA, nrow(lam.boot)),
                          WT_C = rep(NA, nrow(lam.boot)),
                          WT_FP = rep(NA, nrow(lam.boot)),
                          WT_CC = rep(NA, nrow(lam.boot)),
                          BC_C  = rep(NA, nrow(lam.boot)),
                          BC_FP  = rep(NA, nrow(lam.boot)),
                          BC_CC  = rep(NA, nrow(lam.boot)),
                          C_FP = rep(NA, nrow(lam.boot)),
                          C_CC = rep(NA, nrow(lam.boot)),
                          FP_CC = rep(NA, nrow(lam.boot)))

differences$PG_WT <- abs(lam.boot[,5]-lam.boot[,6])
differences$PG_BC <- abs(lam.boot[,5]-lam.boot[,1])
differences$PG_C  <- abs(lam.boot[,5]-lam.boot[,3])
differences$PG_FP <- abs(lam.boot[,5]-lam.boot[,4])
differences$PG_CC <- abs(lam.boot[,5]-lam.boot[,2])

differences$WT_BC <- abs(lam.boot[,6]-lam.boot[,1])
differences$WT_C  <- abs(lam.boot[,6]-lam.boot[,3])
differences$WT_FP <- abs(lam.boot[,6]-lam.boot[,4])
differences$WT_CC <- abs(lam.boot[,6]-lam.boot[,2])

differences$BC_C   <- abs(lam.boot[,1]-lam.boot[,3])
differences$BC_FP  <- abs(lam.boot[,1]-lam.boot[,4])
differences$BC_CC  <- abs(lam.boot[,1]-lam.boot[,2])

differences$C_FP  <- abs(lam.boot[,3]-lam.boot[,4])
differences$C_CC  <- abs(lam.boot[,3]-lam.boot[,2])

differences$FP_CC <- abs(lam.boot[,4]-lam.boot[,2])

#Figure out what xlim should be:
max(summary(differences[,1]), summary(differences[,2]),
    summary(differences[,3]), summary(differences[,4]), 
    summary(differences[,5]), summary(differences[,6]),
    summary(differences[,7]), summary(differences[,8]),
    summary(differences[,9]), summary(differences[,10]),
    summary(differences[,11]), summary(differences[,12]),
    summary(differences[,13]), summary(differences[,14]),
    summary(differences[,15]), PG_WT_obs, PG_BC_obs, PG_C_obs, 
    PG_FP_obs, PG_CC_obs,WT_BC_obs, WT_C_obs, WT_FP_obs, 
    WT_CC_obs, BC_C_obs, BC_FP_obs, BC_CC_obs, C_FP_obs,
    C_CC_obs, FP_CC_obs)

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/randomization.eps", width=width.cm/2.54, 
           height=1.5*width.cm/(2.54), pointsize=pointsize,  encoding = "TeXtext.enc")
#x11(width = width.cm/2.54, height = 1.5*width.cm/(2.54), 
#pointsize = pointsize)
#png(file="./Figures/variability.png")
par(mar = c(3, 3, 2, 3), # Margins
    mgp = c(1.5, .5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    oma=c(0, 1, 4,3.5),
    mai=c(0.3,0.4,0.05,0.1),
    mfrow=c(5, 5))

hist(differences[,1], main = "", xlab="Difference", 
     xlim=c(0, 0.2), breaks=12)
abline(v=PG_WT_obs, col="red")
mtext(side=3, "Wild Turkey", line =2, cex = 1.2)
mtext(side=2, "Punta Gorda", line = 3, cex=1.2)

hist(differences[,2], main = "", xlab="Difference",
     breaks=12, xlim=c(0, 0.2))
abline(v=PG_BC_obs, col="red")
mtext(side=3, "Big Cypress", line =2, cex=1.2)

hist(differences[,3], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=PG_C_obs, col="red")
mtext(side=3, "Chekika", line =2, cex=1.2)

hist(differences[,4], main= "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=PG_FP_obs, col="red")
mtext(side=3, "Fort Pierce", line = 2, cex=1.2)

hist(differences[,5], main="", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=PG_CC_obs, col="red")
mtext(side=3, "Cape Canaveral", line = 2, cex = 1.2)

plot(0,type='n',axes=FALSE,ann=FALSE)
mtext(side=2, "Wild Turkey", line =3, cex=1.2)

hist(differences[,6], main="", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=WT_BC_obs, col="red")

hist(differences[,7], main ="", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=WT_C_obs, col="red")

hist(differences[,8], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=WT_FP_obs, col="red")

hist(differences[,9], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=WT_CC_obs, col="red")

plot(0,type='n',axes=FALSE,ann=FALSE)
mtext(side=2, "Big Cypress", line =3, cex=1.2)
plot(0,type='n',axes=FALSE,ann=FALSE)

hist(differences[,10], main="", xlab="Difference", 
     breaks=12, xlim=c(0,0.2))
abline(v=BC_C_obs, col="red")

hist(differences[,11], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=BC_FP_obs, col="red")

hist(differences[,12], main = "", xlab="Difference", 
     breaks=12, xlim=c(0,0.2))
abline(v=BC_CC_obs, col="red")

plot(0,type='n',axes=FALSE,ann=FALSE)
mtext(side=2, "Chekika", line =3, cex=1.2)
plot(0,type='n',axes=FALSE,ann=FALSE)
plot(0,type='n',axes=FALSE,ann=FALSE)

hist(differences[,13], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=C_FP_obs, col="red")

hist(differences[,14], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=C_CC_obs, col="red")

plot(0,type='n',axes=FALSE,ann=FALSE)
mtext(side=2, "Fort Pierce", line =3, cex=1.2)
plot(0,type='n',axes=FALSE,ann=FALSE)
plot(0,type='n',axes=FALSE,ann=FALSE)
plot(0,type='n',axes=FALSE,ann=FALSE)
hist(differences[,15], main = "", xlab="Difference",
     breaks=12, xlim=c(0,0.2))
abline(v=FP_CC_obs, col="red")
dev.off()

<<<<<<< HEAD
##FIG F.1 TAU #####
load("results_tau_sensitivity.RData")

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/tau.eps", width=width.cm_onepanel/2.54, 
           height=width.cm_onepanel/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
# x11(width = width.cm/2.54, height = height.cm/2.54, 
#     pointsize = pointsize)
=======
#*- Sensitivity analysis of tau#####

setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/randomization.eps", width=width.cm/2.54, 
           height=1.5*width.cm/(2.54), pointsize=pointsize,  encoding = "TeXtext.enc")
x11(width = width.cm/2.54, height = width.cm/(2.54), 
pointsize = pointsize)

load("results_tau_sensitivity.RData")
>>>>>>> 3e0678db558215c0c2ffe614aa1e814956e9ac49

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,0.5,0.5,0.5),
    mfrow=c(1, 1))
plot(log10(results_tau_sensitivity$newparm), 
     results_tau_sensitivity$lambda_PG, col=cols[1], pch=19,
     ylim=c(0.95, 1.75), ylab=expression(lambda),
     xlab=expression(paste("log(", tau[1], ")")))
points(log10(results_tau_sensitivity$newparm),
       results_tau_sensitivity$lambda_WT, col=cols[2], pch=19)
points(log10(results_tau_sensitivity$newparm), 
       results_tau_sensitivity$lambda_BC, col=cols[3], pch=19)
points(log10(results_tau_sensitivity$newparm),
       results_tau_sensitivity$lambda_C, col=cols[4], pch=19)
points(log10(results_tau_sensitivity$newparm),
       results_tau_sensitivity$lambda_FP, col=cols[5], pch=19)
points(log10(results_tau_sensitivity$newparm),
       results_tau_sensitivity$lambda_CC, col=cols[6], pch=19)

legend(-3.75, 1.75, pch=19, col=cols,
       c("PG", "WT", "BC", "C", "FP", "CC"))

points(log10(p.vec_PG[37]), lam.stable_PG, col=cols[1], pch=4)
points(log10(p.vec_WT[37]), lam.stable_WT, col=cols[2], pch=4)
points(log10(p.vec_BC[37]), lam.stable_BC, col=cols[3], pch=4)
points(log10(p.vec_C[37]),  lam.stable_C,  col=cols[4], pch=4)
points(log10(p.vec_FP[37]), lam.stable_FP, col=cols[5], pch=4)
points(log10(p.vec_CC[37]), lam.stable_CC, col=cols[6], pch=4)
dev.off()