# Code to produce figures in Erickson et al. integral projection manuscript. 
# Uses code provided by Ellner and Rees 
# Modifications to Ellner and Rees scripts added by Kelley D. Erickson

library(RColorBrewer)
library(viridis)
library(ggplot2)

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
x1seq_larges <- seq(1.6, 800, length.out = 50 )
x2seq_larges <- seq (16, 800, length.out = 50)


# FIGURE 1: D1 Survival by Biotype #####
load("./BC/p.vec_BC.RData")
load("./CC/p.vec_CC.RData")
load("./C/p.vec_C.RData")
load("./FP/p.vec_FP.RData")
load("./PG/p.vec_PG.RData")
load("./WT/p.vec_WT.RData")

#png(file="D1_survival.png", width=20, height=10, units="in", res=300)
#setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/D1_survival.eps", width=width.cm/2.54, 
          height=(height.cm*2)/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
 #x11(width = width.cm, height = 2*height.cm, 
   #pointsize = pointsize)
par(mfrow=c(2,2)) 
#par(mar = c(7.1,15.1,5.1,5.1))
par(mar=c(7, 11, 2, 1.5 ))
par(oma=c(0,0,1,4)) #2, .5, 2, 4
#Big Cypress, Punta Gorda and Wild Turkey
b0 <- p.vec_BC[1]
b1 <- p.vec_BC[2]
b2 <- p.vec_BC[3]
df$BC_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)
z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) (invlogit(b0 + b1*a + b2*b)))

par(cex.axis = 1, cex.lab = 1.2, cex = 1.2, cex.sub = 2, 
    cex.main = 1.2, lwd = 1.5, bg = "transparent")
image.plot(x1seq_seedlings, x2seq_seedlings, z1, xlab = "", 
           ylab="", main="", col=my.palette, legend.line = 2.5, legend.lab= "P(Survival)", 
           legend.cex=1.2, zlim= c(0,1))
points(seedlings$Diameter_t[seedlings$Surv_tplus1==1 & (seedlings$Site=="Big Cypress" |seedlings$Site == "Punta Gorda" | seedlings$Site == "Wild Turkey")], seedlings$Height_t[seedlings$Surv_tplus1==1& (seedlings$Site=="Big Cypress" |seedlings$Site == "Punta Gorda" | seedlings$Site == "Wild Turkey")],
       pch=19, cex=1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==0 & (seedlings$Site=="Big Cypress" |seedlings$Site == "Punta Gorda" | seedlings$Site == "Wild Turkey")], seedlings$Height_t[seedlings$Surv_tplus1==0&(seedlings$Site=="Big Cypress" |seedlings$Site == "Punta Gorda" | seedlings$Site == "Wild Turkey")],
       pch=1, cex=1.2)
mtext(side=3, "(a) BC, PG, WT", cex = 1.2, line = 1)
mtext(side=1, "Diameter (mm)", line =3, cex=1.2)
mtext(side=2, "Height (cm)", line = 3, cex = 1.2)
#legend(1.63, 15, pch=c(1,19), cex = 0.8, text.width=0.19, c("Died", "Survived"), xpd=T)


#Cape Canaveral
b0 <- p.vec_CC[1]
b1 <- p.vec_CC[2]
b2 <- p.vec_CC[3]
df$CC_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) (invlogit(b0 + b1*a + b2*b)))

par(cex.axis = 1, cex.lab = 1.2, cex = 1.2, cex.sub = 2, 
    cex.main = 1.2, lwd = 1.5, bg = "transparent")
image.plot(x1seq_seedlings, x2seq_seedlings, z1, xlab = "", ylab="", 
           legend.line = 2.5, legend.lab= "P(Survival)", zlim=c(0,1),
           main="", col=my.palette, legend.cex=1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==1 & seedlings$Site=="Cape Canaveral"], seedlings$Height_t[seedlings$Surv_tplus1==1& seedlings$Site=="Cape Canaveral"], 
       pch=19, cex=1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==0 & seedlings$Site=="Cape Canaveral"], seedlings$Height_t[seedlings$Surv_tplus1==0&seedlings$Site=="Cape Canaveral"], 
       pch=1, cex=1.2)
mtext(side=3, "(b) CC", cex = 1.2, line = 1)
mtext(side=1, "Diameter (mm)", line =3, cex=1.2)
mtext(side=2, "Height (cm)", line = 3, cex = 1.2)
#legend(1.63, 15, pch=c(1,19), cex = 0.8, text.width=0.19, c("Died", "Survived"), xpd=T)

#Chekika
b0 <- p.vec_C[1]
b1 <- p.vec_C[2]
b2 <- p.vec_C[3]
df$C_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) (invlogit(b0 + b1*a + b2*b)))

par(cex.axis = 1, cex.lab = 1.2, cex = 1.2, cex.sub = 2, 
    cex.main = 1.2, lwd = 1.5, bg = "transparent")
image.plot(x1seq_seedlings, x2seq_seedlings, z1, xlab = "", ylab="", 
            main="", col=my.palette, zlim=c(0,1), legend.line = 2.5, 
           legend.lab= "P(Survival)", legend.cex=1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==1 & seedlings$Site=="Chekika"], seedlings$Height_t[seedlings$Surv_tplus1==1& seedlings$Site=="Chekika"],
       pch=19, cex = 1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==0 & seedlings$Site=="Chekika"], seedlings$Height_t[seedlings$Surv_tplus1==0&seedlings$Site=="Chekika"],
       pch=1, cex = 1.2)
mtext(side=3, "(c) C", cex = 1.2, line = 1)
mtext(side=1, "Diameter (mm)", line =3, cex=1.2)
mtext(side=2, "Height (cm)", line = 3, cex = 1.2)
#legend(1.63, 15, pch=c(1,19), cex = 0.8, text.width=0.19, c("Died", "Survived"), xpd=T)

#Fort Pierce
b0 <- p.vec_FP[1]
b1 <- p.vec_FP[2]
b2 <- p.vec_FP[3]
df$FP_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) (invlogit(b0 + b1*a + b2*b)))

par(cex.axis = 1, cex.lab = 1.2, cex = 1.2, cex.sub = 2, 
    cex.main = 1.2, lwd = 1.5, bg = "transparent")
image.plot(x1seq_seedlings, x2seq_seedlings, z1, xlab = "", ylab="",
           main="", col=my.palette, legend.line = 2.5, zlim=c(0,1),
           legend.lab= "P(Survival)", legend.cex=1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==1 & seedlings$Site=="Fort Pierce"], seedlings$Height_t[seedlings$Surv_tplus1==1& seedlings$Site=="Fort Pierce"], 
       pch=19, cex=1.2)
points(seedlings$Diameter_t[seedlings$Surv_tplus1==0 & seedlings$Site=="Fort Pierce"], seedlings$Height_t[seedlings$Surv_tplus1==0&seedlings$Site=="Fort Pierce"],
       pch=1, cex = 1.2)
mtext(side=3, "(d) FP", cex = 1.2, line = 1)
mtext(side=1, "Diameter (mm)", line =3, cex=1.2)
mtext(side=2, "Height (cm)", line = 3, cex = 1.2)
#legend(1.63, 15, pch=c(1,19), cex = 0.8, text.width=0.19, c("Died", "Survived"), xpd=T)
#This creates an EPS figure! That LaTeX can read

dev.off()

#*---Figure 1 in ggplot #####
#Big Cypress, Punta Gorda, Wild Turkey
plot_BCs<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= BC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(a) BC, PG, WT") + 
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

#Cape Canaveral
plot_CCs<-ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= CC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(b) CC") + 
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

#Chekika
plot_Cs <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= C_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(c) C ") + 
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

plot_FPs <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= FP_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  ggtitle("(d) FP") + 
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
grid.arrange(plot_BCs2, plot_CCs2,plot_Cs2, plot_FPs2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()




# FIGURE 2 D1 GROWTH  #####
df <- expand.grid(x1seq_seedlings, x2seq_seedlings)

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


#x11(width = width.cm/2.54, height = height.cm/2.54, 
#pointsize = pointsize)




require(gridExtra)
library(grid)

#*-- ggplot #####
#Big Cypress
plot1_BC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_BC)) +
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(a) Big Cypress") + 
 geom_point(data=subset(seedlings2, seedlings2$Site == "Big Cypress"), 
            aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5),
                        limits=c(0.3, 1.5), range= c(0,3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_BC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_BC)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(a) Big Cypress") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Big Cypress"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(breaks = c(3, 6, 9, 12, 15), limits=c(3, 15), 
                        range = c(0, 3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Cape Canaveral
plot1_CC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_CC)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(b) Cape Canaveral") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Cape Canaveral"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5), 
                        limits=c(0.3, 1.5), range= c(0,3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_CC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_CC)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(b) Cape Canaveral") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Cape Canaveral"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(breaks = c(3, 6, 9, 12, 15), limits=c(3, 15),
                        range = c(0, 3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Chekika
plot1_C<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_C)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(c) Chekika") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Chekika"),
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5), 
                        limits=c(0.3, 1.5), range= c(0,3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


plot2_C<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_C)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(c) Chekika") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Chekika"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(breaks = c(3, 6, 9, 12, 15), 
                        limits=c(3, 15), range = c(0, 3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Fort Pierce
plot1_FP<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_FP)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(d) Fort Pierce") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Fort Pierce"),
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5), 
                        limits=c(0.3, 1.5), range= c(0,3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_FP<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_FP)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(d) Fort Pierce") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Fort Pierce"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(breaks = c(3, 6, 9, 12, 15), 
                        limits=c(3, 15), range = c(0, 3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Punta Gorda
plot1_PG<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_PG)) + 
  scale_fill_viridis_c(limits = c(0, 1.7)) +
  ggtitle("(e) Punta Gorda") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Punta Gorda"),
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5), 
                        limits=c(0.3, 1.5), range= c(0,3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_PG<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_PG)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(e) Punta Gorda") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Punta Gorda"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(breaks = c(3, 6, 9, 12, 15), 
                        limits=c(3, 15), range = c(0, 3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n height \n at time \n t + 1 (cm)") +
  labs(size = "Measured \n height \n at time \n t + 1 (cm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order = 1)) +
  guides(fill = guide_colorbar(title.position = "right", order = 2))

#Wild Turkey
plot1_WT<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= diam_WT)) + 
  scale_fill_viridis_c(limits = c(0, 1.6)) +
  ggtitle("(f) Wild Turkey") + 
  geom_point(data=subset(seedlings2, seedlings2$Site == "Wild Turkey"), 
             aes(x=Diameter_t, y=Height_t, size = Diameter_tplus1), shape=1)+
  scale_size_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5), 
                        limits=c(0.3, 1.5), range= c(0,3)) + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=8), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n diameter \n at time \n t + 1 (mm)") +
  labs(size = "Measured \n diameter \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(size = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot2_WT<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= height_WT)) + 
  scale_fill_viridis_c(limits = c(0, 16)) +
  ggtitle("(f) Wild Turkey") + 
  geom_point(data=subset(seedlings2,seedlings2$Site == "Wild Turkey"), 
             aes(x=Diameter_t, y=Height_t, size = Height_tplus1), shape=1 ) +
  scale_size_continuous(breaks = c(3, 6, 9, 12, 15), 
                        limits=c(3, 15), range = c(0, 3)) + 
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

postscript("./Figures/D1_growth_diameter.eps", width=width.cm/2.54, 
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



postscript("./Figures/D1_growth_height.eps", width=width.cm/2.54, 
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


#setEPS(horizontal=F, onefile=F, paper="special")
#postscript("./Figures/D1_growth.eps", width=width.cm/2.54, 
           #height=height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")
 x11(width = width.cm/2.54, height = 2*height.cm/2.54, 
     pointsize = pointsize)##### Image Plots
par(mfrow=c(1,2))
my.palette <- viridis(9)
par(mar=c(5.1,4.1,4.1,6))
coeff_bigger=.5

layout(matrix(c(1,2, 3, 4), 2, 2, byrow = TRUE), 
       widths=c(2,2), heights=c(3, 1))

#Big Cypress IMAGE PLOT
#Future Diameter

image.plot(x1seq_seedlings, x2seq_seedlings, z1_BC, xlab = "Diameter (mm) ", 
           ylab="Height (cm)", col=my.palette, legend.lab = "Diameter at t+1 (mm)", 
           main = "Big Cypress")
points(seedlings2$Diameter_t[seedlings2$Site == "Big Cypress"], 
       seedlings2$Height_t[seedlings2$Site == "Big Cypress"], 
       cex=seedlings2$Diameter_tplus1[seedlings2$Site == "Big Cypress"]*coeff_bigger, 
       pch=1)
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")

#Future Height
image.plot(x1seq_seedlings, x2seq_seedlings, z2_BC, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, main ="Big Cypress",legend.lab = "Height at t+1 (cm)")
points(seedlings2$Diameter_t[seedlings2$Site == "Big Cypress"], seedlings2$Height_t[seedlings2$Site == "Big Cypress"], cex=seedlings2$Height_tplus1[seedlings2$Site == "Big Cypress"]*coeff_bigger,  pch=1, bg = "white")
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(3,8, 13,18, 23, 28, 33, 38, 44)*coeff_bigger,c("3", "8", "13", "23", "28", "33", "38", "44"), xpd=T, horiz = T, title="Height at t +1 (cm)")


plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topright", 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

#Cape Canaveral IMAGE PLOT
#Future Diameter
image.plot(x1seq_seedlings, x2seq_seedlings, z1_CC, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, legend.lab = "Diameter at t+1 (mm)", main = "Cape Canaveral")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Cape Canaveral"], seedlings_g$Height_t[seedlings_g$Site == "Cape Canaveral"], cex=seedlings_g$Diameter_tplus1[seedlings_g$Site == "Cape Canaveral"]*coeff_bigger,  pch=1)
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")

#Future Height
image.plot(x1seq_seedlings, x2seq_seedlings, z2_CC, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, main ="Cape Canaveral",legend.lab = "Height at t+1 (cm)")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Cape Canaveral"], seedlings_g$Height_t[seedlings_g$Site == "Cape Canaveral"], cex=seedlings_g$Height_tplus1[seedlings_g$Site == "Cape Canaveral"]*coeff_bigger,  pch=1, bg = "white")
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(3,8, 13,18, 23, 28, 33, 38, 44)*coeff_bigger,c("3", "8", "13", "23", "28", "33", "38", "44"), xpd=T, horiz = T, title="Height at t +1 (cm)")

#Chekika IMAGE PLOT
#Future Diameter
image.plot(x1seq_seedlings, x2seq_seedlings, z1_C, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, legend.lab = "Diameter at t+1 (mm)", main = "Chekika")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Chekika"], seedlings_g$Height_t[seedlings_g$Site == "Chekika"], cex=seedlings_g$Diameter_tplus1[seedlings_g$Site == "Chekika"]*coeff_bigger,  pch=1)
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")

#Future Height
image.plot(x1seq_seedlings, x2seq_seedlings, z2_C, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, main ="Chekika",legend.lab = "Height at t+1 (cm)")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Chekika"], seedlings_g$Height_t[seedlings_g$Site == "Chekika"], cex=seedlings_g$Height_tplus1[seedlings_g$Site == "Chekika"]*coeff_bigger,  pch=1, bg = "white")
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(3,8, 13,18, 23, 28, 33, 38, 44)*coeff_bigger,c("3", "8", "13", "23", "28", "33", "38", "44"), xpd=T, horiz = T, title="Height at t +1 (cm)")

#Fort Pierce IMAGE PLOT
#Future Diameter
image.plot(x1seq_seedlings, x2seq_seedlings, z1_FP, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, legend.lab = "Diameter at t+1 (mm)", main = "Fort Pierce")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Fort Pierce"], seedlings_g$Height_t[seedlings_g$Site == "Fort Pierce"], cex=seedlings_g$Diameter_tplus1[seedlings_g$Site == "Fort Pierce"]*coeff_bigger,  pch=1)
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")

#Future Height
image.plot(x1seq_seedlings, x2seq_seedlings, z2_FP, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, main ="Fort Pierce",legend.lab = "Height at t+1 (cm)")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Fort Pierce"], seedlings_g$Height_t[seedlings_g$Site == "Fort Pierce"], cex=seedlings_g$Height_tplus1[seedlings_g$Site == "Fort Pierce"]*coeff_bigger,  pch=1, bg = "white")
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(3,8, 13,18, 23, 28, 33, 38, 44)*coeff_bigger,c("3", "8", "13", "23", "28", "33", "38", "44"), xpd=T, horiz = T, title="Height at t +1 (cm)")

#Punta Gorda IMAGE PLOT
#Future Diameter
image.plot(x1seq_seedlings, x2seq_seedlings, z1_PG, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, legend.lab = "Diameter at t+1 (mm)", main = "Punta Gorda")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Punta Gorda"], seedlings_g$Height_t[seedlings_g$Site == "Punta Gorda"], cex=seedlings_g$Diameter_tplus1[seedlings_g$Site == "Punta Gorda"]*coeff_bigger,  pch=1)
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")

#Future Height
image.plot(x1seq_seedlings, x2seq_seedlings, z2_PG, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, main ="Punta Gorda",legend.lab = "Height at t+1 (cm)")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Punta Gorda"], seedlings_g$Height_t[seedlings_g$Site == "Punta Gorda"], cex=seedlings_g$Height_tplus1[seedlings_g$Site == "Punta Gorda"]*coeff_bigger,  pch=1, bg = "white")
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(3,8, 13,18, 23, 28, 33, 38, 44)*coeff_bigger,c("3", "8", "13", "23", "28", "33", "38", "44"), xpd=T, horiz = T, title="Height at t +1 (cm)")

#Wild Turkey IMAGE PLOT
#Future Diameter
image.plot(x1seq_seedlings, x2seq_seedlings, z1_WT, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, legend.lab = "Diameter at t+1 (mm)", main = "Wild Turkey")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Wild Turkey"], seedlings_g$Height_t[seedlings_g$Site == "Wild Turkey"], cex=seedlings_g$Diameter_tplus1[seedlings_g$Site == "Wild Turkey"]*coeff_bigger,  pch=1)
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(0.5,1, 1.5, 2,2.5, 3, 3.5,4), c("0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4"), xpd=T, horiz = T, title="Diameter at t +1 (mm)")

#Future Height
image.plot(x1seq_seedlings, x2seq_seedlings, z2_WT, xlab = "Diameter (mm) ", ylab="Height (cm)", col=my.palette, main ="Wild Turkey",legend.lab = "Height at t+1 (cm)")
points(seedlings_g$Diameter_t[seedlings_g$Site == "Wild Turkey"], seedlings_g$Height_t[seedlings_g$Site == "Wild Turkey"], cex=seedlings_g$Height_tplus1[seedlings_g$Site == "Wild Turkey"]*coeff_bigger,  pch=1, bg = "white")
legend(0, 18.5, pch=1,text.width=0.1, pt.cex = c(3,8, 13,18, 23, 28, 33, 38, 44)*coeff_bigger,c("3", "8", "13", "23", "28", "33", "38", "44"), xpd=T, horiz = T, title="Height at t +1 (cm)")

dev.off()

### FIG 3: D2 Survival ##### 
df2 <- expand.grid(x1seq_larges, x2seq_larges)
#BC, CC, and PG
b0 <- p.vec_BC[19]
b1 <- p.vec_BC[20]
b2 <- p.vec_BC[21]
df2$BC_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

#Chekika
b0 <- p.vec_C[19]
b1 <- p.vec_C[20]
b2 <- p.vec_C[21]
df2$C_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

#Fort Pierce
b0 <- p.vec_FP[19]
b1 <- p.vec_FP[20]
b2 <- p.vec_FP[21]
df2$FP_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

#Wild Turkey
b0 <- p.vec_WT[19]
b1 <- p.vec_WT[20]
b2 <- p.vec_WT[21]
df2$WT_s2 <- invlogit(b0 + b1*df2$Var1 + b2*df2$Var2)

#Figure 3 in ggplot
#Big Cypress, Cape Canaveral and Punta Gorda
plot_BCs2<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= BC_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  geom_point(data=subset(larges, !is.na(larges$Surv_tplus1) &
                           (larges$Site == "Big Cypress" |
                              larges$Site == "Cape Canaveral" |
                              larges$Site == "Punta Gorda")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
             shape=1, size=0.8) +
  scale_color_manual(values=c("red", "black")) +
  ggtitle("(a) BC, CC, PG") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Chekika
plot_Cs2<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= C_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(b) C") + 
  geom_point(data=subset(larges, !is.na(larges$Surv_tplus1) &
                           (larges$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
             shape=1, size=0.8) +
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Fort Pierce
plot_FPs2<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= FP_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(c) FP") + 
  geom_point(data=subset(larges, !is.na(larges$Surv_tplus1) &
                           (larges$Site == "Fort Pierce")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)),
             shape = 1, size=0.8) +
  scale_color_manual(values = c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

#Wild Turkey
plot_WTs2<- ggplot() + xlim(1.6, 800) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= WT_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(d) WT") + 
  geom_point(data=subset(larges, !is.na(larges$Surv_tplus1) &
                           (larges$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)),
             shape = 1, size=0.8) +
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

g <- ggplotGrob(plot_BCs2 + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot_BCs2 <- plot_BCs2 + theme(legend.position ="none")
plot_Cs2 <- plot_Cs2 + theme(legend.position ="none")
plot_FPs2  <- plot_FPs2  + theme(legend.position = "none")
plot_WTs2 <- plot_WTs2 + theme(legend.position = "none")



postscript("./Figures/D2_survival.eps", width=width.cm/2.54, 
           height=(2*height.cm)/2.54, pointsize=pointsize, 
           encoding = "TeXtext.enc")
lay <- rbind(c(1, 1, 2, 2, 5, 5),
             c(1, 1, 2, 2, 5, 5),
             c(3, 3, 4, 4, 5, 5),
             c(3, 3, 4, 4, 5, 5))
grid.arrange(plot_BCs2, plot_Cs2,plot_FPs2, plot_WTs2, legend,
             layout_matrix=lay, widths=c(1, 1, 1, 1, 1, .7))
dev.off()

#*--- Try with separate panels for data points #####
#Big Cypress, Cape Canaveral and Punta Gorda
plot_BCs2_pred<- ggplot() + xlim(1.6, 400) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= BC_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(a) BC, CC, PG") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
 guides(fill = guide_colorbar(title.position = "right", order =2))

plot_BCs2_meas <- ggplot() + xlim(1.6, 400) + ylim(16, 800) + 
geom_point(data=subset(larges,
                !is.na(larges$Surv_tplus1) & (larges$Site == "Big Cypress" |
                            larges$Site == "Cape Canaveral" |
                            larges$Site == "Punta Gorda")),
           aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
           shape=1, size=0.8) +
  ggtitle(" ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Chekika
plot_Cs2_pred<- ggplot() + xlim(1.6, 400) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= C_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(b) C") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot_Cs2_meas <- ggplot()  + xlim(1.6, 400) + ylim(16, 800) + 
  geom_point(data=subset(larges,
              !is.na(larges$Surv_tplus1) & (larges$Site == "Chekika")),
              aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
               shape=1, size=0.8) +
  ggtitle(" ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Fort Pierce
plot_FPs2_pred<- ggplot() + xlim(1.6, 400) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= FP_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(c) FP") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))


plot_FPs2_meas <- ggplot()  + xlim(1.6, 400) + ylim(16, 800) + 
  geom_point(data=subset(larges,
                  !is.na(larges$Surv_tplus1) &
                    (larges$Site == "Fort Pierce")),
                      aes(x=Diameter_t, y=Height_t, 
                      col = factor(Surv_tplus1)), 
                      shape=1, size=0.8) +
  ggtitle(" ") + 
  scale_color_manual(values=c("red", "black")) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(col = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(col = guide_legend(title.position = "right", order =1)) 

#Wild Turkey
plot_WTs2_pred<- ggplot() + xlim(1.6, 400) + ylim(16, 800) + 
  geom_tile(data = df2, aes(x=Var1, y = Var2, fill= WT_s2)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0.6,1)) + 
  ggtitle("(d) WT") + 
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(fill = guide_colorbar(title.position = "right", order =2))

plot_WTs2_meas <- ggplot()  + xlim(1.6, 400) + ylim(16, 800) + 
  geom_point(data=subset(larges,
        !is.na(larges$Surv_tplus1) & (larges$Site == "Wild Turkey")),
        aes(x=Diameter_t, y=Height_t, col = factor(Surv_tplus1)), 
        shape=1, size=0.8) +
  ggtitle(" ") + 
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

###FIG 4: GROWTH OF D2 INDIVIDUALS #####

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


#x11(width = width.cm/2.54, height = height.cm/2.54, 
#pointsize = pointsize)


#Big Cypress
plot1_BC<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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

plot2_BC<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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
plot1_CC<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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

plot2_CC<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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
plot1_C<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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

plot2_C<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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
plot1_FP<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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

plot2_FP<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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
plot1_PG<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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

plot2_PG<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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
plot1_WT<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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

plot2_WT<- ggplot() + xlim(1.6, 500) + ylim(16, 800) + 
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
plot(x2seq_larges, pf_BC, xlab="Height at time t (cm)", ylim=c(0,1), xlim=c(0, 500),
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

plot(x1seq_larges, y_BC, xlab="Diameter at time t (mm)", xlim=c(0,500),
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

plot(x1seq_larges, y_BC, xlab="Diameter at time t (mm)", xlim=c(0,500),
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
barplot(c(lam.stable_BC, lam.stable_CC, lam.stable_C, lam.stable_FP,
          lam.stable_PG, lam.stable_WT), 
        col=cols,
        ylim=c(0.98, 1),
        ylab=expression(paste(lambda)), 
        names.arg=c("BC", "CC", "C", "FP",
                    "PG", "WT"),
         xpd=F)
dev.off()

load("./BC/BC_pool.RData")
load("./CC/CC_pool.RData")
load("./C/C_pool.RData")
load("./FP/FP_pool.RData")
load("./PG/PG_pool.RData")
load("./WT/WT_pool.RData")

#Figure 9: Barplot of contribution of kernel components to differences in lambda #####
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
    mfrow=c(3, 2))

barplot(BC_pool, col=cols[1], ylim=c(-0.13, 0.16),
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), axes=T, cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["BC"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(CC_pool, col=cols[2], ylim=c(-0.13, 0.16), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=0.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["CC"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(C_pool, col=cols[3], ylim=c(-0.13, 0.16), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["C"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(FP_pool, col=cols[4], ylim=c(-0.13, 0.16), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["FP"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(PG_pool, col=cols[5], ylim=c(-0.13, 0.16), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["PG"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

barplot(WT_pool, col=cols[6], ylim=c(-0.13, 0.16), 
        names.arg=c("Seedling", "Maturation", "Fertility", "Adult"), cex.names=.6, cex.axis=1)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["WT"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(3, 1.5), cex=1)
abline(h=0)

dev.off()

#Figure A1: SSD and RV ######

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

load("./BC/y_diam.RData")
load("./BC/y_height.RData")
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

plot(y_diam, ss_diam_BC,xlab="Diameter (mm)",ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=5, lwd=1, ylim=c(0, 0.24), col=cols[1]); 
lines(y_diam, ss_diam_CC, lwd=1, lty=3, col=cols[2])
lines(y_diam, ss_diam_C, lwd=1, lty=1, col=cols[3])
lines(y_diam, ss_diam_FP, lwd=1, lty=4, col=cols[4])
lines(y_diam, ss_diam_PG, lwd=1, lty=6, col=cols[5])
lines(y_diam, ss_diam_WT, lwd=1, lty=7, col=cols[6])

title(main="(a) Stable diameter distribution", cex.main=1);


plot(y_height,ss_height_BC,xlab="Height (cm)",ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=5, lwd=1, ylim=c(0, 0.24), col=cols[1]); 
lines(y_height, ss_height_CC, lwd=1, lty=3, col=cols[2])
lines(y_height, ss_height_C, lwd=1, lty=1, col= cols[3])
lines(y_height, ss_height_FP, lwd=1, lty=4, col=cols[4])
lines(y_height, ss_height_PG, lwd=1, lty=6, col=cols[5])
lines(y_height, ss_height_WT, lwd=1, lty=7, col=cols[6])
title(main="(b) Stable height distribution", cex.main=1); 
legend(350, 0.15, lty=c(5, 1, 3, 4, 6, 7), lwd=1, col=cols,
       c("BC", "CC", "C", "FP", "PG", "WT"), cex=1,  seg.len=3)


plot(y_diam, rv_diam_BC, xlab="Diameter(mm)", ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=5, lwd=1, ylim=c(0, 0.022), col=cols[1])
lines(y_diam, rv_diam_CC, lwd=1, lty=3, col=cols[2])
lines(y_diam, rv_diam_C, lwd=1, lty=1, col=cols[3])
lines(y_diam, rv_diam_FP, lwd=1, lty=4, col=cols[4])
lines(y_diam, rv_diam_PG, lwd=1, lty=6, col=cols[5])
lines(y_diam, rv_diam_WT, lwd=1, lty=7, col=cols[6] )
title(main="(c) Reproductive value: by diameter", cex.main=1);

 
plot(y_height, rv_height_BC, xlab="Height (cm)", ylab="frequency",
     type="l", cex=1, cex.axis=1, cex.lab=1, lty=5, lwd=1, ylim=c(0, 0.022), col=cols[1])
lines(y_height, rv_height_CC, lwd=1, lty=1, col=cols[2])
lines(y_height, rv_height_C, lwd=1, lty=3, col=cols[3])
lines(y_height, rv_height_FP, lwd=1, lty=4, col=cols[4])
lines(y_height, rv_height_PG, lwd=1, lty=6, col=cols[5])
lines(y_height, rv_height_WT, lwd=1, lty=7, col=cols[6])
title(main="(d) Reproductive value: by height", cex.main=1); 


dev.off()




#Figure A2: Marginal elasticity by diameter ######
#setwd("/Users/curculion/Dropbox/matrix/outputs/sites/Eastern")
# load("total.elas_D1_E.RData")
# load("total.elas_F_E.RData")
# load("total.elas_G_E.RData")
# load("total.elas_D2_E.RData")

load("/Users/curculion/Documents/GitHub/Schinus_IPM_code/BC/total.elas.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code//CC/total.elas.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code//C/total.elas.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code//FP/total.elas.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code//PG/total.elas.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code//WT/total.elas.RData")
load("/Users/curculion/Documents/GitHub/Schinus_IPM_code//Overall/total.elas.RData")

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

#
plot(x1_seedlings, rowSums(total.elas_D1_BC) , xlim=c(0, 1.6), xlab="", 
     ylab="", main="",  type='l', lwd=1, lty=5, cex.axis=1.5, col=cols[1])
text(0.2, 0.029, "Seedling", cex=1.5)
lines(x1_seedlings, rowSums(total.elas_D1_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x1_seedlings, rowSums(total.elas_D1_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x1_seedlings, rowSums(total.elas_D1_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x1_seedlings, rowSums(total.elas_D1_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x1_seedlings, rowSums(total.elas_D1_WT), type='l', lwd=1, lty=7, col=cols[6])
mtext(side=3, "(a) Marginal elasticity over diameter", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.elas_F_BC), xlim=c(1.6, 800), xlab="", 
     ylab= "", main="", type='l', lwd=1, lty=5, cex.axis=1.5, col=cols[1])
text(65, 0.029, "Fertility", cex=1.5)
lines(x1_adults, rowSums(total.elas_F_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x1_adults, rowSums(total.elas_F_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x1_adults, rowSums(total.elas_F_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x1_adults, rowSums(total.elas_F_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x1_adults, rowSums(total.elas_F_WT), type='l', lwd=1, lty=7, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_seedlings, rowSums(total.elas_G_BC), xlim=c(0, 1.6), xlab="",
     ylab="", main="", type='l', lwd=1, lty=5, cex.axis=1.5, col=cols[1])
text(0.2, 0.039, "Maturation", cex=1.5)
lines(x1_seedlings, rowSums(total.elas_G_CC), type='l', lty=1, col=cols[2])
lines(x1_seedlings, rowSums(total.elas_G_C), type='l', lty=3, col=cols[3])
lines(x1_seedlings, rowSums(total.elas_G_FP), type='l', lty=4, col=cols[4])
lines(x1_seedlings, rowSums(total.elas_G_PG), type='l', lty=6, col=cols[5])
lines(x1_seedlings, rowSums(total.elas_G_WT), type='l', lty=7, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.elas_D2_BC), xlim=c(1.6, 800), xlab=" ",
     ylab="", main="",  type='l',  lwd=1, lty=5, cex.axis=1.5, col=cols[1])
text(65, 0.029, "Adult", cex=1.5)
lines(x1_adults, rowSums(total.elas_D2_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x1_adults, rowSums(total.elas_D2_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x1_adults, rowSums(total.elas_D2_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x1_adults, rowSums(total.elas_D2_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x1_adults, rowSums(total.elas_D2_WT), type='l', lwd=1, lty=7, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)




plot(x2_seedlings, colSums(total.elas_D1_BC), ylim=c(0, 0.03), xlim=c(0, 16), xlab=" ", 
     ylab="", main="",  type='l', lwd=1, lty=5, cex.axis=1.5, col=cols[1])
text(2, 0.029, "Seedling", cex=1.5)
lines(x2_seedlings, colSums(total.elas_D1_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x2_seedlings, colSums(total.elas_D1_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_seedlings, colSums(total.elas_D1_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_seedlings, colSums(total.elas_D1_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x2_seedlings, colSums(total.elas_D1_WT), type='l', lwd=1, lty=7, col=cols[6])
mtext(side=3, "(b) Marginal elasticity over height", line=1, cex=1.2, adj=0)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.elas_F_BC), ylim=c(0, 0.03), xlim=c(16,800),  xlab=" ",
     ylab= "", main="", type='l', lwd=1, lty=5, cex.axis=1.5, col=cols[1])
text(65, 0.029, "Fertility", cex=1.5)
lines(x2_adults, colSums(total.elas_F_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x2_adults, colSums(total.elas_F_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_adults, colSums(total.elas_F_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_adults, colSums(total.elas_F_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x2_adults, colSums(total.elas_F_WT), type='l', lwd=1, lty=7, col=cols[6])
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)


plot(x2_seedlings, colSums(total.elas_G_BC), ylim= c(0, 0.03), xlim=c(0, 16), xlab="",
     ylab="", main="", type='l', lwd=1, lty=5, cex.axis=1.5, col=cols[1])
lines(x2_seedlings, colSums(total.elas_G_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x2_seedlings, colSums(total.elas_G_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_seedlings, colSums(total.elas_G_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_seedlings, colSums(total.elas_G_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x2_seedlings, colSums(total.elas_G_WT), type='l', lwd=1, lty=7, col=cols[6])
text(2, 0.029, "Maturation", cex=1.5)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.elas_D2_BC), ylim= c(0, 0.03), xlim=c(16, 800),  xlab=" ", 
     ylab="", main="",  type='l',  lwd=1, lty=5, cex.axis=1.5, col=cols[1])
lines(x2_adults, colSums(total.elas_D2_CC), type='l', lwd=1, lty=1, col=cols[2])
lines(x2_adults, colSums(total.elas_D2_C), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_adults, colSums(total.elas_D2_FP), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_adults, colSums(total.elas_D2_PG), type='l', lwd=1, lty=6, col=cols[5])
lines(x2_adults, colSums(total.elas_D2_WT), type='l', lwd=1, lty=7, col=cols[6])
text(75, 0.029, "Adult", cex=1.5)
legend(250, 0.03, lty=c(3, 1, 5, 4, 6, 7), col=cols, lwd=1,
       c("BC", "CC", "C", "FP", "PG", "WT"),
       seg.len=2)
mtext(side=2, "Elasticity", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)
dev.off() #x11(width = width.cm/2.54, height = 1.5*width.cm/(2.54), 
  #   pointsize = pointsize)




###FIGURE B1: Contribution of matrix elements to differences in lambda by diameter #####
load("./BC/total.contrib_BC_pool.RData")
load("./CC/total.contrib_CC_pool.RData")
load("./C/total.contrib_C_pool.RData")
load("./FP/total.contrib_FP_pool.RData")
load("./PG/total.contrib_PG_pool.RData")
load("./WT/total.contrib_WT_pool.RData")

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


plot(x1_seedlings, rowSums(total.contrib_D1_BC_pool), ylim=c(-0.003, 0.01),
     xlim=c(0, 1.6), xlab="", ylab="", main="",  type='l', lwd=1, 
     lty=1, cex.axis=1.5, col=cols[1])
text(0.2, 0.01, "Seedling", cex=1.5)
lines(x1_seedlings, rowSums(total.contrib_D1_CC_pool), type='l', lwd=1, lty=2, col=cols[2])
lines(x1_seedlings, rowSums(total.contrib_D1_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x1_seedlings, rowSums(total.contrib_D1_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x1_seedlings, rowSums(total.contrib_D1_PG_pool), type='l', lwd=1, lty=5, col=cols[5])
lines(x1_seedlings, rowSums(total.contrib_D1_WT_pool), type='l', lwd=1, lty=6, col=cols[6])
mtext(side=3, "(a) Contributions by diameter", line=1, cex=1.2, adj=0)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.contrib_F_BC_pool), ylim=c(-0.003, 0.01), xlim=c(1.6, 800), xlab="", 
     ylab= "", main="", type='l', lwd=1, lty=1, cex.axis=1.5, col=cols[1])
text(65, 0.01, "Fertility", cex=1.5)
lines(x1_adults, rowSums(total.contrib_F_CC_pool), type='l', lwd=1, lty=2, col=cols[2])
lines(x1_adults, rowSums(total.contrib_F_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x1_adults, rowSums(total.contrib_F_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x1_adults, rowSums(total.contrib_F_PG_pool), type='l', lwd=1, lty=5, col=cols[5])
lines(x1_adults, rowSums(total.contrib_F_WT_pool), type='l', lwd=1, lty=6, col=cols[6])
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x1_seedlings, rowSums(total.contrib_G_BC_pool), ylim=c(-0.003, 0.01), xlim=c(0, 1.6), xlab="",
     ylab="", main="", type='l', lwd=1, lty=1, cex.axis=1.5, col=cols[1])
text(0.2, 0.01, "Maturation", cex=1.5)
lines(x1_seedlings, rowSums(total.contrib_G_CC_pool), type='l', lty=2, col=cols[2])
lines(x1_seedlings, rowSums(total.contrib_G_C_pool), type='l', lty=3, col=cols[3])
lines(x1_seedlings, rowSums(total.contrib_G_FP_pool), type='l', lty=4, col=cols[4])
lines(x1_seedlings, rowSums(total.contrib_G_PG_pool), type='l', lty=5, col=cols[5])
lines(x1_seedlings, rowSums(total.contrib_G_WT_pool), type='l', lty=6, col=cols[6])
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)


plot(x1_adults, rowSums(total.contrib_D2_BC_pool), ylim=c(-0.003, 0.01), xlim=c(1.6, 800), xlab=" ",
     ylab="", main="",  type='l',  lwd=1, lty=1, cex.axis=1.5, col=cols[1])
text(65, 0.01, "Adult", cex=1.5)
lines(x1_adults, rowSums(total.contrib_D2_CC_pool), type='l', lwd=1, lty=2, col=cols[2])
lines(x1_adults, rowSums(total.contrib_D2_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x1_adults, rowSums(total.contrib_D2_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x1_adults, rowSums(total.contrib_D2_PG_pool), type='l', lwd=1, lty=5, col=cols[5])
lines(x1_adults, rowSums(total.contrib_D2_WT_pool), type='l', lwd=1, lty=6, col=cols[6])
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Diameter (mm)", line=2,cex=1)

plot(x2_seedlings, colSums(total.contrib_D1_BC_pool), ylim=c(-0.003, 0.01), xlim=c(0, 16), xlab=" ", 
     ylab="", main="",  type='l', lwd=1, lty=1, cex.axis=1.5, col=cols[1])
text(2, 0.01, "Seedling", cex=1.5)
lines(x2_seedlings, colSums(total.contrib_D1_CC_pool), type='l', lwd=1, lty=1, col=cols[2])
lines(x2_seedlings, colSums(total.contrib_D1_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_seedlings, colSums(total.contrib_D1_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_seedlings, colSums(total.contrib_D1_PG_pool), type='l', lwd=1, lty=6, col=cols[5])
lines(x2_seedlings, colSums(total.contrib_D1_WT_pool), type='l', lwd=1, lty=7, col=cols[6])
mtext(side=3, "(b) Contributions by height", line=1, cex=1.2, adj=0)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.contrib_F_BC_pool), ylim=c(-0.003, 0.01), xlim=c(16,800),  xlab=" ",
     ylab= "", main="", type='l', lwd=1, lty=1, cex.axis=1.5, col=cols[1])
text(65, 0.01, "Fertility", cex=1.5)
lines(x2_adults, colSums(total.contrib_F_CC_pool), type='l', lwd=1, lty=2, col=cols[2])
lines(x2_adults, colSums(total.contrib_F_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_adults, colSums(total.contrib_F_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_adults, colSums(total.contrib_F_PG_pool), type='l', lwd=1, lty=5, col=cols[5])
lines(x2_adults, colSums(total.contrib_F_WT_pool), type='l', lwd=1, lty=6, col=cols[6])
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)


plot(x2_seedlings, colSums(total.contrib_G_BC_pool), ylim=c(-0.003, 0.01), xlim=c(0, 16), xlab="",
     ylab="", main="", type='l', lwd=1, lty=1, cex.axis=1.5, col=cols[1])
lines(x2_seedlings, colSums(total.contrib_G_CC_pool), type='l', lwd=1, lty=2, col=cols[2])
lines(x2_seedlings, colSums(total.contrib_G_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_seedlings, colSums(total.contrib_G_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_seedlings, colSums(total.contrib_G_PG_pool), type='l', lwd=1, lty=5, col=cols[5])
lines(x2_seedlings, colSums(total.contrib_G_WT_pool), type='l', lwd=1, lty=6, col=cols[6])
text(2, 0.01, "Maturation", cex=1.5)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)

plot(x2_adults, colSums(total.contrib_D2_BC_pool), ylim=c(-0.003, 0.01), xlim=c(16, 800),  xlab=" ", 
     ylab="", main="",  type='l',  lwd=1, lty=1, cex.axis=1.5, col=cols[1])
lines(x2_adults, colSums(total.contrib_D2_CC_pool), type='l', lwd=1, lty=2, col=cols[2])
lines(x2_adults, colSums(total.contrib_D2_C_pool), type='l', lwd=1, lty=3, col=cols[3])
lines(x2_adults, colSums(total.contrib_D2_FP_pool), type='l', lwd=1, lty=4, col=cols[4])
lines(x2_adults, colSums(total.contrib_D2_PG_pool), type='l', lwd=1, lty=5, col=cols[5])
lines(x2_adults, colSums(total.contrib_D2_WT_pool), type='l', lwd=1, lty=6, col=cols[6])
text(65, 0.01, "Adult", cex=1.5)
legend(250, 0.01, lty=c(1:6), lwd=1, c("BC", "CC", "C", "FP", "PG", "WT"),
       seg.len=2, col=cols)
mtext(side=2, "Contribution", line=2, cex=1)
mtext(side=1, "Height (cm)", line=2,cex=1)
dev.off()



#Figure B3: Contributions of matrix elements to variation in lambda by diameter #####
#With six sites calculating the variability takes too long (more than 50 hours at least)
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


##### Look at sensitivity 
sens_overall <- readRDS("/Users/curculion/Documents/GitHub/Schinus_IPM_code/Overall/sens_overall.rds")

sens_D1_overall <-readRDS("/Users/curculion/Documents/GitHub/Schinus_IPM_code/Overall/sens_D1_overall.rds") 
sens_D2_overall <- readRDS("/Users/curculion/Documents/GitHub/Schinus_IPM_code/Overall/sens_D2_overall.rds") 