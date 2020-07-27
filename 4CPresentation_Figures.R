#4CPresentation_Figures
library(sp)
library(rgdal)
library(ggspatial)
library(rgeos)
library(tmap)
library(sf)
library(raster)
library(dplyr)
library(mapview)
library(ggplot2)
library(shiny)
library(maptools)
library(plyr)
library(ggsn)
library(arm)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(cowplot)
library(gtable)
library(viridis)
cols=viridis(6)

my.palette <- brewer.pal(11, "RdYlBu")
setwd("/Users/curculion/Documents/GitHub/Schinus_IPM_code")
#FIG 1: Map of Study Sites #####
sites <- read.csv('sites.csv', as.is=TRUE)


sitesSpatial <- SpatialPointsDataFrame(
	coords= cbind(sites$Longitude, sites$Latitude),
	data=sites,
	proj4string=CRS('+proj=longlat +datum=WGS84 +no defs +ellps=WGS84 +towgs84=0,0,0'))

load("florida.RData")


xs<-c(-87.63723, -78.5)
ys<-c(31.00211, 24.52042)



png(file="map.png", width=13.5, height=10, units='in', res=50)
par(ps=24, oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
plot(xs, ys, pch=NA)

sp::plot(florida, col='gray80',add=TRUE, xpd=NA, usePolypath=F)
points(sitesSpatial, pch=c(15, 15, 16, 16, 17, 17), cex=4, bg=c('purple', 'purple', 'blue', 'blue', 'red', 'red'))

#points(-87.63723, 31.00211, pch=15, col="red", cex=4)
#points(-78.5, 24.52042, pch=15, col="red", cex=4)
legend(-87, 28, pch=c(15, 16, 17), pt.cex=4, y.intersp=3, x.intersp=2.5, legend=c("Hybrid", "Eastern", "Western"), bty='n')
#text(sitesSpatial, labels=sites$Site)
rect(-87.2, 26.5, -85.5, 28)

#xleft, ybottom, xright, ytop


text(-82, 26, "Big Cypress")
text(-82.8, 26.3, "Wild Turkey")
text(-83.2, 26.9, "Punta Gorda")
text(-79.6, 25.6, "Chekika")
text(-79.5, 27.4, "Fort Pierce")
text(-79.5, 28.4, "Cape Canaveral")
dev.off()

set_Polypath(F)
florida2 <- fortify(florida, region=florida@data$NAME_2)
ggplot() +
  layer_spatial(florida, usePolypath=F)

spplot(florida@data$ID_0, sp.layout=list(scale))

florida@data$id <- rownames(florida@data)
florida2 <- fortify(florida, region="id")
florida2$id <- florida@data$ID_2
florida3 <- join(florida2, florida@data, by="id")
florida2@data$ID_2 <- as.numeric(florida2@data$ID_2)
florida2 <- 
ggplot(data=florida) + 
  geom_sf() + 
  annotation_scale(location = "bl", width_hint=0.5)


ggplot(florida3) + aes(long, lat, group=group) +
  geom_polygon() + geom_path(color="white") + coord_equal() +
  raster::scalebar(lon=-85, lat=26, distance_lon=100, distance_lat = 20,
           distance_legend=40, dist_unit="km", 
           arrow_length=100, arrow_distance = 60, arrow_north_size=6)

#FIG 2: tikz picture #####
#FIG 3: Seedling Survival #####

x1seq_seedlings<-seq(0, 1.6, length.out=50)

x2seq_seedlings<-seq(0, 16, length.out=50)
x1seq_larges <- seq(1.6, 800, length.out = 50 )
x2seq_larges <- seq (16, 800, length.out = 50)

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

png(file="./PowerPoint/seedling_survival_BC.png", width=2.91, 
    height=2.91, units='in', res=300)
plotBC <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= BC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=15))+
  theme(axis.title=element_text(size=15)) +
  theme(title = element_text(size=15)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  ggtitle("Punta Gorda, Wild \n Turkey & Big Cypress") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Big Cypress" |
                              seedlings$Site == "Punta Gorda" |
                              seedlings$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t, size=0.5,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2)) +
theme(plot.margin=unit(c(0.2,1.4,0,0.2),"cm")) + #top, right, bottom, left
  guides(size=F) +
theme(legend.position = "none")
dev.off()

#Chekika
b0 <- p.vec_C[1]
b1 <- p.vec_C[2]
b2 <- p.vec_C[3]
df$C_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_Cs <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= C_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=24))+
  theme(axis.title=element_text(size=36)) +
  theme(title = element_text(size=24)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(legend.text = element_text(size=24)) + 
  ggtitle("Chekika") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t, size=3,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2)) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  guides(size=F)

#Fort Pierce
b0 <- p.vec_FP[1]
b1 <- p.vec_FP[2]
b2 <- p.vec_FP[3]
df$FP_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_FPs <- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= FP_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=24))+
  theme(axis.title=element_text(size=36)) +
  theme(title = element_text(size=24)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(legend.text = element_text(size=24)) + 
  ggtitle("Chekika") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Fort Pierce")),
             aes(x=Diameter_t, y=Height_t, size=3,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2)) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  guides(size=F)

#Cape Canaveral
b0 <- p.vec_CC[1]
b1 <- p.vec_CC[2]
b2 <- p.vec_CC[3]
df$CC_s <- invlogit(b0 + b1*df$Var1 + b2*df$Var2)

plot_CCs<-ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= CC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=24))+
  theme(axis.title=element_text(size=36)) +
  theme(title = element_text(size=24)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(legend.text = element_text(size=24)) + 
  ggtitle("Chekika") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Cape Canaveral")),
             aes(x=Diameter_t, y=Height_t, size=3,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("Diameter at time t (mm)") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  theme(legend.box.margin = margin (20, -5, 20, -11)) +
  guides(shape = guide_legend(title.position = "right", order =1)) +
  guides(fill = guide_colorbar(title.position = "right", order =2)) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  guides(size=F)


g <- ggplotGrob(plot_BCs + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

plot_BCs2 <- plot_BCs + theme(legend.position ="none")
plot_CCs2 <- plot_CCs + theme(legend.position ="none")
plot_Cs2  <- plot_Cs  + theme(legend.position = "none")
plot_FPs2 <- plot_FPs + theme(legend.position = "none")

png(file="./PowerPoint/seedling_survival.png", width=11.64, 
    height=2.91, units='in', res=300)
lay <- c(1,2, 3, 4)
#lay <- rbind(c(1, 1, 2, 2, 5, 5),
             #c(1, 1, 2, 2, 5, 5),
             #c(3, 3, 4, 4, 5, 5),
             #c(3, 3, 4, 4, 5, 5))
grid.arrange(plotBC, plotBC, plotBC, plotBC, legend, ncol=5)
dev.off()






BC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= BC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=15))+
  theme(axis.title=element_text(size=15)) +
  theme(title = element_text(size=15)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(axis.text = element_text(size=15)) +
  ggtitle("Punta Gorda, \n Wild Turkey \n & Big Cypress") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Big Cypress" |
                              seedlings$Site == "Punta Gorda" |
                              seedlings$Site == "Wild Turkey")),
             aes(x=Diameter_t, y=Height_t,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  guides(shape = guide_legend(title.position = "right", order =1, text=15, 
                              label.theme=element_text(size=15), override.aes = list(size=5))) +
  guides(fill = guide_colorbar(title.position = "right", order =2, 
                               label.theme=element_text(size=15), text=15)) +
  guides(size=F)

C<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= C_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=15))+
  theme(axis.title=element_text(size=15)) +
  theme(title = element_text(size=15)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(axis.text = element_text(size=15)) +
  ggtitle("\n \n Chekika") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Chekika")),
             aes(x=Diameter_t, y=Height_t,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  guides(shape = guide_legend(title.position = "right", order =1, text=15, 
                              label.theme=element_text(size=15), override.aes = list(size=5))) +
  guides(fill = guide_colorbar(title.position = "right", order =2, 
                               label.theme=element_text(size=15), text=15)) +
  guides(size=F)

FP<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= FP_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=15))+
  theme(axis.title=element_text(size=15)) +
  theme(title = element_text(size=15)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(axis.text = element_text(size=15)) +
  ggtitle("\n \n Fort Pierce") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Fort Pierce")),
             aes(x=Diameter_t, y=Height_t,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  guides(shape = guide_legend(title.position = "right", order =1, text=15, 
                              label.theme=element_text(size=15), override.aes = list(size=5))) +
  guides(fill = guide_colorbar(title.position = "right", order =2, 
                               label.theme=element_text(size=15), text=15)) +
  guides(size=F)

CC<- ggplot() + xlim(0, 1.6) + ylim(0, 16) + 
  geom_tile(data = df, aes(x=Var1, y = Var2, fill= CC_s)) + 
  scale_fill_gradientn(colors=my.palette, limits=c(0,1)) + 
  theme(axis.text=element_text(size=15))+
  theme(axis.title=element_text(size=15)) +
  theme(title = element_text(size=15)) +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(axis.text = element_text(size=15)) +
  ggtitle("\n \n Cape Canaveral") + 
  geom_point(data=subset(seedlings, !is.na(seedlings$Surv_tplus1) &
                           (seedlings$Site == "Cape Canaveral")),
             aes(x=Diameter_t, y=Height_t,shape = factor(Surv_tplus1))) +
  scale_shape_manual(values=c(1, 16)) +
  xlab("") + ylab("Height at time t (cm)") +
  theme(text = element_text(size=10), plot.margin = margin(c(10,2,2,2))) + 
  labs(fill = "Predicted \n survival \n at time \n t + 1") +
  labs(shape = "Measured \n survival \n at time \n t + 1 (mm)") +
  theme(legend.title.align=0.5) +
  guides(shape = guide_legend(title.position = "right", order =1, text=15, 
                              label.theme=element_text(size=15), override.aes = list(size=5))) +
  guides(fill = guide_colorbar(title.position = "right", order =2, 
                               label.theme=element_text(size=15), text=15)) +
  guides(size=F)

legend <- get_legend(
  # create some space to the left of the legend
  BC + theme(legend.box.margin = margin(0, 0, 0, 0))
) 
prow <- plot_grid(
  BC + theme(aspect.ratio=1, legend.position="none", plot.margin=unit(c(0,0,0,0), "cm")) + labs(x=""),
  C + theme(aspect.ratio=1, legend.position="none", plot.margin=unit(c(0,0,0,0), "in"), 
             axis.text.y=element_blank()) + labs(y="", x=" Diameter at time t (mm)"),
  FP + theme(aspect.ratio=1, legend.position="none", plot.margin=unit(c(0,0,0,0), "cm"), 
              axis.text.y=element_blank()) + labs(y="", x=""),
  CC + theme(aspect.ratio=1, legend.position="none", plot.margin=unit(c(0,0,0,0), "cm"), 
             axis.text.y=element_blank()) + labs(y="", x=""), ncol=4
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
png(file="./PowerPoint/seedling_survival.png", width=11.64, 
    height=2.91, units='in', res=300)
plot_grid(prow, legend, 
          rel_widths=c(3,0.7))
dev.off()

#FIG 7: Population growth rate #####
lam.stable_BC <- readRDS("./BC/lam.stable_BC.rds")
lam.stable_CC <- readRDS("./CC/lam.stable_CC.rds")
lam.stable_C  <- readRDS("./C/lam.stable_C.rds")
lam.stable_FP <- readRDS("./FP/lam.stable_FP.rds")
lam.stable_PG <- readRDS("./PG/lam.stable_PG.rds")
lam.stable_WT <- readRDS("./WT/lam.stable_WT.rds")


x11(width=11, height=5)
png(filename="./Figures/lambdas_for_ppt.png",
    width=11, height=5, res=300)
par(mar = c(5, 3, 2, 1), # Margins,
    oma = c(4, 0, 0, 0),
    mgp = c(1.5, 2, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=F,
    mai=c(0.5,2,0.5,0.5),
    mfrow=c(1, 1))
barCenters <- plot(c(lam.stable_PG, lam.stable_WT, lam.stable_BC,
                     lam.stable_C, lam.stable_FP, lam.stable_CC), 
                   col=cols,
                   cex=3,
                   cex.axis=1.5,
                   ylab="",
                   axes= F,
                   xlab= "",
                   pch=19,
                   xpd=F, 
                   ylim =c(0.98, 1.26))
axis(side=1, labels = c("Punta \n Gorda", "Wild \n Turkey",
                        "Big \n Cypress", "Chekika",
                        "Fort \n Pierce", "Cape \n Canaveral"),
     at = 1:6, line=1.5, cex.axis=1.5)
mtext(side=1, "Site", line=5, cex=2)
axis(side=2, cex.axis=1.5)
abline(h=1, lty=2, lwd=2)
mtext(side=2, expression(paste(lambda)), line=5, cex=2)
segments(3, 1.062, 3, 1.206, lwd=2) #Big Cypress
arrows(3, 1.062, 3, 1.206, angle = 90, 
       code = 3, length = 0.05, lwd=2)

segments(6, 1.055, 6, 1.142, lwd=2) #Cape Canaveral
arrows(6, 1.055, 6, 1.142, angle = 90, 
       code = 3, length = 0.05, lwd=2)

segments(4, 0.986, 4, 1.05, lwd=2) #Chekika
arrows(4, 0.986, 4, 1.05, angle = 90, 
       code = 3, length = 0.05, lwd=2)

segments(5, 1.079, 5, 1.221, lwd=2) #Fort Pierce
arrows(5, 1.079, 5, 1.221, angle = 90, code = 3, length = 0.05, lwd=2)

segments(1, 0.995, 1, 1.140, lwd=2) #Punta Gorda
arrows(1, 0.995, 1, 1.140, angle=90, code = 3, length = 0.05, lwd=2)

segments(2, 1.062, 2, 1.251, lwd=2) #Wild Turkey
arrows(2, 1.062, 2, 1.251, angle = 90, code =3, length = 0.05, lwd=2)
dev.off()

#FIG 9 Barplot Contributions #####
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



BC <- c(BC_pool[1], sum(BC_pool[2:3]), sum(BC_pool[4:5]), sum(BC_pool[6:9]))
CC <- c(CC_pool[1], sum(CC_pool[2:3]), sum(CC_pool[4:5]), sum(CC_pool[6:9]))
C  <- c(C_pool[1],  sum(C_pool[2:3]),  sum(C_pool[4:5]),  sum(C_pool[6:9]))
FP <- c(FP_pool[1], sum(FP_pool[2:3]), sum(FP_pool[4:5]), sum(FP_pool[6:9]))
PG <- c(PG_pool[1], sum(PG_pool[2:3]), sum(PG_pool[4:5]), sum(PG_pool[6:9]))
WT <- c(WT_pool[1], sum(WT_pool[2:3]), sum(WT_pool[4:5]), sum(WT_pool[6:9]))

par(mar = c(3, 1, 1, 1), # Margins
    mgp = c(0, 0, 0), # Distance of axis tickmark labels (second value)
    tcl = 0.3, # Length of axis tickmarks
    xpd=F,
    #mai=c(0.2,.5,0.2,0.1),
    oma=c(1, 5, 1, 0),
    mfrow=c(1, 6))

#1) PG
xx <- barplot(PG, col=cols[1], ylim=c(-0.05, 0.13),
              axes=F, cex.names=0.6, cex.axis=1.5, horiz=F, 
              space=c(0.2,0.2))
axis(side=2, cex.axis=2)
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["site"],"-",lambda["pool"])))
mtext(do.call(expression, Lines),side=2,line=c(4, 1.5), cex=1.5)
mtext(side=1, "Punta \n Gorda", cex=1.5, line=2.5)
text(x=xx[1], y=PG[1]+0.02, label="Seedling", 
     cex=2, srt=90, col=cols[1])
text(x=xx[2], y=PG[2]+0.02, label="Maturation", 
     cex=2, srt=90, col=cols[1])
text(x=xx[3], y=0.015, label="Fertility",
     cex=2, srt=90, col=cols[1])
text(x=xx[4], y=0.01, label="Adult", 
     cex=2, srt=90, col=cols[1])
abline(h=0)


#2) WT
xx <- barplot(WT, col=cols[2], ylim=c(-0.05, 0.13),
              axes=F, cex.names=0.6, cex.axis=1.5, horiz=F, 
              space=c(0.2,0.2))
mtext(side=1, "Wild \n Turkey", cex=1.5, line=2.5)
text(x=xx[1], y=0.015, label="Seedling", 
     cex=2, srt=90, col=cols[2])
text(x=xx[2], y=WT[2]+0.02, label="Maturation", 
     cex=2, srt=90, col=cols[2])
text(x=xx[3], y=0.015, label="Fertility",
     cex=2, srt=90, col=cols[2])
text(x=xx[4], y=WT[4]+0.01, label="Adult", 
     cex=2, srt=90, col=cols[2])
abline(h=0)

#3) BC
xx <- barplot(BC, col=cols[3], ylim=c(-0.05, 0.13),
              axes=F, cex.names=0.6, cex.axis=1.5, horiz=F, 
              space=c(0.2,0.2))
mtext(side=1, "Big \n Cypress", cex=1.5, line=2.5)
text(x=xx[1], y=0.015, label="Seedling", 
     cex=2, srt=90, col=cols[3])
text(x=xx[2], y=0.02, label="Maturation", 
     cex=2, srt=90, col=cols[3])
text(x=xx[3], y=BC[3]+0.015, label="Fertility",
     cex=2, srt=90, col=cols[3])
text(x=xx[4], y=0.01, label="Adult", 
     cex=2, srt=90, col=cols[3])
abline(h=0)

#)4 C
xx <- barplot(C, col=cols[4], ylim=c(-0.05, 0.13),
              axes=F, cex.names=0.6, cex.axis=1.5, horiz=F, 
              space=c(0.2,0.2))
mtext(side=1, "Chekika", cex=1.5, line=2.5)
text(x=xx[1], y=0.015, label="Seedling", 
     cex=2, srt=90, col=cols[4])
text(x=xx[2], y=0.02, label="Maturation", 
     cex=2, srt=90, col=cols[4])
text(x=xx[3], y=0.015, label="Fertility",
     cex=2, srt=90, col=cols[4])
text(x=xx[4], y=0.01, label="Adult", 
     cex=2, srt=90, col=cols[4])
abline(h=0)

#5) FP
xx <- barplot(FP, col=cols[5], ylim=c(-0.05, 0.13),
              axes=F, cex.names=0.6, cex.axis=1.5, horiz=F, 
              space=c(0.2,0.2))
mtext(side=1, "Fort \n Pierce", cex=1.5, line=2.5)
text(x=xx[1], y=FP[1]+0.015, label="Seedling", 
     cex=2, srt=90, col=cols[5])
text(x=xx[2], y=FP[2]+0.02, label="Maturation", 
     cex=2, srt=90, col=cols[5])
text(x=xx[3], y=0.015, label="Fertility",
     cex=2, srt=90, col=cols[5])
text(x=xx[4], y=0.01, label="Adult", 
     cex=2, srt=90, col=cols[5])
abline(h=0)

#6)CC
xx <- barplot(CC, col=cols[6], ylim=c(-0.05, 0.13),
              axes=F, cex.names=0.6, cex.axis=1.5, horiz=F, 
              space=c(0.2,0.2))
mtext(side=1, "Cape \n Canaveral", cex=1.5, line=2.5)
text(x=xx[1], y=CC[1]+0.015, label="Seedling", 
     cex=2, srt=90, col=cols[6])
text(x=xx[2], y=CC[2]+0.02, label="Maturation", 
     cex=2, srt=90, col=cols[6])
text(x=xx[3], y=0.015, label="Fertility",
     cex=2, srt=90, col=cols[6])
text(x=xx[4], y=0.01, label="Adult", 
     cex=2, srt=90, col=cols[6])
abline(h=0)
dev.off()
