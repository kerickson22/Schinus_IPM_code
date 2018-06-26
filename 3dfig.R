#Default width requested by journal: 
width.cm <- 12.9
height.cm <- 6.45
pointsize <- 8
x11(width = width.cm/2.54, height = height.cm/2.54, 
    pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.3,0.3,0.1,0.1))

#####################################################################
# Function for nicely plotting 3D surface 
survival_plot<-function(b0, b1, b2,max.x, max.y, max.z, z, x, y, zlim,
                        xlab, ylab, zlab ) {
  
  
  z.axis<-seq(0, max.z, length.out=4)
  x.axis<-seq(0, max.z, length.out=4)
  y.axis<-seq(0, max.y, length.out=4)
  
  
  pmat<-persp(x, y, z, ticktype="detailed",
              xlab="", ylab="", zlab="", zlim=c(0, zlim),
              box=T, axes=FALSE, 
              mar=c(10, 1, 0, 2), 
              shade=0.25, theta=angle)
  
  
  tick.start <- trans3d(x.axis, 0, 0, pmat)
  tick.end <- trans3d(x.axis, (0 - 0.20), 0, pmat)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
  
  tick.start <- trans3d(0, y.axis, 0, pmat)
  tick.end <- trans3d(0 - 0.02, y.axis, 0, pmat)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
  
  tick.start <- trans3d(0, max.y, z.axis, pmat)
  tick.end <- trans3d(0, (max.y + 1), z.axis, pmat)
  segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
  
  labels <- as.character(round(x.axis, 1))
  label.pos <- trans3d(x.axis, (0 - 0.75), 0, pmat)
  text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1)
  
  xlabel.pos <- trans3d(max.x/2, (min.y-2.5), min.z, pmat)
  text(xlabel.pos$x, xlabel.pos$y, labels=xlab,
       cex=1, srt=20)
  
  labels <- as.character(round(y.axis, 1))
  ylabel.pos <- trans3d(0, y.axis, 0, pmat)
  text(ylabel.pos$x, ylabel.pos$y, labels=labels, adj=c(1.75, NA), cex=1)
  
  ylabel.pos<-trans3d(0-.3, max.y/2, 0, pmat)
  text(ylabel.pos$x, ylabel.pos$y, labels=ylab, 
       cex=1, srt=-60)
  
  labels <- as.character(round(z.axis, 1))
  zlabel.pos <- trans3d(0, max.y+2.5, z.axis, pmat)
  text(zlabel.pos$x, zlabel.pos$y, labels=labels, adj=c(1, NA), cex=1)
  
  zlabel.pos<- trans3d(0-0.3, max.y+2.5, max.z/2, pmat)
  text(zlabel.pos$x, zlabel.pos$y, labels="P(Survival)",
       cex=1, srt=-85)
}
##########################################################################

b0<-p.vec_E[20]
b1<-p.vec_E[21]
b2<-p.vec_E[22]

min.x<-0
max.x<-1.5
min.y<-0
max.y<-15
min.z<-0
max.z<-1
zlim<-1.2

z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))


x11(width = width.cm/2.54, height = height.cm/2.54, 
    pointsize = pointsize)

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.3,0.3,0.1,0.1),
    mfrow=c(1, 2))
survival_plot(b0, b1, b2,max.x, max.y, max.z, z1, x1seq_seedlings,
              x2seq_seedlings, 1.2,
              "Diameter (mm)", "Height (cm)", "P(Survival)")
survival_plot(b0, b1, b2, max.x, max.y, max.z, z1, x1seq_seedlings,
              x2seq_seedlings, 1.2, 
              "Diameter (mm)", "Height (cm)", "P(Survival)")

#This creates an EPS figure! That LaTeX can read
setEPS(horizontal=F, onefile=F, paper="special")
postscript("./Figures/example5.ps", width=width.cm/2.54, 
           height=height.cm/2.54, pointsize=pointsize,  encoding = "TeXtext.enc")

par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.3,0.3,0.1,0.1),
    mfrow=c(1, 2))
survival_plot(b0, b1, b2,max.x, max.y, max.z, z1, x1seq_seedlings,
              x2seq_seedlings, 1.2,
              "Diameter (mm)", "Height (cm)", "P(Survival)")
survival_plot(b0, b1, b2, max.x, max.y, max.z, z1, x1seq_seedlings,
              x2seq_seedlings, 1.2, 
              "Diameter (mm)", "Height (cm)", "P(Survival)")


dev.off()




x11(width=width.cm/2.54, height=height.cm/2.54, pointsize = pointsize )
par(mar = c(3, 3, 2, 1), # Margins
    mgp = c(1.5, 0.5, 0), # Distance of axis tickmark labels (second value)
    tcl = -0.3, # Length of axis tickmarks
    xpd=NA,
    mai=c(0.3,0.3,0.1,0.1),
    mfrow=c(1, 2))


survival_plot(b0, b1, b2, max.x, max.y, max.z, z1, x1seq_seedlings, 
              x2seq_seedlings, 1.2, "Diameter (mm)", "Height (cm)", 
              "P(Survival)") 



















b0<-p.vec_E[20]
b1<-p.vec_E[21]
b2<-p.vec_E[22]

min.x<-0
max.x<-1.5
min.y<-0
max.y<-15
min.z<-0
max.z<-1
z.axis<-seq(0, 1, by=.25)
x.axis<-seq(0, 1.5, length.out=4)
y.axis<-seq(0, 15, length.out=4)

z1<-outer(x1seq_seedlings, x2seq_seedlings, function(a,b) ((exp(b0+b1*a +b2*b)/(1+exp(b0+b1*a +b2*b)))))

#bottom, left, top, right 
pmat<-persp(x1seq_seedlings, x2seq_seedlings, z1, ticktype="detailed",
      xlab="", ylab="", zlab="", zlim=c(0, 1.2),
      box=T, axes=FALSE, 
      mar=c(10, 1, 0, 2), 
      shade=0.25, theta=angle)

#lines(trans3d(x1seq_seedlings, min.y, min.z, pmat) , col="black")
#lines(trans3d(min.x, x2seq_seedlings, min.z, pmat) , col="black")
#lines(trans3d(min.x, min.y, z.axis, pmat) , col="black")
#lines(trans3d(min.x, max.y, z.axis, pmat) , col="black")

tick.start <- trans3d(x.axis, min.y, min.z, pmat)
tick.end <- trans3d(x.axis, (min.y - 0.20), min.z, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
            
tick.start <- trans3d(min.x, y.axis, min.z, pmat)
tick.end <- trans3d(min.x - 0.02, y.axis, min.z, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

tick.start <- trans3d(min.x, max.y, z.axis, pmat)
tick.end <- trans3d(min.x, (max.y + 1), z.axis, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

labels <- as.character(round(x.axis, 1))
label.pos <- trans3d(x.axis, (min.y - 0.75), min.z, pmat)
text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1)

xlabel.pos <- trans3d(max.x/2, (min.y-2.5), min.z, pmat)
text(xlabel.pos$x, xlabel.pos$y, labels="Diameter at time t",
     cex=1, srt=20)

labels <- as.character(round(y.axis, 1))
ylabel.pos <- trans3d(min.x, y.axis, min.z, pmat)
text(ylabel.pos$x, ylabel.pos$y, labels=labels, adj=c(1.75, NA), cex=1)

ylabel.pos<-trans3d(min.x-.3, max.y/2, min.z, pmat)
text(ylabel.pos$x, ylabel.pos$y, labels="Height at time t", 
     cex=1, srt=-60)

labels <- as.character(round(z.axis, 1))
zlabel.pos <- trans3d(min.x, max.y+2.5, z.axis, pmat)
text(zlabel.pos$x, zlabel.pos$y, labels=labels, adj=c(1, NA), cex=1)

zlabel.pos<- trans3d(min.x-0.3, max.y+2.5, max.z/2, pmat)
text(zlabel.pos$x, zlabel.pos$y, labels="P(Survival)",
     cex=1, srt=-85)
#dev.off()

filename <- "example2"
filetype <- "pdf" 
savePlot(filename=filename, type=filetype)
