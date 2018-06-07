rm(list=ls())

#LTRE
m1<-10
m2<-m1+1
m3<-100
m4<-m3+1

h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
h3=(700-1.6)/m3;
y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
h4=(800-16)/m4
 y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2


#####Load objects that were previously run: 
# Have to do in pieces and remove strategically, because it is 12GB of data

################
###OVERALL     #	
################	
	setwd("/Users/kelley/Dropbox/matrix/outputs/m1=10/m3=100")

	#Load A matrix:
	load("D1.RData")	
	load("G.RData")		
	load("F.RData")
	load("D2.RData")
	A_all<-cbind(rbind(D1, G), rbind(F, D2))
		
	#tidy up workspace:
	rm(D1)		
	rm(G)
	rm(F)
	rm(D2)

	#Load lambda:
	load("lam.stable.RData")
		lambda_all<-lam.stable
		rm(lam.stable)

	#Calculate sensitivity matrix: 
	load("stable.dist.RData")
	load("v.RData")

	v.dot.w<-sum(t(v)%*%stable.dist)
	norm_v<-v/v.dot.w
	sens_all<-norm_v%*%t(stable.dist)
	
	elas_all<-sens_all*A_all/lambda_all



	rm(stable.dist)
	rm(v)
	rm(v.dot.w)
	rm(norm_v)


################
###EASTERN     #	
################	
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Eastern")

	#Load A matrix:
	load("D1.RData")	
	load("G.RData")		
	load("F.RData")
	load("D2.RData")
	A_E<-cbind(rbind(D1, G), rbind(F, D2))
		
	#tidy up workspace:
	rm(D1)		
	rm(G)
	rm(F)
	rm(D2)
	
	#Calculate difference matrix: D_E<-A_E - A_all
	D_E<-A_E-A_all
	rm(A_E)
	
	#Calculate matrix of contributions: C_E = D_E %*% sens_all
	C_E<-D_E*sens_all
	rm(D_E)
	

	#Load lambda:
	load("lam.stable.RData")
		lambda_E<-lam.stable
		rm(lam.stable)

lambda_all+sum(C_E)

	#Predict lambda_E<- lambda_all + sum(C_E)
	
C_D1_E<-C_E[1:(m1*m2), 1:(m1*m2)]
C_G_E<-C_E[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
C_F_E<-C_E[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
C_D2_E<-C_E[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

	#Unpack matrix of contributions into pieces: 
	
	#D1#
	plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in 	array
	Plop=outer(1:m1,1:m2,plop); 

	Kvals_contrib_D1_E=array(0,c(m1,m2,m1,m2));  

	for(i in 1:m1){
		for(j in 1:m2){
		for(k in 1:m1){
				kvals= C_D1_E[Plop[k,1:m2],Plop[i,j]]
				Kvals_contrib_D1_E[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_contrib_D2_E=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= C_D2_E[Plop[k,1:m4],Plop[i,j]]
				Kvals_contrib_D2_E[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_contrib_F_E=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=C_F_E[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_contrib_F_E[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_contrib_G_E=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=C_G_E[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_contrib_G_E[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


##############
total.contrib_D1_E<-apply(Kvals_contrib_D1_E, c(3,4), sum);
total.contrib_D2_E<-apply(Kvals_contrib_D2_E, c(3,4), sum);
total.contrib_F_E<-apply(Kvals_contrib_F_E, c(3,4), sum);
total.contrib_G_E<-apply(Kvals_contrib_G_E, c(3,4), sum);
	
	
setwd("/Users/kelley/Dropbox/matrix/outputs/LTRE")	

jpeg('Contributions_to_lambda_height_E.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y2, colSums(total.contrib_D1_E), ylim=c(-0.0015, 0.0001), xlab="Height (cm)", ylab="Contribution", main="D1-E",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.contrib_F_E), ylim=c(-0.0015, 0.0001), xlab="Height (cm)", ylab= "Contribution", main="F-E", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y2, colSums(total.contrib_G_E), ylim=c(-0.0015, 0.0001), xlab="Height (cm)", ylab="Contribution", main="G-E", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.contrib_D2_E), ylim=c(-0.0015, 0.0001), xlab="Height (cm)", ylab="Contribution", main="D2-E",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()

jpeg('Contributions_to_lambda_diam_E.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y1, rowSums(total.contrib_D1_E), ylim=c(-0.0015, 0.0001), xlab="Diameter (mm)", ylab="Contribution", main="D1-E",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.contrib_F_E), ylim=c(-0.0015, 0.0001), xlab="Diameter (mm)", ylab= "Contribution", main="F-E", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y1, rowSums(total.contrib_G_E), ylim=c(-0.0015, 0.0001), xlab="Diameter (mm)", ylab="Contribution", main="G-E", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.contrib_D2_E), ylim=c(-0.0015, 0.0001), xlab="Diameter (mm)", ylab="Contribution", main="D2-E",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()
	
	
	
	
################
###Hybrid     #	
################	
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Hybrid")

	
	#Load A matrix:
	load("D1.RData")	
	load("G.RData")		
	load("F.RData")
	load("D2.RData")
	A_H<-cbind(rbind(D1, G), rbind(F, D2))
		
	#tidy up workspace:
	rm(D1)		
	rm(G)
	rm(F)
	rm(D2)
	
	#Calculate difference matrix: D_H<-A_H - A_all
	D_H<-A_H-A_all
	rm(A_H)
	
	#Calculate matrix of contributions: C_E = D_H * sens_all
	C_H<-D_H*sens_all
	rm(D_H)
	

	#Load lambda:
	load("lam.stable.RData")
		lambda_H<-lam.stable
		rm(lam.stable)

lambda_all+sum(C_H)

	#Predict lambda_H<- lambda_all + sum(C_H)
	
C_D1_H<-C_H[1:(m1*m2), 1:(m1*m2)]
C_G_H<-C_H[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
C_F_H<-C_H[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
C_D2_H<-C_H[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

	#Unpack matrix of contributions into pieces: 
	
	#D1#
	plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in 	array
	Plop=outer(1:m1,1:m2,plop); 

	Kvals_contrib_D1_H=array(0,c(m1,m2,m1,m2));  

	for(i in 1:m1){
		for(j in 1:m2){
		for(k in 1:m1){
				kvals= C_D1_H[Plop[k,1:m2],Plop[i,j]]
				Kvals_contrib_D1_H[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_contrib_D2_H=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= C_D2_H[Plop[k,1:m4],Plop[i,j]]
				Kvals_contrib_D2_H[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_contrib_F_H=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=C_F_H[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_contrib_F_H[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_contrib_G_H=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=C_G_H[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_contrib_G_H[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


##############
total.contrib_D1_H<-apply(Kvals_contrib_D1_H, c(3,4), sum);
total.contrib_D2_H<-apply(Kvals_contrib_D2_H, c(3,4), sum);
total.contrib_F_H<-apply(Kvals_contrib_F_H, c(3,4), sum);
total.contrib_G_H<-apply(Kvals_contrib_G_H, c(3,4), sum);
	
	
setwd("/Users/kelley/Dropbox/matrix/outputs/LTRE")	

jpeg('Contributions_to_lambda_height_H.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y2, colSums(total.contrib_D1_H), ylim=c(-0.0015, 0.002), xlab="Height (cm)", ylab="Contribution", main="D1-H",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.contrib_F_H), ylim=c(-0.0015, 0.002), xlab="Height (cm)", ylab= "Contribution", main="F-H", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y2, colSums(total.contrib_G_H), ylim=c(-0.0015, 0.002), xlab="Height (cm)", ylab="Contribution", main="G-H", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.contrib_D2_H), ylim=c(-0.0015, 0.002), xlab="Height (cm)", ylab="Contribution", main="D2-H",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()

jpeg('Contributions_to_lambda_diam_H.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y1, rowSums(total.contrib_D1_H), ylim=c(-0.0015, 0.002), xlab="Diameter (mm)", ylab="Contribution", main="D1-H",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.contrib_F_H), ylim=c(-0.0015, 0.002), xlab="Diameter (mm)", ylab= "Contribution", main="F-H", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y1, rowSums(total.contrib_G_H), ylim=c(-0.0015, 0.002), xlab="Diameter (mm)", ylab="Contribution", main="G-H", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.contrib_D2_H), ylim=c(-0.0015, 0.002), xlab="Diameter (mm)", ylab="Contribution", main="D2-H",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()
	
	

################
###Western     #	
################	
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Western")

	
	#Load A matrix:
	load("D1.RData")	
	load("G.RData")		
	load("F.RData")
	load("D2.RData")
	A_W<-cbind(rbind(D1, G), rbind(F, D2))
		
	#tidy up workspace:
	rm(D1)		
	rm(G)
	rm(F)
	rm(D2)
	
	#Calculate difference matrix: D_W<-A_W - A_all
	D_W<-A_W-A_all
	rm(A_W)
	
	#Calculate matrix of contributions: C_W = D_W * sens_all
	C_W<-D_W*sens_all
	rm(D_W)
	

	#Load lambda:
	load("lam.stable.RData")
		lambda_W<-lam.stable
		rm(lam.stable)

lambda_all+sum(C_W)

	#Predict lambda_W<- lambda_all + sum(C_W)
	
C_D1_W<-C_W[1:(m1*m2), 1:(m1*m2)]
C_G_W<-C_W[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
C_F_W<-C_W[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
C_D2_W<-C_W[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

	#Unpack matrix of contributions into pieces: 
	
	#D1#
	plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in 	array
	Plop=outer(1:m1,1:m2,plop); 

	Kvals_contrib_D1_W=array(0,c(m1,m2,m1,m2));  

	for(i in 1:m1){
		for(j in 1:m2){
		for(k in 1:m1){
				kvals= C_D1_W[Plop[k,1:m2],Plop[i,j]]
				Kvals_contrib_D1_W[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_contrib_D2_W=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= C_D2_W[Plop[k,1:m4],Plop[i,j]]
				Kvals_contrib_D2_W[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_contrib_F_W=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=C_F_W[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_contrib_F_W[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_contrib_G_W=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=C_G_W[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_contrib_G_W[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


##############
total.contrib_D1_W<-apply(Kvals_contrib_D1_W, c(3,4), sum);
total.contrib_D2_W<-apply(Kvals_contrib_D2_W, c(3,4), sum);
total.contrib_F_W<-apply(Kvals_contrib_F_W, c(3,4), sum);
total.contrib_G_W<-apply(Kvals_contrib_G_W, c(3,4), sum);
	
	
setwd("/Users/kelley/Dropbox/matrix/outputs/LTRE")	

jpeg('Contributions_to_lambda_height_W.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y2, colSums(total.contrib_D1_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab="Contribution", main="D1-W",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.contrib_F_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab= "Contribution", main="F-W", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y2, colSums(total.contrib_G_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab="Contribution", main="G-W", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.contrib_D2_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab="Contribution", main="D2-W",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()

jpeg('Contributions_to_lambda_diam_W.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y1, rowSums(total.contrib_D1_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab="Contribution", main="D1-W",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.contrib_F_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab= "Contribution", main="F-W", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y1, rowSums(total.contrib_G_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab="Contribution", main="G-W", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.contrib_D2_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab="Contribution", main="D2-W",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()


###Putting it all together on one figure...

x11()
jpeg('Contributions_to_lambda_height_all.jpg', quality=500, width=1000, height=1000)
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=24)
plot(y2, colSums(total.contrib_D1_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab="Contribution", main="D1",  col="red", type='l', cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y2, colSums(total.contrib_D1_H), col="purple", type='l', lwd=2)
lines(y2, colSums(total.contrib_D1_E), col="blue", type='l', lwd=2)

plot(y4, colSums(total.contrib_F_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab= "Contribution", main="F", type='l', col="red", cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y4, colSums(total.contrib_F_H), col="purple", type='l', lwd=2)
lines(y4, colSums(total.contrib_F_E), col="blue", type='l', lwd=2)

plot(y2, colSums(total.contrib_G_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab="Contribution", col="red", main="G", type='l', cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y2, colSums(total.contrib_G_H), col="purple", type='l', lwd=2)
lines(y2, colSums(total.contrib_G_E), col="blue", type='l', lwd=2)

plot(y4, colSums(total.contrib_D2_W), ylim=c(-0.0015, 0.005), xlab="Height (cm)", ylab="Contribution", main="D2",  col="red", type='l',  cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y4, colSums(total.contrib_D2_H), col="purple", type='l', lwd=2)
lines(y4, colSums(total.contrib_D2_E), col="blue", type="l", lwd=2)
legend("topright", lty=1, lwd=2, col=c("red", "purple", "blue"), c("Western", "Hybrid", "Eastern"), bty='n', cex=0.8, y.intersp=2)
dev.off()


#For Diameter: 

jpeg('Contributions_to_lambda_diameter_all.jpg', quality=500, width=1000, height=1000)
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=24)
plot(y1, rowSums(total.contrib_D1_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab="Contribution", main="D1",  col="red", type='l', cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y1, rowSums(total.contrib_D1_H), col="purple", type='l', lwd=2)
lines(y1, rowSums(total.contrib_D1_E), col="blue", type='l', lwd=2)

plot(y3, rowSums(total.contrib_F_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab= "Contribution", main="F", type='l', col="red", cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y3, rowSums(total.contrib_F_H), col="purple", type='l', lwd=2)
lines(y3, rowSums(total.contrib_F_E), col="blue", type='l', lwd=2)

plot(y1, rowSums(total.contrib_G_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab="Contribution", col="red", main="G", type='l', cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y1, rowSums(total.contrib_G_H), col="purple", type='l', lwd=2)
lines(y1, rowSums(total.contrib_G_E), col="blue", type='l', lwd=2)

plot(y3, rowSums(total.contrib_D2_W), ylim=c(-0.0015, 0.005), xlab="Diameter (mm)", ylab="Contribution", main="D2",  col="red", type='l',  cex=1, cex.axis=1, cex.lab=1, lwd=2)
lines(y3, rowSums(total.contrib_D2_H), col="purple", type='l', lwd=2)
lines(y3, rowSums(total.contrib_D2_E), col="blue", type="l", lwd=2)
legend("topright", lty=1, lwd=2, col=c("red", "purple", "blue"), c("Western", "Hybrid", "Eastern"), bty='n', cex=0.8, y.intersp=2)
dev.off()



###NEXT STEP: 

lambda_diff_E<-lambda_E-lambda_all
lambda_diff_H<-lambda_H-lambda_all
lambda_diff_W<-lambda_W-lambda_all

barplot(c(lambda_diff_W, lambda_diff_H, lambda_diff_E), col=c("red", "purple", "blue"))

W<-c( 
sum(rowSums(total.contrib_D1_W)), 
sum(rowSums(total.contrib_G_W)), 
sum(rowSums(total.contrib_F_W)), 
sum(rowSums(total.contrib_D2_W)) )

H<-c(
sum(rowSums(total.contrib_D1_H)), 
sum(rowSums(total.contrib_G_H)), 
sum(rowSums(total.contrib_F_H)), 
sum(rowSums(total.contrib_D2_H)) )

E<-c( 
sum(rowSums(total.contrib_D1_E)), 
sum(rowSums(total.contrib_G_E)), 
sum(rowSums(total.contrib_F_E)), 
sum(rowSums(total.contrib_D2_E)) )


x11()

png(file='lambdas_plot.png', width=10, height=10, units="in", res=300)

par(mar=c(5.1, 5.1, 2.1, 0), oma=c(1, 0, 0, 0), ps=24)
barplot(c(lambda_all, lambda_E, lambda_H, lambda_W), col=c("black", "blue", "purple", "red"), ylab=expression(paste(lambda)),  names.arg=c("Overall", "Eastern", "Hybrid", "Western"), ylim=c(1, 1.1),xpd=F)
dev.off()

#text(c(lambda_all, lambda_E, lambda_H, lambda_W), c(1.09, 1.09, 1.09, 1.09), c("1.08", "1.05", "1.08", "1.09"))


png(file='barplot_contributions.png', width=4, height=6, units="in", res=300)
par(mfrow=c(3, 1), mar=c(2.1,7.1,2.1,0), oma = c(1,0,0,0),
ps=24)

barplot(W, col=c("red", "red", "red"), ylim=c(-0.015, 0.03), names.arg=c("D1", "G", "F", "D2"))
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["W"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(5, 3), cex=.6)
abline(h=0)

barplot(H, col=c("purple", "purple", "purple"), ylim=c(-0.015, 0.03), names.arg=c("D1", "G", "F", "D2"))
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["H"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(5, 3), cex=.6)
abline(h=0)

barplot(E, col=c("blue", "blue", "blue"), ylim=c(-0.015, 0.03), names.arg=c("D1", "G", "F", "D2"))
Lines <- list(bquote(paste( "Contribution to" )),
              bquote(paste(lambda["E"],"-", lambda["overall"])))
mtext(do.call(expression, Lines),side=2,line=c(5, 3), cex=.6)
abline(h=0)
dev.off()

x11()


jpeg("contributions_barplot.jpeg", width=1000, height=500, quality=500)
barplot(c(W, H, E), col=c("red", "red", "red", "red", "red", "purple", "purple","purple","purple","purple", "blue","blue","blue","blue","blue"), density=c(100, 10, 10, 10, 10, 100, 10, 10, 10, 10, 100, 10 ,10 ,10, 10), names.arg=c("W", "W:D1", "W:G", "W:F", "W:D2", "H", "H:D1", "H:G", "H:F", "H:D2", 
"E", "E:D1", "E:G", "E:F", "E:D2"), ylab="Lambda_Biotype - Lambda_overall")
dev.off()


#######################################
#
#CV * elas 
#
#######################################

#Eastern
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Eastern")
	load("D1.RData")
	D1_E<-D1
	rm(D1)
	
#Hybrid
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Hybrid")
	load("D1.RData")
	D1_H<-D1
	rm(D1)

#Western
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Western")
	load("D1.RData")
	D1_W<-D1
	rm(D1)
	
	
CV<-function(x) {
	return(sd(x)/mean(x))
	
}



	
cv_D1<-matrix(, nrow=(m1*m2), ncol=(m1*m2))

for(i in 1:(m1*m2)) {
	for(j in 1:(m1*m2)) {
		cv_D1[i,j]<-CV(c(D1_E[i,j], D1_W[i,j], D1_H[i,j]))
		if(is.na(cv_D1[i,j])) {cv_D1[i,j]<-0}
	}
}

D1_avg=(D1_E+D1_H+D1_W)/3
rm(D1_E)
rm(D1_H)
rm(D1_W)
	
#Eastern
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Eastern")
	load("G.RData")
	G_E<-G
	rm(G)
	
#Hybrid
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Hybrid")
	load("G.RData")
	G_H<-G
	rm(G)

#Western
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Western")
	load("G.RData")
	G_W<-G
	rm(G)

cv_G<-matrix(, nrow=(m3*m4), ncol=(m1*m2))


for(i in 1:(m3*m4)) {
	for(j in 1:(m1*m2)) {
		cv_G[i,j]<-CV(c(G_E[i,j], G_W[i,j], G_H[i,j]))
		if(is.na(cv_G[i,j])) {cv_G[i,j]<-0}
	}
}

G_avg<-(G_E+G_H+G_W)/3
rm(G_E)
rm(G_H)
rm(G_W)

#Eastern
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Eastern")
	load("F.RData")
	F_E<-F
	rm(F)
	
#Hybrid
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Hybrid")
	load("F.RData")
	F_H<-F
	rm(F)

#Western
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Western")
	load("F.RData")
	F_W<-F
	rm(F)

cv_F<-matrix(, nrow=(m1*m2), ncol=(m3*m4))
for(i in 1:(m1*m2)) {
	for(j in 1:(m3*m4)) {
		cv_F[i,j]<-CV(c(F_E[i,j], F_W[i,j], F_H[i,j]))
		if(is.na(cv_F[i,j])) {cv_F[i,j]<-0}
	}
}

F_avg<-(F_E + F_W + F_H)/3

rm(F_E)
rm(F_W)
rm(F_H)

#Eastern
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Eastern")
	load("D2.RData")
	D2_E<-D2
	rm(D2)
	
#Hybrid
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Hybrid")
	load("D2.RData")
	D2_H<-D2
	rm(D2)

#Western
	setwd("/Users/kelley/Dropbox/matrix/outputs/sites/Western")
	load("D2.RData")
	D2_W<-D2
	rm(D2)

cv_D2<-matrix(, nrow=(m3*m4), ncol=(m3*m4))
for(i in 1:(m3*m4)) {
	for(j in 1:(m3*m4)) {
		cv_D2[i,j]<-CV(c(D2_E[i,j], D2_W[i,j], D2_H[i,j]))
		if(is.na(cv_D2[i,j])) {cv_D2[i,j]<-0}
	}
}

D2_avg<-(D2_E+D2_H+D2_W)/3
rm(D2_E)
rm(D2_W)
rm(D2_H)


A_cv<-cbind(rbind(cv_D1, cv_G), rbind(cv_F, cv_D2))
rm(cv_D1, cv_G, cv_F, cv_D2)

thing<-A_cv*elas_all

setwd("/Users/kelley/Dropbox/matrix/outputs/LTRE/")
save(thing, file="thing.RData")



thing_D1<-thing[1:(m1*m2), 1:(m1*m2)]
thing_G<-thing[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
thing_F<-thing[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
thing_D2<-thing[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

	#Unpack matrix of contributions into pieces: 
	
	#D1#
	plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in 	array
	Plop=outer(1:m1,1:m2,plop); 

	Kvals_thing_D1=array(0,c(m1,m2,m1,m2));  

	for(i in 1:m1){
		for(j in 1:m2){
		for(k in 1:m1){
				kvals= thing_D1[Plop[k,1:m2],Plop[i,j]]
				Kvals_thing_D1[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_thing_D2=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= thing_D2[Plop[k,1:m4],Plop[i,j]]
				Kvals_thing_D2[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_thing_F=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=thing_F[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_thing_F[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_thing_G=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=thing_G[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_thing_G[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


##############
total.thing_D1<-apply(Kvals_thing_D1, c(3,4), sum);
total.thing_D2<-apply(Kvals_thing_D2, c(3,4), sum);
total.thing_F<-apply(Kvals_thing_F, c(3,4), sum);
total.thing_G<-apply(Kvals_thing_G, c(3,4), sum);
	
	
setwd("/Users/kelley/Dropbox/matrix/outputs/LTRE")	

jpeg('thing_height.jpg')

par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=24)
plot(y2, colSums(total.thing_D1), ylim=c(0, 0.003), xlab="Height (cm)", ylab="CV * elas", main="D1",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.thing_F), ylim=c(0, 0.003), xlab="Height (cm)", ylab= "CV * elas", main="F", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y2, colSums(total.thing_G, na.rm=T), ylim=c(0, 0.003), xlab="Height (cm)", ylab="CV * elas", main="G", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y4, colSums(total.thing_D2), ylim=c(0, 0.003), xlab="Height (cm)", ylab="CV * elas", main="D2",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()

jpeg('thing_diam.jpg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=20)
plot(y1, rowSums(total.thing_D1), ylim=c(0, 0.003), xlab="Diameter (mm)", ylab="CV * elas", main="D1",  type='l', cex=1, cex.axis=0.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.thing_F), ylim=c(0, 0.003), xlab="Diameter (mm)", ylab= "CV * elas", main="F", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y1, rowSums(total.thing_G, na.rm=T), ylim=c(0, 0.003), xlab="Diameter (mm)", ylab="CV * elas", main="G", type='l', cex=1, cex.axis=.5, cex.lab=1, lwd=2)
plot(y3, rowSums(total.thing_D2), ylim=c(0, 0.003), xlab="Diameter (mm)", ylab="CV * sens", main="D2",  type='l',  cex=1, cex.axis=.5, cex.lab=1, lwd=2)
dev.off()


######################################################
#
###Some mysteries to solve: 
# SD and mean matrix (of A_E, A_H, and A_W) both have unique entries, but when 
# divided SD by mean, suddenly have columns that are identical in blocks. Why?
##


A_avg<-(A_E+A_H+A_W)/3
rm(A_E)
rm(A_H)
rm(A_W)
A_diff<-A_avg-A_all


persp(A_diff, theta=-40, xlab="Time t", ylab="Time t + 1", zlab="A_avg-A_all)")

x11()
persp(A_diff, theta=-90, xlab="Time t", ylab="Time t + 1", zlab="A_avg-A_all)")


#Plot just the F component: 

x11()

png(file="A_avg-A_overall.png", width=10, height=10, units="in", res=300)
par(ps=24)
persp(A_diff[1:110, 111:10210], theta=150, phi=20, ylab="Size at time t", xlab="Size at time t + 1", zlab="A_average - A_overall")
dev.off()


#persp(x1seq_seedlings, x2seq_seedlings, z_sdling_grad, ticktype="detailed", theta=-40, zlim=c(0, 1.02), xlab="Diameter (mm)", ylab="Height (cm)", zlab="P(graduating)", col = color[facetcol])



