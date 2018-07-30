# Constructs a integral projection model using diameter and height as two continuous
# structuring state variables, where population transition probabilities are
# modeled on two size domains. Modified code from Ellner and Rees, 2006. 



#clear everything
rm(list=ls(all=TRUE))

require(Matrix);
library(MASS)
library(mvtnorm)
library(msm)

 

# Set matrix size (to show up errors) and convergence tolerance. 
# Make m1 and m2 smaller to make this run faster.  
m1=10
m2=m1+1
m3=100
m4=m3+1
tol=1.e-8; 




# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(200),...) {
	nx=length(x); ny=length(y); 
	x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
	y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}
#============================================================================================#
# (I) Parameters and demographic functions for computing the kernel 
#============================================================================================#

#probability of seedling establishment
    

# rounded parameter vector, see Table 1 for more details
p.vec<-rep(0,39);
p.vec[1]<- surv_mod_larges$coeff[1]			#Intercept survival for large domain
p.vec[2]<- surv_mod_larges$coeff[2]			#diam slope survival for large domain
p.vec[3]<- surv_mod_larges$coeff[3]			#height slope survival for large domain
p.vec[4]<- growth_mod_larges$coeff[1,1]		#intercept growth diam for large domain
p.vec[5]<- growth_mod_larges$coeff[2,1]		#diam slope growth diam for large domain
p.vec[6]<- growth_mod_larges$coeff[3,1]		#height slope growth diam for large domain
p.vec[7]<- sigma_growth_diam_larges			#sd growth diam for large domain
p.vec[8]<- growth_mod_larges$coeff[1,2]		#intercept growth height for large domain
p.vec[9]<- growth_mod_larges$coeff[2,2]		#diam slope growth height for large domain 
p.vec[10]<- growth_mod_larges$coeff[3,2]		#height slope growth height for large domain
p.vec[11]<- sigma_growth_height_larges		#sd growth height for large domain 
p.vec[12]<-	mod_repro$coeff[1]				#intercept repro 
p.vec[13]<-	mod_repro$coeff[2]				#height slope repro
p.vec[14]<-	4.8894							#slope fruit production (diam)
p.vec[15]<- mu_diam							#mean recruit diam
p.vec[16]<-sd_diam							#sd recruit diam
p.vec[17]<-mu_height							#mean recruit height
p.vec[18]<-sd_height							#sd recruit height
p.vec[19]<-p_surv_seedlings					#probability of seedling establishment
p.vec[20]<-surv_mod_seedlings$coeff[1]		#intercept of seedling survival mod 
p.vec[21]<-surv_mod_seedlings$coeff[2]		#diam slope survival for seedlings
p.vec[22]<-surv_mod_seedlings$coeff[3]		#height slope survival for seedlings
p.vec[23]<-growth_mod_seedlings$coeff[1,1]	#intercept growth diameter for seedlings
p.vec[24]<-growth_mod_seedlings$coeff[2,1]	#diam slope growth diam for seedlings
p.vec[25]<-growth_mod_seedlings$coeff[3,1]	#height slope growth diam for seedlings
p.vec[26]<-sigma_growth_diam_seedlings		#sd growth diam for seedlings
p.vec[27]<-growth_mod_seedlings$coeff[1,2]	#intercept growth height for seedlings
p.vec[28]<-growth_mod_seedlings$coeff[2,2]	#diam slope growth height for seedlings
p.vec[29]<-growth_mod_seedlings$coeff[3,2]	#height slope growth height for seedlings
p.vec[30]<-sigma_growth_height_seedlings		#sd growth height for seedlings
p.vec[31]<-sdlng_grad_mod$coeff[1]			#intercept graduation 
p.vec[32]<-sdlng_grad_mod$coeff[2]			#diam slope graduation
p.vec[33]<-sdlng_grad_mod$coeff[3]			#height slope graduation
p.vec[34]<-mu_grad_diam						#mean graduate diam
p.vec[35]<-sd_grad_diam						#sd graduate diam
p.vec[36]<-mu_grad_height					#mean graduate height
p.vec[37]<-sd_grad_height					#sd graduate height
p.vec[38]<-0.005								#Probability fruit removed from tree
p.vec[39]<-0.002								#Probability seed survives journey

#p.vecs are defined externally in script 1demography_models.R
load('p.vec_overall.RData')
load('p.vec_E.RData')
load('p.vec_H.RData')
load('p.vec_W.RData')
#============================================================================# 
#  Define the kernels and iteration matrix:
#Domain 1: Seedling Domain
#	for diam=diameter in range [0, 1.6]
#	for height in range [0, 16]
#Domain 2: Larger Domain
#	 for diam=diameter in range [1.6,700]
#and for height=height in range [16, 800]
#============================================================================# 




#Survival-Growth of Seedlings (in D1)
pyx1=function(diamp,heightp,diam, height, params) { (1-((exp(params[31]+params[32]*diam + params[33]*height)/(1+exp(params[31]+params[32]*diam + params[33]*height))
)))*(exp(params[20]+params[21]*diam+params[22]*height)/(1+exp(params[20]+params[21]*diam+params[22]*height)))*dtnorm(diamp,mean=params[23] + params[24]*diam + params[25]*height,sd=params[26], lower=0, upper=1.6)*dtnorm(heightp,mean=params[27]+params[28]*diam + params[29]*height,sd=params[30], lower=0, upper=16)}


#Survival-Growth of Larger Plants (in D2)

pyx2=function(diamp,heightp,diam, height, params) { 0.997*(exp(params[1]+params[2]*diam+params[3]*height)/(1+exp(params[1]+params[2]*diam+params[3]*height)))*dtnorm(diamp,mean=params[4] + params[5]*diam + params[6]*height,sd=params[7], lower=1.6, upper=700)*dtnorm(heightp,mean=params[8]+params[9]*diam + params[10]*height,sd=params[11], lower=16, upper=800)}


#Fecundity=P(fruiting)*# of Fruits Produced*P(survival of seeds)*Distribution of Seedling Diameters*Distribution of Seedling Heights
fyx=function(diamp,heightp,diam,height, params) {(exp(params[12]+params[13]*height)/(1+exp(params[12]+params[13]*height)))*(params[14]*diam*diam)*params[38]*params[39]*params[19]*dtnorm(diamp,mean=params[15],sd=params[16], lower=0, upper=1.6)*dtnorm(heightp,mean=params[17],sd=params[18], lower=0, upper=16)}

#Graduation = P(graduation)*Distribution of Graduates Diams*Distribution of Graduates Heights
gyx=function(diamp,heightp,diam,height, params) {((exp(params[20]+params[21]*diam+params[22]*height)/(1+exp(params[20]+params[21]*diam+params[22]*height))))*(exp(params[31]+params[32]*diam + params[33]*height)/(1+exp(params[31]+params[32]*diam + params[33]*height))
)*dtnorm(diamp,mean=params[34],sd=params[35], lower=1.6, upper=700)*dtnorm(heightp,mean=params[36],sd=params[37], lower=16, upper=800)}

#kyx=function(diamp,heightp,diam,height, params) {pyx(diamp,heightp, diam,height, params)+fyx(diamp,heightp,diam,height, params)}

# Compute meshpoints  

###Fix these to match domains 
h1=1.6/m1; 
y1=(h1/2)*((0:(m1-1))+(1:m1)); #for diameter in D1
h2=16/m2
y2=(h2/2)*((0:(m2-1))+(1:m2)); #for height in D1
h3=(700-1.6)/m3;
y3=(h3/2)*((0:(m3-1))+(1:m3))+1.6; #for diameter in D2
h4=(800-16)/m4
 y4=(h4/2)*((0:(m4-1))+(1:m4))+16; #for height in D2

# Compute the iteration matrix. With a bit of vectorizing it's not too slow,
# though you can probably do better if you need to. The shortcuts here have 
# been checked against the results from code that uses loops for everything.

#####Part One: Overall model 
###Construct D1 (Seedling Domain):
plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in A 
Plop=outer(1:m1,1:m2,plop); 

D1=matrix(0,m1*m2,m1*m2); 
Kvals_D1=array(0,c(m1,m2,m1,m2));  



for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals=pyx1(y1[k],y2[1:m2],y1[i],y2[j], p.vec_overall)
				D1[Plop[k,1:m2],Plop[i,j]]=kvals
				Kvals_D1[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}
D1=D1*h1*h2 #multiply D1 by widths
		
###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 

D2=matrix(0,m3*m4,m3*m4); 
Kvals_D2=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals=pyx2(y3[k],y4[1:m4],y3[i],y4[j], p.vec_overall)
				D2[Plop[k,1:m4],Plop[i,j]]=kvals
				Kvals_D2[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		
D2=D2*h3*h4 #Multiply D2 by widths

###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);
F=matrix(0,m1*m2,m3*m4); 
Kvals_F=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=fyx(y1[k], y2[1:m2], y3[i], y4[j], p.vec_overall)
			F[Plop1[k, 1:m2], Plop2[i,j]]=kvals
			Kvals_F[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}
F=F*h1*h2

###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);
G=matrix(0,m3*m4,m1*m2); 
Kvals_G=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=gyx(y3[k], y4[1:m4], y1[i], y2[j], p.vec_overall)
			G[Plop1[k, 1:m4], Plop2[i,j]]=kvals
			Kvals_G[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}
G=G*h3*h4

#Putting the matrix together:


yellow_white<-rbind(D1, G)
white_green<-rbind(F, D2)
new_A<-cbind(yellow_white, white_green)

A<-new_A

save(A, file="A_matrix_overall.RData")
#rm(list=ls(all=TRUE))
#load("A_matrix.RData")

#lambda<-eigen(A)$values[1]
#lambda






#============================================================================# 
#  Find lambda, w by iteration 
#
#  Note: the Matrix package is used to iterate more quickly. The check  
#  for convergence requires extracting matrix entries via the @x slot
#  of a Matrix object. Matrix is S4-style -- see ?Matrix. 
#============================================================================# 
A2=Matrix(A); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	nt1=A2%*%nt;
	qmax=sum(abs((nt1-lam*nt)@x));  
	lam=sum(nt1@x); 
	nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
	cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
stable.dist=nt
lam.stable=lam;

# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 



##Alternate calculation of rv:

A3=t(A2); nt=Matrix(1,m1*m2+m3*m4,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	nt1=A3%*%nt;
	qmax=sum(abs((nt1-lam*nt)@x));  
	lam=sum(nt1@x); 
	nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
	cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2+m3*m4,1); 
#stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
#stable.dist=nt
lam.stable.t=lam; 


# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 






#============================================================================# 
#  Find reproductive value function by iteration.  
#  Note that we use left-multiplication to do the transpose iteration,
#  which is equivalent to using the transpose kernel. 
#============================================================================# 
vt=Matrix(1,1,m1*m2+m3*m4); vt1=vt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	vt1=vt%*%A2;
	qmax=sum(abs((vt1-lam*vt)@x));  
	lam=sum(vt1@x); 
	vt@x=(vt1@x)/lam;   
	cat(lam,qmax,"\n");
} 
v=t(matrix(vt@x,1,m1*m2+m3*m4)); 
lam.stable.t=lam; 



#####################################
####
#### ELASTICITY AND SENSITIVITY
####
#####################################




# Compute elasticity matrix 
#Ellner: 
#repro.val=matrix(v, m1, m2);
#stable.state=matrix(stable.dist, m1, m2);
#v.dot.w=sum(stable.state_smalls*repro.val_smalls)
#sens=outer(repro.val_smalls,stable.state_smalls)/v.dot.w
#elas=sens*Kvals/lam.stable;
# # Compute matrix of total(=integrated) elasticities for all transitions (x_i,q_j) -> anywhere 
#total.elas=h1*h2*apply(elas,c(3,4),sum);
# Checks
#cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
#cat("Integrated elasticity=",sum(h1*h2*h1*h2*elas_D1)," Should =1","\n"); 


####Method: Compute elasticity matrix first then separate:


v.dot.w<-sum(t(v)%*%stable.dist)
norm_v<-v/v.dot.w
check<-t(norm_v)%*%stable.dist #should be 1

rv<-norm_v

sens<-norm_v%*%t(stable.dist)
elas<-sens*A/lam.stable

plot(colSums(elas[1:m1*m2,]))

#Break the elasticity matrix back into its component parts: 
elas_D1<-elas[1:(m1*m2), 1:(m1*m2)]
elas_G<-elas[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
elas_F<-elas[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
elas_D2<-elas[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]

sens_D1<-sens[1:(m1*m2), 1:(m1*m2)]
sens_G<-sens[((m1*m2)+1):(m1*m2+m3*m4), 1:(m1*m2)]
sens_F<-sens[1:(m1*m2), ((m1*m2)+1):(m1*m2 + m3*m4)]
sens_D2<-sens[((m1*m2)+1):(m1*m2+m3*m4), ((m1*m2)+1):(m1*m2 + m3*m4)]





########
###Now go back to unpack the elasticity matrices into components:
########

plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in array
Plop=outer(1:m1,1:m2,plop); 

Kvals_elas_D1=array(0,c(m1,m2,m1,m2));  

for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals= elas_D1[Plop[k,1:m2],Plop[i,j]]
				Kvals_elas_D1[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}


###Construct D2 (Large Domain):
plop=function(i,j) {(j-1)*m3+i} # for putting values in proper place in A 
Plop=outer(1:m3,1:m4,plop); 


Kvals_elas_D2=array(0,c(m3,m4,m3,m4));  



for(i in 1:m3){
	for(j in 1:m4){
		for(k in 1:m3){
				kvals= elas_D2[Plop[k,1:m4],Plop[i,j]]
				Kvals_elas_D2[k,1:m4,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}		


###Construct F (Fecundity):
plop1=function(i, j) {(j-1)*m1 + i}
plop2=function(i, j) {(j-1)*m3 + i}
Plop1=outer(1:m1,1:m2,plop1); 
Plop2=outer(1:m3, 1:m4, plop2);

Kvals_elas_F=array(0, c(m1, m2, m3, m4))

for(i in 1:m3) {
	for (j in 1:m4) {
		for (k in 1:m1) {
			kvals=elas_F[Plop1[k, 1:m2], Plop2[i,j]]
			Kvals_elas_F[k, 1:m2, i, j]=kvals
		}}
		cat(i, "\n");
}


###Construct G (Graduation):
plop1=function(i, j) {(j-1)*m3 + i}
plop2=function(i, j) {(j-1)*m1 + i}
Plop1=outer(1:m3,1:m4,plop1); 
Plop2=outer(1:m1, 1:m2, plop2);

Kvals_elas_G=array(0, c(m3, m4, m1, m2))

for(i in 1:m1) {
	for (j in 1:m2) {
		for (k in 1:m3) {
			kvals=elas_G[Plop1[k, 1:m4], Plop2[i,j]]
			Kvals_elas_G[k, 1:m4, i, j]=kvals
		}}
		cat(i, "\n");
}


##############
total.elas_D1<-apply(Kvals_elas_D1, c(3,4), sum);
total.elas_D2<-apply(Kvals_elas_D2, c(3,4), sum);
total.elas_F<-apply(Kvals_elas_F, c(3,4), sum);
total.elas_G<-apply(Kvals_elas_G, c(3,4), sum);

#Separate:
repro.val_smalls=matrix(v[1:(m1*m2)], m1,m2);
stable.state_smalls=matrix(stable.dist[1:(m1*m2)], m1, m2);

repro.val_larges=matrix(v[((m1*m2)+1):((m1*m2)+(m3*m4))], m3, m4);
stable.state_larges=matrix(stable.dist[((m1*m2) +1):((m1*m2) + (m3*m4))], m3, m4);


#============================================================================# 
#  Calculate stable state distribution and state-dependent total elasticity 
#============================================================================# 

ss_diam1<-apply(stable.state_smalls, 1, sum)
ss_diam2<-apply(stable.state_larges, 1, sum)
ss_diam<-c(ss_diam1, ss_diam2)
y_diam<-c(y1, y3)



ss_height1<-apply(stable.state_smalls, 2, sum)
ss_height2<-apply(stable.state_larges, 2, sum)
ss_height<-c(ss_height1, ss_height2)
y_height<-c(y2, y4)




rv_diam1<-apply(repro.val_smalls, 1, sum)
rv_diam2<-apply(repro.val_larges, 1, sum)
rv_diam<-c(rv_diam1, rv_diam2)
y_diam<-c(y1, y3)





rv_height1<-apply(repro.val_smalls, 2, sum)
rv_height2<-apply(repro.val_larges, 2, sum)
rv_height<-c(rv_height1, rv_height2)
y_height<-c(y2, y4)



total.elas_D1_overall<-total.elas_D1
total.elas_D2_overall<-total.elas_D2
total.elas_G_overall<-total.elas_G
total.elas_F_overall<-total.elas_F

save(total.elas_D1_all, file="total.elas_D1_all.RData")
save(total.elas_D2_all, file="total.elas_D2_all.RData")
save(total.elas_G_all, file="total.elas_G_all.RData")
save(total.elas_F_all, file="total.elas_F_all.RData")



setwd("/Users/ke2/Dropbox/matrix/outputs/sites/")


##Redo these with comparisons between overall and Eastern 
elas_D1_ref_E<-total.elas_D1_all-total.elas_D1_E
elas_D2_ref_E<-total.elas_D2_all-total.elas_D2_E
elas_G_ref_E<-total.elas_G_all-total.elas_G_E
elas_F_ref_E<-total.elas_F_all-total.elas_F_E

elas_D1_ref_H<-total.elas_D1_all-total.elas_D1_H
elas_D2_ref_H<-total.elas_D2_all-total.elas_D2_H
elas_G_ref_H<-total.elas_G_all-total.elas_G_H
elas_F_ref_H<-total.elas_F_all-total.elas_F_H

elas_D1_ref_W<-total.elas_D1_all-total.elas_D1_W
elas_D2_ref_W<-total.elas_D2_all-total.elas_D2_W
elas_G_ref_W<-total.elas_G_all-total.elas_G_W
elas_F_ref_W<-total.elas_F_all-total.elas_F_W


jpeg('marginal_elasticity_height_ref.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(colSums(elas_D1_ref_W), xlab="Height (cm)", ylim=c(-.01, .01), ylab="Change in Elasticity", main="D1",  type='l', col=col_W_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_D1_ref_H), col=col_H_dark, type='l', lwd=2)
lines(colSums(elas_D1_ref_E), col=col_E_dark, type='l', lwd=2)


plot(colSums(elas_F_ref_H), ylim=c(-.01, 0.01),  xlab="Height (cm)", ylab= "Change in Elasticity", main="F", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_F_ref_E), type='l', col=col_E_dark, lwd=2)
lines(colSums(elas_F_ref_W), type='l', col=col_W_dark, lwd=2)

plot(colSums(elas_G_ref_H), ylim=c(-0.01, 0.01), xlab="Height (cm)", ylab="Change in Elasticity", main="G", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_G_ref_E), type='l', col=col_E_dark, lwd=2)
lines(colSums(elas_G_ref_W), type='l', col=col_W_dark, lwd=2)

plot(colSums(elas_D2_ref_H), ylim=c(-0.01, 0.01), xlab="Height (cm)", ylab="Change in Elasticity", main="D2", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(colSums(elas_D2_ref_E), type='l', col=col_E_dark, lwd=2)
lines(colSums(elas_D2_ref_W), type='l', col=col_W_dark, lwd=2)

dev.off()

###
jpeg('marginal_elasticity_diam_ref.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(rowSums(elas_D1_ref_W), xlab="Diameter (mm) ", ylim=c(-.01, .01), ylab="Change in Elasticity", main="D1",  type='l', col=col_W_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_D1_ref_H), col=col_H_dark, type='l', lwd=2)
lines(rowSums(elas_D1_ref_E), col=col_E_dark, type='l', lwd=2)


plot(rowSums(elas_F_ref_H), ylim=c(-.01, 0.01),  xlab="Diameter (mm)", ylab= "Change in Elasticity", main="F", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_F_ref_E), type='l', col=col_E_dark, lwd=2)
lines(rowSums(elas_F_ref_W), type='l', col=col_W_dark, lwd=2)

plot(rowSums(elas_G_ref_H), ylim=c(-0.01, 0.01), xlab="Diameter (mm)", ylab="Change in Elasticity", main="G", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_G_ref_E), type='l', col=col_E_dark, lwd=2)
lines(rowSums(elas_G_ref_W), type='l', col=col_W_dark, lwd=2)

plot(rowSums(elas_D2_ref_H), ylim=c(-0.01, 0.01), xlab="Diameter (mm)", ylab="Change in Elasticity", main="D2", type='l', col=col_H_dark, cex=2, cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(rowSums(elas_D2_ref_E), type='l', col=col_E_dark, lwd=2)
lines(rowSums(elas_D2_ref_W), type='l', col=col_W_dark, lwd=2)

dev.off()



###


jpeg('marginal_elasticity_diam_all.jpeg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=36)
plot(rowSums(total.elas_D1_all), ylim=c(0, 0.07),xlab="Size", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
plot(rowSums(total.elas_F_all),  ylim=c(0, 0.07), xlab="Size", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
plot(rowSums(total.elas_G_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
plot(rowSums(total.elas_D2_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
dev.off()

jpeg('marginal_elasticity_diam_height_all.jpeg')
par(mfrow=c(2,2), mar=c(5.1,6.1,4.1,2.1))
par(ps=36)
plot(rowSums(total.elas_D1_all), ylim=c(0, 0.07),xlab="Size", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_D1_all), type='l', lwd=2, col="red")
plot(rowSums(total.elas_F_all),  ylim=c(0, 0.07), xlab="Size", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_F_all), type='l', lwd=2, col="red")
plot(rowSums(total.elas_G_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_G_all), type='l', lwd=2, col="red")
plot(rowSums(total.elas_D2_all), ylim=c(0, 0.07),  xlab="Size", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=.5, cex.lab=1, xaxt='n')
lines(colSums(total.elas_D2_all), type='l', lwd=2, col="red")
dev.off()







###

jpeg('marginal_elasticity_diam_lineage.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(rowSums(total.elas_D1_E), ylim=c(0, 0.07),xlab="Diameter (mm)", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_D1_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_D1_H), type='l', lwd=2, col=col_H_dark)

plot(rowSums(total.elas_F_E),  ylim=c(0, 0.07), xlab="Diameter (mm)", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_F_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_F_H), type='l', lwd=2, col=col_H_dark)

plot(rowSums(total.elas_G_E), ylim=c(0, 0.07),  xlab="Diameter (mm)", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_G_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_G_H), type='l', lwd=2, col=col_H_dark)

plot(rowSums(total.elas_D2_E), ylim=c(0, 0.07),  xlab="Diameter (mm)", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(rowSums(total.elas_D2_W), type='l', lwd=2, col=col_W_dark)
lines(rowSums(total.elas_D2_H), type='l', lwd=2, col=col_H_dark)
dev.off()

###

jpeg('marginal_elasticity_height_lineage.jpeg')
par(mfrow=c(2,2), mar=c(5.1,5.1,4.1,2.1))
plot(colSums(total.elas_D1_E), ylim=c(0, 0.07),xlab="Height (cm)", ylab="Elasticity", main="D1",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_D1_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_D1_H), type='l', lwd=2, col=col_H_dark)

plot(colSums(total.elas_F_E),  ylim=c(0, 0.07), xlab="Height (cm)", ylab= "Elasticity", main="F", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_F_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_F_H), type='l', lwd=2, col=col_H_dark)

plot(colSums(total.elas_G_E), ylim=c(0, 0.07),  xlab="Height (cm)", ylab="Elasticity", main="G", type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_G_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_G_H), type='l', lwd=2, col=col_H_dark)

plot(colSums(total.elas_D2_E), ylim=c(0, 0.07),  xlab="Height (cm)", ylab="Elasticity", main="D2",  type='l', lwd=2, cex.axis=1.5, cex.lab=1.5, col=col_E_dark)
lines(colSums(total.elas_D2_W), type='l', lwd=2, col=col_W_dark)
lines(colSums(total.elas_D2_H), type='l', lwd=2, col=col_H_dark)
dev.off()

################################################################
#####
#####	Investigating the effects of varying the seed parameters
#####
#################################################################

setwd("/Users/kelley/Dropbox/March 2016 Documents/Documents/Grad/dissertation/Pratt_demographic_data") 

library(calibrate)

data3<-read.csv("seed_params.csv", header=T)


png(file="varying_tau.png", width=10, height=10, units="in", res=300)
par(ps=24, mar=c(5.1, 5.1, 4.1, 2.1))
plot(log10(data3$tau_postdispersal[1:4]), data3$lambda[1:4], type='b', ylab=expression(lambda), xlab=expression(log(tau[1])))
dev.off()

textxy(data3$param, data3$lambda, paste("(", data3$param, ", ", data3$lambda, ")"), cex=.7)
dev.off()

lambdas<-c(1.045574, 1.079985, 1.086125)

plot(lambdas, col=c(E_col_dark, H_col_dark, W_col_dark), 
