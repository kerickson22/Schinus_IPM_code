# Transient dynamics analysis using the IPM constructed from 
# demography_sites.R.  
# Kelley D. Erickson 

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
setwd("/Users/curculion/Dropbox/matrix/outputs/m1=10/m3=100")

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

eigenvalue_all <- eigen(A_all)
save(eigenvalue_all, file = "eigenvalue_all.RData")

rm(A_all)
rm(eigenvalue_all)

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

eigenvalue_E <- eigen(A_E)
save(eigenvalue_E, file = "eigenvalue_E.RData")

rm(A_E)
rm(eigenvalue_E)

