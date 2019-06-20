
CV<-function(x) {
  return(sd(x)/mean(x))
  
}
#D2
D2_BC <- readRDS("./BC/D2_BC.rds")
D2_CC <- readRDS("./CC/D2_CC.rds")
D2_C  <- readRDS("./C/D2_C.rds")
D2_FP <- readRDS("./FP/D2_FP.rds")
D2_PG <- readRDS("./PG/D2_PG.rds")
D2_WT <- readRDS("./WT/D2_WT.rds")


cv_D2<-matrix(nrow=(m3*m4), ncol=(m3*m4))
for(i in 1:(m3*m4)) {
  for(j in 1:(m3*m4)) {
    cv_D2[i,j]<-CV(c(D2_BC[i,j], D2_CC[i,j], D2_C[i,j], D2_FP,
                     D2_PG, D2_WT))
    if(is.na(cv_D2[i,j])) {cv_D2[i,j]<-0}
  }
  cat("i=", i, "\n")
}
save(cv_D2, file="cv_D2.RData")