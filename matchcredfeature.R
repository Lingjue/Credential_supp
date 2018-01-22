


matchcredfeature = function(dt1T1, dt1T2, ppm = 15, drt = 1) {
  
match_combined=data.table()
nomatch1_combined=data.table()
nomatch2_combined=data.table()
quipu_match=data.table()

quipu_m1 = numeric()
quipu_m2 = numeric()
rtmean_m1 = numeric()
rtmean_m2 = numeric()
basemz_m1 = numeric()
basemz_m2 = numeric()
quipu_nm1 = numeric()
quipu_nm2 = numeric()
quipu_dm1 = numeric()
quipu_dm2 = numeric()
match_idx2 =numeric()
fratio = numeric()

for (tmp.quipu1 in unique(dt1T1[,quipu])){

  idx1 <- which(dt1T1[,"quipu"]==tmp.quipu1)
  
  basemz1 = min(dt1T1[idx1,mz])
  minmz1 = basemz1-basemz1*ppm/1e6
  maxmz1 = basemz1+basemz1*ppm/1e6
  basert1 = mean(dt1T1[idx1,rt])
  minrt1 = basert1-drt
  maxrt1 = basert1+drt
  
  idx2 <- which(dt1T2[,"mz"]>=minmz1 & dt1T2[,"mz"]<=maxmz1 & dt1T2[,"rt"]>=minrt1 & dt1T2[,"rt"]<=maxrt1)
  
  if(length(idx2)<1){
    cat("No match w/ quipu#:",tmp.quipu1,"\n")
    nomatch1_combined = rbind(nomatch1_combined,dt1T1[idx1])
  }
  if(length(idx2)>1){
    tmp.quipu2 = unique(dt1T2[idx2,quipu])
    cat("Duplicate match between quipu#1:",tmp.quipu1,"and quipu#2:",tmp.quipu2, "\n")
    quipu_dm1 = c(quipu_dm1,tmp.quipu1)
    quipu_dm2 = c(quipu_dm2,tmp.quipu2)
    
  }
  if(length(idx2)==1){
    
    tmp.quipu2 <- dt1T2[idx2,quipu]
    basemz2 = dt1T2[idx2,mz]
    idx2 <- which(dt1T2[,"quipu"]==tmp.quipu2)  #extract all quipu in set 2
    merged = cbind.fill(dt1T1[idx1,],dt1T2[idx2,])
    ratio = dt1T1[idx1[1],ratio]/dt1T2[idx2[1],ratio]
    basert2 = mean(dt1T2[idx2,rt])
    fratio = c(fratio,ratio)
    quipu_m1 = c(quipu_m1,tmp.quipu1)
    quipu_m2 = c(quipu_m2,tmp.quipu2)
    basemz_m1 = c(basemz_m1,basemz1)
    basemz_m2 = c(basemz_m2,basemz2)
    rtmean_m1 = c(rtmean_m1,basert1)
    rtmean_m2 = c(rtmean_m2,basert2)
    match_combined = rbind(match_combined,merged)
    match_idx2 = c(match_idx2,idx2)
   }
}



quipu_match = data.table(cbind(quipu_m1,rtmean_m1,basemz_m1,quipu_m2,rtmean_m2,basemz_m2,fratio))
nomatch2_combined = dt1T2[!match_idx2]
nomatch1 = length(unique(dt1T1[,quipu]))-length(quipu_match[,quipu_m1])
nomatch2 = length(unique(dt1T2[,quipu]))-length(quipu_match[,quipu_m2])

cat(length(quipu_match[,quipu_m1]), "matches are found.\n")
cat(nomatch1, "unmatched quipu#1 and",nomatch2,"unmatched quipu#2 are found.\n")
cat(length(quipu_dm1),"quipu#1 and",length(quipu_dm2),"quipu#2 are duplicated match.\n")

list(quipu_match,match_combined,nomatch1_combined,nomatch2_combined)
}
