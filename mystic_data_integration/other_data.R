pfs <- function(start_date,pd_date,fu_date,eot,cens_os){
  x<-start_date
  y<-pd_date
  z<-fu_date
  a<-eot
  b<-cens_os
  pfs1 <- as.numeric(y-x)/30.41
  pfs2 <- as.numeric(z-x)/30.41
  pfs3 <- as.numeric(a-x)/30.41
  
  pfsf <- ifelse(!is.na(pfs1), pfs1,
                 ifelse(!is.na(pfs2), pfs2,
                        ifelse(!is.na(pfs3), pfs3, "error")))
  pfsc <- ifelse(!is.na(pfs1), 1,
                 ifelse(!is.na(pfs2), b,
                        ifelse(!is.na(pfs3), 0, "error")))
  ts <- ifelse(pfsf=="error"&is.na(x), "cons", "re")
  ts <- ifelse(ts=="re"&pfsf!="error", "cons", ts)
  
  if(sum(ts=="re")!=0){
    print("There's an error. One date missing for one or more patients! Check dates")
  }else{
    print("Troubleshooting successful. All suitable patients have a value for the PFS. Be aware of other errors")
  }
  if(length(unique(is.negative(pfsf)))==2){
    print("There's an error! Negative survival times, check dates")
  }else{
    print ("No other errors found!")
  }
  pfs_time <<- as.numeric(pfsf)
  pfs_cens <<- as.numeric(pfsc)
}

ntable <- function(x){
  lvl <- levels(as.factor(x))
  nass <- ifelse(sum(is.na(x))!=0,"NA",0)
  #lvl <- c(lvl,nas)
  sum(is.na(x))
  a <- NULL
  b <- NULL
  nas <- NULL
  nasp <- NULL
  
  for(i in 1:length(lvl)){
    a[i] <- sum(x==lvl[i], na.rm=T)
    nas <- sum(is.na(x))
    b[i] <- sum(x==lvl[i], na.rm=T)/length(x)*100
    nasp <- sum(is.na(x))/length(x)*100
    
  }
  if(nass!=0){
    dafra <- data.frame("Levels" = c(lvl,nass),
                        "values" = c(a,nas),
                        "percentages" = c(b,nasp))
  }else{
    dafra <- data.frame("Levels" = c(lvl),
                        "values" = c(a),
                        "percentages" = c(b))
  }
  dafra
}


ntable2 <- function(character, groups){
  vectt <- NULL
  combsf <- NULL
  props <- NULL
  props2 <- NULL
  ns <- NULL
  y <- character
  x <- groups
  xlev <- levels(as.factor(x))
  ylev <- levels(as.factor(y))
  combs <- combn(c(xlev,ylev),2, simplify=F)
  for(i in 1:length(combs)){
    if((combs[i][[1]][[1]] %in% xlev & combs[i][[1]][[2]] %in% xlev) | (combs[i][[1]][[1]] %in% ylev & combs[i][[1]][[2]] %in% ylev) ){
    }else{
      combsf[i] <- combs[i]
    }
  }
  combsf <- combsf[!sapply(combsf,is.null)]
  t1 <- table(x,y)
  trow <- nrow(t1)
  tcol <- ncol(t1)
  combt <- trow*tcol
  combs2 <- paste(sapply(combsf, paste, collapse='"'))
  combs3 <- stringr::str_replace(combs2, '\"', ' ')
  for(i in 1:length(xlev)){
    ns[i] <- sum(x==xlev[i], na.rm=T)
  }
  for (i in 1:combt){
    vectt[i] <- t(t1)[i]
  }
  n <- length(unlist(combsf))/2
  result <- unlist(combsf)[c(2*(1:n)-1,2*(n:1))]
  result <- result[1:n]
  for(i in 1:combt){
    props2[i] <- sum(x==result[i], na.rm=T)
  }
  props2 <- (vectt/props2)*100
  vect <- data.frame("Levels"=combs3, "Values"=vectt,"Proportions 1st Level"=props2)
  vect
}

print_progress <- function(i, n, initial_iter){
  when_is_10perc <-( n-initial_iter)*0.1
  perc_points <- c(round(when_is_10perc)+initial_iter, 2*round(when_is_10perc)+initial_iter,3*round(when_is_10perc)+initial_iter,
                   4*round(when_is_10perc)+initial_iter,5*round(when_is_10perc)+initial_iter,6*round(when_is_10perc)+initial_iter,
                   7*round(when_is_10perc)+initial_iter,8*round(when_is_10perc)+initial_iter,9*round(when_is_10perc)+initial_iter)
  if(i==initial_iter){
    ini_time <<- Sys.time()
  }
  if(i==n){
    final_time <<- Sys.time()
    print(paste0("Complete! ",round(final_time-ini_time,2)))
  }
  if(i %in% perc_points){
    pct <- which(perc_points==i)*10
    meanw_time <<- Sys.time()
    if(i==perc_points[1]){
      print(paste0(pct, "%", " Elapsed: ", round(meanw_time-ini_time,2)))
      print(paste0("Approximate finish time ", format(strptime(Sys.time()+9*(meanw_time-ini_time), "%Y-%m-%d %H:%M:%S"), '%H:%M')))
    }else{
      print(paste0(format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), '%H:%M')," " ,pct, "%", " Elapsed: ", round(meanw_time-ini_time, 2)))
      print(paste0("Expect ", round((100-pct)*(meanw_time-ini_time)/pct,2)," to finish"))
    }
  }
}

len <- function(x){
  length(x)
}
