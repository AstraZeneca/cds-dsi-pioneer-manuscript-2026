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

len <- function(x){
  length(x)
}
