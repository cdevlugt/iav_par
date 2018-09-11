require(compiler)
enableJIT(3)

require(doParallel)

subunitDecompress <- function(df, sub=1, drops=c(1:8), drop=T, name="Subunit"){
  # decompresses the data by replicating a column based on read count; can also remove columns
  #
  # Args:
  #   df: data frame containing data
  #   sub:  column containing data to be replicated
  #   drops:  columns to be dropped
  #   drop: should the other subunits be dropped
  #
  # Returns:
  #   df: decompress data frame for a given subunit
  if(drop){
    drops <- drops[drops != sub]#retains subunit column if in drops range
    if(length(drops)!=0){
      df <- df[,-drops]
      sub <- sub - sum(drops < sub)#one more opperation versus replicating than dropping, but stores less stuff in ram during replication
    }
  }
  df <- df[(rep(row.names(df), df[,sub])),]
  df[,sub] <- colnames(df)[sub]
  colnames(df)[sub] <- name
  return (df)
}
subunitSeperate <- function(df, subs=c(1:8), decompress=F, name="Subunit"){
  #Takes a data frame with read counts in columns and then clones the rows seperating the subunits so that they are on seperate rows
  #
  #Args:
  # df: data frame with subunit read counts
  # subs: column numbers of subunits
  # decompress: flag indicating if the subunits should be decompressed
  # name: name for collumn if decompress is T
  #
  #Return:
  # dfFinal: data frame with columns cloned based on subunit
  
  for(i in 1:length(subs)){
    dfTemp <- df[df[,subs[i]]>0,]
    if(decompress){
      dfTemp <- subunitDecompress(dfTemp, sub=subs[i], drops=subs, drop=T, name=name)
    } else {
      dfTemp[,subs[-i]] <- 0
    }
    if(i==1){
      dfFinal <- dfTemp
    } else {
      dfFinal <- rbind(dfFinal, dfTemp)
    }
    rm(dfTemp)
  }
  return(dfFinal)
}
seqComp <- function(df, chars=c('A', 'C', 'G', 'T'), Sequences=10, otherCols=NA, ignore.case=F, Occur=T, Freq=T, All=F, purines=F, strong=F){
  #calculates the relative or absolute portion of a string that is made up of a substring; works to determine the portion of a string made up of
  # RNA or DNA bases/motifs
  #
  #Args:
  # df: dataframe or vector containg the strings
  # chars:  substrings to search for
  # Sequences:  columns containing strings to analyze
  # otherCols:  allows for retention of other data frame columns; you may also cbind results and ignore this feature
  # ignore.case:  ignore case when searching for substrings
  # Occur:  count occurances
  # Freq: count frequency
  # All:  return all data rather than just occurances and/or frequency
  # purines:  determine A and G content of string
  # strong: determine G and C content of string
  #
  #return:
  # cbind(df, dfR): data frame containing sequences and other columns
  # dfR:  dataframe containing sequence composition as specified in Occur and Freq
  
  
  if(!Occur && !Freq){
    stop("Specify either Occur of Freq as true; Occur returns occurances, Freq returns frequency; both may be true")
  }
  
  cols <- Sequences
  Sequences <- as.vector(1:length(Sequences))
  
  if(!is.na(otherCols[1])){
    cols <- c(otherCols, cols)
    otherCols <- 1:length(otherCols)
    Sequences <- as.vector((length(cols)-length(Sequences)+1):length(cols))
  }
  df <- data.frame(df[,cols], stringsAsFactors = F)
  cols <- ncol(df)
  if(Occur){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(chars[j])
        colnames(df)[ncol(df)] <- paste("Abs",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    dfR <- df[,c((cols+1):ncol(df))]
    df <- df[,1:cols]
  }
  if(Freq){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    
    if(Occur){
      dfR <- cbind(dfR, df[,c((cols+1):ncol(df))])
    } else {
      dfR <- df[,c((cols+1):ncol(df))]
    }
    dfR[is.na(dfR)] <-0
    df <- df[,1:cols]
  }
  
  if(All){
    return(cbind(df, dfR))
  } else {
    return(dfR)
  }
}

#setwd("") #if not default directory
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

dfPreProcess <- function(df, subs=F, sep=F, strain=''){
  #prepares the data for downstream analysis
  # reads | values | series | id | strain
  #
  #args:
  # df  data.frame
  #
  #return:
  # df  data.frame with columns and values needed for downstream processing
  
  #for this program length, series, and reads are needed
  
  df <- df[,c(1:8, grep('Sequence', colnames(df))[1], grep('Trim_Sequence', colnames(df))[1])] #delimit to needed columns
  df$series <- 'No Realignment'
  df$series[df$Sequence!=df$Trim_Sequence] <- 'Realigned'
  
  #gets global read
  df$reads <- rowSums(df[,c(1:8)])
  
  #get G+C 
  df <- cbind(df, seqComp(df, chars = 'G', Sequences =c(9,10), otherCols = T, Occur = F, Freq=T, All = F, purines=F, strong = T)[,c(2,4)])
  
  if(!subs){
    df <- df[,-c(1:8)]
  }
  
  #split data
  df2 <- df[df$Sequence!=df$Trim_Sequence,] #create second df for 'Trimmed'
  df2$series <- 'Trimmed'
  
  #decrease column number
  df <- df[,-c(grep('Sequence', colnames(df))[1], grep('Trim_Sequence', colnames(df))[1], grep('GC_Trim_Sequence', colnames(df))[1])]
  df2 <- df2[,-c(grep('Sequence', colnames(df2))[1], grep('Trim_Sequence', colnames(df2))[1], grep('GC_Sequence', colnames(df2))[1])]
  
  #set colname as 'values'
  colnames(df)[grep('Sequence', colnames(df))] <- 'values'
  colnames(df2)[grep('Trim_Sequence', colnames(df2))] <- 'values'
  
  #combine data frames without rbind (faster)
  df[(nrow(df)+1):(nrow(df)+nrow(df2)),] <- df2
  rm(df2) #remove df2
  gc()
  
  df$id <- 1:nrow(df) #adds id for monte carlo randomization
  
  # reads | values | series | id
  if(!subs){
    df <- df[,c(2,3,1,4)] # orders columns
  } else {
    df <- df[,c(1:8, 10, 11, 9, 12)]
    if(sep){
      df <- df[,-9]
      df <- subunitSeperate(df,1:8,T)
    }
  }
  
  if(strain!=''){
    df$strain <- strain
  }
  
  return(df)
}

setwd("./Sequence Composition")
hk <- dfPreProcess(hk, T, T, 'Hong Kong')
pr8 <- dfPreProcess(pr8, T, T, 'Puerto Rico')
wsn <- dfPreProcess(wsn, T, T, 'WSN')
bri <- dfPreProcess(bri, T, T, 'Brisbane')

#construct table

statsTable <- function(df, by=1){
  #constructs table as
  #subunit|series|mean|sd|ci1|ci2|pval vs s1|pval vs s2|pval vs s3
  
  series <- unique(df$series) #series names
  targets <- unique(df[,by]) #get targets (subunits of strains)
  targetName <- colnames(df)[by]
  
  st <- data.frame(target='All', series=series, mean=0, sd=0, ci1=0, ci2=0, p1=0, p2=0, p3=0, stringsAsFactors = F)
  
  statVals <- function(df, series){
    #faster to do each opperation individually rather than as a group
    
    cores <- 3
    cl <- makeCluster(cores) 
    registerDoParallel(cl)
    
    comb <- function(x, ...) {
      #combines things, stolen from: 
      # https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    
    stVals <- foreach(i=1:3, .combine='comb', .multicombine=TRUE,
                      .init=list(list(), list(), list(), list(), list())) %dopar% {
                        list(mean(df$values[df$series==series[i]]), 
                             sd(df$values[df$series==series[i]]), 
                             t.test(df$values[df$series==series[i]], df$values[df$series==series[1]], var.equal = T)$p.value, 
                             t.test(df$values[df$series==series[i]], df$values[df$series==series[2]], var.equal = T)$p.value, 
                             t.test(df$values[df$series==series[i]], df$values[df$series==series[3]], var.equal = T)$p.value)
                      }
    
    stopCluster(cl)
    
    st1 <- data.frame(
      mean=unlist(stVals[[1]]), 
      sd=unlist(stVals[[2]]), 
      ci1=0, 
      ci2=0, 
      p1=unlist(stVals[[3]]), 
      p2=unlist(stVals[[4]]), 
      p3=unlist(stVals[[5]]),
      stringsAsFactors = F
    )
    
    st1$ci1 <- st1$mean - 2*st1$sd
    st1$ci2 <- st1$mean + 2*st1$sd
    
    return(st1)
  }
  
  st[,c(3:9)] <- statVals(df, series)
  
  for(i in 1:length(targets)){
    st[(nrow(st)+1):(nrow(st)+3),] <- data.frame(
      taget=targets[i],
      series=series,
      mean=0, 
      sd=0, 
      ci1=0, 
      ci2=0, 
      p1=0, 
      p2=0, 
      p3=0,
      stringsAsFactors = F)
    
    st[(3*i + 1):(3*i + 3),c(3:9)] <- statVals(df[df[,by]==targets[i],], series)
  }
  
  return(st)
}

hks <- statsTable(hk, 1)
pr8s <- statsTable(pr8, 1)
wsns <- statsTable(wsn, 1)
bris <- statsTable(bri, 1)

global <- rbind(hk, pr8, wsn, bri)
gs <- statsTable(global, 5)

write.csv(file='Hong Kong G+C.csv', hks, row.names = F)
write.csv(file='Puerto Rico 8 G+C.csv', pr8s, row.names = F)
write.csv(file='WSN G+C.csv', wsns, row.names = F)
write.csv(file='Brisbane G+C.csv', bris, row.names = F)
write.csv(file='All G+C.csv', gs, row.names = F)

#####Genome levels
subunitDecompress <- function(df, sub=1, drops=c(1:8), drop=T, name="Subunit"){
  # decompresses the data by replicating a column based on read count; can also remove columns
  #
  # Args:
  #   df: data frame containing data
  #   sub:  column containing data to be replicated
  #   drops:  columns to be dropped
  #   drop: should the other subunits be dropped
  #
  # Returns:
  #   df: decompress data frame for a given subunit
  if(drop){
    drops <- drops[drops != sub]#retains subunit column if in drops range
    if(length(drops)!=0){
      df <- df[,-drops]
      sub <- sub - sum(drops < sub)#one more opperation versus replicating than dropping, but stores less stuff in ram during replication
    }
  }
  df <- df[(rep(row.names(df), df[,sub])),]
  df[,sub] <- colnames(df)[sub]
  colnames(df)[sub] <- name
  return (df)
}
subunitSeperate <- function(df, subs=c(1:8), decompress=F, name="Subunit"){
  #Takes a data frame with read counts in columns and then clones the rows seperating the subunits so that they are on seperate rows
  #
  #Args:
  # df: data frame with subunit read counts
  # subs: column numbers of subunits
  # decompress: flag indicating if the subunits should be decompressed
  # name: name for collumn if decompress is T
  #
  #Return:
  # dfFinal: data frame with columns cloned based on subunit
  
  for(i in 1:length(subs)){
    dfTemp <- df[df[,subs[i]]>0,]
    if(decompress){
      dfTemp <- subunitDecompress(dfTemp, sub=subs[i], drops=subs, drop=T, name=name)
    } else {
      dfTemp[,subs[-i]] <- 0
    }
    if(i==1){
      dfFinal <- dfTemp
    } else {
      dfFinal <- rbind(dfFinal, dfTemp)
    }
    rm(dfTemp)
  }
  return(dfFinal)
}
seqComp <- function(df, chars=c('A', 'C', 'G', 'T'), Sequences=10, otherCols=NA, ignore.case=F, Occur=T, Freq=T, All=F, purines=F, strong=F){
  #calculates the relative or absolute portion of a string that is made up of a substring; works to determine the portion of a string made up of
  # RNA or DNA bases/motifs
  #
  #Args:
  # df: dataframe or vector containg the strings
  # chars:  substrings to search for
  # Sequences:  columns containing strings to analyze
  # otherCols:  allows for retention of other data frame columns; you may also cbind results and ignore this feature
  # ignore.case:  ignore case when searching for substrings
  # Occur:  count occurances
  # Freq: count frequency
  # All:  return all data rather than just occurances and/or frequency
  # purines:  determine A and G content of string
  # strong: determine G and C content of string
  #
  #return:
  # cbind(df, dfR): data frame containing sequences and other columns
  # dfR:  dataframe containing sequence composition as specified in Occur and Freq
  
  
  if(!Occur && !Freq){
    stop("Specify either Occur of Freq as true; Occur returns occurances, Freq returns frequency; both may be true")
  }
  
  cols <- Sequences
  Sequences <- as.vector(1:length(Sequences))
  
  if(!is.na(otherCols[1])){
    cols <- c(otherCols, cols)
    otherCols <- 1:length(otherCols)
    Sequences <- as.vector((length(cols)-length(Sequences)+1):length(cols))
  }
  df <- data.frame(df[,cols], stringsAsFactors = F)
  cols <- ncol(df)
  if(Occur){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(chars[j])
        colnames(df)[ncol(df)] <- paste("Abs",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    dfR <- df[,c((cols+1):ncol(df))]
    df <- df[,1:cols]
  }
  if(Freq){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    
    if(Occur){
      dfR <- cbind(dfR, df[,c((cols+1):ncol(df))])
    } else {
      dfR <- df[,c((cols+1):ncol(df))]
    }
    dfR[is.na(dfR)] <-0
    df <- df[,1:cols]
  }
  
  if(All){
    return(cbind(df, dfR))
  } else {
    return(dfR)
  }
}
dfPreProcess <- function(df, subs=c(1:8)){
  #calculate the mean and standard deviation of G+C content in the genome
  
  #extracts the mini genome and the G+C content
  df <- df[,c(subs, grep('mini_genome', colnames(df)))] #delimit data to subs and genome
  df <- seqComp(df, 'C', max(subs)+1, subs, F, F, T, T, F, T) #get genome composition
  df <- df[,c(subs, grep('Rel_GC_', colnames(df)))] #delimit data
  
  #decompress
  df <- subunitSeperate(df, subs, T)
  
  vals <- c(mean(df$Rel_GC_mini_genome), sd(df$Rel_GC_mini_genome))
  return(vals)
}

setwd("..")
hk <- dfPreProcess(read.csv("HK_All_Match.csv", stringsAsFactors = F))
pr8 <- dfPreProcess(read.csv("PR8_All_Match.csv", stringsAsFactors = F))
wsn <- dfPreProcess(read.csv("WSN_All_Match.csv", stringsAsFactors = F))
bri <- dfPreProcess(read.csv("BRI_All_Match.csv", stringsAsFactors = F))

df <- data.frame(pr8 = pr8, hk = hk, wsn = wsn, bri = bri, stringsAsFactors = F)
row.names(df) <- c('mean', 'standard deviation')
setwd("./Sequence Composition")
write.csv(df, 'genome.csv')
