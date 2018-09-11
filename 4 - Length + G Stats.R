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

#data directory
#setwd("")
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
  # subs  return counts of subunits (T) or just total (F)
  #
  #return:
  # df  data.frame with columns and values needed for downstream processing
  
  #for this program length, series, and reads are needed
  
  df <- df[,c(1:8, grep('Length', colnames(df)), grep('Trim_Len', colnames(df)), grep('NT_coupure', colnames(df)))] #delimit to needed columns
  df$series <- 'No Realignment'
  df$series[df$Length!=df$Trim_Len] <- 'Realigned'
  
  #gets global read
  df$reads <- rowSums(df[,c(1:8)])
  
  if(!subs){
    df <- df[,-c(1:8)]
  }
  
  df2 <- df[df$Length!=df$Trim_Len,] #create second df for 'Trimmed'
  df2$series <- 'Trimmed'
  
  #decrease column number
  df <- df[,-grep('Trim_Len', colnames(df))]
  df2 <- df2[,-grep('Length', colnames(df2))]
  
  #set colname as 'values'
  colnames(df)[grep('Length', colnames(df))] <- 'values'
  colnames(df2)[grep('Trim_Len', colnames(df2))] <- 'values'
  
  #factor in G+1
  df$values[grepl('G',substr(df$NT_coupure,1,1))] <- df$values[grepl('G',substr(df$NT_coupure,1,1))] +1
  df2$values[grepl('G',substr(df2$NT_coupure,1,1))] <- df2$values[grepl('G',substr(df2$NT_coupure,1,1))] +1
  
  #remove NT_coupure
  df <- df[,-grep('NT_coupure', colnames(df))]
  df2 <- df2[,-grep('NT_coupure', colnames(df2))]
  
  #combine data frames without rbind (faster)
  df[(nrow(df)+1):(nrow(df)+nrow(df2)),] <- df2
  rm(df2) #remove df2
  
  
  df$id <- 1:nrow(df) #adds id for monte carlo randomization
  
  if(!subs){
    df <- df[,c(3,1,2,4)] # orders columns
  } else {
    df <- df[,c(1:8, 11, 9, 10, 12)]
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

#you will need to create a stats folder in the length folder
setwd("./Length/stats")
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

write.csv(file='Hong Kong Lengths.csv', hks, row.names = F)
write.csv(file='Puerto Rico 8 Lengths.csv', pr8s, row.names = F)
write.csv(file='WSN Lengths.csv', wsns, row.names = F)
write.csv(file='Brisbane Lengths.csv', bris, row.names = F)
write.csv(file='All Lengths.csv', gs, row.names = F)