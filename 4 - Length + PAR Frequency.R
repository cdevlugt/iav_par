require(compiler)

enableJIT(3)

require(ggplot2)
require(cowplot)

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

#setwd("") #data directory
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

dfPreProcess <- function(df, strain='') {
  #drop rows unneeded downstream
  df <- df[,c(1:8, grep('Length', colnames(df)), grep('Trim_Len', colnames(df)), grep('NT_c', colnames(df)))]
  
  #label series
  df$series <- 'No Realignment'
  df$series[df$Length!=df$Trim_Len] <- 'Trimmed'
  
  #assign L+1 to gcomp series
  df$Trim_Len[startsWith(df$NT_coupure, 'G')] <- df$Trim_Len[startsWith(df$NT_coupure, 'G')]  +1
  
  #drop length, rename trim_length
  df <- df[,c(1:8,10, 12)]
  colnames(df)[9] <- 'Length'
  
  #name strain
  df$strain <- strain
  
  #decompress data
  df <- subunitSeperate(df, decompress = T)
  
  return(df)
}

hk <- dfPreProcess(hk, 'Hong Kong')
pr8 <- dfPreProcess(pr8, 'Puerto Rico')
wsn <- dfPreProcess(wsn, "WSN")
bri <- dfPreProcess(bri, 'Brisbane')


realignmentFrequency <- function(df=NA, col=2, limits=c(8:17), quant=F, freq=T, name='') {
  
  #spoof the graph and return legend if nothing is passed in df
  if(is.null(nrow(df))){
    re <- data.frame(series=c('No_Realignment', 'Trimmed'), x=1, y=1)
    
    reGraph <- ggplot() +
      geom_bar(data=re, aes(group=series, y=y, x=x, fill=series), stat='identity') + 
      theme_bw() +
      scale_fill_manual(name='Series', breaks=c('No_Realignment', 'Trimmed'), 
                        values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                        labels=c(No_Realignment='No Realignment', Trimmed='Realigned'))
    return(get_legend(reGraph))
  } else if(nrow(df)==0) {
    return(NULL)
  }
  
  if(quant!=F){
    limits <- quantile(df[,col], prob = seq(0, 1, length=(quant+1)))
    iter <- length(limits) - 1
  } else {
    iter <- length(limits)
  }
  
  #seperate data into temporary data frames
  for(i in 1:iter){
    assign(paste('df',i,sep=''), df[df[,col]==limits[i],])
  }
  
  #create dataframe for plot
  if(quant!=F){
    re <- data.frame(x1=limits[1:(length(limits)-1)], x2=limits[2:length(limits)], NoR=0, Trim=0, all=0, stringsAsFactors = F)
  } else {
    re <- data.frame(x1=limits[1:(length(limits))], x2=limits[1:length(limits)], NoR=0, Trim=0, all=0, stringsAsFactors = F)
  }
  
  #get number of rows in each temporary data frame for each of the 2 series
  for(i in 1:iter){
    re[i,3] <- nrow(get(paste('df',i,sep=''))[get(paste('df',i,sep=''))$series=='No Realignment',])
    re[i,4] <- nrow(get(paste('df',i,sep=''))[get(paste('df',i,sep=''))$series=='Trimmed',])
    rm(list=paste('df', i, sep='')) #removes the i-th temporary df used for sorting
  }
  
  re$all <- rowSums(re[,c(3,4)])
  if(quant!=F){
    xlabs <- paste(signif(re[,1],3), '% to ', signif(re[,2], 3), '%', sep='')
  } else {
    xlabs <- limits
  }
  re$PAR <- re$Trim/re$all * 100
  write.csv(re, paste(name, '.csv', sep=''), row.names = F) #write for fold change
  re$x <- 1:nrow(re)
  
  if(freq){
    re[,c(3,4)] <- re[,c(3,4)]/re$all *100
    
    #get global realignment frequency
    reFreq <- nrow(df[df$series=='Trimmed',])/nrow(df) * 100
  } else {
    #global realignment frequency as count
    reFreq <- data.frame(y=(re$all * nrow(df[df$series=='Trimmed',])/nrow(df)), x=0)
    reFreq$x <- 1:nrow(reFreq)
  }
  
  re <- data.frame(
    series=c(rep('No_Realignment', iter), rep('Trimmed', iter)), x=re$x, y=c(re$NoR, re$Trim)
  )
  
  reGraph <- ggplot() +
    geom_bar(data=re, aes(group=series, y=y, x=x, fill=series), stat='identity') + 
    theme_bw() +
    scale_x_continuous(breaks=1:(length(xlabs)), labels = xlabs) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab('Length (nt)') +
    scale_fill_manual(breaks=c('No_Realignment', 'Trimmed'), 
                      values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                      labels=c(No_Realignment='No Realignment', Trimmed='Trimmed'), guide=F) 
  
  if(freq){
    reGraph <- reGraph +
      geom_hline(yintercept = reFreq) +
      ylab('Frequency (%)')
  } else {
    reGraph <- reGraph +
      stat_summary(fun.y='mean', data=reFreq, colour="#000000", geom="errorbar", aes(x=x, y=y,ymax=..y.., ymin=..y..), width=1, linetype="solid", show.legend=F) +
      ylab('Count')
  }
  
  rm(df, col, limits, quant, re)
  
  return(reGraph)
}
getGraphs <- function(df, target='Subunit', seqLen=NA, col=2, limits=c(8:17), quant=F, freq=T){
  #generates and saves graphs from a df from dfPreProcess()
  #
  #Args
  # df: dataframe from dfPreProcess
  # col:  target column
  # limits: limits for groups
  # quant:  number of quantiles
  # target: column to create graphs around (ex., Strain, Subunit)
  # seqLen: vector of sequences to examine; if NA analyzes all lengths
  
  colnames(df)[grep(target, colnames(df))] <- 'target'
  subs <- unique(df$target)
  if(target=="Strain" || target=="strain" ){
    name <- 'Global'
  } else {
    name <- df$strain[1]
  }
  
  if(length(seqLen)==1 && is.na(seqLen)){
    seqLen <- c(min(df$Length):max(df$Length))
  }
  
  LoN <- realignmentFrequency() #graph legend
  
  print(target)
  
  for(i in 1:ceiling(length(subs)/4)){
    print(subs[4*(i-1)+1])
    titleText <- gsub("\\.","",paste(subs[4*(i-1)+1],target, sep=" "))
    title <- ggdraw() + draw_label(titleText, fontface='bold')
    graph1 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+1],], col, limits, quant, freq = freq, name=paste(name ,titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    
    if(length(subs)<(4*(i-1)+2)){
      graph2 <- NULL
    } else {
      print(subs[4*(i-1)+2])
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+2],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph2 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+2],], col, limits, quant, freq = freq, name=paste(name ,titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(length(subs)<(4*(i-1)+3)){
      graph3 <- NULL
    } else {
      print(subs[4*(i-1)+3])
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+3],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph3 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+3],], col, limits, quant, freq = freq, name=paste(name ,titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(length(subs)<(4*(i-1)+4)){
      graph4 <- NULL
    } else {
      print(subs[4*(i-1)+4])
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+4],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph4 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+4],], col, limits, quant, freq = freq, name=paste(name ,titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(i==1){
      graph <- plot_grid(graph1, graph2, graph3, graph4, LoN, ncol=5, nrow=1, rel_widths = c(5,5,5,5,1.6), labels = c(LETTERS[(4*(i-1)+1):(4*(i-1)+4)], NULL))
    } else {
      graph <- plot_grid(graph, plot_grid(graph1, graph2, graph3, graph4, NULL, ncol=5, nrow=1, rel_widths = c(5,5,5,5,1.6), labels = c(LETTERS[(4*(i-1)+1):(4*(i-1)+4)], NULL)), rel_widths=c(1,1), rel_heights=c(i-1, 1), ncol=1, nrow=2)
    }
    
    rm(graph1, graph2, graph3, graph4)
    gc()
    
    if(i==ceiling(length(subs)/4)){
      print('plotting')
      save_plot(paste(gsub("\\.","",paste(name, "_by_", target, sep="")),'.png',sep=""), graph, base_height = 5.5*i, base_width = 21.6, dpi=900)
      rm(graph)
      gc()
    }
  }
}

#you will need to create a freq_G file in the Length file
setwd("./Length/freq_G")

getGraphs(hk, 'Subunit', seqLen = 9:16, freq=T)
getGraphs(pr8, 'Subunit', seqLen = 9:16, freq=T)
getGraphs(wsn, 'Subunit', seqLen = 9:16, freq=T)
getGraphs(bri, 'Subunit', seqLen = 9:16, freq=T)

getGraphs(rbind(pr8, hk, wsn, bri), 'strain', seqLen = 9:16, freq=T)