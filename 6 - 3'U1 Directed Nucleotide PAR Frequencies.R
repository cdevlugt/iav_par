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

#setwd("") #if not using default directory, set it here
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

dfPreProcess <- function(df, strain='') {
  #takes the influenza data.frames and extracts relevant data
  #also breaks data into series and subunits
  
  #drop rows unneeded downstream
  df <- df[,c(1:8, 11, 31, 16)]
  
  #label series
  df$series <- 'No Realignment'
  df$series[df$Sequence!=df$Trim_Sequence] <- 'Trimmed'
  
  #select data that matches series (trimmed for trimmed, final for no realignment and realigned)
  df <- df[,c(1:8,10,11,12)]
  df$end <- substr(df$Trim_Sequence, nchar(df$Trim_Sequence), nchar(df$Trim_Sequence))
  df$NT_coupure <- substr(df$NT_coupure,1,1)
  df$NT_coupure <- gsub('[^G]','', df$NT_coupure)
  df$end <- paste(df$end, df$NT_coupure, sep='')
  df$Length <- nchar(df$Trim_Sequence)
  df <- df[,c(1:8,11,12,13)]
  
  df$end <- gsub('T', 'U', df$end)
  df$strain <- strain
  
  #decompress data
  df <- subunitSeperate(df, decompress = T)
  
  return(df)
}

hk <- dfPreProcess(hk, 'Hong Kong')
pr8 <- dfPreProcess(pr8, 'Puerto Rico')
wsn <- dfPreProcess(wsn, "WSN")
bri <- dfPreProcess(bri, 'Brisbane')

realignmentFrequency <- function(df=NULL, freq=T, name=''){
  
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
  
  #generate dataframe with realignment frequencies
  re <- data.frame(
    series = c(rep('No_Realignment',8), rep('Trimmed',8), rep('all', 8)),
    end = c('U', 'C', 'G', 'A', 'UG', 'CG', 'GG', 'AG'),
    y=0, 
    x=c(1:8),
    stringsAsFactors = F
  )
  
  #get values
  for(i in 1:8){
    re[i,3] <- nrow(df[df$series=='No Realignment' & df$end==re[i,2],])
    re[(i+8),3] <- nrow(df[df$series=='Trimmed' & df$end==re[(i+8),2],])
    re[(i+16),3] <- nrow(df[df$end==re[(i+16),2],])
  }
  
  #write fold changes
  output <- data.frame(x1 = re[,2][re[,1]=='No_Realignment'],
                       x2 = re[,2][re[,1]=='Trimmed'],
                       NoR = re$y[re[,1]=='No_Realignment'],
                       Trim = re$y[re[,1]=='Trimmed'],
                       all = re$y[re[,1]=='all'],
                       PAR = re$y[re[,1]=='Trimmed']/re$y[re[,1]=='all']*100
                       )
  write.csv(output, paste(name, '.csv', sep=''))
  
  #get frequencies
  if(freq){
    re[1:8,3] <- re[1:8,3] / re[17:24,3] *100
    re[9:16, 3] <- re[9:16,3] / re[17:24,3] *100
  }
  
  if(freq){
    #get global realignment frequency
    reFreq <- nrow(df[df$series=='Trimmed',])/nrow(df) * 100
  } else {
    #global realignment frequency as count
    reFreq <- data.frame(y=(re[17:24,3] * nrow(df[df$series=='Trimmed',])/nrow(df)), x=0)
    reFreq$x <- 1:nrow(reFreq)
  }
  
  re <- re[-c(17:24),] #drop useless data
  
  #factor for order
  re[,1] <- factor(re[,1], levels=c('No_Realignment', 'Trimmed'))
  
  reGraph <- ggplot() +
    geom_bar(data=re, aes(group=series, y=y, x=x, fill=series), stat='identity') + 
    theme_bw() +
    scale_x_continuous(breaks=c(1:8), labels = c('U (-1.2 kcal/mol)', 'C (-1.4 kcal/mol)', 'G (-1.57 kcal/mol)', 'A (-1.73 kcal/mol)',
                                                 'UG (-1.2 kcal/mol)', 'CG (-1.4 kcal/mol)', 'GG (-1.57 kcal/mol)', 'AG (-1.73 kcal/mol)')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab('mRNA Sequence End') +
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
  gc()
  return(reGraph)
}
getGraphs <- function(df, target='Subunit', seqLen=NA, freq=T){
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
  if(target=="Strain" || target=='strain'){
    name <- 'Global'
  } else {
    name <- df$strain[1]
  }
  
  if(length(seqLen)==1 && is.na(seqLen)){
    seqLen <- c(min(df$Length):max(df$Length))
  }
  
  LoN <- realignmentFrequency() #graph legend
  
  print(target)
  
  lenGraph <- function(df, var, sub, target, name, freq){
    for(j in 1:ceiling(length(var)/4)){
      LoN <- realignmentFrequency() #graph legend
      
      print(paste("Length =",var[4*(j-1)+1]))
      titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+1], " Nucleotides", sep=""))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph1 <- plot_grid(title, realignmentFrequency(df[df$Length==var[4*(j-1)+1],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
      
      if(length(var)<(4*(j-1)+2)){
        graph2 <- NULL
      } else {
        print(paste('length =',var[4*(j-1)+2]))
        titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+2], " Nucleotides", sep=""))
        title <- ggdraw() + draw_label(titleText, fontface='bold')
        graph2 <- plot_grid(title, realignmentFrequency(df[df$Length==var[4*(j-1)+2],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
      }
      
      if(length(var)<(4*(j-1)+3)){
        graph3 <- NULL
      } else {
        print(paste('length =',var[4*(j-1)+3]))
        titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+3], " Nucleotides", sep=""))
        title <- ggdraw() + draw_label(titleText, fontface='bold')
        graph3 <- plot_grid(title, realignmentFrequency(df[df$Length==var[4*(j-1)+3],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
      }
      
      if(length(var)<(4*(j-1)+4)){
        graph4 <- NULL
      } else {
        print(paste('length =',var[4*(j-1)+4]))
        titleText <- gsub("\\.","",paste(sub," ", target, "- ", var[4*(j-1)+4], " Nucleotides", sep=""))
        title <- ggdraw() + draw_label(titleText, fontface='bold')
        graph4 <- plot_grid(title, realignmentFrequency(df[df$Length==var[4*(j-1)+4],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
      }
      
      if(j==1){
        graph <- plot_grid(graph1, graph2, graph3, graph4, LoN, ncol=5, nrow=1, rel_widths = c(5,5,5,5,1.6), labels = c(LETTERS[(4*(j-1)+1):(4*(j-1)+4)], NULL))
      } else {
        graph <- plot_grid(graph, plot_grid(graph1, graph2, graph3, graph4, NULL, ncol=5, nrow=1, rel_widths = c(5,5,5,5,1.6), labels = c(LETTERS[(4*(j-1)+1):(4*(j-1)+4)], NULL)), rel_widths=c(1,1), rel_heights=c(j-1, 1), ncol=1, nrow=2)
      }
      
      rm(graph1, graph2, graph3, graph4)
      gc()
      
      if(j==ceiling(length(var)/4)){
        print('plotting')
        save_plot(paste(gsub("\\.","",paste(name,"_", sub, "_by_Length_",sep="")),'.png',sep=""), graph, base_height = 5.5*j, base_width = 21.6, dpi=900)
      }
    }
  }
  gc()
  
  for(i in 1:ceiling(length(subs)/4)){
    print(subs[4*(i-1)+1])
    lenGraph(df[df$target==subs[4*(i-1)+1],], seqLen, sub=subs[4*(i-1)+1], name=name, target=target, freq = freq)
    titleText <- gsub("\\.","",paste(subs[4*(i-1)+1],target, sep=" "))
    title <- ggdraw() + draw_label(titleText, fontface='bold')
    graph1 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+1],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    
    if(length(subs)<(4*(i-1)+2)){
      graph2 <- NULL
    } else {
      print(subs[4*(i-1)+2])
      lenGraph(df[df$target==subs[4*(i-1)+2],], seqLen, sub=subs[4*(i-1)+2], name=name, target=target, freq = freq)
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+2],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph2 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+2],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(length(subs)<(4*(i-1)+3)){
      graph3 <- NULL
    } else {
      print(subs[4*(i-1)+3])
      lenGraph(df[df$target==subs[4*(i-1)+3],], seqLen, sub=subs[4*(i-1)+3], name=name, target=target, freq = freq)
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+3],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph3 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+3],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(length(subs)<(4*(i-1)+4)){
      graph4 <- NULL
    } else {
      print(subs[4*(i-1)+4])
      lenGraph(df[df$target==subs[4*(i-1)+4],], seqLen, sub=subs[4*(i-1)+4], name=name, target=target, freq = freq)
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+4],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph4 <- plot_grid(title, realignmentFrequency(df[df$target==subs[4*(i-1)+4],], freq = freq, paste(name, titleText, sep=' ')),ncol=1, nrow=2, rel_heights = c(0.5,5))
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

#the subunits of each strain
setwd("./Helical Stability/Sequence Ends/Hong Kong")
getGraphs(hk, 'Subunit', seqLen = 9:16, freq=T)
gc()
setwd('..')
setwd("./Puerto Rico")
getGraphs(pr8, 'Subunit', seqLen = 9:16, freq=T)
gc()
setwd('..')
setwd("./WSN")
getGraphs(wsn, 'Subunit', seqLen = 9:16, freq=T)
gc()
setwd('..')
setwd("./Brisbane")
getGraphs(bri, 'Subunit', seqLen = 9:16, freq=T)
gc()

#all strains
setwd("..")
getGraphs(rbind(pr8, hk, wsn, bri), 'strain', seqLen = 9:16, freq=T)
gc()
