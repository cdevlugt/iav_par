require(compiler)
enableJIT(3)

require(ggplot2)
require(ggpubr)
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
    df <- df[,-drops]
    sub <- sub - sum(drops < sub)#one more opperation versus replicating than dropping, but stores less stuff in ram during replication
  }
  df <- df[(rep(row.names(df), df[,sub])),]
  df[,sub] <- colnames(df)[sub]
  colnames(df)[sub] <- name
  return (df)
}
adjustByLength <- function(df, sub=NA, startLen=26, additions=c(29:36), rounds=NA, minAdjust=1, end=F, decompress=F, byAddition=T, byLength=T){
  #takes a data frame containing a start point and lengths of adjustments, and returns a data frame containing
  # the average adjustment either to or from that length. both strings and intger additions may be passed
  #
  #Args:
  # df: the dataframe
  # sub:  the column containing read counts; if unspecified it will take all data
  # startLen: column containing start length
  # additions:  columns containing additions
  # rounds: column containing rounds of addition
  # minAdjust:  minimum adjustment for tabulated data
  # end:  should data be returned as end lengths and average addition to that end length (F)
  #   or as the start length plus average addition from that length(T)
  # decompress: decompresses based on reads in subunit... this requires a large amount of RAM
  # byAddition:  average uncompressed data by addition; if both byAddition and byLength, averages are calculated by global
  # byLength: average uncompressed data by length; if both byAddition and byLength, averages are calculated by global
  #
  #return:
  # dfAdj:  data frame containing lengths and average additions to or from that point depending on specification
  # dfABL:  decompressed file
  
  #determines addition type
  if(is.numeric(df[1,additions[1]])){
    char <- F
  } else if(is.character(df[1,additions[1]])) {
    char <- T
    seqs <- unique(df[,additions[1]])
    if(length(additions) > 1){
      for(i in 2:length(additions)){
        seqs <- c(seqs, df[,additions[i]])
      }
    }
    seqs <- unique(seqs)
    seqs <- seqs[seqs!=""]
    seqs <- seqs[order(nchar(seqs), seqs)]
  } else {
    stop("Invalid type in additions argument")
  }
  
  #if specified removes all rows with rounds of 0
  if(!(is.na(rounds))){
    df <- df[df[,rounds]>0,]
  } else {
    df[,(ncol(df)+1)] <- 0
    rounds <- ncol(df)
  }
  
  #removes rows with less than the minimum adjustment
  if(!char){
    df <- df[rowSums(df[,additions])>=minAdjust,]
  } else {
    df$tempcol <- df[,additions[1]]
    if(length(additions)>1){
      for(i in 2:length(additions)){
        df$tempcol <- paste(df$tempcol, df[,additions[i]], sep="")
      }
    }
    df <- df[nchar(df$tempcol)>=minAdjust,]
    df <- df[,-ncol(df)]
  }
  
  #takes all if sub unspecified, drops all with sub reads == 0
  if(is.na(sub)){
    df[,(ncol(df)+1)] <- 1
    sub <- ncol(df)
  }
  df <- df[df[,sub]>0,]
  
  #removes useless data
  cols <- c(sub, rounds, startLen, additions)
  df <- df[,cols]
  gc()
  
  #gets length adjustments 
  dfAdj <- df[, c(1, 3, 4)]
  if(!char){
    dfAdj[,(ncol(dfAdj)+1)] <- dfAdj[,2] + dfAdj[,3]
  } else {
    dfAdj[,(ncol(dfAdj)+1)] <- dfAdj[,2] + nchar(dfAdj[,3])
  }
  
  #obtain remainder of additions
  if(length(additions)>2){
    for(i in 2:length(additions)){
      dfAdj[,(ncol(dfAdj)+1)] <- df[,(3+i)]
      if(!char){
        dfAdj[,(ncol(dfAdj)+1)] <- dfAdj[,ncol(dfAdj)] + dfAdj[,(ncol(dfAdj)-1)]
        dfAdj[,ncol(dfAdj)][dfAdj[,(ncol(dfAdj)-1)]==0] <- 0
      } else {
        dfAdj[,(ncol(dfAdj)+1)] <- nchar(dfAdj[,ncol(dfAdj)]) + dfAdj[,(ncol(dfAdj)-1)]
        dfAdj[,ncol(dfAdj)][dfAdj[,(ncol(dfAdj)-1)]==""] <- 0
      }
    }
  }
  rm(df) #removes df
  gc()
  
  #creates data frame as: length, addition, end
  dfABL <- dfAdj[,c(1,2,3)]
  colnames(dfABL)[2:3] <- c("Length", "Addition")
  if(!char){
    dfABL$End <- dfABL$Length + dfABL$Addition
  } else {
    dfABL$End <- dfABL$Length + nchar(dfABL$Addition)
  }
  
  #add in remainder of rounds via rbind
  if(length(additions)>1){
    for(i in 2:length(additions)){
      dfABLtemp <- dfAdj[,c(1,(i*2),(i*2+1))]
      colnames(dfABLtemp)[2:3] <- c("Length", "Addition")
      if(!char){
        dfABLtemp$End <- dfABLtemp$Length + dfABLtemp$Addition
        dfABLtemp <- dfABLtemp[dfABLtemp$Addition > 0,]
      } else {
        dfABLtemp$End <- dfABLtemp$Length + nchar(dfABLtemp$Addition)
        dfABLtemp <- dfABLtemp[nchar(dfABLtemp$Addition) > 0,]
      }
      dfABL <- rbind(dfABL, dfABLtemp)
    }
  }
  rm(dfABLtemp)
  gc()
  
  if(!decompress){ #average by multiplying by read count
    if(!char){
      dfABL[,3] <- dfABL[,3] * dfABL[,1]
      
      if(end){
        dfAdj <- data.frame(table(dfABL$End))
        col <- 4
      } else {
        dfAdj <- data.frame(table(dfABL$Length))
        col <- 2
      }
      dfAdj[,1] <- as.numeric(as.character(dfAdj[,1]))
      
      for(i in 1:nrow(dfAdj)){
        dftemp <- dfABL[dfABL[,col]==dfAdj[i,1],]
        dfAdj[i,2] <- sum(dftemp[,3])/sum(dftemp[,1]) 
      }
    } else {
      if(end){
        col <- 4
      } else {
        col <- 2
      }
      
      #lengths to create data frame
      #lengths <- unique(dfABL[,col])
      #lengths <- lengths[order(lengths)]
      lengths <- 9:16 #data range for all analysis
      dfABL <- dfABL[dfABL$Length>=min(lengths) & dfABL$Length<=max(lengths),]
      
      dfAdj <- data.frame(Length=rep(lengths, length(seqs)), stringsAsFactors = F) #reassigns dfAdj for return
      dfAdj$Addition <- ""
      dfAdj <- dfAdj[order(dfAdj[,1]),]
      dfAdj$Frequency <- 0
      dfAdj$Seq_Reads <- 0
      dfAdj$Reads <- 0
      row.names(dfAdj) <- 1:nrow(dfAdj)
      
      globalCount <- sum(dfABL[,1])
      
      if(byAddition == byLength){
        count <- globalCount
      }
      
      dfAdj$LenCount <- 0
      dfAdj$AddCount <- 0
      dfAdj$GlobalCount <- globalCount
      
      rm(globalCount)
      
      for(i in 1:length(lengths)){
        dfi <- dfABL[dfABL[,col]==lengths[i],]
        lengthCount <- sum(dfi[,1])
        if(!byAddition && byLength){
          count <- lengthCount
        }
        for(j in 1:length(seqs)){
          dfj <- dfi[dfi[,3]==seqs[j],] #data frame with specific length (start or end) for a given addition
          addCount <- sum(dfABL[,1][dfABL[,3]==seqs[j]])
          if(byAddition && !byLength){
            count <- addCount
          }
          dfAdj[((i-1)*length(seqs)+j),2] <- seqs[j] #assigns sequence to Addition column
          dfAdj[((i-1)*length(seqs)+j),3] <- sum(dfj[,1])/count #assigns frequency to frequency column
          dfAdj[((i-1)*length(seqs)+j),4] <- sum(dfj[,1]) #total number of reads
          dfAdj[((i-1)*length(seqs)+j),5] <- count #assigns number of reads of that length to Reads
          dfAdj[((i-1)*length(seqs)+j),6] <- lengthCount
          dfAdj[((i-1)*length(seqs)+j),7] <- addCount
          rm(dfj) #close to save RAM
          gc()
        }
        rm(dfi) #close to save RAM
        gc()
      }
    }
    dfAdj[is.na(dfAdj)] <- 0
    
    return(dfAdj) #heatmap return
  } else { #decompressed table
    if(end){
      col <- 4
    } else {
      col <- 2
    }
    rm(dfAdj)
    gc()
    
    dfABL <- dfABL[,c(1,col,3)] #get reads, start or end, adjustment
    dfABL <- dfABL[order(dfABL[,2], nchar(dfABL[,3]), dfABL[,3]),]
    
    #length delimiting
    dfABL <- dfABL[dfABL[,2] >= 8 & dfABL[,2] <= 17,]
    
    dfABL <- subunitDecompress(dfABL, sub=1, drop=F)
    dfABL <- dfABL[,c(2,3)]
    colnames(dfABL) <- c("length", "addition")
    
    gc()
    return(dfABL)
  }
}

#setwd("") #directory if not default
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

#functions to simplify data analysis
seqAdditionHeatmap <- function(df=NA, end=F, threshold=0){
  #Works for data generated from adjustByLength(decompress=F)
  #generates a heatmap from the passed data frame; returns heatmap
  #
  #Args:
  # df: dataframe generated in adjustByLength(decompress=F)
  # end:  as in adjustByLength; specifies if sequence was added after (F) or before (T) length
  # treshold: drops reads below this value
  #
  #return:
  # (get_legend(hm)): the legend for the graph, if df is NA
  # graph: heatmap
  
  if(!end){
    xName <- "Length Before Addition (nt)"
  } else {
    xName <- "Length After Addition (nt)"
  }
  
  #legend
  if(all(is.na(df))){
    df <- data.frame(Length=1, Addition="G", Frequency=1, stringsAsFactors = F)
    returnLegend <- T
  } else {
    returnLegend <- F
    
    xdf <- data.frame(x=min(df[,1]):max(df[,1]), stringsAsFactors = F)
    xdf$y <- 0
    for(i in 1:nrow(xdf)){
      xdf[i,2] <- df[df[,1]==(min(df[,1])+i-1),][1,6]
    }
    xdf[is.na(xdf)] <- 0
    xdf$y <- xdf$y/df[1,8]
    xdf$x <- factor(xdf$x, levels = unique(xdf$x))
    
    ydf <- data.frame(x = unique(df[,2]), y=df[1:(length(unique(df[,2]))),7]/df[1,8], stringsAsFactors = F)
    ydf$x <- factor(ydf$x, levels = unique(ydf$x))
    
    #blocks ggplot from reordering your data
    dfLevels <- min(df[,1]):max(df[,1])
    lens <- unique(df[,1])
    umlen <- dfLevels[!(dfLevels %in% lens)]
    seqs <- unique(df[,2])
    rm(lens)
    if(length(umlen)>0){
      for( i in 1:length(umlen)){
        dftemp <- data.frame(Length=rep(umlen[i]), Addition=seqs, Frequency=0, Seq_Reads=0, Reads=0, LenCount=0, AddCount=0, GlobalCount=0, stringsAsFactors = F)
        df <- rbind(df, dftemp)
        rm(dftemp)
      }
    }
    
    df[,1] <- as.character(df[,1])
    df[,1] <- factor(df[,1], levels=dfLevels)
    rm(dfLevels)
    df[,2] <- factor(df[,2], levels=seqs)
    rm(seqs)
    df <- df[df[,6] >= threshold,]
  }
  
  #graph edges
  t <- -.89
  b <- -.89
  l <- -1.71
  r <- -1.71
  
  hm <- ggplot(data=df, aes(x=Length, y=Addition, fill=Frequency)) +
    geom_tile(colour="white",size=0.25) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme_bw() +
    xlab(xName) +
    ylab("Sequence Added") +
    theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
    theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
    theme(legend.title =element_blank(), legend.text=element_text(size=7)) +
    guides(fill = guide_colorbar(title.vjust = -20)) +
    theme(plot.margin=unit(c(t, r,0,0), 'cm')) #top, right, bottom, left
  
  #if returning legend; i.e., no passed df
  if(returnLegend){
    return(get_legend(hm))
  }
  
  #if not returning legend
  hm <- hm + theme(legend.position='none')
  
  xfill <- 'black'
  yfill <- 'black'
  
  xbar <- ggplot(data=xdf, aes(x=x, y=y, fill=y)) + 
    geom_bar(stat='identity') +
    scale_y_continuous(limits=c(0,1)) +
    #scale_x_continuous(breaks=seq(min(xdf$x), max(xdf$x),  by=1), limits=c(min(xdf$x), max(xdf$x))) +
    ylab("") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    clean_theme() +
    theme(plot.margin = unit(c(t,r,b,l), 'cm')) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme(legend.position='none')
  
  ybar <- ggplot(data=ydf, aes(x=x, y=y, fill=y)) + 
    geom_bar(stat='identity') +
    scale_y_continuous(limits=c(0,1)) +
    #scale_x_discrete(breaks=NULL) +
    ylab("") +
    #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    clean_theme() +
    theme(plot.margin = unit(c(t,r,b,l), 'cm')) + 
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme(legend.position='none') +
    coord_flip()
  
  graph <- plot_grid(xbar, NULL, hm, ybar, ncol=2, nrow=2, rel_widths=c(2,.5,2,.5), rel_heights=c(1,1,2,2), align = 'hv')
  graph <- graph +
    theme(plot.margin = unit(c(0,-r,0,0), 'cm')) +
    theme(aspect.ratio=8/24)
  
  return(graph)
}
generateHeatmaps <- function(df, threshold=12500, subs=c(1:8), startLen=26, additions=c(17:24), rounds=14, minAdjust=1, global=T, byAddition=F, byLength=F, label="Subunit"){
  #Analyzes data in the program and generates heatmaps
  #
  #args:
  # df: a data frame containing data
  # threshold:  threshold for the data (readcounts lower than this will be dropped from the heatmap)
  # subs: columns containing IAV subuunits
  # startLen: column containing the length of the trimmed sequence
  # additions:  columns containing the sequences added in each round
  # minAdjust:  minimum length of adjustments a sequence must have to be considered
  # global: should the global values also be plotted
  # byAddition: needed for adjustByLength
  # byLength: needed for adjustByLength
  # label:  follows the subunit name as specified in the column; i.e., if column name is "HA" and label="Subunit" title of the graph is "HA Subunit"
  #
  #return:
  # graph:  a cowplot of heatmaps for all selected subunits with a legend and titles
  
  if(global){
    df$All <- rowSums(df[,c(subs)])
    subs <- c(ncol(df), subs)
  }
  
  letterNames <- LETTERS[1:(2*length(subs))]
  graphNames <- vector(mode="character", length=0)
  plotLegend <- seqAdditionHeatmap()
  
  #generate data.frames for all subunits and global
  for(i in 1:length(subs)){
    if(global && i==1){
      graphNames[i] <- paste(colnames(df)[subs[i]],"combined", sep="_")
      title <- paste("All ",label,"s", sep="")
    } else {
      if(grepl("\\.", colnames(df)[subs[i]], perl=T)){
        graphNames[i] <- paste(substr(colnames(df)[subs[i]],1,regexpr("\\.",colnames(df)[subs[i]], perl=T)-1),"combined", sep="_")
        title <- paste(substr(colnames(df)[subs[i]],1,regexpr("\\.",colnames(df)[subs[i]], perl=T)-1),label, sep=" ")
      } else {
        graphNames[i] <- paste(colnames(df)[subs[i]],"combined", sep="_")
        title <- paste(colnames(df)[subs[i]],label, sep=" ")
      }
    }
    title <- ggdraw() + draw_label(title, fontface='bold', size = 12)
    assign(graphNames[i], plot_grid(
      seqAdditionHeatmap(
        adjustByLength(df, sub=subs[i], startLen=startLen, additions = additions, rounds = rounds, minAdjust = minAdjust, end=F, decompress = F, byAddition=byAddition, byLength=byLength),
        end = F, threshold = threshold), 
      labels=c(letterNames[((i-1)*1+1)],letterNames[((i-1)*1+1)]), nrow=1, ncol=1, rel_widths = c(1,1))) +
      theme(aspect.ratio=8/24)
    assign(graphNames[i], plot_grid(title, get(graphNames[i]), ncol=1, rel_heights=c(1, 8))) +
      theme(aspect.ratio=9/24)
  }
  
  graph <- get(graphNames[1])
  
  if(length(graphNames)>1){
    for(i in 2:length(graphNames)){
      if(i==2){
        graph <- plot_grid(plot_grid(graph, get(graphNames[i]), ncol=1, rel_heights = c((i-1),1)), plotLegend, ncol=2, rel_widths=c(24,4))
      } else {
        graph <- plot_grid(graph, plot_grid(get(graphNames[i]), plotLegend + theme(legend.position='none'), ncol=2, rel_widths=c(24,4)),ncol=1, rel_heights = c((i-1),1))
      }
    }
  } else {
    graph <- plot_grid(graph, plotLegend, ncol=2, rel_widths=c(24,4))
  }
  return(graph)
}

height <- (9 * 4)/4.5
width <- 28/4.5
dpi <- 900

maxRounds <- max(hk$rounds, pr8$rounds, bri$rounds, wsn$rounds)

#make a large data frame with all strains
pr8$PuertoRico8 <- rowSums(pr8[,c(1:8)])
pr8$HongKong <- 0
pr8$WSN <- 0
pr8$Brisbane <- 0
pr8 <- pr8[,c((ncol(pr8)-3):(ncol(pr8)), grep("Trim_Len", colnames(pr8)), grep("rounds", colnames(pr8)), grep("Seq_Trim", colnames(pr8)))]
dfRounds <- max(pr8$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    pr8[,(ncol(pr8)+1)] <- ""
    colnames(pr8)[ncol(pr8)] <- paste("Seq_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

hk$PuertoRico8 <- 0
hk$HongKong <- rowSums(hk[,c(1:8)])
hk$WSN <- 0
hk$Brisbane <- 0
hk <- hk[,c((ncol(hk)-3):(ncol(hk)), grep("Trim_Len", colnames(hk)), grep("rounds", colnames(hk)), grep("Seq_Trim", colnames(hk)))]
dfRounds <- max(hk$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    hk[,(ncol(hk)+1)] <- ""
    colnames(hk)[ncol(hk)] <- paste("Seq_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

wsn$PuertoRico8 <- 0
wsn$HongKong <- 0
wsn$WSN <- rowSums(wsn[,c(1:8)])
wsn$Brisbane <- 0
wsn <- wsn[,c((ncol(wsn)-3):(ncol(wsn)), grep("Trim_Len", colnames(wsn)), grep("rounds", colnames(wsn)), grep("Seq_Trim", colnames(wsn)))]
dfRounds <- max(wsn$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    wsn[,(ncol(wsn)+1)] <- ""
    colnames(wsn)[ncol(wsn)] <- paste("Seq_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

bri$PuertoRico8 <- 0
bri$HongKong <- 0
bri$WSN <- 0
bri$Brisbane <- rowSums(bri[,c(1:8)])
bri <- bri[,c((ncol(bri)-3):(ncol(bri)), grep("Trim_Len", colnames(bri)), grep("rounds", colnames(bri)), grep("Seq_Trim", colnames(bri)))]
dfRounds <- max(bri$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    bri[,(ncol(bri)+1)] <- ""
    colnames(bri)[ncol(bri)] <- paste("Seq_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

df <- rbind(pr8, hk, wsn, bri)
rm(pr8, bri, hk, wsn)
gc()

colnames(df)[1:4] <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')

setwd("./Sequence Additions") #folder must be created
save_plot('All_Strain_Length.png',generateHeatmaps(df, 0, subs=c(1:4), startLen = grep("Trim_Len", colnames(df)), additions =  grep("Seq_Trim", colnames(df)), rounds = grep("rounds", colnames(df)), minAdjust = 1, global = F, byAddition = F, byLength = T, label="Strain"), base_height = height, base_width = width, dpi=dpi)
save_plot('All_Strain_Addition.png',generateHeatmaps(df, 0, subs=c(1:4), startLen = grep("Trim_Len", colnames(df)), additions =  grep("Seq_Trim", colnames(df)), rounds = grep("rounds", colnames(df)), minAdjust = 1, global = F, byAddition = T, byLength = F, label="Strain"), base_height = height, base_width = width, dpi=dpi)
save_plot('All_Strain_Total_Realignment.png',generateHeatmaps(df, 0, subs=c(1:4), startLen = grep("Trim_Len", colnames(df)), additions =  grep("Seq_Trim", colnames(df)), rounds = grep("rounds", colnames(df)), minAdjust = 1, global = F, byAddition = F, byLength = F, label="Strain"), base_height = height, base_width = width, dpi=dpi)
