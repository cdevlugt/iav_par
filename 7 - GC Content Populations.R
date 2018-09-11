require(ggplot2)
require(cowplot)
require(compiler)

enableJIT(3)

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
ntConGraph <- function(df, seqName="Purine", Freq=T, genome=NA){
  #Graphs relative or absolute base content of a sequence as calculated by seqComp
  # Box plot with geom_point for data points and mean displayed as a dashed line
  #
  #Args:
  # df: data frame created by seqComp with only the read counts, Sequence, Trim_Sequence, and content counts columns
  # seqName:  y-axis label name 
  # Freq: is the data frequency or occurance
  # genome: genome data to add
  #
  #Return:
  # graph:  box plot with points (alpha values represent read counts) and mean values added
  
  #Df as
  #subunit#Sequence#Trim_Sequence#Rel_GC#Rel_GC (trimmed)
  
  #need 3 series: changes, unchanged, full trimmed
  #force pass subunit, sequence, Trim_Sequence, the 2 Rel or Abs
  df <- df[df[,1]>0,]
  df$Series <- ""
  df$Series[df$Sequence==df$Trim_Sequence] <- "No_Realignment"
  df$Series[df$Series==""] <- "Realigned"
  df <- df[,c(1,ncol(df), 4,5)]
  df <- df[order(df[,2]),]
  
  dfTemp <- df[df$Series=='Realigned',]
  dfTemp$Series <- 'Trimmed'
  
  dfTemp <- dfTemp[,c(1,2,4)]
  df <- df[,c(1,2,3)]
  
  colnames(dfTemp)[3] <- 'y'
  colnames(df)[3] <- 'y'
  
  df[(nrow(df)+1):(nrow(df)+nrow(dfTemp)),] <- dfTemp
  rm(dfTemp)
  gc()
  
  if(!is.na(genome)[1]){gAdd <- T}
  
  if(gAdd){
    genome <- genome[genome[,1]>0,]
    genome$Series <- 'Genome'
    genome <- genome[,c(1,4,3)] #drop sequence
    colnames(genome)[3] <- 'y'
    
    #add to df
    df[(nrow(df)+1):(nrow(df)+nrow(genome)),] <- genome
    rm(genome)
    gc()
  }
  
  #points
  contVals <- unique(df$y)
  
  if(gAdd){
    dfPoints <- rbind(
      data.frame(Series='Genome', y=contVals, x=1, n=0, stringsAsFactors = F),
      data.frame(Series='No_Realignment', y=contVals, x=2, n=0, stringsAsFactors = F),
      data.frame(Series='Trimmed', y=contVals, x=4, n=0, stringsAsFactors = F),
      data.frame(Series='Realigned', y=contVals, x=3, n=0, stringsAsFactors = F)
    )
    
    for(j in 1:length(contVals)){
      dfPoints[(1+(j-1)),4] <- sum(df[,1][df$Series=='Genome' & df$y==contVals[j]])
      dfPoints[(1*length(contVals)+1+(j-1)),4] <- sum(df[,1][df$Series=='No_Realignment' & df$y==contVals[j]])
      dfPoints[(2*length(contVals)+1+(j-1)),4] <- sum(df[,1][df$Series=='Trimmed' & df$y==contVals[j]])
      dfPoints[(3*length(contVals)+1+(j-1)),4] <- sum(df[,1][df$Series=='Realigned' & df$y==contVals[j]])
    }
    
    dfPoints <- dfPoints[dfPoints$n!=0,]
    
    dfPoints$Series <- factor(dfPoints$Series, levels =c('Genome', 'No_Realignment', 'Trimmed', 'Realigned'))
    dfPoints <- dfPoints[order(dfPoints$Series),]
    
    rm(df)
    gc()
    
    #means and SD
    df <- subunitDecompress(dfPoints, sub=4, drop=F)
    
    stats <- data.frame(
      Series=c('Genome', 'No_Realignment', 'Trimmed', 'Realigned'),
      mean=c(mean(df$y[df$Series=='Genome']), mean(df$y[df$Series=='No_Realignment']), mean(df$y[df$Series=='Trimmed']),mean(df$y[df$Series=='Realigned'])),
      sd=c(sd(df$y[df$Series=='Genome']), sd(df$y[df$Series=='No_Realignment']), sd(df$y[df$Series=='Trimmed']), sd(df$y[df$Series=='Realigned'])),
      x=c(1,2,3,4)
    )
    stats$Series <- factor(stats$Series, levels =c('Genome', 'No_Realignment', 'Trimmed', 'Realigned'))
    stats <- stats[order(stats$Series),]
    
    rm(df)
    gc()
    
    #alpha scales
    dfPoints$n[dfPoints$Series=='Genome'] <- dfPoints$n[dfPoints$Series=='Genome'] /max(dfPoints$n[dfPoints$Series=='Genome'])
    dfPoints$n[dfPoints$Series=='No_Realignment'] <- dfPoints$n[dfPoints$Series=='No_Realignment'] /max(dfPoints$n[dfPoints$Series=='No_Realignment'])
    dfPoints$n[dfPoints$Series=='Trimmed'] <- dfPoints$n[dfPoints$Series=='Trimmed'] /max(dfPoints$n[dfPoints$Series=='Trimmed'])
    dfPoints$n[dfPoints$Series=='Realigned'] <- dfPoints$n[dfPoints$Series=='Realigned'] /max(dfPoints$n[dfPoints$Series=='Realigned'])
    
    if(Freq){
      y_name <- paste(seqName," Content (%)", sep="")
    } else {
      y_name <- paste(seqName," Content (nt)", sep="")
    }
    
    graph <-ggplot() +
      geom_point(data=dfPoints, aes(y=y, x=x, group=Series, colour=Series, alpha=n), show.legend=F, position = position_jitter(width=.2)) +
      scale_colour_manual(name='legend',
                          values=c("black", #colours are applied to alphabetical rather than factor order
                                   hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), 
                                   hcl(h=seq(15,375, length=(3+1))[3], c=100, l=65),
                                   hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65))) + 
      geom_errorbar(data=stats, aes(ymin=mean-2*sd, ymax=mean+2*sd, x=x)) +
      stat_summary(fun.y='mean', data=stats, colour="#000000", geom="errorbar", aes(group=Series,x=x, y=mean,ymax=..y.., ymin=..y..), width=.75, linetype="dashed", show.legend=F) +
      ylab(y_name) +
      xlab("") +
      scale_x_continuous(labels=c('Genome', 'No Realignment', 'Pre-Realignment', 'Post-Realignment'), breaks=c(1,2,3,4)) + 
      theme_bw()
    
    
    #add on y-axis numbering
    if(Freq){
      graph <- graph + scale_y_continuous(breaks=seq(0, 100, by=10), limits=c(-1, 101))
    } else {
      ymax <- ceiling(max(dfPoints$y))
      if((ymax %% 2)>0){
        ymax <- ymax + 1
      }
      graph <- graph + scale_y_continuous(breaks=seq(0, ymax, by=2), limits=c(0, ymax))
    }
    
  } else {
    dfPoints <- rbind(
      data.frame(Series='No_Realignment', y=contVals, x=1, n=0, stringsAsFactors = F),
      data.frame(Series='Realigned', y=contVals, x=2, n=0, stringsAsFactors = F),
      data.frame(Series='Trimmed', y=contVals, x=3, n=0, stringsAsFactors = F)
    )
    
    for(j in 1:length(contVals)){
      dfPoints[(1+(j-1)),4] <- sum(df[,1][df$Series=='No_Realignment' & df$y==contVals[j]])
      dfPoints[(length(contVals)+1+(j-1)),4] <- sum(df[,1][df$Series=='Realigned' & df$y==contVals[j]])
      dfPoints[(2*length(contVals)+1+(j-1)),4] <- sum(df[,1][df$Series=='Trimmed' & df$y==contVals[j]])
    }
    
    dfPoints <- dfPoints[dfPoints$n!=0,]
    
    dfPoints$Series <- factor(dfPoints$Series, levels =c('No_Realignment', 'Realigned', 'Trimmed'))
    dfPoints <- dfPoints[order(dfPoints$Series),]
    
    #means and SD
    df <- subunitDecompress(dfPoints, sub=4, drop=F)
    
    stats <- data.frame(
      Series=c('No_Realignment', 'Realigned', 'Trimmed'),
      mean=c(mean(df$y[df$Series=='No_Realignment']), mean(df$y[df$Series=='Realigned']), mean(df$y[df$Series=='Trimmed'])),
      sd=c(sd(df$y[df$Series=='No_Realignment']), sd(df$y[df$Series=='Realigned']), sd(df$y[df$Series=='Trimmed'])),
      x=c(1,2,3)
    )
    
    rm(df)
    
    #alpha scales
    dfPoints$n[dfPoints$Series=='No_Realignment'] <- dfPoints$n[dfPoints$Series=='No_Realignment'] /max(dfPoints$n[dfPoints$Series=='No_Realignment'])
    dfPoints$n[dfPoints$Series=='Trimmed'] <- dfPoints$n[dfPoints$Series=='Trimmed'] /max(dfPoints$n[dfPoints$Series=='Trimmed'])
    dfPoints$n[dfPoints$Series=='Realigned'] <- dfPoints$n[dfPoints$Series=='Realigned'] /max(dfPoints$n[dfPoints$Series=='Realigned'])
    
    if(Freq){
      y_name <- paste(seqName," Content (%)", sep="")
    } else {
      y_name <- paste(seqName," Content (nt)", sep="")
    }
    
    graph <-ggplot() +
      geom_point(data=dfPoints, aes(y=y, x=x, group=Series, colour=Series, alpha=n), show.legend=F, position = position_jitter(width=.2)) +
      geom_errorbar(data=stats, aes(ymin=mean-2*sd, ymax=mean+2*sd, x=x)) +
      stat_summary(fun.y='mean', data=stats, colour="#000000", geom="errorbar", aes(group=Series,x=x, y=mean,ymax=..y.., ymin=..y..), width=.75, linetype="dashed", show.legend=F) +
      ylab(y_name) +
      xlab("") +
      scale_x_continuous(labels=c("No Realignment","Post-Realignment","Pre-Realignment"), breaks=c(1,2,3)) + 
      theme_bw() +
      scale_colour_manual(name='legend', breaks=c('No_Realignment', 'Realigned', 'Trimmed'), 
                          values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Realigned=hcl(h=seq(15,375, length=(3+1))[3], c=100, l=65), Trimmed=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                          labels=c(No_Realignment='No Realignment', Realigned='Realigned', Trimmed='Trimmed'))
    
    #add on y-axis numbering
    if(Freq){
      graph <- graph + scale_y_continuous(breaks=seq(0, 100, by=10), limits=c(-1, 101))
    } else {
      ymax <- ceiling(max(dfPoints$y))
      if((ymax %% 2)>0){
        ymax <- ymax + 1
      }
      graph <- graph + scale_y_continuous(breaks=seq(0, ymax, by=2), limits=c(0, ymax))
    }
  }
  return(graph)
}
generateConGraphs <- function(df, strain="", subs=c(1:8), global=T, seq=10, trim=15, label='Subunit', Freq=T, GConly=F){
  #Generates and saves graphs of A, C, G, T, GC, AT, purine, and pyrimidine content of sequences using seq comp and
  # ntCompGraph
  #
  #Args:
  # df: the data frame with subunits, Sequence, and Trim_Sequence
  # strain: name of the strain; used for title purposes
  # subs: subunits of the influenza virus
  # Global: calculate values for all subunits
  # seq:  column containing sequences
  # trim: column containing trim lengths
  # Freq: frequency (T) analysis or absolute (F) analysis
  # GConly: only analyze GC content
  
  if(global){
    df$All <- rowSums(df[,c(subs)])
    subs <- c(ncol(df), subs)
  }
  
  letterNames <- LETTERS[1:length(subs)]
  graphNames <- vector(mode="character", length=0)
  seqName <- c('A', 'C', 'G', 'T', 'Purine', 'GC', 'Pyrimidine', 'AT', 'GCA')
  
  #extract mini_genome
  genome <- df[,c(subs, grep("mini_genome", colnames(df)))]
  df <- df[,-grep("mini_genome", colnames(df))]
  
  df <- df[,c(subs, seq, trim)]
  subs <- 1:length(subs)
  seq <- max(subs) + 1
  trim <- seq + 1
  
  if(Freq){
    df <- cbind(df, seqComp(df, Sequences = seq, ignore.case = T, Occur = F, Freq = T, purines = T, strong = T))
    df$Rel_Pyrimidine_Sequence <- 100 - df[,(ncol(df)-1)]
    df$Rel_AT_Sequence <- 100 - df[,(ncol(df))]
    df$Rel_GCA_Sequence <- 100 - df[,grep('Rel_T_', colnames(df))[1]]
    df <- cbind(df, seqComp(df, Sequences = trim, ignore.case = T, Occur = F, Freq = T, purines = T, strong = T))
    df$Rel_Pyrimidine_Trim_Sequence <- 100 - df[,(ncol(df)-1)]
    df$Rel_AT_Trim_Sequence <- 100 - df[,(ncol(df))]
    df$Rel_GCA_Trim_Sequence <- 100 - df[,grep('Rel_T_', colnames(df))[2]]
    
    genome <- cbind(genome, seqComp(genome, Sequences = length(subs)+1, ignore.case = T, Occur = F, Freq = T, purines = T, strong = T))
    genome$Rel_Pyrimidine_Sequence <- 100 - genome[,(ncol(genome)-1)]
    genome$Rel_AT_Sequence <- 100 - genome[,(ncol(genome))]
    genome$Rel_GCA_Sequence <- 100 - genome[,grep('Rel_T_', colnames(genome))[1]]
  } else {
    df <- cbind(df, seqComp(df, Sequences = seq, ignore.case = T, Occur = T, Freq = F, purines = T, strong = T))
    df$Abs_Pyrimidine_Sequence <- nchar(df$Sequence) - df[,(ncol(df)-1)]
    df$Abs_AT_Sequence <- nchar(df$Sequence) - df[,(ncol(df))]
    df$Abs_GCA_Sequence <- nchar(df$Sequence) - df[,grep('Abs_T_', colnames(df))[1]]
    df <- cbind(df, seqComp(df, Sequences = trim, ignore.case = T, Occur = T, Freq = F, purines = T, strong = T))
    df$Abs_Pyrimidine_Trim_Sequence <- nchar(df$Trim_Sequence) - df[,(ncol(df)-1)]
    df$Abs_AT_Trim_Sequence <- nchar(df$Trim_Sequence) - df[,(ncol(df))]
    df$Abs_GCA_Trim_Sequence <- nchar(df$Sequence) - df[,grep('Abs_T_', colnames(df))[2]]
    
    genome <- cbind(genome, seqComp(genome, Sequences = length(subs)+1, ignore.case = T, Occur = T, Freq = F, purines = T, strong = T))
    genome$Abs_Pyrimidine_Sequence <- nchar(genome$mini_genome) - genome[,(ncol(genome)-1)]
    genome$Abs_AT_Sequence <- nchar(genome$mini_genome) - genome[,(ncol(genome))]
    genome$Abs_GCA_Sequence <- nchar(genome$mini_genome) - genome[,grep('Abs_T_', colnames(genome))[1]]
  }
  
  if(GConly){
    df <- df[,c(subs,seq,trim, grep("_GC_",colnames(df)))]
    genome <- genome[,c(subs, (length(subs)+1), grep("_GC_",colnames(genome)))]
    seqName <- 'GC'
  }
  
  gc()
  
  #generate graphs for all subunits and global
  for(k in 1:length(seqName)){
    print(seqName[k])
    for(i in 1:length(subs)){
      
      if(global && i==1){
        graphNames[i] <- paste(colnames(df)[subs[i]],"combined", sep="_")
        title <- "All Subunits"
      } else {
        if(grepl("\\.", colnames(df)[subs[i]], perl=T)){
          graphNames[i] <- paste(substr(colnames(df)[subs[i]],1,regexpr("\\.",colnames(df)[subs[i]], perl=T)-1),"combined", sep="_")
          title <- paste(substr(colnames(df)[subs[i]],1,regexpr("\\.",colnames(df)[subs[i]], perl=T)-1),label, sep=" ")
        } else {
          graphNames[i] <- paste(colnames(df)[subs[i]],"combined", sep="_")
          title <- paste(colnames(df)[subs[i]],label, sep=" ")
        }
      }
      print(title)
      title <- ggdraw() + draw_label(title, fontface='bold')
      assign(graphNames[i], plot_grid(title, 
                                      ntConGraph(df[,c(subs[i], seq, trim, grep(paste('_', seqName[k], "_", sep=''), colnames(df)))],
                                                 seqName = seqName[k], 
                                                 Freq=Freq, 
                                                 genome = genome[,c(subs[i], (length(subs)+1), grep(paste('_', seqName[k], "_", sep=''), colnames(genome)))][genome[,subs[i]]>0,]) +
                                        theme(legend.position='none'),
                                      labels=c("",letterNames[i]), ncol=1, rel_heights=c(0.1, 1)))
    }
    graph <- plot_grid(get(graphNames[1]), get(graphNames[2]),get(graphNames[3]),get(graphNames[4]),ncol=4, nrow=1, rel_widths = c(1,1,1,1))
    
    if((length(graphNames)/4)>1){
      for(i in 2:ceiling(length(graphNames)/4)){
        graph <- plot_grid(graph, plot_grid(get(graphNames[(4*(i-1)+1)]),get(graphNames[(4*(i-1)+2)]),get(graphNames[(4*(i-1)+3)]),get(graphNames[(4*(i-1)+4)]),ncol=4, nrow=1, rel_widths = c(1,1,1,1)), ncol=1, rel_heights = c((i-1),1))
      }
    }
    save_plot(paste(strain, " ", seqName[k]," content",".png", sep=""),graph, base_height = 6*ceiling(length(graphNames)/4), base_width = 24, dpi=900)
    rm(graph)
    gc()
  }
}

#setwd("") #dir if not default
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

#set up data to be combined and graphed
pr8[,ncol(pr8)+1] <- rowSums(pr8[,c(1:8)])
pr8[,ncol(pr8)+1] <- 0
pr8[,ncol(pr8)+1] <- 0
pr8[,ncol(pr8)+1] <- 0
pr8 <- pr8[,c((ncol(pr8)-3):ncol(pr8), 11,31,27)]
colnames(pr8)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

hk[,ncol(hk)+1] <- 0
hk[,ncol(hk)+1] <- rowSums(hk[,c(1:8)])
hk[,ncol(hk)+1] <- 0
hk[,ncol(hk)+1] <- 0
hk <- hk[,c((ncol(hk)-3):ncol(hk), 11,31,27)]
colnames(hk)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

wsn[,ncol(wsn)+1] <- 0
wsn[,ncol(wsn)+1] <- 0
wsn[,ncol(wsn)+1] <- rowSums(wsn[,c(1:8)])
wsn[,ncol(wsn)+1] <- 0
wsn <- wsn[,c((ncol(wsn)-3):ncol(wsn), 11,31,27)]
colnames(wsn)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

bri[,ncol(bri)+1] <- 0
bri[,ncol(bri)+1] <- 0
bri[,ncol(bri)+1] <- 0
bri[,ncol(bri)+1] <- rowSums(bri[,c(1:8)])
bri <- bri[,c((ncol(bri)-3):ncol(bri), 11,31,27)]
colnames(bri)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

setwd("./Sequence Composition")
All_strain <- rbind(pr8, hk, wsn, bri)
rm(pr8, hk, wsn, bri)
gc()
generateConGraphs(All_strain, strain = "All_strain rel", subs=c(1:4), seq=5, trim=6, global = F, GConly = T, label = "Strain")
