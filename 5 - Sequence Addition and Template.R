#requires 5 - Sequence Addition Table
require(compiler)
enableJIT(3)

require(ggplot2)
require(ggpubr)
require(cowplot)
require(doParallel)

#####All data
#get data
#setwd() #if using a different directory
setwd("./Sequence Additions")
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
templates <- c("All", "3'U4", "3'C4")
cluster <- 4
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[j],' nucleotide sequence additions.csv', sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

for(i in 1:length(strains)){
  dfs[[i]]$strain <- strains[i] #label strains
  #limit dfs[[i]] to target | additions | total | percent | strain
  dfs[[i]] <- dfs[[i]][,c(grep("target|additions|total|percent|strain", colnames(dfs[[i]])))]
}

#get the unique additions
additions <- unique(dfs[[1]]$additions)

#merge data for polt
all <- dfs[[1]][grep(paste(templates, collapse = "|"), dfs[[1]]$target),]
if(length(strains)>1){
  for(i in 2:length(strains)){
    all <- rbind(all, dfs[[i]][grep(paste(templates, collapse = "|"), dfs[[i]]$target),])
  }
}

#factor for plot
all$additions <- factor(all$additions, levels=additions)
all$target <- factor(all$target, levels = templates)
write.csv(all, 'Sequence addition by template.csv', row.names = F)
all$percent <- all$percent/100

#all
for(i in 1:length(strains)){
  hm <- ggplot(data=all[all$strain==strains[i],], aes(x=target, y= additions, fill = percent)) +
    geom_tile(colour="white",size=0.25) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme_bw() +
    ggtitle(strains[i]) +
    xlab('Conserved Sequence') +
    ylab('Sequence Added') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=3),
          axis.text = element_text(size=3), 
          axis.title = element_text(size=4),
          plot.title = element_text(size = 5, hjust = 0.5)) +
    guides(fill= guide_colorbar(barheight=4, barwidth = 0.5))
  #steal legend
  hml <- get_legend(hm)
  #remove legend
  hm <- hm + theme(legend.position = 'none')
  
  if(i==1){
    graph <- hm
  } else if(i==length(strains)){
    graph <- plot_grid(graph, hm, hml, ncol = 3, nrow = 1, rel_widths = c(length(strains)-1, 1, 0.5))
  } else {
    graph <- plot_grid(graph, hm, ncol = 2, nrow = 1, rel_widths = c(i-1, 1))
  }
  rm(hm, hml)
}

save_plot('Template and Addition Heatmap.png', graph, base_height = 1, base_width = 4)

#####all 3+
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
templates <- c("All", "3'U4", "3'C4")
cluster <- 4
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[j],' nucleotide sequence additions.csv', sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

for(i in 1:length(strains)){
  dfs[[i]]$strain <- strains[i] #label strains
  #limit dfs[[i]] to target | additions | total | percent | strain
  dfs[[i]] <- dfs[[i]][,c(grep("target|additions|total|percent|strain", colnames(dfs[[i]])))]
  
  #remove all shorter than 3 nt
  dfs[[i]] <- dfs[[i]][nchar(dfs[[i]]$additions)>=3,]
  
  #convert percentages for 3nt +
  for(j in 1:length(templates)){
    dfs[[i]]$percent[dfs[[i]]$target==templates[j]] <- dfs[[i]]$total[dfs[[i]]$target==templates[j]] / sum(dfs[[i]]$total[dfs[[i]]$target==templates[j]]) *100
  }
}

#get the unique additions
additions <- unique(dfs[[1]]$additions)

#merge data for polt
all <- dfs[[1]][grep(paste(templates, collapse = "|"), dfs[[1]]$target),]
if(length(strains)>1){
  for(i in 2:length(strains)){
    all <- rbind(all, dfs[[i]][grep(paste(templates, collapse = "|"), dfs[[i]]$target),])
  }
}

#factor for plot
all$additions <- factor(all$additions, levels=additions)
all$target <- factor(all$target, levels = templates)
write.csv(all, 'Sequence addition by template 3nt.csv', row.names = F)
all$percent <- all$percent/100

#3nt +
for(i in 1:length(strains)){
  hm <- ggplot(data=all[all$strain==strains[i],], aes(x=target, y= additions, fill = percent)) +
    geom_tile(colour="white",size=0.25) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme_bw() +
    ggtitle(strains[i]) +
    xlab('Conserved Sequence') +
    ylab('Sequence Added') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=3),
          axis.text = element_text(size=3), 
          axis.title = element_text(size=4),
          plot.title = element_text(size = 5, hjust = 0.5)) +
    guides(fill= guide_colorbar(barheight=4, barwidth = 0.5))
  #steal legend
  hml <- get_legend(hm)
  #remove legend
  hm <- hm + theme(legend.position = 'none')
  
  if(i==1){
    graph2 <- hm
  } else if(i==length(strains)){
    save_plot('Template and Addition All.png', plot_grid(graph, 
                                                         plot_grid(graph2, hm, NULL, ncol = 3, nrow = 1, rel_widths = c(length(strains)-1, 1, 0.5)),
                                                         ncol=1, nrow=2),
              base_height = 2, base_width = 4)
    graph2 <- plot_grid(graph2, hm, hml, ncol = 3, nrow = 1, rel_widths = c(length(strains)-1, 1, 0.5))
  } else {
    graph2 <- plot_grid(graph2, hm, ncol = 2, nrow = 1, rel_widths = c(i-1, 1))
  }
  rm(hm, hml)
}

save_plot('Template and Addition Heatmap 3nt+.png', graph, base_height = 1, base_width = 4)

