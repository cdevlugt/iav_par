#requires the script: 6 - 3'U1 Directed Nucleotide PAR Frequencies
#requires the script: 3 - Nucleotide Enrichment Monte Carlo

require(ggplot2)
require(cowplot)

#import strain PAR
#setwd() #if not the default directory, set it here
setwd('./Overview')
parav <- rowMeans(read.csv('PAR rates.csv', stringsAsFactors = F, row.names = 1)[2,])

#####PAR frequencies and sequence end
#import PAR frequencies for sequence end data
dfPreProcess <- function(df, strain=''){
  
  #combine with and without G
  df$Trim[1:4] <- df$Trim[1:4] + df$Trim[5:8]
  df$all[1:4] <- df$all[1:4] + df$all[5:8]
  df$PAR <- df$Trim / df$all * 100
  
  #strain
  df$strain <- strain
  
  #delimit
  df <- df[,grep('PAR|strain', colnames(df))]
  df <- df[1:4,]
  
  #fix x-units
  nrgs <- c(-1.2, -1.4, -1.57, -1.73)
  df$x1 <- nrgs
  
  #remove NA
  df <- df[!is.na(df$PAR),]
  
  return(df)
}
setwd('..')
setwd('./Helical Stability/Sequence Ends') #requires the script 6 - 3'U1 Directed Nucleotide PAR Frequencies
hk <- dfPreProcess(read.csv('Global Hong Kong strain.csv', stringsAsFactors = F), 'Hong Kong')
pr8 <- dfPreProcess(read.csv('Global Puerto Rico strain.csv', stringsAsFactors = F), 'Puerto Rico')
wsn <- dfPreProcess(read.csv('Global WSN strain.csv', stringsAsFactors = F), 'WSN')
bri <- dfPreProcess(read.csv('Global Brisbane strain.csv', stringsAsFactors = F), "Brisbane")

#generate linear models
hklm <- lm(data = hk, PAR ~ x1)
pr8lm <- lm(data = pr8, PAR ~ x1)
wsnlm <- lm(data = wsn, PAR ~ x1)
brilm <- lm(data = bri, PAR ~ x1)

#average and merge data
aav <- rowMeans(cbind(hk$PAR, pr8$PAR, wsn$PAR, bri$PAR))
aav <- data.frame(PAR = aav, strain = 'Average', x1 = c(-1.2, -1.4, -1.57, -1.73))
aavlm <- lm(data = aav, PAR ~ x1)

all <- rbind(pr8, hk, wsn, bri, aav)
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane', 'Average')
all$strain <- factor(all$strain, levels = strains)

#get R^2
r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = all[all$strain==strains[i],], PAR ~ x1))$r.squared
}

corPlot <- ggplot(all, aes(x=x1, y=PAR, group=strain, colour=strain)) +
  geom_point(size=0.7) +
  geom_smooth(method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=3),
        axis.text = element_text(size=3), 
        axis.title = element_text(size=4))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65),
                              "dark grey"),
                     labels=c(do.call('expression', #add R^2 values with supersctipts to legend
                                      list(
                                        bquote(.(strains[1])~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#')))),
                                        bquote(.(strains[2])~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#')))),
                                        bquote(.(strains[3])~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#')))),
                                        bquote(.(strains[4])~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#')))),
                                        bquote(.(strains[5])~(R^{2}*~'='~.(formatC(r2val[5], digits = 3, flag='#'))))
                                      )))) +
  labs(x=bquote(Delta~G[37]^{o}*' (kcal/mol)'), y='Prime-and-Realign Frequency (%)') +
  scale_x_continuous(limits = c(min(all$x1), max(all$x1)),
                     breaks=unique(all$x1),
                     labels =paste(c('U', 'C', 'G', 'A'), ' (', unique(all$x1), 'kcal/mol)', sep='')) +
  scale_y_continuous(limits = c(0, 40), breaks=seq(0, 40, by=10)) +
  geom_hline(yintercept = parav)

#steal legend
corLegend <- get_legend(corPlot)
#remove legend
corPlot <- corPlot + theme(legend.position = 'none')

save_plot('correlation_plot.png', plot=plot_grid(corPlot, corLegend, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)

#####G+1 and normal distinction
#import PAR frequencies for sequence end data
dfPreProcess <- function(df, strain=''){
  
  #strain
  df$strain <- strain
  
  #delimit
  df <- df[,grep('x2|PAR|strain', colnames(df))]
  
  #fix x-units
  nrgs <- c(-1.2, -1.4, -1.57, -1.73, -1.2, -1.4, -1.57, -1.73)
  df$x1 <- nrgs
  
  #remove NA
  df <- df[!is.na(df$PAR),]
  
  return(df)
}

#setwd('./Helical Stability/Sequence Ends') #should already be set above
hk <- dfPreProcess(read.csv('Global Hong Kong strain.csv', stringsAsFactors = F), 'Hong Kong')
pr8 <- dfPreProcess(read.csv('Global Puerto Rico strain.csv', stringsAsFactors = F), 'Puerto Rico')
wsn <- dfPreProcess(read.csv('Global WSN strain.csv', stringsAsFactors = F), 'WSN')
bri <- dfPreProcess(read.csv('Global Brisbane strain.csv', stringsAsFactors = F), "Brisbane")

#average and merge data
aav <- rowMeans(cbind(hk$PAR, pr8$PAR, wsn$PAR, bri$PAR))
aav <- data.frame(x2 = hk$x2, PAR = aav, strain = 'Average', x1 = c(-1.2, -1.4, -1.57, -1.73, -1.2, -1.4, -1.57, -1.73))

#generate linear models without G+1
hklmN <- lm(data = hk[1:4,], PAR ~ x1)
pr8lmN <- lm(data = pr8[1:4,], PAR ~ x1)
wsnlmN <- lm(data = wsn[1:4,], PAR ~ x1)
brilmN <- lm(data = bri[1:4,], PAR ~ x1)
aavlmN <- lm(data = aav[1:4,], PAR ~ x1)

#generate linear models with G+1
hklmG <- lm(data = hk[5:8,], PAR ~ x1)
pr8lmG <- lm(data = pr8[5:8,], PAR ~ x1)
wsnlmG <- lm(data = wsn[5:8,], PAR ~ x1)
brilmG <- lm(data = bri[5:8,], PAR ~ x1)
aavlmG <- lm(data = aav[5:8,], PAR ~ x1)

all <- rbind(pr8, hk, wsn, bri, aav)
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane', 'Average')
all$strain <- factor(all$strain, levels = strains)

all$shape <- 16
all$shape[nchar(all$x2)==2] <- 71

allN <- all[nchar(all$x2)==1,]
allG <- all[nchar(all$x2)==2,]

#get R^2
r2val <- vector(mode = 'numeric', length=2*length(strains))
for(i in 1:length(strains)){
  r2val[2*(i-1)+1] <- summary(lm(data = all[all$strain==strains[i],][1:4,], PAR ~ x1))$r.squared
  r2val[2*(i-1)+2] <- summary(lm(data = all[all$strain==strains[i],][5:8,], PAR ~ x1))$r.squared
}

corPlotNG <- ggplot(all, aes(x=x1, y=PAR, group=strain, colour=strain, shape=shape)) +
  geom_point(size=0.7) +
  scale_shape_identity() +
  geom_smooth(data = allN, aes(x=x1, y=PAR, group=strain, colour=strain), method='lm', se=F, formula = y~x, size=0.5) +
  geom_smooth(data = allG, aes(x=x1, y=PAR, group=strain, colour=strain), method='lm', se=F, formula = y~x, size=0.5, linetype='dashed') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=3),
        axis.text = element_text(size=3), 
        axis.title = element_text(size=4))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65),
                              "dark grey"),
                     labels = c(bquote(.(strains[1])~atop('No G'~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#'))),
                                                          ~'With G'~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#'))))),
                                bquote(.(strains[2])~atop('No G'~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#'))),
                                                          ~'With G'~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#'))))),
                                bquote(.(strains[3])~atop('No G'~(R^{2}*~'='~.(formatC(r2val[5], digits = 3, flag='#'))),
                                                          ~'With G'~(R^{2}*~'='~.(formatC(r2val[6], digits = 3, flag='#'))))),
                                bquote(.(strains[4])~atop('No G'~(R^{2}*~'='~.(formatC(r2val[7], digits = 3, flag='#'))),
                                                          ~'With G'~(R^{2}*~'='~.(formatC(r2val[8], digits = 3, flag='#'))))),
                                bquote(.(strains[5])~atop('No G'~(R^{2}*~'='~.(formatC(r2val[9], digits = 3, flag='#'))),
                                                          ~'With G'~(R^{2}*~'='~.(formatC(r2val[10], digits = 3, flag='#')))))
                     )) +
  labs(x=bquote(Delta~G[37]^{o}*' (kcal/mol)'), y='Prime-and-Realign Frequency (%)') +
  scale_x_continuous(limits = c(min(all$x1), max(all$x1)),
                     breaks=unique(all$x1),
                     labels =paste(c('U', 'C', 'G', 'A'), ' (', unique(all$x1), 'kcal/mol)', sep='')) +
  scale_y_continuous(limits = c(0, 40), breaks=seq(0, 40, by=10)) +
  geom_hline(yintercept = parav)

#steal legend
corLegendNG <- get_legend(corPlotNG)
#remove legend
corPlotNG <- corPlotNG + theme(legend.position = 'none')

save_plot('correlation_plot_NG.png', plot=plot_grid(corPlotNG, corLegendNG, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)

#save both
save_plot('Free Energy and PAR.png', plot= plot_grid(corPlot, corLegend, corPlotNG, corLegendNG, rel_widths = c(2.5, 1, 2.5, 1), ncol=4, nrow=1), base_height = 2, base_width = 7, dpi = 600)

#####Save all correlations
cors <- data.frame(`Puerto Roco`=c(summary(pr8lm)$r.squared, summary(pr8lmN)$r.squared, summary(pr8lmG)$r.squared),
                   `Hong Kong`=c(summary(hklm)$r.squared, summary(hklmN)$r.squared, summary(hklmG)$r.squared),
                   WSN=c(summary(wsnlm)$r.squared, summary(wsnlmN)$r.squared, summary(wsnlmG)$r.squared),
                   Brisbane=c(summary(brilm)$r.squared, summary(brilmN)$r.squared, summary(brilmG)$r.squared),
                   Average=c(summary(aavlm)$r.squared, summary(aavlmN)$r.squared, summary(aavlmG)$r.squared))
row.names(cors) <- c('All', 'No G', 'With G')

write.csv(cors, 'correlations.csv', row.names = T)

#####Difference between the 2 linear models
#generate linear models without G+1
hklmN <- lm(data = hk[1:4,], x1 ~ PAR)
pr8lmN <- lm(data = pr8[1:4,], x1 ~ PAR)
wsnlmN <- lm(data = wsn[1:4,], x1 ~ PAR)
brilmN <- lm(data = bri[1:4,], x1 ~ PAR)
aavlmN <- lm(data = aav[1:4,], x1 ~ PAR)

#generate linear models with G+1
hklmG <- lm(data = hk[5:8,], x1 ~ PAR)
pr8lmG <- lm(data = pr8[5:8,], x1 ~ PAR)
wsnlmG <- lm(data = wsn[5:8,], x1 ~ PAR)
brilmG <- lm(data = bri[5:8,], x1 ~ PAR)
aavlmG <- lm(data = aav[5:8,], x1 ~ PAR)

#generate linear models with G+1
hklm2 <- lm(data = hk[5:8,], PAR ~ x1)
pr8lm2 <- lm(data = pr8[5:8,], PAR ~ x1)
wsnlm2 <- lm(data = wsn[5:8,], PAR ~ x1)
brilm2 <- lm(data = bri[5:8,], PAR ~ x1)
aavlm2 <- lm(data = aav[5:8,], PAR ~ x1)

hm <- mean((predict(hklm2) * hklmN$coefficients[2] + hklmN$coefficients[1]) - predict(hklmG))
pm <- mean((predict(pr8lm2) * pr8lmN$coefficients[2] + pr8lmN$coefficients[1]) - predict(pr8lmG))
wm <- mean((predict(wsnlm2) * wsnlmN$coefficients[2] + wsnlmN$coefficients[1]) - predict(wsnlmG))
bm <- mean((predict(brilm2) * brilmN$coefficients[2] + brilmN$coefficients[1]) - predict(brilmG))
aam <- mean((predict(aavlm2) * aavlmN$coefficients[2] + aavlmN$coefficients[1]) - predict(aavlmG))

am <- mean(c(hm, pm, wm, bm))
ams <- sd(c(hm, pm, wm, bm))


##### Nucleotide Enrichment and Free energy
dfPreProcess <- function(df, strain=''){
  
  #combine with and without G
  df$Trim[1:4] <- df$Trim[1:4] + df$Trim[5:8]
  df$all[1:4] <- df$all[1:4] + df$all[5:8]
  df$PAR <- df$Trim / df$all * 100
  
  #strain
  df$strain <- strain
  
  #delimit
  df <- df[,grep('PAR|strain', colnames(df))]
  df <- df[1:4,]
  
  #fix x-units
  nrgs <- c(-1.2, -1.4, -1.57, -1.73)
  df$x1 <- nrgs
  
  #remove NA
  df <- df[!is.na(df$PAR),]
  
  return(df)
}

#setwd('./Helical Stability/Sequence Ends') #should already be set above
#import and reorder to A, C, G, U
hk <- dfPreProcess(read.csv('Global Hong Kong strain.csv', stringsAsFactors = F), 'Hong Kong')[c(4,2,3,1),]
pr8 <- dfPreProcess(read.csv('Global Puerto Rico strain.csv', stringsAsFactors = F), 'Puerto Rico')[c(4,2,3,1),]
wsn <- dfPreProcess(read.csv('Global WSN strain.csv', stringsAsFactors = F), 'WSN')[c(4,2,3,1),]
bri <- dfPreProcess(read.csv('Global Brisbane strain.csv', stringsAsFactors = F), "Brisbane")[c(4,2,3,1),]

aav <- rowMeans(cbind(hk$PAR, pr8$PAR, wsn$PAR, bri$PAR))

#all nucleotide enrichments
#requires 3 - Nucleotide Enrichment Monte Carlo
setwd('..')
setwd('..')
setwd('./G+1/p-values')
allNT <- data.frame(pr8 = c(t(read.csv('Puerto Rico global obtained z and p.csv')[1,6:9])),
           hk = c(t(read.csv('Hong Kong global obtained z and p.csv')[1,6:9])),
           wsn = c(t(read.csv('WSN global obtained z and p.csv')[1,6:9])),
           bri = c(t(read.csv('Brisbane global obtained z and p.csv')[1,6:9])),
           energy = c(-1.73, -1.4, -1.57, -1.2)
           )
ntav <- rowMeans(allNT[,1:4])

hk$end <- allNT$hk
pr8$end <- allNT$pr8
wsn$end <- allNT$wsn
bri$end <- allNT$bri

aav <- data.frame(PAR = aav,
                  strain = 'Average',
                  x1 = c(-1.73, -1.4, -1.57, -1.2),
                  end = ntav
                  )

#rm and gc
rm(allNT, ntav)
gc()

#merge data
all <- rbind(pr8, hk, wsn, bri, aav)
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane', 'Average')
all$strain <- factor(all$strain, levels = strains)
all$shape <- c(65, 67, 71, 85,
               65, 67, 71, 85,
               65, 67, 71, 85,
               65, 67, 71, 85,
               65, 67, 71, 85)

setwd('..')
setwd('..')
setwd('./Helical Stability/Sequence Ends')

#selection and energy
#get R^2
r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = all[all$strain==strains[i],], end ~ x1))$r.squared
}

corPlotES <- ggplot(all, aes(x=x1, y=end, group=strain, colour=strain)) +
  geom_point(size=0.7) +
  geom_smooth(method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=3),
        axis.text = element_text(size=3), 
        axis.title = element_text(size=4))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65),
                              "dark grey"),
                     labels=c(do.call('expression', #add R^2 values with supersctipts to legend
                                      list(
                                        bquote(.(strains[1])~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#')))),
                                        bquote(.(strains[2])~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#')))),
                                        bquote(.(strains[3])~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#')))),
                                        bquote(.(strains[4])~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#')))),
                                        bquote(.(strains[5])~(R^{2}*~'='~.(formatC(r2val[5], digits = 3, flag='#'))))
                                      )))) +
  labs(x=bquote(Delta~G[37]^{o}*' (kcal/mol)'), y="mRNA (%)") +
  scale_x_continuous(limits = c(min(all$x1), max(all$x1)),
                     breaks=unique(all$x1),
                     labels =paste(c('A', 'C', 'G', 'U'), ' (', unique(all$x1), 'kcal/mol)', sep='')) +
  scale_y_continuous(limits = c(0, 50), breaks=seq(0, 50, by=10))

#steal legend
corLegendES <- get_legend(corPlotES)
#remove legend
corPlotES <- corPlotES + theme(legend.position = 'none')

save_plot('cleavage and energy correlation plot.png', plot=plot_grid(corPlotES, corLegendES, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)


#selection and PAR
#get R^2
r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = all[all$strain==strains[i],], end ~ PAR))$r.squared
}

corPlotSP <- ggplot(all, aes(x=end, y=PAR, group=strain, colour=strain, shape=shape)) +
  geom_point(size=0.7) +
  scale_shape_identity() +
  geom_smooth(method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=3),
        axis.text = element_text(size=3), 
        axis.title = element_text(size=4))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65),
                              "dark grey"),
                     labels=c(do.call('expression', #add R^2 values with supersctipts to legend
                                      list(
                                        bquote(.(strains[1])~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#')))),
                                        bquote(.(strains[2])~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#')))),
                                        bquote(.(strains[3])~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#')))),
                                        bquote(.(strains[4])~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#')))),
                                        bquote(.(strains[5])~(R^{2}*~'='~.(formatC(r2val[5], digits = 3, flag='#'))))
                                      )))) +
  labs(x="mRNA (%)", y="Prime-and-Realign Frequency (%)") +
  scale_x_continuous(limits = c(0, 50), breaks=seq(0, 50, by=10)) +
  scale_y_continuous(limits = c(0, 40), breaks=seq(0, 40, by=10)) +
  geom_hline(yintercept = parav)

#steal legend
corLegendSP <- get_legend(corPlotSP)
#remove legend
corPlotSP <- corPlotSP + theme(legend.position = 'none')

save_plot('cleavage and PAR correlation plot.png', plot=plot_grid(corPlotSP, corLegendSP, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)

#cobine and save plots
save_plot('Selection correlations.png', plot= plot_grid(corPlotSP, corLegendSP, corPlotES, corLegendES, rel_widths = c(2.5, 1, 2.5, 1), ncol=4, nrow=1), base_height = 2, base_width = 7, dpi = 600)