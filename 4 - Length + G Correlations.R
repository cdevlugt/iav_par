#requires the script 4 - Length + PAR frequency 
require(ggplot2)
require(cowplot)

#import PAR frequencies for length data
setwd('./Length/freq_G/')
hk <- read.csv('Global Hong Kong strain.csv', stringsAsFactors = F)
pr8 <- read.csv('Global Puerto Rico strain.csv', stringsAsFactors = F)
wsn <- read.csv('Global WSN strain.csv', stringsAsFactors = F)
bri <- read.csv('Global Brisbane strain.csv', stringsAsFactors = F)

#add strain and drop useless columns
hk$strain <- 'Hong Kong'
hk <- hk[,grep('x1|PAR|strain', colnames(hk))]

pr8$strain <- 'Puerto Rico'
pr8 <- pr8[,grep('x1|PAR|strain', colnames(pr8))]

wsn$strain <- 'WSN'
wsn <- wsn[,grep('x1|PAR|strain', colnames(wsn))]

bri$strain <- 'Brisbane'
bri <- bri[,grep('x1|PAR|strain', colnames(bri))]

#log transform
hk$PAR <- log(hk$PAR)
pr8$PAR <- log(pr8$PAR)
wsn$PAR <- log(wsn$PAR)
bri$PAR <- log(bri$PAR)

#remove inf
hk <- hk[!is.infinite(hk$PAR),]
pr8 <- pr8[!is.infinite(pr8$PAR),]
wsn <- wsn[!is.infinite(wsn$PAR),]
bri <- bri[!is.infinite(bri$PAR),]

#remove NA
hk <- hk[!is.na(hk$PAR),]
pr8 <- pr8[!is.na(pr8$PAR),]
wsn <- wsn[!is.na(wsn$PAR),]
bri <- bri[!is.na(bri$PAR),]

#generate linear models
hklm <- lm(data = hk, PAR ~ x1)
pr8lm <- lm(data = pr8, PAR ~ x1)
wsnlm <- lm(data = wsn, PAR ~ x1)
brilm <- lm(data = bri, PAR ~ x1)

all <- rbind(pr8, hk, wsn, bri)
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
all$strain <- factor(all$strain, levels = strains)

#get R^2
r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = all[all$strain==strains[i],], PAR ~ x1))$r.squared
}

corPlot <- ggplot(all, aes(x=x1, y=PAR, group=strain, colour=strain)) +
  geom_point(size=0.5) +
  geom_smooth(method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=4),
        axis.text = element_text(size=4), 
        axis.title = element_text(size=6))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65)),
                     labels=c(do.call('expression', #add R^2 values with supersctipts to legend
                                      list(
                                        bquote(.(strains[1])~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#')))),
                                        bquote(.(strains[2])~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#')))),
                                        bquote(.(strains[3])~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#')))),
                                        bquote(.(strains[4])~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#'))))
                                      )))) +
  labs(x='Length (nt)', y='Prime-and-Realign Frequency (ln(%))') +
  scale_x_continuous(limits = c(min(all$x1), max(all$x1)), breaks=unique(all$x1), labels =unique(all$x1))

#steal legend
corLegend <- get_legend(corPlot)
#remove legend
corPlot <- corPlot + theme(legend.position = 'none')


save_plot('correlation_plot.png', plot=plot_grid(corPlot, corLegend, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)

#all up to 13
all <- all[all$x1<=13,]

#get R^2
r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = all[all$strain==strains[i],], PAR ~ x1))$r.squared
}

corPlot13 <- ggplot(all, aes(x=x1, y=PAR, group=strain, colour=strain)) +
  geom_point(size=0.5) +
  geom_smooth(method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=4),
        axis.text = element_text(size=4), 
        axis.title = element_text(size=6))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65)),
                     labels=c(do.call('expression', #add R^2 values with supersctipts to legend
                                      list(
                                        bquote(.(strains[1])~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#')))),
                                        bquote(.(strains[2])~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#')))),
                                        bquote(.(strains[3])~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#')))),
                                        bquote(.(strains[4])~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#'))))
                                      )))) +
  labs(x='Length (nt)', y='Prime-and-Realign Frequency (ln(%))') +
  scale_x_continuous(limits = c(min(all$x1), max(all$x1)), breaks=unique(all$x1), labels =unique(all$x1))

#steal legend
corLegend13 <- get_legend(corPlot13)
#remove legend
corPlot13 <- corPlot13 + theme(legend.position = 'none')


save_plot('correlation_plot_13.png', plot=plot_grid(corPlot13, corLegend13, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)