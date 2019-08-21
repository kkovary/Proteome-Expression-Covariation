library(tidyverse)
library(ggpubr)

#############################
#### Import Protein Data ####
#############################

files = list.files('TMTdata/', pattern = 'TargetProtein', full.names = T)

for(i in seq_along(files)){
  if(i == 1){
    rm(protRaw)
  }
  if(!exists('protRaw')){
    protRaw = read_tsv(file = files[i]) %>% mutate(file = files[i])
  } else{
    temp = read_tsv(file = files[i]) %>% mutate(file = files[i])
    colnames(temp) = colnames(protRaw)
    protRaw = rbind(protRaw, temp)
    rm(temp)
  }
}
protRaw = protRaw %>% select(file, everything())



#########################
### Tidy Protein Data ###
#########################

protRaw = protRaw %>% mutate(run = c('ctrl','a','b','c','d','e')[match(protRaw$file, files)]) %>% 
  select(file, run, everything())

protTidy = protRaw %>% gather(key = label, value = intensity, `Ratios: (127) / (126)` : `Ratios: (131) / (126)`) %>% 
  select(run, label, intensity, everything(), -file) %>% 
  mutate(condition = ifelse(run == 'ctrl', 
                            paste0('ctrl',127:131)[match(as.character(label),paste0('Ratios: (',c(127:131), ') / (126)'))],
                            c(0,20,40,60,80)[match(as.character(label),paste0('Ratios: (',c(127:131), ') / (126)'))]),
         intensity = as.numeric(intensity)) %>%
  select(run, label, condition, everything())


###################################
#### Protein Variance Analysis ####
###################################

runVar = protTidy %>% group_by(Accession, run) %>% summarise(cv = sd(intensity, na.rm = T) / mean(intensity, na.rm = T)) %>% mutate(group = 'run')
condVar = protTidy %>% group_by(Accession, condition) %>% summarise(cv = sd(intensity, na.rm = T) / mean(intensity, na.rm = T)) %>% mutate(group = 'cond')

allVar = full_join(runVar, condVar, by = 'Accession') %>% gather(key = group, value = A)

ggplot(runVar, aes(cv, colour = run)) + geom_freqpoly(bins = 100) + theme_bw() + facet_grid(run~.)
ggplot(condVar, aes(cv, colour = condition)) + geom_freqpoly(bins = 100) + theme_bw() + facet_grid(condition~.)

############################
#### Test for Normality ####
############################

ggplot(protTidy, aes(log2(intNorm))) + geom_histogram(bins = 100) + facet_grid(run~label) + theme_bw()
ggqqplot(data = protTidy, x = 'intensity', facet.by = c('run', 'label')) + ylim(0,3)
shapiro.test()

protTidy = protTidy %>% group_by(run, label) %>% mutate(intNorm = intensity / median(intensity, na.rm = T)) %>%
  select(run, label, condition, intensity, intNorm, everything())

ggplot(protTidy, aes(log2(intNorm), colour = run)) + geom_density() + facet_grid(label~.) + theme_bw()
ggplot(protTidy, aes(log2(intNorm), colour = label)) + geom_density() + facet_grid(run~.) + theme_bw()

##############################
#### Protein Correlations ####
##############################

protDF = protTidy %>% filter(run != 'ctrl', condition != 'ctrl') %>% group_by(run) %>%
  unite('cell', c('run','condition')) %>% select(cell, Accession, intensity) %>% 
  filter(!is.na(log2(intensity)), !is.na(intensity)) %>%
  spread(key = cell, value = intensity)

protMat = as.matrix(protDF %>% select(-Accession))
rownames(protMat) = protDF$Accession

corMat = cor(t(protMat), use = 'pairwise.complete.obs')
corMat[which(is.na(corMat))] = 0
pheatmap(corMat, show_colnames = F, show_rownames = F)

pca = prcomp(formula = ~., data = protDF[,2:ncol(protDF)], center = T, scale = T, na.action=na.exclude)
pcaDF = tibble(sample = rownames(pca$rotation), PC1 = pca$rotation[,1], PC2 = pca$rotation[,2])
pcaDF = pcaDF %>% separate(sample, into = c('run','time'))
ggplot(filter(pcaDF, run != 'ctrl'), aes(x = PC1, y = PC2, colour = time)) + geom_point(size = 5, alpha = 0.5)

