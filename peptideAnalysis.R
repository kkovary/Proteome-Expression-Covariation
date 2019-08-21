library(tidyverse)

###################
### Import Data ###
###################

files = list.files('TMTdata/', pattern = 'TargetPeptideSpectrumMatch', full.names = T)

for(i in seq_along(files)){
  if(i == 1){
    rm(pepRaw)
  }
  if(!exists('pepRaw')){
    pepRaw = read_tsv(file = files[i]) %>% mutate(file = files[i], massID = 1:nrow(.))
  } else{
    temp = read_tsv(file = files[i]) %>% mutate(file = files[i], massID = 1:nrow(.))
    colnames(temp) = colnames(pepRaw)
    pepRaw = rbind(pepRaw, temp)
    rm(temp)
  }
}
pepRaw = pepRaw %>% select(file, massID, everything()) %>% filter(!is.na(`126`), 
                                                                  `Quan Usage` == 'Use',
                                                                  Confidence != 'Low')

#################
### Tidy data ###
#################

pepRaw = pepRaw %>% mutate(run = c('ctrl','a','b','c','d','e')[match(pepRaw$file, files)]) %>% 
  group_by(file, `Annotated Sequence`) %>% mutate(pepID = letters[1:n()]) %>% ungroup() %>%
  select(file, run, pepID, everything())

pepTidy = pepRaw %>% gather(key = label, value = intensity, '126':'131') %>% 
  select(run, label, intensity, everything(), -file) %>% 
  mutate(condition = ifelse(run == 'ctrl', 'ctrl',c('ctrl',0,20,40,60,80)[match(label,c(126:131))]),
         intensity = as.numeric(intensity)) %>%
  select(run, label, condition, everything())

################################
### Control Peptide Analysis ###
################################

ctrlPeps = pepTidy %>% group_by(massID, run) %>%
  select(run, label, massID, `Annotated Sequence`, intensity, Confidence) %>%
  group_by(massID, `Annotated Sequence`,run,Confidence) %>%
  summarise(cv = sd(intensity, na.rm = T) / mean(intensity, na.rm = T))

ggplot(ctrlPeps, aes(x = cv, colour = run)) + geom_density() + theme_bw()
ggplot(ctrlPeps, aes(x = Confidence, y = cv, colour = Confidence)) + geom_boxplot() + theme_bw() +
  facet_wrap(~run)

###################################
### Combine Duplicated Peptides ###
###################################

pepTidy %>% group_by(run, label) %>% summarise(dup = sum(duplicated(`Annotated Sequence`)) / n()) %>% 
  unite(.,'Cell', c('run','label')) %>% 
  ggplot(aes(y = 100*dup, x = Cell)) + geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + ylim(0,25)

pepNorm = pepTidy %>% group_by(run, massID) %>% 
  mutate(pepNorm = intensity / intensity[label == 126]) %>%
  select(run, label, condition, pepNorm, everything()) %>% group_by(run, label, `Annotated Sequence`) %>%
  mutate(pepNorm = mean(pepNorm)) %>% filter(!duplicated(`Annotated Sequence`)) %>% ungroup()


test = pepTidy %>% group_by(run, label) %>% 
  ungroup() %>% group_by(`Annotated Sequence`, pepID, run) %>%
  mutate(pepNorm = intensity / intensity[label == 126])

test2 = filter(test, `Annotated Sequence` == 'aAAEGPmk', pepID == 'a') %>%
  group_by(run) %>% mutate(norm = intensity / intensity[label == 126]) %>%
  select(run, label, condition, pepID, intensity, norm)

model = lm(norm ~ run, data = test2)
model = lm(intensity~run+label, data = test2)
test2$resid = model$residuals
ggplot(test2, aes(x = label, y = norm, colour = run)) + 
  geom_point()
ggplot(test2, aes(x = label, y = resid, colour = run)) + 
  geom_point()

testProt = test %>% group_by(run, condition, `Master Protein Accessions`) %>%
  summarise(intensity = median(intensity, na.rm = T)) %>%
  #group_by(run, `Master Protein Accessions`) %>% 
  #mutate(intensity = intensity / intensity[condition == 'ctrl']) %>%
  rename(Accession = `Master Protein Accessions`)

#### Correlation between loading controls
pepCtrl = pepTidy %>% filter(label == 126) %>% select(run, `Annotated Sequence`, pepID, intensity) %>%
  unite('peptide', c(`Annotated Sequence`, 'pepID')) %>% 
  group_by(run) %>% mutate(intensity = intensity / median(intensity)) %>%
  group_by(peptide) %>% mutate(intensity = intensity / median(intensity)) %>% spread(run, intensity)

pepMat = matrix(data = as.vector(unlist(pepCtrl[,2:7])), ncol = 6, nrow = nrow(pepCtrl),dimnames = list(rows = pepCtrl$peptide, cols = names(pepCtrl)[2:7]))
remove = which(is.na(pepMat), arr.ind = T)[,1] %>% as.vector()
pepMat = pepMat[-remove,]

pepCor = cor(pepMat)
mypal = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(201))
pheatmap(pepCor, breaks = seq(-1,1, by = 0.01), color = mypal)

ggpairs(as.data.frame(log2(pepMat)))
ggplot(gather(pepCtrl, 'run','intensity',2:ncol(pepCtrl)), aes(x = log2(intensity), colour = run)) + geom_density()

goldPeps = pepCtrl %>% gather('run','intensity',2:ncol(.)) %>% group_by(run) %>% 
  mutate(sd = sd(log2(intensity), na.rm = T)) %>% rowwise() %>% 
  mutate(intensity = ifelse(abs(log2(intensity)) > 2*sd, NA, intensity)) %>% select(-sd) %>% ungroup()
ggplot(goldPeps, aes(x = log2(intensity), colour = run)) + geom_density()
ggpairs(log2(spread(goldPeps, run, intensity)[,2:7]))

goldPeps = goldPeps %>% group_by(peptide) %>% mutate(keep = ifelse(sum(is.na(intensity)) > 0, F, T)) %>%
  ungroup() %>% filter(keep == T) %>% select(-keep)
goldPeps = goldPeps$peptide %>% unique()

pepTidy %>% group_by(run, label, `Master Protein Accessions`) %>% mutate(keep = ifelse(n()<2, F, T)) %>% select(keep)

protNorm = pepTidy %>% group_by(run, label, `Master Protein Accessions`) %>% 
  mutate(keep = ifelse(n()<2, F, T)) %>% ungroup() %>% filter(keep == T) %>%
  unite(peptide, c(`Annotated Sequence`,'pepID')) %>% filter(peptide %in% goldPeps) %>%
  group_by(run, label) %>% mutate(intensity = intensity / median(intensity, na.rm = T)) %>% ungroup() %>%
  group_by(run, peptide) %>% mutate(intensity = intensity / intensity[label == 126]) %>% ungroup %>%
  filter(label != 126, run != 'ctrl') %>% group_by(run, label, `Master Protein Accessions`) %>%
  summarise(intensity = log2(median(intensity))) %>% rename(Uniprot = `Master Protein Accessions`)

protMat = protNorm %>% unite(cell, c('run', 'label')) %>% spread(cell, intensity)
protMat = matrix(data = as.vector(unlist(protMat[,2:ncol(protMat)])), 
                 ncol = (ncol(protMat)-1), nrow = nrow(protMat), 
                 dimnames = list(rows = protMat$Uniprot, cols = names(protMat)[2:ncol(protMat)]))
remove = which(is.na(protMat), arr.ind = T)[,1] %>% as.vector()
protMat = protMat[-remove,]

protCor = cor(t(protMat))

heat = pheatmap(protCor, breaks = seq(-1,1, by = 0.01), color = mypal, cutree_rows = 10, cutree_cols = 10)
write.csv(cutree(heat$tree_row, k = 10)[cutree(heat$tree_row, k = 10) == 1] %>% as.data.frame(),'/Users/kylekovary/Downloads/cluster.csv')
write.csv(cutree(heat$tree_row, k = 10) %>% as.data.frame(),'/Users/kylekovary/Downloads/background.csv')

pca = prcomp(as.data.frame(protMat))
pcaDF = tibble(cell = rownames(pca$rotation), PC1 = pca$rotation[,1], PC2 = pca$rotation[,2]) %>%
  separate(cell, into = c('run','channel'))
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = channel)) + geom_point(size = 5)
ggplot(tibble(prot = rownames(pca$x),PC1 = pca$x[,1], PC2 = pca$x[,2]), aes(x = PC1, y = PC2, label = prot)) + geom_point()

