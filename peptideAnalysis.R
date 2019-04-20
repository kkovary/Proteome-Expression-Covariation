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
  select(file, run, everything())

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




### Combine Peptides to Proteins ###









