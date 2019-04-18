library(tidyverse)

### Import Data ###
files = list.files('TMTdata/', pattern = 'TargetPeptideSpectrumMatch', full.names = T)
rm(pepRaw)
for(i in seq_along(files)){
  if(!exists('pepRaw')){
    pepRaw = read_tsv(file = files[i]) %>% mutate(file = files[i]) %>% select(file, everything())
  } else{
    pepRaw = rbind(pepRaw, 
                   read_tsv(file = files[i], col_names = colnames(pepRaw)) %>% 
                     mutate(file = files[i]) %>% 
                     select(file, everything()))
  }
}
files = tibble(files = files, names = c('ctrl','a','b','c','d','e'))


