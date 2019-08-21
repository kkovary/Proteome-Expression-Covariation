library(tidyverse)
library(readxl)
library(stringr)

# Load in and combine data
high <- read_xlsx(path = '/Users/kylekovary/Documents/GitHub/Proteome-Expression-Covariation/high_vs_low_var_GO.xlsx',
                  sheet = 1) %>%
  mutate(Var = 'high')
low <- read_xlsx(path = '/Users/kylekovary/Documents/GitHub/Proteome-Expression-Covariation/high_vs_low_var_GO.xlsx',
                 sheet = 2) %>%
  mutate(Var = 'low')

highLow <- full_join(high, low, by = c('Definitions')) %>% 
  dplyr::select(Definitions, Counts.x, Counts.y) %>%
  dplyr::rename(High_Var = Counts.x, Low_Var = Counts.y) %>%
  gather(key = Var, value = Counts, High_Var:Low_Var)

highLow[which(is.na(highLow$Counts)), 'Counts'] = 0

# Tidy and calculate percentages
highLow <- highLow %>% filter(!Definitions %in% c('biological_process',
                                                 'molecular_function',
                                                 'cellular_component',
                                                 'cell',
                                                 'cytoplasm',
                                                 'intracellular')) %>% 
  group_by(Var) %>% 
  dplyr::mutate(per = Counts / sum(Counts)) %>%
  ungroup() #%>% mutate(Definitions = fct_reorder(Definitions, per))

# "Other" categories
notOther <- highLow %>% group_by(Definitions) %>% summarise(sumPerc = sum(per)) %>%
  filter(sumPerc > 0.04)
notOther <- as.character(notOther$Definitions)

other <- highLow %>% filter(!Definitions %in% notOther) %>% group_by(Var) %>%
  summarise(Counts = sum(Counts), per = sum(per)) %>% mutate(Definitions = 'other') %>%
  dplyr::select(Definitions, Var, Counts, per)

highLow <- highLow %>% filter(Definitions %in% notOther)
highLow <- rbind(highLow, other)
highLow$Definitions = str_wrap(highLow$Definitions, width = 20)

# Plots
ggplot(highLow, aes(x = factor(1), y = per, fill = Definitions)) +
  geom_bar(stat = 'identity') + facet_grid(~Var) + 
  guides(fill=FALSE)

ggplot(highLow, aes(x = factor(1), y = per, fill = Definitions)) +
  geom_bar(width = 1) + coord_polar('y') + facet_grid(~Var, scales = 'free') + 
  guides(fill=FALSE)

ggplot(highLow, aes(x = 2, y = per, fill = Definitions)) +
  geom_bar(width = 1, stat = "identity", color = "white") + facet_wrap(~Var) + 
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) + theme_void()# + guides(fill=FALSE)
