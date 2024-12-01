# require(devtools)
# install_github('saezlab/OmnipathR')
library(tidyverse)
library(OmnipathR)
library(liana)


curated <- OmnipathR::curated_ligand_receptor_interactions() %>%
  # if no need for complexes use liana decomplexify function
  liana:::decomplexify() %>% 
  # select consensus stimulation
  filter(consensus_stimulation == 1) %>%
  select(source_genesymbol, target_genesymbol, sources, curation_effort, n_resources) %>%
  distinct()

# review resource
curated<-curated %>% filter(source_genesymbol != 'IL10RB') %>% filter(target_genesymbol != 'CD274')

write.csv(curated, '../data/receptor_ligand_association/liana_omni_receptor_ligand_interactions_curated.csv')
