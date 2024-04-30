
library(tidyverse)

#I really ought to finally get around to replacing the incorrect 3 layer network that I have put off doing anything about for a while now. 

#Intro: read in all the data and lists of edges

#Read in the list of TFs
TFs_full <- read.csv("../data_for_networks/Effector_TF_data/Nematode_TFs_with_cyst.csv", stringsAsFactors = FALSE)
TFs <- TFs_full %>%  select(Gene_ID)

#Read in the list of effectors
Effectors_full <- read.csv("../data_for_networks/Effector_TF_data/Nematode_Effectors_with_cyst.csv", stringsAsFactors = FALSE) 
Effectors<- Effectors_full %>% select(Gene_ID)

#Read in the list of plant genes
Plant_genes_full <- read.csv("../data_for_networks/Plant_data_with_attributes/Annotated_difference_of_means.csv", stringsAsFactors = FALSE)
Plant_genes <- Plant_genes_full %>%  select(Gene_ID)

#Read in the edges between TFs and effectors

TF_edges_large <- read.table("Output_files/output_coefficient_file_effector_TF_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".") %>% 
  rename("TFs"=V1,"Effectors"=V2) %>% select(-V3)

#Keep only the edges involving a TF
TF_edges <- inner_join(TF_edges_large,TFs, by=join_by(TFs==Gene_ID))

#Read in the edges between effectors and plant genes
Plant_edges_large <- read.table("Output_files/output_coefficient_file_RBS_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".")  %>% 
  rename("Plant_genes"=V1,"Nematode_genes"=V2) %>% select(-V3)

#Keep only the edges involving a plant gene
Plant_edges <- inner_join(Plant_edges_large,Plant_genes, by=join_by(Plant_genes==Gene_ID))

#SUGR should be given a value of 3, connected effectors given a value of 2, and plant-genes connected to those given a value of 1

sugr <- "Hsc_gene_14352"

TFs_full_new <- TFs_full %>%  mutate(SUGR_proximity=case_match(Gene_ID,sugr~3,.default = 0)) %>% select(Gene_ID,SUGR_proximity)

#Now we need to find the names of the set of effectors to which sugr is connected

Effector_set <- TF_edges %>% filter(TFs==sugr) %>% select(Effectors)

Effectors_full_new <- Effectors_full %>% mutate(SUGR_proximity=0)

for(i in Effector_set){
  Effectors_full_new <- Effectors_full_new %>% mutate(SUGR_proximity=case_match(Gene_ID,i~2,.default = SUGR_proximity))
}

Effectors_full_new <- Effectors_full_new %>% select(Gene_ID,SUGR_proximity)

#Now we need to find the names of the genes in the set of plant genes to which those effectors are connected

#Function to take an effector name and return the relevant plant gene set
plants_by_effector <- function(eff_gene){
  plant_set = Plant_edges %>% filter(Nematode_genes==eff_gene) %>% select(Plant_genes)
  return(plant_set)
}

all_plant <- lapply(Effector_set$Effectors, plants_by_effector)

all_plant_un <- unlist(all_plant) %>% unique()

Plant_genes_full_new <- Plant_genes_full %>% mutate(SUGR_proximity=0)

for(i in all_plant_un){
  Plant_genes_full_new <- Plant_genes_full_new %>% mutate(SUGR_proximity=case_match(Gene_ID,i~1,.default = SUGR_proximity))
}

Plant_genes_full_new <- Plant_genes_full_new %>% select(Gene_ID,SUGR_proximity)

#Write to csv
TF_csv <- read.csv("For_writeup/3_layer_network_work/Nematode_TFs_cyst.csv", stringsAsFactors = FALSE) %>% select(-SUGR_proximity)
TF_csv <- TF_csv %>% right_join(TFs_full_new)

Effector_csv_1 <-  read.csv("For_writeup/3_layer_network_work/Nematode_effectors_cyst.csv", stringsAsFactors = FALSE) %>% select(-SUGR_proximity)
Effector_csv_1 <-  Effector_csv_1 %>% right_join(Effectors_full_new)

Effector_csv_2 <-  read.csv("For_writeup/3_layer_network_work/Nematode_effectors_non_cyst.csv", stringsAsFactors = FALSE) %>% select(-SUGR_proximity)
Effector_csv_2 <-  Effector_csv_2 %>% right_join(Effectors_full_new)

Plant_csv <-  read.csv("For_writeup/3_layer_network_work/Arabidopsis_data.csv", stringsAsFactors = FALSE) %>% select(-SUGR_proximity)
Plant_csv <-  Plant_csv %>% right_join(Plant_genes_full_new)

#Actually write it back
write.csv(TF_csv,"For_writeup/3_layer_network_work/Nematode_TFs_cyst.csv",row.names = F)
write.csv(Effector_csv_1,"For_writeup/3_layer_network_work/Nematode_effectors_cyst.csv",row.names = F)
write.csv(Effector_csv_2,"For_writeup/3_layer_network_work/Nematode_effectors_non_cyst.csv",row.names = F)
write.csv(Plant_csv,"For_writeup/3_layer_network_work/Arabidopsis_data.csv",row.names = F)
