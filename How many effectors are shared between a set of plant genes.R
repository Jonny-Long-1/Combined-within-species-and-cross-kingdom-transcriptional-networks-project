#This code should take a particular set of plant genes and return a data frame containing 2 columns
#The first column should contain the names of every effector to which a gene in that set of plant genes is connected
#The second column should indicate how many of the genes in the original plant gene set that effector is connected to 

library(tidyverse)

edges <- read.table("For_writeup/output_coefficient_file_659_plant_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".") %>% 
  select(-V3) %>%  rename("Plant_genes"=V1, "Effectors"=V2)

Input_1_wide <- read.csv("../data_for_networks/Plant_data_with_attributes/Nematode_effectors.csv", stringsAsFactors = FALSE) %>% 
  select(Gene_ID,Cluster)

edges_cluster <- left_join(edges,Input_1_wide, by = join_by(Effectors==Gene_ID))

plant_gene_set <- c("AT1G26420",
                    "AT1G16370",
                    "AT1G47560",
                    "AT2G18750",
                    "AT2G31751",
                    "AT2G34350",
                    "AT3G19710",
                    "AT4G30340",
                    "AT4G36820",
                    "AT5G23820")

#Finds the total set of effectors to which those plant genes are connected

effector_set <- c()

for(gene in plant_gene_set){
  temp <- edges %>% filter(Plant_genes==gene)
  effector_set<- c(effector_set,temp$Effectors)
}

effector_set <- effector_set %>% unique()

#This is a function that takes the name of an effector and a set of plant genes, and returns how many of those plant genes that effector is connected to
count_edges <- function(effector, plant_set){
  n= 0
  for(gene in plant_set){
    temp_2 = edges %>% filter(Plant_genes == gene)
    if(effector %in% temp_2$Effectors){
      n = n+1
    }
  }
  return(n)
}

effector_edge_count <- lapply(effector_set, count_edges, plant_set = plant_gene_set) %>% unlist()

final <- cbind(effector_set,effector_edge_count)

table(final[,2])

#Of the 10 most connected plant genes, the vast majority of their connections are shared between all 10 (105 effectors)
#There are only 29 effectors connected to one or more of these plant genes but not all 10


set <- edges_cluster %>%  filter(Plant_genes == "AT1G26420")

table(set[,3])
