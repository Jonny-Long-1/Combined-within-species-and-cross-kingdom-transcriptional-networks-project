#3 layer connectivity analysis plan

#Read in the coefficient files for the effector-TF and effector-plant networks
  #Requires me remaking the network with this data, but I think that's fine

#Iterate through the effectors. Create some kind of data structure to link each effector to the plant genes to which it is connected.
#Iterate through the TFs. 
  #For each effector to which that TF is connected, add all of that effector's assosciated plant genes to a data structure for the TF
  #Once this process is finished, call the unique() function on the set of plant genes.

#You should now have a unique set of plant genes assosciated with each TF on which to perform data analysis

library(tidyverse)

#Read in the list of TFs
TFs <- read.csv("../data_for_networks/Effector_TF_data/Nematode_TFs_with_cyst.csv", stringsAsFactors = FALSE) %>% 
  select(Gene_ID)

#Read in the list of effectors
Effectors <- read.csv("../data_for_networks/Effector_TF_data/Nematode_Effectors_with_cyst.csv", stringsAsFactors = FALSE) %>% 
  select(Gene_ID)

#Read in the list of plant genes
Plant_genes <- read.csv("../data_for_networks/Plant_data_with_attributes/Annotated_difference_of_means.csv", stringsAsFactors = FALSE) %>% 
  select(Gene_ID)

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

#Iterate through all the effectors in the network 
plant_effector_edges <- list() #The list in which the vectors of plant genes connected to each effector will be kept
index <- 1 #The first index of the list
unique_effectors <- unique(TF_edges$Effectors) #The unique set of effectors connected to a TF
for(i in unique_effectors){
  plant_connections <- Plant_edges %>% 
    filter(Nematode_genes==i) %>% #For each effector, generate a vector containing all plant genes to which it is connected.
    unique() #Should be redundant
  plant_effector_edges[[index]] <- plant_connections[,1]  #Add that vector to the list
  index <- index+1 
}

names(plant_effector_edges) <- unique_effectors #name each vector in the list for later identification

plant_TF_edges_df <- data.frame() #The data frame in which the final indirect connections will be stored
unique_TFs <- unique(TF_edges$TFs) #The unique set of all transcription factors in the network

#Iterate through all the TFs 
for(i in unique_TFs){
  to_effectors <- TF_edges %>% 
    filter(TFs == i) %>% #The set of all effectors connected to a TF
    unique() #Should be redundant
  TF_plant_temp <- c() #A vector to contain all the plant genes indirectly connected to that TF
  for(j in to_effectors[,2]){ #Iterate through every effector to which that TF is connected
    plant_temp <- plant_effector_edges[[j]] #The vector containing all plant genes to which that effector is connected
    TF_plant_temp <- c(TF_plant_temp,plant_temp) # Add to the vector for the set of TF-plant connections
  }
  TF_plant_temp <- unique(TF_plant_temp) #duplicates not needed
  TF_df_temp <- rep(i,length(TF_plant_temp)) 
  df_temp <- cbind(TF_df_temp, TF_plant_temp) #A matrix containing the TF of interest in every cell in column 1, and the set of plant genes connected to it in the cells in column 2 
  plant_TF_edges_df <- rbind(plant_TF_edges_df, df_temp) #Add this to the bottom of the dataframe of TF-plant gene connections
}

colnames(plant_TF_edges_df)<- c("TFs","Plant_genes")


TF_degree <- data.frame(table(plant_TF_edges_df$TFs))
#The TFs indirectly connected to the most plants are:
  #Hsc_gene_24270
  #Hsc_gene_7118
  #Hsc_gene_24269 - 3rd most connected to effectors
  #Hsc_gene_11646
  #Hsc_gene_25239
  #Hsc_gene_17197 - Interesting and identified by Priya and Alexis
  #Hsc_gene_13584
  #Hsc_gene_8926
  #Hsc_gene_18314
  #Hsc_gene_11570

#Interestingly, many of the TFs most connected to effectors do not appear in this list
  #I hypothesise that this list contains more J2_10 transcription factors (yet to fully investigate)

ggplot(data=TF_degree, aes(x=Freq))+geom_histogram(bins = 15) + labs(x= "Number of connections to plant genes")
#Shows a heavy right skew like the other degree distributions plotted, but potentially more fat-tailed

print(mean(TF_degree$Freq))
#The mean TF is indirectly connected to 564 plant genes
#The range of connections is 1-1581

print(quantile(TF_degree$Freq,0.2))
#50% of TFs are indirectly connected to 433 plant genes or fewer
#20% of TFs are indirectly connected to 99 plant genes or fewer

TF_degree_distribution <- data.frame(table(TF_degree$Freq)) %>% rename(degree=Var1, number=Freq) %>%  
  mutate("Cumulative_count" = cumsum(number)) %>% 
  mutate("Edges_accounted_for" = as.numeric(as.character(degree))*as.integer(number)) %>% 
  mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

ggplot(data=TF_degree_distribution, aes(x=degree, y=Edges_accounted_for)) + geom_point() +theme(axis.text.x = element_text(angle=90))
#50% of connections come from the 31 (21.4%) most connected TFs
#The most highly connected TF is connected to 18.9% of all plant genes (1581)
  #Or 37.4% of all plant genes with an indirect connection to some TF


#Plant gene statistics
plant_degree <- data.frame(table(plant_TF_edges_df$Plant_genes))

#Since not every effector has an edge to a TF, not every plant gene has an indirect connection to any nematode TFs
  #Only 4229 plant genes have an indirect connection to a TF (72.9% of plant genes with an edge to at least one effector)

ggplot(data=plant_degree, aes(x=Freq))+geom_histogram(bins=20) +labs(x="Number of connections to nematode TFs")
#The degree distribution is once again right skewed, but with an apparent peak at higher degrees and a depression for medium values
print(mean(plant_degree$Freq))
#The average plant gene with connections is connected to 19.3 TFs (13.3%)
#When the plant genes with no connections are taken into account, the mean falls to 9.53 TFs


temp <- data.frame(table(plant_degree$Freq)) %>% rename(degree=Var1, number=Freq) %>% 
  mutate("degree"=as.numeric(as.character(degree))) %>% 
  mutate("number"=as.integer(number))

line_1 <- matrix(c(0,4129),nrow=1,ncol=2) %>% data.frame()
colnames(line_1) <- c("degree","number")

plant_degree_distribution_zeros <- rbind(line_1,temp)

plant_degree_distribution_zeros <- plant_degree_distribution_zeros %>% mutate("Cumulative_count"= cumsum(number)) %>% 
  mutate("Edges_accounted_for"=degree*number) %>% 
  mutate("cumulative_edges"=cumsum(Edges_accounted_for))

#plant_degree_distribution <- data.frame(table(plant_degree$Freq)) %>% rename(degree=Var1, number=Freq) %>%  
  #mutate("Cumulative_count" = cumsum(number)) %>% 
  #mutate("Edges_accounted_for" = as.numeric(as.character(degree))*as.integer(number)) %>% 
  #mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

ggplot(data=plant_degree_distribution_zeros, aes(x=degree, y=number)) +geom_point()
ggplot(data=plant_degree_distribution_zeros, aes(x=degree, y=Edges_accounted_for)) + geom_point()
#50% if all edges come from the 853(20.2%) most highly connected plant genes with connections (or 9.96% most highly connected of all plant genes in the network)
 #Connected to 42 TFs or more



