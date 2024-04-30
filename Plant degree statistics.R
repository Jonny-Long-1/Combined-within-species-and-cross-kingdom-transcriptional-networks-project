#This code will compare the connectivity of genes from different groups
library(tidyverse)

#Take all the edges in the network, rename the columns appropriately and drop weights which are not needed
edges <- read.table("Output_files/output_coefficient_file_RBS_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".")  %>% 
  rename("Plant_genes"=V1,"Nematode_genes"=V2) %>% select(-V3)
print(length(edges))

#Extract only the edges between a plant gene and an effector
edges_efef<- edges
edges<- edges[grep("AT", edges$Plant_genes),] 
print(length(edges))

plant_degree <- data.frame(table(edges$Plant_genes))
ggplot(data=plant_degree, aes(x=Freq)) + geom_histogram(binwidth = 1) +labs(x="Plant gene degree", y="Frequency") #+ geom_vline(xintercept=6, color="red")# +scale_y_continuous(trans="log10")
#More plant genes are correlated with one effector than any other value
print(mean(plant_degree$Freq))
#The mean degree for plant genes is 23.7

print(quantile(plant_degree$Freq,0.5))
#Results would be different if disconnected nodes were included

#50% of plant genes are connected to 6 effectors or fewer
#18.1% of plant genes in the network are connected to only a single effector

plant_degree_distribution <- data.frame(table(plant_degree$Freq)) %>% rename(degree=Var1, number=Freq) %>%  
  mutate("Cumulative_count" = cumsum(number)) %>% 
  mutate("Edges_accounted_for" = as.numeric(as.character(degree))*as.integer(number)) %>% 
  mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

ggplot(data=plant_degree_distribution, aes(x=degree, y=Edges_accounted_for)) + geom_point()

#50% of edges come from only the 12% most highly connected plant genes (connected to 78 effectors or more)
#The most highly connected plant gene is connected to 18.4% of all effectors (121)


effector_degree <- data.frame(table(edges$Nematode_genes)) %>% rename(Gene_ID=Var1)
ggplot(data=effector_degree, aes(x=Freq)) + geom_histogram(binwidth = 10) + geom_vline(xintercept=85, color="red")# + scale_y_continuous(trans="log2")
print(mean(effector_degree$Freq))
print(quantile(effector_degree$Freq, 0.5))

effector_degree_efef <- data.frame(table(edges_efef$Nematode_genes)) %>% rename(Gene_ID=Var1)
ggplot(data=effector_degree_efef, aes(x=Freq)) + geom_histogram(binwidth = 10) +labs(x="Effector degree", y="Frequency") #+ geom_vline(xintercept=85, color="red")# + scale_y_continuous(trans="log2")
print(mean(effector_degree_efef$Freq))
print(quantile(effector_degree_efef$Freq, 0.5))

#The mean degree of effectors is 209
#19% of effectors are connected to 13 plant genes or fewer 
#50% of effectors are connected to 85 plant genes or fewer 
#2 effectors (0.3%) are connected to only a single plant gene 

effector_degree_distribution <- data.frame(table(effector_degree$Freq), stringsAsFactors = FALSE) %>% rename(degree=Var1, number=Freq) %>%  
  mutate("Cumulative_count" = cumsum(number)) %>% 
  mutate("Edges_accounted_for" = as.numeric(as.character(degree))*number) %>% 
  mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

#50% of the edges come from the top 80 (12.1%) of the most highly connected effectors (connected to 630 plant genes or more)
#The most highly connected effectors are each connected to 1011 (17.4% of all) plant genes

ggplot(data=effector_degree_distribution, aes(x=degree, y=Edges_accounted_for)) + geom_point()

effector_degree_distribution_efef <- data.frame(table(effector_degree_efef$Freq), stringsAsFactors = FALSE) %>% rename(degree=Var1, number=Freq) %>%  
  mutate("Cumulative_count" = cumsum(number)) %>% 
  mutate("Edges_accounted_for" = as.numeric(as.character(degree))*number) %>% 
  mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

#By effector cluster
effectors_in_network <- read_csv("../data_for_networks/Plant_data_with_attributes/Nematode_effectors.csv") %>% select(Gene_ID,Cluster) 
effectors_in_network$Cluster[is.na(effectors_in_network$Cluster)] <- "none"
effector_degree <- effector_degree %>% left_join(effectors_in_network, by="Gene_ID")
clusters <- unique(effector_degree$Cluster)

mean_degree <- matrix(ncol=length(clusters)+1, nrow=2, dimnames = list(c("Number of nodes in cluster","Average degree"), c("overall", clusters)))
mean_degree[1,1]<- length(effector_degree$Gene_ID)
mean_degree[2,1]<- mean(effector_degree$Freq)


for (i in c(1:length(clusters))){
  specific <- clusters[i]
  cluster_degrees <- effector_degree %>% filter(Cluster==specific)
  mean_degree[1,1+i]<- length(cluster_degrees$Freq)
  mean_degree[2,1+i] <- mean(cluster_degrees$Freq)
}

mean_degree_df<- data.frame(t(mean_degree[,-1]))

#ggplot(data=mean_degree_df, aes(x=Number.of.nodes.in.cluster, y=Average.degree)) + geom_point() +geom_text_repel(aes(label=clusters))
#10hpi (yellow) and J2_10hpi(orange) have average degree almost twice as high as any other cluster. 
#48 hpi(dark green) has a very low average degree, especially for its size. 
#12dpi_male and 12dpi_female(blue) also both have relatively low degree. 
#10hpi_48hpi(light green) is relatively intermediate between 10hpi and 48hpi in both size and degree. 
