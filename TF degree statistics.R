library(tidyverse)
library(ggrepel)

#Take all the edges in the network, rename the columns appropriately and drop weights which are not needed
edges <- read.table("Output_files/output_coefficient_file_effector_TF_0.975.txt", sep = "\t", header = FALSE, quote = "", dec=".") %>% 
  rename("TFs"=V1,"Effectors"=V2) %>% select(-V3)

#Keep only the edges between TFs and effectors
TFs <- read.csv("../data_for_networks/Effector_TF_data/Nematode_TFs_with_cyst.csv", stringsAsFactors = FALSE) %>% 
  select(Gene_ID)
Effectors <- read.csv("../data_for_networks/Effector_TF_data/Nematode_Effectors_with_cyst.csv", stringsAsFactors = FALSE) %>% 
  select(Gene_ID)

#How many unique effectors are in the network
v1_effectors <- inner_join(edges,Effectors,by=join_by(TFs == Gene_ID))
a<- unique(v1_effectors$TFs)
b<- unique(v1_effectors$Effectors)
c<- c(a,b) %>% unique()

edges_filtered <- inner_join(edges,TFs,by=join_by(TFs == Gene_ID)) 

edges_filtered$TFs %>% unique() %>% length()
edges_filtered$Effectors %>% unique() %>% length()

edges_filtered$TFs %>% unique() %>% length() %>% print()
#The network has 145 connected TFs
edges_filtered$Effectors %>% unique() %>% length()
#The network has 397 connected effectors

#The network has 10,857 edges, of which 1,410 are effector-TF and 9,447 are effector-effector

TF_degree <- data.frame(table(edges_filtered$TFs)) %>% rename(Gene_ID=Var1)
ggplot(data=TF_degree, aes(x=Freq)) + geom_histogram(binwidth = 1)+labs(x="TF degree", y="Frequency") #+ geom_vline(xintercept=5, color="red")# +scale_y_continuous(trans="log10")
ggsave("./degree_dist_figs/TF_degree.png")
print(mean(TF_degree$Freq))
#The mean TF is connected to 9.72 effectors

print(quantile(TF_degree$Freq,0.5))
#50% of TFs are connected to 5 effectors or fewer

TF_degree_distribution <- data.frame(table(TF_degree$Freq)) %>% rename(degree=Var1, number=Freq) %>%  
  mutate("Cumulative_count" = cumsum(number)) %>% 
  mutate("Edges_accounted_for" = as.numeric(as.character(degree))*as.integer(number)) %>% 
  mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

ggplot(data=TF_degree_distribution, aes(x=degree, y=Edges_accounted_for)) + geom_point()

#24.1% of TFs are connected to only a single effector
#50% of all TF-effector edges originate from only 15.2% of TFs
#The most highly connected TF is connected to 13.4% of effectors in the network with TF connections, or 8.04% of 659 effectors (53)

effector_degree <- data.frame(table(edges_filtered$Effectors)) %>% rename(Gene_ID=Var1)
ggplot(data=effector_degree, aes(x=Freq)) + geom_histogram() +labs(x="Effector degree", y="Frequency")#+ geom_vline(xintercept=85, color="red")# + scale_y_continuous(trans="log2")
ggsave("./degree_dist_figs/eff_TF_degree.png")

print(mean(effector_degree$Freq))
print(quantile(effector_degree$Freq, 0.5))

#The mean effector is connected to 2.44 TFs (1.68%)

effector_degree_distribution <- data.frame(table(effector_degree$Freq), stringsAsFactors = FALSE) %>% rename(degree=Var1, number=Freq) 

effector_degree_distribution$degree <- as.numeric(as.character(effector_degree_distribution$degree)) 
effector_degree_distribution<- effector_degree_distribution %>% add_row(degree=0, number=262, .before=1)
 
effector_degree_distribution <- effector_degree_distribution%>%  
  mutate("Cumulative_count" = cumsum(number)) %>% 
  mutate("Edges_accounted_for" = degree*number) %>% 
  mutate("Cumulative_edges" = cumsum(Edges_accounted_for))

#31.2% of effectors are connected to no TFs
#12.7% of effectors are connected to only a single TF
#50% of edges to TFs are acounted for by the top 12.3% most highly connected effectors

ggplot(data=effector_degree_distribution, aes(x=degree, y=Edges_accounted_for)) + geom_point()

#Statistics differing by cluster

#Effector cluster
effectors_in_network <- read_csv("../data_for_networks/Effector_TF_data/Nematode_Effectors_with_cyst.csv") %>% select(Gene_ID,Cluster) 
effectors_in_network$Cluster[is.na(effectors_in_network$Cluster)] <- "none"
effector_degree <- effector_degree %>% right_join(effectors_in_network, by="Gene_ID")
effector_degree$Freq[is.na(effector_degree$Freq)] <- 0


clusters_effector <- unique(effector_degree$Cluster)

#TF cluster
TFs_in_network <- read_csv("../data_for_networks/Effector_TF_data/Nematode_TFs_with_cyst.csv") %>% select(Gene_ID,Cluster) 
TFs_in_network$Cluster[is.na(TFs_in_network$Cluster)] <- "none"
TF_degree <- TF_degree %>% left_join(TFs_in_network, by="Gene_ID")
clusters_TF <- unique(TF_degree$Cluster)

TF_cluster_count <- table(select(TF_degree,Cluster)) %>% data.frame()
#The 5 largest TF clusters are J2_10 with 20 genes (orange), Cyst_J2 with 19 genes (grey), 24 dpi with 15 (dark purple),
# 12 dpi_female 12dpi_male with 11 (blue) and 10hpi and J2 with 9 (yellow and red respectively)


#Now to find the average number of connections to TFs for effectors in each cluster
mean_degree_effector <- matrix(ncol=length(clusters_effector)+1, nrow=2, dimnames = list(c("Number of nodes in cluster","Average degree"), c("overall", clusters_effector)))
mean_degree_effector[1,1]<- length(effector_degree$Gene_ID)
mean_degree_effector[2,1]<- mean(effector_degree$Freq)


for (i in c(1:length(clusters_effector))){
  specific <- clusters_effector[i]
  cluster_degrees <- effector_degree %>% filter(Cluster==specific)
  mean_degree_effector[1,1+i]<- length(cluster_degrees$Freq)
  mean_degree_effector[2,1+i] <- mean(cluster_degrees$Freq)
}

mean_degree_effector <- mean_degree_effector[,-1] %>% t() %>% data.frame()

ggplot(data=mean_degree_effector, aes(x=Number.of.nodes.in.cluster, y=Average.degree)) + geom_point() +geom_text_repel(aes(label=clusters_effector))

#Of the clusters containing more than 10 nodes, J2_10(orange) is the most connected at 6.03, 12_fem_12_male is second at 3.80, 
# J2_10_48 (brown) third at 3.24, and then 10_hpi and 10hpi_48hpi (yellow and green) at 2.67 and 2.49 respectively.

