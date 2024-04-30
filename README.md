
## Programs used
R 4.2.2

Gephi 0.10

## OS
Windows 11

## Workflow
### 2-class Network creation
Requires 2 files, containing the expression data across all time points of the genes to be considered in both gene classes

Load these files into 2_class_network_creator.R, then run 2_class_gexf_creator on the output

This will create a gexf file containing a visualisation of the network that can be opened in Gephi. To create the visualisations seen in the paper, the layout applied was first Fruchterman Reingold with an area of 100,000, followed by Network Splitter 3d (https://gephi.org/plugins/#/plugin/network-splitter-3d). 

For the effector-TF network, [z] was set equal to the degree of the node for TFs, and to 0 for effectors. 
The paramaters applied were then:

Z-Maximum Level: 50 (equal to the largest [z] in the network)

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80

For all effector-plant network visualisations, [z] was set to 80 for all effectors, and to 0 for all plant genes. 
The parameters applied were then:

Z-Maximum Level: 80

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80

### Adding GO Slim attributes
The list of GO terms and GO slim terms assosciated with plant genes in the network was downloaded from TAIR, and the subset of genes with assosciated terms of interest identified by GO_Slim_attribute_applier_1.R. Presence or absence in this subset was then applied to the data as a binary attribute by GO_Slim_attribute_applier_2.R. 

### Analysis of stress or defence genes
Bootstrapping to analyse an over or underabundance of edges between effector superclusters and genes involved in immunity or defense, as well as presentation of this data was done using Comprehensive Go histogram plotter.R. 

### 3-class network creation
First, load the expression data for nematode TFs and nematode effectors into 2 class network creator.R, and load the resulting file into 3 class TF_effector gefx creator.R.
Next load the expression data for nematode effectors and plant genes into 2 class network creator.R, and load both the resulting file and the outputs of 3 class TF_effector gefx creator.R into 3 class effector_plant gene gefx creator.R.

This will create a gexf file containing a visualisation of the network that can be opened in Gephi. To create the visualisations seen in the paper, the layout applied was first Fruchterman Reingold with an area of 100,000, followed by Network Splitter 3d (https://gephi.org/plugins/#/plugin/network-splitter-3d). 

For the network, [z] was set equal to (degree+80) for TFs, to 80 for effectors, and to 0 for plant genes. 
The paramaters applied were then:

Z-Maximum Level: 133 (equal to the largest [z] in the network)

Z-Distance Factor: 10

Z-Scale: 100

Alfa: 80


