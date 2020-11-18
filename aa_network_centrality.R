

# Build the protein structure networks based on atomic correlation data derived from NMA
# Explore centrality measures

# Identifies important aminoacids in the structure based on their closeness centrality, betweenness centrality
# and belonging to the inner core in a graph derived from the protein (domain) structure.

# -- aurelio.moya@ucl.ac.uk

##### Protein (domain) contact network (PCN) derived from the correlation of normal modes #######
# Normal mode analysis (NMA) is one of the major simulation techniques used to probe large-scale motions in biomolecules.
# The PCN is built based on the cross correlation matrix obtained by NMA.
# The correlation in the fluctuation between two residues i,j measures the dynamic information flow through the edge connecting i and j
# We draw a link between i and j if their fluctuations are correlated (cutoff on correlation matrix = 0.3).
# Therefore, the network edges are weighted by the degree of strength in the correlated motions of contacting residues: a strong correlation
# in the motion between contacting residues implies that knowing how one residue moves better enables one to predict the motion of the other, 
# thereby suggesting a strong information flow between the two residues. 
# The weights can be transformed into ‘‘effective distances’’ between connecting nodes, with strong correlations resulting in shorter effective amino acid - amino acid distances.
# Clarke, D. et al. Identifying Allosteric Hotspots with Dynamics: Application to Inter- and Intra-species Conservation. Structure/Folding and Design 1–13 (2016). doi:10.1016/j.str.2016.03.008

##### Betweenness centrality ########

# is strongly linked to the centrality of nodes (residues) in terms of their capability to transfer signals throughout the protein.
# the depletion of residues having high between-ness centrality values is supposed to interrupt the allosteric communication among regions of the proteins that lie far apart.
# Di Paola, L., De Ruvo, M., Paci, P., Santoni, D. & Giuliani, A. Protein contact networks: an emerging paradigm in chemistry. Chem. Rev. 113, 1598–1613 (2013).
# A node with a large betweenness centrality B controls the flow of information in a network
# Karain, W. I. & Qaraeen, N. I. Weighted protein residue networks based on joint recurrences between residues. BMC Bioinformatics 16, 173 (2015).


##### Closeness centrality ########
# The closeness centrality is connected to the aptitude of a node to participate in the signal transmission
# throughout the protein structure. High closeness centrality nodes were
# demonstrated to correspond to residues located in the active
# site of ligand-binding proteins or to evolutionary conserved residues.
# Di Paola, L., De Ruvo, M., Paci, P., Santoni, D. & Giuliani, A. Protein contact networks: an emerging paradigm in chemistry. Chem. Rev. 113, 1598–1613 (2013).
# A node with a high closeness centrality C plays a principal part in the transmission of information to all other residues in the network 
# Karain, W. I. & Qaraeen, N. I. Weighted protein residue networks based on joint recurrences between residues. BMC Bioinformatics 16, 173 (2015).

# Nodes with large B and C values have been shown to lie in critical regions in proteins, and are usually binding free energy hotspots, or are located in
# the vicinity of hotspots [9, 2, 10]. Protein hotspot residues play a key role in protein-protein interactions.
# Karain, W. I. & Qaraeen, N. I. Weighted protein residue networks based on joint recurrences between residues. BMC Bioinformatics 16, 173 (2015).


##### K-cores ########
# K-cores are a way to measure the core-periphery structure of a network. Protein contact networks have a densely intraconnected core and a relatively sparse connected periphery.
# A k-core in a network is the subnwtwork where nodes are connected with k other nodes.
# Residues in the inner core (the highest k core):
# - are evolutionarily more conserved than those belonging to the lower order (outer) cores, i.e. the periphery.
# - receptor sites for known ligand molecules of most proteins occur in the innermost core.
# - are associated with structural pockets and cavities in binding or active sites.
# - are more probale to hold deleterious or intolerant mutations.
# - conform the stabilization centre of the protein.

library(bio3d)
library(igraph)
library(dplyr)

############## Functions ################


node_in_list <- function(node, lista) {
	
	if(node %in% lista) {
                return(1)
                }
        else {
                return(0)
                }

}

##################

betweenness2pdb <- function(pdb, network) {
        pdb2<- pdb
        pdb2$atom$b <- 0
        for (i in V(network)$name){
                betweenness <- V(network)$betweenness[V(network)$name == i]
                cat(i, betweenness,"\n")
                pdb2$atom$b[pdb2$atom$resno == i] <- betweenness
        }
        
        return(pdb2)
}

##################

kcore2pdb <- function(pdb, network) {
        pdb2<- pdb
        pdb2$atom$o <- 0
        for (i in V(network)$name){
                kcore <- coreness(network)[V(network)$name == i]
                cat(i, kcore,"\n")
                pdb2$atom$o[pdb2$atom$resno == i] <- kcore
        }
        
        return(pdb2)
}

##################



### Load the PDB and the aa interaction network as arguments

cat("specify input files and output name as argument. Example: ./path/to/pdb/files/1t5wA01.pdb 1t5wA01\n")

# arguments to pass from the command line

args <-commandArgs(trailingOnly = TRUE)

pdbfile <- args[1]
outputname <- args[2]
chainselect <- args[3]

print(args)

# Load the pdb file

pdb <- read.pdb(pdbfile)

chainA <-atom.select(pdb,chain=chainselect)
#nowater <-atom.select(pdb,"water",inverse=TRUE)
#ligand <-atom.select(pdb,"ligand",inverse=TRUE)
#nuc <-atom.select(pdb,"nucleic",inverse=TRUE)

# clean pdb (remove ligands and waters). 2 ways, either produce a pdb with just the protein
# 1 select "protein"
prot <- atom.select(pdb,"protein", verbose = T)
pdb <- trim.pdb(pdb, prot)

# 2 clean pdb
pdb <- clean.pdb(pdb=pdb,rm.wat = T, rm.lig = T, rm.h = T, fix.chain = T) ## fix.chain = T removes Warning message: 
# In clean.pdb(pdb = pdb, rm.wat = T, rm.lig = T, rm.h = T) : 
#        PDB is still not clean. Try fix.chain=TRUE and/or fix.aa=TRUE 


# trim the PDB file -> CHAIN, .

# pdb <- trim.pdb(pdb, chainA)
# pdb <- trim.pdb(pdb, nowater)
# pdb <- trim.pdb(pdb, ligand)
# pdb <- trim.pdb(pdb, nuc)

###### NMA and dynamic cross correlation network. Not used now

# Compute the normal modes
modes <- nma(pdb)

# dynamic cross-correlation matrix. This determines the cross-correlations of atomic displacements.
cij <- dccm(modes)
net <- cna(cij, cutoff.cij=0.3)
red <- net$network


### Network analysis

# Generate a data frame with node names as row names

aa_net.df <- as.data.frame(V(red)$name)
names(aa_net.df) <- c("Node")


## nonhub_bottlenecks and inner core nodes
# Bottleneck nodes are those that communicate modules. They can be also hubs, but I'm interested here in those that are not hubs

bet75 <- quantile(betweenness(red, normalized = T),probs = 0.75)
deg75 <- quantile(degree(red, mode = "all"),probs = 0.75)
deg50 <- quantile(degree(red, mode = "all"),probs = 0.5)
close75 <- quantile(closeness(red, normalized = T),probs = 0.75)
core85 <- quantile(coreness(red), probs = 0.85)

betlist <- betweenness(red,normalized = T)
deglist <- degree(red, mode = "all")
closelist <- closeness(red,normalized = T)
corelist <- coreness(red)

betwe.cen <- names(betlist)[betlist >= bet75]
lowdeg <- names(deglist)[deglist <= deg50]
bridges <- names(betwe.cen)[names(betwe.cen) %in% names(lowdeg)]
close.cen <- names(closelist)[closelist >= close75]
cores <- names(corelist)[corelist >= core85]



### Network centrality measures



## Betweeness centrality

aa_net.df$Betweenness <- betweenness(red,normalized = T)


## Closeness centrality

aa_net.df$Closeness <- closeness(red,normalized = T)



## K-core

aa_net.df$k_core <- coreness(red)




aa_net.df$betweenness.cent <- sapply(aa_net.df$Node, node_in_list, betwe.cen)
aa_net.df$closeness.cent <- sapply(aa_net.df$Node, node_in_list, close.cen)
aa_net.df$inner.core <- sapply(aa_net.df$Node, node_in_list, cores)


## Write output

write.table(aa_net.df, file = paste(outputname, ".tsv", sep =""), sep = "\t", row.names = F,col.names = T, quote = F)
write.table(as.character((filter(aa_net.df, inner.core==1))$Node), file= paste(outputname, "_innercore.txt", sep =""), sep = "\t", row.names = F,col.names = F, quote = F)
write.table(as.character((filter(aa_net.df, betweenness.cent==1))$Node), file= paste(outputname, "_bet-central.txt", sep =""), sep = "\t", row.names = F,col.names = F, quote = F)
write.table(as.character((filter(aa_net.df, closeness.cent==1))$Node), file= paste(outputname, "_close-central.txt", sep =""), sep = "\t", row.names = F,col.names = F, quote = F)


## add the centrality to vertex attribute

vertex_attr(red,"betweenness") <- aa_net.df$Betweenness
vertex_attr(red,"closeness") <- aa_net.df$Closenees
vertex_attr(red,"kcore") <- aa_net.df$k_core 

write_graph(red,file = paste(outputname, ".gml", sep =""), format = "gml")

# Write the betweenness in the B-Factor and the coreness in the Occupancy fields of the pdb
pdb.scores <- betweenness2pdb(pdb,red)
pdb.scores <- kcore2pdb(pdb.scores,red)
write.pdb(pdb.scores, file = paste(outputname, "_annotated.pdb", sep =""))
