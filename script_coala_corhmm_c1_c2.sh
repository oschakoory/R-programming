#!/usr/bin/Rscript

###########################################################################

# Data simulation to classify genotypes using hidden markov model (HMM)
# 7 February 2020
# Author : Oshma Chakoory (Internship M2 Bioinformatics)

###########################################################################

#The aim of data simulation is to optimse the method of classification so as to minimise % error before tetsting on real dataset
#First we will use coala library to simulate different dataset (varying in number of leaves) 
#All these simulated dataset will be entered in corhmm library to test its efficiency 
#In corhmm, we will vary the parameters so that the number of errors between simulated dataset and corhmm datset remains relatively low 

args <- commandArgs() # get arguments

print ("setting working directory...")
#path <-"/mnt/grille/2020_oshma_corhmm/"
path <-"/home/sge/2020_oshma_corhmm/"
setwd(path)

print ("loading libraries...")
library(coala) #library to simulate dataset
library(ggtree) #library to generate phylogenetic tree
library(ape) #library to generate phylogenetic tree
library(corHMM) #library to classify genotypes

#start time
print ("start time...")
ptm <- proc.time()

###########################################################################
#
#             Data simulation with R library 'coala' 
#
###########################################################################


#we provide coala sample size (k), estimated number of mutation (mu), rate of mutation (theta).
#coala simuale the tree using the model of coalescent.
#we also provide emission probabilities of 2 genotypes to determine 2 phenotypes
#expected results: Phylogenetic tree (format phylo to be used in corHMM), phenotypes and genotypes

# sample size of population N
# also the number of leaves on the phylogenetic tree
print ("starting coala...")
print ("initialise parameters...")
sim <- as.numeric(args[6])
k <- as.numeric(args[7])

# estimated number of mutation per generations on 1 locus
mu <- 7e-6 #Number of mutation found in Borrelia

# Size of the population
N <- 1e5

# rate of mutation 
theta <- 4.1

# numbre of sites on locus
nb_site <- 1000000

# initialisation
m1 <- coal_model(sample_size=k,loci_number=1,loci_length=nb_site,ploidy=1)
m1 <- m1+sumstat_trees(name="arbre")
m1 <- m1+feat_mutation(rate=theta,model="IFS") # ajout de mutations selon un modèle de sites infinis
m1 <- m1+sumstat_seg_sites()

# simulation
print ("tree simulation...")
simulation_m1 <- simulate(m1,nsim=1,seed=NULL)
arbre <- read.tree(text=simulation_m1$arbre[[1]]) # récupérer l'arbre
arbre$node.label <- c(1:arbre$Nnode)

#retrieve all sites of mutation
snp <- simulation_m1$seg_sites[[1]]

nb_snp_row <- apply(snp,1,sum) #nb of SNP par row
nb_snp_col <- apply(snp,2,sum) #nb of SNP per column

#nb_snp_row%%2 #composed of 0 and 1

## create data.frame to annotate tree with haplotypes
haplotype <- data.frame("id"=1:k,"haplotype"=nb_snp_row%%2) #nb_snp_row%%2 donne des 0 et 1

#create a df to store all information from corHMM
#data_class2 <- data.frame(simulation=character(), loglik=numeric(), AIC=numeric(), AICc=numeric(), nerror=numeric(), percent_error=numeric(), sum_unmatchedP=numeric())

#create a directory to save all results for class 1 & 2 genotypes
dir.create(paste0("snp_plus1"),showWarnings = FALSE)
dir.create(file.path(paste0("snp_plus1"),paste0("simulation_", sim)),showWarnings = FALSE)
dir.create(file.path(paste0("snp_plus1"),paste0("simulation_", sim),paste0("class_1")),showWarnings = FALSE)
dir.create(file.path(paste0("snp_plus1"),paste0("simulation_", sim),paste0("class_2")),showWarnings = FALSE)
#######################
#
#Emission probabilities
#
#######################

#define probability emission for class 2 genotype
print("defining emission probabilities...")
EP <- load(args[8])

#add G in from of 0 and 1 to match EP table
genotype <- paste("G",nb_snp_row%%2, sep="")

#view new results
#df[genotype,]

#function to draw phenotype
print("drawing phenotypes based on emission probabilities...") 
tirage_phenotype <- function(p) { return(sample(x=c(0,1),size=1,prob=p)) }
#tirage_phenotype( df[genotype,][1,])
#tirage_phenotype( df[genotype,][2,])

#list of phenotype for each genotype
phenotype <- apply(df[genotype,],1,tirage_phenotype)
phenotype_df <-data.frame("id"=1:k,phenotype)
geno_pheno <- cbind(haplotype,phenotype=phenotype_df$phenotype)

########################################################################
#
#         Classification of genotypes using R library 'corHMM' 
#
########################################################################

print("starting corHMM...")

#from coala we have:
#arbre format phylo
#phenotype_df: list of 0 and 1
#corhmm takes these 2 parameters to calculate transition probabilities of hidden states
#expected results: Transition probabilities of hidden states ($tip.states) and probabilities of internal nodes ($states)

###################
#
#corHMM class 1
#
###################

#run corhmm with arbre and phenotype_df
#rate.cat: 1 class of hidden states
#aim of running corhmm with rat.cat = 1 is to compare $AICc with class 2's $AICc
#expected results: class 2's AICc is greater than the $AICc of class 1 showing that the estimation of hidden sates is better if we have more than 1 hidden states
print("running corHMM class 1...")
class1_corhmm <- corHMM(arbre,phenotype_df,rate.cat=1,node.states="marginal",optim.method="twoStep", nstarts=100, sann.its=5000)
#class1_corhmm <- corHMM(arbre,phenotype_df,rate.cat=1,node.states="marginal",optim.method="subplex", nstarts=10)

#after running corhmm, we get
# $solution representing the Q matrix 
# $solution.se representing standard deviation of the matrix Q (diagn=T for corhmm to calculate)
# $tip.states representing the transition proabilities of hidden states of leaves
# $states representing the internal nodes probabilities

#define file name for class 1 results
file_name <- as.character(args[9])

my_results_class1 <- list(Emission_probability=df,corHMM_class1 = class1_corhmm,coala_phenotype= phenotype_df,arbre_coala=arbre)

print("saving results for corHMM class 1...")
save(my_results_class1,file=file.path(paste0("snp_plus1"),paste0("simulation_", sim),paste0("class_1"),paste0("res_",file_name,".RData")))


##################
#
#corHMM class 2
#
##################

#run corhmm with arbre and phenotype_df
#rate.cat: 2 classes of hidden states
print("running corHMM class 2...")
class2_corhmm <- corHMM(arbre,phenotype_df,rate.cat=2,node.states="marginal",optim.method="twoStep", nstarts=100, sann.its=5000)
#class2_corhmm <- corHMM(arbre,phenotype_df,rate.cat=2,node.states="marginal",optim.method="subplex", nstarts=10)

#get EP of the leaves of the phylogenetic tree ($tip.states)
#from $tip.states, we select the column with the highest TP and retrieve its associated genotypes
print("working on transition probabilities of leaves...")
colnames(class2_corhmm$tip.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
TP_geno_corHMM <- data.frame(apply(class2_corhmm$tip.states,1,max))
etat_feuilles <- data.frame(colnames(class2_corhmm$tip.states)[apply(class2_corhmm$tip.states,1,which.max)])

#since we do not know which of G0 and G1 is R1 or R2, we do both condition and select the best one
#associate G0 = R1 and G1= R2
genotype_feuille_G0R1 <- data.frame(apply(etat_feuilles, 1, function(r) any(r %in% c("(0,R2)", "(1,R2)")))*1)
#associate G0 = R2 and G1 = R1
genotype_feuille_G0R2 <- data.frame(apply(etat_feuilles, 1, function(r) any(r %in% c("(0,R1)", "(1,R1)")))*1)

#convert phenotype to 0 and 1
phenotype_feuille <- data.frame(apply(etat_feuilles, 1, function(r) any(r %in% c("(1,R1)", "(1,R2)")))*1)

#keep same order of tips/leaves 
tiplabel <- data.frame(as.numeric(arbre$tip.label))
simulation_df<-data.frame(tiplabel,etat_feuilles,TP_geno_corHMM,genotype_feuille_G0R1,genotype_feuille_G0R2, phenotype_feuille)
colnames(simulation_df) <- c("tip_order", "leaf_state", "Transition_prob", "leaf_genotype_G0R1", "leaf_genotype_G0R2","leaf_phenotype")
simulation_df<-simulation_df[order(simulation_df$tip_order),]

#get EP of all the internal noeuds of the phylogenetic tree
#from $states, we select the column with the highest prob and retrieve its associated genotypes
#since $state represents the internal nodes we have k-1 states, 1 being the root with prob 1/4 at each geno-pheno
print("working on probabilities of internal nodes...")
etat_noeuds <- data.frame(colnames(class2_corhmm$states)[apply(class2_corhmm$states,1,which.max)])

genotype_noeud_G0R1 <- data.frame(apply(etat_noeuds, 1, function(r) any(r %in% c("(0,R2)", "(1,R2)")))*1)
genotype_noeud_G0R2 <- data.frame(apply(etat_noeuds, 1, function(r) any(r %in% c("(0,R1)", "(1,R1)")))*1)

#we fuse the column genotype_coala and genotype_corhmm in a df to count errors
#count number of unmatched genotype between coala and corHMM
print("comparing genotypes of coala and corHMM...")
fused_df <- data.frame(geno_pheno$haplotype,simulation_df$leaf_genotype_G0R1,simulation_df$leaf_genotype_G0R2, simulation_df$Transition_prob)
fused_df[,c('match_G0R1')] <- apply(fused_df[,c(1,2)], 1, sum)
fused_df[,c('match_G0R2')] <- apply(fused_df[,c(1,3)], 1, sum)

#we count number of errors for both conditions and select the one with minimum error
print("calculating errors...")
count_error_G0R1 <- nrow(fused_df[fused_df$match_G0R1 == 1,])
count_error_G0R2 <- nrow(fused_df[fused_df$match_G0R2 == 1,])

if (count_error_G0R1 < count_error_G0R2) {
  count_error=count_error_G0R1
  compare_df <- data.frame(simulation_df$tip_order,geno_pheno$haplotype,simulation_df$leaf_genotype_G0R1,simulation_df$Transition_prob,fused_df$match_G0R1)
  node_df<-data.frame("id"=1:nrow(genotype_noeud_G0R1),etat_noeuds,genotype_noeud_G0R1)
} else {
  count_error=count_error_G0R2
  compare_df <- data.frame(simulation_df$tip_order,geno_pheno$haplotype,simulation_df$leaf_genotype_G0R2,simulation_df$Transition_prob,fused_df$match_G0R2)
  node_df<-data.frame("id"=1:nrow(genotype_noeud_G0R2),etat_noeuds,genotype_noeud_G0R2)
}

colnames(compare_df) <- c("tip_order","coala_genotype", "corHMM_genotype", "Transition_prob", "match")
colnames(node_df) <- c("id", "node_state","node_genotype")

percent_error <- (count_error/k)*100
unmatched_probability <- (sum(compare_df$Transition_prob[compare_df$match==1]))/k

#define file_name for class 2 results
file_name <- as.character(args[10])

#save phylogeneotic tree as pdf
print("Drawing phylogenetic tree...")
if(T)
{
  #pdf(file.path(paste0("simulation_", sim),paste0("class_2"),file = paste0("tree_",file_name,".pdf")))
  ## tree annotation
  nb_noeud<-arbre$Nnode+k # (k-1)+k = 2k-1
  couleur_branche<-rep("blue",time=nb_noeud)
  
  tree_CC<-ggtree(tr=arbre,color=couleur_branche)+theme_tree2()
  #tree_corhmm<-ggtree(tr=arbre)
  tree_CC<-tree_CC %<+% geno_pheno # contains phenotype and genotype of coala
  tree_CC <- tree_CC %<+% compare_df #contains phenotype and genotype of corHMM
  tree_CC <- tree_CC %<+% node_df #contains internal node
  tree_CC<-tree_CC + geom_tippoint(mapping=aes(color=factor(haplotype),x=x+.005),size=4)
  tree_CC <- tree_CC + geom_tippoint(mapping=aes(color=as.factor(corHMM_genotype),x=x+.04),size=4, shape=1)
  tree_CC<-tree_CC + geom_tippoint(aes(color=factor(phenotype),x=x+.08),size=4,shape=15)
  #tree_CC <- tree_CC + geom_tippoint(aes(color=as.factor(leaf_phenotype),x=x+.08),size=4,shape=15)
  tree_CC <- tree_CC + geom_nodepoint(aes(color=as.factor(node_genotype)), size=1.5)
  plot(tree_CC)
  dev.off()
}
#save all results in my_list
my_results_class2 <- list(seg_sites=snp, arbre_coala=arbre, emission_probability=df, coala_simulation = geno_pheno, coala_phenotype= phenotype_df, 
                   corHMM_class2 = class2_corhmm, transmission_prob_leaf = simulation_df, node_probability= node_df, compare_table = compare_df, final_tree=tree_CC,
                   count_error=count_error, percentage_error=percent_error, unmatched_probability=unmatched_probability)

#end time of script
print("end time...") 
my_results_class2$time <- proc.time() - ptm

#save all results in .RData so that the results can be easily loaded in R to be analysed
print("saving results for corHMM class2...")
save(my_results_class2,file=file.path(paste0("snp_plus1"),paste0("simulation_", sim),paste0("class_2"),paste0("res_",file_name,".RData")))

print("End of script")
