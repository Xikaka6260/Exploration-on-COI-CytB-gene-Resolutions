#### Assignment 2 ------ Version 7 ----- BOLD + NCBI 

### Zizhen Zhong 

### OCT 24, 2024

#Packages Used ---- 
library(tidyverse)
library(Biostrings)
library(rentrez)
library(dplyr)
library(DECIPHER)
library(ape)
library(phytools)
library(dendextend)
library(gplots)


#### ---- Part1 a : Obtain Dataset (COI gene) ----

#Obtained data from BOLD system on OCT 14, 2024 
BOLD.data <- read_tsv("../data/Canidae_data_BOLD.txt")

#To get an overview about this dataset 
summary(BOLD.data)
names(BOLD.data)
dim(BOLD.data)
#==> The dimension of this original data frame is 2886 by 80 variables 
unique(BOLD.data$markercode)
#==> Various markercodes data are available, however, I am only interested in COI-5P gene. 

#Minimize dataframe size to obtain only columns wanted, and then polish the data to obtain rows with only COI gene sequences available. Columns of latitude and longitude is not being used at the moment, so I will not filter and remove the missing values from these columns yet. 
COI.BOLD <- BOLD.data %>% 
  select(processid, genus_name, species_name, markercode, nucleotides, lat, lon) %>% 
  filter(!is.na(markercode)) %>% 
  filter(!is.na(genus_name)) %>% 
  filter(!is.na(species_name)) %>%
  filter(markercode == "COI-5P")

# Filtration check : 
unique(COI.BOLD$markercode)
sum(is.na(COI.BOLD$genus_name))
sum(is.na(COI.BOLD$species_name))
dim(COI.BOLD)
#==> dimension of this new data frame is now 2309 by 7 variables, filtering done successfully only samples with COI-5P markercode remained and all rows with missing values under genus and species name columns are removed. 

#Get an idea of what are the available genus
unique(COI.BOLD$genus_name)
#==> There are 11 unique genus groups available. 

#I filtered the dataframe into a smaller dimension, because I am only interested in species of fox. These are 6 genus groups are wolf and dog species representations I will be filtering out. 
df_Fox_COI <- COI.BOLD %>% 
  filter(!genus_name == "Canis") %>% 
  filter(!genus_name == "Nyctereutes") %>% 
  filter(!genus_name == "Lycaon") %>% 
  filter(!genus_name == "Speothos") %>% 
  filter(!genus_name == "Chrysocyon") %>% 
  filter(!genus_name == "Cuon") %>% 
  filter(!is.na(nucleotides))


#Filter Checks : 
dim(df_Fox_COI)
#==> Dimension of this fox data frame is 469 rows by 7 variables, filtration done successfully. 

unique(df_Fox_COI$genus_name)
unique(df_Fox_COI$species_name)
#==> There are 5 Genus of foxes and 18 species of foxes.The other genus groups had been successfully removed  

#Get an overview of the sequence structure 
head(df_Fox_COI$nucleotides)
#==> Some sequences contain a lot of "N" at the end of the sequence 


#Re-filter the data frame and set more parameters. Created a new sequence column that has removed all "N" and "-" at the beginning and end of the sequence. Also, remove sequences that have invalid bases other than "ATGC", and remove sequences that contains "Y" bases. Lastly, I set some parameters to remove poor quality sequences. The threshold for gap and missing value tolerance is 1% prior to sequence alignment 
df_Fox_COI <- df_Fox_COI %>% 
  mutate(nucleotides2 = str_remove_all(nucleotides, "^N+|N+$")) %>% 
  filter(nucleotides2 == str_remove_all(nucleotides2, "^-+|-+$")) %>% 
  filter(!nucleotides2 == str_remove_all(nucleotides2, "[ATGC]")) %>% 
  filter(nucleotides2 == str_remove_all(nucleotides2, "Y+")) %>% 
  filter(str_count(nucleotides2, "-") <= (0.01 * str_count(nucleotides2))) %>% 
  filter(str_count(nucleotides2, "N") <= (0.01 * str_count(nucleotides2)))

#Filter Checks: 
dim(df_Fox_COI)
view(df_Fox_COI)
#==> dimension is 440 rows by 8 variables, filtering happened successfully

#Confirm gaps and missing datas have been removed 
sum(str_count(df_Fox_COI$nucleotides2, "-"))
sum(str_count(df_Fox_COI$nucleotides2, "N"))
sum(str_count(df_Fox_COI$nucleotides2, "Y"))
#==> There are 0 gaps "-" or "Y" and 4 "N" remained after filtering.


#Check distribution of sequence length 
hist(nchar(df_Fox_COI$nucleotides2), main = "Distribution of Fox Sequence Length (COI) ", xlab = "Length of Basepairs per Sequence", ylab = "Number of Species", col = "gold")
 #==> Base on the histogram distribution, there are many data with very short nucleotide length, and majority of the data have length of 1500bps. 


summary(nchar(df_Fox_COI$nucleotides2))
#==>Some sequences are quite short compare to the average, where the minimum length is 178bp, first quantile is 657bp in length, and the median is 1542bp. 

table(str_count(df_Fox_COI$nucleotides2))
#==> Based on the sequence count table, length <= 310 are very different from the rest of the data. Since some of the basepair length are so dissimilar compare to the majority of the data, these short sequences need to be removed to enhance the data quality.

#Made a new dataframe that stores the short sequences samples 
df_short_sequence <- df_Fox_COI %>% 
  filter(nchar(nucleotides2) <= 310)

dim(df_short_sequence)
#==> Base on the dimension I know there are 37 data samples with shorter than 310bp in sequence length

#checking what are the different species with short sequence length 
unique(df_short_sequence$species_name)
table(df_Fox_COI$species_name)
#==> Base on the result, it seems like all of the Lycalopex genus have only a short sequence length in the BOLD data system. So I decided to filter out the Lycalopex genus and the other 2 species prior to the alignment. This is because the sequence length difference is too big compare to other species samples. It will cause a lot gaps introduced during multiple sequence alignment, and results will be less reliable. Also, based on the table, I know that there are 16 samples of Cerdocyon thous and 54 samples of Vulpes vulpes samples; thus removing the short data sequence samples from these 2 species will not affect my analysis. 

#filter out short sequences 
df_Fox_COI <- df_Fox_COI %>% 
  filter(!nchar(nucleotides2) <= 310)

dim(df_Fox_COI)
#==> dimension of this fox data frame is now 403 by 8, filter has been done 

#visualize again to see the distribution of sequence length 
hist(nchar(df_Fox_COI$nucleotides2), main = "Distribution of Fox Sequence Length (COI) ", xlab = "Length of Basepairs per Sequence", ylab = "Number of Species", col = "gold")
#==> The sequence length differences is less diverse and confirmed the removal of sequence length <= 310bps. 

#To check again in terms of numerical sequence distribution 
summary(nchar(df_Fox_COI$nucleotides2))
table(nchar(df_Fox_COI$nucleotides2))
#==> The distribution of sequence length is less varied and the shortest sequence is now 516bp in length.


unique(df_Fox_COI$species_name)
unique(df_Fox_COI$genus_name)
#==> There are 5 genuses with a total of 18 species left for the following alignment step 


rm(df_short_sequence)



#### ---- Part1 b : Random Sample Selection Per Species (COI) ----

#Prior to performing multiple sequence alignment, and constructing a phylogenetic tree, I randomly selected a sample from each species as a representation. 
set.seed(123)
selected_COIsample <- df_Fox_COI %>%
  group_by(species_name) %>% 
  sample_n(1) %>% 
  ungroup()

#Quality Check : 
dim(selected_COIsample)
unique(selected_COIsample$species_name)
summary(str_count(selected_COIsample$nucleotides2))
#==> The dimension of the data frame has shown that there are 18 samples selected, and are one from each species. The mean sequence length is 910bp, with minimum being 551bp and maximum being 1542bp.

#setting the 18 sampled species names into a new vector for downstream analysis
BOLD.names <- selected_COIsample$species_name



#### ---- Part1 c : Multiple Sequence Alignment (COI) ----

#Assigning corresponding species names to the sequences 
names(selected_COIsample$nucleotides2) <- selected_COIsample$species_name

#Set the nucleotide sequences into DNAStringSet format for downstream sequence alignment and then re-assign this DNAStringSet format back
selected_COIsample <- DNAStringSet(selected_COIsample$nucleotides2)

#Check if it has successfully changed into Stringset 
class(selected_COIsample)

#To perform multiple sequence alignment of those 18 species using muscle package, with default iterations, gap opening and gap extension penalty of 15 and 6.66. These parameters were selected based on this research (https://doi.org/10.1371/journal.pone.0014156). 
COI_alignment <- DNAStringSet(muscle::muscle(selected_COIsample, gapopening = 15, gapextend = 6.66, use.names = TRUE))

BrowseSeqs(COI_alignment)
#==> Majority of gaps were introduced at the beginning and ends of the sequences since some of the sequences have greater length than others, but the overall alignment looks acceptable.  



#### ---- Part1 d : Phylogenetic Tree Construction (COI) ----

#this code helps to turn sequences in DNAStringSet format 
COIbin <- as.DNAbin(COI_alignment)

class(COIbin) 
#Confirm success in translating to DNAbin format 

# For model, the reason why I picked TN93 is because the model accounts for different base frequencies and transition/transversion biases. These parameters are very important, because it is now commonly known that same nucleotide base group more structurally similar. Thus it is easier to occur transition than transversion. Also, a research done by Bohlin & Pettersson in 2019 supported the idea that in eukaryotic genome, the regions with more recombination events tend to show higher GC content (doi: 10.1016/j.csbj.2019.03.001)

#GTR model was ruled out although it has similar parameters, because most recent research has questioned the validity of the model (https://doi.org/10.1093/bioinformatics/btk001)
distanceMatrix_COI <- dist.dna(COIbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

head(distanceMatrix_COI)
#==> Based on the distance matrix tables, I can see that distance matrix is 0 when both axis are the same specie, this confirms that the matrix has performed successfully. I noticed species in the same genus usually exhibit distance below 0.02, with one exception which is Lycalopex sechurae. The distance matrix for Lycalopex sechurae with any species under the same genus is relative high, at the range of 0.5 - 0.7. 


#Better visualize the distance matrix, for quality check 
heatmap.2(distanceMatrix_COI, key = TRUE, key.title = "Dissimilarity", main = "Distance Matrix (COI)", cexRow = 0.5, cexCol = 0.5)
#==> Species within the same genus appear redder on the map compared to those of different genus. This indicates that these species are more genetically similar based on COI gene. Whereas Vulpes macrotis stands out as particularly distinct from Lycalopex griseus and Lycalopex vetulus. Similarly, Vulpes chama shows a significant higher distance value from Otocyon meglotis. This suggests that, based on the COI gene, the genetic composition of Vulpes chama and Otocyon meglotis is quite distinct, indicating they may be more distantly related. Overall, the distance matrix generated using the TN93 model seem to have performed as expected, clustering species from the same genus closer together. 


#Clustering method I prefer is neighbor joining, because it constructs unrooted phylogenetic tree. It is a good fit for my project because it only consider how species are related in terms of similarity in sequences without having to specify ancestor root. The cutoff value for a specie to be put into one cluster is 0.02. This cutoff value was decided based on the distance matrix table. 
clusters.COI <- DECIPHER::TreeLine(myDistMatrix = distanceMatrix_COI,
                                   method = "NJ",
                                   cutoff = 0.02,
                                   showPlot = TRUE,
                                   type = "dendrogram",
                                   verbose = TRUE)


#===> Base on the phylogeny tree, most of the species are arranged in accordance to the their genus, however with 3 exceptions. Which is Vulpes marcrotis, Otocyon megalotis, and Lycalopex sechurae. 

#Checking back at the alignment 
BrowseSeqs(COI_alignment)
#==> Vulpes marcrotis obtains 585bp, Otocyon megalotis obtains 588bp, Lycalopex sechurae 582bp. These 3 sequences all have relatively short sequence in length, compare to the rest samples selected. However, they are not the shortesting sequence. Thus, I think the reason why these 3 species were not grouped accurately. 



#### ---- Part2 a : Obtain Dataset (CytB gene) ----

#Here is getting an overall idea of how many data are in NCBI nuccore database with these 5 genus of Fox I am working with and also contain cytochrome B gene.
entrez_search( db = "nuccore", term = "Vulpes[ORGN] OR Lycalopex[ORGN] OR Urocyon[ORGN] or Cerdocyon[ORGN] or Otocyon[ORGN] AND CytB[Gene]",  retmax = 100)
#==> There are about 1162 data samples in nuccore database

#I went to the NCBI webpage and searched for the similar terms, and noticed majority of the sequence length is between 300 - 880bps.
cytB_search <- entrez_search( db = "nuccore", term = "Vulpes[ORGN] OR Lycalopex[ORGN] OR Urocyon[ORGN] or Cerdocyon[ORGN] or Otocyon[ORGN] AND CytB[Gene] AND 300:880 [SLEN]", retmax = 300, use_history = TRUE)

cytB_search$count
#==> There are 414 samples in the nuccore database that meets my criteria but we only obtained 300 IDs from the database. 

#I want to know how many unique species I have obtained from the search. Here I set a new vector to store all the unique IDs from the search
ID.cytB <- unique(cytB_search$ids)

class(ID.cytB)
#This is a character vector that contains all 145 ID numbers from the search. 

#Since the vector is a super long character vector, I want to turn every ID as a unique value, thus set this vector into list format. 
ID.cytB <- as.list(ID.cytB)

class(ID.cytB)
#This is now a list vector 

#In order to get organism names from each IDs, first I created a new empty character vector to store all of the organisms names after my loop 
organism.names <- vector("character", length(ID.cytB))

#This loop function is used to put every single fox IDs I obtained from the previous step and check for its corresponding organism name. Here I am doing for all 145 IDs I have obtained 
for (i in seq_along(ID.cytB)) {
  check.ID <- entrez_summary(db = "nuccore", id = ID.cytB[i])
  organism.names[i] <- check.ID$organism
}


unique(organism.names)
#==> There are 17 unique species obtained from the search. 

#assigning these 17 species into a new vector for downstream analysis 
NCBI_names <- unique(organism.names)


#### ---- Part2 b : Fetching Dataset (CytB gene) ----

#Fetching nucleotide sequence data from NCBI, for each of the IDs obtained previously. 
#cytB_fetch <- entrez_fetch(db = "nuccore", id = ID.cytB, rettype = "fasta")

#view(cytB_fetch)
#class(cytB_fetch)
#head(cytB_fetch)
#Noticing this is a super long character vector. 

#Saved this fasta file into working directory, every sequence is separated by \n\n. 
#write(cytB_fetch, "CytB_fetch.fasta", sep = "\n\n")




#### ---- Part2 c : Upload and Filter fetched (CytB gene) ----

#Read and load the fetched file to countable stringsets 
cytBstringset <- readDNAStringSet("CytB_fetch.fasta")

class(cytBstringset)
#the fasta file has been successfully read in and transformed into biostring mode

# What are the names for each corresponding sequences 
names(cytBstringset)
#==> Noticed the names are too long, and included many unrelated information to my project

#Create a dataframe to store the names and sequences into different columns 
df_cytB <- data.frame(CytB_names = names(cytBstringset), CytBsequence = paste(cytBstringset))

dim(df_cytB)
#==> There are 300 rows by 2 columns 

#Create a new column and input names obtaining only the second and third words from the original CytB_names columns 
df_cytB$species_names <- word(df_cytB$CytB_names, 2L, 3L)

dim(df_cytB)
view(df_cytB)
#==> There are 300 rows by 3 columns, the new column has been created successfully

#Quality check
unique(df_cytB$species_names)
#==> The species names are looking perfect with only containing words needed, no need further filtering. 

#Comparing the species names searched and species names from fetched 
setdiff(unique(df_cytB$species_names), NCBI_names)
#==> The different species names are Pseudalopex genus group and Dusicyon genus. Based on literature, the Pseudalopex genus is another name for Lycalopex (https://www.researchgate.net/publication/307512612_Lycalopex_culpaeus); Dusicyon is another name for Cerdocyon (https://www.researchgate.net/publication/283815197_Notes_on_crab-eating_fox_Dusicyon_thous_seed_dispersal_and_food_habits_in_southeastern_Brazil). 

#Editing the species name to another synonyms name
df_cytB$species_names[df_cytB$species_names == "Pseudalopex sechurae" ] <- "Lycalopex sechurae"
df_cytB$species_names[df_cytB$species_names == "Pseudalopex griseus" ] <- "Lycalopex griseus"
df_cytB$species_names[df_cytB$species_names == "Pseudalopex vetulus" ] <- "Lycalopex vetulus"
df_cytB$species_names[df_cytB$species_names == "Pseudalopex gymnocercus" ] <- "Lycalopex gymnocercus"
df_cytB$species_names[df_cytB$species_names == "Pseudalopex culpaeus" ] <- "Lycalopex culpaeus"
df_cytB$species_names[df_cytB$species_names == "Dusicyon thous" ] <- "Cerdocyon thous"

#renaming check 
setdiff(unique(df_cytB$species_names), NCBI_names)
#==> There are 0 different species names between the fetched data and the searched data, thus it also indicated that the renaming process was done successfully. 

unique(df_cytB$species_names)
#==> There are 16 unique species with CytB gene sequences


#To confirm if the fetched sequences have similar range of sequence length as the search parameter 300:880bp. 
hist(str_count(df_cytB$CytBsequence), main = "Fox Species sequence Length Fetched (CytB)", xlab = "Basepair Length", ylab = "Species Count", col = "pink")
#==> The histogram looks skewed to both ends, where majority of the sequence length are either 400 bp long or 880bp long. 




#### ---- Part2 d : Random Sample Selection per Species (CytB gene) ----

#Prior to performing multiple sequence alignment, and constructing a phylogenetic tree, I randomly selected a sample from each species as a representation. 
set.seed(432)
selected_cytBsample <- df_cytB %>%
  group_by(species_names) %>% 
  sample_n(1) %>% 
  ungroup()

#Quality Check : 
dim(selected_cytBsample)
unique(selected_cytBsample$species_names)
summary(str_count(selected_cytBsample$CytBsequence))
#==> The dimension of the data frame has shown that there are 16 samples selected, and are one from each species. The mean sequence length is 373bp, with minimum being 317bp and maximum being 476bp.




#### ---- Part2 e : Multiple Sequence Alignment (CytB) ----

#Assigning corresponding species names to the sequences 
names(selected_cytBsample$CytBsequence) <- selected_cytBsample$species_names

#Set the nucleotide sequences into DNAStringSet format for downstream sequence alignment and then re-assign this DNAStringSet format back
selected_cytBsample <- DNAStringSet(selected_cytBsample$CytBsequence)

#Check if it has successfully changed into Stringset 
class(selected_cytBsample)

#Perform multiple sequence alignment of the 16 species with same parameters as COI gene alignment. 
cytB_alignment <- DNAStringSet(muscle::muscle(selected_cytBsample, gapopening = 15, gapextend = 6.66, use.names = TRUE))

BrowseSeqs(cytB_alignment)
#==> Notice the gaps introduced in the CytB gene alignment was way less than COI gene alignment. This is mostly because the range of sequence length in CytB gene sample is less diversed than COI gene. 




#### ---- Part2 f : Phylogenetic Tree Construction (CytB gene) ----

#this code helps to turn sequences in DNAStringSet format 
cytBbin <- as.DNAbin(cytB_alignment)

class(cytBbin) 
#Confirm success in translating to DNAbin format 

#I used the same model and same parameter done for COI gene for cytB gene 
distanceMatrix_cytB <- dist.dna(cytBbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

head(distanceMatrix_cytB)
#==> The values in distance matrix in CytB gene between these species are way larger than COI gene. However, comparing the distance values between same genus, the values are still relatively small similar to the COI gene distance. Thus the cutoff values 0.02 for generating the phylogeny tree can also be use for CytB gene. 


#quick check for the matrix values 
heatmap.2(distanceMatrix_cytB, key = TRUE, key.title = "Dissimilarity", main = "Distance Matrix (CytB)", cexRow = 0.5, cexCol = 0.5)
#==> Otocyon megalotis shows a more distinct dissimilarity in CytB gene composition with Lycalopex sechurae and Vulpes macrotis. Lycalopex sechurae also shows a distinct dissimilarity with Urocyon cinereoargenteus. Overall, the distance matrix generated looks reasonable where species with the same genus name are more genetically similar. 


#constructing the phylogeny tree with the same parameters as COI gene.  
clusters.cytB <- DECIPHER::TreeLine(myDistMatrix = distanceMatrix_cytB,
                                   method = "NJ",
                                   cutoff = 0.02,
                                   showPlot = TRUE,
                                   type = "dendrogram",
                                   verbose = TRUE, processors = 1)


#===> Noticing the species from the same genus are placed close together, with Vulpes macrotis and Lycalopex sechurae being the exceptions. These are the same 2 exceptions observed from COI gene dendrogram.  




#### ---- Part3 : Comparing the Trees ---- 

#To see what are the different names that 2 gene samples dont have in common 
setdiff(BOLD.names, NCBI_names)
#==> Lycalopex fulvipes, Urocyon littoralis, Vulpes chama, and Vulpes velox are species only obtained from BOLD database. 

setdiff(NCBI_names, BOLD.names)
#==> Vulpes rueppellii and Vulpes pallida are species only obtained from NCBI database. Vulpes macrotis zinseri is the same as as Vulpes macrotis, it does exist in BOLD sample. 

#==> Thus in order to interpret the dendrogram more accurately, I would start by disregarding the placement and cluster group of these uncommon species


#This code checks if the tree I have constructed previously is in phylo format, if it not in phylo format, then change it to phylo format. 
if (!inherits(clusters.cytB, "phylo")) {
  clusters.cytB <- as.phylo(clusters.cytB)
}

if (!inherits(clusters.COI, "phylo")) {
  clusters.COI <- as.phylo(clusters.COI)
}

#Here I put both the COI phylogeny tree and CytB phylogeny tree side by side 
par(mfrow = c(1, 2))
plotTree(clusters.cytB, main = "Tree 1 : CytB gene" )
plotTree(clusters.COI, main = "Tree 2 : COI gene", direction = "leftwards")
#==> Despite using identical methods, parameters, models, and cutoff values to analyze sequences from both genes, the resulting dendrograms have quite different resolutions. The COI dendrogram (right) has Cerdocyon thous grouped closer to the Lycalopex genus group; In CytB dendrogram (left) Cerdocyon thous is grouped far away from the Lycalopex group. 

# ==> Also, noticing COI dendrogram has all Vulpes genus group close together, whereas CytB dendrogram has Otocyon megalotis separating the Vulpes group. But in both cases, Otocyon megalotis are clustered close to Vulpes lagopus and Vulpes vulpes.  

#==> Lastly, the outgroup for both dendrograms have Vulpes macrotis and Lycalopex sechurae as previously mentioned. 

