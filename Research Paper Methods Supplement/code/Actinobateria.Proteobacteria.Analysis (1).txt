setwd("c:/school/coding/Bio Lab")

#Read in rarefied table
otu_table <- read.table('hmp_rarefied_table.txt', 
                        comment="",
                        header=TRUE,
                        sep="\t",
                        skip=1,
                        as.is=TRUE,
                        check.names=F,
                        row=1)
#Read in metadata
metadata <- read.table("HMP_mapping_file.txt",
                       comment="",
                       header=TRUE,
                       sep="\t",
                       check.names=F,
                       row=1)

# Subsetting
metadata <- metadata[metadata$BODY_SITE=="UBERON:feces" ,]

#metadata <- metadata[metadata$CHRONIC_CONDITION == "y",]
#metadata <- metadata[metadata$CHRONIC_CONDITION == "n",]

IDs_Keep <- intersect(colnames(otu_table), rownames(metadata))
metadata2 <- metadata[IDs_Keep,]
otu_table <- otu_table[IDs_Keep]

# Add a taxa column to the OTU table
otu_table$OTU_ID <- row.names(otu_table)

#Read in the taxa table
taxa_table <- read.table('97_otu_taxonomy.txt', sep = "\t")

colnames(taxa_table) <- c("OTU_ID","Taxa")

#Merge the two tables
otu_taxa <- merge(otu_table, taxa_table)

#Keep only the entries in the OTU table
rownames(otu_taxa) <- otu_taxa$OTU_ID

#Removing the OTU_ID column. No idea why, just being told to do it.
Rows <- dim(otu_taxa)[2]
drop_IDs <- c(2:(Rows-1),Rows)

otu_taxa2 <- otu_taxa[,drop_IDs]



#Output the table
#write.table(otu_taxa2,"hmp_taxa_otu.tsv", sep = "\t", row.names = F)

# taxa <- read.table('hmp_taxa_otu.tsv',
#                    comment="",
#                    header=TRUE,
#                    sep="\t",
#                    skip=F,
#                    as.is=TRUE,
#                    check.names=F)

level <- 2

#List of names
names_split <- array(dim = c(length(otu_taxa2$Taxa),level))

taxa_names <- as.character(otu_taxa2$Taxa)


for (i in 1:length(taxa_names)){
  names_split[i,] <- head(strsplit(taxa_names[i], ";", fixed = T)[[1]], n = level)
}
otu_taxa2$SampleID <- rownames(otu_taxa2)


taxa_names <- apply(names_split, 1, function(x) paste(x[1:level],sep ="", collapse = ";"))

otu_taxa2$taxonomy <- taxa_names

#We actually added 3 columns, not 1, so we need to subtract 3 here
sample_num <- ncol(otu_taxa2) -3 
otu_taxa2$taxonomy

otu_taxa3 <- aggregate(otu_taxa2[,1:sample_num], by = list(otu_taxa2$taxonomy), FUN = sum)

#Sets the first column to taxonomy, since this is what we aggreated by
names(otu_taxa3)[1] <- "taxonomy"
nrow(otu_taxa3)

rownames(otu_taxa3) <- otu_taxa3$taxonomy

#Removes the taxonomy column
otu_taxa3 <- otu_taxa3[,!names(otu_taxa3) == "taxonomy"]

otu_taxa3[otu_taxa3 < sum(colSums(otu_taxa3))/1000000] <- 0

otu_taxa3[otu_taxa3 < 2] <- 0

otu_taxa3 <- otu_taxa3[rowSums(otu_taxa3 > 0) > (0.05*ncol(otu_taxa3)),]


#Calculating relative abundance
for(i in 1:ncol(otu_taxa3)){
  otu_taxa3[,i] <- otu_taxa3[,i]/sum(otu_taxa3[,i])
}

#formatting for plotting
otu_taxa3 <- data.frame(t(otu_taxa3))

otu_taxa3$SampleID <- rownames(otu_taxa3)



library(reshape2)
library(plyr)
library(ggplot2)
otu_taxa3 <- melt(otu_taxa3, id.vars = "SampleID",variable.name = "Taxa", value.name = "RelativeAbundance")
metadata2$SampleID <- rownames(metadata2)
otu_taxa3 <- merge(otu_taxa3, metadata2, by = "SampleID")

ggplot(otu_taxa3, aes(x = CHRONIC_CONDITION, y = RelativeAbundance, fill = Taxa)) + geom_bar(stat = "identity", position = "fill") +
  labs(x= "Chronic Condition")

taxa_list <-"k__Bacteria..p__Actinobacteria"
filtered <- subset(otu_taxa3, is.element(otu_taxa3$Taxa,taxa_list))

color <- c("cadetblue1","pink")
ggplot(filtered, aes(x = CHRONIC_CONDITION, y = RelativeAbundance, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(y = "Relative Abundance", x = "Chronic Condition")+
  scale_fill_manual(labels = c("Proteobacteria", "Actinobacteria"), values = color) +
  scale_x_discrete(labels = c("n","y"))




shapiro.test(as.numeric(filtered$CHRONIC_CONDITION=='n'))
wilcox.test(as.numeric(filtered$CHRONIC_CONDITION=='n'), as.numeric(filtered$CHRONIC_CONDITION=='y'), na.rm=TRUE)
