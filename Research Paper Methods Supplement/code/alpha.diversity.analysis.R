#This is the working directory
setwd("c:/school/coding/Bio Lab")

#Rarefied table. We had an issue with this before. 
otu_table <- read.table("table_from_HMP.txt", 
                        comment="",
                        header=TRUE,
                        sep="\t",
                        skip=1,
                        as.is=TRUE,
                        check.names=F,
                        row=1)

metadata <- read.table("HMP_mapping_file.txt", 
                       sep = "\t", 
                       comment="", 
                       header=T,
                       check.names=F, 
                       row=1)


#Converted from qza -> biom -> tsv
alpha <- as.data.frame(read.table("alpha-diversity.txt",
                                  sep='\t',
                                  header=TRUE,
                                  as.is=TRUE,
                                  check.names=FALSE,
                                  row=1))

#Section where we subswt to data to remove unwanted categories. Also makes the dataset smaller so our laptops can run it.
metadata <- metadata[metadata$BODY_SITE=="UBERON:feces" ,]
metadata.no.cc <- metadata[metadata$CHRONIC_CONDITION=="n", ]
metadata.cc <-metadata[metadata$CHRONIC_CONDITION=="y", ]
metadata

#Our OTU table doesn't have a taxonomy column, so we don't remove the last column


IDs_Keep_cc <- intersect(colnames(otu_table), rownames(metadata.cc))
IDs_Keep_no_cc <- intersect(colnames(otu_table), rownames(metadata.no.cc))

# Now let's filter the metadata to keep only those samples
metadata2.cc <- metadata[IDs_Keep_cc,]
metadata2.no.cc <- metadata[IDs_Keep_no_cc,]



# Alpha diversity has the samples as row names
alpha.cc <- alpha[IDs_Keep_cc, ]
alpha.no.cc <- alpha[IDs_Keep_no_cc, ]

#Add an alpha column to the metadata
#metadata.alpha <- as.data.frame(metadata2$alpha)

#metadata2.vf <- metadata2[metadata2$BODY_SITE=="UBERON:vaginal fornix" ,]

#no.cc.x <- rownames(metadata2.vf$CHRONIC_CONDITION=='n')
#alpha.no.cc <- alpha[no.cc.x, ]
hist(alpha.no.cc, xlab="Alpha Diversity", main='NCC')
hist(alpha.cc, xlab="Alpha Diversity", main='CC')


#hist(metadata.alpha.no.cc$alpha, xlab="Alpha Diversity", main='Females')

shapiro.test(alpha.no.cc) #p-value = 4.352e-05
shapiro.test(alpha.cc) #p-value = 0.1853
wilcox.test(alpha.no.cc, alpha.cc, na.rm=TRUE) #p-value = 0.5421

## PLOTTING

# loading needed package
library(ggplot2)

x_labels <- c('Gut')

# plotting the data
# geom_boxplot plots the shannon data in boxplots, separating the values by body site
# scale_fill_manuel fills the boxplots with colors representing chronic condition status
# labs adds labels and a caption (our names)
ggplot() +
  geom_boxplot(data=metadata, aes(x= BODY_SITE, y= alpha, fill= CHRONIC_CONDITION)) +
  scale_fill_manual(values= c("cadetblue3", "white"), labels = c("No Chronic Condition", "Chronic Condition"), name = "") + labs(x = "Chronic Condition Status", y = "Shannon's Diversity Index") +
  scale_x_discrete(labels = x_labels) + theme(legend.position="bottom")
