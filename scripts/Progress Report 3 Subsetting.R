#This is the working directory
setwd("/Users/keliansoh/Comp Bio/")

#What the data directory is - this is referenced to get files
data_dir <- "Comp Bio/"

#Get path to OTU table and map
otu_fp <- paste(data_dir, "table_from_HMP.txt.zmdownload", sep='')
map_fp <- paste(data_dir, "HMP_mapping_file.txt", sep='')

#Load map 
metadata <- read.table("HMP_mapping_file.txt", 
                       sep = "\t", 
                       comment="", 
                       header=T,
                       check.names=F, 
                       row=1)

#The hypothesis: Antibiotic use is associated with a significant reduction in microbiome diversity across various body sites. A significant reduction in alpha diversity will be observed in male patients undergoing an antibiotic treatment; obesity is expected to correlate with additional reductions in male patients.

#What groups are we interested in:
#Note: May decide to go the more general route and focus on the supersite which has less variables.
levels(metadata$BODY_SITE) #UBERON:tongue, UBERON:fossa, UBERON:skin, UBERON:nostril, UBERON:mouth, UBERON:vaginal fornix, UBERON:feces, UBERON:oropharynx, UBERON:vagina
levels(metadata$SEX) #male, female
levels(metadata$HMPBODYSUPERSITE) #Oral, Skin, Airways, Urogenital_tract, Gastrointestinal_tract
levels(metadata$OBESITY) #n, y, None
levels(metadata$CHRONIC_CONDITION) #n, y, None

#Keep only the body sites we're interested in:
metadata_sub <- metadata[metadata$BODY_SITE=="UBERON:tongue" | metadata$BODY_SITE=="UBERON:fossa" | metadata$BODY_SITE=="UBERON:skin" | metadata$BODY_SITE=="UBERON:nostril" | metadata$BODY_SITE=="UBERON:mouth" | metadata$BODY_SITE=="UBERON:vaginal fornix" | metadata$BODY_SITE=="UBERON:feces" | metadata$BODY_SITE=="UBERON:oropharynx" | metadata$BODY_SITE=="UBERON:vagina",] 

#Keep only the sex we're interested in:
metadata_sub <- metadata[metadata$SEX=="male" | metadata$SEX=="female",] 

#Check out number of samples for both genders:
nrow(metadata_sub[metadata_sub$SEX=="male",])
nrow(metadata_sub[metadata_sub$SEX=="female",])

#Read in the otu table:
otu_table <- read.table(otu_fp, sep="\t", 
                        comment="", 
                        header=T, 
                        skip=1, 
                        as.is=T, 
                        check.names=F, 
                        row=1)

RStudio.Version()$version
