
rm(list=ls())
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

#package needed to run this code
library(stringr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(car)
library(patchwork)
library(ggvenn)

"%not%" <- Negate("%in%") #to create a "not in" sign, the opposite of %in%

#________________________________

#--------Filter dataset 0.01%-----------------------------------------------
df <- read.csv("~/Downloads/Dataset/ASV-table-metabarcoding.csv", sep =";", row.names = 1)

df1 <- t(df)

otudf <- df1

df1 <- as.data.frame(t(apply(df1, 1, function(x) x/sum(x))))
rowSums(df1)

#remove asv < 0.01% of reads of the total reads but keeping reads number
for (i in 1:ncol(otudf)) {
  
  for (x in 1:nrow(otudf)) {
    
    if (df1[x,i] < 0.0001) { #to account for cross contamination and barcode switching during sequencing
      
      otudf[x,i] <- 0
    }
  }
}

rowSums(otudf)
colSums(otudf)

otudf <- otudf[,colSums(otudf) > 0] #from 1072 asv to 1070 asv

datafilt <- t(otudf)

write.csv(datafilt, "~/Downloads/Dataset/ASV-table-metabarcoding-filt.csv")


#--------create dataset with blast hits and metadata---------------------------------------------------------------

#load data
unoise4 <- read.csv("~/Downloads/Dataset/ASV-table-metabarcoding-filt.csv", row.names = 1)
metadata <- read.csv("~/Downloads/Dataset/MetadataCOImocksedimetoh20220106.csv", sep = ";")

otudf <- as.data.frame(t(unoise4))
otudf$baseclear_id <- row.names(otudf)
otudf <- otudf[,c(ncol(otudf), 1:(ncol(otudf)-1))]
row.names(otudf) <- NULL

#join metadata with asv data
df1 <- merge(metadata[,c(1,7:8)], otudf, by = "baseclear_id")

write.csv(df1, "~/Downloads/Dataset/workingdataset-metabarcoding.csv", row.names=FALSE)


#add taxonomy blast results - with a melt dataset
df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding.csv")
hits <- read.csv("~/Downloads/Dataset/otuhit75-all.csv", sep = ";")
taxo <- read.csv("~/Downloads/Dataset/otutaxo75-all.csv", sep = ";")
seq <- read.csv("~/Downloads/Dataset/sequences-315-325.csv", sep = ";")

taxo$Species <- as.character(as.factor(taxo$Species))
hits$Species <- as.character(as.factor(hits$Species))
taxo$Phylum <- as.character(as.factor(taxo$Phylum))
taxo$Class <- as.character(as.factor(taxo$Class))
taxo$Order <- as.character(as.factor(taxo$Order))
taxo$Family <- as.character(as.factor(taxo$Family))

df_m <- melt(df, id.vars = c("baseclear_id","registr.code","field.nmbr"), 
             value.name = "reads", variable.name = "Otus")

#attach taxo to hits
df2 <- as.data.frame(merge(hits, taxo[,c(1,3:7)], by = "Species"))
#attached sequences to ASVs
df3 <- merge(df_m, df2[,c(1,2,3,7:11)], by = "Otus", all = TRUE)
df4 <- merge(df3, seq, by = "Otus", all = TRUE)

df4 <- df4[df4$reads > 0,]
df4 <- df4[!is.na(df4$field.nmbr),]

write.csv(df4, "~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv", row.names=FALSE)


#--------create stats filtering criteria--------------------------------------

df_pro <- read.csv("~/Downloads/Dataset/ASV-table-metabarcoding.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/metabarcoding-Raw-reads-table.csv", sep = ";")
df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding.csv")

#get sum of reads/sample from processed data - ASV table without cutoff
df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[,c(1234,1233)]
#get sum of reads/sample from raw data
df_raw1 <- df_raw[1:nrow(df_raw),c(1,2)]
#get sum of reads/sample from ASV table with cutoff 0.1%
colnames(df)[2] <- "Sample_name"
df_cut <- df %>% group_by(Sample_name) %>% summarise(readscutoff = sum(reads))

df_all <- merge(df_raw1, df_pro1, by = "Sample_name", all = TRUE)
df_all$ratioretainedpercent <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100
df_all <- merge(df_all, df_cut, by = "Sample_name")

write.csv(df_all, "~/Downloads/Dataset/metabarcoding-readpersamplestats.csv", row.names=FALSE)


#--------Basic stats on the datasets-------------------------
df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv")

df$replicate <- substr(df$registr.code,13,13)
selected <- c("Mock1", "Mock2", "Mock3", "Mock4", "Mock5", "Mock6")

#how many ASVs not assigned to forams
notfor <- unique(df[df$Phylum %not% "Foraminifera",])
notfor <- na.omit(notfor)
notforottu <- unique(notfor$Otus)
notforphyl <- unique(notfor$Phylum)

#how many ASVs assigned to forams
foram <- unique(df[df$Phylum %in% "Foraminifera",])
foramdf <- unique(foram$Otus)

#subset for mocks
mocks <- df[df$field.nmbr %in% selected,]

#subset for sedim and etoh
sedimentetoh <- df[df$field.nmbr %not% selected,]
sedimentetoh$medium <- substr(sedimentetoh$field.nmbr,1,1)
sedimentetoh$location <- substr(sedimentetoh$field.nmbr,3,20)
sedim <- subset(sedimentetoh, subset = medium %in% "S")
etoh <- subset(sedimentetoh, subset = medium %in% "E")

mocks$counts <- 1
sedim$counts <- 1
etoh$counts <- 1

#basic stats on mock communities
allm <- mocks
mocks <- subset(mocks, mocks$Phylum == "Foraminifera")
phylm <- mocks[mocks$Identity > 75,]
clam <-  mocks[mocks$Identity > 80,]
ordm <-  mocks[mocks$Identity > 84,]
famm <-  mocks[mocks$Identity > 96,]
spm <-  mocks[mocks$Identity > 99.4,]

sum(famm$reads)/sum(phylm$reads)*100 #proportion of reads at the family
sum(spm$reads)/sum(phylm$reads)*100 #proportion of reads at the species

spmu <- levels(droplevels(mocks$Species[mocks$Identity > 99.4]))
notm <- levels(droplevels(mocks$Otus[is.na(mocks$Identity)]))

#basic stats on bulk-DNA from sediment samples
alls <- sedim
sedimna <- sedim[is.na(sedim$Phylum),]
sedim <- subset(sedim, sedim$Phylum == "Foraminifera")
phyls <- sedim[sedim$Identity > 75,]
clas <- sedim[sedim$Identity > 80,]
ords <- sedim[sedim$Identity > 84,]
fams <- sedim[sedim$Identity > 96,]
sps <- sedim[sedim$Identity > 99.4,]

sum(sedimna$reads)/sum(alls$reads)*100 #proportion of reads for non-assigned
sum(phyls$reads)/sum(alls$reads)*100 #proportion of reads at the phylum
sum(clas$reads)/sum(alls$reads)*100 #proportion of reads at the class
sum(ords$reads)/sum(alls$reads)*100 #proportion of reads at the order
sum(fams$reads)/sum(alls$reads)*100 #proportion of reads at the family
sum(sps$reads)/sum(alls$reads)*100 #proportion of reads at the species

spsu <- levels(droplevels(sedim$Species[sedim$Identity > 99.4]))
nots <- levels(droplevels(sedim$Otus[is.na(sedim$Identity)]))

#basic stats on eDNA from ethanol samples
alle <- etoh
etohna <- etoh[is.na(etoh$Phylum),]
etoh <- subset(etoh, etoh$Phylum == "Foraminifera")
phyle <- etoh[etoh$Identity > 75,]
clae <- etoh[etoh$Identity > 80,]
orde <- etoh[etoh$Identity > 84,]
fame <- etoh[etoh$Identity > 96,]
spe <- etoh[etoh$Identity > 99.4,]

sum(etohna$reads)/sum(alle$reads)*100 #proportion of reads not assigned
sum(phyle$reads)/sum(alle$reads)*100 #proportion of reads at the phylum
sum(clae$reads)/sum(alle$reads)*100 #proportion of reads at the class
sum(orde$reads)/sum(alle$reads)*100 #proportion of reads at the order
sum(fame$reads)/sum(alle$reads)*100 #proportion of reads at the family
sum(spe$reads)/sum(alle$reads)*100 #proportion of reads at the species

#--------Supplementary Fig. S4 - bubble chart for mocks---------------
df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv")

#subset
selected <- c("Mock1", "Mock2", "Mock3", "Mock4", "Mock5", "Mock6")
mocks <- df[df$field.nmbr %in% selected,]
mocks <- mocks[mocks$Identity > 99.4,1:ncol(mocks)]

mocks$replicate <- as.factor(substr(mocks$registr.code,13,13))
mocks$sample <- as.factor(substr(mocks$registr.code,1,12))


mocks <- mocks %>% group_by(field.nmbr, sample, registr.code,replicate, Species) %>% summarise(reads = sum(reads))

#include expected species in mocks (exp)
expected <- read.csv("C:/Temp/01-PhD-project/Chapter2-TestingCOImarker/Sequences-analysis/Project-21076-545/expected-mock-species.csv", sep = ";")
mocks <- rbind.data.frame(mocks, expected)

mocks <- mocks[order(mocks$Species, decreasing = TRUE),]
orderlevel <- unique(mocks$Species)

ggplot(mocks, aes(x=replicate , y=factor(Species, level = orderlevel), size = reads)) + 
  geom_point(aes(shape = replicate, color = replicate)) + 
  scale_shape_manual(values = c(16,16,16,16,16,4)) +
  scale_color_manual(values = c("black","black","black","black","black","red")) +
  facet_grid(~field.nmbr, scales="free_x") +
  theme_minimal() + xlab("Replicates") + ylab("Species (> 99.4% ID)")





#--------Supplementary Fig. S1 - finding best identity threshold using database------------------

dfm <- read.csv("~/Downloads/Dataset/Foram_COI_refs_25-01-2022-idmatrix.csv", sep = ";", row.names = 1)

#________CLASS: find optimum identity threshold___________
df_heatmap <- dfm
df_heatmap[is.na(df_heatmap)] <- 100

for (i in 1:ncol(df_heatmap)) {
  colnames(df_heatmap)[i] <- paste(str_split_fixed(colnames(df_heatmap)[i], "_", 4)[1], i)
}
for (i in 1:nrow(df_heatmap)) {
  rownames(df_heatmap)[i] <- paste(str_split_fixed(rownames(df_heatmap)[i], "_", 4)[1], i)
}

Breaks <- c(60, 80, 100)

myColor <- colorRampPalette(c("gray", "black"))(2)

#to see where the NA are
pheatmap(df_heatmap, border_color=NA,color = myColor, breaks = Breaks,
         cluster_cols = TRUE, cluster_rows = TRUE, main = "Identity threshold at 80% or Classes")

#________ORDER: find optimum identity threshold___________

df_heatmap <- dfm
df_heatmap[is.na(df_heatmap)] <- 100

for (i in 1:ncol(df_heatmap)) {
  colnames(df_heatmap)[i] <- paste(str_split_fixed(colnames(df_heatmap)[i], "_", 4)[2], i)
}
for (i in 1:nrow(df_heatmap)) {
  rownames(df_heatmap)[i] <- paste(str_split_fixed(rownames(df_heatmap)[i], "_", 4)[2], i)
}

Breaks <- c(60, 84, 100)

myColor <- colorRampPalette(c("gray", "black"))(2)

#to see where the NA are
pheatmap(df_heatmap, border_color=NA,color = myColor, breaks = Breaks,
         cluster_cols = TRUE, cluster_rows = TRUE, main = "Identity threshold at 84% for Orders")

#________Family: find optimum identity threshold___________

df_heatmap <- dfm
df_heatmap[is.na(df_heatmap)] <- 100

for (i in 1:ncol(df_heatmap)) {
  colnames(df_heatmap)[i] <- paste(str_split_fixed(colnames(df_heatmap)[i], "_", 4)[3], i)
}
for (i in 1:nrow(df_heatmap)) {
  rownames(df_heatmap)[i] <- paste(str_split_fixed(rownames(df_heatmap)[i], "_", 4)[3], i)
}

Breaks <- c(60, 96, 100)

myColor <- colorRampPalette(c("gray","black"))(2)

#to see where the NA are
pheatmap(df_heatmap, border_color=NA, color = myColor, breaks = Breaks,
         cluster_cols = TRUE, cluster_rows = TRUE, main = "Identity threshold at 96% for Families")

#________Species: find optimum identity threshold___________

df_heatmap <- dfm
df_heatmap[is.na(df_heatmap)] <- 100

for (i in 1:ncol(df_heatmap)) {
  colnames(df_heatmap)[i] <- paste(str_split_fixed(colnames(df_heatmap)[i], "_", 4)[4], i)
}
for (i in 1:nrow(df_heatmap)) {
  rownames(df_heatmap)[i] <- paste(str_split_fixed(rownames(df_heatmap)[i], "_", 4)[4], i)
}

Breaks <- c(60, 99.4, 100)

myColor <- colorRampPalette(c("gray","black"))(2)

#to see where the NA are
pheatmap(df_heatmap, border_color=NA, color = myColor, breaks = Breaks,
         cluster_cols = TRUE, cluster_rows = TRUE, main = "Identity threshold at 99.4% for Species")

#______Final heatmap____________________________________________

df_heatmap <- dfm
df_heatmap[is.na(df_heatmap)] <- 100

Breaks <- c(60, 80, 84, 96, 99.4, 100)

myColor <- colorRampPalette(c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))(5)

#to see where the NA are
pheatmap(df_heatmap, border_color=NA,color = myColor, breaks = Breaks,
         cluster_cols = TRUE, cluster_rows = TRUE)
         #main = "Sequence identity threshold (%) fitting Foraminifera taxonomy best")


#--------Figure 2 - pie charts to assess ratio-forams vs non-forams--------

df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv")

df$Species <- as.character(as.factor(df$Species))
df$Phylum <- as.character(as.factor(df$Phylum))
df$Class <- as.character(as.factor(df$Class))
df$Order <- as.character(as.factor(df$Order))
df$Family <- as.character(as.factor(df$Family))
df$Genus <- as.character(as.factor(df$Genus))
df$Phylum[is.na(df$Phylum)] <- "not assigned"
df$Species[is.na(df$Species)] <- "not assigned"
df$Class[is.na(df$Class)] <- "not assigned"
df$Order[is.na(df$Order)] <- "not assigned"
df$Family[is.na(df$Family)] <- "not assigned"
df$Genus[is.na(df$Genus)] <- "not assigned"


df$replicate <- substr(df$registr.code,13,13)
selected <- c("Mock1", "Mock2", "Mock3", "Mock4", "Mock5", "Mock6")

#subset for sedim and etoh
sedimentetoh <- df[df$field.nmbr %not% selected,]
sedimentetoh$medium <- substr(sedimentetoh$field.nmbr,1,1)
sedimentetoh$location <- substr(sedimentetoh$field.nmbr,3,20)
sedim <- subset(sedimentetoh, subset = medium %in% "S")
etoh <- subset(sedimentetoh, subset = medium %in% "E")

sedim$counts <- 1
sedim <- sedim[order(sedim$Species),]
sedim <- sedim[order(sedim$Genus),]
sedim <- sedim[order(sedim$Family),]
sedim <- sedim[order(sedim$Order),]
sedim <- sedim[order(sedim$Class),]
sedim <- sedim[order(sedim$Phylum),]
etoh$counts <- 1
etoh <- etoh[order(etoh$Species),]
etoh <- etoh[order(etoh$Genus),]
etoh <- etoh[order(etoh$Family),]
etoh <- etoh[order(etoh$Order),]
etoh <- etoh[order(etoh$Class),]
etoh <- etoh[order(etoh$Phylum),]

unique(sedim$Order)

#set colors:
Colors <- c('#E6194B','#3399FF',
            '#E6194B','#3399FF','#3CB44B',
            '#3CB44B','#3399FF','#CCFFFF',
            '#FF9999','#E6194B','#FF9999',
            '#99FF99','#3CB44B','#006633',
            '#CCFFFF','#99FFFF','#66FFFF','#99CCFF', 
            '#3399FF','#0066CC', '#0033CC', '#003366','#A9A9A9')
Names <- c('Foraminifera','not assigned',
           'Monothalamids','Globothalamea','Tubothalamea',
           'Miliolida','Rotaliida','Textulariida',
           'Group-3','Clade-E','Lacrogromiidae',
           'Alveolinidae','Peneroplidae','Soritidae',
           'Amphisteginidae','Calcarinidae','Glabratellidae','Murrayinellidae', 
           'Nummulitidae','Rosalinidae', 'undefined-clade', 'Uvigerinidae', 'other taxa')

pcolors <- as.data.frame(cbind(Names, Colors))
colnames(pcolors) <- c("Phylum", "Colors")
ccolors <- as.data.frame(cbind(Names, Colors))
colnames(ccolors) <- c("Class", "Colors")
ocolors <- as.data.frame(cbind(Names, Colors))
colnames(ocolors) <- c("Order", "Colors")
fcolors <- as.data.frame(cbind(Names, Colors))
colnames(fcolors) <- c("Family", "Colors")


#transform dataframe to get: Species code, and relative abundance for all years as column
sedim1 <- sedim %>% group_by(Phylum, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(sedim1)){
  if (sedim1$totalcounts[i]/sum(sedim1$totalcounts) < 0.02) {
    sedim1[i,1] <- "other taxa"
  }}
sedim1 <- sedim1 %>% group_by(Phylum) %>% summarise(totalcounts = sum(totalcounts))
sedim1 <- merge(sedim1, pcolors, by ="Phylum")
sedim1 <- sedim1[order(match(sedim1$Phylum,Names)),]
sedim2 <- sedim %>% subset(Phylum == "Foraminifera") %>% group_by(Class, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(sedim2)){
  if (sedim2$totalcounts[i]/sum(sedim2$totalcounts) < 0.02) {
    sedim2[i,1] <- "other taxa"
  }}
sedim2 <- sedim2 %>% group_by(Class) %>% summarise(totalcounts = sum(totalcounts))
sedim2 <- merge(sedim2, ccolors, by ="Class")
sedim2 <- sedim2[order(match(sedim2$Class,Names)),]
sedim3 <- sedim %>% subset(Phylum == "Foraminifera") %>% group_by(Order, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(sedim3)){
  if (sedim3$totalcounts[i]/sum(sedim3$totalcounts) < 0.02) {
    sedim3[i,1] <- "other taxa"
  }}
sedim3 <- sedim3 %>% group_by(Order) %>% summarise(totalcounts = sum(totalcounts))
sedim3 <- merge(sedim3, ocolors, by ="Order")
sedim3 <- sedim3[order(match(sedim3$Order,Names)),]
sedim4 <- sedim %>% subset(Phylum == "Foraminifera") %>% group_by(Family, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(sedim4)){
  if (sedim4$totalcounts[i]/sum(sedim4$totalcounts) < 0.02) {
    sedim4[i,1] <- "other taxa"
  }}
sedim4 <- sedim4 %>% group_by(Family) %>% summarise(totalcounts = sum(totalcounts))
sedim4 <- merge(sedim4, fcolors, by ="Family")
sedim4 <- sedim4[order(match(sedim4$Family,Names)),]

etoh1 <- etoh %>% group_by(Phylum, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(etoh1)){
  if (etoh1$totalcounts[i]/sum(etoh1$totalcounts) < 0.02) {
    etoh1[i,1] <- "other taxa"
  }}
etoh1 <- etoh1 %>% group_by(Phylum) %>% summarise(totalcounts = sum(totalcounts))
etoh1 <- merge(etoh1, pcolors, by ="Phylum")
etoh1 <- etoh1[order(match(etoh1$Phylum,Names)),]
etoh2 <- etoh %>% subset(Phylum == "Foraminifera") %>% group_by(Class, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(etoh2)){
  if (etoh2$totalcounts[i]/sum(etoh2$totalcounts) < 0.02) {
    etoh2[i,1] <- "other taxa"
  }}
etoh2 <- etoh2 %>% group_by(Class) %>% summarise(totalcounts = sum(totalcounts))
etoh2 <- merge(etoh2, ccolors, by ="Class")
etoh2 <- etoh2[order(match(etoh2$Class,Names)),]
etoh3 <- etoh %>% subset(Phylum == "Foraminifera") %>% group_by(Order, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(etoh3)){
  if (etoh3$totalcounts[i]/sum(etoh3$totalcounts) < 0.02) {
    etoh3[i,1] <- "other taxa"
  }}
etoh3 <- etoh3 %>% group_by(Order) %>% summarise(totalcounts = sum(totalcounts))
etoh3 <- merge(etoh3, ocolors, by ="Order")
etoh3 <- etoh3[order(match(etoh3$Order,Names)),]
etoh4 <- etoh %>% subset(Phylum == "Foraminifera") %>% group_by(Family, Otus) %>% summarise(Mean = mean(counts)) %>% summarise(totalcounts = sum(Mean))
for (i in 1:nrow(etoh4)){
  if (etoh4$totalcounts[i]/sum(etoh4$totalcounts) < 0.02) {
    etoh4[i,1] <- "other taxa"
  }}
etoh4 <- etoh4 %>% group_by(Family) %>% summarise(totalcounts = sum(totalcounts))
etoh4 <- merge(etoh4, fcolors, by ="Family")
etoh4 <- etoh4[order(match(etoh4$Family,Names)),]

#print 2 row of 4 piecharts
par(mfrow=c(2,4))
#to make a piechart slope 
lbls <- sedim1$Phylum
slices <- sedim1$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2,border=NA,col=as.character(sedim1$Colors),
    main="Sediment Phylum")
lbls <- sedim2$Class
slices <- sedim2$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2,border=NA, col=as.character(sedim2$Colors),
    main="Sediment Class")
lbls <- sedim3$Order
slices <- sedim3$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2, border=NA,col=as.character(sedim3$Colors),
    main="Sediment Order")
lbls <- sedim4$Family
slices <- sedim4$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2, border=NA,col=as.character(sedim4$Colors),
    main="Sediment Family")

lbls <- etoh1$Phylum
slices <- etoh1$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2,border=NA, col=as.character(etoh1$Colors),
    main="EtOH Phylum")
lbls <- etoh2$Class
slices <- etoh2$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2, border=NA,col=as.character(etoh2$Colors),
    main="EtOH Class")
lbls <- etoh3$Order
slices <- etoh3$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2, border=NA,col=as.character(etoh3$Colors),
    main="EtOH Order")
lbls <- etoh4$Family
slices <- etoh4$totalcounts
pct <- round(slices/sum(slices)*100)
lbls1 <- paste(lbls, pct) # add percents to labels
lbls1 <- paste(lbls1,"%",sep="") # ad % to labels
pie(slices,labels = lbls1, radius = 2,border=NA, col=as.character(etoh4$Colors),
    main="EtOH Family")

#--------Figure 3 - detected species vs expected in mocks---------------------------------------
df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv")

#subset
selected <- c("Mock1", "Mock2", "Mock3", "Mock4", "Mock5", "Mock6")
mocks <- df[df$field.nmbr %in% selected,]
mocks <- mocks[mocks$Identity > 99.4,1:ncol(mocks)]

mocks$replicate <- as.factor(substr(mocks$registr.code,13,13))
mocks$sample <- as.factor(substr(mocks$registr.code,1,12))

mocks <- mocks %>% group_by(field.nmbr, sample, registr.code,replicate, Species) %>% summarise(reads = sum(reads))

#include expected species in mocks (exp)
expected <- read.csv("C:/Temp/01-PhD-project/Chapter2-TestingCOImarker/Sequences-analysis/Project-21076-545/expected-mock-species.csv", sep = ";")
mocks <- rbind.data.frame(mocks, expected)

prabs <- mocks
prabs$reads <- 1

prabs <- dcast(prabs, field.nmbr + sample + Species ~ replicate, fun.aggregate = sum, value.var = "reads")

prabs$det <- 0

for (i in 1:nrow(prabs)) {
  if (sum(prabs[i,c("a","b","c","d","e")]) > ((ncol(select_if(prabs, is.numeric))-2)/2)) {
    prabs$det[i] <- 1
  } else {
    prabs$det[i] <- 0
  }
}

expdet <- prabs[,c(1:3,9,10)]
expdet <- melt(expdet, id.vars = c("field.nmbr","sample","Species"), value.name = "reads", variable.name = "replicate")
expdet$reads[expdet$reads == 0] <- NA
#expdet <- expdet[!(expdet$replicate == "exp" & expdet$reads == 0),]


expdet <- expdet[order(expdet$Species, decreasing = TRUE),]
orderlevel <- unique(expdet$Species)

ggplot(expdet, aes(x=replicate , y=factor(Species, level = orderlevel), size = reads)) + 
  geom_point(aes(shape = replicate, color = replicate)) + 
  scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(values = c("black","green4")) +
  facet_grid(~field.nmbr, scales="free_x") +
  theme_minimal() + xlab("Replicates") + ylab("Species (> 99.4% ID)") + theme(legend.position = "none")



#--------Figure 6 - detected species in bulk-DNA vs morphological samples---------------------------------------
df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv")

#subset
selected <- c("Mock1", "Mock2", "Mock3", "Mock4", "Mock5", "Mock6")
sedimentetoh <- df[df$field.nmbr %not% selected,]
sedimentetoh <- sedimentetoh[sedimentetoh$Identity > 99.4,1:ncol(sedimentetoh)]
sedimentetoh$replicate <- as.factor(substr(sedimentetoh$registr.code,13,13))
sedimentetoh$sample <- as.factor(substr(sedimentetoh$registr.code,1,12))

sedimentetoh <- sedimentetoh %>% group_by(field.nmbr, sample, registr.code, replicate, Species) %>% summarise(reads = sum(reads))
sedimentetoh$location <- substr(sedimentetoh$field.nmbr,3,20)

#include expected species in mocks (exp)
expected <- read.csv("~/Downloads/Dataset/expected-sedimentetoh-speciesv2.csv", sep = ";")
type <- read.csv("~/Downloads/Dataset/morpho-vs-edna-sedimentetoh.csv", sep = ";")
type <- unique(type[,c(3,11)])
new <- c("Calcarina_sp._NBCLAB5165", "LBF")
type <- rbind(type, new)

sedimentetoh <- rbind.data.frame(sedimentetoh, expected)

sedimentetoh <- merge(sedimentetoh, type, by = "Species")

sedimentetoh <- sedimentetoh[!is.na(sedimentetoh$location),]
sedimentetoh$medium <- substr(sedimentetoh$field.nmbr,1,1)
sedimentetoh$island <- substr(sedimentetoh$field.nmbr,3,7)
sedimentetoh$location[sedimentetoh$location == "UPG90-6"] <- "UPG90-06"

sedimentetoh$Species <- as.factor(sedimentetoh$Species)

#align names between morphological and referenced morphospecies in our database
levels(sedimentetoh$Species)[match("Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS",levels(sedimentetoh$Species))] <- "Amphisorus_spp."
levels(sedimentetoh$Species)[match("Amphisorus_SpS",levels(sedimentetoh$Species))] <- "Amphisorus_spp."
levels(sedimentetoh$Species)[match("Neorotalia_gaimardi_&_Baculogypsina_sphaerulata",levels(sedimentetoh$Species))] <- "Neorotalia_gaimardi"
levels(sedimentetoh$Species)[match("Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._NBCLAB5247",levels(sedimentetoh$Species))] <- "Calcarina_hispida_spengleri"
levels(sedimentetoh$Species)[match("Calcarina_spengleri",levels(sedimentetoh$Species))] <- "Calcarina_hispida_spengleri"
levels(sedimentetoh$Species)[match("Calcarina_hispida_&_Calcarina_sp_NBCLAB5163",levels(sedimentetoh$Species))] <- "Calcarina_hispida"
levels(sedimentetoh$Species)[match("Calcarina_sp._NBCLAB5164",levels(sedimentetoh$Species))] <- "Calcarina_mayori"
levels(sedimentetoh$Species)[match("Peneroplis_sp2_&_Peneroplis_pertusus_NBCLAB5117_&_Dendritina_ambigua",levels(sedimentetoh$Species))] <- "Peneroplis_spp."
levels(sedimentetoh$Species)[match("Peneroplis_sp1",levels(sedimentetoh$Species))] <- "Peneroplis_spp."
levels(sedimentetoh$Species)[match("Peneroplis_planatus",levels(sedimentetoh$Species))] <- "Peneroplis_spp."
levels(sedimentetoh$Species)[match("Peneroplis_pertusus",levels(sedimentetoh$Species))] <- "Peneroplis_spp."
levels(sedimentetoh$Species)[match("Heterostegina_depressa_sp1",levels(sedimentetoh$Species))] <- "Heterostegina_depressa"
levels(sedimentetoh$Species)[match("Heterostegina_depressa_sp2",levels(sedimentetoh$Species))] <- "Heterostegina_depressa"
levels(sedimentetoh$Species)[match("Sorites_sp1",levels(sedimentetoh$Species))] <- "Sorites_spp."
levels(sedimentetoh$Species)[match("Sorites_sp2",levels(sedimentetoh$Species))] <- "Sorites_spp."
levels(sedimentetoh$Species)[match("Sorites_orbiculus",levels(sedimentetoh$Species))] <- "Sorites_spp."
levels(sedimentetoh$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(sedimentetoh$Species))] <- "Amphistegina_lessonii"
levels(sedimentetoh$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(sedimentetoh$Species))] <- "Amphistegina_lessonii"
levels(sedimentetoh$Species)[match("Amphistegina_papillosa_Spermonde_NBCLAB3741",levels(sedimentetoh$Species))] <- "Amphistegina_papillosa"
levels(sedimentetoh$Species)[match("Parasorites_sp1",levels(sedimentetoh$Species))] <- "Parasorites_sp."

prabs <- sedimentetoh
prabs$reads <- 1

prabs1 <- dcast(prabs, field.nmbr + sample + medium + location + Type + Species ~ replicate, fun.aggregate = sum, value.var = "reads")

prabs1$det <- 0

sample <- unique(prabs$sample)

#make a loop to calculate the detected species 
#(present in more than half of the replicates)
for (i in 1:length(sample)) {
  
  x <- unique(prabs$replicate[prabs$sample == sample[[i]]])
  d <- as.character(x[-c(length(x))])
  
for (z in 1:nrow(prabs1)) {
  
  if (prabs1$sample[z] == sample[[i]] & sum(prabs1[z,d]) > ((length(d))/2)) {
    
    prabs1$det[z] <- 1
  }
}
}


expdet <- prabs1[,c(1:6,14,15)]
expdet$exp[expdet$exp > 1] <- 1
expdet <- melt(expdet, id.vars = c("field.nmbr","sample", "medium", "location","Species", "Type"), value.name = "reads", variable.name = "replicate")
expdet$reads[expdet$reads == 0] <- NA
#expdet <- expdet[!(expdet$replicate == "exp" & expdet$reads == 0),]

expdet <- expdet[!expdet$medium == "E",]

expdet$island <- substr(expdet$field.nmbr,3,7)

expdet <- expdet[order(expdet$Species, decreasing = TRUE),]
expdet <- expdet[order(expdet$Type, decreasing =TRUE),]
orderlevel <- unique(expdet$Species)
expdet <- expdet[order(expdet$medium, decreasing = TRUE),]
orderlevel1 <- unique(expdet$medium)



ggplot(expdet, aes(x=factor(medium, level = orderlevel1) , y=factor(Species, level = orderlevel), size = reads)) + 
  geom_point(aes(shape = replicate, color = replicate)) + 
  scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(values = c("black","green4")) +
  facet_grid(~location, scales="free_x") +
  theme_minimal() + xlab("Sample type") + ylab("Species (> 99.4% ID)") + theme(legend.position = "none")



#--------Figure 5 - bar charts for comparison bulk-DNA vs morphological samples---------------------------------------
df <- read.csv("~/Downloads/Dataset/morpho-vs-edna-sedimentetoh.csv", sep = ";")

df$medium <- substr(df$field.nmbr,1,1)
df$location <- substr(df$field.nmbr,3,20)
df$island <- substr(df$field.nmbr,3,7)
df$replicate <- substr(df$registr.code,13,13)
df$Species <- as.factor(df$Species)
df$Genus <- as.factor(df$Genus)
df$Family <- as.factor(df$Family)

df1 <- df %>% subset(Identity > 99.4) %>% group_by(island, location, medium, field.nmbr, Type, Phylum, Class, Order, Family, Genus, Species) %>% summarise(medreads = median(reads))
df2 <- df %>% subset(Identity > 99.4) %>% group_by(island, location, medium, field.nmbr, Type) %>% summarise(sumreads = sum(reads))
df3 <- df1 %>% subset(Type == "LBF")
df4 <- df1 %>% subset(Type == "other_foram")

#make the LBF from sedim-etoh comparable with traditional by grouping the species accordingly
df5 <- as.data.frame(df3)
levels(df5$Species)[match("Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS",levels(df5$Species))] <- "Amphisorus_spp."
levels(df5$Species)[match("Amphisorus_SpS",levels(df5$Species))] <- "Amphisorus_spp."
levels(df5$Species)[match("Neorotalia_gaimardi_&_Baculogypsina_sphaerulata",levels(df5$Species))] <- "Neorotalia_gaimardi"
levels(df5$Species)[match("Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._NBCLAB5247",levels(df5$Species))] <- "Calcarina_hispida_spengleri"
levels(df5$Species)[match("Calcarina_spengleri",levels(df5$Species))] <- "Calcarina_hispida_spengleri"
levels(df5$Species)[match("Calcarina_hispida_&_Calcarina_sp_NBCLAB5163",levels(df5$Species))] <- "Calcarina_hispida"
levels(df5$Species)[match("Calcarina_sp._NBCLAB5164",levels(df5$Species))] <- "Calcarina_mayori"
levels(df5$Species)[match("Peneroplis_sp2_&_Peneroplis_pertusus_NBCLAB5117_&_Dendritina_ambigua",levels(df5$Species))] <- "Peneroplis_spp."
levels(df5$Species)[match("Peneroplis_sp1",levels(df5$Species))] <- "Peneroplis_spp."
levels(df5$Species)[match("Peneroplis_planatus",levels(df5$Species))] <- "Peneroplis_spp."
levels(df5$Species)[match("Peneroplis_pertusus",levels(df5$Species))] <- "Peneroplis_spp."
levels(df5$Species)[match("Heterostegina_depressa_sp1",levels(df5$Species))] <- "Heterostegina_depressa"
levels(df5$Species)[match("Heterostegina_depressa_sp2",levels(df5$Species))] <- "Heterostegina_depressa"
levels(df5$Species)[match("Sorites_sp1",levels(df5$Species))] <- "Sorites_spp."
levels(df5$Species)[match("Sorites_sp2",levels(df5$Species))] <- "Sorites_spp."
levels(df5$Species)[match("Sorites_orbiculus",levels(df5$Species))] <- "Sorites_spp."
levels(df5$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(df5$Species))] <- "Amphistegina_lessonii"
levels(df5$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(df5$Species))] <- "Amphistegina_lessonii"
levels(df5$Species)[match("Amphistegina_papillosa_Spermonde_NBCLAB3741",levels(df5$Species))] <- "Amphistegina_papillosa"
levels(df5$Species)[match("Parasorites_sp1",levels(df5$Species))] <- "Parasorites_sp."
levels(df5$Genus)[match("Peneroplis_Dendritina",levels(df5$Genus))] <- "Peneroplis"
levels(df5$Genus)[match("Neorotalia_Baculogypsina",levels(df5$Genus))] <- "Neorotalia"

levels(droplevels(df5$Species))
levels(droplevels(df5$Genus))
levels(droplevels(df5$Family))

df5 <- df5[order(df5$Species),]
df5 <- df5[order(df5$Genus),]
df5 <- df5[order(df5$Family),]

fam <- unique(droplevels(df5$Family))
gen <- unique(droplevels(df5$Genus))
spe <- unique(droplevels(df5$Species))

#remove ethanol eDNA samples to just compare Bulk vs morphological
df5 <- df5[!df5$medium == "E",]

colspe <- c('#FFE119','#F58231','#FF9999', '#FF6666','#E6194B', '#990033',
            '#CCFFFF','#99CCFF', '#3399FF','#0066CC', '#0033CC', '#003366',
            '#CC0099', '#660066','#A9A9A9','#CCFFCC','#99FF99','#3CB44B','#006633')

ggplot(df5, aes(fill=factor(Species, level = spe), y=medreads, x=field.nmbr)) + 
  geom_bar(position="fill", stat="identity") + facet_grid(location~., scales="free_y") +
  scale_fill_manual(name = "Species",values = colspe) + theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank(), axis.title.y = element_blank()) +
  ggtitle("Relative abundance of LBF species") + coord_flip() +
  ylab("Relative abundance of species (> 99.4% ID)")


#--------Figure 1 and Sup. Fig. S2 - NMDS plots--------------

df <- read.csv("~/Downloads/Dataset/workingdataset-metabarcoding-taxo-seq.csv")

df$Species <- as.character(as.factor(df$Species))
df$Phylum <- as.character(as.factor(df$Phylum))
df$Class <- as.character(as.factor(df$Class))
df$Order <- as.character(as.factor(df$Order))
df$Family <- as.character(as.factor(df$Family))
df$Genus <- as.character(as.factor(df$Genus))
df$Phylum[is.na(df$Phylum)] <- "not assigned"
df$Species[is.na(df$Species)] <- "not assigned"
df$Class[is.na(df$Class)] <- "not assigned"
df$Order[is.na(df$Order)] <- "not assigned"
df$Family[is.na(df$Family)] <- "not assigned"
df$Genus[is.na(df$Genus)] <- "not assigned"

df$medium <- as.factor(substr(df$field.nmbr,1,1))
df$location <- substr(df$field.nmbr,3,20)
df$location[df$location == "UPG90-6"] <- "UPG90-06"
df$location <- as.factor(df$location)
df$island <- substr(df$field.nmbr,3,7)
df$island[df$island == "ck1"] <- "Mock1"
df$island[df$island == "ck2"] <- "Mock2"
df$island[df$island == "ck3"] <- "Mock3"
df$island[df$island == "ck4"] <- "Mock4"
df$island[df$island == "ck5"] <- "Mock5"
df$island[df$island == "ck6"] <- "Mock6"
df$island <- as.factor(df$island)
df$replicate <- as.factor(substr(df$registr.code,13,13))


#_________________________MOCKS_______________________________________________
mocks <- df[df$medium %in% "M",]

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(mocks, registr.code + medium + location + replicate + island ~ Otus, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(mocks, registr.code ~ Otus, fun.aggregate = sum, value.var = "reads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$registr.code
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,6:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(df_NMDS$island)
location <- droplevels(df_NMDS$location)
medium <- droplevels(df_NMDS$medium)
replicate <- df_NMDS$replicate

colmed <- c(15) 
colils <- c('#FFE119','#F58231','#E6194B','#3399FF','#660066','#3CB44B')

#Significance is the p value
ano <- anosim(df_NMDS_vegan, location)

#for plot
plot(ord, disp="sites", type = "n", main = paste ("Mock communities - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
points(ord, disp="sites", pch = colmed[medium], cex=1.5, col = colils[island])
with(ord, legend(x = "topright", legend = levels(island), col = colils, pch = c(15)))



#________________________Sedimetoh_______________________________________________
sedet <- df[df$medium %in% c("E","S"),]

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(sedet, registr.code + medium + location + replicate + island ~ Otus, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(sedet, registr.code ~ Otus, fun.aggregate = sum, value.var = "reads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$registr.code
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,6:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(df_NMDS$island)
location <- droplevels(df_NMDS$location)
medium <- droplevels(df_NMDS$medium)
replicate <- droplevels(df_NMDS$replicate)

colmed <- c(17,19) 
colils <- c('#FF9999','#E6194B', '#990033','#99CCFF','#3399FF','#0033CC',
            '#A9A9A9','#99FF99','#3CB44B','#006633')

#Significance is the p value
ano <- anosim(df_NMDS_vegan, location) 
ano <- anosim(df_NMDS_vegan, medium)

#for plot
plot(ord, disp="sites", type = "n", main = paste ("Sediment & EtOH Samples - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3))) 
points(ord, disp="sites", pch = colmed[medium], cex=1.5, col = colils[location]) 
with(ord, legend(x = "topright", legend = levels(location), col = colils, pch = c(15))) 
with(ord, legend(x = "bottomright", legend = levels(medium), pch = colmed))


#________________________EtoH_______________________________________________
sedet <- df[df$medium %in% "E",]

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(sedet, registr.code + medium + location + replicate + island ~ Otus, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(sedet, registr.code ~ Otus, fun.aggregate = sum, value.var = "reads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$registr.code
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,6:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(df_NMDS$island)
location <- droplevels(df_NMDS$location)
medium <- droplevels(df_NMDS$medium)
replicate <- droplevels(df_NMDS$replicate)

colmed <- c(17) 
colils <- c('#FF9999','#E6194B', '#990033','#99CCFF','#3399FF','#0033CC',
            '#A9A9A9','#99FF99','#3CB44B','#006633')

#Significance is the p value
ano <- anosim(df_NMDS_vegan, location) 
#for plot
plot(ord, disp="sites", type = "n", main = paste ("EtOH Samples - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3))) 
points(ord, disp="sites", pch = colmed[medium], cex=1.5, col = colils[location]) 
with(ord, legend(x = "topright", legend = levels(location), col = colils, pch = c(17))) 

#________________________Sedim_______________________________________________
sedet <- df[df$medium %in% "S",]

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(sedet, registr.code + medium + location + replicate + island ~ Otus, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(sedet, registr.code ~ Otus, fun.aggregate = sum, value.var = "reads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$registr.code
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,6:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(df_NMDS$island)
location <- droplevels(df_NMDS$location)
medium <- droplevels(df_NMDS$medium)
replicate <- droplevels(df_NMDS$replicate)

colmed <- c(19) 
colils <- c('#FF9999','#E6194B', '#990033','#99CCFF','#3399FF','#0033CC',
            '#A9A9A9','#99FF99','#3CB44B','#006633')

#Significance is the p value
ano <- anosim(df_NMDS_vegan, location) 
#for plot
plot(ord, disp="sites", type = "n", main = paste ("Sediment Samples - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3))) 
points(ord, disp="sites", pch = colmed[medium], cex=1.5, col = colils[location]) 
with(ord, legend(x = "topright", legend = levels(location), col = colils, pch = c(19))) 

#--------Figure 4 - NMDS morphological vs bulk-DNA LBF--------------------
df <- read.csv("~/Downloads/Dataset/morpho-vs-edna-sedimentetoh.csv", sep = ";")

df$Phylum <- as.factor(df$Phylum)
df$Class <- as.factor(df$Class)
df$Order <- as.factor(df$Order)
df$Family <- as.factor(df$Family)
df$Genus <- as.factor(df$Genus)
df$Species <- as.factor(df$Species)

levels(df$Species)[match("Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS",levels(df$Species))] <- "Amphisorus_spp."
levels(df$Species)[match("Amphisorus_SpS",levels(df$Species))] <- "Amphisorus_spp."
levels(df$Species)[match("Neorotalia_gaimardi_&_Baculogypsina_sphaerulata",levels(df$Species))] <- "Neorotalia_gaimardi"
levels(df$Species)[match("Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._NBCLAB5247",levels(df$Species))] <- "Calcarina_hispida_spengleri"
levels(df$Species)[match("Calcarina_spengleri",levels(df$Species))] <- "Calcarina_hispida_spengleri"
levels(df$Species)[match("Calcarina_hispida_&_Calcarina_sp_NBCLAB5163",levels(df$Species))] <- "Calcarina_hispida"
levels(df$Species)[match("Calcarina_sp._NBCLAB5164",levels(df$Species))] <- "Calcarina_mayori"
levels(df$Species)[match("Peneroplis_sp2_&_Peneroplis_pertusus_NBCLAB5117_&_Dendritina_ambigua",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Peneroplis_sp1",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Peneroplis_planatus",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Peneroplis_pertusus",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Heterostegina_depressa_sp1",levels(df$Species))] <- "Heterostegina_depressa"
levels(df$Species)[match("Heterostegina_depressa_sp2",levels(df$Species))] <- "Heterostegina_depressa"
levels(df$Species)[match("Sorites_sp1",levels(df$Species))] <- "Sorites_spp."
levels(df$Species)[match("Sorites_sp2",levels(df$Species))] <- "Sorites_spp."
levels(df$Species)[match("Sorites_orbiculus",levels(df$Species))] <- "Sorites_spp."
levels(df$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(df$Species))] <- "Amphistegina_lessonii"
levels(df$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(df$Species))] <- "Amphistegina_lessonii"
levels(df$Species)[match("Amphistegina_papillosa_Spermonde_NBCLAB3741",levels(df$Species))] <- "Amphistegina_papillosa"
levels(df$Species)[match("Parasorites_sp1",levels(df$Species))] <- "Parasorites_sp."
levels(df$Genus)[match("Peneroplis_Dendritina",levels(df$Genus))] <- "Peneroplis"
levels(df$Genus)[match("Neorotalia_Baculogypsina",levels(df$Genus))] <- "Neorotalia"

df$Species <- as.character(as.factor(df$Species))

df$medium <- as.factor(substr(df$field.nmbr,1,1))
df$location <- substr(df$field.nmbr,3,20)
df$location[df$location == "UPG90-6"] <- "UPG90-06"
df$location <- as.factor(df$location)
df$island <- substr(df$field.nmbr,3,7)
df$island[df$island == "ck1"] <- "Mock1"
df$island[df$island == "ck2"] <- "Mock2"
df$island[df$island == "ck3"] <- "Mock3"
df$island[df$island == "ck4"] <- "Mock4"
df$island[df$island == "ck5"] <- "Mock5"
df$island[df$island == "ck6"] <- "Mock6"
df$island <- as.factor(df$island)
df$replicate <- as.factor(substr(df$registr.code,13,13))

sedet <- df[df$Type %in% "LBF",]
sedet <- sedet[sedet$Identity > 99.4,]

sedet <- sedet[!sedet$medium == "E",]

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(sedet, registr.code + medium + location + replicate + island ~ Species, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(sedet, registr.code ~ Species, fun.aggregate = sum, value.var = "reads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$registr.code
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
#df_NMDS_vegan[df_NMDS_vegan > 0] <- 1
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,6:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(df_NMDS$island)
location <- droplevels(df_NMDS$location)
medium <- droplevels(df_NMDS$medium)
replicate <- droplevels(df_NMDS$replicate)

colmed <- c(19, 8) 
colils <- c('#FF9999','#E6194B', '#990033','#99CCFF','#3399FF','#0033CC',
            '#A9A9A9','#99FF99','#3CB44B','#006633')

#Significance is the p value
ano <- anosim(df_NMDS_vegan, location) 
#for plot
plot(ord, disp="sites", type = "n", main = paste ("Traditional & sediment - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3))) 
points(ord, disp="sites", pch = colmed[medium], cex=2, col = colils[location]) 
with(ord, legend(x = "topright", legend = levels(location), col = colils, pch = c(15))) 
with(ord, legend(x = "bottomright", legend = levels(medium), pch = colmed))
plot(en, col = "black")

#--------Supplementary Fig. S3 - Diversity indexes bulk-DNA, eDNA, morphological-------------------------

df <- read.csv("~/Downloads/Dataset/morpho-vs-edna-sedimentetoh.csv", sep = ";")
df1 <- df %>% subset(Identity > 99.4)
LBF <-  df1 %>% subset(Type == "LBF")

df1 <- dcast(data = df1, field.nmbr + registr.code ~ Species, fun.aggregate = sum, value.var = "reads")
LBF1 <- dcast(data = LBF, field.nmbr + registr.code ~ Species, fun.aggregate = sum, value.var = "reads")

#Community information on the samples
df2 <- df1
df2$Richness_s <- apply(df1[,3:ncol(df1)]>0,1,sum)
df2$Abundance <- apply(df1[,3:ncol(df1)],1,sum)
#Diversity indexes
df2$Shannon_H <- vegan::diversity(df1[,3:ncol(df1)], index="shannon")
df2$Simpson_D <- vegan::diversity(df1[,3:ncol(df1)], index="simpson")
df2$Evenness_J <- vegan::diversity(df1[,3:ncol(df1)], index="simpson")/log(df2$Richness)
#True diversiy indexes
df2$True_shannon <- exp(df2$Shannon_H)
df2$True_simpson <- 1/df2$Simpson_D

#melt dataset to plot diversity index per habitat, per year, per index
df3 <- melt(df2, id.vars=c("field.nmbr", "registr.code", "Richness_s", "Abundance", "Shannon_H", 
                           "Simpson_D", "Evenness_J", "True_shannon", "True_simpson"),
            variable.name = "Species_Code", value.name = "value")

df4 <- unique(melt(df3, id.vars=c("field.nmbr", "registr.code"),
            variable.name = "Diversity_Index", value.name = "Value"))

df5 <- subset(df4, subset = Diversity_Index %in% c("Richness_s", "Shannon_H", 
                                                   "Simpson_D", "Evenness_J", "True_shannon", "True_simpson"))
df5$medium <- as.factor(substr(df5$field.nmbr,1,1))
df5$Value <- as.numeric(df5$Value)
df5$location <- substr(df5$field.nmbr,3,20)

#only ethanol and sediment
df6 <- subset(df5, subset = medium %in% c("S", "E"))
df6$type <- "all"          

LBF2 <- LBF1
LBF2$Richness_s <- apply(LBF1[,3:ncol(LBF1)]>0,1,sum)
LBF2$Abundance <- apply(LBF1[,3:ncol(LBF1)],1,sum)
#Diversity indexes
LBF2$Shannon_H <- vegan::diversity(LBF1[,3:ncol(LBF1)], index="shannon")
LBF2$Simpson_D <- vegan::diversity(LBF1[,3:ncol(LBF1)], index="simpson")
LBF2$Evenness_J <- vegan::diversity(LBF1[,3:ncol(LBF1)], index="simpson")/log(LBF2$Richness)
#True diversiy indexes
LBF2$True_shannon <- exp(LBF2$Shannon_H)
LBF2$True_simpson <- 1/LBF2$Simpson_D

#melt dataset to plot diversity index per habitat, per year, per index
LBF3 <- melt(LBF2, id.vars=c("field.nmbr", "registr.code", "Richness_s", "Abundance", "Shannon_H", 
                           "Simpson_D", "Evenness_J", "True_shannon", "True_simpson"),
            variable.name = "Species_Code", value.name = "value")

LBF4 <- unique(melt(LBF3, id.vars=c("field.nmbr", "registr.code"),
                   variable.name = "Diversity_Index", value.name = "Value"))

LBF5 <- subset(LBF4, subset = Diversity_Index %in% c("Richness_s", "Shannon_H", 
                                                   "Simpson_D", "Evenness_J", "True_shannon", "True_simpson"))
LBF5$medium <- as.factor(substr(LBF5$field.nmbr,1,1))
LBF5$Value <- as.numeric(LBF5$Value)
LBF5$location <- substr(LBF5$field.nmbr,3,20)

LBF5$type <- "LBF"

#put both LBF and all together
df7 <- rbind(df6, LBF5)
df8 <- subset(df7, subset = Diversity_Index %in% c("Richness_s"))

ggplot(df8, aes(x=Value, y=field.nmbr)) + 
  geom_boxplot(aes(col = medium)) +
  facet_grid(location~type, scales="free") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Species richness")

stats <- df8 %>% group_by(medium, type) %>% summarise(average = mean(Value))


#--------Add-on to Figure 4 and 5 - Venn diagram LBF morphological, bulk-DNA------------------

df <- read.csv("~/Downloads/Dataset/morpho-vs-edna-sedimentetoh.csv", sep = ";")
df1 <- df %>% subset(Identity > 99.4)
LBF <-  df1 %>% subset(Type == "LBF")
LBF$medium <- as.factor(substr(LBF$field.nmbr,1,1))
LBF$island <- substr(LBF$field.nmbr,3,7)

LBF$Species <- as.factor(LBF$Species)

levels(LBF$Species)[match("Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS",levels(LBF$Species))] <- "Amphisorus_spp."
levels(LBF$Species)[match("Amphisorus_SpS",levels(LBF$Species))] <- "Amphisorus_spp."
levels(LBF$Species)[match("Neorotalia_gaimardi_&_Baculogypsina_sphaerulata",levels(LBF$Species))] <- "Neorotalia_gaimardi"
levels(LBF$Species)[match("Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._NBCLAB5247",levels(LBF$Species))] <- "Calcarina_hispida_spengleri"
levels(LBF$Species)[match("Calcarina_spengleri",levels(LBF$Species))] <- "Calcarina_hispida_spengleri"
levels(LBF$Species)[match("Calcarina_hispida_&_Calcarina_sp_NBCLAB5163",levels(LBF$Species))] <- "Calcarina_hispida"
levels(LBF$Species)[match("Calcarina_sp._NBCLAB5164",levels(LBF$Species))] <- "Calcarina_mayori"
levels(LBF$Species)[match("Peneroplis_sp2_&_Peneroplis_pertusus_NBCLAB5117_&_Dendritina_ambigua",levels(LBF$Species))] <- "Peneroplis_spp."
levels(LBF$Species)[match("Peneroplis_sp1",levels(LBF$Species))] <- "Peneroplis_spp."
levels(LBF$Species)[match("Peneroplis_planatus",levels(LBF$Species))] <- "Peneroplis_spp."
levels(LBF$Species)[match("Peneroplis_pertusus",levels(LBF$Species))] <- "Peneroplis_spp."
levels(LBF$Species)[match("Heterostegina_depressa_sp1",levels(LBF$Species))] <- "Heterostegina_depressa"
levels(LBF$Species)[match("Heterostegina_depressa_sp2",levels(LBF$Species))] <- "Heterostegina_depressa"
levels(LBF$Species)[match("Sorites_sp1",levels(LBF$Species))] <- "Sorites_spp."
levels(LBF$Species)[match("Sorites_sp2",levels(LBF$Species))] <- "Sorites_spp."
levels(LBF$Species)[match("Sorites_orbiculus",levels(LBF$Species))] <- "Sorites_spp."
levels(LBF$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(LBF$Species))] <- "Amphistegina_lessonii"
levels(LBF$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(LBF$Species))] <- "Amphistegina_lessonii"
levels(LBF$Species)[match("Amphistegina_papillosa_Spermonde_NBCLAB3741",levels(LBF$Species))] <- "Amphistegina_papillosa"
levels(LBF$Species)[match("Parasorites_sp1",levels(LBF$Species))] <- "Parasorites_sp."

LBF$Species <- as.character(as.factor(LBF$Species))

LBF <- LBF[!LBF$medium == "E",] #we do not take ethanol into account.

#per medium
x <- list(
  Morphological = unique(LBF$Species[LBF$medium == "T"]), 
  Sediment = unique(LBF$Species[LBF$medium == "S"])
  
)

ggvenn(
 x, 
  fill_color = c('#E6194B', '#3399FF', '#3CB44B'),
   stroke_color = NA,
  stroke_size = 0, set_name_size = 6,
 show_percentage = FALSE
)

#per island
v <- list(
  Barangbaringan = unique(LBF$Species[LBF$island == "UPG82"]), 
  Langkadea = unique(LBF$Species[LBF$island == "UPG90"]),
  Pajenekang = unique(LBF$Species[LBF$island == "UPG91"]),
  "Bone Lola" = unique(LBF$Species[LBF$island == "UPG92"])
)

ggvenn(
  v, 
  fill_color = c('#E6194B', '#3399FF','#A9A9A9', '#3CB44B'),
  stroke_color = NA,
  stroke_size = 0, set_name_size = 6, 
  show_percentage = FALSE
)



#--------PermANOVA analysis to check differences between samples types and sites-----------------------

#________A. bulk-DNA and eDNA different at a given site_______________

df <- read.csv("~/Downloads/Dataset/morpho-vs-edna-sedimentetoh.csv", sep = ";")

df$Phylum <- as.factor(df$Phylum)
df$Class <- as.factor(df$Class)
df$Order <- as.factor(df$Order)
df$Family <- as.factor(df$Family)
df$Genus <- as.factor(df$Genus)
df$Species <- as.factor(df$Species)

levels(df$Species)[match("Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS",levels(df$Species))] <- "Amphisorus_spp."
levels(df$Species)[match("Amphisorus_SpS",levels(df$Species))] <- "Amphisorus_spp."
levels(df$Species)[match("Neorotalia_gaimardi_&_Baculogypsina_sphaerulata",levels(df$Species))] <- "Neorotalia_gaimardi"
levels(df$Species)[match("Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._NBCLAB5247",levels(df$Species))] <- "Calcarina_hispida_spengleri"
levels(df$Species)[match("Calcarina_spengleri",levels(df$Species))] <- "Calcarina_hispida_spengleri"
levels(df$Species)[match("Calcarina_hispida_&_Calcarina_sp_NBCLAB5163",levels(df$Species))] <- "Calcarina_hispida"
levels(df$Species)[match("Calcarina_sp._NBCLAB5164",levels(df$Species))] <- "Calcarina_mayori"
levels(df$Species)[match("Peneroplis_sp2_&_Peneroplis_pertusus_NBCLAB5117_&_Dendritina_ambigua",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Peneroplis_sp1",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Peneroplis_planatus",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Peneroplis_pertusus",levels(df$Species))] <- "Peneroplis_spp."
levels(df$Species)[match("Heterostegina_depressa_sp1",levels(df$Species))] <- "Heterostegina_depressa"
levels(df$Species)[match("Heterostegina_depressa_sp2",levels(df$Species))] <- "Heterostegina_depressa"
levels(df$Species)[match("Sorites_sp1",levels(df$Species))] <- "Sorites_spp."
levels(df$Species)[match("Sorites_sp2",levels(df$Species))] <- "Sorites_spp."
levels(df$Species)[match("Sorites_orbiculus",levels(df$Species))] <- "Sorites_spp."
levels(df$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(df$Species))] <- "Amphistegina_lessonii"
levels(df$Species)[match("Amphistegina_lessonii_NBCLAB5174",levels(df$Species))] <- "Amphistegina_lessonii"
levels(df$Species)[match("Amphistegina_papillosa_Spermonde_NBCLAB3741",levels(df$Species))] <- "Amphistegina_papillosa"
levels(df$Species)[match("Parasorites_sp1",levels(df$Species))] <- "Parasorites_sp."
levels(df$Genus)[match("Peneroplis_Dendritina",levels(df$Genus))] <- "Peneroplis"
levels(df$Genus)[match("Neorotalia_Baculogypsina",levels(df$Genus))] <- "Neorotalia"

df$Species <- as.character(as.factor(df$Species))

df$medium <- as.factor(substr(df$field.nmbr,1,1))
df$location <- substr(df$field.nmbr,3,20)
df$location[df$location == "UPG90-6"] <- "UPG90-06"
df$location <- as.factor(df$location)
df$island <- substr(df$field.nmbr,3,7)
df$island[df$island == "ck1"] <- "Mock1"
df$island[df$island == "ck2"] <- "Mock2"
df$island[df$island == "ck3"] <- "Mock3"
df$island[df$island == "ck4"] <- "Mock4"
df$island[df$island == "ck5"] <- "Mock5"
df$island[df$island == "ck6"] <- "Mock6"
df$island <- as.factor(df$island)
df$replicate <- as.factor(substr(df$registr.code,13,13))

#sedet <- df[df$Identity > 99.4,]
sedet <- df
sedet <- sedet[!sedet$medium == "T",]
sedet <- na.omit(sedet)

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(sedet, registr.code + medium + location + replicate + island ~ Species, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- na.omit(dcast(sedet, registr.code ~ Species, fun.aggregate = sum, value.var = "reads"))
rownames(df_NMDS_vegan) <- na.omit(df_NMDS_vegan$registr.code)
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
#df_NMDS_vegan[df_NMDS_vegan > 0] <- 1

data <- df_NMDS_vegan
env <- df_NMDS[,1:5]

#permanova analysis that there is a significant difference between eDNA and bulk-DNA communities
adonis2(data ~ medium, data = env)

#permanova analysis comparing each sample site between eDNA and bulk-DNA communities
all <- cbind(env, data)

interest_S <- all[all$medium %in% "S",]
interest_E <- all[all$medium %in% "E",]

#bulk-DNA samples
data_S <- interest_S[,6:ncol(interest_S)]
env_S <- interest_S[,1:5]
adonis2(data_S ~ location, data = env_S) #permanova analysis

#morphological samples
data_E <- interest_E[, 6:ncol(interest_E)]
env_E <- interest_E[,1:5]
adonis2(data_E ~ location, data = env_E) #permanova analysis

locations <- as.character(unique(all$location))

grouping <- c(1:3)
location <- c(1:3)
r2 <- c(1:3)
f_value <- c(1:3)
p_value <- c(1:3)
res_ado <- data.frame(cbind(grouping, location, r2, f_value, p_value))


for (i in 1: length(locations)) {
  interest_a <- all[all$location %in% locations[[i]],]
  
  if (length(unique(interest_a$medium)) < 2) {next}
  else {
  
  data_a <- interest_a[, 6:ncol(interest_a)]
  env_a <- interest_a[,1:5]
  ado <- adonis2(data_a ~ medium, data = env_a) #permanova analysis
  
  res_ado[i, 1] <- "medium"
  res_ado[i, 2] <- locations[[i]]
  res_ado[i, 3] <- ado[1,3] #R2
  res_ado[i, 4] <- ado[1,4] #F-value
  res_ado[i, 5] <- ado[1,5] #p-value
   
}
}


#________B. difference between UPG82-UPG92 for molecular vs morphological_______________

sedet <- df[df$Identity > 99.4,]
sedet <- sedet[!sedet$medium == "E",]
sedet <- na.omit(sedet)

#To build the NMDS, I have to "unmelt" the dataframe with reshape
df_NMDS <- dcast(sedet, registr.code + medium + location + replicate + island ~ Species, fun.aggregate = sum, value.var = "reads")

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- na.omit(dcast(sedet, registr.code ~ Species, fun.aggregate = sum, value.var = "reads"))
rownames(df_NMDS_vegan) <- na.omit(df_NMDS_vegan$registr.code)
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
#df_NMDS_vegan[df_NMDS_vegan > 0] <- 1

data <- df_NMDS_vegan
env <- df_NMDS[,1:5]

all <- cbind(env, data)

interest_S <- all[all$island %in% c("UPG82", "UPG92") & all$medium %in% "S",]
interest_T <- all[all$island %in% c("UPG82", "UPG92") & all$medium %in% "T",]

#bulk-DNA samples
data_S <- interest_S[,6:ncol(interest_S)]
env_S <- interest_S[,1:5]
adonis2(data_S ~ island, data = env_S) #permanova analysis

#morphological samples
data_T <- interest_T[, 6:ncol(interest_T)]
env_T <- interest_T[,1:5]
adonis2(data_T ~ island, data = env_T) #permanova analysis




