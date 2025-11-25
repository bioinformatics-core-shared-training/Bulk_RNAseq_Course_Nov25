###This is demo for day3 
##Annotation for DESeq2 results
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(ggvenn)
library(pheatmap)

results.d11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
head(results.d11)

#querying the database 
ah <- AnnotationHub()
ah
ah[1]
dbResult <- query(ah,c("EnsDb","Mus musculus","102"))
dbResult
dbResult[1]

#download the result 
MouseEnsDb <- ah[[dbResult$ah_id]]
MouseEnsDb

#extraction genes
annotations <- genes(MouseEnsDb,return.type="data.frame")
view(annotations)
head(annotations)
colnames(annotations)
t(annotations[1, ])

#subset the table 
annot <- annotations %>% 
  dplyr::select(gene_id, gene_name, description, entrezid) %>%
  dplyr::filter(gene_id %in% rownames(results.d11))
view(annot)
length(annot$entrezid)
length(unique(annot$entrezid))
sum(is.na(annot$entrezid))
dupEntrez <- annot %>% 
  dplyr::filter(!is.na(entrezid)) %>%
  add_count(entrezid) %>%
  arrange(entrezid) %>%
  dplyr::filter(n>1)
head(dupEntrez)

annot <- readRDS("RObjects/Ensembl_annotations.rds")
view(annot)

#add annotation to our result 
results.d11 <- as.data.frame(results.d11) %>%
  rownames_to_column("GeneID") %>% 
  left_join(annot,"GeneID")
head(results.d11)
write_tsv(results.d11, "results.d11_Results_Annotated.txt") 

###Visualization of DE results 
#p-values sanity check 
hist(results.d11$pvalue)

#MA plot
results.d11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
ddsObj <- readRDS("RObjects/DESeqDataSet.interaction.rds")
annot <- readRDS("RObjects/Ensembl_annotations.rds")
shrink.11 <- lfcShrink(ddsObj, 
                       res = results.d11,
                       type = "ashr")
head(shrink.11)
shrinkTab.11 <- as.data.frame(shrink.11) %>% 
  rownames_to_column("GeneID") %>%
  left_join(annot,"GeneID")
head(shrinkTab.11)
saveRDS(shrinkTab.33, file="Shrunk_Results.d33.rds")
#MA plot 
plotMA(results.d11,alpha=0.05)
par(mfrow=c(1,2))
plotMA(results.d11,alpha=0.05)
plotMA(shrink.11,alpha=0.05)

#Volcano plots 
ggplot(shrinkTab.11,aes( x= log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = padj < 0.05),size=1) + 
  geom_text(data= ~top_n(.x, 1, wt = -padj), aes(label=Symbol)) +
  labs(x = "log2(fold change)", y = "-log10(p-value)", colour = "FDR<5%",
       title= "Infected vs Uninfected (day11)")
shrinkTab.11[shrinkTab.11$Symbol=="Irgm1",]


#exercise one 
results.d33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")
annot <- readRDS("RObjects/Ensembl_annotations.rds")
shrink.d33 <- lfcShrink(ddsObj, 
                        res = results.d33,
                        type = "ashr")
shrinkTab.33 <- as.data.frame(shrink.33) %>% 
  rownames_to_column("GeneID") %>%
  left_join(annot,"GeneID")
head(shrinkTab.33)
ggplot(shrinkTab.33,aes( x= log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = padj < 0.05),size=1) + 
  geom_text(data= ~top_n(.x, 1, wt = -padj), aes(label=Symbol)) +
  labs(x = "log2(fold change)", y = "-log10(p-value)", colour = "FDR<5%",
       title= "Infected vs Uninfected (day33)")
shrinkTab.33[shrinkTab.33$Symbol=="H2-Eb1",]

ggplot(shrinkTab.33,aes( x= log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = padj < 0.05),size=1) + 
  geom_text(data= ~top_n(.x, 10, wt = -padj), aes(label=Symbol)) +
  labs(x = "log2(fold change)", y = "-log10(p-value)", colour = "FDR<5%",
       title= "Infected vs Uninfected (day33)")

label_genes <- shrinkTab.33 %>% 
  filter(!is.na(padj)) %>%
  slice_min(padj,n=10)

ggplot(shrinkTab.33,aes( x= log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = padj < 0.05),size=1) + 
  geom_text_repel(data=label_genes,aes(label=Symbol),size=3,max.overlaps = Inf) +
  labs(x = "log2(fold change)", y = "-log10(p-value)", colour = "FDR<5%",
       title= "Infected vs Uninfected (day33)")

#MAplot with ggplot
# dev.off()
plotMA(shrink.33,alpha=0.05)
head(shrinkTab.33)
ggplot(shrinkTab.33, aes(x=log2(baseMean),y=log2FoldChange)) +
  geom_point(aes(colour=padj<0.05),size=1) + 
  scale_y_continuous(limit=c(-4,4),oob=scales::squish) +
  labs(x = "log2(Mean Expression)", 
       y = "log2(Fold Change)",
       colour = "FDR < 5%",
       title = "Infected vs Uninfected (day 33)") + theme_bw()


#strip chart to check expression of specific genes 
geneID <- filter(shrinkTab.11,Symbol=="Il10ra") %>%
  pull(GeneID)
geneID
plotCounts(ddsObj,
           gene=geneID,
           intgroup=c("TimePoint", "Status", "Replicate"),
           returnData = T) %>% 
  ggplot(aes(x=Status,y=log2(count))) +
  geom_point(aes(fill=Replicate),shape=21,size=2) +
  facet_wrap(~TimePoint) + 
  expand_limits(y=0) +
  labs(title = "Normalised counts - Interleukin 10 receptor, alpha") +theme_bw()



#exercise 3
geneID <- filter(shrinkTab.11,Symbol=="Jchain") %>%
  pull(GeneID)
geneID
plotCounts(ddsObj,
           gene=geneID,
           intgroup=c("TimePoint", "Status", "Replicate"),
           returnData = T) %>% 
  ggplot(aes(x=Status,y=log2(count))) +
  geom_point(aes(fill=Replicate),shape=21,size=2) +
  facet_wrap(~TimePoint) + 
  expand_limits(y=0) +
  labs(title = "Normalised counts - Jchain") +theme_bw()

shrinkTab.33[shrinkTab.33$Symbol=="Jchain",]

#Venn diagram of DE genes 
Upregulated_d11 <- shrinkTab.11 %>%
  filter(padj < 0.05) %>% 
  filter(log2FoldChange>0) %>% 
  pull("GeneID")
head(Upregulated_d11)

getGenes <- function(shrTab,direction){
  sign <- ifelse(direction=="up",1,-1)
  shrTab %>%
    filter(padj < 0.05) %>% 
    filter(log2FoldChange>0) %>% 
    pull("GeneID")
}
vennList <- list(Upregulated_d11 = getGenes(shrinkTab.11, "up"),
                 Downregulated_d11 = getGenes(shrinkTab.11, "down"),
                 Upregulated_d33 = getGenes(shrinkTab.33, "up"),
                 Downregulated_d33 = getGenes(shrinkTab.33, "down"))
vennList
str(vennList)
ggvenn(vennList,set_name_size=4)
#upset plot 

#heatmap of gene expression 
#filtering
sigGenes <- shrinkTab.11 %>%
  top_n(300,wt=-padj) %>%
  pull("GeneID")
#plot heatmap 
plotDat<-vst(ddsObj)[sigGenes,] %>%
  assay()
head(plotDat)
colours <- c("royalblue3", "ivory", "orangered3") 
hmPalette <- colorRampPalette(colours)(100)
hmPalette 
pheatmap(plotDat,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         scale="row",
         show_rownames = FALSE,
         color=hmPalette
)
#add annotation
annot_df<-colData(ddsObj) %>%
  as.data.frame() %>%
  select(Status, TimePoint)
annot_df
pheatmap(plotDat,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         scale="row",
         show_rownames = FALSE,
         show_colnames=FALSE,
         annotation_col=annot_df,
         color=hmPalette,
         cutree_rows=5
)

#change colours for annotation 
annot_col <- list(Status = c("Uninfected" = "darkgreen",
                             "Infected" = "palegreen"),
                  TimePoint = c("d11" = "lightblue",
                                "d33" = "darkblue"))


pheatmap(plotDat,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         scale="row",
         show_rownames = FALSE,
         show_colnames=FALSE,
         annotation_col=annot_df,
         annotation_colors=annot_col,
         color=hmPalette,
         cutree_rows=5
)























































































