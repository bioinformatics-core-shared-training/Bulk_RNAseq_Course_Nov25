library(tidyverse)
library(clusterProfiler)


# KEGG

search_kegg_organism('mouse', by='common_name') # mmu

shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")

sigGenes <- shrink.d11 %>%
  drop_na(Entrez, padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(Entrez)

keggRes <- enrichKEGG(gene = sigGenes, organism = 'mmu')
tab <- as_tibble(keggRes)

library(pathview)
logFC <- shrink.d11$log2FoldChange
names(logFC) <- shrink.d11$Entrez

pathview(gene.data = logFC,
         pathway.id = "mmu04612",
         species = "mmu",
         limit = list(gene = 20, cpd = 1))

# Exercise 1

logFC <- shrink.d11 %>%
  drop_na(padj, Entrez) %>%
  filter(padj < 0.01) %>%
  pull(log2FoldChange, Entrez)

pathview(gene.data = logFC,
         pathway.id = "mmu04659",
         species = "mmu",
         limit = list(gene = 5, cpd = 1))

# GO terms

library(org.Mm.eg.db)

sigGenes_GO <- shrink.d11 %>%
  drop_na(padj) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  pull(GeneID)

universe <- shrink.d11$GeneID

ego <- enrichGO(gene = sigGenes_GO,
                universe = universe,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pvalueCutoff = 0.01,
                readable = TRUE)

tab2 <- as_tibble(ego)

barplot(ego, showCategory = 20)
dotplot(ego, font.size = 14)

library(enrichplot)
ego_pt <- pairwise_termsim(ego)
emapplot(ego_pt)

# GSEA

library(msigdb)
library(ExperimentHub)

eh <- ExperimentHub()
query(eh, c('msigdb', 'mm', '2023'))
msigdb.mm <- getMsigdb(org = 'mm', id = 'EZID', version = '2023.1')

listCollections(msigdb.mm)

rankedGenes <- shrink.d11 %>%
  drop_na(GeneID, padj, log2FoldChange) %>%
  mutate(rank = log2FoldChange) %>%
  filter(!is.na(Entrez)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)
head(rankedGenes)

hallmarks <- subsetCollection(msigdb.mm, 'h')
msigdb_ids <- geneIds(hallmarks)

term2gene <- enframe(msigdb_ids, name = "gs_name", value = "entrez") %>%
  unnest(entrez)
term2gene

gseaRes <- GSEA(rankedGenes,
                TERM2GENE = term2gene,
                pvalueCutoff = 1,
                minGSSize = 15,
                maxGSSize = 500)

tab3 <- as_tibble(gseaRes) %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, wt=-p.adjust) %>% 
  dplyr::select(-core_enrichment) %>%
  mutate(across(c("enrichmentScore", "NES"), ~round(.x, digits=3))) %>% 
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))

gseaplot(gseaRes,
         geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE",
         title = "HALLMARK_INFLAMMATORY_RESPONSE")

# Exercise 2

rankedGenes.e11 <- shrink.d11 %>%
  drop_na(GeneID, padj, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  filter(!is.na(Entrez)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)

gseaRes.e11 <- GSEA(rankedGenes.e11,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1,
                    minGSSize = 15,
                    maxGSSize = 500)  
  
tab4 <- as_tibble(gseaRes.e11) %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, wt=-p.adjust) %>% 
  dplyr::select(-core_enrichment) %>%
  mutate(across(c("enrichmentScore", "NES"), ~round(.x, digits=3))) %>% 
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))
  
# day 33

shrink.d33 <- readRDS("RObjects/Shrunk_Results.d33.rds") 

rankedGenes.e33 <- shrink.d33 %>%
  drop_na(GeneID, padj, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  filter(!is.na(Entrez)) %>%
  arrange(desc(rank)) %>%
  pull(rank, Entrez)

gseaRes.e33 <- GSEA(rankedGenes.e33,
                    TERM2GENE = term2gene,
                    pvalueCutoff = 1,
                    minGSSize = 15,
                    maxGSSize = 500)  

tab5 <- as_tibble(gseaRes.e33) %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, wt=-p.adjust) %>% 
  dplyr::select(-core_enrichment) %>%
  mutate(across(c("enrichmentScore", "NES"), ~round(.x, digits=3))) %>% 
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))


  