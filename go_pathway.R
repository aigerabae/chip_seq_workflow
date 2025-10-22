## 3.	Gene Ontology (GO) category-based enrichment.
Cellular pathways define a broad function within a cell, whereas Gene Ontology categories (Biological process, BP and Molecular function, MF) define a more precise function. Hence GO categories can give a finer resolution of the functional phenotype for a list of genes. You would test enrichment for the MF and BP classes of GO categories. There is also a third class of GO category: Cellular component (CC), enrichment analysis for whom makes sense when you suspect your genes to belong to specific cellular compartment such as lysosome or cell membrane.


  ## For mouse genome annotation
library(org.Mm.eg.db)
## For enrichment calculation and plotting
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggnewscale)
##
my.gene.list <- read.table("shared_genes.txt",header=FALSE,stringsAsFactors = FALSE)
## Check the structure
str(my.gene.list)
# 'data.frame':	528 obs. of  1 variable:
#  $ V1: chr  "0610040F04Rik" "1010001N08Rik" "1110002J07Rik" "1600002D24Rik" ...
##
## There are 528 gene symbols read
## Are there any duplicates?
length(unique(my.gene.list$V1))
# [1] 528
## None. So, list is ready to be processed


## GO enrichment for Molecular function (MF) class
enrchmnt.GO.mf <- enrichGO(gene = my.gene.list$V1, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
## Note that by default, the above function test for molecular function (MF) ontology class
## You can access the enrichment results like -
dim(enrchmnt.GO.mf@result)
# [1] 635   9
## What are the column names?
colnames(enrchmnt.GO.mf@result)
# [1] "ID"          "Description" "GeneRatio"   "BgRatio"     "pvalue"      "p.adjust"    "qvalue"      "geneID"     
# [9] "Count"      

##
## You can look at the first 5 entries by printing them out
## Skipping the 8th column as it has long lists of gene symbols
(enrchmnt.GO.mf@result)[1:5,c(1:7,9)]
#                    ID                                                              Description GeneRatio   BgRatio
# GO:0001227 GO:0001227 DNA-binding transcription repressor activity, RNA polymerase II-specific    19/463 328/28438
# GO:0001217 GO:0001217                             DNA-binding transcription repressor activity    19/463 330/28438
# GO:0001664 GO:0001664                                       G protein-coupled receptor binding    18/463 318/28438
# GO:0017124 GO:0017124                                                       SH3 domain binding    11/463 133/28438
# GO:0005096 GO:0005096                                                GTPase activator activity    20/463 405/28438
#                  pvalue     p.adjust       qvalue Count
# GO:0001227 2.209745e-06 0.0007391321 0.0006458189    19
# GO:0001217 2.415464e-06 0.0007391321 0.0006458189    19
# GO:0001664 5.612370e-06 0.0011449236 0.0010003804    18
# GO:0017124 1.202033e-05 0.0015933358 0.0013921820    11
# GO:0005096 1.301745e-05 0.0015933358 0.0013921820    20
##
## How many enriched molecular function (MF) categories?
length(which(enrchmnt.GO.mf@result$p.adjust < 0.05))
# [1] 17
##
## Now you can make and save plots
## Dotplot to show the desired number of top categories (based on adjusted p-val)
dotplot(enrchmnt.GO.mf,showCategory=5)
##
## Save the plot
pdf("sharedGenes.enrchmnt.GO.mf.dotplot.pdf")
dotplot(enrchmnt.GO.mf,showCategory=5)
dev.off()
##
## You can also plot an interaction network of the enriched categories and their neighbours (or related categories)
## For that first calculate the pairwise similarities
enrchmnt.GO.mf.sm <- pairwise_termsim(enrchmnt.GO.mf)
## Now plot using the emapplot function and save it
pdf("sharedGenes.enrchmnt.GO.mf.emapplot.pdf")
emapplot(enrchmnt.GO.mf.sm,cex_category=0.5,cex_label_category = 0.5)
dev.off()
##
## Like you did for MF class, you can also test enrichment in the Biological process (BP) class
enrchmnt.GO.bp <- enrichGO(gene = my.gene.list$V1, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
## Note that here you are explicitly specifying the ontology class to test by specifying *ont = "BP"*
##
## How many biological processes are found significant?
length(which(enrchmnt.GO.bp@result$p.adjust < 0.05))
# [1] 270
##
## Checking out the top BP categories
(enrchmnt.GO.bp@result)[1:5,c(1:7,9)]
#                    ID                   Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue Count
# GO:0060562 GO:0060562 epithelial tube morphogenesis    27/463 394/29008 2.907943e-10 1.188767e-06 9.173794e-07    27
# GO:0001655 GO:0001655 urogenital system development    26/463 385/29008 8.525857e-10 1.742685e-06 1.344842e-06    26
# GO:0072001 GO:0072001      renal system development    23/463 342/29008 8.771659e-09 1.195285e-05 9.224092e-06    23
# GO:0003002 GO:0003002               regionalization    24/463 379/29008 1.334567e-08 1.278911e-05 9.869442e-06    24
# GO:0001822 GO:0001822            kidney development    22/463 327/29008 1.837151e-08 1.278911e-05 9.869442e-06    22
##
## Create and save the dotplot and ontology category interaction network plot
pdf("sharedGenes.enrchmnt.GO.bp.dotplot.pdf")
dotplot(enrchmnt.GO.bp,showCategory=5)
dev.off()
##
enrchmnt.GO.bp.sm <- pairwise_termsim(enrchmnt.GO.bp)
pdf("sharedGenes.enrchmnt.GO.bp.emapplot.pdf")

# cex options don't work so i commented out that line and made one without those options
#emapplot(enrchmnt.GO.bp.sm,cex_category=0.5,cex_label_category = 0.5)
emapplot(enrchmnt.GO.bp.sm)
dev.off()


##
## Now, you would perform enrichment for KEGG pathways
## In this case the gene identifiers need to be in Entrez ID format. You have gene symbols instead
## So, you would use the Id translator function to retrieve matching Entrez IDs for the symbols you have
my.entrez <- bitr(my.gene.list$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# 'select()' returned 1:1 mapping between keys and columns
# Warning message:
# In bitr(my.gene.list$V1, fromType = "SYMBOL", toType = "ENTREZID",  :
#   1.72% of input gene IDs are fail to map...
##
## Note the warning message. Its informing us that some of the symbols do not have an associated Entrez ID
##
## How many Entrez IDs you have retrieved?
dim(my.entrez)
# [1] 515   2
##
## Check few entries
head(my.entrez,3)
#          SYMBOL ENTREZID
# 1 0610040F04Rik    75394
# 3 1110002J07Rik    68488
# 4 1600002D24Rik    69776
##
## Now perform enrichment over KEGG pathways
enrchmnt.KEGG <- enrichKEGG(gene = my.entrez$ENTREZID, organism = 'mmu', minGSSize = 50)
##
## The way to access the results would be as before
dim(enrchmnt.KEGG@result)
# [1] 221 14
colnames(enrchmnt.KEGG@result)
# [1] "ID"          "Description" "GeneRatio"   "BgRatio"     "pvalue"      "p.adjust"    "qvalue"      "geneID"     
# [9] "Count"      
(enrchmnt.KEGG@result)[1:5,c(1:7,9)]
#                ID                Description GeneRatio  BgRatio       pvalue     p.adjust       qvalue Count
# mmu04015 mmu04015     Rap1 signaling pathway    18/184 214/9006 3.398395e-07 7.340533e-05 6.045566e-05    18
# mmu05224 mmu05224              Breast cancer    12/184 147/9006 4.484261e-05 4.843002e-03 3.988632e-03    12
# mmu04151 mmu04151 PI3K-Akt signaling pathway    18/184 359/9006 3.829257e-04 2.757065e-02 2.270682e-02    18
# mmu05221 mmu05221     Acute myeloid leukemia     7/184  70/9006 5.365152e-04 2.897182e-02 2.386081e-02     7
# mmu04010 mmu04010     MAPK signaling pathway    15/184 294/9006 9.919582e-04 3.734703e-02 3.075852e-02    15
##
## How many KEGG pathways with significant enrichment?
length(which(enrchmnt.KEGG@result$p.adjust < 0.05))
# [1] 17
##
## You can visualise the pathway that have come significant by using the *browseKEGG* function
## which will open the particular pathway on your web browser.
## And the genes which were present in your input would be highlighted by red colour font (as compared to standard black colour font) symbols
##
## You can see the KEGG pathway ID from the *enrchmnt.KEGG@result* print out
## The ID for the top pathway "Rap1 signaling" is mmu04015
## So -
browseKEGG(enrchmnt.KEGG,'mmu04015')
## That should open the mouse "Rap1 signaling" pathway with your input genes highlighted
## You can click the "Image file" button on the top of the web page and save the image with an appropriate file name
## 

