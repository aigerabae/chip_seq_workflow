# chip_seq_workflow
msc bioinformatics, m6, chip-seq tutorial  


The purpose of ChIP-seq is to find regions of DNA that are associated with certain proteins. For that, proteins are immunoprecipitated with specific antibodies along with the DNA fragments bound by it. That mixture is then cleaned of antibodies and proteins and sequenced. The resulting sequences are aligned to a reference genome and peaks are called (ie the areas where there are substantially more fragments) which indicate that that are was targeted by that protein that was immunoprecipitated. After QC, the genes in those regions are identified and GO ontology or KEGG analysis can be performed to analyze what pathways are affected or linked to that protein.  

tutorial.md highlights command line work done to perform alignment, QC, and peak calling. go_pathway.R highlights the work done in R that does enrichment analysis of found genes
