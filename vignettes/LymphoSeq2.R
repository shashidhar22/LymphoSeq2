## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data---------------------------------------------------------------------
system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2") #For the TCRB files
system.file("extdata", "IGH_sequencing", package = "LymphoSeq2") #For the IGH files.

## ----readImmuno---------------------------------------------------------------
library(LymphoSeq2)
library(tidyverse)
study_files <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2") 
study_table <- LymphoSeq2::readImmunoSeq(study_files)

## ----showStudyTable-----------------------------------------------------------
study_table

## ----getSamplesNames----------------------------------------------------------
study_table %>% 
dplyr::pull(repertoire_id) %>% 
unique()

## ----subSample----------------------------------------------------------------
TRB_Unsorted_0 <- study_table %>% filter(repertoire_id == "TRB_Unsorted_0")
TRB_Unsorted_0

## ----subSamplePattern---------------------------------------------------------
CMV <- study_table %>% 
       dplyr::filter(str_detect(repertoire_id, "CMV"))
CMV

## ----loadMeta-----------------------------------------------------------------
TCRB_metadata <- readr::read_csv(system.file("extdata", "TCRB_metadata.csv", package = "LymphoSeq2"))
TCRB_metadata

## ----mergeMeta----------------------------------------------------------------
study_table <- dplyr::left_join(study_table, TCRB_metadata, by = c("repertoire_id" = "samples"))
study_table

## ----metaFilter---------------------------------------------------------------
unsorted_300 <- study_table %>% 
                dplyr::filter(day > 300 & phenotype == "Unsorted")
unsorted_300

## ----readAmino----------------------------------------------------------------
aa_table <- LymphoSeq2::productiveSeq(study_table = study_table, aggregate = "junction_aa", 
                          prevalence = FALSE)
  aa_table

## ----readNucleotide-----------------------------------------------------------
nuc_table <- LymphoSeq2::productiveSeq(study_table = study_table, aggregate = "junction", 
                          prevalence = FALSE)
nuc_table

## ----showUnprod---------------------------------------------------------------
study_table %>% 
dplyr::filter(repertoire_id == "TRB_Unsorted_0") 

## ----showProdAA---------------------------------------------------------------
aa_table %>% 
dplyr::filter(repertoire_id == "TRB_Unsorted_0")

## ----showSummary--------------------------------------------------------------
LymphoSeq2::clonality(study_table = study_table)

## ----readIGH------------------------------------------------------------------
IGH_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq2")
IGH_table <- LymphoSeq2::readImmunoSeq(path = IGH_path)
LymphoSeq2::clonalRelatedness(study_table = IGH_table, editDistance = 10)

## ----getIGHProd, fig.width = 7, fig.height = 10, comment = ""-----------------
nuc_IGH_table <- LymphoSeq2::productiveSeq(study_table = IGH_table, aggregate = "junction")
LymphoSeq2::phyloTree(study_table = nuc_IGH_table, 
                      repertoire_ids = "IGH_MVQ92552A_BL", 
                      type = "junction", 
                      layout = "rectangular")

## ----alignSeq, results='asis'-------------------------------------------------
alignment <- LymphoSeq2::alignSeq(study_table = nuc_IGH_table, 
                                  repertoire_ids = "IGH_MVQ92552A_BL", 
                                  type = "junction_aa", 
                                  method = "ClustalW")
msa::msaPrettyPrint(alignment, output = "asis")

## ----searchSeq----------------------------------------------------------------
LymphoSeq2::searchSeq(study_table = aa_table, 
                      sequence = "CASSPVSNEQFF", 
                      seq_type = "junction_aa", 
                      match = "global", 
                      edit_distance = 0)

## ----searchPublished----------------------------------------------------------
LymphoSeq2::searchPublished(study_table = aa_table) %>% 
dplyr::filter(!is.na(PMID))

## ----plotLorenz---------------------------------------------------------------
samples <- aa_table %>% 
           dplyr::pull(repertoire_id) %>% unique()
LymphoSeq2::lorenzCurve(repertoire_ids = samples, study_table = aa_table)

## ----topSeqPlot---------------------------------------------------------------
LymphoSeq2::topSeqsPlot(study_table = aa_table, top = 10)

## ----bhattacharyya------------------------------------------------------------
bhattacharyya_matrix <- LymphoSeq2::scoringMatrix(aa_table, mode = "Bhattacharyya")
LymphoSeq2::pairwisePlot(bhattacharyya_matrix)

## ----commonSeq----------------------------------------------------------------
common <- LymphoSeq2::commonSeqs(study_table = aa_table, 
                                 repertoire_ids =  c("TRB_Unsorted_0", "TRB_Unsorted_32"))
common

## ----commonVenn, fig.height=5, fig.width=5------------------------------------
LymphoSeq2::commonSeqsVenn(repertoire_ids = c("TRB_Unsorted_32", "TRB_Unsorted_83"), 
                           productive_aa = aa_table)

## ----commonVennT, fig.height=5, fig.width=5-----------------------------------
LymphoSeq2::commonSeqsVenn(repertoire_ids = c("TRB_Unsorted_0", "TRB_Unsorted_32", "TRB_Unsorted_83"), 
                           productive_aa = aa_table)

## ----commonSeqPlot------------------------------------------------------------
LymphoSeq2::commonSeqsPlot("TRB_Unsorted_32", "TRB_Unsorted_83", 
                           productive_aa = aa_table, show = "common")

## ----commonBar----------------------------------------------------------------
LymphoSeq2::commonSeqsBar(productive_aa = aa_table, 
                          repertoire_ids = c("TRB_CD4_949", "TRB_CD8_949", 
                                             "TRB_Unsorted_949", "TRB_Unsorted_1320"), 
                          color_sample = "TRB_CD8_949",
                          labels = "no")

## ----diffAbund----------------------------------------------------------------
LymphoSeq2::differentialAbundance(study_table = aa_table, 
                                  repertoire_ids =c("TRB_Unsorted_949", 
                                                    "TRB_Unsorted_1320"), 
                                  type = "junction_aa", q = 0.01)

## ----uniqueSeq----------------------------------------------------------------
unique_seqs <- LymphoSeq2::uniqueSeqs(productive_table = aa_table)
unique_seqs

## ----seqMatrix----------------------------------------------------------------
sequence_matrix <- LymphoSeq2::seqMatrix(productive_aa = aa_table, sequences = unique_seqs$junction_aa, )
sequence_matrix

## ----topFreq------------------------------------------------------------------
top_freq <- LymphoSeq2::topFreq(productive_table = aa_table, frequency = 0.001)
top_freq

## ----mergeData----------------------------------------------------------------
top_freq_matrix <- dplyr::full_join(top_freq, sequence_matrix)
top_freq_matrix

## ----cloneTrack---------------------------------------------------------------
ctable <- LymphoSeq2::cloneTrack(study_table = aa_table,
                                 sample_list = c("TRB_CD8_949", "TRB_CD8_CMV_369"))
LymphoSeq2::plotTrack(ctable)

## ----plotCloneSet-------------------------------------------------------------
ttable <- LymphoSeq2::topSeqs(aa_table, top = 10)
ctable <- LymphoSeq2::cloneTrack(ttable)
LymphoSeq2::plotTrack(ctable, alist = c("CASSESAGSTGELFF", "CASSLAGDSQETQYF")) + theme(legend.position = "bottom")

## ----plotSingular-------------------------------------------------------------
lalluvial <- LymphoSeq2::plotTrackSingular(ctable)
lalluvial[[1]]

## ----geneFreq-----------------------------------------------------------------
vGenes <- LymphoSeq2::geneFreq(productive_nt = nuc_table, locus = "V", family = TRUE)
vGenes

## ----chordDiagram-------------------------------------------------------------
top_seqs <- LymphoSeq2::topSeqs(productive_table = nuc_table, top = 1)
LymphoSeq2::chordDiagramVDJ(study_table = top_seqs, 
                            association = "VJ", 
                            colors = c("darkred", "navyblue"))

## ----wordcloud----------------------------------------------------------------
vGenes <- LymphoSeq2::geneFreq(productive_nt = nuc_table, locus = "V", family = TRUE)
library(RColorBrewer)
library(grDevices)
RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
library(wordcloud)
vtable <- vGenes %>% dplyr::filter(repertoire_id == "TRB_Unsorted_83")
wordcloud::wordcloud(words = vtable$gene_name, 
                     freq = vtable$gene_frequency, 
                     colors = RedBlue)

## ----heatmap------------------------------------------------------------------
vGenes <- LymphoSeq2::geneFreq(nuc_table, locus = "V", family = TRUE) %>% 
           tidyr::pivot_wider(id_cols = gene_name, 
                              names_from = repertoire_id, 
                              values_from = gene_frequency, 
                              values_fn = sum,
                              values_fill = 0)
gene_names <- vGenes %>% 
              dplyr::pull(gene_name)
vGenes <- vGenes %>% 
          dplyr::select(-gene_name) %>% 
          as.matrix()
rownames(vGenes) <- gene_names
pheatmap::pheatmap(vGenes, scale = "row")

## ----geneFreqPlot-------------------------------------------------------------
vGenes <- LymphoSeq2::geneFreq(productive_nt = nuc_table, locus = "V", family = TRUE)
multicolors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Set1")))(28)
ggplot2::ggplot(vGenes, aes(x = repertoire_id, y = gene_frequency, fill = gene_name)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_y_continuous(expand = c(0, 0)) + 
  ggplot2::guides(fill = guide_legend(ncol = 2)) +
  ggplot2::scale_fill_manual(values = multicolors) + 
  ggplot2::labs(y = "Frequency (%)", x = "", fill = "") +
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## ----searchSequences----------------------------------------------------------
LymphoSeq2::searchSeq(study_table = aa_table, sequence = "CASSESAGSTGELFF", seq_type = "junction_aa")

## ----removeSeq----------------------------------------------------------------
cleansed <- LymphoSeq2::removeSeq(study_table = aa_table, sequence = "CASSESAGSTGELFF")
LymphoSeq2::searchSeq(study_table = cleansed, sequence = "CASSESAGSTGELFF", seq_type = "junction_aa")

## ----rarefation, warning=FALSE------------------------------------------------
LymphoSeq2::plotRarefactionCurve(study_table = aa_table)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

