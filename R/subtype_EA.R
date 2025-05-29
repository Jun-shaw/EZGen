
#' @export
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import GOplot
#' @import ggplot2
#' @import dplyr
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @references
#' Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5):284-287.
#' @examples # There is no example and please refer to vignette.

subtype_EA <- function(DEGs.res = DEGs.res,
                         data.base='GO',#c('GO','KEGG')
                         method   ='basic',#c('basic','circos','dendrogram','Taiji','bubble')
                         lfc.col  =c('#CE5C69','#5770A6'),
                         zsc.col  =c('#CE5C69', '#EDF5F7', '#5770A6'),
                         cat.col  =c('#CE5C69','#A281B1','#5770A6'))
{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOplot)
  library(ggplot2)
  library(dplyr)
  if (!dir.exists('EA_plot')) {
    dir.create('EA_plot')
  }
  DEGs.dt <- DEGs.res$DEGs
  colnames(DEGs.dt)[c(1,3)] <- c("SYMBOL","logFC") # Modify column names to prepare for the next merge
  DEGs.use <- DEGs.dt[,c("SYMBOL","logFC")] # Keep only the gene names and fold change columns

  Gene <- bitr(DEGs.use$SYMBOL,
               fromType = "SYMBOL", # Type of input data
               toType = c("ENTREZID"), # Type of data to convert to
               OrgDb = org.Hs.eg.db) # Species

 if (data.base=='GO') {
  # GO enrichment analysis --------------------------------------------------------------------
  GO <- enrichGO(gene = Gene$ENTREZID, # Input gene "ENTREZID"
                 OrgDb = org.Hs.eg.db,# Annotation information
                 keyType = "ENTREZID",
                 ont = "all",     # Optional items BP/CC/MF
                 pAdjustMethod = "BH", # Method for p-value adjustment
                 pvalueCutoff = 1,   # Threshold for p-value
                 qvalueCutoff = 1, # Threshold for q-value
                 minGSSize = 5,
                 maxGSSize = 5000,
                 readable = TRUE)   # Whether to convert entrez id to symbol
  GO_result <- as.data.frame(GO) # Convert result

  if (method=='basic') {
    # Beautification version --------------------------------------------------------------------
    # Select the top 10 significant items from each category
    BP_top10 <- GO_result %>%
      filter(ONTOLOGY == "BP") %>%  # Filter sub-data where ONTOLOGY is BP
      arrange(p.adjust) %>%        # Arrange in ascending order by p.adjust
      head(5)                      # Extract top 10 rows

    CC_top10 <- GO_result %>%
      filter(ONTOLOGY == "CC") %>%  # Filter sub-data where ONTOLOGY is CC
      arrange(p.adjust) %>%        # Arrange in ascending order by p.adjust
      head(5)                      # Extract top 10 rows

    MF_top10 <- GO_result %>%
      filter(ONTOLOGY == "MF") %>%  # Filter sub-data where ONTOLOGY is MF
      arrange(p.adjust) %>%        # Arrange in ascending order by p.adjust
      head(5)                      # Extract top 10 rows

    merge_data <- rbind(BP_top10,CC_top10,MF_top10) # Merge the three datasets
    merge_data$ONTOLOGY <- factor(merge_data$ONTOLOGY, levels = c("BP", "CC","MF"))
    merge_data$logPvalue <- (-log(merge_data$p.adjust)) # Data processing

    colnames(merge_data) # Check column names
    merge_data$Description <- factor(merge_data$Description, levels = unique(merge_data$Description[order(merge_data$ONTOLOGY)]))

    # Horizontal bar plot - plot based on p-value #
    pdf(file=paste0('EA_plot/',"GO_bubble.pdf"),width = 10,height = 8)
    p1 <- ggplot(merge_data,
                 aes(x=Description,y=Count, fill=pvalue)) +  # Define x and y axes; fill color according to ONTOLOGY
      geom_bar(stat="identity", width=0.8) +  # Width of the bar plot
      scale_fill_gradient(low = lfc.col[1],high =lfc.col[2] ) + # Fill color of bar plot
      xlab("GO term") + # x-axis label
      ylab("Gene Number") +  # y-axis label
      labs(title = "GO Terms Enrich")+ # Set title
      theme_bw() +
      theme(axis.text.x=element_text(family="sans", color="black",angle = 70,vjust = 1, hjust = 1 ,size=11)) # Font style, color, and x-axis angle adjustment
    # Add grouping boxes based on GO enrichment analysis classification information #
    p1 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')
    print(p1)
    dev.off()
    p1
  }
  else
  {
    GO_result <- GO_result[(GO_result$pvalue < 0.05 & GO_result$p.adjust < 0.05),]
    go_result <- data.frame(Category = GO_result$ONTOLOGY,
                            ID = GO_result$ID, # Contains pathway ID of enrichment analysis results
                            Term = GO_result$Description, # Contains pathway descriptions of enrichment analysis results
                            Genes = gsub("/", ", ", GO_result$geneID), # Contains gene IDs from enrichment analysis results, separated by commas if multiple
                            adj_pval = GO_result$p.adjust) # Contains adjusted p-values from enrichment analysis results
    # Create a data frame containing gene IDs and differential expression values
    genelist <- data.frame(ID = DEGs.use$SYMBOL, logFC = DEGs.use$logFC)
    genelist <- genelist %>% distinct(ID, .keep_all = TRUE)
    # Set the row names of the genelist data frame to ID
    row.names(genelist)=genelist[,1]
    # Generate a visual chart of the pathway enrichment analysis results
    library(GOplot)
    circ <- circle_dat(go_result, genelist)
    head(circ)

    # Visualization preparation and plotting ----------------------------------------------------------------
    nrow(go_result) # Check the total number of pathways
    term_num = 3 # Set the number of GO entries to display in the graph

    nrow(genelist) # Check the total number of genes
    gene_num = nrow(genelist)  # Set the upper limit for the number of genes to display

    # Generate a chord diagram of the pathway enrichment analysis results
    chord <- chord_dat(circ, # A data frame containing pathway enrichment analysis results
                       genelist[1:gene_num,], # A data frame containing gene IDs and differential expression values, including only the top gene_num genes
                       go_result$Term[1:term_num] # A vector containing pathway descriptions, including only the top term_num pathways
    )

    if (method=='circos') {
      # I. GO enrichment circle diagram
      pdf(file=paste0('EA_plot/',"GO_circos.pdf"),width = 12,height = 12)
      p1 <- GOChord(chord,
              title = "GOcircos", # Title
              space = 0.01, # Set the distance between genes to 0.01
              gene.order = 'logFC',    # Sort genes based on logFC values
              gene.size = 2, # Size of gene labels
              nlfc = 1, # Define the number of logFC columns (default = 1)
              lfc.col=zsc.col, # Specified logFC fill colors: c (color for low values, color for midpoint, color for high values)
              lfc.min = -3, # Specify the minimum value for logFC ratio (default = -3)
              lfc.max = 3, # Specify the maximum value for logFC ratio (default = 3)
              gene.space = 0.2,       # Distance of genes from the circle
              border.size = 0.2, # Define the size of functional area borders
              process.label = 12, # Size of legend
              limit = c(1,5)
              # Controls the number of GO items related to genes in the chord diagram, ensuring clear and readable visualization.
              # The first value is 1, indicating that each gene should be assigned at least one GO item.
              # The second value is 5, indicating that each displayed GO item should be associated with at least 5 genes.
      )
      print(p1)
      dev.off()
      p1
    }
    if (method=='dendrogram') {
      pdf(file=paste0('EA_plot/',"GO_dendrogram.pdf"),width = 10,height = 8)
      p1 <-GOCluster(circ,
                metric = "euclidean", # Set distance metric method to Euclidean distance
                clust.by='logFC', # Cluster based on logFC values
                lfc.col=zsc.col, # Specified logFC fill colors: c (color for low values, color for midpoint, color for high values)
                as.character(go_result[1:term_num,3]) # Select the top termNum GO terms from go_result
      )
      print(p1)
      dev.off()
      p1
    }
    if (method=='Taiji') {
      pdf(file=paste0('EA_plot/',"GO_Taiji.pdf"),width = 16,height = 8)
      p1 <-GOCircle(circ,nsub = 8, lfc.col = lfc.col,zsc.col = zsc.col)
      go_result$ID[1:20]
      print(p1)
      dev.off()
      p1
    }
    if (method=='bubble') {
      pdf(file=paste0('EA_plot/',"GO_bubble.pdf"),width = 10,height = 8)
      p1 <-GOBubble(circ,
               labels = 20, # The threshold for displaying labels in the graph (GO: ...), labels will be marked for terms with -log(p.adj) > 6, the larger the value, the fewer displays
               display= 'multiple', # Display bubbles for multiple variables on the same graph
               colour = cat.col, # Set colors
               bg.col = T # Whether to add background color
      )
      print(p1)
      dev.off()
      p1
    }
  }
 }
  if (data.base=='KEGG') {
  KEGG <- enrichKEGG(gene = Gene$ENTREZID,
                     organism = "hsa",
                     keyType = "kegg", #KEGG数据库
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)

  df <- as.data.frame(KEGG)
  df_sorted <- df[order(df$p.adjust), ][1:10, ]
  colnames(df_sorted)
  df_sorted$Description
  df_sorted$Description <- gsub(" - Homo sapiens \\(human\\)", "", df_sorted$Description)
  pdf(file=paste0('EA_plot/',"KEGG_bubble.pdf"),width = 10,height = 8)
  p1 <-ggplot(df_sorted, aes(x=GeneRatio, y=Description, size=-log10(p.adjust), color=p.adjust)) +
    geom_point(alpha=0.7) +
    scale_color_gradient(low=lfc.col[1], high=lfc.col[2]) +
    labs(title="The most enrichment KEGG pathway", x="GeneRatio", y="KEGG pathway", size="-log10(P-value)", color="P-value") +
    theme_minimal() +
    theme(axis.text.y = element_text(face = 'bold',size = 12),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 18, face = "bold"))
  print(p1)
  dev.off()
  p1
  }
}





