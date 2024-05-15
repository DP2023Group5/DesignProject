library("limma")
library("AnnotationDbi")
library(clusterProfiler)
library(msigdbr)
# Homo sapiens annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(data.table)
genes_up = read.csv("~/Desktop/design_project/average_chrono_ORA/genes_up_average.csv")

# Split the 'EntrezGeneSymbol' column and keep only the first part
genes_up$EntrezGeneSymbol <- sapply(strsplit(genes_up$EntrezGeneSymbol, "\\|"), `[`, 1)

# This returns a named vector which we can convert to a data frame, where
# the keys (Ensembl IDs) are the names
entrez_vector <- mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    # The vector of gene identifiers we want to map
    keys = genes_up$EntrezGeneSymbol,
    # Replace with the type of gene identifiers in your data
    keytype = "SYMBOL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # In the case of 1:many mappings, return the
    # first one. This is default behavior!
    multiVals = "first"
)



proteins <- as.data.frame(fread("~/Desktop/design_project/belgium_somalogic_annotation_20211124.csv", na = ""))
universe <- proteins %>%
    filter(UniProt != "Q99LC4") %>%
    dplyr::select(EntrezGeneID) 
universe <- na.omit(universe)
universe <- universe$EntrezGeneID
universe <- unique(universe)

library(tibble)

go <- goana(entrez_vector, universe = universe, species = "Hs") %>% 
    as_tibble() %>%
    filter(P.DE <= 0.05)

library(GO.db)

go <- go %>%
    mutate(
        GOID = mapIds(
            GO.db, .$Term, "GOID", "TERM"
        ) %>% unname()
    ) %>%
    dplyr::select(GOID, everything()) %>%
    arrange(P.DE)

library(ggplot2)
library(forcats)

# Define custom colors for the legend
legend_colors <- c("BP" = "#1f77b4",  # Blue color for BP
                   "MF" = "#ff7f0e",  # Orange color for MF
                   "CC" = "#2ca02c")  # Green color for CC

# Create the ggplot object with data and aesthetics
graph <- ggplot(data = go, aes(x = fct_rev(fct_reorder(Term, P.DE)), y = P.DE, color = Ont)) +
    
    # Add the column geom
    geom_col(fill = "#69b3a2", alpha = 0.7, width = 0.5) +
    
    # Adjust the scale and labels
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    
    # Add a title and subtitle
    labs(title = "Distribution of P-values by Term",
         subtitle = "Ordered by P-value",
         x = NULL,
         y = "P-value",
         color = "Ontology") +
    
    # Apply a theme with custom appearance
    theme_minimal() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          plot.subtitle = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.text.x = element_blank(),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10)) +
    
    # Specify custom colors for the legend
    scale_color_manual(values = legend_colors) +
    
    # Fill legend symbols with custom colors
    guides(color = guide_legend(override.aes = list(fill = legend_colors,color = "black")))

# Print the graph
go_enrich <- enrichGO(gene = entrez_vector,
                      universe = universe,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.1)

library(enrichplot)
library(ggupset)

upsetplot(go_enrich)

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 20, 
        title = "GO Biological Pathways",
        font.size = 8)


dotplot(go_enrich)
emapplot(go_enrich)
cnetplot(go_enrich, categorySize="pvalue")

## KEGG 

kegg_organism = "hsa"
kk <- enrichKEGG(gene=entrez_vector, universe=universe,organism=kegg_organism, pvalueCutoff = 0.2)

barplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)


 cdotplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)

cnetplot(kk, categorySize="pvalue")



