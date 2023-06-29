
pheatmap(exp.ALL[,order(all.cl, decreasing=TRUE)],
         annotation_col=ann,
         cluster_cols=FALSE, 
         cluster_rows=FALSE, 
         scale=heatmap.var, 
         breaks=bks,
         show_colnames=FALSE,
         show_rownames=show.rownames, 
         na_col = na_col, 
         color = color, 
         legend = legend)

pheatmap(mat,
         angle_col = c("270", "0", "45", "90", "315"),
         annotation_row = NA, 
         annotation_col = NA,
         annotation = NA, 
         annotation_colors = NA, 
         annotation_legend = TRUE,
         annotation_names_row = TRUE, 
         annotation_names_col = TRUE, 
         breaks = NA, 
         border_color = "grey60",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
         cellwidth = NA, 
         cellheight = NA,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         clustering_callback = identity2, 
         cutree_rows = NA, 
         cutree_cols = NA,
         display_numbers = F,
         drop_levels = TRUE,
         filename = NA, 
         fontsize = 10, 
         fontsize_row = fontsize, 
         fontsize_col = fontsize,
         fontsize_number = 0.8 * fontsize, 
         gaps_row = NULL, 
         gaps_col = NULL, 
         height = NA,
         scale = "none", 
         show_rownames = T, 
         show_colnames = T,
         silent = FALSE,
         treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows, 50, 0), 
         treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0), 
         kmeans_k = NA,
         labels_row = NULL,
         labels_col = NULL,
         legend = TRUE, 
         legend_breaks = NA,
         legend_labels = NA, 
         main = NA,
         na_col = "#DDDDDD",
         number_format = "%.2f", 
         number_color = "grey30", 
         width = NA, 
         ...)

         
         
         
         
         

         
         
         
         
         