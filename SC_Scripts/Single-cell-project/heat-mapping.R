#df_bulks<- read_csv("~/Desktop/Course_Materials/week6/single_cell/stages.csv")
df_bulks<- read_csv("~/Desktop/Course_Materials/week6/single_cell/sig_genes_table.csv")
df_bulks <- df_bulks[df_bulks$...1 %in% rownames(sce),]
stage_columns <- c("Early_stage", "Intermediate1_stage", "Intermediate2_stage", "Intermediate3_stage", "Late_stage")
df_bulks$Max_Stage <- apply(df_bulks[, stage_columns], 1, function(x) stage_columns[which.max(x)])

plotGroupedHeatmap(sce, 
                   features = df_bulks$...1,
                   group = "label",
                   block = "Sample", 
                   scale = TRUE, center = TRUE, 
                   zlim = c(-3, 3),
                   ) 

######################################

# Create a data frame for column annotation
annotation_df <- data.frame(Stage = df_bulks$Max_Stage)

rownames(annotation_df) <- df_bulks$...1
df_bulks$Max_Stage <- factor(df_bulks$Max_Stage)

colData(sce)
culster_names$label <- culster_names$X
culster_names$X <- NULL

####################################################################################
table(sce$timepoint)
timepoint_clusters <-table(sce$label, sce$timepoint)
timepoint_clusters[,1] <-  timepoint_clusters[,1] / 23429
timepoint_clusters[,2] <-  timepoint_clusters[,2] / 34562
timepoint_clusters[,3] <-  timepoint_clusters[,3] / 12064
max_column_names <- apply(timepoint_clusters, 1, function(x) names(x)[which.max(x)])
culster_names$time_point <- max_column_names

culster_names$Group_Time <- paste(culster_names$time_point, culster_names$Major_group)
library(RColorBrewer)

time_point <- c(
  "300" = "#654F6F",
  "400" = "#5C5D8D",
  "500" = "#99a1a6"
)
Major_group <- c(
  "body_wall" = "#654F6F",
  "cuticle" = "#5C5D8D",
  "somatic_gonad" = "#99a1a6",
  "excretory_duct" = "#A8C69F",
  "intestine" = "#F9FBB2",
  "excretory_canal" = "#1A281F",
  "basal_lamina" = "#E6EBE0",
  "pharynx" = "#c2948a"
  
)
Stage <- c(
  "Early_stage" = "#A8C69F",
  "Intermediate1_stage" = "#1A281F",
  "Intermediate2_stage" = "#E6EBE0",
  "Intermediate3_stage" = "#c2948a",
  "Late_stage" = "#654F6F"
)

annotation_colors <- list(time_point = time_point,
                          Major_group = Major_group,
                          Stage = Stage)
plotGroupedHeatmap(sce,
                   features = df_bulks$...1,
                   group = "label",
                   scale = TRUE,
                   center = TRUE,
                   zlim = c(-3, 3),
                   annotation_row = annotation_df[1],
                   annotation_col = culster_names[c(2,4)],
                   cutree_cols = 10,
                   cutree_rows = 9,
                   angle_col = 45,
                   fontsize = 11,
                   show_rownames = F,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   treeheight_row = 0,
                   annotation_colors = annotation_colors
                   
)


table(culster_names[3])
culster_names
annotation_df %>% dplyr::filter()
keep <- annotation_df[,1] != "Early_stage" & annotation_df[,1] != "Intermediate1_stage"
new_df <- annotation_df[keep, , drop = FALSE]
new_df_bulks <- df_bulks[keep, , drop = FALSE]




annotation_colors <- list(time_point = time_point)
plotGroupedHeatmap(sce,
                   features = new_df_bulks$...1,
                   group = "label",
                   scale = TRUE,
                   center = TRUE,
                   zlim = c(-3, 3),
                   annotation_row = new_df,
                   annotation_col = culster_names[2:3],
                   cutree_cols = 14,
                   cutree_rows = 9,
                   angle_col = 45,
                   fontsize = 11,
                   show_rownames = F,
                   color = colorRampPalette(rev(brewer.pal(11, "PRGn")))(100),
                   treeheight_row = 0
)

temp <- sce


culster_names$label <- row.names(culster_names)
merged_data <- left_join(as.data.frame(colData(temp)), culster_names, by = "label")
colData(temp) <- DataFrame(merged_data)


plotReducedDim(sce, 
               dimred = "UMAP",
               colour_by = "timepoint", 
               text_by = "timepoint") +
  ggtitle("Leiden k=25 clusters")

time_umap <- plotReducedDim(temp, 
               dimred = "UMAP",
               colour_by = "time_point") 

groups_umap <- plotReducedDim(temp, 
                       dimred = "UMAP",
                       colour_by = "Major_group") 

time_umap + groups_umap

annotation_df$Symbol <- row.names(annotation_df)
row_stages <- left_join(as.data.frame(rowData(temp)), annotation_df, by="Symbol")
rowData(temp) <- DataFrame(row_stages)

valid_rows <- !apply(is.na(rowData(temp)), 1, any)
temp <- temp[valid_rows, ]
temp$Sample



