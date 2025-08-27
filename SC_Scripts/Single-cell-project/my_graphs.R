time_point <- c(
  "300" = "#BC4749",
  "400" = "#80AB82",
  "500" = "#5C5D8D"
)
Major_group <- c(
  "body_wall" = "#80AB82",
  "cuticle" = "#5C5D8D",
  "somatic_gonad" = "#BC474972",
  "excretory_duct" = "#247BA0",
  "intestine" = "#99a1a6",
  "excretory_canal" = "#1A281F",
  "basal_lamina" = "#F06C9B",
  "pharynx" = "#c2948a"
  
)
Stage <- c(
  "Early_stage" = "#BC4749",
  "Intermediate1_stage" = "#1A281F",
  "Intermediate2_stage" = "#E6EBE0",
  "Intermediate3_stage" = "#80AB82",
  "Late_stage" = "#654F6F"
)


plotReducedDim(temp, 
               dimred = "UMAP",
               colour_by = "timepoint",)+
  theme(
    legend.text = element_text(size = 30),
    legend.key.size = unit(3, "lines")
  )



plotReducedDim(temp, 
               dimred = "UMAP",
               colour_by = "Major_group",
               text_size = 11) +
  scale_color_manual(values = Major_group)+
  theme(
    legend.text = element_text(size = 30),
    legend.key.size = unit(3, "lines")
  )




############################################


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
                   fontsize = 17,
                   show_rownames = F,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   treeheight_row = 0,
                   annotation_colors = annotation_colors
                   
)
########################################
plotReducedDim(sce, dimred = "PCA", colour_by = "timepoint")



