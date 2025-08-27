sce <- readRDS("/home/participant/Course_Materials/week6/single_cell/sce-corrected-clustered")
markers <- scoreMarkers(sce, 
                        groups = sce$label, 
                        block = sce$Sample)


# loop through list of marker genes and extract top-ranked gene names
top_markers_all <- lapply(markers, function(x){
  x %>% 
    as.data.frame() %>% 
    dplyr::filter(rank.logFC.cohen < 10) %>% 
    rownames()
})

# examining this list reveals several known markers of immune cells


sce$timepoint
rowData(sce)


top_markers_all[10]
# violin plot
Neuronal <- c("unc-119", "rab-3", "unc-4", "eat-4", "dat-1", "tph-1", "cho-1")
Muscle <- c("myo-3", "unc-54", "hlh-1")
Pharynx <- c("myo-2", "pha-4")
Embryo_Developmental <- c("pal-1", "mex-5", "skn-1")

plotExpression(sce, x = "label", features = Embryo_Developmental)
plotReducedDim(sce, 
               dimred = "UMAP",
               by_exprs_values = "logcounts",
               colour_by = "mex-5",
               text_by = "label")

known_genes <- c(
  "unc-119", # Neuronal
  "myo-3", # Muscle
  "myo-2", # Pharynx
  "skn-1", # Early Embryo/Developmental Markers
  "pie-1",  # Germline Markers
  "elt-2", # entestinal Markers
  "dpy-7"#Hypodermal Markers
)

bulk_genes <- c("act-5","aqp-4","bre-1","cah-4","cki-1","cnb-1",
                "col-117","col-130","cpn-3","deb-1","dpy-14","dsl-3",
                "eps-8","fkb-3","gln-3","gpd-3","gta-1","hil-4","hot-5",
                "ifb-1","ifb-2","lbp-6","lec-6","lec-9","mlc-2","mup-2",
                "nlp-24","nlp-29","nlt-1","pqn-94","rpl-3","rpl-7",
                "rpl-9","rpl-11.2","rpl-15","rpl-17","rpl-18","rpl-22",
                "rps-3","rps-4","rps-5","rps-7","rps-11","rps-19",
                "rps-21","sqt-3","sup-12","sym-1","tnc-2","tni-1","tts-2"
                ,"unc-22","vha-12","wrt-10","C35C5.10","F07A11.4","F08G5.6"
                ,"F09B12.3","F10G8.8","clec-196","ccch-1","F53F8.4",
                "dct-18","rack-1","K07C5.9","pitp-1","M176.5","R07E3.1"
                ,"ttr-15","dod-6","hsp-12.1","sipa-1","mlt-11","W04A8.4"
                ,"W06F12.2","Y37D8A.6","Y37D8A.16","Y57A10A.23","glb-1"
                ,"ZK1053.4","C08F1.10","C16D9.1","C23H5.8","C37A2.7","test-1"
                ,"acdh-1","F14B8.6","F28B4.3","drd-2","F56C9.7","F56F10.1",
                "K11H12.7","M60.4","pigv-1","T09B4.5","T19D12.1",
                "T23E7.2","Y110A2AL.4","ttr-36", "M60.2") ##could not find M60.2

bulk_genes_50 <- c("act-5","cah-4","cki-1","cnb-1","dpy-14","eps-8","fkb-3","gln-3","hil-4","ifb-2","lec-9","mlc-2","mup-2","nlp-24","nlp-29","rpl-7","rpl-11.2","rpl-17","rpl-18","rps-4","rps-7","rps-11","rps-19","sqt-3","sup-12","sym-1","tni-1","tts-2","vha-12","C35C5.10","clec-196","ccch-1","dct-18","rack-1","M176.5","R07E3.1","dod-6","hsp-12.1","W04A8.4","Y37D8A.6","Y37D8A.16","C37A2.7","test-1","K11H12.7","M60.2","M60.4","pigv-1","T09B4.5","T19D12.1","Y110A2AL.4")

bulk_genes_50 <- bulk_genes_50[bulk_genes_50  %in% rownames(sce)]

unique(colData(sce)$Sample)


plotGroupedHeatmap(sce, 
                   features = bulk_genes_50,
                   group = "label",
                   block = "Sample", 
                   scale = TRUE, center = TRUE, 
                   zlim = c(-3, 3)
                   )

culster_names <- read.csv("/home/participant/Course_Materials/week6/single_cell/cluster_names.csv")

