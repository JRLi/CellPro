rm(list=ls())
library(ComplexHeatmap)
library(circlize)
library(marray)
load('IMC_LUAD.RData')
outTiff = 'Cluster.tiff'
NA_less = 0 ######

sigs = CoxUniProximity169[CoxUniProximity169$Pu<0.01,'name']
dfx = IMC_proximity[row.names(IMC_clinical),sigs]
check = rowSums(is.na(dfx))
se = names(check[check <= NA_less])
length(se) # NA: 0: 218, 1: 259, 2: 319, 3: 347

info2 = IMC_clinical[se,]
dfx2 = dfx[se,]
info2[is.na(info2$Stage.Late1),'Stage.Late1']=1
info2$Stage = ifelse(info2$Stage.Late1==0, 'early','late')
info2$Smoking = ifelse(info2$Smoking==0, 'non-smoking','smoking')
info2$Sex = ifelse(info2$Sex.F==0, 'male','female')
info2$Age = ifelse(info2$Age.75==0, '≤75','>75')
info2$Stage = factor(info2$Stage, levels = c('early','late'))
info2$Smoking = factor(info2$Smoking, levels = c('non-smoking','smoking'))
info2$Sex = factor(info2$Sex, levels = c('female','male'))
info2$Age = factor(info2$Age, levels = c('≤75','>75'))

#heatmap_col = colorRamp2(c(min(dfx2), 0, max(dfx2)), c("#377EB8", "black", "#E41A1C"))
heatmap_col <- maPalette(low="green", mid= "black", high="red", k=40)
hist_col = c("lepidic"="lightskyblue1","papillary"="olivedrab2",'acinar'='khaki',
             'micropapillary'='darkorange1','solid'='#CC0000')
Stage_col = c("early" = "#CCEBC5", "late" = "#CC0000")
Sex_col = c('female'='lightskyblue1','male'='darkorange1')
Smoking_col = c('non-smoking'='#CCEBC5','smoking'='darkorange1')
ha = rowAnnotation(
  Subtype = info2$hist,
  Stage = info2$Stage,
  Sex = info2$Sex,
  Smoking = info2$Smoking,
  col = list(Subtype = hist_col,
             Stage = Stage_col,
             Sex = Sex_col,
             Smoking = Smoking_col),
  na_col = "white", border = TRUE,
  show_legend = c(TRUE, TRUE, TRUE, TRUE),
  show_annotation_name = F,
  annotation_legend_param = list(
    Subtype = list(title = "A. Subtype"),
    Stage = list(title = "B. Stage"),
    Sex = list(title = "C. Sex"),
    Smoking = list(title = "D. Smoking"))
)


htmp <- Heatmap(
  dfx2, 
  heatmap_legend_param = list(border = "black", title = "Proximity", # direction = "horizontal"
                              title_position = "leftcenter-rot",title_gp = gpar(fontsize = 11.5, fontface = "bold"), # title_position = "topcenter"
                              legend_width = unit(3.6, "cm")),
  col = heatmap_col,
  show_heatmap_legend = T,# column_title = "Cell Y", col = heatmap_col, 
  column_title_side = "bottom", column_title_gp = gpar(fontsize = 0),
  row_title_gp = gpar(fontsize = 6),row_km = 2, column_km = 3,
  show_column_names = T, show_row_names =  F,border = TRUE, 
  show_column_dend = T, show_row_dend = T, column_names_rot = 45,
  cluster_columns = T, cluster_rows = T, #row_split = sps,row_gap = unit(0, "mm"),
  column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 4),
  left_annotation = ha, clustering_method_rows = "complete" # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
)

# cairo_pdf(paste0('plot/Cluster_Z_',cells,surfix1,'.pdf'),width = wth,height = 15)
tiff(outTiff,width = 7,height = 6, units = 'in',res=300)
ht = draw(htmp, heatmap_legend_side="left", annotation_legend_side="left",
          legend_grouping = "original")
annotation_titles = c(Subtype = "A",
                      Stage = "B",
                      Sex = 'C',
                      Smoking = "D")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], y = unit(1, "npc") + unit(1, "mm"), just = "bottom")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
dev.off()

