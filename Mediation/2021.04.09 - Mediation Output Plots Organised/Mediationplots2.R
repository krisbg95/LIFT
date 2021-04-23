## Libraries
library(tidyverse) 
library(readxl) 
library(colorspace) 
library (circlize)
library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(aspace) 
library(tiff) 
library(devEMF) 

## Function 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X,na="NA"))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
} 

all = read_excel_allsheets(paste0(getwd(),"/medranked_second_analysis.xlsx"))   
Regions_order <- read.csv(paste0(getwd(),"/Regions_order.csv"))


## Structural connectivity
mframe = all$sc_total[,1:3]
colnames(mframe) = c("Region1","Region2","Estimate")
mframe$Estimate <- as.numeric(mframe$Estimate)

# LH = t(unique(str_subset(t(mframe[,1:2]), "LH-"))) 
# RH = t(unique(str_subset(t(mframe[,1:2]), "RH-")))  

morder = as.data.frame(unique(c(mframe$Region1,mframe$Region2)));names(morder) = "Region" 
morder = merge.data.frame(Regions_order,morder,by = "Region");morder = morder[order(morder$Order),]; morder$Lobe <- as.factor(morder$Lobe) 
levels(morder$Lobe) = c(2:8);mcolour = as.integer(morder[,2]); names(mcolour) = morder[,1] 

  
col_fun = colorRamp2(c(-0.5,0.06), c("blue","red"), transparency = 0.7) 
# lgd_links = Legend(at = round(range(mframe$Estimate),2),col_fun = col_fun) ## legend options using ComplexHeatmap 
lgd_links = Legend(at = c(-0.5,0.06),col_fun = col_fun,legend_height = unit(4, "cm"),labels_gp = gpar(fontsize = 20)) ## legend options using ComplexHeatmap 
par(mfrow = c(1, 2),cex=1.2,font=2,mar = c(0, 0, 0, 0),oma = c(0,0,0,0), mgp = c(0, 0, 0),xpd = NA)
circos.par(canvas.xlim= c(0,0),canvas.ylim=c(-1.8,1.8),start.degree = 90)
chordDiagram(mframe, order = morder$Region,grid.col = mcolour,scale=T,link.decreasing = TRUE,link.sort= T,directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",diffHeight = mm_h(0),
             big.gap = 1,preAllocateTracks = 1,annotationTrack = "grid",annotationTrackHeight = 0.03,col=col_fun) 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important 
draw(lgd_links, x = unit(0.03, "npc") - unit(2, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear() 

# dev.print(tiff, "MediationSC.tiff", res=600, height=9.07, width=9.3, units="in",compression = "lzw")
# dev.print(emf, "MediationSC.emf",height=9.07, width=9.3) 


## Functional connectivity
mframe = all$rsfc_total[,1:3]
colnames(mframe) = c("Region1","Region2","Estimate")

# LH = t(unique(str_subset(t(mframe[,1:2]), "LH-"))) 
# RH = t(unique(str_subset(t(mframe[,1:2]), "RH-")))  

morder = as.data.frame(unique(c(mframe$Region1,mframe$Region2)));names(morder) = "Region" 
morder = merge.data.frame(Regions_order,morder,by = "Region");morder = morder[order(morder$Order),]; morder$Lobe <- as.factor(morder$Lobe) 
levels(morder$Lobe) = c(2:8);mcolour = as.integer(morder[,2]); names(mcolour) = morder[,1] 


col_fun = colorRamp2(range(mframe$Estimate), c("blue","red"), transparency = 0.7) 
lgd_links = Legend(at = round(c(-0.5,0.06),2),col_fun = col_fun) ## legend options using ComplexHeatmap 
# par(cex=0.65,font=1,mar = c(0, 0, 0, 0),oma = c(0,0,0,0), mgp = c(0, 0, 0),xpd = NA)
circos.par(canvas.xlim= c(0,0),canvas.ylim=c(-1.8,1.8),start.degree = 90)
chordDiagram(mframe,scale=T,order = morder$Region,grid.col = mcolour,link.decreasing = TRUE,link.sort= T,
             big.gap = 1,preAllocateTracks = 1,annotationTrack = "grid",annotationTrackHeight = 0.03,col=col_fun) 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important 
#draw(lgd_links, x = unit(0.05, "npc") - unit(2, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear() 

# dev.print(tiff, "MediationrsFC.tiff", res=600, height=9.07, width=9.3, units="in",compression = "lzw")
# dev.print(emf, "MediationrsFC.emf",height=9.07, width=9.3) 

dev.print(tiff, "Mediation2.tiff",res=600,height=10.5, width=23,units="in",compression = "lzw")
dev.print(emf, "Mediation2.emf",height=10.5, width=23)


