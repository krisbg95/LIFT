## Libraries
library(tidyverse) 
library(readxl) 
library(colorspace) 
library (circlize)
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(aspace) 
library(tiff) 
library(devEMF) 

## Structural connectivity
mnames <- read_csv(paste0(getwd(),"/Outputs/Region_names.csv"),col_names = F)   
allc <- read_csv(paste0(getwd(),"/Outputs/SCtotal.csv"),col_names = F)   
allce <- read_csv(paste0(getwd(),"/Outputs/Significant Mediators/SCtotal.csv"),col_names=F) 
allce[1,2:ncol(allce)] <- t(allc[2:nrow(allc),])  
rm(allc)
allce[1,1] <- "Connection" 

connames <- as.character(allce[1,2:ncol(allce)])   

mframe <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Region1", "Region2", "Estimate"))

for(i in 1:length(connames)){
  before = sub("-.*", "", connames[i]) 
  after = sub(".*-", "", connames[i]) 
  before = mnames[as.integer(before),1] 
  after = mnames[as.integer(after),1] 
  n = 1 + i
  estimate = allce[3,n] 
  mframe[i,] <- cbind(before,after,estimate)
} 

mframe$Estimate <- as.numeric(mframe$Estimate)
LH = t(unique(str_subset(t(mframe[,1:2]), "LH-"))) 
RH = t(unique(str_subset(t(mframe[,1:2]), "RH-"))) 

col_fun = colorRamp2(range(mframe$Estimate), c("blue","red"), transparency = 0.7) 
# lgd_links = Legend(at = round(range(mframe$Estimate),2),col_fun = col_fun) ## legend options using ComplexHeatmap 
lgd_links = Legend(at = c(-0.5,0.04),col_fun = col_fun) ## legend options using ComplexHeatmap 
par(mfrow = c(1, 2),cex=0.65,font=2,mar = c(0, 1, 0, 0),oma = c(0,0,0,0), mgp = c(0, 0, 0),xpd = NA)
circos.par(canvas.xlim= c(0,0),canvas.ylim=c(-1.3,1.3),start.degree = 90)
chordDiagram(mframe, order = c(RH,LH),grid.col = "black",scale=T,link.decreasing = TRUE,link.sort= T,directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",diffHeight = mm_h(0),
             big.gap = 1,preAllocateTracks = 1,annotationTrack = "grid",annotationTrackHeight = 0.02,col=col_fun) 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important 
draw(lgd_links, x = unit(0.05, "npc") - unit(2, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear() 

# dev.print(tiff, "MediationSC.tiff", res=600, height=9.07, width=9.3, units="in",compression = "lzw")
# dev.print(emf, "MediationSC.emf",height=9.07, width=9.3) 


## Functional connectivity
allc <- read_csv(paste0(getwd(),"/Outputs/rsFCtotal.csv"),col_names = F)   
allce <- read_csv(paste0(getwd(),"/Outputs/Significant Mediators/rsFCsigtotal.csv"),col_names=F) 
allce[1,2:ncol(allce)] <- t(allc[2:nrow(allc),])  
rm(allc)
allce[1,1] <- "Connection" 

connames <- as.character(allce[1,2:ncol(allce)])   

mframe <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Region1", "Region2", "Estimate"))

for(i in 1:length(connames)){
  before = sub("-.*", "", connames[i]) 
  after = sub(".*-", "", connames[i]) 
  before = mnames[as.integer(before),1] 
  after = mnames[as.integer(after),1] 
  n = 1 + i
  estimate = allce[3,n] 
  mframe[i,] <- cbind(before,after,estimate)
} 

mframe$Estimate <- as.numeric(mframe$Estimate)
LH = t(unique(str_subset(t(mframe[,1:2]), "LH-"))) 
RH = t(unique(str_subset(t(mframe[,1:2]), "RH-"))) 
col_fun = colorRamp2(range(mframe$Estimate), c("blue","red"), transparency = 0.7) 
lgd_links = Legend(at = round(range(mframe$Estimate),2),col_fun = col_fun) ## legend options using ComplexHeatmap 
# par(cex=0.65,font=1,mar = c(0, 0, 0, 0),oma = c(0,0,0,0), mgp = c(0, 0, 0),xpd = NA)
circos.par(canvas.xlim= c(0,0),canvas.ylim=c(-1.3,1.3),start.degree = 90)
chordDiagram(mframe,scale=T,grid.col = "black",order = c(RH,LH),link.decreasing = TRUE,link.sort= T,
             big.gap = 1,preAllocateTracks = 1,annotationTrack = "grid",annotationTrackHeight = 0.02,col=col_fun) 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important 
#draw(lgd_links, x = unit(0.05, "npc") - unit(2, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear() 

# dev.print(tiff, "MediationrsFC.tiff", res=600, height=9.07, width=9.3, units="in",compression = "lzw")
# dev.print(emf, "MediationrsFC.emf",height=9.07, width=9.3) 

dev.print(tiff, "Mediation.tiff",res=600,height=11, width=18.8,units="in",compression = "lzw")
dev.print(emf, "Mediation.emf",height=11, width=18.8)