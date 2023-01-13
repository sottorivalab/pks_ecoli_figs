####################################################
####################################################
###########FigureS1 phylogenetic tree################
####################################################
####################################################


#####Step1##########################################
#####run BEAST with AlvesParameter##################
#####the name of Summary Tree---C***_WGS_clean.summary.tree
#####the path of Summary Tree---'/Users/bchen/Desktop/Projects/PKSFiles/Data/All_summary_tree/'


#####Step2##########################################
#####plot megerd BEAST trees
#####Load packages
#if(!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")
#BiocManager::install("treeio")
#BiocManager::install("phylobase")
library("ggtree")
library("treeio")
library(ape)
library(ggplot2)
library(phylobase)
library(dplyr)
####################################################


######define own operator not in ###################
`%ni%` <- Negate(`%in%`) 
####################################################


#######define multiplot function ###################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
      
    }
  }
}
####################################################


#######ID of EPICC Patients#########################

Pids<-c('C516','C519','C524','C527',
        'C528','C530','C531',
        'C537','C543','C544','C547','C548',
        'C549','C550','C551','C552',
        'C555','C560','C561','C562')

Pids<-as.vector(t(matrix(Pids,nrow=5,ncol=4)))
######name of outgroup for each pid ###################
outgroupsamples<-c('Z1_B1','E1_B1','E1_B1','E1_B1',
                   'Z1_B1','E1_B1','Z1_B1',
                   'Z1_B1','E1_B1','Z1_B1','Z1_B1','Z1_B1',
                   'Z1_B1','Z1_B1','E1_B1','Z1_B1',
                   'Z1_B1','E1_B1','E1_B1','E1_B1')
outgroupsamples<-as.vector(t(matrix(outgroupsamples,nrow=5,ncol=4)))
####################################################





######list of potential normals
Pnormal<-list(
  'C519'=c('A1_G10','A1_G1','A1_G5','B1_G2','B1_G3',
           'D1_G2','D1_G9','E1_G1'),
  'C524'=c('A1_G3','A1_G6'),
  'C527'=c('A1_G10','A1_G8','C1_G7','C1_G9','D1_B1'),
  'C528'=c('B1_G6','D1_G6','D1_G8','D1_G9'),
  'C530'=c('B1_G2','B1_G7','C1_G10','C1_G7','E1_G1',
           'E1_G2'),
  'C531'=c('D1_G1','D1_G2','D1_G7'),
  'C537'=c('B1_G3','B1_G4','D1_G7','D1_G10','E1_G3'),
  'C543'=c('A1_G2','A1_G10','C1_G1','C1_G8','D1_G9',
           'D1_G10'),
  'C544'=c('B1_G1','B1_G2','C1_G3','E1_G1','E1_G3'),
  'C547'=c('A1_G4','B1_G2','B1_G3','B1_G4','C1_G7',
           'D1_G6','D1_G8','E1_G1'),
  'C549'=c('A1_G2','A1_G8','C1_G4','C1_G8'),
  'C550'=c('A1_G2','B1_G2','B1_G3','C1_G3','C1_G4',
           'D1_G2','D1_G5'),
  'C551'=c('D1_G5'),
  'C552'=c('A1_G1','B1_G10','B1_G8','B1_G9','C1_G1',
           'C1_G8','E1_G3'),
  'C555'=c('B1_G10','C1_G1','C1_G7','C1_G9'),
  'C560'=c('A1_G8','D1_G5','D1_G9'),
  'C562'=c('B1_G3','B1_G4','C1_G6','D1_G2','E1_G2')
)


######name of outgroup for each pid ###################
outgroupsamples<-c('Z1_B1','E1_B1','E1_B1','E1_B1','E1_B1','E1_B2',
                   'E1_B1','Z1_B1','E1_B1','Z1_B1','Z1_B1','Z1_B1',
                   'Z1_B1','E1_B1','E1_B1','Z1_B1','E1_B1','Z1_B1',
                   'Z1_B1','Z1_B1','Z1_B1','Z1_B1','E1_B1','Z1_B1',
                   'Z1_B1','Z1_B1','Z1_B1','E1_B1','E1_B1','E1_B1')
outgroupsamples<-as.vector(t(matrix(outgroupsamples,nrow=6,ncol=5)))
####################################################


######save tree to mylist###################
path='/Users/bchen/Desktop/Projects/PKSFiles/Data/All_summary_tree/All_summary_tree/'
mylist<-list()
for(i in (1:length(Pids))){
  pid=Pids[i]
  pidPnormal=Pnormal[[pid]]
  outgroupsample=outgroupsamples[i]
  file=paste0(path,pid,"_WGS_clean.summary.tree")
  if (file.exists(file)) {
    tree <- read.beast(file)
    tree2<-root.phylo(as.phylo(tree),tree@phylo$tip.label[tree@phylo$Nnode])
    if (pid=="C519"){
      tree2<-drop.tip(tree2,"A1_G5_D1")
    }else{
      tree2=tree2 
    }
    tree2$tip.label<-gsub("\\_D1$", "", tree2$tip.label)
    tree2<-root.phylo(tree2,outgroup = outgroupsample)
    tree2$tip.label<-gsub(paste0("\\EPICC_",pid,"_"), "", tree2$tip.label)
    tree3<-root(tree2,outgroupsample)
    tree4<-drop.tip(tree3, outgroupsample)
    sampledata<-data.frame(sample=tree2$tip.label,
                           #region=gsub("[0-9]_\\w+", "", tree2$tip.label),
                           region=sub("^(\\w).*$", "\\1", tree2$tip.label),
                           color=sample(c('black'), length(tree2$tip.label), replace=T))
    sampledata <-sampledata  %>% mutate(color = ifelse(sample %ni% pidPnormal & as.character(region) == "A", "#E41A1C",
                                                       ifelse(sample %ni% pidPnormal & as.character(region) == "B", "#377EB8",
                                                              ifelse(sample %ni% pidPnormal & as.character(region) == "C", "#4DAF4A", 
                                                                     ifelse(sample %ni% pidPnormal & as.character(region) == "D", "#984EA3",
                                                                            ifelse(sample %ni% pidPnormal & as.character(region) == "F", "black",
                                                                                   ifelse(sample %ni% pidPnormal & as.character(region) == "G", "black",
                                                                                          ifelse(sample %ni% pidPnormal & as.character(region) == "H", "black",
                                                                                                 ifelse(sample %ni% pidPnormal & as.character(region) == "Z", "grey",
                                                                                                        ifelse(sample %in% pidPnormal,"#00CDCD",
                                                                                                               "#FF7F00"))))))))))
    #tree4<-remove_root_tip(tree3,outgroupsample) 
    tree4<-tree3
    myplot<-ggtree(tree4,branch.length = "branch.length",size = 0.8)  %<+% sampledata+
      geom_tiplab(size=3,aes(color=I(color)),show.legend = NA) +
      geom_rootedge()+geom_tippoint(size=1,aes(color=I(color)),show.legend = NA) + 
      theme(plot.title = element_text(size=10,face="bold"),legend.position = c(0.8, 1),
            legend.justification = c(0,1),
            legend.title=element_text(size=10),legend.text=element_text(size=10)) + 
      ggtitle(pid)+xlim(0,2.5*max(tree@data$height)) + 
      ggtree::geom_treescale(y = -1,x = 0,fontsize = 3,linesize = 0.5,offset = 0.2) 
    #+ scale_color_manual(values=c("#E41A1C","#377EB8", "#4DAF4A", "#984EA3","grey","black","cyan3", "#FF7F00"), label=c("region A","region B","region C","region D","blood","polyps","adjacent normal","distant normal"))
    mylist[[pid]]<-myplot
    
  } else {
    mylist[[pid]]<-0
  }
}
############################################
#plot_legend 
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:4, ylim=0:4)
legend("topleft", ncol=1,legend =c("region A","region B","region C","region D","polyps","adjacent normal","distant normal"), 
       pch=16, pt.cex=1, cex=1, bty='n',
       col = c("#E41A1C","#377EB8", "#4DAF4A", "#984EA3","black","cyan3", "#FF7F00"))
mtext("Samples", at=0.1, cex=1.2)

######generate pdf files ###################
pdf(paste0("/Users/bchen/Desktop/Projects/PKSFiles/Figures/FigureS1.pdf"), width = 18, height =10)
multiplot(plotlist=c(mylist), cols=6)
dev.off()
############################################



