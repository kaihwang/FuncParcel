setwd('/Volumes/neuro/bin/FuncParcel/Data/')

### import libraries
library(scales)
library(ggplot2)
library(reshape2)
library(plyr)
library("grid")


###load data
Thalamus_Data = read.csv('Thalamus_nodal_WTA.csv', header=TRUE)
Cortical_Data = read.csv('Cortical_nodal_WTA.csv', header=TRUE)
Cortical_Data <- Cortical_Data[Cortical_Data$Functional.Network!='Other',] #take out unassigned ROIs

Thalamus_plus_cortical_data <- read.csv('Cortical_plus_thalamus_nodal_WTA.csv', header=TRUE)
Thalamus_plus_cortical_data <- Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification!='Unclassified',]  #take out voxels not within Morel atlas
Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',]  #take out voxels not within Morel atlas
Thalamus_Data <- Thalamus_Data[Thalamus_Data$Functional.Network!='Other',] #take out unassigned ROIs

### plot summary graphs
X_order <- c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nConnector Hubs", "Cortical \nProvincial Hubs", "Cortical \nNon Hubs")
Thalamus_plus_cortical_data$Classification <-factor(Thalamus_plus_cortical_data$Classification, levels=X_order)

#PC, try boxplot for Mark
Thalamus_plus_cortical_data_1<- Thalamus_plus_cortical_data
levels(Thalamus_plus_cortical_data_1$Classification) <-c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nConnector Hubs", "Cortical \nNon Connector Hubs", "Cortical \nNon Connector Hubs")

Thalamus_boxplot <- ggplot(data = Thalamus_plus_cortical_data_1, aes(x=Classification, y=PC)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
Thalamus_boxplot <- Thalamus_boxplot + theme_grey(base_size = 8)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
plot(Thalamus_boxplot)
#ggsave(filename = 'PC_classification_box.pdf', plot = Thalamus_boxplot, units = c("in"),width=4, height=2,dpi=300) 

#PC kernal density
CI_colors <- c("blue", "red","yellow")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = PC)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none") #+ theme(legend.title = element_blank())
kplot <-kplot + coord_cartesian(ylim=c(0,7),xlim=c(0,1)) +ylab("Kernel Density") + xlab("PC (partial correlations)")
plot(kplot)
ggsave(filename = 'PC_classification_kernal.pdf', plot = kplot, units = c("in"),width=2.3, height=2,dpi=300) 

#PC kernal contrast full and partial
CI_colors <- c("blue", "red","yellow")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = fPC)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none") #+ theme(legend.title = element_blank())
kplot <-kplot + coord_cartesian(ylim=c(0,7),xlim=c(0,1)) + ylab("Kernel Density") + xlab("PC (full correlations)")
plot(kplot)
ggsave(filename = 'fPC_classification_kernal.pdf', plot = kplot, units = c("in"),width=2.3, height=2,dpi=300) 

#bPC kernal density
#CI_colors <- c("green", "blue","red")
#kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = bPC)) 
#kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
#kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none")#+ theme(legend.title = element_blank())
#plot(kplot)
#ggsave(filename = 'bPC_classification_kernal.pdf', plot = kplot, units = c("in"),width=3, height=2,dpi=300) 

#PC kernal density for each func network, plus "stick" for thalamic parcel
X_order <- c('CO','DM','T','latO','mO','mT','sFP','SM','FP')
Cortical_Data$Functional.Network <-factor(Cortical_Data$Functional.Network, levels=X_order)
CI_colors <- c("#6B006C", "red","pink","cyan","blue", "brown", "#0E6E6C", "#008100", "yellow")
kplot <- ggplot(data = Cortical_Data, aes(x = PC)) 
kplot <-kplot + stat_density(aes( group = Functional.Network,fill = Functional.Network ), alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_color_manual(values=CI_colors )+scale_fill_manual(values=CI_colors ) + theme(legend.position="none")
kplot <-kplot + geom_vline(xintercept=0.60, color = "#6B006C", size=1, linetype="F1") 
kplot <-kplot + geom_vline(xintercept=0.62, color = "red", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.94, color = "yellow", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.78, color = "cyan", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.79, color = "blue", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.76, color = "brown", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.91, color = "#0E6E6C", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.68, color = "#008100", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.83, color = "pink", size=1, linetype="F1")
kplot <- kplot + coord_cartesian(ylim=c(0,5))
plot(kplot)
#ggsave(filename = 'PC_pacel_kernal.pdf', plot = kplot, units = c("in"),width=2.8, height=2,dpi=300)

#WMD boxplot
Thalamus_plus_cortical_data_2<- Thalamus_plus_cortical_data
levels(Thalamus_plus_cortical_data_2$Classification) <-c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nNon Provincial Hubs", "Cortical \nProvincial Hubs", "Cortical \nNon Provincial Hubs")
X_order <- c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nProvincial Hubs", "Cortical \nNon Provincial Hubs")
Thalamus_plus_cortical_data_2$Classification <-factor(Thalamus_plus_cortical_data_2$Classification, levels=X_order)
Thalamus_boxplot <- ggplot(data = Thalamus_plus_cortical_data_2, aes(x=Classification, y=WMD)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
Thalamus_boxplot <- Thalamus_boxplot + theme_grey(base_size = 8)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
Thalamus_boxplot <- Thalamus_boxplot + ylim( -4, 4) 
plot(Thalamus_boxplot)
#ggsave(filename = 'WMD_classification_box.pdf', plot = Thalamus_boxplot, units = c("in"),width=4, height=2,dpi=300) 

#WMD kernal density
CI_colors <-c("blue", "red","yellow")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = WMD)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none")
kplot <-kplot + coord_cartesian(ylim=c(0,.6), xlim=c( -4.2, 8.2))   + ylab("Kernel Density") + xlab("WMD (partial correlations)")
plot(kplot)
ggsave(filename = 'WMD_classification_kernal.pdf', plot = kplot, units = c("in"),width=2.3, height=2,dpi=300) 

CI_colors <-c("blue", "red","yellow")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = fWMD)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none")
kplot <-kplot + coord_cartesian(ylim=c(0,.6), xlim=c( -4.2, 8.2))   + ylab("Kernel Density") + xlab("WMD (full correlations)")
plot(kplot)
ggsave(filename = 'fWMD_classification_kernal.pdf', plot = kplot, units = c("in"),width=2.3, height=2,dpi=300) 

#WMD kernal density for each func network, plus "stick" for thalamic parcel
CI_colors <- c("#6B006C", "red","pink","cyan","blue", "brown", "#0E6E6C", "#008100", "yellow")
kplot <- ggplot(data = Cortical_Data, aes(x = WMD)) 
kplot <-kplot + stat_density(aes( group = Functional.Network,fill = Functional.Network ), alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none")
kplot <-kplot + xlim(-4,4)# +ylim(0,1)
kplot <-kplot + geom_vline(xintercept=3.1, color = "#6B006C", size=1, linetype="F1")  #
kplot <-kplot + geom_vline(xintercept=2.93, color = "red", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=2.95, color = "yellow", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=-0.36, color = "cyan", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.88, color = "blue", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=1.23, color = "brown", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=4.04, color = "#0E6E6C", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=0.04, color = "#008100", size=1, linetype="F1")
kplot <-kplot + geom_vline(xintercept=-.35, color = "pink", size=1, linetype="F1")
kplot <- kplot + coord_cartesian(ylim=c(0,.75))
plot(kplot)
ggsave(filename = 'WMD_parcel_kernal.pdf', plot = kplot, units = c("in"),width=2.8, height=2,dpi=300) 

#CogFlex boxplot
Thalamus_boxplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x=Classification, y=cog)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
Thalamus_boxplot <- Thalamus_boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 7.5) +ylab("Cognitive Flexibility")
plot(Thalamus_boxplot)
ggsave(filename = 'Cog_classification_box.pdf', plot = Thalamus_boxplot, units = c("in"),width=5.5, height=2,dpi=300) 

#CogFlex kernal plot
CI_colors <-c("blue", "red","yellow")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = cog)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity")  + ylab("Kernel Density")
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme( legend.position="none",axis.text = element_text(colour = "black")) + xlab("Cognitive Flexibility")
plot(kplot)
ggsave(filename = 'Cog_kernal.pdf', plot = kplot, units = c("in"),width=3, height=2,dpi=300) 


##2D wmd-pc
newrow <- data.frame(PC=c(.60, .62, .94, .785, .786, .76, .91, .68, .84), WMD=c(3.1, 2.94, 2.94, -.36, .82, 1.23, 4.04, .04, -.35), nodetype="Thalamic Parcels", Functional.Network=c("CO","DM", "FP", "latO", "mO", "mT", "sFP", "SM", "T"))
Data = rbind.fill(Cortical_Data, newrow)
CI_colors <- c("#6B006C", "red","yellow","cyan","blue", "brown", "#0E6E6C", "#008100", "pink")
wmdpcplt<-ggplot(data = Data, aes(x=PC, y=WMD, color=Functional.Network, shape = nodetype, size=nodetype)) + geom_point(size=2.5 ) + scale_color_manual(values=CI_colors ) + facet_wrap(~ Functional.Network) +xlim(0,1) + ylim(-2.25, 4.1)
wmdpcplt<- wmdpcplt + theme_grey(base_size = 8) 
plot(wmdpcplt)
ggsave(filename = '2d_pc_by_wmd.pdf', plot = wmdpcplt, units = c("in"),width=6, height=5,dpi=300) 

###WMD kernal plot facet
CI_colors <- c("#6B006C", "red","yellow","cyan","blue", "brown", "#0E6E6C", "#008100", "pink")
kplot <- ggplot(data = Cortical_Data, aes(x = WMD)) 
kplot <-kplot + stat_density(aes( group = Functional.Network,fill = Functional.Network ), alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none")
kplot <-kplot + xlim(-4,4.1) + facet_wrap(~ Functional.Network)
parcelmean<-c(3.1, 2.94, 2.94, -.36, .82, 1.23, 4.04, .04, -.35)
Functional.Network<-c("CO","DM", "FP", "latO", "mO", "mT", "sFP", "SM", "T")
parceldf<- data.frame(parcelmean,Functional.Network)
kplot <-kplot+ geom_vline(data =parceldf, aes(xintercept=parcelmean), size=1)  + coord_cartesian(ylim=c(0,1)) + facet_wrap(~ Functional.Network)  + ylab("Kernel Density")
plot(kplot)
ggsave(filename = 'kernal_wmd_by_network.pdf', plot = kplot, units = c("in"),width=3, height=3,dpi=300) 

###PC kernal plot facet
CI_colors <- c("#6B006C", "red","yellow","cyan","blue", "brown", "#0E6E6C", "#008100", "pink")
kplot <- ggplot(data = Cortical_Data, aes(x = PC)) 
kplot <-kplot + stat_density(aes( group = Functional.Network,fill = Functional.Network ), alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.position="none")
kplot <-kplot + xlim(0,1) + facet_wrap(~ Functional.Network) + scale_x_continuous(breaks=c(0,0.4,0.8))
parcelmean<-c(.60, .62, .94, .785, .786, .76, .91, .68, .84)
Functional.Network<-c("CO","DM", "FP", "latO", "mO", "mT", "sFP", "SM", "T")
parceldf<- data.frame(parcelmean,Functional.Network)
kplot <-kplot+ geom_vline(data =parceldf, aes(xintercept=parcelmean), size=1) + facet_wrap(~ Functional.Network) + coord_cartesian(ylim=c(0,5))  + ylab("Kernel Density")
plot(kplot)
ggsave(filename = 'kernal_pc_by_network.pdf', plot = kplot, units = c("in"),width=3, height=3,dpi=300) 


### box plot for atlases
Nuclei_order <-c('AN', 'LD', 'MD', 'CL', 'CeM', 'CM', 'Pf', 'Li', 'PuA', 'PuI', 'PuL','PuM','LP','Po','SG','MGN','LGN','VA','VL','VM','VPI','VPL','VPM' )
Thalamus_Data$Morel.Parcellations_f = factor(Thalamus_Data$Morel.Parcellations, levels=Nuclei_order) 
#CI_colors <- c("#008080", "purple", "green", "red", "yellow", "magenta", "cyan", "pink", "blue", "pink")

Variables_to_plot <- c('WMD'  ) #'NNC', 'BNWR', 'bcc' 'WMD'
for (v in Variables_to_plot){
  parcelmean<-c(3.1, 2.94, 2.94, -.36, .82, 1.23, 4.04, .04, -.35)
  category<-c("CO","DM", "FP", "latO", "mO", "mT", "sFP", "SM", "T")
  parceldf<- data.frame(parcelmean,category)
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, outlier.shape = NA, fill = "grey80", coef=3) 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 5) +ggtitle("Thalamus Functional Atlas") + geom_hline(aes(yintercept=1.04), colour="blue", linetype="dashed")
  boxplot <- boxplot + geom_boxplot(data=parceldf,aes(factor(category), y=parcelmean),inherit.aes=FALSE,color="gold",size=0.5)
  plot(boxplot)
  ggsave(filename = paste(v,'_tha_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.3) 
  
  boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 5) +ggtitle("Cortical ROIs") + geom_hline(aes(yintercept=1.04), colour="blue", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_cortical_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.3) 
  
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 5) +ggtitle("Oxford-FSL Atlas") + geom_hline(aes(yintercept=1.04), colour="blue", linetype="dashed")  
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_an_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.85, height=1.3) 
  
  FO_Data <-Thalamus_Data[Thalamus_Data$Classification=='First Order \nThalamic Nuclei',] 
  boxplot <- ggplot(data = FO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 5) +ggtitle("Morel Atlas (First Order Nuclei)") + geom_hline(aes(yintercept=1.04), colour="blue", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_fo_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.7, height=1.3) 
  
  HO_Data <-Thalamus_Data[Thalamus_Data$Classification=='Higher Order \nThalamic Nuclei',] 
  boxplot <- ggplot(data = HO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 5) +ggtitle("Morel Atlas (Higher Order Nuclei)") + geom_hline(aes(yintercept=1.04), colour="blue", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_ho_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.3) 
  
}  

Variables_to_plot <- c('PC'  ) #'NNC', 'BNWR', 'bcc' 'WMD'
for (v in Variables_to_plot){
  parcelmean<-c(.60, .62, .94, .785, .786, .76, .91, .68, .84)
  category<-c("CO","DM", "FP", "latO", "mO", "mT", "sFP", "SM", "T")
  parceldf<- data.frame(parcelmean,category)
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, outlier.shape = NA, fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( 0,1) +ggtitle("Thalamus Functional Atlas") + geom_hline(aes(yintercept=.63), colour="blue", linetype="dashed")
  boxplot <- boxplot + geom_boxplot(data=parceldf,aes(factor(category), y=parcelmean),inherit.aes=FALSE,color="gold",size=0.5)
  plot(boxplot)
  ggsave(filename = paste(v,'_tha_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.3) 
  
  boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( 0,1) +ggtitle("Cortical ROIs") + geom_hline(aes(yintercept=.63), colour="blue", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_cortical_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.3) 
  
  parcelmean<-c(.54, .84, .78, .88, .71, .54, .83 )
  category<-c("M","O", "PFC", "PL", "pM", "S", "T")
  parceldf<- data.frame(parcelmean,category)
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( 0,1) +ggtitle("Oxford-FSL Atlas") + geom_hline(aes(yintercept=.63), colour="blue", linetype="dashed")  
  boxplot <- boxplot + geom_boxplot(data=parceldf,aes(factor(category), y=parcelmean),inherit.aes=FALSE,color="gold",size=0.5)
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_an_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.85, height=1.3) 
  
  FO_Data <-Thalamus_Data[Thalamus_Data$Classification=='First Order \nThalamic Nuclei',] 
  parcelmean<-c(.74, .82, .62, .8, .69 )
  category<-c("AN","LGN", "MGN", "VL", "VP")
  parceldf<- data.frame(parcelmean,category)
  boxplot <- ggplot(data = FO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( 0,1) +ggtitle("Morel Atlas (First Order Nuclei)") + geom_hline(aes(yintercept=.63), colour="blue", linetype="dashed")
  boxplot <- boxplot + geom_boxplot(data=parceldf,aes(factor(category), y=parcelmean),inherit.aes=FALSE,color="gold",size=0.5)  
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_fo_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.7, height=1.3) 
  
  HO_Data <-Thalamus_Data[Thalamus_Data$Classification=='Higher Order \nThalamic Nuclei',] 
  parcelmean<-c(.79, .86, .84, .83, .63, .77, .67, .82, .57, .56 )
  category<-c("IL","LP", "MD", "Po", "PuA", "PuI", "PuL", "PuM", "VA", "VM")
  parceldf<- data.frame(parcelmean,category)
  boxplot <- ggplot(data = HO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( 0,1) +ggtitle("Morel Atlas (Higher Order Nuclei)") + geom_hline(aes(yintercept=.63), colour="blue", linetype="dashed")
  boxplot <- boxplot + geom_boxplot(data=parceldf,aes(factor(category), y=parcelmean),inherit.aes=FALSE,color="gold",size=0.5)    
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_ho_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.3) 
  
}  


### patient stuff
setwd('/Volumes/neuro/bin/FuncParcel/Data/')
Data = read.csv('Q_df.csv', header=TRUE)
mData <-melt(Data,id=c('SubjID','Density','PC.Damage.Score'))
levels(mData$variable) <-c('Whole \nBrain','Lesioned \nHemisphere',"Intact \nHemisphere")
qplot <- ggplot(data = mData, aes(x=Density, y=value, color =SubjID)) + geom_line(size=1.5) + ylim(-6.5,0) + facet_wrap(~variable)
qplot <- qplot + ylab('Q (z-score)') + scale_colour_discrete(name  ="Patient") + theme_grey(base_size = 10)                                                              
plot(qplot)
ggsave(filename ='Q.pdf', plot = qplot, units = c("in"),width=4.5, height=1.7) 

setwd('~/Google Drive/Projects/Thalamus-Rest/')
# plot patient lesion extent
Variables_to_plot <- c('S1', 'S2', 'S3', 'S4')
for (v in Variables_to_plot){
  PT_Data <-read.csv('Patient_locaition.csv', header=TRUE)
  
  #pData <- PT_Data[PT_Data$Patient==v,]
  #pData <- pData[pData$Atlas == 'Functional',]
  #pt_plot <-ggplot(data=pData, aes(x=Location, y=Size))
  #pt_plot <- pt_plot  + geom_bar(stat="identity") + theme_classic(base_size = 10) + theme( axis.title.x=element_blank()) + ylab(bquote('Lesion Size ('~mm^3*')'))
  #plot(pt_plot)
  #ggsave(filename = paste(v,'_pt_loc_functional.pdf', sep=''), plot = pt_plot, units = c("in"),width=2, height=1.5) 
  
  pData <- PT_Data[PT_Data$Patient==v,]
  pData <- pData[pData$Atlas == 'Morel',]
  pt_plot <-ggplot(data=pData, aes(x=Location, y=Size))
  pt_plot <- pt_plot  + geom_bar(stat="identity") + theme_classic(base_size = 10) + theme( axis.title.x=element_blank())+ ylab(bquote('Lesion Size ('~mm^3*')'))
  plot(pt_plot)
  ggsave(filename = paste(v,"_pt_loc_morel.pdf", sep=''), plot = pt_plot,units = c("in"),width=2,  height=1.5)
}

#plot patient v control NMI
NMIDATA <- read.csv('NMI.csv', header = TRUE)
NMIplot <- ggplot(data=NMIDATA, aes(x=Group, y=NMI, color = Group)) +geom_point(size=2) + scale_colour_manual(values=c("Black","DarkGreen", "Blue","#888888")) 
NMIplot <- NMIplot + theme_grey(base_size = 8) + labs(y="Adjusted NMI") 
NMIplot <- NMIplot +theme(plot.title=element_text( size=28)) + scale_x_discrete(breaks=NULL)
NMIplot
ggsave(filename = "pt_nmi_plot.pdf", plot = NMIplot, units = c("in"),width=2.5, height=2)



### AN, MD, PuM, Intra network stregnth 
setwd('~/Google Drive/Projects/Thalamus-Rest/')
CI_colors <- c("#6B006C", "red","yellow","cyan","blue", "brown", "#0E6E6C", "#008100", "pink")
Nuclei_Data <- read.csv('NucleiNetworkStrength.csv')
Variables_to_plot <- c('MD', 'PuM', 'VL', 'An', 'LGN', 'LP', 'IL', 'VA','VM'  )
for (v in Variables_to_plot){
  plot_Data <- Nuclei_Data[Nuclei_Data$Nuclei==v,] 
  n_plot <- ggplot(data = plot_Data, aes(x=factor(Network), y=Connectivity.Porportion))
  n_plot <- n_plot + geom_bar(stat = "identity", aes(fill=Network)) + labs(y = "% of Total \nConnectivity Weight") + theme_classic(base_size = 8)
  n_plot <- n_plot +scale_fill_manual(values=CI_colors ) + geom_hline(aes(yintercept=11), colour="#990000", linetype="dashed")
  n_plot <- n_plot + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text = element_text(colour = "black"), legend.position="none")
  ggsave(filename = paste(v,'_plot.pdf', sep=''), plot = n_plot, units = c("in"),width=1.25, height=1) 
  plot(n_plot)
}


### plot to compare with cog flexibility
cog_plot <- ggplot(data = Thalamus_Data, aes(x=PC, y=cog))
cog_plot <- cog_plot + geom_point()
cog_plot <- cog_plot + geom_smooth(method=lm)
plot(cog_plot)



### look at behav correlate with TRSE
setwd('/Volumes/neuro/bin/FuncParcel/Data/')
RT_Data = read.csv('bseries_trse_RT.csv', header=TRUE)
Accu_Data = read.csv('bseries_trse_Accu.csv', header=TRUE)
Data<-merge(RT_Data,Accu_Data,by ='sub')
plt<-ggplot(Data, aes(relev_scene_CO, y=relev_scene_acc)) + geom_point(shape=19, size=3) + stat_smooth(method=lm, fullrange=TRUE) + theme_grey(base_size=24) 
plot(plt)
summary(lm(relev_scene_acc ~ HF_sFP , Data))


### behav correlate of TDSigEI
setwd('/Volumes/neuro/bin/FuncParcel/Data/')
Data = read.csv('TDSigEI_behav.csv', header=TRUE)
plt<-ggplot(Data, aes(FH_sFP, y=FH_Accu)) + geom_point(shape=19, size=3) + stat_smooth(method=lm, fullrange=TRUE) + theme_grey(base_size=24) 
plot(plt)
summary(lm(HF_RT~HF_VL , Data))


### HCP behav
setwd('/Volumes/neuro/bin/FuncParcel/Data/')
Data = read.csv('HCP_behav.csv', header=TRUE)
plt<-ggplot(Data, aes(WM_MD, y=WM)) + geom_point(shape=19, size=3) + stat_smooth(method=lm, fullrange=TRUE) + theme_grey(base_size=24) 
plot(plt)
summary(lm(RELATIONAL ~ (RELATIONAL_MD) , Data))

