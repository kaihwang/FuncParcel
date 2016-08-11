setwd('/Volumes/neuro/bin/FuncParcel/Data/')

# import libraries
library(ggplot2)
library(reshape2)
library(plyr)
library("grid")
#load data

Thalamus_Data = read.csv('Thalamus_nodal_WTA.csv', header=TRUE)
Cortical_Data = read.csv('Cortical_nodal_WTA.csv', header=TRUE)
Cortical_Data <- Cortical_Data[Cortical_Data$Functional.Network!='Other',] 

Thalamus_plus_cortical_data <- read.csv('Cortical_plus_thalamus_nodal_WTA.csv', header=TRUE)
Thalamus_plus_cortical_data <- Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification!='Unclassified',] 
Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',] 
Thalamus_Data <- Thalamus_Data[Thalamus_Data$Functional.Network!='Other',] 

### plot summary graph
X_order <- c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nConnector Hubs", "Cortical \nProvincial Hubs", "Cortical \nNon Hubs")
Thalamus_plus_cortical_data$Classification <-factor(Thalamus_plus_cortical_data$Classification, levels=X_order)

#PC, try boxplot for Mark
Thalamus_plus_cortical_data_1<- Thalamus_plus_cortical_data
levels(Thalamus_plus_cortical_data_1$Classification) <-c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nConnector Hubs", "Cortical \nNon Connector Hubs", "Cortical \nNon Connector Hubs")


Thalamus_boxplot <- ggplot(data = Thalamus_plus_cortical_data_1, aes(x=Classification, y=PC)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
Thalamus_boxplot <- Thalamus_boxplot + theme_grey(base_size = 8)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
plot(Thalamus_boxplot)
ggsave(filename = 'PC_classification_box.pdf', plot = Thalamus_boxplot, units = c("in"),width=4, height=2,dpi=300) 

#PC kernal density
CI_colors <- c("green", "blue","red")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = PC)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.title = element_blank())
plot(kplot)
ggsave(filename = 'PC_classification_kernal.pdf', plot = kplot, units = c("in"),width=4, height=2,dpi=300) 

#bPC kernal density
CI_colors <- c("green", "blue","red")
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = bPC)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.title = element_blank())
plot(kplot)
ggsave(filename = 'bPC_classification_kernal.pdf', plot = kplot, units = c("in"),width=4, height=2,dpi=300) 



#WMD
Thalamus_plus_cortical_data_2<- Thalamus_plus_cortical_data
levels(Thalamus_plus_cortical_data_2$Classification) <-c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nNon Provincial Hubs", "Cortical \nProvincial Hubs", "Cortical \nNon Provincial Hubs")
X_order <- c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei',"Cortical \nProvincial Hubs", "Cortical \nNon Provincial Hubs")
Thalamus_plus_cortical_data_2$Classification <-factor(Thalamus_plus_cortical_data_2$Classification, levels=X_order)
Thalamus_boxplot <- ggplot(data = Thalamus_plus_cortical_data_2, aes(x=Classification, y=WMD)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
Thalamus_boxplot <- Thalamus_boxplot + theme_grey(base_size = 8)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
Thalamus_boxplot <- Thalamus_boxplot + ylim( -4, 4) 
plot(Thalamus_boxplot)
ggsave(filename = 'WMD_classification_box.pdf', plot = Thalamus_boxplot, units = c("in"),width=4, height=2,dpi=300) 

#WMD kernal density
kplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x = WMD)) 
kplot <-kplot + stat_density(aes( group = nodetype, fill = nodetype ),size=2, alpha=0.7, position="identity") 
kplot <-kplot + theme_grey(base_size = 8) +scale_fill_manual(values=CI_colors ) + theme(legend.title = element_blank())
plot(kplot)
ggsave(filename = 'WMD_classification_kernal.pdf', plot = kplot, units = c("in"),width=4, height=2,dpi=300) 

#CogFlex
Thalamus_boxplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x=Classification, y=cog)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
Thalamus_boxplot <- Thalamus_boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 7.5) +ylab("Cognitive Flexibility")
plot(Thalamus_boxplot)
ggsave(filename = 'Cog_classification_box.pdf', plot = Thalamus_boxplot, units = c("in"),width=5.5, height=2,dpi=300) 


Nuclei_order <-c('AN', 'LD', 'MD', 'CL', 'CeM', 'CM', 'Pf', 'Li', 'PuA', 'PuI', 'PuL','PuM','LP','Po','SG','MGN','LGN','VA','VL','VM','VPI','VPL','VPM' )
Thalamus_Data$Morel.Parcellations_f = factor(Thalamus_Data$Morel.Parcellations, levels=Nuclei_order)
#CI_colors <- c("#008080", "purple", "green", "red", "yellow", "magenta", "cyan", "pink", "blue", "pink")

### volcano plot
Variables_to_plot <- c('WMD'  ) #'NNC', 'BNWR', 'bcc' 'WMD'

for (v in Variables_to_plot){
  
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, outlier.shape = NA, fill = "grey80", coef=3) 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Thalamus Functional Atlas") + geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_tha_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.5) 
  
  boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Cortical ROIs") + geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_cortical_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.5) 
  
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Oxford-FSL Atlas") + geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed")  
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_an_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.85, height=1.5) 
  
  FO_Data <-Thalamus_Data[Thalamus_Data$Classification=='First Order \nThalamic Nuclei',] 
  boxplot <- ggplot(data = FO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Morel Atlas (First Order Nuclei)") + geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_fo_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.7, height=1.5) 
  
  HO_Data <-Thalamus_Data[Thalamus_Data$Classification=='Higher Order \nThalamic Nuclei',] 
  boxplot <- ggplot(data = HO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Morel Atlas (Higher Order Nuclei)") + geom_hline(aes(yintercept=.8), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_ho_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.5) 
  
}  

Variables_to_plot <- c('PC'  ) #'NNC', 'BNWR', 'bcc' 'WMD'

for (v in Variables_to_plot){
  
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, outlier.shape = NA, fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Thalamus Functional Atlas") + geom_hline(aes(yintercept=.63), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_tha_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.5) 
  
  boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Cortical ROIs") + geom_hline(aes(yintercept=.63), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_cortical_fn_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.5) 
  
  boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Oxford-FSL Atlas") + geom_hline(aes(yintercept=.63), colour="#990000", linetype="dashed")  
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_an_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.85, height=1.5) 
  
  FO_Data <-Thalamus_Data[Thalamus_Data$Classification=='First Order \nThalamic Nuclei',] 
  boxplot <- ggplot(data = FO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Morel Atlas (First Order Nuclei)") + geom_hline(aes(yintercept=.63), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_fo_box.pdf', sep=''), plot = boxplot, units = c("in"),width=2.7, height=1.5) 
  
  HO_Data <-Thalamus_Data[Thalamus_Data$Classification=='Higher Order \nThalamic Nuclei',] 
  boxplot <- ggplot(data = HO_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,outlier.shape = NA,fill = "grey80") 
  boxplot <- boxplot + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  boxplot <- boxplot + ylim( -2, 4) +ggtitle("Morel Atlas (Higher Order Nuclei)") + geom_hline(aes(yintercept=.63), colour="#990000", linetype="dashed")
  plot(boxplot)
  ggsave(filename = paste(v,'_thalamus_ho_box.pdf', sep=''), plot = boxplot, units = c("in"),width=3.4, height=1.5) 
  
}  


#### patient stuff
setwd('~/Google Drive/Projects/Thalamus-Rest/')
# plot patient lesion extent
Variables_to_plot <- c('S1', 'S2', 'S3', 'S4')

for (v in Variables_to_plot){
  PT_Data <-read.csv('Patient_locaition.csv', header=TRUE)
  pData <- PT_Data[PT_Data$Patient==v,]
  pData <- pData[pData$Atlas == 'Functional',]
  pt_plot <-ggplot(data=pData, aes(x=Location, y=Size))
  pt_plot <- pt_plot  + geom_bar(stat="identity") + theme_classic(base_size = 10) + theme( axis.title.x=element_blank()) + ylab(bquote('Lesion Size ('~mm^3*')'))
  plot(pt_plot)
  ggsave(filename = paste(v,'_pt_loc_functional.pdf', sep=''), plot = pt_plot, units = c("in"),width=2, height=1.5) 
  
  pData <- PT_Data[PT_Data$Patient==v,]
  pData <- pData[pData$Atlas == 'Morel',]
  pt_plot <-ggplot(data=pData, aes(x=Location, y=Size))
  pt_plot <- pt_plot  + geom_bar(stat="identity") + theme_classic(base_size = 10) + theme( axis.title.x=element_blank())+ ylab(bquote('Lesion Size ('~mm^3*')'))
  plot(pt_plot)
  ggsave(filename = paste(v,"_pt_loc_morel.pdf", sep=''), plot = pt_plot,units = c("in"),width=2,  height=1.5)
}

#plot patient
PT_Data = read.csv('Patient_Q_v_PC.csv', header=TRUE)
pt_q_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Q.Diff))
pt_q_plot <- pt_q_plot + geom_bar(stat = "identity")
pt_q_plot <- pt_q_plot + labs(x= "Patients", y = "Hemispheric Difference \n in Modularity") + theme_classic(base_size = 10)+ theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
ggsave(filename = "pt_q_plot.pdf", plot = pt_q_plot, units = c("in"),width=1.5, height=2) 
plot(pt_q_plot)

pt_pc_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Lesioned.PC))
pt_pc_plot <- pt_pc_plot + geom_bar(stat = "identity") 
pt_pc_plot <- pt_pc_plot + ylim( 0, 1) 
pt_pc_plot <- pt_pc_plot + labs(x= "Patients", y = "Mean PC (lesion)") + theme_classic(base_size = 10)+ theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black")) 
ggsave(filename = "pt_pc_plot.pdf", plot = pt_pc_plot, units = c("in"),width=1.5, height=2) 
plot(pt_pc_plot)


pt_pc_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Lesioned.WMD))
pt_pc_plot <- pt_pc_plot + geom_bar(stat = "identity") 
pt_pc_plot <- pt_pc_plot + ylim( 0, 2) 
pt_pc_plot <- pt_pc_plot + labs(x= "Patients", y = "Mean WMD (lesion)") + theme_classic(base_size = 10)+ theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black")) 
ggsave(filename = "pt_WMD_plot.pdf", plot = pt_pc_plot, units = c("in"),width=1.5, height=2) 
plot(pt_pc_plot)

#### AN, MD, PuM, Intra network stregnth 
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
  ggsave(filename = paste(v,'_plot.pdf', sep=''), plot = n_plot, units = c("in"),width=1.5, height=1.25) 
  plot(n_plot)
}


### plot to compare with cog flexibility
cog_plot <- ggplot(data = Thalamus_Data, aes(x=PC, y=cog))
cog_plot <- cog_plot + geom_point()
cog_plot <- cog_plot + geom_smooth(method=lm)
plot(cog_plot)


## do permutation statistics
DV<-c(Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='First Order \nThalamic Nuclei',14], Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='Cortical \nConnector Hubs',14])

IV <- factor(rep(c("FO", "CC"), c(length(Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='First Order \nThalamic Nuclei',14]), length(Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='Cortical \nConnector Hubs',14]))))  
library(coin)  
pvalue(oneway_test(DV ~ IV, distribution=approximate(B=9999)))
mean(Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='First Order \nThalamic Nuclei',14])
mean(Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='Higher Order \nThalamic Nuclei',14])
mean(Thalamus_plus_cortical_data[Thalamus_plus_cortical_data$Classification=='Cortical \nConnector Hubs',14])
