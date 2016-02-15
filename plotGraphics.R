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

Nuclei_order <-c('AN', 'LD', 'MD', 'CL', 'CeM', 'CM', 'Pf', 'Li', 'PuA', 'PuI', 'PuL','PuM','LP','Po','SG','MGN','LGN','VA','VL','VM','VPI','VPL','VPM' )
Thalamus_Data$Morel.Parcellations_f = factor(Thalamus_Data$Morel.Parcellations, levels=Nuclei_order)

#CI_colors <- c("#008080", "purple", "green", "red", "yellow", "magenta", "cyan", "pink", "blue", "pink")

### volcano plot
Variables_to_plot <- c('WMD'  ) #'NNC', 'BNWR', 'bcc' 'WMD'

for (v in Variables_to_plot){
  
  volplot <- ggplot(data = Thalamus_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_grid(. ~ Functional.Network) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim(-2,2) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Thalamus Functional Atlas")
  plot(volplot)
  ggsave(filename = paste(v,'_tha_fn_density.pdf', sep=''), plot = volplot, units = c("in"),width=3.4, height=1.5) 
  
  volplot <- ggplot(data = Cortical_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_grid(. ~ Functional.Network) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim( 0,9) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Cortical ROI")
  plot(volplot)
  ggsave(filename = paste(v,'_cortical_fn_density.pdf', sep=''), plot = volplot, units = c("in"),width=3.4, height=1.5) 
  
  volplot <- ggplot(data = Thalamus_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_grid(. ~ Anatomical.Parcellations) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim( 0,9) +  theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Thalamus Oxford-FSL Atlas")
  plot(volplot)
  ggsave(filename = paste(v,'_thalamus_an_density.pdf', sep=''), plot = volplot, units = c("in"),width=2.85, height=1.5) 
  
  Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',] 
  volplot <- ggplot(data = Thalamus_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_wrap(~Morel.Parcellations_f, ncol = 8) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim( 0,9) +theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Thalamus Morel Atlas") + theme(panel.margin.y = unit(1.25, "cm"))
  plot(volplot)
  ggsave(filename = paste(v,'_thalamus_morel_density.pdf', sep=''), plot = volplot, units = c("in"),width=3.4, height=4.25) 
  
}  

### boxplot to compare nodal roles between each partition, for thalamus and cortex


for (v in Variables_to_plot){
  
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Functional Parcellation", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_tha_fn.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=2.5, height=2.5) 

  Cortical_boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, ,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Cortical_boxplot <- Cortical_boxplot + ggtitle(paste("Cortical ROI", v)) + theme_grey(base_size = 10) +  theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  Cortical_boxplot <- Cortical_boxplot + ylim( 0, 1)
  plot(Cortical_boxplot)
  ggsave(filename = paste(v,'_cortical_box_consensus.pdf', sep=''), plot = Cortical_boxplot, units = c("in"),width=3, height=2.5) 
  
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Oxford-FSL Atlas", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_fsl_fn.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=2, height=2.5) 
  
  Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',] 
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Morel Atlas", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_morel_fn.png', sep=''), plot = Thalamus_boxplot, units = c("in"),width=6, height=2.5) 
}

Variables_to_plot <- c('BNWR'  )
for (v in Variables_to_plot){
  
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Functional Parcellation", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_tha_fn.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=2.5, height=2.5) 
  
  Cortical_boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, ,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Cortical_boxplot <- Cortical_boxplot + ggtitle(paste("Cortical ROI", v)) + theme_grey(base_size = 10) +  theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  Cortical_boxplot <- Cortical_boxplot + ylim( 0, 1)
  plot(Cortical_boxplot)
  ggsave(filename = paste(v,'_cortical_box_consensus.pdf', sep=''), plot = Cortical_boxplot, units = c("in"),width=3, height=2.5) 
  
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Oxford-FSL Atlas", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_fsl_fn.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=2, height=2.5) 
  
  Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',] 
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Morel Atlas", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_morel_fn.png', sep=''), plot = Thalamus_boxplot, units = c("in"),width=6, height=2.5) 
}

Variables_to_plot <- c('WMD'  )
for (v in Variables_to_plot){
  
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Functional Parcellation", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( -2, 2) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_tha_fn.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=2.5, height=2.5) 
  
  Cortical_boxplot <- ggplot(data = Cortical_Data, aes_string(x="Functional.Network", y=v)) + geom_boxplot(outlier.colour = NULL, ,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Cortical_boxplot <- Cortical_boxplot + ggtitle(paste("Cortical ROI", v)) + theme_grey(base_size = 10) +  theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  Cortical_boxplot <- Cortical_boxplot + ylim( -2, 2)
  plot(Cortical_boxplot)
  ggsave(filename = paste(v,'_cortical_box_consensus.pdf', sep=''), plot = Cortical_boxplot, units = c("in"),width=3, height=2.5) 
  
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Anatomical.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Oxford-FSL Atlas", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank(), axis.text.y=element_blank(),legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( -2, 2) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_fsl_fn.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=2, height=2.5) 
  
  Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',] 
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Morel.Parcellations", y=v)) + geom_boxplot(outlier.colour = NULL,fill = "grey80") #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Morel Atlas", v))  + theme_grey(base_size = 10)  + theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
  Thalamus_boxplot <- Thalamus_boxplot + ylim( -2, 2) 
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_morel_fn.png', sep=''), plot = Thalamus_boxplot, units = c("in"),width=7, height=2.5) 
}


#### patient stuff
setwd('~/Google Drive/Projects/Thalamus-Rest/')
#plot patient
PT_Data = read.csv('Patient_Q_v_PC.csv', header=TRUE)
pt_q_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Q.Diff))
pt_q_plot <- pt_q_plot + geom_bar(stat = "identity")
pt_q_plot <- pt_q_plot + labs(x= "Patients", y = "Hemispheric Difference in Modularity") + theme_grey(base_size = 10)+ theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black"))
ggsave(filename = "pt_q_plot.pdf", plot = pt_q_plot, units = c("in"),width=3, height=3) 
plot(pt_q_plot)

pt_pc_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Lesioned.PC))
pt_pc_plot <- pt_pc_plot + geom_bar(stat = "identity") 
pt_pc_plot <- pt_pc_plot + ylim( 0, 80) 
pt_pc_plot <- pt_pc_plot + labs(x= "Patients", y = "Lesioned voxel's mean PC") + theme_grey(base_size = 10)+ theme( axis.title.x=element_blank(), legend.position="none", axis.text = element_text(colour = "black")) 
ggsave(filename = "pt_pc_plot.pdf", plot = pt_pc_plot, units = c("in"),width=3, height=3) 
plot(pt_pc_plot)


#### AN, MD, PuM, Intra network stregnth 
setwd('~/Google Drive/Projects/Thalamus-Rest/')
CI_colors <- c("#0E6E6C", "#6B006C", "red", "yellow", "cyan", "brown", "#008100", "pink", "blue")
Nuclei_Data <- read.csv('NucleiNetworkStrength.csv')
plot_Data <- Nuclei_Data[Nuclei_Data$Nuclei=='MD',] 
n_plot <- ggplot(data = plot_Data, aes(x=factor(Network), y=Connectivity.Strength))
n_plot <- n_plot + geom_bar(stat = "identity", aes(fill=Network)) + labs(y = "Z-Score") + theme_classic(base_size = 10)
n_plot <- n_plot +scale_fill_manual(values=CI_colors ) 
n_plot <- n_plot + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text = element_text(colour = "black"), legend.position="none")
ggsave(filename = "MD_plot.pdf", plot = n_plot, units = c("in"),width=2, height=1.25) 
plot(n_plot)


### plot to compare within thalamus nodal role
Variables_to_plot <- c('PC','within_PC','WMD', 'within_WMD', 'target_PC')

Thalamus_pc_plot <-ggplot(data = Thalamus_Data, aes(x=target_PC, y=PC))+geom_point()
plot(Thalamus_pc_plot)

### look at correlation between dataset
Thalamus_pc_plot <-ggplot(data = Thalamus_Data, aes(x=MGH_NNC, y=NKI_NNC))+geom_point()
plot(Thalamus_pc_plot)

### look at correlation between thalamic voxel nodal property and its cortical target
Thalamus_Target_Data = read.csv('Thalamus_target.csv', header=TRUE)
total <- merge(Thalamus_Target_Data,Thalamus_Data,by="Voxel")

Thalamus_pc_plot <-ggplot(data = total, aes(x=PC , y=Target.PC))+geom_point() + geom_jitter(position = position_jitter(width = .1)) +  scale_colour_manual(values=CI_colors) 
Thalamus_pc_plot <- Thalamus_pc_plot + stat_smooth(method = 'lm', fill="blue", colour="darkblue", size=2) + theme_grey(base_size = 32) + ylab("Cortical ROIs' PC")+ xlab("Thalamic Voxels' PC")    
plot(Thalamus_pc_plot)
ggsave(filename = "Thalamocortical_PC.pdf", plot = Thalamus_pc_plot, units = c("in"),width=15, height=8) 

Thalamus_wmd_plot <-ggplot(data = total, aes(x=WMD , y=Target.WMD))+geom_point() + geom_jitter(position = position_jitter(width = .1)) +  scale_colour_manual(values=CI_colors) 
Thalamus_wmd_plot <- Thalamus_wmd_plot + stat_smooth(method = 'lm', fill="blue", colour="darkblue", size=2) + theme_grey(base_size = 32)+ ylab("Cortical ROIs' WMD")+ xlab("Thalamic Voxels' WMD") 
plot(Thalamus_wmd_plot)
ggsave(filename = "Thalamocortical_WMD.pdf", plot = Thalamus_wmd_plot, units = c("in"),width=15, height=8) 

## look at tsnr
Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Associated.System", y="TSNR", fill="Associated.System")) + geom_boxplot() #+geom_point() + geom_jitter(position = position_jitter(width = .1))
plot(Thalamus_boxplot)


  
                      