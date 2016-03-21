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

### plot summary graph
X_order <- c('First Order \nThalamic Nuclei','Higher Order \nThalamic Nuclei','Nonspecific \nThalamic Nuclei',"Cortical \nConnector Hubs", "Cortical \nProvincial Hubs", "Cortical \nNone Hubs")
Thalamus_plus_cortical_data$Classification <-factor(Thalamus_plus_cortical_data$Classification, levels=X_order)

#PC
volplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x=PC))
volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
volplot <- volplot + facet_grid(. ~ Classification) + coord_flip() + theme_grey(base_size = 10)
volplot <- volplot  + xlim(0,1) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
plot(volplot)
ggsave(filename = 'PC_classification.pdf', plot = volplot, units = c("in"),width=6.5, height=1.7,dpi=300) 

#WMD
volplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x=WMD))
volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
volplot <- volplot + facet_grid(. ~ Classification) + coord_flip() + theme_grey(base_size = 10) 
volplot <- volplot  + xlim(-3,4.5) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
plot(volplot)
ggsave(filename = 'WMD_classification.pdf', plot = volplot, units = c("in"),width=6.5, height=1.7,dpi=300) 

#CogFlex
volplot <- ggplot(data = Thalamus_plus_cortical_data, aes(x=cog))
volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
volplot <- volplot + facet_grid(. ~ Classification) + coord_flip() + theme_grey(base_size = 10) 
volplot <- volplot  + xlim(0,8) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
plot(volplot)
ggsave(filename = 'Cog_classification.pdf', plot = volplot, units = c("in"),width=6.5, height=1.7,dpi=300) 


Nuclei_order <-c('AN', 'LD', 'MD', 'CL', 'CeM', 'CM', 'Pf', 'Li', 'PuA', 'PuI', 'PuL','PuM','LP','Po','SG','MGN','LGN','VA','VL','VM','VPI','VPL','VPM' )
Thalamus_Data$Morel.Parcellations_f = factor(Thalamus_Data$Morel.Parcellations, levels=Nuclei_order)
#CI_colors <- c("#008080", "purple", "green", "red", "yellow", "magenta", "cyan", "pink", "blue", "pink")

### volcano plot
Variables_to_plot <- c('PC'  ) #'NNC', 'BNWR', 'bcc' 'WMD'

for (v in Variables_to_plot){
  
  volplot <- ggplot(data = Thalamus_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_grid(. ~ Functional.Network) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim(0,1) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Thalamus Functional Atlas")
  plot(volplot)
  ggsave(filename = paste(v,'_tha_fn_density.pdf', sep=''), plot = volplot, units = c("in"),width=3.4, height=1.5) 
  
  volplot <- ggplot(data = Cortical_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_grid(. ~ Functional.Network) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim( 0,1) + theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Cortical ROI")
  plot(volplot)
  ggsave(filename = paste(v,'_cortical_fn_density.pdf', sep=''), plot = volplot, units = c("in"),width=3.4, height=1.5) 
  
  volplot <- ggplot(data = Thalamus_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_grid(. ~ Anatomical.Parcellations) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim( 0,1) +  theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Thalamus Oxford-FSL Atlas")
  plot(volplot)
  ggsave(filename = paste(v,'_thalamus_an_density.pdf', sep=''), plot = volplot, units = c("in"),width=2.85, height=1.5) 
  
  Thalamus_Data <- Thalamus_Data[Thalamus_Data$Morel.Parcellations!='Unclassified',] 
  volplot <- ggplot(data = Thalamus_Data, aes_string(x=v))
  volplot <- volplot + stat_density(aes(ymax = ..density..,  ymin = -..density.. ),geom = "ribbon", position = "identity" )
  volplot <- volplot + facet_wrap(~Morel.Parcellations_f, ncol = 8) + coord_flip() + theme_grey(base_size = 10) 
  volplot <- volplot  + xlim( 0,1) +theme( axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.text = element_text(colour = "black")) 
  volplot <- volplot + ggtitle("Thalamus Morel Atlas") + theme(panel.margin.y = unit(1.25, "cm"))
  plot(volplot)
  ggsave(filename = paste(v,'_thalamus_morel_density.pdf', sep=''), plot = volplot, units = c("in"),width=3.4, height=4.25) 
  
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
CI_colors <- c("#6B006C", "red","yellow","cyan","blue", "brown", "#0E6E6C", "#008100", "pink")
Nuclei_Data <- read.csv('NucleiNetworkStrength.csv')
Variables_to_plot <- c('MD', 'PuM', 'VL', 'An'  )
for (v in Variables_to_plot){
  plot_Data <- Nuclei_Data[Nuclei_Data$Nuclei==v,] 
  n_plot <- ggplot(data = plot_Data, aes(x=factor(Network), y=Connectivity.Strength))
  n_plot <- n_plot + geom_bar(stat = "identity", aes(fill=Network)) + labs(y = "Z-Score") + theme_classic(base_size = 10)
  n_plot <- n_plot +scale_fill_manual(values=CI_colors ) 
  n_plot <- n_plot + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text = element_text(colour = "black"), legend.position="none")
  ggsave(filename = paste(v,'_plot.pdf', sep=''), plot = n_plot, units = c("in"),width=2, height=1.25) 
  plot(n_plot)
}


### plot to compare with cog flexibility
cog_plot <- ggplot(data = Thalamus_Data, aes(x=PC, y=cog))
cog_plot <- cog_plot + geom_point()
cog_plot <- cog_plot + geom_smooth(method=lm)
plot(cog_plot)



  
                      