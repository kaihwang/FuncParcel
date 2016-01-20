setwd('/Volumes/neuro/bin/FuncParcel/Data/')

# import libraries
library(ggplot2)
library(reshape2)
library(plyr)
#load data

Thalamus_Data = read.csv('Thalamus_nodal_WTA.csv', header=TRUE)
Cortical_Data = read.csv('Cortical_nodal_WTA.csv', header=TRUE)
Cortical_Data <- Cortical_Data[Cortical_Data$Associated.System!='Other',] 


CI_colors <- c("#008080", "purple", "green", "red", "yellow", "magenta", "cyan", "pink", "blue", "pink")

### boxplot to compare nodal roles between each partition, for thalamus and cortex
Variables_to_plot <- c('PC'  ) #'NNC', 'BNWR', 'bcc' 'WMD'


for (v in Variables_to_plot){
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Associated.System", y=v, fill="Associated.System", colour="Associated.System")) + geom_boxplot(outlier.colour = NULL) #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Thalamus", v)) + scale_fill_manual(values=CI_colors)  + theme_grey(base_size = 32)  + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position="none")
  Thalamus_boxplot <- Thalamus_boxplot + ylim( 0, 1) 
  plot(Thalamus_boxplot)


  #Cortical_boxplot <- ggplot(data = Cortical_Data, aes_string(x="Associated.System", y=v, fill="Associated.System", colour="Associated.System")) + geom_boxplot(outlier.colour = NULL) #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  #Cortical_boxplot <- Cortical_boxplot + ggtitle(paste("Cortical", v))+ scale_fill_manual(values=CI_colors) + theme_grey(base_size = 32)  + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position="none")
  #Cortical_boxplot <- Cortical_boxplot + ylim( -1, 5)
  #plot(Cortical_boxplot)
  #ggsave(filename = paste(v,'_cortical_box_consensus.pdf', sep=''), plot = Cortical_boxplot, units = c("in"),width=6, height=4) 
}


#plot patient
PT_Data = read.csv('Patient_Q_v_PC.csv', header=TRUE)
pt_q_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Hemispheric_Difference_in_Modularity))
pt_q_plot <- pt_q_plot + geom_bar(stat = "identity")
pt_q_plot <- pt_q_plot + labs(x= "Patients", y = "Hemispheric Difference in Modularity") + theme_grey(base_size = 20)
ggsave(filename = "pt_q_plot.pdf", plot = pt_q_plot, units = c("in"),width=6, height=6) 
plot(pt_q_plot)

pt_pc_plot <- ggplot(data = PT_Data, aes(x=factor(SubjID), y=Lesioned_site_PC))
pt_pc_plot <- pt_pc_plot + geom_bar(stat = "identity")
pt_pc_plot <- pt_pc_plot + labs(x= "Patients", y = "Lesioned voxel's PC") + theme_grey(base_size = 20)
ggsave(filename = "pt_pc_plot.pdf", plot = pt_pc_plot, units = c("in"),width=6, height=6) 
plot(pt_q_plot)


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


  
                      