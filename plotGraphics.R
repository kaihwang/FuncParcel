setwd('/Volumes/neuro/bin/FuncParcel/Data/')


# import libraries
library(ggplot2)
library(reshape2)
library(plyr)
#load data

Thalamus_Data = read.csv('Thalamus_nodal.csv', header=TRUE)
Thalamus_Data <- Thalamus_Data[Thalamus_Data$PC<1,]  
#Thalamus_Data <- Thalamus_Data[Thalamus_Data$target_PC<20,]

Cortical_Data = read.csv('Cortical_nodal.csv', header=TRUE)
Cortical_Data <- Cortical_Data[Cortical_Data$PC>0,]  
Cortical_Data <- Cortical_Data[Cortical_Data$Associated.System!='Other',] 

CI_colors <- c("#008080", "purple", "green", "red", "yellow", "magenta", "cyan", "pink", "blue")

### boxplot to compare nodal roles between each partition, for thalamus and cortex
Variables_to_plot <- c('WMD'  ) #'NNC', 'BNWR', 'bcc' 'WMD'
for (v in Variables_to_plot){
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Associated.System", y=v, fill="Associated.System", colour="Associated.System")) + geom_boxplot(outlier.colour = NULL) #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Thalamus", v)) + scale_fill_manual(values=CI_colors)  + theme_grey(base_size = 32)  + theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.x=element_blank(), legend.position="none")
  Thalamus_boxplot <- Thalamus_boxplot + ylim( -1, 5)
  plot(Thalamus_boxplot)
  ggsave(filename = paste(v,'_tha_box.pdf', sep=''), plot = Thalamus_boxplot, units = c("in"),width=8, height=8) 

  Cortical_boxplot <- ggplot(data = Cortical_Data, aes_string(x="Associated.System", y=v, fill="Associated.System", colour="Associated.System")) + geom_boxplot(outlier.colour = NULL) #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Cortical_boxplot <- Cortical_boxplot + ggtitle(paste("Cortical", v))+ scale_fill_manual(values=CI_colors) + theme_grey(base_size = 32)  + theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.x=element_blank(), legend.position="none")
  Cortical_boxplot <- Cortical_boxplot + ylim( -1, 5)
  plot(Cortical_boxplot)
  ggsave(filename = paste(v,'_cortical_box.pdf', sep=''), plot = Cortical_boxplot, units = c("in"),width=8, height=8) 
}

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




## look at patient data
WarmScale<-c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
BlueScale<-c("#9ecae1", "#6baed6",  "#4292c6", "#2171b5", "#084594")
cScale <- c("#e7298a",
"#beaed4",
"#fdc086",
"#ffff99",
"#386cb0")
patient_data = read.csv('patient_df.csv', header=TRUE)
melt_data <- melt(patient_data, id_vars = c("SubjID", "Voxel", "CI", "PC"), measure.vars = c("Target_total_weight_bn", "nonTarget_total_weight_bn"))
plot_data <- ddply(melt_data, c("SubjID", "CI", "variable"), summarise, mean = mean(value, na.rm=TRUE))
patient_plot <-ggplot(plot_data , aes(x=SubjID , y=mean, colour=variable))+geom_point(size=3, aes(shape = factor(SubjID)))+ scale_colour_manual(values=c("Red", "Black")) + facet_grid(.~CI)
patient_plot <- patient_plot + theme_grey(base_size = 16) # + theme(axis.text.x=element_text(angle=90, hjust=1))                                                                                                                                            
patient_plot


patient_data = read.csv('patient_nework_df.csv', header=TRUE)
patient_data$SubjID <- factor(patient_data$SubjID)
patient_plot <-ggplot(patient_data, aes(x=SubjID, y=Between_network_connectivity_weight))+geom_bar(stat="identity")
patient_plot                      
                      
                      