setwd('/Volumes/neuro/bin/FuncParcel/Data/')


# import libraries
library(ggplot2)
library(reshape2)

#load data

Thalamus_Data = read.csv('Thalamus_nodal.csv', header=TRUE)
Thalamus_Data <- Thalamus_Data[Thalamus_Data$PC<1,]  
#Thalamus_Data <- Thalamus_Data[Thalamus_Data$target_PC<20,]

Cortical_Data = read.csv('Cortical_nodal.csv', header=TRUE)
Cortical_Data <- Cortical_Data[Cortical_Data$PC>0,]  
Cortical_Data <- Cortical_Data[Cortical_Data$Associated.System!='Other',] 

CI_colors <- c("#008080", "purple", "green", "red", "yellow", "black", "cyan", "pink", "blue")

### boxplot to compare nodal roles between each partition, for thalamus and cortex
Variables_to_plot <- c('PC', 'WMD', 'NNC', 'BNWR', 'bcc')
for (v in Variables_to_plot){
  Thalamus_boxplot <- ggplot(data = Thalamus_Data, aes_string(x="Associated.System", y=v, fill="Associated.System", colour="Associated.System")) + geom_boxplot(outlier.colour = NULL) #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Thalamus_boxplot <- Thalamus_boxplot + ggtitle(paste("Thalamus", v)) + scale_fill_manual(values=CI_colors)  + theme_grey(base_size = 32)  + theme(axis.text.x=element_text(angle=90, hjust=1))
  plot(Thalamus_boxplot)

  Cortical_boxplot <- ggplot(data = Cortical_Data, aes_string(x="Associated.System", y=v, fill="Associated.System", colour="Associated.System")) + geom_boxplot(outlier.colour = NULL) #+geom_point() + geom_jitter(position = position_jitter(width = .1))
  Cortical_boxplot <- Cortical_boxplot + ggtitle(paste("Cortical", v))+ scale_fill_manual(values=CI_colors) + theme_grey(base_size = 32)  + theme(axis.text.x=element_text(angle=90, hjust=1))
  plot(Cortical_boxplot)
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

Thalamus_pc_plot <-ggplot(data = total, aes(x=Associated.System.y , y=Target.PC, color =Associated.System.y ))+geom_point() + geom_jitter(position = position_jitter(width = .1)) +  scale_colour_manual(values=CI_colors) 
plot(Thalamus_pc_plot)


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

plot_data = melt(patient_data, id_vars = c("SubjID", "Voxel", "CI", "PC"), measure.vars = c("Target_connected_weight", "nonTarget_connected_weight"))

patient_plot <-ggplot(plot_data , aes(x=PC , y=value, colour=variable))+geom_point(size=3, aes(shape = factor(SubjID)))+ scale_colour_manual(values=c("Red", "Black")) + facet_grid(.~CI)
                                                                                                                                            
patient_plot
