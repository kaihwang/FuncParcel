# script for plotting
setwd('/Volumes/neuro/bin/FuncParcel/Data/')

# import libraries
library(ggplot2)
library(plyr)
library(reshape2)
library(scales)
library(grid)

BlueScale<-c("#9ecae1", "#6baed6",  "#4292c6", "#2171b5", "#084594")
WarmScale<-c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
### global graph metrics
# read data
Data = read.csv('GraphGlobalData.csv', header=TRUE)
PatientData <- subset(Data, Group !='Control')
levels(PatientData$Group) <- c("", "Striatal Patients", "Thalamic Patients")

Variables_to_plot <- c('Q_zscore')
for (v in Variables_to_plot){
  
  fig_global <- ggplot(data=PatientData, aes_string(x='Density', y=v)) + facet_wrap(~Group) + ggtitle('Modularity(Q)') + theme_grey(base_size = 32) + theme(panel.margin = unit(5, "lines"))
  fig_global <- fig_global + geom_line(aes(color=Subject), size = 2) + xlim( 0.05, 0.2) + scale_colour_manual(values=c(BlueScale, WarmScale), labels = c('T1', 'T2', 'T3', 'T4', 'T5', 'S1', 'S2', 'S3', 'S4', 'S5')) + labs(y ="Z Score")
  plot(fig_global)

}
ggsave(filename = "IndivModularity.pdf", plot = fig_global, width=10, height=8) 

### nodal properties
Data <- read.csv('PatientsNodalZscoreData.csv', header=TRUE)
PatientData <- subset(Data, node_type !='all') 
PatientData$row.names <-NULL
PatientData$X <-NULL
Variables_to_plot <- c("Between_Module_Weight", "Within_Module_Weight")

for (v in Variables_to_plot){

  plotData<-melt(data=Data,value.name = "z_score",variable.name = "Metric", id.vars=c("Group","Subject","Density", "node_type"), measure.vars = c(paste("Target_",v, sep=""), paste("nonTarget_", v, sep="")))
  levels(plotData$node_type) <- c("All Cortical ROIs", "Cortical Connector Hubs", "Non Hubs", "Cortical Provincial Hubs")
  levels(plotData$Group) <- c("Striatal Patients", "Thalamic Patients")
  fig_nodal <- ggplot(data=plotData, aes(x=Density, y=z_score, linetype = Metric )) + ggtitle(gsub("_"," ",v)) + theme_grey(base_size = 32)
  fig_nodal <- fig_nodal + facet_grid(Group ~ node_type) + scale_colour_manual(name = "Patients", values=c(BlueScale, WarmScale),labels = c('T1', 'T2', 'T3', 'T4', 'T5', 'S1', 'S2', 'S3', 'S4', 'S5') ) + theme(panel.margin = unit(0.5, "lines"))
  fig_nodal <- fig_nodal + geom_line(aes(color=Subject), size = 2) + xlim( 0.05, 0.2) + ylim( -2.5, 2.5) + labs(y ="Z Score") + scale_linetype_discrete(name = "", labels = c("Target", "Non-Target"))

  plot(fig_nodal)
  ggsave(filename = paste(v,'.pdf', sep=''), plot = fig_nodal, units = c("in"),width=24, height=8) 
}


# plot NMI

NMIDATA <- read.csv('NMI_dataframe.csv', header = TRUE)
levels(NMIDATA$Group) <- c("Controls","Striatal Patients", "Thalamic Patients")
NMIplot <- ggplot(data=NMIDATA, aes(x=Group, y=NMI, color = Group)) + geom_boxplot(size=0.5)+geom_point(size=4) + scale_colour_manual(values=c("Black","DarkGreen", "Blue","#888888")) 
NMIplot <- NMIplot + theme_grey(base_size = 32) + labs(y="NMI") + ggtitle('Normalized mutual information (NMI)\n of modularity partition \n compared to Power/Gordon Atlases')
NMIplot <- NMIplot +theme(plot.title=element_text( size=28)) + scale_x_discrete(breaks=NULL)
NMIplot

ggsave(filename = "NMI.pdf", plot = NMIplot, width=10, height=8) 
