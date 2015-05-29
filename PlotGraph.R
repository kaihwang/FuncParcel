# script for plotting
setwd('/Volumes/neuro/bin/FuncParcel/Data/')

# import libraries
library(ggplot2)
library(plyr)
library(reshape2)
library(scales)

BlueScale<-c("#9ecae1", "#6baed6",  "#4292c6", "#2171b5", "#084594")
WarmScale<-c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
### global graph metrics
# read data
Data = read.csv('GraphGlobalData.csv', header=TRUE)
PatientData <- subset(Data, Group !='Control')
Variables_to_plot <- c('Q_zscore', 'CC_zscore')
for (v in Variables_to_plot){
  
  fig_global <- ggplot(data=PatientData, aes_string(x='Density', y=v)) + facet_wrap(~Group)
  fig_global <- fig_global + geom_line(aes(color=Subject), size = 2) + xlim( 0.05, 0.2) + scale_colour_manual(values=c(BlueScale, WarmScale))
  plot(fig_global)

}

### nodal properties
Data = read.csv('PatientsNodalZscoreData.csv', header=TRUE)
PatientData <- subset(Data, node_type !='all') 
PatientData$row.names <-NULL
PatientData$X <-NULL
Variables_to_plot <- c("PC", "Between_Module_Weight", "WMD", "Within_Module_Weight")

for (v in Variables_to_plot){

  plotData<-melt(data=Data,value.name = "z_score",variable.name = "Metric", id.vars=c("Group","Subject","Density", "node_type"), measure.vars = c(paste("Target_",v, sep=""), paste("nonTarget_", v, sep="")))
  levels(plotData$node_type) <- c("All Cortical ROIs", "Cortical Connector Hubs", "Non Hubs", "Cortical Provincial Hubs")
  levels(plotData$Group) <- c("Striatal Patients", "Thalamic Patients")
  fig_nodal <- ggplot(data=plotData, aes(x=Density, y=z_score, linetype = Metric )) + ggtitle(gsub("_"," ",v))
  fig_nodal <- fig_nodal + facet_grid(Group ~ node_type) + scale_colour_manual(name = "Patients", values=c(BlueScale, WarmScale))
  fig_nodal <- fig_nodal + geom_line(aes(color=Subject), size = 1) + xlim( 0.05, 0.2) +labs(y ="Z Score") + scale_linetype_discrete(name = "", labels = c("Target", "Non-Target"))

  plot(fig_nodal)
}


# plot NMI

NMIDATA <- read.csv('nmi.csv', header = TRUE)
NMIplot <- ggplot(data=NMIDATA, aes(x=Group, y=Value, color = Group)) + geom_boxplot(size=0.5)+geom_point(size=4) + scale_colour_manual(values=c("Black","Yellow", "Blue","#888888")) 
NMIplot <- NMIplot + theme_grey(base_size = 16) + labs(y="Mutual Information")
NMIplot
