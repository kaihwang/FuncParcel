# script for plotting
setwd('/Volumes/neuro/bin/FuncParcel/Data/')

# import libraries
library(ggplot2)
library(plyr)
library(reshape2)
library(scales)

### global graph metrics
# read data
Data = read.csv('GraphGlobalData.csv', header=TRUE)
PatientData <- subset(Data, Group !='Control')
Variables_to_plot <- c('Q_zscore', 'CC_zscore')
for (v in Variables_to_plot){
  
  fig_global <- ggplot(data=PatientData, aes_string(x='Density', y=v)) + facet_wrap(~Group)
  fig_global <- fig_global + geom_line(aes(color=Subject), size = 2) + xlim( 0.05, 0.2)
  plot(fig_global)

}

### nodal properties
Data = read.csv('PatientsNodalZscoreData.csv', header=TRUE)
PatientData <- subset(Data, node_type !='all') 

PatientData$row.names <-NULL
PatientData$X <-NULL
rm(plotData)
plotData<-melt(data=PatientData,value.name = "z_score",variable.name = "Metric", id.vars=c("Group","Subject","Density", "node_type"), measure.vars = c("Target_PC", "nonTarget_PC"))

fig_nodal <- ggplot(data=plotData, aes(x=Density, y=z_score, linetype = Metric )) 
fig_nodal <- fig_nodal + facet_grid(node_type ~ Group)
fig_nodal <- fig_nodal + geom_line(aes(color=Subject), size = 1) + xlim( 0.05, 0.2) 
plot(fig_nodal)
