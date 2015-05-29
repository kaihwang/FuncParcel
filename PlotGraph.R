# script for plotting
setwd('~/bin//FuncParcel/Data/')

# import libraries
library(ggplot2)
library(plyr)
library(reshape2)
library(scales)

### global graph metrics
# read data
Data = read.csv('~/bin/FuncParcel/Data/GraphGlobalData.csv', header=TRUE)
PatientData <- subset(Data, Group !='Control')
Variables_to_plot <- c('Q_zscore', 'CC_zscore')
for (v in Variables_to_plot){
  
  fig_global <- ggplot(data=PatientData, aes_string(x='Density', y=v)) + facet_wrap(~Group)
  fig_global <- fig_global + geom_line(aes(color=Subject), size = 2) + xlim( 0.05, 0.2)
  plot(fig_global)

}

### nodal properties
Data = read.csv('~/bin/FuncParcel/Data/PatientsNodalZscoreData.csv', header=TRUE)
PatientData <- subset(Data, node_type !='all') 
#plotData<-melt(data=NodalDATA, id.vars=c("Subject","Density"), 
               measure.vars = c("Cortical_Target_Within_Weight", "Cortical_nonTarget_Within_Weight"))
PatientData$row.names <-NULL
PatientData$X <-NULL
rm(plotData)
plotData<-melt(data=PatientData, id.vars=c("Group","Subject","Density", "node_type"), measure.vars = c("Target_Between_Module_Weight", "nonTarget_Between_Module_Weight"))

fig_nodal <- ggplot(data=plotData, aes(x=Density, y=value, linetype = variable )) 
fig_nodal <- fig_nodal + facet_grid(Group ~ node_type)
fig_nodal <- fig_nodal + geom_line(aes(color=Subject), size = 1) + xlim( 0.05, 0.2) 
plot(fig_nodal)
