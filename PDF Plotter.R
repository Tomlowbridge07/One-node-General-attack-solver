#This file has the pdf plotter

library(ggplot2)
library(reshape)
library(tikzDevice)

PlotPDF<-function(PDF,from,to,stepsize)
{
  XCoordinates=seq(from,to,stepsize)
  YCoordinates=vector(length=length(XCoordinates))
  for(i in 1:length(YCoordinates))
  {
    YCoordinates[i]=PDF(XCoordinates[i])
  }
  XYCoordiantes=data.frame(XCoordinates)
  XYCoordiantes=cbind(XYCoordiantes,YCoordinates)
  MeltedDataFrame=melt(XYCoordiantes,id="XCoordinates")
  
  Plot<-ggplot(MeltedDataFrame,aes(x=XCoordinates,y=value),show.legend='True')+geom_line()
  print(Plot)
}