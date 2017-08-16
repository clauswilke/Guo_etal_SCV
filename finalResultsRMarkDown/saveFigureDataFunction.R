# save figure data

saveFigure<-function(figureObj,name, dataSetName, folderName)
{
  temp<-print(figureObj)
  temp<-as.data.frame(temp$data)
  temp%>%dplyr::select(time=x, density=density, group, colour)->temp
  write.csv(x = temp, file = paste0("./",folderName,"/","SCV_",dataSetName,"_",name,".csv"))
}

saveFigure2<-function(figureObj,name, dataSetName, folderName)
{
  temp<-print(figureObj)
  temp<-as.data.frame(temp$data)
  temp%>%dplyr::select(x, y, group, colour)->temp
  write.csv(x = temp, file = paste0("./",folderName,"/","SCV_",dataSetName,"_",name,".csv"))
}