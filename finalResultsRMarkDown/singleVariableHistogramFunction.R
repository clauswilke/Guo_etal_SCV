# Single variable histogram function.

singleVariableHistogramFunction<-function(DF, x_axis, x_limits=NULL, colVec, maximumTime, Adjust=1)
{
  if(x_axis=="COMB_maximum_y")
  {x_axis_label="Maximum"
    title_part1="Maximum Distribution"}
  
  if(x_axis=="COMB_slope1")
  {x_axis_label="Slope"
    title_part1="Slope Distribution"}
  
  if(x_axis=="COMB_midPoint1_x")
  {x_axis_label="Midpoint"
    title_part1="Midpoint Distribution"}
  
  if(x_axis=="COMB_incrementTime")
  {x_axis_label="Infection Time"
    title_part1="Infection Time Distribution"}
  
  if(x_axis=="COMB_startPoint_x")
  {x_axis_label="Start Point"
   title_part1="Start Point Distribution"}
  
  elements=as.vector(unique(DF$decision_bio))
  if(setequal(elements,c("infection&lysis"))){title_part2="(InfLys)"}
  if(setequal(elements,c("infection"))){title_part2="(Inf)"}
  if(setequal(elements,c("infection", "infection&lysis"))){title_part2="(Inf + InfLys)"}
  
  
  fig<-ggplot2::ggplot(DF, aes_string(x = x_axis, color = "dataSet_bio")) +
    geom_line(aes(y=..density..), stat="density",size=2, adjust=Adjust) +
    xlab(x_axis_label)+
    ggtitle(paste(title_part1,title_part2,sep = "  "))+
    theme_bw()+
    scale_colour_manual(values=colVec,
                        name=paste0("Conditions\n", maximumTime, " hours"))+
    theme( axis.text.x=element_text(size=16),
           axis.text.y=element_text(size=16),
           axis.title.x=element_text(size=16),
           axis.title.y=element_text(size=16),
           legend.title=element_text(size=14),
           legend.text=element_text(size=14),
           legend.position="bottom",
           legend.direction="vertical")+
    scale_y_continuous(expand = c(0,0)) 
  
  yRange<-ggplot_build(fig)$layout$panel_ranges[[1]]$y.range
  yRange2<-1.1*yRange[2]
  
  xRange<-ggplot_build(fig)$layout$panel_ranges[[1]]$x.range
  xRange2<-1.1*xRange[2]
  
  fig<-ggplot2::ggplot(DF, aes_string(x = x_axis, color = "dataSet_bio")) +
    geom_line(aes(y=..density..), stat="density",size=2, alpha=0.8, adjust=Adjust) +
    xlab(x_axis_label)+
    ggtitle(paste(title_part1,title_part2,sep = "  "))+
    theme_bw()+
    scale_colour_manual(values=colVec,
                        name=paste0("Conditions\n",maximumTime, " hours"))+
    theme( axis.text.x=element_text(size=16),
           axis.text.y=element_text(size=16),
           axis.title.x=element_text(size=16),
           axis.title.y=element_text(size=16),
           legend.title=element_text(size=14),
           legend.text=element_text(size=14),
           panel.grid.major = element_blank(),
           legend.position="bottom",
           legend.direction="vertical")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, yRange2))
  
      
  if(!is.null(x_limits))
  {fig=fig+scale_x_continuous(expand = c(0, 0), limits = x_limits)}
  
  if(is.null(x_limits))
  {fig=fig+scale_x_continuous(expand = c(0, 0))}
  
  return(fig)
}