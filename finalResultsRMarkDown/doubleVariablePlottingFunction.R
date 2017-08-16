# Single variable histogram function.

doubleVariablePlottingFunction<-function(DF, line_vs_point,
                                         x_axis, y_axis, 
                                         x_limits=NULL, y_limits=NULL,
                                         colVec, maximumTime,
                                         sizeForLines=2, sizeForPoints=2)
{
  if(x_axis=="COMB_maximum_y")
  {x_axis_label="Maximum"}
  
  if(x_axis=="COMB_slope1")
  {x_axis_label="Slope"}
  
  if(x_axis=="COMB_midPoint1_x")
  {x_axis_label="Midpoint"}
  
  if(x_axis=="COMB_incrementTime")
  {x_axis_label="Infection Time"}
  
  if(x_axis=="COMB_startPoint_x")
  {x_axis_label="Start Point"}
  
  if(y_axis=="COMB_maximum_y")
  {y_axis_label="Maximum"}
  
  if(y_axis=="COMB_slope1")
  {y_axis_label="Slope"}
  
  if(y_axis=="COMB_midPoint1_x")
  {y_axis_label="Midpoint"}
  
  if(y_axis=="COMB_incrementTime")
  {y_axis_label="Infection Time"}
  
  if(y_axis=="COMB_startPoint_x")
  {y_axis_label="Start Point"}
  
  elements=as.vector(unique(DF$decision_bio))
  if(setequal(elements,c("infection"))){title_part3="(Inf)"}
  if(setequal(elements,c("infection&lysis"))){title_part3="(InfLys)"}
  if(setequal(elements,c("infection", "infection&lysis"))){title_part3="(Inf + InfLys)"}
  
  
  # Start for figure
  fig<-ggplot(DF, aes_string(x = x_axis, y = y_axis, color = "dataSet_bio"))
  
  # Adding geoms to figure
  if(line_vs_point=="line"){fig<-fig+geom_density2d(size = sizeForLines, alpha=0.7)}
  if(line_vs_point=="point"){fig<-fig+geom_point(size = sizeForPoints, alpha=0.7) }
  
  # Figure style
  fig= fig +
    theme_bw()+
    scale_colour_manual(values=colVec,
                        name=paste0("Conditions\n",maximumTime, " hours"))+
    xlab(x_axis_label)+
    ylab(y_axis_label)+
    ggtitle(paste(x_axis_label, "vs", y_axis_label, 
                  "Distribution",title_part3, sep = " "))+
    theme( axis.text.x=element_text(size=16),
           axis.text.y=element_text(size=16),
           axis.title.x=element_text(size=16),
           axis.title.y=element_text(size=16),
           legend.title=element_text(size=14),
           legend.text=element_text(size=14) )
  
  # Figure axis limits
  if(!is.null(x_limits))
  {fig=fig+scale_x_continuous(expand = c(0, 0), limits = x_limits)}
  
  if(is.null(x_limits))
  {fig=fig+scale_x_continuous(expand = c(0, 0))}
  
  if(!is.null(y_limits))
  {fig=fig+scale_y_continuous(expand = c(0, 0), limits = y_limits)}
  
  if(is.null(x_limits))
  {fig=fig+scale_y_continuous(expand = c(0, 0))}
  
  
  return(fig)
}