# ---- Generate data sets that will be used for plots ----
#********************************************
# Load sicegar related files
for (counter01 in 1:length(listFileNames) ) {
  dataFileName=listFileNames[counter01];
  
  load(paste0('../tidyData/tidy',dataFileName,'.RDA'))
  assign(x = "tidyDf", value = get(paste0(dataFileName,"TidyGrouped")))
  tidyDf %>%
    dplyr::group_by(Cell,Mutation,MOI,dataSet) %>%
    dplyr::select(time = hpi, intensity = nRat_a)->tidyDf
  
  currentDf=read.csv(file = paste0("../detailedSicegarResults/sicegar_output_currentDf_",
                                   dataFileName,
                                   ".csv"))
  
  currentDf %>%
    dplyr::group_by(Cell,Mutation,MOI,dataSet) -> currentDf
  
  
  load(file=paste0("../tidyData/","tidy",dataFileName,"_DDOM",".RDA"))
  assign(x = "DDOM", value = get(paste0(dataFileName,"TidyGrouped_DDOM")))
  DDOM %>%
    dplyr::mutate(dataSet=dataFileName) %>%
    dplyr::mutate(DDOM=ifelse(Dirt==1|Dried_wells==1|Out_of_focus==1|Multiple_infections==1,1,0))->DDOM
  
  dplyr::left_join(currentDf,DDOM)->currentDf
  
  currentDf$Dirt=as.numeric(currentDf$Dirt)
  currentDf$Dried_wells=as.numeric(currentDf$Dried_wells)
  currentDf$Out_of_focus=as.numeric(currentDf$Out_of_focus)
  currentDf$Multiple_infections=as.numeric(currentDf$Multiple_infections)
  
  if(counter01==1)
  {
    currentDfList = currentDf
    tidyDfList = tidyDf
  }
  if(counter01!=1)
  {
    currentDfList = dplyr::bind_rows(currentDfList,currentDf)
    tidyDfList = dplyr::bind_rows(tidyDfList,tidyDf)
  }
  
  remove(currentDf, DDOM, tidyDf)
}
#******************************************


#******************************************
# Calculate Maximum time
maximumTime=max(tidyDfList$time)
#******************************************


#******************************************
# make a summary data frame that will be used as main frame
currentDfList -> currentDf
tidyDfList -> tidyDf
remove(currentDfList, tidyDfList)
#******************************************


#******************************************
# Generate biological decision // refactor currentDf
currentDf %>% dplyr::mutate(decision_bio = DEC_decision) -> currentDf
currentDf$decision_bio <- replace_fun(input_vector = currentDf$decision_bio,
                                      initialVal = c("sigmoidal", "double_sigmoidal", "no_signal"),
                                      finalVal = c("infection", "infection&lysis", "no signal"))
currentDf$decision_bio <- factor(currentDf$decision_bio,
                                 levels=c("no signal", "ambiguous", "infection", "infection&lysis"))

currentDf %>% dplyr::mutate(dataSet_bio = dataSet) -> currentDf
currentDf$dataSet_bio <- replace_fun(input_vector = currentDf$dataSet_bio,
                                     initialVal= listFileNames,
                                     finalVal= listConditionNames)
currentDf$dataSet_bio <- factor(currentDf$dataSet_bio,
                                levels = listOrderNames)

currentDf$Cell=factor(currentDf$Cell)
currentDf$MOI=factor(currentDf$MOI)
currentDf$dataSet=factor(currentDf$dataSet)
#******************************************


#******************************************
# Filter out DDOM==1
currentDf -> currentDfwithDDOM
currentDf %>% dplyr::filter(DDOM==0) -> currentDf
#******************************************


#******************************************
# Get rid of unnecessary columns
currentDf %>%
  dplyr::mutate(decision = DEC_decision) -> currentDf

currentDf[grep("SCALING_*.*", colnames(currentDf))] <- NULL
currentDf[grep("SM_*.*", colnames(currentDf))] <- NULL
currentDf[grep("DSM_*.*", colnames(currentDf))] <- NULL
currentDf[grep("DEC_*.*", colnames(currentDf))] <- NULL
#******************************************


#******************************************
# generate 3 data sets
currentDf %>%
  dplyr::filter(decision_bio == "infection" | decision_bio == "infection&lysis") %>%
  dplyr::group_by(Cell, Mutation, MOI, dataSet, dataSet_bio) -> temp_inf_infLys

currentDf %>%
  dplyr::filter(decision_bio == "infection&lysis") %>%
  dplyr::group_by(Cell, Mutation, MOI, dataSet, dataSet_bio) -> temp_infLys

currentDf %>%
  dplyr::filter(decision_bio == "infection") %>%
  dplyr::group_by(Cell, Mutation, MOI, dataSet, dataSet_bio) -> temp_inf
#******************************************


#******************************************
# save 4 data sets

dataSetName=paste(gsub("SCV","",listFileNames),collapse = "_")

folderName=paste0("variables_",dataSetName);
if(!dir.exists(folderName)){dir.create(folderName)}

fileName_csv=paste0("./",folderName,"/","inf_infLys_",dataSetName,".csv")
write.csv(x = temp_inf_infLys,file = fileName_csv)

fileName_csv=paste0("./",folderName,"/","infLys_",dataSetName,".csv")
write.csv(x = temp_infLys, file = fileName_csv)

fileName_csv=paste0("./",folderName,"/","inf_",dataSetName,".csv")
write.csv(x = temp_inf, file = fileName_csv)

fileName_csv=paste0("./",folderName,"/","detailed_",dataSetName,".csv")
write.csv(x = currentDfwithDDOM, file = fileName_csv)
# End of data set preperation
#******************************************


# ---- Distribution of labels ----
fig00A<-ggplot(currentDfwithDDOM, aes_string("dataSet_bio", fill="decision_bio")) +
  geom_bar(position="fill", width = .8) +
  xlab("conditions")+
  ylab("Relative Frequency")+
  ggtitle("Unfiltered")+
  theme_bw()+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  scale_fill_manual(values=c("red", "orange", "lightgreen", "darkgreen"),
                    name=paste0("conditions\n",maximumTime," hours"))
print(fig00A)


fig00B<-ggplot(currentDf, aes_string("dataSet_bio", fill="decision_bio")) +
  geom_bar(position="fill", width = .8) +
  xlab("conditions")+
  ylab("Relative Frequency")+
  ggtitle("Filtered")+
  theme_bw()+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  scale_fill_manual(values=c("red", "orange", "lightgreen", "darkgreen"),
                    name=paste0("conditions\n",maximumTime," hours"))
print(fig00B)

# ---- Distribution of labels WIDE ----
fig_distributionOfLabels <- cowplot::plot_grid(fig00A, fig00B, labels = c("A", "B"))
print(fig_distributionOfLabels)



# ---- Distribution of DDOM data ----
##******************************************
currentDfwithDDOM %>%
  dplyr::mutate(cellProblem=ifelse(Dirt==1,"Dirt","OK"),
                cellProblem=ifelse(Dried_wells==1, "Dried wells", cellProblem),
                cellProblem=ifelse(Out_of_focus==1, "Out of focus", cellProblem),
                cellProblem=ifelse(Multiple_infections==1, "Multiple infections", cellProblem)) %>%
  dplyr::group_by(cellProblem,decision_bio)->currentDfwithDDOM


currentDfwithDDOM %>%
  dplyr::summarise(numSamples=length(Cell)) %>%
  dplyr::group_by(decision_bio)%>%
  dplyr::mutate(percentage=100*numSamples/sum(numSamples))->currentDfwithDDOM_summary


template = as.data.frame(expand.grid(cellProblem=c("Dirt","Dried wells", "Out of focus", "Multiple infections", "OK"),
                                     decision_bio=c("no signal","ambiguous","infection","infection&lysis")))


dplyr::left_join(template, currentDfwithDDOM_summary) -> currentDfwithDDOM_summary
currentDfwithDDOM_summary$numSamples[is.na(currentDfwithDDOM_summary$numSamples)] <- 0
currentDfwithDDOM_summary$percentage[is.na(currentDfwithDDOM_summary$percentage)] <- 0

currentDfwithDDOM_summary %>%
  dplyr::group_by(cellProblem) %>%
  dplyr::summarise(numSamples=sum(numSamples))%>%
  dplyr::mutate(percentage=100*numSamples/sum(numSamples),
                decision_bio="Total")->currentDfwithDDOM_summary2

dplyr::bind_rows(currentDfwithDDOM_summary,currentDfwithDDOM_summary2)->currentDfwithDDOM_summary

currentDfwithDDOM_summary$cellProblem<-factor(currentDfwithDDOM_summary$cellProblem,
                                              levels=c("Dirt", "Dried wells", "Out of focus", "Multiple infections", "OK"))
##******************************************

# ---- Distribution of DDOM ----
##******************************************
fig00C<-ggplot(currentDfwithDDOM_summary, aes_string(x="decision_bio", y="percentage", group="cellProblem", fill="cellProblem")) +
  #geom_point(size=3)+
  geom_bar(stat="identity", position = "dodge", width = .8) +
  scale_x_discrete()+
  scale_fill_manual(values=c("#e41a1c", "#377eb8", "#984ea3", "#ff7f00", "green"),
                    name=paste0("Problem\n",maximumTime," hours"))+
  theme_bw()+
  xlab("data clusters")+
  ylab("Problematic Percentage")+
  ggtitle("Problematic Well Distribution")+
  theme( axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=12),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=12),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14),
         panel.grid.major.y= element_line(colour = "grey30"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,110))

print(fig00C)
##******************************************


# ---- Maximum Associated figures and tables ----

#***********************************************
# maxima figure limits
fig1a<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                       x_axis="COMB_maximum_y",
                                       colVec=colVec01, maximumTime=maximumTime)

fig2a<-singleVariableHistogramFunction(DF=temp_infLys,
                                       x_axis="COMB_maximum_y",
                                       colVec=colVec01, maximumTime=maximumTime)

fig3a<-singleVariableHistogramFunction(DF=temp_inf,
                                       x_axis="COMB_maximum_y",
                                       colVec=colVec01, maximumTime=maximumTime)

xRange<-ggplot_build(fig1a)$layout$panel_ranges[[1]]$x.range
xRange_1a<-1.0*xRange[2]
xRange<-ggplot_build(fig2a)$layout$panel_ranges[[1]]$x.range
xRange_2a<-1.0*xRange[2]
xRange<-ggplot_build(fig3a)$layout$panel_ranges[[1]]$x.range
xRange_3a<-1.0*xRange[2]

xRange=max(xRange_1a, xRange_2a, xRange_3a)
#***********************************************


#***********************************************
fig1<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                      x_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange),
                                      colVec=colVec01, maximumTime = maximumTime)
#print(fig1)
saveFigure(figureObj = fig1,name = "maximum_Estimate_inf_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj01_table_ttest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_maximum_y","t.test")
obj01_table_utest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_maximum_y","wilcox.test")
obj01_table_difMean <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_maximum_y","difMean")
obj01_table_effectSize <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_maximum_y","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj01_table_ttest[[1]],"**\n\n",obj01_table_ttest[[2]],"\n"))
knitr::kable(obj01_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj01_table_utest[[1]],"**\n\n",obj01_table_utest[[2]],"\n"))
knitr::kable(obj01_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj01_table_difMean[[1]],"**\n\n",obj01_table_difMean[[2]],"\n"))
knitr::kable(obj01_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj01_table_effectSize[[1]],"**\n\n",obj01_table_effectSize[[2]],"\n"))
knitr::kable(obj01_table_effectSize[[3]])
#***********************************************


#*****************************************
fig2<-singleVariableHistogramFunction(DF=temp_infLys,
                                      x_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange),
                                      colVec=colVec01, maximumTime=maximumTime)
#print(fig2)
saveFigure(figureObj = fig2, name = "maximum_Estimate_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj02_table_ttest <- tableFunction(temp_infLys,"dataSet_bio","COMB_maximum_y","t.test")
obj02_table_utest <- tableFunction(temp_infLys,"dataSet_bio","COMB_maximum_y","wilcox.test")
obj02_table_difMean <- tableFunction(temp_infLys,"dataSet_bio","COMB_maximum_y","difMean")
obj02_table_effectSize <- tableFunction(temp_infLys,"dataSet_bio","COMB_maximum_y","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj02_table_ttest[[1]],"**\n\n",obj02_table_ttest[[2]],"\n"))
knitr::kable(obj02_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj02_table_utest[[1]],"**\n\n",obj02_table_utest[[2]],"\n"))
knitr::kable(obj02_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj02_table_difMean[[1]],"**\n\n",obj02_table_difMean[[2]],"\n"))
knitr::kable(obj02_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj02_table_effectSize[[1]],"**\n\n",obj02_table_effectSize[[2]],"\n"))
knitr::kable(obj02_table_effectSize[[3]])
#*****************************************


#*****************************************
fig3<-singleVariableHistogramFunction(DF=temp_inf,
                                      x_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange),
                                      colVec=colVec01, maximumTime=maximumTime)
#print(fig1)
saveFigure(figureObj = fig3,name = "maximum_Estimate_inf", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj03_table_ttest <- tableFunction(temp_inf,"dataSet_bio","COMB_maximum_y","t.test")
obj03_table_utest <- tableFunction(temp_inf,"dataSet_bio","COMB_maximum_y","wilcox.test")
obj03_table_difMean <- tableFunction(temp_inf,"dataSet_bio","COMB_maximum_y","difMean")
obj03_table_effectSize <- tableFunction(temp_inf,"dataSet_bio","COMB_maximum_y","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj03_table_ttest[[1]],"**\n\n",obj03_table_ttest[[2]],"\n"))
knitr::kable(obj03_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj03_table_utest[[1]],"**\n\n",obj03_table_utest[[2]],"\n"))
knitr::kable(obj03_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj03_table_difMean[[1]],"**\n\n",obj03_table_difMean[[2]],"\n"))
knitr::kable(obj03_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj03_table_effectSize[[1]],"**\n\n",obj03_table_effectSize[[2]],"\n"))
knitr::kable(obj03_table_effectSize[[3]])
#***********************************************


# ---- Maximum Associated figures and tables WIDE ----
fig_MaxValue<-cowplot::plot_grid(fig1,fig2, fig3, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_MaxValue)


# ---- Maximum Associated violin ----
fig4a<-ggplot(temp_inf_infLys, aes(x=dataSet_bio, y=COMB_maximum_y)) +
  geom_violin(fill="grey50", scale = "area") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Maximum")+
  ggtitle("Maximum Distribution Violin (Inf+InfLys)")+
  #scale_y_log10()+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))+
  ylim(0,xRange)

#print(fig4a)

fig4b<-ggplot(temp_infLys, aes(x=dataSet_bio, y=COMB_maximum_y)) +
  geom_violin(fill="grey50", scale = "area") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Maximum")+
  ggtitle("Maximum Distribution Violin (InfLys)")+
  #scale_y_log10()+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))+
  ylim(0,xRange)

#print(fig4b)

fig4c<-ggplot(temp_inf, aes(x=dataSet_bio, y=COMB_maximum_y)) +
  geom_violin(fill="grey50", scale = "area") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Maximum")+
  ggtitle("Maximum Distribution Violin (Inf)")+
  #scale_y_log10()+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))+
  ylim(0,xRange)

#print(fig4c)

fig_MaxValue_Violin<-cowplot::plot_grid(fig4a, fig4b, fig4c, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_MaxValue_Violin)


# ---- Slope Associated figures and tables ----


#***********************************************
# slope figure limits
fig5a<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                       x_axis="COMB_slope1",
                                       colVec=colVec01, maximumTime=maximumTime)

fig6a<-singleVariableHistogramFunction(DF=temp_infLys,
                                       x_axis="COMB_slope1",
                                       colVec=colVec01, maximumTime=maximumTime)

fig7a<-singleVariableHistogramFunction(DF=temp_inf,
                                       x_axis="COMB_slope1",
                                       colVec=colVec01, maximumTime=maximumTime)

xRange<-ggplot_build(fig5a)$layout$panel_ranges[[1]]$x.range
xRange_5a<-1.0*xRange[2]
xRange<-ggplot_build(fig6a)$layout$panel_ranges[[1]]$x.range
xRange_6a<-1.0*xRange[2]
xRange<-ggplot_build(fig7a)$layout$panel_ranges[[1]]$x.range
xRange_7a<-1.0*xRange[2]

xRange=max(xRange_5a, xRange_6a, xRange_7a)
#***********************************************


#*****************************************
fig5<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                      x_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange),
                                      colVec=colVec01, maximumTime=maximumTime)
#print(fig5)
saveFigure(figureObj = fig5,name = "slope_Estimate_inf_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj05_table_ttest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_slope1","t.test")
obj05_table_utest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_slope1","wilcox.test")
obj05_table_difMean <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_slope1","difMean")
obj05_table_effectSize <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_slope1","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj05_table_ttest[[1]],"**\n\n",obj05_table_ttest[[2]],"\n"))
knitr::kable(obj05_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj05_table_utest[[1]],"**\n\n",obj05_table_utest[[2]],"\n"))
knitr::kable(obj05_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj05_table_difMean[[1]],"**\n\n",obj05_table_difMean[[2]],"\n"))
knitr::kable(obj05_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj05_table_effectSize[[1]],"**\n\n",obj05_table_effectSize[[2]],"\n"))
knitr::kable(obj05_table_effectSize[[3]])
#*****************************************


#*****************************************
fig6<-singleVariableHistogramFunction(DF=temp_infLys,
                                      x_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange),
                                      colVec=colVec01, maximumTime=maximumTime)

#print(fig6)
saveFigure(figureObj = fig6,name = "slope_Estimate_infLys", dataSetName = dataSetName, folderName=folderName)


# TABLE
obj06_table_ttest <- tableFunction(temp_infLys,"dataSet_bio","COMB_slope1","t.test")
obj06_table_utest <- tableFunction(temp_infLys,"dataSet_bio","COMB_slope1","wilcox.test")
obj06_table_difMean <- tableFunction(temp_infLys,"dataSet_bio","COMB_slope1","difMean")
obj06_table_effectSize <- tableFunction(temp_infLys,"dataSet_bio","COMB_slope1","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj06_table_ttest[[1]],"**\n\n",obj06_table_ttest[[2]],"\n"))
knitr::kable(obj06_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj06_table_utest[[1]],"**\n\n",obj06_table_utest[[2]],"\n"))
knitr::kable(obj06_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj06_table_difMean[[1]],"**\n\n",obj06_table_difMean[[2]],"\n"))
knitr::kable(obj06_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj06_table_effectSize[[1]],"**\n\n",obj06_table_effectSize[[2]],"\n"))
knitr::kable(obj06_table_effectSize[[3]])
#*****************************************


#*****************************************
fig7<-singleVariableHistogramFunction(DF=temp_inf,
                                      x_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange),
                                      colVec=colVec01, maximumTime=maximumTime)

#print(fig7)
saveFigure(figureObj = fig7,name = "slope_Estimate_inf", dataSetName = dataSetName, folderName=folderName)


# TABLE
obj07_table_ttest <- tableFunction(temp_inf,"dataSet_bio","COMB_slope1","t.test")
obj07_table_utest <- tableFunction(temp_inf,"dataSet_bio","COMB_slope1","wilcox.test")
obj07_table_difMean <- tableFunction(temp_inf,"dataSet_bio","COMB_slope1","difMean")
obj07_table_effectSize <- tableFunction(temp_inf,"dataSet_bio","COMB_slope1","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj07_table_ttest[[1]],"**\n\n",obj07_table_ttest[[2]],"\n"))
knitr::kable(obj07_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj07_table_utest[[1]],"**\n\n",obj07_table_utest[[2]],"\n"))
knitr::kable(obj07_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj07_table_difMean[[1]],"**\n\n",obj07_table_difMean[[2]],"\n"))
knitr::kable(obj07_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj07_table_effectSize[[1]],"**\n\n",obj07_table_effectSize[[2]],"\n"))
knitr::kable(obj07_table_effectSize[[3]])
#*****************************************



# ---- Slope Associated figures and tables WIDE ----
fig_Slope<-cowplot::plot_grid(fig5,fig6, fig7, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_Slope)


# ---- Slope Associated violin ----
fig8a<-ggplot(temp_inf_infLys, aes(x=dataSet_bio, y=COMB_slope1)) +
  geom_violin(fill="grey50", scale = "area") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Slope")+
  ggtitle("Slope Distribution Violin (Inf+InfLys)")+
  #scale_y_log10()+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))+
  ylim(0,xRange)

#print(fig8a)

fig8b<-ggplot(temp_infLys, aes(x=dataSet_bio, y=COMB_slope1)) +
  geom_violin(fill="grey50", scale = "area") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Slope")+
  ggtitle("Slope Distribution Violin (InfLys)")+
  #scale_y_log10()+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)

#print(fig8b)

fig8c<-ggplot(temp_inf, aes(x=dataSet_bio, y=COMB_slope1)) +
  geom_violin(fill="grey50", scale = "area") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Slope")+
  ggtitle("Slope Distribution Violin (Inf)")+
  #scale_y_log10()+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)

#print(fig8c)

fig_Slope_Violin<-cowplot::plot_grid(fig8a,fig8b,fig8c, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_Slope_Violin)


# ---- Midpoint Associated figures and tables ----

#*****************************************
# Scaling Related
fig11a<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                        x_axis="COMB_midPoint1_x",
                                        colVec=colVec01, maximumTime=maximumTime)

fig12a<-singleVariableHistogramFunction(DF=temp_infLys,
                                        x_axis="COMB_midPoint1_x",
                                        colVec=colVec01, maximumTime=maximumTime)

fig13a<-singleVariableHistogramFunction(DF=temp_inf,
                                        x_axis="COMB_midPoint1_x",
                                        colVec=colVec01, maximumTime=maximumTime)

xRange<-ggplot_build(fig11a)$layout$panel_ranges[[1]]$x.range
xRange_11a<-1.0*xRange[2]
xRange<-ggplot_build(fig12a)$layout$panel_ranges[[1]]$x.range
xRange_12a<-1.0*xRange[2]
xRange<-ggplot_build(fig13a)$layout$panel_ranges[[1]]$x.range
xRange_13a<-1.0*xRange[2]

xRange=max(xRange_11a, xRange_12a, xRange_13a)
#*****************************************


#*****************************************
fig11<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                       x_axis="COMB_midPoint1_x",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime)
#print(fig11)
saveFigure(figureObj = fig11,name = "midPoint_Estimate_inf_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj11_table_ttest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_midPoint1_x","t.test")
obj11_table_utest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_midPoint1_x","wilcox.test")
obj11_table_difMean <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_midPoint1_x","difMean")
obj11_table_effectSize <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_midPoint1_x","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj11_table_ttest[[1]],"**\n\n",obj11_table_ttest[[2]],"\n"))
knitr::kable(obj11_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj11_table_utest[[1]],"**\n\n",obj11_table_utest[[2]],"\n"))
knitr::kable(obj11_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj11_table_difMean[[1]],"**\n\n",obj11_table_difMean[[2]],"\n"))
knitr::kable(obj11_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj11_table_effectSize[[1]],"**\n\n",obj11_table_effectSize[[2]],"\n"))
knitr::kable(obj11_table_effectSize[[3]])
#*****************************************


#*****************************************
fig12<-singleVariableHistogramFunction(DF=temp_infLys,
                                       x_axis="COMB_midPoint1_x",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime)

#print(fig12)
saveFigure(figureObj = fig12,name = "midPoint_Estimate_infLys", dataSetName = dataSetName, folderName=folderName)


# TABLE
obj12_table_ttest <- tableFunction(temp_infLys,"dataSet_bio","COMB_midPoint1_x","t.test")
obj12_table_utest <- tableFunction(temp_infLys,"dataSet_bio","COMB_midPoint1_x","wilcox.test")
obj12_table_difMean <- tableFunction(temp_infLys,"dataSet_bio","COMB_midPoint1_x","difMean")
obj12_table_effectSize <- tableFunction(temp_infLys,"dataSet_bio","COMB_midPoint1_x","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj12_table_ttest[[1]],"**\n\n",obj12_table_ttest[[2]],"\n"))
knitr::kable(obj12_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj12_table_utest[[1]],"**\n\n",obj12_table_utest[[2]],"\n"))
knitr::kable(obj12_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj12_table_difMean[[1]],"**\n\n",obj12_table_difMean[[2]],"\n"))
knitr::kable(obj12_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj12_table_effectSize[[1]],"**\n\n",obj12_table_effectSize[[2]],"\n"))
knitr::kable(obj12_table_effectSize[[3]])
#*****************************************


#*****************************************
fig13<-singleVariableHistogramFunction(DF=temp_inf,
                                       x_axis="COMB_midPoint1_x",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime)

#print(fig13)
saveFigure(figureObj = fig13,name = "midPoint_Estimate_inf", dataSetName = dataSetName, folderName=folderName)


# TABLE
obj13_table_ttest <- tableFunction(temp_inf,"dataSet_bio","COMB_midPoint1_x","t.test")
obj13_table_utest <- tableFunction(temp_inf,"dataSet_bio","COMB_midPoint1_x","wilcox.test")
obj13_table_difMean <- tableFunction(temp_inf,"dataSet_bio","COMB_midPoint1_x","difMean")
obj13_table_effectSize <- tableFunction(temp_inf,"dataSet_bio","COMB_midPoint1_x","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj13_table_ttest[[1]],"**\n\n",obj13_table_ttest[[2]],"\n"))
knitr::kable(obj13_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj13_table_utest[[1]],"**\n\n",obj13_table_utest[[2]],"\n"))
knitr::kable(obj13_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj13_table_difMean[[1]],"**\n\n",obj13_table_difMean[[2]],"\n"))
knitr::kable(obj13_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj13_table_effectSize[[1]],"**\n\n",obj13_table_effectSize[[2]],"\n"))
knitr::kable(obj13_table_effectSize[[3]])
#*****************************************


# ---- Midpoint Associated figures and tables WIDE ----
fig_Midpoint<-cowplot::plot_grid(fig11,fig12, fig13, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_Midpoint)


# ---- Midpoint Associated violin ----
fig15a<-ggplot(temp_inf_infLys, aes(x=dataSet_bio, y=COMB_midPoint1_x)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Midpoint")+
  ggtitle("Midpoint Distribution Violin (inf + infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=12),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig15a)


fig15b<-ggplot(temp_infLys, aes(x=dataSet_bio, y=COMB_midPoint1_x)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Midpoint")+
  ggtitle("Midpoint Distribution Violin (infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig15b)

fig15c<-ggplot(temp_inf, aes(x=dataSet_bio, y=COMB_midPoint1_x)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Midpoint")+
  ggtitle("Midpoint Distribution Violin (infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig15c)

fig_midPoint_Violin<-cowplot::plot_grid(fig15a,fig15b, fig15c, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_midPoint_Violin)









# ---- Infection time Associated figures and tables ----


#*****************************************
# Scaling Related
fig16a<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                        x_axis="COMB_incrementTime",
                                        colVec=colVec01, maximumTime=maximumTime)

fig17a<-singleVariableHistogramFunction(DF=temp_infLys,
                                        x_axis="COMB_incrementTime",
                                        colVec=colVec01, maximumTime=maximumTime)

fig18a<-singleVariableHistogramFunction(DF=temp_inf,
                                        x_axis="COMB_incrementTime",
                                        colVec=colVec01, maximumTime=maximumTime)

xRange<-ggplot_build(fig16a)$layout$panel_ranges[[1]]$x.range
xRange_16a<-1.0*xRange[2]
xRange<-ggplot_build(fig17a)$layout$panel_ranges[[1]]$x.range
xRange_17a<-1.0*xRange[2]
xRange<-ggplot_build(fig18a)$layout$panel_ranges[[1]]$x.range
xRange_18a<-1.0*xRange[2]

xRange=max(xRange_16a, xRange_17a, xRange_18a)
#*****************************************


#*****************************************
fig16<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                       x_axis="COMB_incrementTime",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime,
                                       Adjust = 1)

#print(fig16)
saveFigure(figureObj = fig16,name = "infectiontime_Estimate_inf_infLys", dataSetName = dataSetName, folderName=folderName)


# TABLE
obj16_table_ttest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_incrementTime","t.test")
obj16_table_utest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_incrementTime","wilcox.test")
obj16_table_difMean <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_incrementTime","difMean")
obj16_table_effectSize <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_incrementTime","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj16_table_ttest[[1]],"**\n\n",obj16_table_ttest[[2]],"\n"))
knitr::kable(obj16_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj16_table_utest[[1]],"**\n\n",obj16_table_utest[[2]],"\n"))
knitr::kable(obj16_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj16_table_difMean[[1]],"**\n\n",obj16_table_difMean[[2]],"\n"))
knitr::kable(obj16_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj16_table_effectSize[[1]],"**\n\n",obj16_table_effectSize[[2]],"\n"))
knitr::kable(obj16_table_effectSize[[3]])
#*****************************************


#*****************************************
fig17<-singleVariableHistogramFunction(DF=temp_infLys,
                                       x_axis="COMB_incrementTime",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime,
                                       Adjust = 1)
#print(fig17)
saveFigure(figureObj = fig17,name = "infectiontime_Estimate_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj17_table_ttest <- tableFunction(temp_infLys,"dataSet_bio","COMB_incrementTime","t.test")
obj17_table_utest <- tableFunction(temp_infLys,"dataSet_bio","COMB_incrementTime","wilcox.test")
obj17_table_difMean <- tableFunction(temp_infLys,"dataSet_bio","COMB_incrementTime","difMean")
obj17_table_effectSize <- tableFunction(temp_infLys,"dataSet_bio","COMB_incrementTime","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj17_table_ttest[[1]],"**\n\n",obj17_table_ttest[[2]],"\n"))
knitr::kable(obj17_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj17_table_utest[[1]],"**\n\n",obj17_table_utest[[2]],"\n"))
knitr::kable(obj17_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj17_table_difMean[[1]],"**\n\n",obj17_table_difMean[[2]],"\n"))
knitr::kable(obj17_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj17_table_effectSize[[1]],"**\n\n",obj17_table_effectSize[[2]],"\n"))
knitr::kable(obj17_table_effectSize[[3]])
#*****************************************


#*****************************************
fig18<-singleVariableHistogramFunction(DF=temp_inf,
                                       x_axis="COMB_incrementTime",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime,
                                       Adjust = 1)
#print(fig18)
saveFigure(figureObj = fig18,name = "infectiontime_Estimate_inf", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj18_table_ttest <- tableFunction(temp_inf,"dataSet_bio","COMB_incrementTime","t.test")
obj18_table_utest <- tableFunction(temp_inf,"dataSet_bio","COMB_incrementTime","wilcox.test")
obj18_table_difMean <- tableFunction(temp_inf,"dataSet_bio","COMB_incrementTime","difMean")
obj18_table_effectSize <- tableFunction(temp_inf,"dataSet_bio","COMB_incrementTime","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj18_table_ttest[[1]],"**\n\n",obj18_table_ttest[[2]],"\n"))
knitr::kable(obj18_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj18_table_utest[[1]],"**\n\n",obj18_table_utest[[2]],"\n"))
knitr::kable(obj18_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj18_table_difMean[[1]],"**\n\n",obj18_table_difMean[[2]],"\n"))
knitr::kable(obj18_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj18_table_effectSize[[1]],"**\n\n",obj18_table_effectSize[[2]],"\n"))
knitr::kable(obj18_table_effectSize[[3]])
#*****************************************


# ---- Infection time Associated figures and tables WIDE ----
fig_InfectionTime<-cowplot::plot_grid(fig16,fig17, fig18, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_InfectionTime)

# ---- Infection time associated violin ----
fig16b<-ggplot(temp_inf_infLys, aes(x=dataSet_bio, y=COMB_incrementTime)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Infection Time")+
  ggtitle("Infection Time Distribution Violin (inf + infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig16b)


fig17b<-ggplot(temp_infLys, aes(x=dataSet_bio, y=COMB_incrementTime)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Infection Time")+
  ggtitle("Infection Time Distribution Violin (infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig17b)

fig18b<-ggplot(temp_inf, aes(x=dataSet_bio, y=COMB_incrementTime)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Infection Time")+
  ggtitle("Infection Time Distribution Violin (infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig18b)

fig_infectionTime_Violin<-cowplot::plot_grid(fig16b, fig17b, fig18b, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_infectionTime_Violin)





# ---- Start point Associated figures and tables ----

# Scaling Related
fig20a<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                        x_axis="COMB_startPoint_x",
                                        colVec=colVec01, maximumTime=maximumTime)

fig21a<-singleVariableHistogramFunction(DF=temp_infLys,
                                        x_axis="COMB_startPoint_x",
                                        colVec=colVec01, maximumTime=maximumTime)

fig22a<-singleVariableHistogramFunction(DF=temp_inf,
                                        x_axis="COMB_startPoint_x",
                                        colVec=colVec01, maximumTime=maximumTime)

xRange<-ggplot_build(fig20a)$layout$panel_ranges[[1]]$x.range
xRange_20a<-1.0*xRange[2]
xRange<-ggplot_build(fig21a)$layout$panel_ranges[[1]]$x.range
xRange_21a<-1.0*xRange[2]
xRange<-ggplot_build(fig22a)$layout$panel_ranges[[1]]$x.range
xRange_22a<-1.0*xRange[2]

xRange=max(xRange_20a, xRange_21a, xRange_22a)
#*****************************************


#*****************************************
# Generating start point column
#*****************************************
fig20<-singleVariableHistogramFunction(DF=temp_inf_infLys,
                                       x_axis="COMB_startPoint_x",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime,
                                       Adjust = 1)

#print(fig20)
saveFigure(figureObj = fig20, name = "startpoint_Estimate_inf_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj20_table_ttest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_startPoint_x","t.test")
obj20_table_utest <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_startPoint_x","wilcox.test")
obj20_table_difMean <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_startPoint_x","difMean")
obj20_table_effectSize <- tableFunction(temp_inf_infLys,"dataSet_bio","COMB_startPoint_x","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj20_table_ttest[[1]],"**\n\n",obj20_table_ttest[[2]],"\n"))
knitr::kable(obj20_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj20_table_utest[[1]],"**\n\n",obj20_table_utest[[2]],"\n"))
knitr::kable(obj20_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj20_table_difMean[[1]],"**\n\n",obj20_table_difMean[[2]],"\n"))
knitr::kable(obj20_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj20_table_effectSize[[1]],"**\n\n",obj20_table_effectSize[[2]],"\n"))
knitr::kable(obj20_table_effectSize[[3]])
#*****************************************


#*****************************************
fig21<-singleVariableHistogramFunction(DF=temp_infLys,
                                       x_axis="COMB_startPoint_x",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime,
                                       Adjust = 1)

#print(fig21)
saveFigure(figureObj = fig21, name = "startpoint_Estimate_infLys", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj21_table_ttest <- tableFunction(temp_infLys,"dataSet_bio","COMB_startPoint_x","t.test")
obj21_table_utest <- tableFunction(temp_infLys,"dataSet_bio","COMB_startPoint_x","wilcox.test")
obj21_table_difMean <- tableFunction(temp_infLys,"dataSet_bio","COMB_startPoint_x","difMean")
obj21_table_effectSize <- tableFunction(temp_infLys,"dataSet_bio","COMB_startPoint_x","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj21_table_ttest[[1]],"**\n\n",obj21_table_ttest[[2]],"\n"))
knitr::kable(obj21_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj21_table_utest[[1]],"**\n\n",obj21_table_utest[[2]],"\n"))
knitr::kable(obj21_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj21_table_difMean[[1]],"**\n\n",obj21_table_difMean[[2]],"\n"))
knitr::kable(obj21_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj21_table_effectSize[[1]],"**\n\n",obj21_table_effectSize[[2]],"\n"))
knitr::kable(obj21_table_effectSize[[3]])
#*****************************************


#*****************************************
fig22<-singleVariableHistogramFunction(DF=temp_inf,
                                       x_axis="COMB_startPoint_x",
                                       x_limits=c(-0.06*xRange,xRange),
                                       colVec=colVec01, maximumTime=maximumTime,
                                       Adjust = 1)

#print(fig22)
saveFigure(figureObj = fig22, name = "startpoint_Estimate_inf", dataSetName = dataSetName, folderName=folderName)

# TABLE
obj22_table_ttest <- tableFunction(temp_inf,"dataSet_bio","COMB_startPoint_x","t.test")
obj22_table_utest <- tableFunction(temp_inf,"dataSet_bio","COMB_startPoint_x","wilcox.test")
obj22_table_difMean <- tableFunction(temp_inf,"dataSet_bio","COMB_startPoint_x","difMean")
obj22_table_effectSize <- tableFunction(temp_inf,"dataSet_bio","COMB_startPoint_x","t.test.effectSize")


knitr::asis_output(paste0("\n\n\n**",obj22_table_ttest[[1]],"**\n\n",obj22_table_ttest[[2]],"\n"))
knitr::kable(obj22_table_ttest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj22_table_utest[[1]],"**\n\n",obj22_table_utest[[2]],"\n"))
knitr::kable(obj22_table_utest[[3]])

knitr::asis_output(paste0("\n\n\n**",obj22_table_difMean[[1]],"**\n\n",obj22_table_difMean[[2]],"\n"))
knitr::kable(obj22_table_difMean[[3]])

knitr::asis_output(paste0("\n\n\n**",obj22_table_effectSize[[1]],"**\n\n",obj22_table_effectSize[[2]],"\n"))
knitr::kable(obj22_table_effectSize[[3]])
#*****************************************


# ---- start point Associated figures and tables WIDE ----
fig_StartPoint<-cowplot::plot_grid(fig20, fig21, fig22, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_StartPoint)

# ---- Start point associated violin ----
fig20b<-ggplot(temp_inf_infLys, aes(x=dataSet_bio, y=COMB_startPoint_x)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Start Point")+
  ggtitle("Start Point Distribution Violin (inf + infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig20b)


fig21b<-ggplot(temp_infLys, aes(x=dataSet_bio, y=COMB_startPoint_x)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Start Point")+
  ggtitle("Start Point Distribution Violin (infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig21b)

fig22b<-ggplot(temp_inf, aes(x=dataSet_bio, y=COMB_startPoint_x)) +
  geom_violin(fill="grey50") +
  #geom_point(aes(x=as.factor(MOI), y=mean(M1_Estimate))) +
  theme_bw()+
  ylab("Start Point")+
  ggtitle("Start Point Distribution Violin (infLys)")+
  theme( axis.text.x=element_text(angle = 45, hjust = 1, vjust=1, size=12),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )+
  ylim(0,xRange)
#print(fig22b)

fig_startPoint_Violin<-cowplot::plot_grid(fig20b, fig21b, fig22b, labels = c("A", "B", "C"), ncol = 3, nrow = 1.5, scale = .95)
print(fig_startPoint_Violin)









#### 2D

# ---- Midpoint vs Slope Associated figures and tables ----

#***********************************************
# Find Figure Limits
fig11a<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="point",
                                       x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                       x_limits=NULL, y_limits=NULL,
                                       colVec=colVec01, maximumTime=maximumTime,
                                       sizeForLines=2, sizeForPoints=2)

fig12a<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="point",
                                       x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                       x_limits=NULL, y_limits=NULL,
                                       colVec=colVec01, maximumTime=maximumTime,
                                       sizeForLines=2, sizeForPoints=2)


xRange<-ggplot_build(fig11a)$layout$panel_ranges[[1]]$x.range
xRange_11a<-1.0*xRange[2]
xRange<-ggplot_build(fig12a)$layout$panel_ranges[[1]]$x.range
xRange_12a<-1.0*xRange[2]

xRange=max(xRange_11a,xRange_12a)


yRange<-ggplot_build(fig11a)$layout$panel_ranges[[1]]$y.range
yRange_11a<-1.0*yRange[2]
yRange<-ggplot_build(fig12a)$layout$panel_ranges[[1]]$y.range
yRange_12a<-1.0*yRange[2]

yRange=max(yRange_11a,yRange_12a)
#***********************************************


#***********************************************
# Generate Figures
fig11<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="line",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig11)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","midpointVSslope_inf_infLys_line.svg"),
                   plot = fig11, ncol = 2)

fig12<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="line",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig12)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","midpointVSslope_infLys_line.svg"),
                   plot = fig12, ncol = 2)

fig13<-doubleVariablePlottingFunction(DF=temp_inf, line_vs_point="line",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig13)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","midpointVSslope_inf_line.svg"),
                   plot = fig13, ncol = 2)

fig14<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="point",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig14,name = "midPoint_slope_inf_infLys", dataSetName = dataSetName, folderName=folderName)

fig15<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="point",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig15,name = "midPoint_slope_infLys", dataSetName = dataSetName, folderName=folderName)

fig16<-doubleVariablePlottingFunction(DF=temp_inf, line_vs_point="point",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_slope1",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig16,name = "midPoint_slope_inf", dataSetName = dataSetName, folderName=folderName)
#***********************************************


#***********************************************
# ---- Midpoint vs Slope Associated figures and tables WIDE ----
fig_midpoint_slope<-cowplot::plot_grid(fig11, fig14, fig12, fig15, fig13, fig16,
                                       labels = c("A", "B", "C", "D", "E", "F"), nrow = 3)
print(fig_midpoint_slope)
#***********************************************






# ---- Slope vs Maximum Associated figures and tables ----

#***********************************************
# Find Figure Limits
fig21a<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="point",
                                       x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                       x_limits=NULL, y_limits=NULL,
                                       colVec=colVec01, maximumTime=maximumTime,
                                       sizeForLines=2, sizeForPoints=2)

fig22a<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="point",
                                       x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                       x_limits=NULL, y_limits=NULL,
                                       colVec=colVec01, maximumTime=maximumTime,
                                       sizeForLines=2, sizeForPoints=2)


xRange<-ggplot_build(fig21a)$layout$panel_ranges[[1]]$x.range
xRange_21a<-1.0*xRange[2]
xRange<-ggplot_build(fig22a)$layout$panel_ranges[[1]]$x.range
xRange_22a<-1.0*xRange[2]

xRange=max(xRange_21a,xRange_22a)


yRange<-ggplot_build(fig21a)$layout$panel_ranges[[1]]$y.range
yRange_21a<-1.0*yRange[2]
yRange<-ggplot_build(fig22a)$layout$panel_ranges[[1]]$y.range
yRange_22a<-1.0*yRange[2]

yRange=max(yRange_21a,yRange_22a)
#***********************************************


#***********************************************
# Generate Figures
fig21<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="line",
                                      x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig21)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","slopeVSmaximum_inf_infLys_line.svg"),
                   plot = fig21, ncol = 2)

fig22<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="line",
                                      x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
#print(fig22)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","slopeVSmaximum_infLys_line.svg"),
                   plot = fig22, ncol = 2)

fig23<-doubleVariablePlottingFunction(DF=temp_inf, line_vs_point="line",
                                      x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
#print(fig23)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","slopeVSmaximum_inf_line.svg"),
                   plot = fig23, ncol = 2)

fig24<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="point",
                                      x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig24,name = "slope_maximum_inf_infLys", dataSetName = dataSetName, folderName=folderName)

fig25<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="point",
                                      x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig25,name = "slope_maximum_infLys", dataSetName = dataSetName, folderName=folderName)

fig26<-doubleVariablePlottingFunction(DF=temp_inf, line_vs_point="point",
                                      x_axis="COMB_slope1", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig26,name = "slope_maximum_inf", dataSetName = dataSetName, folderName=folderName)
#***********************************************


#***********************************************
# ---- Slope vs Maximum Associated figures and tables WIDE----
fig_slope_maximum<-cowplot::plot_grid(fig21, fig24, fig22, fig25, fig23, fig26,
                                      labels = c("A", "B", "C", "D", "E", "F"), nrow = 3)
print(fig_slope_maximum)
#***********************************************




# ---- Midpoint vs Maximum Associated figures and tables ----
#***********************************************
# Find Figure Limits
fig31a<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="point",
                                       x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                       x_limits=NULL, y_limits=NULL,
                                       colVec=colVec01, maximumTime=maximumTime,
                                       sizeForLines=2, sizeForPoints=2)

fig32a<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="point",
                                       x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                       x_limits=NULL, y_limits=NULL,
                                       colVec=colVec01, maximumTime=maximumTime,
                                       sizeForLines=2, sizeForPoints=2)


xRange<-ggplot_build(fig31a)$layout$panel_ranges[[1]]$x.range
xRange_31a<-1.0*xRange[2]
xRange<-ggplot_build(fig32a)$layout$panel_ranges[[1]]$x.range
xRange_32a<-1.0*xRange[2]

xRange=max(xRange_31a,xRange_32a)


yRange<-ggplot_build(fig31a)$layout$panel_ranges[[1]]$y.range
yRange_31a<-1.0*yRange[2]
yRange<-ggplot_build(fig32a)$layout$panel_ranges[[1]]$y.range
yRange_32a<-1.0*yRange[2]

yRange=max(yRange_31a,yRange_32a)
#***********************************************


#***********************************************
# Generate Figures
fig31<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="line",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig31)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","midPointVSmaximum_inf_infLys_line.svg"),
                   plot = fig31, ncol = 2)

fig32<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="line",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig32)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","midPointVSmaximum_infLys_line.svg"),
                   plot = fig32, ncol = 2)

fig33<-doubleVariablePlottingFunction(DF=temp_inf, line_vs_point="line",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig33)
cowplot::save_plot(filename = paste0("./",folderName,"/","SCV_",dataSetName,"_","midPointVSmaximum_inf_line.svg"),
                   plot = fig33, ncol = 2)

fig34<-doubleVariablePlottingFunction(DF=temp_inf_infLys, line_vs_point="point",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)

saveFigure2(figureObj = fig34,name = "midPoint_maximum_inf_infLys", dataSetName = dataSetName, folderName=folderName)

fig35<-doubleVariablePlottingFunction(DF=temp_infLys, line_vs_point="point",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig35)

saveFigure2(figureObj = fig35,name = "midPoint_maximum_infLys", dataSetName = dataSetName, folderName=folderName)

fig36<-doubleVariablePlottingFunction(DF=temp_inf, line_vs_point="point",
                                      x_axis="COMB_midPoint1_x", y_axis="COMB_maximum_y",
                                      x_limits=c(-0.06*xRange,xRange*1.06), y_limits=c(-0.06*yRange,yRange*1.06),
                                      colVec=colVec01, maximumTime=maximumTime,
                                      sizeForLines=2, sizeForPoints=2)
print(fig36)

saveFigure2(figureObj = fig36,name = "midPoint_maximum_inf", dataSetName = dataSetName, folderName=folderName)
#***********************************************


# ---- Midpoint vs Maximum Associated figures and tables WIDE ----
fig_midpoint_maximum<-cowplot::plot_grid(fig31,fig34, fig32, fig35, fig33, fig36,
                                         labels = c("A", "B", "C", "D", "E", "F"), nrow = 3)
print(fig_midpoint_maximum)
