## Generating Data Table By Using Tidy_R

# Initial Command to Reset the System
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
# setwd('C:/Users/kelvin/Desktop/Main Git Folder/single_cell_virology')    #Laptop
# setwd('C:/Users/Eda Umut/Desktop/GIt Projects/single_cell_virology') # home computer
setwd('/Users/umut/GitHub/Guo_etal_SCV/codeSicegar/') # mac computer


#Load Libraries
require(ggplot2)
require(plyr)
require(dplyr)
require(reshape2)
require(tidyr)

# Data Set and data properties selection

#listFileNames=c("SCV028", "SCV029", "SCV030") # the correct data sets
#listMOI_level=c("50","500", "5000") # write the correct MOI's (must be numbers)
#listMutation_Value=c("WT", "WT", "WT") # write the mutant type

#listFileNames=c("SCV031", "SCV032") # the correct data sets
#listMOI_level=c("500", "500") # write the correct MOI's (must be numbers)
#listMutation_Value=c("WT", "WT") # write the mutant type

#listFileNames=c("SCV035", "SCV036", "SCV037") # the correct data sets
#listMOI_level=c("500","500", "500") # write the correct MOI's (must be numbers)
#listMutation_Value=c("WT", "WT", "WT") # write the mutant type

#listFileNames=c("SCV038", "SCV039", "SCV040", "SCV041") # the correct data sets
#listMOI_level=c("50","5000", "50", "5000") # write the correct MOI's (must be numbers)
#listMutation_Value=c("WT", "WT", "HR", "HR") # write the mutant type

listFileNames=c("SCV096", "SCV097", "SCV098", "SCV099") # the correct data sets
listMOI_level=c("2000","2000", "2000", "2000") # write the correct MOI's (must be numbers)
listMutation_Value=c("WT", "WT", "HR", "HR") # write the mutant type

for (i in 1:length(listFileNames) ) {
  dataFileName=listFileNames[i]
  MOI_level=listMOI_level[i]
  Mutation_Value=listMutation_Value[i]

  Initialdata <- read.csv(paste0('../rawData/',dataFileName,'.csv'))

  newColNames <-c("hpi",paste0("X",as.character(c(1:(length(colnames(Initialdata))-1)))))
  colnames(Initialdata)<-newColNames
  finalPosition<-length(colnames(Initialdata));
  Initialdata_nRat_a <-Initialdata[1:43,]
  Initialdata_Infection <-Initialdata[68,];
  Initialdata_Infection <- replace(Initialdata_Infection,
                                   1:length(Initialdata_Infection),
                                   NA) # additional Line
  row.names(Initialdata_Infection) <- NULL

  InitialdataDDOM <-as.data.frame(t(Initialdata[44:47,2:finalPosition]));
  colnames(InitialdataDDOM) <- c("Dirt", "Dried_wells", "Out_of_focus", "Multiple_infections")
  InitialdataDDOM %>%
    dplyr::mutate(Cell=as.numeric(sub("X","",rownames(InitialdataDDOM)))) ->InitialdataDDOM
  rownames(InitialdataDDOM)<-NULL

  Initialdata_nRat_a<-gather(Initialdata_nRat_a,Cell,nRat_a,2:finalPosition)
  Initialdata_nRat_a$nRat_a=as.numeric(Initialdata_nRat_a$nRat_a)+runif(nrow(Initialdata_nRat_a), min=0, max=1)*10^-16 # add noise
  Initialdata_Infection<-gather(Initialdata_Infection,Cell,Infection,2:finalPosition)
  Initialdata_Infection<-Initialdata_Infection[,2:3]
  Initialdata<-Initialdata_nRat_a[,2:1]
  Initialdata$GFP<-NA
  Initialdata$BG<-NA
  Initialdata<-inner_join(Initialdata,Initialdata_Infection)
  Initialdata$Mutation<-as.factor(rep(c(Mutation_Value),each=nrow(Initialdata)))
  Initialdata$MOI<-as.factor(rep(c(MOI_level),each=nrow(Initialdata)))
  Initialdata$dataSet<-as.factor(rep(c(dataFileName),each=nrow(Initialdata)))
  Initialdata$nGFP<-NA
  Initialdata<-inner_join(Initialdata,Initialdata_nRat_a)

  Initialdata$Cell<-sub("X", "", Initialdata$Cell)

  Initialdata$Cell <- as.numeric(Initialdata$Cell)
  Initialdata$hpi <- as.numeric(as.character(Initialdata$hpi))
  Initialdata$MOI <- factor(Initialdata$MOI)
  Initialdata$Infection <- factor(Initialdata$Infection)
  Initialdata$nRat_a <- as.numeric(Initialdata$nRat_a)

  InitialdataGrouped <- dplyr::group_by(Initialdata, MOI, Mutation, Cell, dataSet) # Grouping Calculations
  InitialdataGrouped <- dplyr::mutate(InitialdataGrouped,
                                      nRat_b=(nRat_a-min(nRat_a))/(max(nRat_a)-min(nRat_a))
  )



  variableName01 <-paste0(dataFileName,"TidyGrouped")
  assign(variableName01,InitialdataGrouped)
  variableName02 <-paste0(dataFileName,"TidyGrouped_DDOM")
  assign(variableName02,InitialdataDDOM)

  remove(Initialdata,
         Initialdata_nRat_a,
         Initialdata_Infection,
         InitialdataGrouped,
         InitialdataDDOM)

  fileName=paste0("tidy",dataFileName,".RDA")
  save(list=c(variableName01),file=paste0("../tidyData/",fileName))

  fileName=paste0("tidy",dataFileName,"_DDOM.RDA")
  save(list=c(variableName02),file=paste0("../tidyData/",fileName))
}
