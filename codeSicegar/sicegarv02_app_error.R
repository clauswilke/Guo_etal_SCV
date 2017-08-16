# Sicegar V0.2 application


###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
seedNo=14159
set.seed(seedNo)
###*****************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/single_cell_virology/codeSicegarV02/"))} # mac computer
###*****************************


###*****************************
# REQUIRED LIBRARIES
# Data tracking
require("sicegar")
require("DESeq2")
require("dplyr")
require("tidyr")

# Graphing
require("ggplot2")
require("cowplot")

# Text manipulation
require("stringr")
###*****************************


###*****************************
#Load Functions
###*****************************


###*****************************
#Load data
# Load data and nicely format it
dataFileName="SCV114";
load(paste0('../tidyData/tidy',dataFileName,'.RDA'))
assign(x = "currentDf", value = get(paste0(dataFileName,"TidyGrouped")))
###*****************************


###*****************************
# Sicegar v0.2
sicegar_fitFunction<-function(data)
{
  dataInput1=data.frame(time=data$time, intensity=data$intensity)
  dataInputName=unique(data$Cell)
  #print(dataInputName)
  sicegar::fitAndCategorize(dataInput = dataInput1,
                            dataInputName = dataInputName) ->fitObj
  output = fitObj$summaryVector
  return(output)
}

currentDf %>%
  dplyr::group_by(Cell,Mutation,MOI,dataSet) %>%
  dplyr::select(time = hpi, intensity = nRat_a) %>%
  dplyr::do(resultsSicegar=sicegar_fitFunction(.))->currentDF_analyzeResults

currentDF_analyzeResults %>%
  dplyr::group_by(Cell,Mutation,MOI,dataSet) %>%
  dplyr::do(.,data.frame(t(unlist(.$resultsSicegar,use.names = TRUE)))) ->currentDF_analyzeResults2
###*****************************


###*****************************
# Combine Results
currentDf %>% dplyr::filter(max(nRat_a)>0.3)%>%dplyr::summarize()

number=10

currentDf%>%
  dplyr::filter(Cell==number) %>%
  dplyr::group_by()%>%
  dplyr::select(time=hpi, intensity=nRat_a)->data

data%>%
  sicegar::fitAndCategorize(dataInput = .,
                            n_runs_min_sm = 100, n_runs_max_sm = 1000,
                            n_runs_min_dsm = 100,n_runs_max_dsm = 1000) ->fitObj

fitObj%>%
  .$summaryVector %>%
  utils::str(.)

decision<-fitObj$summaryVector$decision

if(decision=="sigmoidal")
{
  sicegar::figureModelCurves(dataInput = data,
                             sigmoidalFitVector = fitObj$sigmoidalModel,
                             showParameterRelatedLines = T)
}

if(decision=="double_sigmoidal")
{
  sicegar::figureModelCurves(dataInput = data,
                             doubleSigmoidalFitVector = fitObj$doubleSigmoidalModel,
                             showParameterRelatedLines = T)
}

if(decision %in% c("ambiguous", "no_signal"))
{
  sicegar::figureModelCurves(dataInput = data)
}
###*****************************
