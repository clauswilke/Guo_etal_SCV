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
{setwd(paste0("/Users/umut/GitHub/Guo_etal_SCV/codeSicegar/"))} # mac computer
###*****************************


###*****************************
# REQUIRED LIBRARIES
# Data tracking
require("sicegar")
require("DESeq2")
require("tidyverse")

# Graphing
require("ggplot2")
require("cowplot")

# Text manipulation
require("stringr")
###*****************************


###*****************************
#Load Functions
###*****************************

#********************************************
# Chose data name
# simply change this to pick up correct files to analyse

listFileNames=c("SCV099")
min_max_val = 0.4
#********************************************

for (counter02 in 1:length(listFileNames) )
{
  ###*****************************
  #Load data
  # Load data and nicely format it
  dataFileName = listFileNames[counter02];
  load(paste0('../tidyData/tidy',dataFileName,'.RDA'))
  assign(x = "currentDf", value = get(paste0(dataFileName,"TidyGrouped")))
  ###*****************************
  
  # currentDf %>%
  #   dplyr::filter(Cell %in% seq(1,100))->currentDf
  
  ###*****************************
  # Sicegar v0.2
  sicegar_fitFunction<-function(data)
  {
    dataInput1=data.frame(time=data$time, intensity=data$intensity)
    dataInputName=unique(data$Cell)
    print(dataInputName)
    sicegar::fitAndCategorize(dataInput = dataInput1,
                              dataInputName = dataInputName,
                              threshold_minimum_for_intensity_maximum = min_max_val) ->fitObj
    
    # Data Scaling Parameters
    dataScalingParametersDF <- data.frame(t(fitObj$normalizedInput$dataScalingParameters))
    colnames(dataScalingParametersDF) <- paste0("SCALING_", colnames(dataScalingParametersDF))
    
    # Sigmoidal
    if(typeof(fitObj$sigmoidalModel)!="logical")
    {
      sigmoidalDF <- data.frame(fitObj$sigmoidalModel)
      sigmoidalDF[grep(pattern = "dataScalingParameters", colnames(sigmoidalDF))]<-NULL
      colnames(sigmoidalDF) <- paste0("SM_", colnames(sigmoidalDF))
    }
    
    # Double Sigmoidal
    if(typeof(fitObj$doubleSigmoidalModel)!="logical")
    {
      doublesigmoidalDF <- data.frame(fitObj$doubleSigmoidalModel)
      doublesigmoidalDF[grep(pattern = "dataScalingParameters", colnames(doublesigmoidalDF))]<-NULL
      colnames(doublesigmoidalDF) <- paste0("DSM_", colnames(doublesigmoidalDF))
    }
    
    #Decision Process
    if(typeof(fitObj$decisionProcess)!="logical")
    {
      decisionProcessDF <- as.data.frame(t(fitObj$decisionProcess))
      colnames(decisionProcessDF) <- paste0("DEC_", colnames(decisionProcessDF))
    }
    
    # Summary File
    summaryDF <- data.frame(t(fitObj$summaryVector))
    summaryDF[,c(1,2)] <- NULL
    if(decisionProcessDF["DEC_decision"]=="sigmoidal")
    {
      summaryDF %>%
        dplyr::rename(midPoint1_x = midPoint_x,
                      midPoint1_y = midPoint_y,
                      slope1 = slope)->summaryDF
    }
    
    if(ncol(summaryDF)!=0)
    {
      colnames(summaryDF) <- paste0("COMB_", colnames(summaryDF)) 
    }
    
    if(exists("doublesigmoidalDF") & exists("sigmoidalDF"))
    {
      output <- dplyr::bind_cols(dataScalingParametersDF, 
                                 sigmoidalDF, 
                                 doublesigmoidalDF, 
                                 decisionProcessDF,
                                 summaryDF)
    }
    else
    {
      output <- dplyr::bind_cols(dataScalingParametersDF, 
                                 decisionProcessDF,
                                 summaryDF)
    }
    
    names(output) <- sub("\\.", "_",names(output))
    output2 <- purrr::flatten(as.vector(output))
    names(output2) <- names(output)
    return(output2)
  }
  
  
  
  currentDf %>%
    dplyr::group_by(Cell,Mutation,MOI,dataSet) %>%
    dplyr::select(time = hpi, intensity = nRat_a) %>%
    dplyr::do(resultsSicegar=sicegar_fitFunction(.))->currentDF_analyzeResults
  
  currentDF_analyzeResults %>%
    dplyr::group_by(Cell,Mutation,MOI,dataSet) %>%
    dplyr::do(.,as.data.frame(t(unlist((.$resultsSicegar))))) ->currentDF_analyzeResults2
  ###*****************************
  
  
  #******************************************
  # Save Necessary Files
  write.csv(x = currentDF_analyzeResults2,
            file = paste0("../detailedSicegarResults/",
                          "sicegar_output_",
                          "currentDf_"
                          ,dataFileName,
                          ".csv"),
            quote = FALSE, row.names = FALSE)
  
  #******************************************

}
