# the file will generate the tables used in the experiments


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
{setwd(paste0("/Users/umut/GitHub/Guo_etal_SCV/tables/"))} # mac computer
###*****************************


###*****************************
# Necessary libraries
require(tidyverse)
require(lazyeval)
###*****************************

###*****************************
# Install table function
source("../finalResultsRMarkDown/sicegarTableFunction.R")
source("../finalResultsRMarkDown/sicegarTableFunctionWODiag.R")
###*****************************


###*****************************
# Install data
# Read Samples
samp_028_029_030 <- read.csv("../finalResultsRMarkDown/variables_028_029_030/inf_infLys_028_029_030.csv")
samp_029_031_035 <- read.csv("../finalResultsRMarkDown/variables_029_031_035/inf_infLys_029_031_035.csv")
samp_031_032 <- read.csv("../finalResultsRMarkDown/variables_031_032/inf_infLys_031_032.csv")
samp_035_036_037 <- read.csv("../finalResultsRMarkDown/variables_035_036_037/inf_infLys_035_036_037.csv")
samp_038_039_040_041 <- read.csv("../finalResultsRMarkDown/variables_038_039_040_041/inf_infLys_038_039_040_041.csv")
samp_096_097_098_099 <- read.csv("../finalResultsRMarkDown/variables_096_097_098_099/inf_infLys_096_097_098_099.csv")


samp <- dplyr::bind_rows(samp_028_029_030, 
                         samp_031_032, 
                         samp_035_036_037, 
                         samp_038_039_040_041, 
                         samp_096_097_098_099)
###*****************************


###*****************************
# Table_02
#  Mean and standard deviation values for maximum, midpoint, slope, 
# infection time and start point.

samp%>%
  dplyr::group_by(dataSet, dataSet_bio, Mutation, MOI) %>%
  dplyr::summarise("Maximum mean" = mean(COMB_maximum_y),
                   "Maximum sd" = sd(COMB_maximum_y),
                   "Midpoint mean" = mean(COMB_midPoint1_x),
                   "Midpoint sd" = sd(COMB_midPoint1_x),
                   "Slope mean" = mean(COMB_slope1),
                   "Slope sd" = sd(COMB_slope1),
                   "Infection time mean" = mean(COMB_incrementTime),
                   "Infection time sd" = sd(COMB_incrementTime),
                   "Startpoint mean" = mean(COMB_startPoint_x),
                   "Startpoint sd" = sd(COMB_startPoint_x)) -> table_02

ini_df = table_02[setdiff(seq(1,ncol(table_02)),grep(" sd| mean" ,x = colnames(table_02)))]
rounded_df = signif(table_02[grep(" sd| mean" ,x = colnames(table_02))],2)
table_02_round = bind_cols(x = ini_df, y = rounded_df)

write.csv(x = table_02_round, file = "table_02.csv", quote = F, row.names = F, col.names = T)
###*****************************


###*****************************
# P value Table Function
table_p_value <- function(df, classifier)
{
  
  df %>%
    dplyr::group_by_(classifier) %>%
    dplyr::summarise(numSample = n())-> df_summary
  
  print.data.frame(df_summary)
  
  cat("\nAll tests use pairwise comparisons using t tests with pooled SD \n")
  cat("All test use *fdr* as P value adjustment method \n")
  
  # Table_maximum_new
  temp = pairwise.t.test(df$COMB_maximum_y, 
                         df[[classifier]], 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Maximum"
  
  theMatrix <- matrix(nrow = length(df_summary[[classifier]]), ncol = length(df_summary[[classifier]]))
  colnames(theMatrix) <- df_summary[[classifier]]
  rownames(theMatrix) <- df_summary[[classifier]]
  theMatrix[upper.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  theMatrix[lower.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  diag(theMatrix) <- "1.00"
  
  theMatrix %>%
    rbind(.,"Sample size" = as.character(df_summary[["numSample"]])) %>%
    cbind(.,"Sample size" = c(as.character(df_summary[["numSample"]]),"")) -> theMatrix
  
  theMatrix %>%
    cbind(rn=rownames(theMatrix),.) %>%
    cbind("Test parameter" = temp$data.name,.) %>%
    rbind(colnames(.),.) -> theMaximumMatrix
  
  cat("\n",paste0(temp$data.name, "\n"))
  print(temp$p.value)
  
  # # Table_maximum_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_maximum_y", 
  #               testType = "t.test")
  
  # Table_slope_new
  temp = pairwise.t.test(df$COMB_slope1, 
                         df[[classifier]], 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Slope"
  
  theMatrix <- matrix(nrow = length(df_summary[[classifier]]), ncol = length(df_summary[[classifier]]))
  colnames(theMatrix) <- df_summary[[classifier]]
  rownames(theMatrix) <- df_summary[[classifier]]
  theMatrix[upper.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  theMatrix[lower.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  diag(theMatrix) <- "1.00"
  
  theMatrix %>%
    rbind(.,"Sample size" = as.character(df_summary[["numSample"]])) %>%
    cbind(.,"Sample size" = c(as.character(df_summary[["numSample"]]),"")) -> theMatrix
  
  theMatrix %>%
    cbind(rn=rownames(theMatrix),.) %>%
    cbind("Test parameter" = temp$data.name,.) %>%
    rbind(colnames(.),.) -> theSlopeMatrix
  
  cat("\n",paste0(temp$data.name, "\n"))
  print(temp$p.value)
  
  # # Table_slope_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_slope1", 
  #               testType = "t.test")
  
  # Table_midPoint_new
  temp = pairwise.t.test(df$COMB_midPoint1_x, 
                         df[[classifier]], 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Midpoint"
  
  theMatrix <- matrix(nrow = length(df_summary[[classifier]]), ncol = length(df_summary[[classifier]]))
  colnames(theMatrix) <- df_summary[[classifier]]
  rownames(theMatrix) <- df_summary[[classifier]]
  theMatrix[upper.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  theMatrix[lower.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  diag(theMatrix) <- "1.00"
  
  theMatrix %>%
    rbind(.,"Sample size" = as.character(df_summary[["numSample"]])) %>%
    cbind(.,"Sample size" = c(as.character(df_summary[["numSample"]]),"")) -> theMatrix
  
  theMatrix %>%
    cbind(rn=rownames(theMatrix),.) %>%
    cbind("Test parameter" = temp$data.name,.) %>%
    rbind(colnames(.),.) -> theMidpointMatrix
  
  cat("\n",paste0(temp$data.name, "\n"))
  print(temp$p.value)
  
  # # Table_midPoint_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_midPoint1_x", 
  #               testType = "t.test")
  
  # Table_startPoint_new
  temp = pairwise.t.test(df$COMB_startPoint_x, 
                         df[[classifier]], 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Start point"
  
  theMatrix <- matrix(nrow = length(df_summary[[classifier]]), ncol = length(df_summary[[classifier]]))
  colnames(theMatrix) <- df_summary[[classifier]]
  rownames(theMatrix) <- df_summary[[classifier]]
  theMatrix[upper.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  theMatrix[lower.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  diag(theMatrix) <- "1.00"
  
  theMatrix %>%
    rbind(.,"Sample size" = as.character(df_summary[["numSample"]])) %>%
    cbind(.,"Sample size" = c(as.character(df_summary[["numSample"]]),"")) -> theMatrix
  
  theMatrix %>%
    cbind(rn=rownames(theMatrix),.) %>%
    cbind("Test parameter" = temp$data.name,.) %>%
    rbind(colnames(.),.) -> theStartpointMatrix
  
  cat("\n",paste0(temp$data.name, "\n"))
  print(temp$p.value)
  
  # # Table_startPoint_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_startPoint_x", 
  #               testType = "t.test")
  
  # Table_incrementTime_new
  temp = pairwise.t.test(df$COMB_incrementTime, 
                         df[[classifier]], 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Infection time"
  
  theMatrix <- matrix(nrow = length(df_summary[[classifier]]), ncol = length(df_summary[[classifier]]))
  colnames(theMatrix) <- df_summary[[classifier]]
  rownames(theMatrix) <- df_summary[[classifier]]
  theMatrix[upper.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  theMatrix[lower.tri(x = theMatrix, diag = FALSE)] <- as.character(signif(temp$p.value[lower.tri(temp$p.value, diag=TRUE)],2))
  diag(theMatrix) <- "1.00"
  
  theMatrix %>%
    rbind(.,"Sample size" = as.character(df_summary[["numSample"]])) %>%
    cbind(.,"Sample size" = c(as.character(df_summary[["numSample"]]),"")) -> theMatrix
  
  theMatrix %>%
    cbind(rn=rownames(theMatrix),.) %>%
    cbind("Test parameter" = temp$data.name,.) %>%
    rbind(colnames(.),.) -> theInfectiontimeMatrix
  
  cat("\n",paste0(temp$data.name, "\n"))
  print(temp$p.value)
  
  
  # # Table_infectionTime_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_incrementTime", 
  #               testType = "t.test")
  
  
  theMatrix = rbind(theMaximumMatrix,
                    theSlopeMatrix,
                    theMidpointMatrix,
                    theStartpointMatrix,
                    theInfectiontimeMatrix)
  
  theMatrix[theMatrix == "rn"] <- ""
  theMatrix[theMatrix == "Test parameter"] <- ""
  
  return(theMatrix)
}
###*****************************



###*****************************
# Convert tables to nice outputs function
###*****************************



###*****************************
# Table 03
df <- samp_028_029_030
matrix_028_029_030 = table_p_value(df = df,
                                   classifier = "dataSet_bio")
write.csv(x = matrix_028_029_030, file = "table_03.csv", quote = F, row.names = F, col.names = F)
###*****************************


###*****************************
# Table 04
df <- samp_035_036_037
matrix_035_036_037 = table_p_value(df = df,
                                   classifier = "dataSet_bio")
write.csv(x = matrix_035_036_037, file = "table_04.csv", quote = F, row.names = F, col.names = F)
###*****************************


###*****************************
# Table 05
samp_035_036_037 %>%
  dplyr::filter(dataSet == "SCV035") -> df

matrix_035 = table_p_value(df = df,
                           classifier = "decision_bio")

write.csv(x = matrix_035, file = "table_05.csv", quote = F, row.names = F, col.names = F)
###*****************************


###*****************************
# Table 06
df <- samp_038_039_040_041
matrix_038_039_040_041 = table_p_value(df = df,
                                       classifier = "dataSet_bio")
write.csv(x = matrix_038_039_040_041, file = "table_06.csv", quote = F, row.names = F, col.names = F)
###*****************************


###*****************************
# Table 07
df <- samp_031_032
matrix_031_032 = table_p_value(df = df,
                               classifier = "dataSet_bio")
write.csv(x = matrix_031_032, file = "table_07.csv", quote = F, row.names = F, col.names = F)
###*****************************


###*****************************
# Table 08
df <- samp_029_031_035
matrix_029_031_035 = table_p_value(df = df,
                                   classifier = "dataSet_bio")
write.csv(x = matrix_029_031_035, file = "table_08.csv", quote = F, row.names = F, col.names = F)
###*****************************