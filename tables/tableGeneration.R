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
  dplyr::summarise(maximum_mean = mean(COMB_maximum_y),
                   maximum_sd = sd(COMB_maximum_y),
                   midPoint_mean = mean(COMB_midPoint1_x),
                   midPoint_sd = sd(COMB_midPoint1_x),
                   slope_mean = mean(COMB_slope1),
                   slope_sd = sd(COMB_slope1),
                   infectionTime_mean = mean(COMB_incrementTime),
                   infectionTime_sd = sd(COMB_incrementTime),
                   startPoint_mean = mean(COMB_startPoint_x),
                   startPoint_sd = sd(COMB_startPoint_x)) -> table_02

###*****************************


###*****************************
# P value Table Function
table_p_value <- function(df)
{
  df %>%
    dplyr::group_by(dataSet_bio) %>%
    dplyr::summarise(numSample = n())-> df_summary
  
  print.data.frame(df_summary)
  
  # Table_maximum_new
  temp = pairwise.t.test(df$COMB_maximum_y, 
                         df$dataSet_bio, 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Maximum"
  print(temp)
  
  # # Table_maximum_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_maximum_y", 
  #               testType = "t.test")
  
  # Table_slope_new
  temp = pairwise.t.test(df$COMB_slope1, 
                         df$dataSet_bio, 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Slope"
  print(temp$data.name)
  print(temp$p.value)
  
  # # Table_slope_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_slope1", 
  #               testType = "t.test")
  
  # Table_midPoint_new
  temp = pairwise.t.test(df$COMB_midPoint1_x, 
                         df$dataSet_bio, 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Midpoint"
  print(temp$data.name)
  print(temp$p.value)
  
  # # Table_midPoint_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_midPoint1_x", 
  #               testType = "t.test")
  
  # Table_startPoint_new
  temp = pairwise.t.test(df$COMB_startPoint_x, 
                         df$dataSet_bio, 
                         p.adjust.method = "fdr", 
                         pool.sd = TRUE)
  temp$data.name <- "Start point"
  print(temp$data.name)
  print(temp$p.value)
  
  # # Table_startPoint_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_startPoint_x", 
  #               testType = "t.test")
  
  # Table_incrementTime_new
  temp = pairwise.t.test(df$COMB_incrementTime, 
                  df$dataSet_bio, 
                  p.adjust.method = "fdr", 
                  pool.sd = TRUE)
  temp$data.name <- "Infection time"
  
  print(temp$data.name)
  print(temp$p.value)
  
  
  # # Table_infectionTime_old
  # tableFunction(dataTable = df, 
  #               seperationVariable = "dataSet_bio", 
  #               condition = "COMB_incrementTime", 
  #               testType = "t.test")
}
###*****************************




###*****************************
# Table 03
df <- samp_028_029_030
table_p_value(df)
###*****************************


###*****************************
# Table 04
df <- samp_035_036_037
table_p_value(df)
###*****************************


###*****************************
# Table 05
# df <- samp_035_036_037
# table_p_value(df)
###*****************************


###*****************************
# Table 06
df <- samp_038_039_040_041
table_p_value(df)
###*****************************


###*****************************
# Table 07
df <- samp_031_032
table_p_value(df)
###*****************************


###*****************************
# Table 08
df <- samp_029_031_035
table_p_value(df)
###*****************************