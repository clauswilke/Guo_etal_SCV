# TABLE FUNCTION 

tableFunctionWODiag <-function(dataTable, seperationVariable, condition, testType){
  
  lengthFunction <-function(vec){u=length(vec); return(u)}
  dots <- list(interp(~lengthFunction(var), var = as.name(condition)))
  
  sig_level_val=0.05
  power_val=0.95
  
  dataTable%>%
    dplyr::group_by_("dataSet",seperationVariable) %>%
    dplyr::summarize_(.dots = setNames(dots, c("usedData")))-> dataSize
  
  fig_table=data.frame(matrix(vector(), nrow(dataSize)+1, nrow(dataSize)+1))
  colnames(fig_table) <- c(as.vector(dataSize[[seperationVariable]]),"numSamples")
  row.names(fig_table) <- c(as.vector(dataSize[[seperationVariable]]),"numSamples")
  fig_table[ncol(fig_table),1:(nrow(dataSize))] <- as.vector(dataSize$usedData)
  fig_table[1:(nrow(dataSize)),nrow(fig_table)] <- as.vector(dataSize$usedData)
  
  
  for (counter01 in 1:(nrow(dataSize))){
    for (counter02 in 1:(nrow(dataSize))){
      
      
      var_inp=seperationVariable
      var_out="tempCol"
      mutateFunction <-function(name){u=name; return(u)}
      dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
      dataTable %>%  mutate_(.dots = setNames(dots, c(var_out))) -> dataTable
      
      temp01=as.vector(dataSize[[seperationVariable]])[counter01]
      dataTable %>% dplyr::filter(tempCol %in% temp01)->firstDataSet
      firstDataSet=firstDataSet[[condition]]
      
      temp02=as.vector(dataSize[[seperationVariable]])[counter02]
      dataTable %>% dplyr::filter(tempCol %in% temp02) -> secondDataSet
      secondDataSet=secondDataSet[[condition]]
      
      
      if(length(unique(firstDataSet))==1 & length(unique(secondDataSet))==1)
      {
        firstDataSet=firstDataSet+runif(length(firstDataSet), 0, 1)*mean(firstDataSet)/100000
        secondDataSet=secondDataSet+runif(length(secondDataSet), 0, 1)*mean(secondDataSet)/100000
      }
      
      if (testType=="t.test")
      {
        test_Value <- "NA"
        try(test_Value <- t.test(firstDataSet, secondDataSet)$p.value, silent = TRUE)
      }
      if (testType=="wilcox.test")
      {
        test_Value <- "NA"
        try(test_Value <- wilcox.test(firstDataSet, secondDataSet)$p.value, silent = TRUE)
      }
      if (testType=="difMean")
      {
        test_Value=abs(mean(firstDataSet)-mean(secondDataSet))
        #if(is.na(test_Value)){browser()}
      }
      if (testType=="t.test.effectSize")
      {
        d_Value <- "NA"
        try(d_Value <- pwr.t2n.test(n1 = length(firstDataSet),
                                n2 = length(secondDataSet),
                                d = NULL,
                                sig.level = sig_level_val,
                                power = power_val)$d, silent = TRUE)
        test_Value <- "NA"
        try(test_Value <- d_Value * 
          sqrt((sd(firstDataSet,na.rm = T)^2 + sd(secondDataSet,na.rm = T)^2)/2),  silent = TRUE) #need to check this part
      } 
      #test_Value=d_Value*1}
      if(testType=="t.test.effectSize")
      {
        fig_table[counter01,counter02] <- "NA"
        try(fig_table[counter01,counter02] <- sprintf("%.3e",test_Value), silent = TRUE)
      }
      if(testType=="t.test" | testType=="wilcox.test" | testType=="difMean")
      {
        fig_table[counter01,counter02] <- "NA"
        try(fig_table[counter01,counter02] <- sprintf("%.5f",test_Value), silent = TRUE)
      }
    }
  }
  fig_table_notCorrected=fig_table
  
  if (testType %in% c("t.test","wilcox.test")){
    dataPortion=fig_table[1:nrow(fig_table)-1,1:ncol(fig_table)-1]
    getUpper=dataPortion[upper.tri(dataPortion, diag=TRUE)]
    getUpperAdj=p.adjust(as.numeric(getUpper),method = "fdr")
    #getUpperAdj=as.numeric(getUpper)
    
    adj.matrix=matrix(nrow = (nrow(fig_table)-1), ncol = (ncol(fig_table)-1))
    adj.matrix2=matrix(nrow = (nrow(fig_table)-1), ncol = (ncol(fig_table)-1))
    adj.matrix[upper.tri(adj.matrix, diag=TRUE)] <- getUpperAdj
    adj.matrix2=adj.matrix
    diag(adj.matrix2) = 1 # new line
    adj.matrix2[lower.tri(adj.matrix2)] <- t(adj.matrix)[lower.tri(adj.matrix)]
    
    fig_table[1:(nrow(fig_table)-1),1:(ncol(fig_table)-1)]<-sprintf("%.5f",adj.matrix2)
  }
  
  
  if(setequal( as.vector(unique(dataTable[["decision_bio"]])) , "infection" )){filterTypeLength=0}
  if(setequal( as.vector(unique(dataTable[["decision_bio"]])) , "infection&lysis" )){filterTypeLength=1}
  if(setequal( as.vector(unique(dataTable[["decision_bio"]])) , c("infection" ,"infection&lysis") )){filterTypeLength=2}
  
  # filterTypeLength=length(as.vector(unique(dataTable[["decision_bio"]])))
  # # if this length is 1 it is infLys
  # # if this length is 2 it is inf + infLys
  
  if (condition=="COMB_maximum_y" & filterTypeLength==0){
    additionalText02="Maximum Distribution (inf)"
  }
  if (condition=="COMB_maximum_y" & filterTypeLength==1){
    additionalText02="Maximum Distribution (infLys)"
  }
  if (condition=="COMB_maximum_y" & filterTypeLength==2){
    additionalText02="Maximum Distribution (inf + infLys)"
  }
  if (condition=="maximum_Estimate_norm" & filterTypeLength==0){
    additionalText02="Maximum Distribution Median Normalized (inf)"
  }
  if (condition=="maximum_Estimate_norm" & filterTypeLength==1){
    additionalText02="Maximum Distribution Median Normalized (infLys)"
  }
  if (condition=="maximum_Estimate_norm" & filterTypeLength==2){
    additionalText02="Maximum Distribution Median Normalized (inf + infLys)"
  }
  if (condition=="COMB_slope1" & filterTypeLength==0){
    additionalText02="Slope Distribution (inf)"
  }
  if (condition=="COMB_slope1" & filterTypeLength==1){
    additionalText02="Slope Distribution (infLys)"
  }
  if (condition=="COMB_slope1" & filterTypeLength==2){
    additionalText02="Slope Distribution (inf + infLys)"
  }
  if (condition=="slope_Estimate_norm" & filterTypeLength==0){
    additionalText02="Slope Distribution Median Normalized (inf)"
  }
  if (condition=="slope_Estimate_norm" & filterTypeLength==1){
    additionalText02="Slope Distribution Median Normalized (infLys)"
  }
  if (condition=="slope_Estimate_norm" & filterTypeLength==2){
    additionalText02="Slope Distribution Median Normalized (inf + infLys)"
  }
  if (condition=="COMB_midPoint1_x" & filterTypeLength==0){
    additionalText02="Midpoint Distribution (inf)"
  }
  if (condition=="COMB_midPoint1_x" & filterTypeLength==1){
    additionalText02="Midpoint Distribution (infLys)"
  }
  if (condition=="COMB_midPoint1_x" & filterTypeLength==2){
    additionalText02="Midpoint Distribution (inf + infLys)"
  }
  if (condition=="midPoint_Estimate_norm" & filterTypeLength==0){
    additionalText02="Midpoint Distribution Median Normalized (inf)"
  }
  if (condition=="midPoint_Estimate_norm" & filterTypeLength==1){
    additionalText02="Midpoint Distribution Median Normalized (infLys)"
  }
  if (condition=="midPoint_Estimate_norm" & filterTypeLength==2){
    additionalText02="Midpoint Distribution Median Normalized (inf + infLys)"
  }
  if (condition=="COMB_incrementTime" & filterTypeLength==0){
    additionalText02="infection time Distribution (inf)"
  }
  if (condition=="COMB_incrementTime" & filterTypeLength==1){
    additionalText02="infection time Distribution (infLys)"
  }
  if (condition=="COMB_incrementTime" & filterTypeLength==2){
    additionalText02="infection time Distribution (inf + infLys)"
  }
  if (condition=="infectiontime_Estimate_norm" & filterTypeLength==0){
    additionalText02="infection time Distribution Median Normalized (inf)"
  }
  if (condition=="infectiontime_Estimate_norm" & filterTypeLength==1){
    additionalText02="infection time Distribution Median Normalized (infLys)"
  }
  if (condition=="infectiontime_Estimate_norm" & filterTypeLength==2){
    additionalText02="infection time Distribution Median Normalized (inf + infLys)"
  }
  if (condition=="COMB_startPoint_x" & filterTypeLength==0){
    additionalText02="start point Distribution (inf)"
  }
  if (condition=="COMB_startPoint_x" & filterTypeLength==1){
    additionalText02="start point Distribution (infLys)"
  }
  if (condition=="COMB_startPoint_x" & filterTypeLength==2){
    additionalText02="Start point (inf + infLys)"
  }  
  
  if (testType=="t.test"){
    additionalText01="ttest. adjusted p values"
    #additionalText01="ttest"
  }
  if (testType=="wilcox.test"){
    additionalText01="Mann-Whitney U Test, adjusted p values"
    #additionalText01="Mann-Whitney U Test (wilcox.test)"
  }
  if (testType=="difMean"){
    additionalText01="Diffrence between means"
  }
  if (testType=="t.test.effectSize"){
    additionalText01=paste0("effect size needed for Power = ",power_val," and sig.level = ",sig_level_val)
  }
  
  
  fig_table[is.na(fig_table)] <- c("")
  returnObj=list(additionalText01,additionalText02,fig_table)
  
  return(returnObj)
}