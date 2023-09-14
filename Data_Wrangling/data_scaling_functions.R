library(mixOmics) #masks select!!!
library(dplyr)
library(ggpubr)
library(qs)
library(janitor)
library(tidyverse)
library(FSA, exclude = "select")
library(paletteer)
library(rstatix , exclude = "select")
#library(ggbeeswarm)
library(ggprism)
library(openxlsx)
#library(coin)
library(recipes)
library(textclean)
library(ggplot2)
library(ggrepel)
#pivot wider
library("pheatmap")
library("paletteer")
library(docstring)
library("gt")


long_to_wide_df <- function(long_df,  Col_to_spread, Col_with_values ) {
  
  
  #' Formats long to wide style data columns. Only works on one pair of "Variable-Data" columns at a time!
  #'
  #' @param long_df the input long df
  #' @param Col_to_spread the col which gets divided into multible new cols 
  #' @param Col_with_values the col which contains the relevant values
  #'
  #' @return 
  #' @export
  #'
  #' @examples long_to_wide_df(data, "Variable", "Measurement")
  
  
  #only include colums which together uniquely identify samples per row in combination!!
  
  wide_df <- long_df  %>%
    pivot_wider(names_from =  Col_to_spread, values_from =  Col_with_values)
  print("longtowide")
  return(wide_df)
}




remove_NA_Cols <- function (wide_df) {
  
  #' Removes columns which only contain missing values
  #'
  #' @param wide_df a df in wide format
  #' @return 
  #' @export
  #'
  #' @examples
  
  not_all_na <- function(x) {!all(is.na(x))}
  wide_df <- wide_df  %>% select_if(not_all_na)
  return (wide_df)
}


OptimizeVariableNames <- function(df, dataformat = "long", datastructureinput, list_of_Names = NULL) {
  
  #' Removes non ASCII Characters from Variable Col or names of cols or replaces names of cols with given vector
  #'
  #' @param df 
  #' @param dataformat select "long" for cleaning 1 column or "wide" for all names of the df 
  #' @param datastructureinput the list of Column names which define the datastructure
 
 
  #' @param list_of_Names a vector with the same length and order as the current col names of a long df
  #'
  #' @return
  #' @export
  #'
  #' @examples 
  

  Variable = datastructureinput$Variable
  
  if (is.null(list_of_Names) == FALSE) {
    if (len(list_of_Names) != len(names(df))) {
      print("vector has a different length than there are col names")
    }
      
    names(df) <- list_of_Names
    return (df)
  }
  
if (dataformat == "long") {
  df[Variable] <- replace_non_ascii(df[[Variable]])
  df[Variable] <-gsub("\\s", "",df[[Variable]])
  df[Variable] <-gsub(",", ".",df[[Variable]])
  return (df)
  }

  if (dataformat == "wide") {
    names(df) <- replace_non_ascii(names(df))
    names(df) <-gsub("\\s", "",names(df))
    names(df)<-gsub(",", ".",names(df))
    return (df)
  }
  
}




RemoveZeroVariance <- function (long_df, datastructureinput = NULL) {
  
  
  #' remove 0 variance variables, can handle 2 variable-columns, adds label, removes variables which are only in 1 cohort. 
  #'
  #' @param long_df, unfisnished doc
  #' @param datastructureinput 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  
  cohort = datastructureinput$cohort
  Variable = datastructureinput$Variable
  Variable2 = datastructureinput$Variable2
  measured = datastructureinput$measured
  
  if (is.null(Variable2) == TRUE) {
     Variable2 = ""
   }
  

    prepo_data_long <- long_df
    
    
    prepo_data_long$label<- paste0(prepo_data_long[[Variable]],prepo_data_long[[Variable2]])
    
    prepo_data_long[cohort]<- as.factor(prepo_data_long[[cohort]])
    
# removes metadata columns..... needs fixing
    
    
    #only includes features which appear atleast twice
    prepo_data_long <- prepo_data_long %>% dplyr::group_by(label,.data[[cohort]] ) %>%  dplyr::filter(n()>2) %>% ungroup()  
    
    #checks that every measured value appears in all cohorts################
    x <-list()
    for (i in as.character((unique(prepo_data_long[[cohort]])))){
      print(i)
      temp <-unique(prepo_data_long%>%dplyr::filter(.data[[cohort]] == i )%>%dplyr::select(label)%>% pull())
      x[[i]] <- temp
    }
    cols_intersection <- c(Reduce(intersect,x))
    
    #only includes features which are in both groups 
    prepo_data_long <- prepo_data_long %>%dplyr::filter(label %in% c(cols_intersection))
    
    #eliminates 0 standart deviation
    temp2 <- prepo_data_long %>% dplyr::group_by(label) %>% dplyr::summarise(sd_var1 = sd(.data[[measured]], na.rm=TRUE))  %>% dplyr::filter(sd_var1 !=0)
    variables_without_Zero_std <-c(temp2$label)
    df_g_backup <- prepo_data_long %>% dplyr::filter(label %in% variables_without_Zero_std)
    df_g <- df_g_backup
    
    print("0 variance removed")   
return(df_g)
}



identify_outliers_new <- function(data, ..., variable = NULL){
  
  is.outlier_new  <- NULL
  if(is_grouped_df(data)){
    results <- data %>%
      doo(identify_outliers_new , ..., variable = variable)
    if(nrow(results) == 0) results <- as.data.frame(results)
    return(results)
  }
  
  if(!inherits(data, "data.frame"))
    stop("data should be a data frame")
  variable <- data %>% get_selected_vars(..., vars = variable)
  n.vars <- length(variable)
  if(n.vars > 1)
    stop("Specify only one variable")
  values <- data %>% pull(!!variable)
  results <- data %>%
    mutate(
      is.outlier_new  = is_outlier_new (values),
      is.extreme_new = is_extreme_new (values)
    ) %>%
    filter(is.outlier_new  == TRUE)
  if(nrow(results) == 0) results <- as.data.frame(results)
  results
}


#' @describeIn outliers detect outliers in a numeric vector. Returns logical vector.
#' @export
is_outlier_new  <- function(x, coef = 1.5){
  res <- x
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  .IQR <- IQR(x, na.rm = TRUE)
  upper.limit <- Q3 + (coef*.IQR)
  lower.limit <- Q1 - (coef*.IQR)
  outlier <- ifelse(x < lower.limit | x > upper.limit, TRUE, FALSE )
  outlier
}

#' @describeIn outliers detect extreme points in a numeric vector. An alias of
#'   \code{is_outlier()}, where coef = 3. Returns logical vector.
#' @export
is_extreme_new  <- function(x){
  is_outlier_new (x, coef = 30)
}



get_selected_vars <- function(x, ..., vars = NULL){
  
  if(is_grouped_df(x))
    x <- x %>% dplyr::ungroup()
  dot.vars <- rlang::quos(...)
  
  if(length(vars) > 0){
    return(vars)
  }
  if (length(dot.vars) == 0) selected <- colnames(x)
  else selected <- tidyselect::vars_select(names(x), !!! dot.vars)
  selected %>% as.character()
}



# Identify outliers by groups
#x <- df_g %>%
#  group_by(Variable,Celltype,outcome) %>%
#  identify_outliers_new("Measurement")
# }

######################
# replace missing values by 1/4 of the minimum per variable

# Limit of detection (1/4 of min for each var), value sets divisor
.replace.by.lod <- function(x, value){
  lod <- min(x[x>0], na.rm=T)/value;
  x[x==0|is.na(x)] <- lod;
  return(x);
}

# replace missing values by 1/4 of the minimum per variable
ReplaceMissingByLoD <- function(wide_df, value = 4){
  
  
  # Limit of detection (1/4 of min for each var), value sets divisor
  .replace.by.lod <- function(x){
    lod <- min(x[x>0], na.rm=T)/value;
    x[x==0|is.na(x)] <- lod;
    return(x);
  }
  
  
  metadata <- wide_df[metadata_columns]
  int.mat <-  wide_df %>% select (-any_of(c(metadata_columns)))
  #int.mat <- as.matrix(matrix);
  
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  int.mat <- apply(int.mat, 2, .replace.by.lod);
  rownames(int.mat) <- rowNms;
  colnames(int.mat) <- colNms;
  result <- bind_cols(metadata, int.mat)
  
  return (result);
}

#ReplaceMissingByLoD(input_heatmap,4)






########### create df for calculations #########
#join matrix and metadata by 1 or more shared ID Cols and create unique label col

metadata_columns <- (c("Alias" ,"Celltype","Stimulation","outcome","cohort"))

Join.Metadata <- function(matrix,metadata, ID1, ID2 = NULL, ID3 = NULL){
  
  result <-merge(matrix, metadata, 
                 by = c(ID1,ID2,ID3))
  #result$sampleID <- paste0(matrix[,ID1],"",matrix[,ID2],"",matrix[,ID3])
  
  
  return(result);
}

#Join.Metadata(matrix3,metadata3, ID1 = "Alias", ID2 ="Celltype")




#from wide to long df, metadata gets put to own columns
create.long.Df <- function(x, metadata_columns){
  long_df<- tidyr::gather( x ,"Variable","Measurement",- c(metadata_columns))
  return(long_df);
}

#create.long.Df(input_heatmap)

#filter df by metadata
`%notin%` <- Negate(`%in%`)
#filtered_df <- raw_df %>% filter (Celltype == "CD14") %>% filter (Alias %notin%  c("PHOA", "H8LB"))
  


#log2transform data 
log2transform <- function(x){
  metadata <- x[metadata_columns]
  matrix <-  x %>% dplyr::select (-any_of(c(metadata_columns)))
  
  matrix <- matrix +1 
  matrix <- log2((abs(matrix))^(sign(matrix)))
  result <- bind_cols(metadata, matrix)
  
  return(result);
}

#normalize
Normalizer <- function(x, scaleNorm, ref.smpl = NULL, ref = NULL){
  metadata <- x[metadata_columns]
  matrix <-  x %>% dplyr::select (-any_of(c(metadata_columns)))
  
  #paretonorm
  ParetoNorm<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T));
  }
  
  
  # normalize to zero mean and unit variance
  AutoNorm<-function(x){
    (x - mean(x))/sd(x, na.rm=T);
  }
  
  # normalize to zero mean but variance/SE
  RangeNorm<-function(x){
    if(max(x) == min(x)){
      x;
    }else{
      (x - mean(x))/(max(x)-min(x));
    }
  }
  
  
  # normalize to zero mean but variance/SE
  MeanCenter<-function(x){
    x - mean(x);
  }
  
  
  # normalize by median
  MedianNorm<-function(x){
    x/median(x, na.rm=T);
  }
  
  # normalize by a reference sample (probability quotient normalization)
  # ref should be the name of the reference sample
  ProbNorm<-function(x, ref.smpl){
    x/median(as.numeric(x/ref.smpl), na.rm=T)
  }
  
  # normalize by a reference reference (i.e. creatinine)
  # ref should be the name of the cmpd
  CompNorm<-function(x, ref){
    1000*x/x[ref];
  }
  
  # perform quantile normalization on the raw data (can be log transformed later by user)
  # https://stat.ethz.ch/pipermail/bioconductor/2005-April/008348.html
  QuantileNormalize <- function(data){
    return(t(preprocessCore::normalize.quantiles(t(data), copy=FALSE)));
  }
  
  #apply function
  if(scaleNorm=='MeanCenter'){
    data<-apply(matrix, 2, MeanCenter);
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'){
    data<-apply(matrix, 2, AutoNorm);
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(matrix, 2, ParetoNorm);
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(matrix, 2, RangeNorm);
    scalenm<-"Range Scaling";
  }else if(scaleNorm=='MedianNorm'){
    data<-apply(matrix, 2, MedianNorm);
    scalenm<-"MedianNorm";
  }else if(scaleNorm=='ProbNorm'){
    data<-apply(matrix, 2, ProbNorm);
    scalenm<-"ProbNorm";
  }else if(scaleNorm=='CompNorm'){
    data<-apply(matrix, 2, CompNorm);
    scalenm<-"CompNorm";
  }else if(scaleNorm=='QuantileNormalize'){
    data<-apply(matrix, 2, QuantileNormalize);
    scalenm<-"QuantileNormalize";
  }else{
    scalenm<-"N/A";
  }
  result <- bind_cols(metadata, data)
  comment(result) <- scalenm
  return(result);
}

#CompNorm <- "QuantileNormalize"
#Normalizer(input_heatmap, CompNorm)




#remove columns with more than 70% 0 or NAN values
RemoveMissingData <- function(wide_df, limit = 0.5){
  
  metadata <- wide_df[metadata_columns]
  matrix <-  wide_df %>% select (-any_of(c(metadata_columns)))
  
  zero_pct <- matrix %>% summarize(across(everything(), ~ mean(. == 0|is.na(.))))
  # select columns with less than or equal to 70% 0s
  tb_filtered <- matrix %>% select(which(zero_pct <= limit))
  selected_cols<- colnames(tb_filtered)
  matrix <- matrix %>% select( c(colnames(tb_filtered)))
  result <- bind_cols(metadata, matrix)
  return(result);
}


#RemoveMissingData(input_heatmap, limit = 0.7)



#converts data cols to numeric
ToNumeric <- function(wide_df, metadata_columns){
  metadata <- wide_df[metadata_columns]
  matrix <-  wide_df %>% select (-any_of(c(metadata_columns)))
  
  StrScrub <- function(x) {
    gsub("[^0-9.-]","", x)
  }
  
  commaScrub <- function(x) {
    gsub(",",".", x)
  }
  
  matrix <-  mutate_all(matrix, funs(commaScrub)) %>%
    mutate_all( funs(StrScrub))
    
  matrix <- matrix %>% mutate_if(is.character, as.numeric)
  result <- bind_cols(metadata, matrix)
  return(result);
}

#ToNumeric(input_heatmap)
  
  
RemoveColswithMissingData <- function(wide_df){
  
  metadata <- wide_df[metadata_columns]
  matrix <-  wide_df %>% select (-any_of(c(metadata_columns)))
  
  matrix <- matrix %>% select(-c(names(which(apply(is.na(.), 2, any)))))
  
 
  result <- bind_cols(metadata, matrix)
  return(result);
}




#####

CalculateFoldChange <- function (long_df, groupingCount, datastructureinput, used_stat = "mean") {
  cohort = datastructureinput$cohort
  Variable = datastructureinput$Variable
  Variable2 = datastructureinput$Variable2
  measured = datastructureinput$measured
  firstcohort =datastructureinput$firstcohort
  secondcohort = datastructureinput$secondcohort
  
  
  
  if (groupingCount == 1) {
    
    if (used_stat == "mean") {
    mean_df <- long_df %>% dplyr::group_by(.data[[cohort]], .data[[Variable]] )%>%    dplyr::summarize(Mean= mean(.data[[measured]], na.rm=TRUE))
    }
    if (used_stat == "median") {
    mean_df <- long_df %>% dplyr::group_by(.data[[cohort]], .data[[Variable]] )%>%    dplyr::summarize(Mean= median(.data[[measured]], na.rm=TRUE))
    }
    
    #splitting by cohort and joining for calculations
    high_mean <- mean_df %>% dplyr::filter( .data[[cohort]] == firstcohort)%>% dplyr::rename( Firstcohort_mean = Mean)
    low_mean <- mean_df %>% dplyr::filter( .data[[cohort]] == secondcohort)%>% dplyr::rename( Secondcohort_mean = Mean)
    highandlow <- high_mean %>% left_join( low_mean, by = c(Variable ))
    #keeps negatives intact log(a/b) == log(a)-log(b))
    log2foldchanges <-  log2((1+abs(highandlow$Secondcohort_mean ))^(sign(highandlow$Secondcohort_mean))) - log2((1+abs(highandlow$Firstcohort_mean ))^(sign(highandlow$Firstcohort_mean)))
    log2foldchanges<- cbind(log2foldchanges)
    log2foldchanges <- as.data.frame(log2foldchanges)
    #log2foldchanges <- tibble::rownames_to_column(log2foldchanges)
    log2foldchanges[c(Variable)]<- unique(highandlow[c(Variable)])  #vector of unique measurements in the correct order
    log2foldchanges <- log2foldchanges %>% left_join( highandlow[,c(Variable,"Firstcohort_mean","Secondcohort_mean")],  by = c(Variable ))
    
    return (log2foldchanges)
    
  }
  if (groupingCount == 2) {
 
    if (used_stat == "mean") {
  mean_df <- long_df %>% dplyr::group_by(.data[[cohort]], .data[[Variable]], .data[[Variable2]]  )%>%    dplyr::summarize(Mean= mean(.data[[measured]], na.rm=TRUE))
    }
    
    if (used_stat == "median") {
      mean_df <- long_df %>% dplyr::group_by(.data[[cohort]], .data[[Variable]], .data[[Variable2]]  )%>%    dplyr::summarize(Mean= median(.data[[measured]], na.rm=TRUE))
    }
    
    
   #splitting by cohort and joining for calculations
  high_mean <- mean_df %>% dplyr::filter( .data[[cohort]] == firstcohort)%>% dplyr::rename( Firstcohort_mean = Mean)
  low_mean <- mean_df %>% dplyr::filter( .data[[cohort]] == secondcohort)%>% dplyr::rename( Secondcohort_mean = Mean)
  highandlow <- high_mean %>% left_join( low_mean, by = c(Variable,Variable2 ))
  #keeps negatives intact log(a/b) == log(a)-log(b))
  log2foldchanges <-  log2((1+abs(highandlow$Secondcohort_mean ))^(sign(highandlow$Secondcohort_mean))) - log2((1+abs(highandlow$Firstcohort_mean ))^(sign(highandlow$Firstcohort_mean)))
  log2foldchanges<- cbind(log2foldchanges)
  log2foldchanges <- as.data.frame(log2foldchanges)
  #log2foldchanges <- tibble::rownames_to_column(log2foldchanges)
  log2foldchanges[c(Variable,Variable2)]<- unique(highandlow[c(Variable,Variable2)]) #vector of unique measurements in the correct order
  log2foldchanges <- log2foldchanges %>% left_join( highandlow[,c(Variable,Variable2,"Firstcohort_mean","Secondcohort_mean")],  by = c(Variable,Variable2 ))
  
  return (log2foldchanges)
  }
}




kruskal_testing <- function (long_df, groupingCount, posthoc, datastructureinput ) {
  cohort = datastructureinput$cohort
  Variable = datastructureinput$Variable
  Variable2 = datastructureinput$Variable2
  measured = datastructureinput$measured
  
  
  
  
  #first global testing###################
  
  if (groupingCount == 1 ) {
  long_df$rowID <- paste0(long_df[[Variable]])
  stat.test <-long_df %>%
  group_by(!!sym(Variable)) %>%
  rstatix::kruskal_test(as.formula(paste(measured, '~', cohort))) %>%
  add_significance()
  stat.test$rowID <- paste0(stat.test[[Variable]])
  stat.test$pothoc <- "not_possible"
  #filters out nonsignificant results
  filtered_P <-  dplyr::filter(stat.test, p < 0.06)
  colselect <-  filtered_P$rowID
  significant_MCs<-dplyr::filter(long_df,  rowID %in% colselect)
  
  if ((posthoc == TRUE) && (nrow(significant_MCs) > 0))  {
    print("significant differences in Kruskal testing detected!")
    stat.test.adj <-significant_MCs %>%
      group_by(!!sym(Variable)) %>%
      rstatix::dunn_test(as.formula(paste(measured, '~', cohort))) %>%
      adjust_pvalue(method = 'BH')%>%
      add_significance()
    stat.test.adj$rowID <- paste0(stat.test.adj[[Variable]])
    stat.test.adj$posthoc <- "Dunn"
  }
  else{
    return(stat.test)
  }
  
  
  }
  
  if (groupingCount == 2 ) {
    long_df$rowID<- paste0(long_df[[Variable]],".",long_df[[Variable2]])
    stat.test <-long_df %>%
      group_by(!!sym(Variable), !!sym(Variable2)) %>%
      rstatix::kruskal_test(as.formula(paste(measured, '~', cohort))) %>%
      add_significance() 
    stat.test$rowID <- paste0(stat.test[[Variable]],".",stat.test[[Variable2]])
    stat.test$posthoc <- "not_possible"
    
    
    
    #filters out nonsignificant results
    filtered_P <-  dplyr::filter(stat.test, p < 0.06)
    colselect <-  filtered_P$rowID
    significant_MCs<-dplyr::filter(long_df, rowID %in% colselect)
    
    if ((posthoc == TRUE) && nrow(significant_MCs) > 0)  {
      print("significant differences in Kruskal Testing  detected!")
      stat.test.adj <-significant_MCs %>%
        group_by(!!sym(Variable), !!sym(Variable2)) %>%
        rstatix::dunn_test(as.formula(paste(measured, '~', cohort))) %>%
        adjust_pvalue(method = 'BH')%>%
        add_significance()
      stat.test.adj$rowID <- paste0(stat.test.adj[[Variable]],".",stat.test.adj[[Variable2]])
      stat.test.adj$posthoc <- "Dunn"
    }
    else{
      return(stat.test)
    }
    
  }
  return (stat.test.adj)
}


?t_test

Students_T_Test <- function (long_df, groupingCount, datastructureinput) {
  
  
  
  cohort = datastructureinput$cohort
  Variable = datastructureinput$Variable
  Variable2 = datastructureinput$Variable2
  measured = datastructureinput$measured
  
 if (nrow(long_df)==0) {
   print("no data detected")
   return(NULL)
 }
  
  
  if (groupingCount == 1 ) {
    long_df$rowID <- paste0(long_df[[Variable]])
    stat.test <-long_df %>%
      group_by(!!sym(Variable)) %>%
      rstatix::t_test(as.formula(paste(measured, '~', cohort)),  alternative = "two.sided", detailed = TRUE) %>%
      adjust_pvalue(method = "BH")%>% 
      add_significance()
    stat.test$rowID <- paste0(stat.test[[Variable]])
    
  }
  
  if (groupingCount == 2 ) {
    long_df$rowID<- paste0(long_df[[Variable]],".",long_df[[Variable2]])
    stat.test <-long_df %>%
      group_by(!!sym(Variable), !!sym(Variable2)) %>%
      rstatix::t_test(as.formula(paste(measured, '~', cohort)),  alternative = "two.sided") %>%
      adjust_pvalue(method = "BH")%>% 
      add_significance() 
    stat.test$rowID <- paste0(stat.test[[Variable]],".",stat.test[[Variable2]])
  }

return (stat.test)
}

Sig.For.Plotting <- function (stat.test, scaling_of_plot, filter_of_sig, x, use_posthoc = T, datastructureinput){
  #places significance brackets correctly and filters out nonsign. stats rows. Scaling MUST be identical as plot!
  
  dim_df <- stat.test %>% dplyr::filter(p < filter_of_sig)
  
  if (nrow(dim_df) == 0 ){
    print("no significant events detected")
    return (NULL)
  }
  
  
  
  stat.test.adj <- stat.test %>% add_xy_position(scales = scaling_of_plot, x = x, 
                                                 dodge = 0.8, fun = "max")  

  
  
    if (use_posthoc == T){
      stat.test.adj <- stat.test.adj %>% dplyr::filter(p.adj < filter_of_sig) #y.trans = Log2
      
      dim_df <- stat.test.adj %>% dplyr::filter(p.adj < filter_of_sig)
      
      if (nrow(dim_df) == 0 ){
        print("no significant events in posthoc test detected")
        return (NULL)
      }
      
      
    }
  
  
  
  if (use_posthoc == F){
    stat.test.adj <- stat.test.adj %>% dplyr::filter(p < filter_of_sig) #y.trans = Log2
  }
  cohort = datastructureinput$cohort
  Variable2 = datastructureinput$Variable2
  
  if (is.null(stat.test.adj$cohort) == T){
  stat.test.adj[cohort] <- "A"
 # stat.test.adj[Variable2] <- "A"
  }
  return (stat.test.adj)
  
}


#coloring function 

create_color_mapping <- function(inputDF, datastructure, selection, palette_name, palette_style) {
  if (palette_style == "continuous") {
    colorrange = paletteer_c
  }
  if (palette_style == "dynamic") {
    colorrange = paletteer_dynamic
  }
  if (palette_style == "discrete") {
    colorrange = paletteer_d
  }
  
  
  
  cohort_unique <- unique(inputDF[[datastructure[[selection]]]])
  print(paste("amount of variabels:", length(cohort_unique)))
  colors_selected <- c(colorrange (palette_name, n = length(cohort_unique)))
  color_p <- setNames(colors_selected, c(cohort_unique))
  return(color_p)

  }
  


  
setPlottingFormat <- function(annotation, ...) {
    #' set the structure of your plot. Function creates a named list. 
    #'@param annotation name the list for its function with a str vector
    #' @param cohort defines what is being plotted on x axis, input equals to column names of long df
    #' @param Variable defines the facetting parameter
    #' @param Variable2 defines a subparameter that is being compared on the x axis
    #' @param Measurement defines the column storing the numeric values 
    #'
    #' @return a list object which is fed into plotting functions
    #' @export
    #'
    #' @examples setPlottingFormat(annotation, cohort = "Celltype", Variable = "Variable", Variable2 = "cohort", Measurement = "Measurement")
    
  
      
      args <- list(...)
      attr(args, "annotation") <- annotation
      print(attr(args, "annotation"))
      
      
      return(args)
    }





#plotting functions

#1 boxplot
Boxplots_simple<- function (long_df,datastructure, stat.test = NULL, group.colors,  dotplot = T, P_values = T, posthoc = F,show_n = F,
                            title = "Test", headersize = 10, legendsize = 30, alpha_dots = 0.5, alpha_value = 0.5,x_axis_labelsize = 5,  Facetting_by_Variable = T) { 
  #for nice stats labels, rounded and styled
  if (is.null(stat.test) == FALSE) {
    
    if (posthoc == T) {
      p_value <- "p.adj"
    } else{ p_value <- "p"}
    
  stat.test$p.adj.ital <-round(stat.test[[p_value]], 5)
  print("using stats data")
  
  }
  
  cohort <- datastructure$cohort
  
  #data in
  data <-long_df
  data[[cohort]] <- as.factor(data[[cohort]])
 
  
  #plotting
  p<-ggplot(data, aes(x = .data[[datastructure$cohort]], y = .data[[datastructure$measured]],  fill =.data[[datastructure$Variable2]] ))+#
    ggtitle(title) +
    #p<-ggplot(significant_MCs, aes(x = cohort, y = MFIs)) +
    #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
    
    #annotation_logticks(long = unit(0.3, "cm"),sides = "l") +
    #coord_fixed() +
    stat_boxplot( aes(colour = .data[[datastructure$cohort]], fill =.data[[datastructure$Variable2]])) +
    
    
    
    theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
    theme(panel.grid.minor = element_blank())+
    scale_color_manual(values=group.colors)+
    scale_fill_manual( values= alpha(c(group.colors),alpha_value))+
    #facet_wrap(as.formula(paste('~', Variable)), scales = "fixed", ncol = 4)+
    
   
    
    #theme_minimal()+
    theme(strip.text.x=element_text(size= headersize),  legend.text=element_text(size=legendsize)) +
    theme(axis.title.x=element_blank(),
         axis.text.x = element_text(size=x_axis_labelsize)
          # axis.text.x=element_blank(),
          #axis.ticks.x=element_blank()
    ) 

  if (Facetting_by_Variable == T) {
   p <- p + facet_wrap(as.formula(paste('~', datastructure$Variable)), scales = "free", ncol = 4)
  }
  
  if (dotplot == T) {
  #  p = p + geom_point(position=position_jitterdodge(),alpha = alpha_dots, ,show.legend = F)
     # Boxplot width can be adjusted
      
     p = p+  geom_jitter(show.legend=FALSE, aes(color=.data[[datastructure$Variable2]]), 
                  position=position_dodge(width=0.8), 
                  size=2,
                  alpha = alpha_dots
      )
    

    
  }
  
  
  if (P_values == T) {
    print("mapping stats data")
    p = p + add_pvalue(stat.test, 
                       xmin = "xmin", 
                       xmax = "xmax",
                       y.position = "y.position",
                       label = "p = {p.adj.ital}",
                       tip.length = 0)
  }
  
  
  if(show_n == T) {
   
     #calculate stats
    sumstats <- data %>% 
      
      group_by(.data[[datastructure$Variable]], 
               #.data[[datastructure$Variable2]],
               .data[[datastructure$cohort]]) %>% 
      get_summary_stats(Measurement, show = c("iqr","mean","median", "n", "q1", "q3"))
    sumstats$upper_limit <- sumstats$q3+ 1.5 * sumstats$iqr
    sumstats$lower_limit <- sumstats$q1- 1.5 * sumstats$iqr
    
    
    
    p <- p + geom_text(data = sumstats, aes(x = .data[[datastructure$cohort]], y =upper_limit+1, color = .data[[datastructure$Variable2]], label = paste0("n:", n)), size= 2,check_overlap = F,
                       vjust = 0,  position=position_dodge(width=0.8)) 
    
    
    
    
  }
  
  return (p)
}


Heatmap.simple <- function(input,metadata_columns, breaks = seq(-2, 2, length.out = 100), 
                           sampleID = "sampleID", cohort = "cohort", Variable2, my_colour, title_heatmap) 
  {
  
  #input: wide df
  #metadata_columns: cols with metadata, also nonnumeric
  #sampleID: unique identifier per sample
  #cohort: the groups you want to compare
  #Variable2: another group to compare
  
  
  #todo: add datastructure
 
    print(sampleID)
    input_heatmap <- input
    
    #annotation row
    print ( input_heatmap %>% group_by(sampleID) %>% filter(n()>1) %>% summarize(n=n()))
    rownames(input_heatmap) <- paste0(input_heatmap[[sampleID]],"_", input_heatmap[[cohort]])
    input_heatmap <- input_heatmap%>% dplyr::arrange(.data[[cohort]])
    
    if ( is.null(Variable2) == FALSE) {
      input_heatmap <- input_heatmap%>% dplyr::arrange(.data[[cohort]], .data[[Variable2]] )
    }
    
    input_heatmap_numbers <- input_heatmap %>% dplyr :: select (-c(metadata_columns))
    
    heatmap_final <- (as.matrix(input_heatmap_numbers )) 
    
    
    row.names(heatmap_final) <- paste0(input_heatmap[[sampleID]],"_", input_heatmap[[cohort]])
    
    
    
    p <- data.frame(Cohort1 = input_heatmap[[cohort]])
    rownames(p) <- paste0(input_heatmap[[sampleID]],"_", input_heatmap[[cohort]])
    p$Cohort1<- as.character(p$Cohort1)
    
    names(p)[1] <- cohort
    
    if ( is.null(Variable2) == FALSE) {
    p[[Variable2]] <- input_heatmap[[Variable2]]
    }
    
    print(p)
    print(rownames(heatmap_final))
    #my_colour = list(
     # Patient_Clusters= c("1" = "red" , "2" = "green", "3" = "blue", "4"= "violet", "5" = "orange" ),
     # outcome = c(survivor = "grey" , nonsurvivor = "black")
   # )
    
    
    plotheat = pheatmap::pheatmap(heatmap_final,
                                  #input_heatmap[selected_cols] %>% select (-c(Alias ,Celltype,Stimulation,outcome,cohort)),
                                  #input_heatmap[-c(55:61)],
                                  #correlation,
                                  #color = paletteer_c('pals::ocean.balance', n = 100),
                                  scale = "column",
                                  cellwidth = 15, cellheight = 6,
                                  annotation_colors = my_colour,
                                  annotation_row = p,
                                  # annotation_col = col,
                                  fontsize = 7,
                                  breaks = breaks,
                                  # display_numbers = TRUE,
                                  # cutree_rows = 6,
                                  cluster_rows = F,
                                  cluster_cols=T,
                                  main = title_heatmap,
                                  silent = F)
    #plotheat <- plotheat + theme_minimal()
    
    
    return(plotheat)
}
  

#calculates VIP scores of LPS_DA and log2fold changes
VIP_Log2fold <- function (wide_DF_normalised, wide_df_not_normalized = NULL, metadata_columns, datastructure = datastructure, VIPtitle,  VIPmin = 1, used_stat = "mean", ncomp = ncomp,groupingCount = 2) {
  wide_DF_normalised <- wide_DF_normalised
  
  test <- as.matrix(wide_DF_normalised %>% dplyr::select(-c(metadata_columns)))
  rownames(test) <- wide_DF_normalised[[datastructure$sampleID]]
  
  X <- test
  Y <- as.factor(wide_DF_normalised[[datastructure$cohort]])
  #imputate NAs
  X <- impute.nipals(X = X, ncomp = ncomp)
  
  #pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 
  srbct.splsda <- splsda(X, Y, ncomp = ncomp) 
  
  final.vip <- vip(srbct.splsda)
  
  print(final.vip[1:5])
  
  VIP<- as_tibble(final.vip, rownames = "features")%>% dplyr::select ( features,comp1) %>% arrange (comp1) %>% 
    dplyr::filter (comp1 > VIPmin)
  
  VIP$Variable <- factor(VIP$features, levels = VIP$features)
  
  input_2 <- wide_DF_normalised  %>% dplyr::select (VIP$features, metadata_columns)
  
   long_df <- create.long.Df(input_2, metadata_columns)
   
   if (is.null(wide_df_not_normalized)== F){
     wide_df_not_normalized <- wide_df_not_normalized %>% filter(sampleID %in% input_2$sampleID)
     wide_df_not_normalized <-  wide_df_not_normalized%>% dplyr::select (VIP$features, metadata_columns)
     long_df <-  create.long.Df(wide_df_not_normalized, metadata_columns)
     print("using raw data for logfold")
   }
  
  
  #print(long_df[1:6,])
  
 
  Log2foldchanges <- CalculateFoldChange ( long_df = long_df, groupingCount = groupingCount, datastructure= datastructure,used_stat)
  
  print(Log2foldchanges[1:5,])
  
  Log2foldchanges <- Log2foldchanges %>% left_join(VIP[c("Variable", "comp1")], by = "Variable" ) %>% rename(VIP= comp1)
  
  #VIP vs LOG2fold changes as a vlcano plot

  
  Log2foldchanges$color_breaks <- cut(Log2foldchanges$log2foldchanges, 
                                      breaks =  c(-Inf, -0.5, 0.5, Inf),
                                      labels = c("Low", "Mid", "High"))
  
  
  
  p<- ggplot(Log2foldchanges, aes(log2foldchanges, VIP, color=Log2foldchanges$color_breaks)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("Low" = "blue", "Mid" = "grey", "High" = "red"),
                       name = "Value Range",
                       breaks = c("Low", "Mid", "High"),
                       labels = c(paste("increased in",datastructure$firstcohort), "no change", paste("increased in", datastructure$secondcohort))) +
    
    theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
    theme(panel.grid.minor = element_blank())+
    
    geom_text_repel(hjust = 0.2, nudge_x = 0.75, size = 3, data=subset(Log2foldchanges, ((VIP  > 1 )& (log2foldchanges > 0.2 | log2foldchanges < -0.2))),
                    aes(log2foldchanges ,VIP,label=Variable, colour = "black"))+
    ggtitle(VIPtitle) +
    expand_limits(y = max(Log2foldchanges$VIP), x = Log2foldchanges$log2foldchanges)
  
  result <- list()
  result$plot <- p
  result$VIP <- VIP
  result$logfold <- Log2foldchanges
  return(result)
  
}

#PLot PCA

plot_PCA <- function(input_cluster_2, datastructure, metadata_columns,imputation = TRUE, ind.names = F, titlePCA = "PCA", components = c(1:2), ncomp =ncomp) {
  
  test <- as.matrix(input_cluster_2 %>% dplyr::select(-c(metadata_columns)))
  rownames(test) <- input_cluster_2[[datastructure$sampleID]]
  
  X <- test
  Y <- as.factor(input_cluster_2[[datastructure$cohort]])
  
  if (imputation == T) {
  #imputate NAs
  X <- impute.nipals(X = X, ncomp = ncomp)
  }
  pca.srbct = pca(X, ncomp = ncomp, center = TRUE, scale = TRUE) 
  p <- plot(pca.srbct)
  p2 <- plotIndiv(pca.srbct, group = Y, ind.names = ind.names,
            comp = components,# plot the samples projected
            legend = TRUE, 
            ellipse = TRUE,
           # style = 'graphics',
            title = titlePCA) # onto the PCA subspace
  plotlist <- list()
  plotlist$variance <- p
  plotlist$PCA <- p2
  
  return(plotlist)
}
  
#plot Least partial squares DA
plot_LPS_DA <- function(input_cluster_2, datastructure,metadata_columns, imputation = TRUE, ind.names = F, titlePCA = "PCA", components = c(1:2), ncomp = ncomp) {
  
  #datastructure$cohort is the grouping factor
  
  test <- as.matrix(input_cluster_2 %>% dplyr::select(-c(metadata_columns)))
  rownames(test) <- input_cluster_2[[datastructure$sampleID]]
  
  X <- test
  Y <- as.factor(input_cluster_2[[datastructure$cohort]])
  
  if (imputation == T) {
    #imputate NAs
    X <- impute.nipals(X = X, ncomp = ncomp)
  }
  srbct.splsda <- splsda(X, Y, ncomp = ncomp)
  
  p2 <- plotIndiv(srbct.splsda, group = Y, ind.names = ind.names,
                  comp = components,# plot the samples projected
                  legend = TRUE, 
                  ellipse = T,
                  style = "ggplot2",
                  title = titlePCA)  # onto the PCA subspace
  
  ?plotIndiv
  
  plotlist <- list()
  
  plotlist$LPSDA <- p2
  
  return(plotlist)
}


#reorders factor of long df as defined
reorder_long_df_factors <- function(long_df, order, datastructure) {
  
  result <- long_df %>%  mutate(cohort_as_factor = factor(long_df[[datastructure$cohort]], levels=order)) %>%ungroup()%>%
    select(-c(datastructure$cohort))
  result <- result %>% rename (!!datastructure$cohort  := cohort_as_factor)
  return(result)
}

?get_summary_stats
#calculates sumstats

calculate_Sumstats <- function (long_df,datastructure, grouping_vector, selected_stats = c("ci","mean", "n", "q1", "q3") ) {
  
  #' calculatis summary statistics with statix package (?get_summary_stats)
  #'
  #' @param long_df a df in long format, filter to selected parameters first
  #' @param datastructure list with datastructure
  #' @param grouping_vector vector containing the columns by which the data is getting grouped (values in datastructure)
  #' @param selected_stats vector with list of selected stats to calculate (values include: "full", "common", "robust", "five_number", "mean_sd", "mean_se", "mean_ci", "median_iqr", "median_mad", "quantile", "mean", "median", "min", "max")
  #'
  #' @return a df with sumstats 
  #' @export
  #'
  #' @examples 
  #' grouping_vector <- c(datastructure$Variable, datastructure$Variable2, datastructure$Variable3)
  #'  selected_stats <- c("ci","mean", "n", "q1", "q3")
  #'  calculate_Sumstats(long_df, datastructure, grouping_vector, selected_stats )
  #'  
  calculate_Sumstats(long_df, datastructure, grouping_vector, selected_stats )
  sumstats <- long_df %>% group_by(.dots = grouping_vector)%>% 
    get_summary_stats(datastructure[["measured"]], show = selected_stats)
  return (sumstats)
}





Overview.Histogram<- function (long_df,datastructure, group.colors,  
                               title = "Test", x_axis_label = "Measurement", headersize = 10, legendsize = 30,  alpha_value = 0.5) { 
  
  
  
  cohort <- datastructure$cohort
  
  #data in
  data <-long_df
  data[[cohort]] <- as.factor(data[[cohort]])
  
  
  #plotting
  p<-ggdensity(data,  x = datastructure$measured,  fill =datastructure$cohort, add = "mean", rug = TRUE )+#
    ggtitle(title) +
    
    
    theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
    
    xlab(x_axis_label)+
    theme(panel.grid.minor = element_blank())+
    scale_color_manual(values=group.colors)+
    scale_fill_manual( values= alpha(c(group.colors),alpha_value)
                       
    ) 
  return (p)
}



Overview.Plot<- function (long_df,datastructure, group.colors,  
                          title = "Test", headersize = 10, legendsize = 30,  alpha_value = 0.5) { 
  #for nice stats labels, rounded and styled
  
  
  cohort <- datastructure$cohort
  
  #data in
  data <-long_df
  data[[cohort]] <- as.factor(data[[cohort]])
  
  
  #plotting
  p<-ggplot(data, aes(y = .data[[datastructure$Variable]], x = .data[[datastructure$measured]],  fill =.data[[datastructure$cohort]] ))+#
    ggtitle(title) +
    ggtitle(title) +
    stat_boxplot() +
    
    
    theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
    theme(panel.grid.minor = element_blank())+
    scale_color_manual(values=group.colors)+
    scale_fill_manual( values= alpha(c(group.colors),alpha_value))+
    
    theme(strip.text.x=element_text(size= headersize),  legend.text=element_text(size=legendsize)) +
    theme(axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          #axis.ticks.x=element_blank()
          
    ) 
  
  
  return (p)
}



longitudinal.Plot<- function (long_df,datastructure, group.colors, sumstats, use_aggregated_data = F, 
                              title = "Test",y_label = "Measurement", headersize = 10, legendsize = 30, x_axis_labelsize = 10,
                              alpha_value = 0.5, pointsize = 3, dodgevalue = 0.5, stats_used = "median", show_n = F) { 
  
  
  #' plotting tool to generate a lineplot, with optional stats
  #'
  #' @param long_df the input df
  #' @param datastructure predefined datastructure: $cohort is x axis,  $measured is y axis, $Variable2 is color, $Variable3 is what is connected
  #' @param group.colors colors defined by function 
  #' @param sumstats stats defined by sumstats function
  #' @param use_aggregated_data plot mean with errorbar or multible line plots or "Box" = Boxplot
  #' @param title title 
  #' @param y_label label of y axis
  #' @param headersize size of text
  #' @param legendsize size of legend
  #' @param alpha_value alpha of line
  #' @param pointsize size of dot
  #' @param dodgevalue dodge od plots
  #' 
  #'
  #' @return ggplot object
  #' @export
  #'
  #' @examples longitudinal.Plot(cleaned, plateletdata,group.colors, sumstats, use_aggregated_data = F, title = "Test", y_label = "MFI", headersize = 10, legendsize = 30,  alpha_value = 0.8, pointsize = 4) 

  
  
  pd <- position_dodge(dodgevalue)
  
  if (use_aggregated_data == T)
  {
    
    plot <- ggplot(sumstats, aes(x=.data[[datastructure$cohort]], y=.data[[stats_used]],  colour=.data[[datastructure$Variable2]], group=.data[[datastructure$Variable3]])) + 
      
      geom_errorbar(aes(ymin= q1, ymax= q3, colour=.data[[datastructure$Variable2]]), alpha = 0.3, width=0.5, position=pd) +
      #geom_errorbar(aes(ymin= (lowerCI), ymax= (upperCI), colour=.data[[datastructure$Variable2]]), alpha = 0.3, width=0.5, position=pd) +
      
      geom_line(position=pd) +
      geom_point(position=pd, size=pointsize)+ 
      scale_color_manual(values=group.colors)+
      facet_wrap(as.formula(paste('~', datastructure$Variable)), scales = "free", ncol = 6) +
      theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
      theme(strip.text.x=element_text(size= headersize),  legend.text=element_text(size=legendsize),
            axis.text.x = element_text(size=x_axis_labelsize)) +
      theme(panel.grid.minor = element_blank())+ 
      ylab(y_label)+
      ggtitle(title)
    
    return (plot)
  }
  
  if (use_aggregated_data == "Box")
  {
  
  p<-ggplot(long_df, aes(x = .data[[datastructure$cohort]], y = .data[[datastructure$measured]],  fill =.data[[datastructure$Variable2]] ))+#
    ggtitle(title) +
   
    stat_boxplot(aes(colour = .data[[datastructure$cohort]], fill =.data[[datastructure$Variable2]])) +
    
    
    
    theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
    theme(panel.grid.minor = element_blank())+
    scale_color_manual(values=group.colors)+
    scale_fill_manual( values= alpha(c(group.colors),alpha_value))+

    theme(strip.text.x=element_text(size= headersize),  legend.text=element_text(size=legendsize), 
          axis.text.x = element_text(size=x_axis_labelsize)) +
    theme(axis.title.x=element_blank(),
        
    ) +
  
   facet_wrap(as.formula(paste('~', datastructure$Variable)), scales = "free", ncol = 4) +  # Boxplot width can be adjusted
    
    geom_jitter(show.legend=FALSE, aes(color=.data[[datastructure$Variable2]]), 
                position=position_dodge(width=0.8), 
                size=pointsize,
                alpha = alpha_value
                )+
    stat_summary(show.legend=FALSE, aes(group=.data[[datastructure$Variable2]], 
                     color=.data[[datastructure$Variable2]]), 
                 fun=stats_used, geom="line", 
                 position=position_dodge(width=0.8),
                 alpha = alpha_value
                 ) 
    # guides(color=FALSE, fill=FALSE)
  
  #display n on top of boxplot
  if(show_n == T) {
 p <- p + geom_text(data = sumstats, aes(x = .data[[datastructure$cohort]], y =upper_limit+1, color = .data[[datastructure$Variable2]], label = paste0("n:", n)), size= 2,check_overlap = F,
              vjust = 0,  position=position_dodge(width=0.8)) 
}
   return(p)
  
  }
  
  
  
  
  
  
  
  plot <- ggplot(long_df, aes(x=.data[[datastructure$cohort]], y=.data[[datastructure$measured]],  colour=.data[[datastructure$Variable2]], group=interaction(.data[[datastructure$Alias]], .data[[datastructure$Variable3]]))) + 
  
    
   geom_line(position=pd, alpha = alpha_value) +
    geom_point(position=pd, size=pointsize)+ 
    scale_color_manual(values=group.colors)+
    facet_wrap(as.formula(paste('~', datastructure$Variable)), scales = "free", ncol = 6) +
    theme_prism(base_size = 11, palette = "candy_bright", base_line_size = 0.5,axis_text_angle = 45, base_family = "sans")+
    theme(strip.text.x=element_text(size= headersize),  legend.text=element_text(size=legendsize), 
          axis.text.x = element_text(size=x_axis_labelsize)) +
    theme(panel.grid.minor = element_blank())+ 
    ggtitle(title)+
    ylab(y_label)
  
  return(plot)
  
}

#removes data which cannot be interpreted by stats functions
clean_data_for_stats <- function(long_df,datastructure, selected_cohort, groupingvariables)  {
  result <- long_df %>%
    filter(!!sym(datastructure$cohort) %in% selected_cohort) %>%
    group_by(across(all_of(groupingvariables))) %>%            # Group by feature and cohort
    filter(n() >= 2) %>%                     # Ensure there are at least 2 observations
    ungroup() %>%
    group_by(!!sym(datastructure$Variable)) %>%
    filter(n_distinct(!!sym(datastructure$cohort)) > 1) %>%       # Ensure feature exists in more than one cohort
    ungroup()
  if(nrow(result)== 0) {
    print("empty dataframe")
  }
  return(result)
}


