

setwd("/Volumes/Krell_SSD/Projects/Research/Main_Projects/2022/COVID") 
source("code/Data_cleaning/data_scaling_functions.R")

##############load MS data #########################

metadata_location <- "metadata"
rawdata <-"data/rawdata/MassSpec"
MS_preprocessed <-"data/processed_data/MassSpec"
Experiment_files <- "/QS/"
filname_combined_data_MS <- "final_combined_data_lipidomics_16_06_23.xlsx"
data_location_metadata <- "data/clinical_data/PDMS_data_cleaned"

lipidomics <- read.xlsx(paste0(MS_preprocessed, "/", filname_combined_data_MS)) %>% distinct(Alias,Celltype, Stimulation, .keep_all = TRUE) 
lipidomics <- lipidomics %>%   mutate(Alias = str_replace_all(Alias, "PGLB", "PQLB")) #missspelled!

MS_experiment_lipid$metadata$clinical_data

#############
#load experiment#
MS_experiment_lipid <- qread(paste0(MS_preprocessed, Experiment_files,"lipid_MS_experiment.qs"))

#save experiment#
qsave(MS_experiment_lipid, paste0(MS_preprocessed, Experiment_files,"lipid_MS_experiment.qs"))
############

#colnames
colnames(lipidomics) <-  MS_experiment_lipid$metadata$col_names


#metadata
medlabsobs_covid <- read.xlsx( paste0(data_location_metadata,'/medlabsobs_covid_13_04_23.xlsx')) %>% full_join( metadata_covid, by = "ADM_ADMISSIONNR")

#### metadata####
updated_metadata <- read.xlsx( paste0(data_location_metadata,"/metadata_updated.xlsx"))

 #join the metadata and data matrix
raw_input <- Join.Metadata(lipidomics,updated_metadata , ID1 = "Alias", ID2 ="Celltype") 


#define the metadata cols
metadata_columns <- (c("Alias" ,"Celltype","Stimulation","outcome","cohort", "ADM_ADMISSIONNR", "sampleID"))


#######. ##########
`%notin%`<- Negate(`%in%`)
input <- raw_input
input <- input %>% filter (Celltype %notin%  c("SUP"))  %>% filter (Stimulation == "UN")


#create unique row-ID 
input$sampleID <- paste0(input$Alias, input$Celltype)
#check for duplicates
input  %>% group_by(sampleID) %>% filter(n()>1) %>% summarize(n=n())

input %>% filter(Alias == "PY3M")



# **********perform processing****************#
projectfile <- list()

projectfile$rawdata <- input

for (i in unique(input$Celltype) ){
#for (i in c("CD16", "CD14") ){
  
# 1)  filter df by metadata
  `%notin%` <- Negate(`%in%`)
  
  projectfile[[i]]$filtered_df <- input %>% 
    #filter (Celltype %in% c("CD16", "CD14")) %>% 
    filter (Alias %notin%  c("H8LB")) %>% 
    filter (Stimulation == "UN") %>% 
    filter (Celltype == i)
    
   

  
    
 
  #convert to numeric
  projectfile[[i]]$filtered_df <- ToNumeric(projectfile[[i]]$filtered_df)

  #remove missing data
  projectfile[[i]]$selected_df <- RemoveMissingData(projectfile[[i]]$filtered_df, limit = 0.8)

  #replace missing data by 1/4th of col minimum
  projectfile[[i]]$replaced_df <- ReplaceMissingByLoD(projectfile[[i]]$selected_df,4)
  
  #make long df of data 
  projectfile[[i]]$replaced_long_df <- create.long.Df(projectfile[[i]]$replaced_df, metadata_columns)
  
  
  #log2 transform 
  projectfile[[i]]$transformed_df <- log2transform(projectfile[[i]]$replaced_df)
  
  #normalize data
  projectfile[[i]]$Normalized_df<- Normalizer(projectfile[[i]]$transformed_df, scaleNorm = "ParetoNorm")
  
  #make long df
  projectfile[[i]]$combined_transformed_Normalized_long_df <- create.long.Df(projectfile[[i]]$Normalized_df, metadata_columns)
  
  #make one big df out of the other preselected columns
  
  projectfile$combined_replaced_total <- bind_rows(projectfile[[i]]$replaced_df, projectfile$combined_replaced)  %>% distinct( .keep_all = TRUE)
  projectfile$combined_replaced <- RemoveColswithMissingData(projectfile$combined_replaced_total)
  projectfile$combined_transformed_df <- log2transform(projectfile$combined_replaced)
  projectfile$combined_transformed_Normalized_df<- Normalizer(projectfile$combined_transformed_df, scaleNorm = "ParetoNorm")
  projectfile$combined_transformed_Normalized_df_long_df <- create.long.Df(projectfile$combined_transformed_Normalized_df,metadata_columns)
  
} 

#stats

datastructureinput_stats <- list()
datastructureinput_stats$cohort <- "outcome"
datastructureinput_stats$Variable <- "Variable"
datastructureinput_stats$Variable2 <- "cohort"
datastructureinput_stats$measured <- "Measurement"


#compare total lipids per celltype
#input_all <- projectfile$combined_transformed_Normalized_df_long_df #%>% dplyr::group_by(Alias ,.data[[cohort]], ,.data[[Variable]]) %>%   
  #dplyr::summarize(Measurement= mean(.data[[measured]], na.rm=TRUE)) %>% ungroup()


input_all <- create.long.Df( projectfile$combined_transformed_Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")
order2 <-  c("CD14", "CD16", "CD19", "CD8", "CD4", "CD66B") 

input_all <- reorder_long_df_factors(input_all, order2,datastructureinput_stats )

stats <- kruskal_testing(input_all, groupingCount = 2, posthoc = T, datastructureinput_stats) 



stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.05, x = "Celltype",  use_posthoc = T,  datastructureinput_stats)









#CD14
input_CD14 <- create.long.Df(projectfile$CD14$Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")

input_CD14 <- reorder_long_df_factors(input_CD14, order,datastructureinput_stats )

stats <- kruskal_testing(input_CD14, groupingCount = 1, posthoc = T, datastructureinput_stats) 

 

stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.05, x = "outcome",  use_posthoc = T,  datastructureinput_stats)





#CD16

input_CD16 <- create.long.Df(projectfile$CD16$Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")

input_CD16 <- reorder_long_df_factors(input_CD16, order,datastructureinput_stats )

#by cohort
stats <- kruskal_testing(input_CD16, groupingCount = 1, posthoc = T, datastructureinput_stats) 



stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.05, x = "outcome",  use_posthoc = T,  datastructureinput_stats)





#CD4

input_CD4 <- create.long.Df(projectfile$CD4$Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")

input_CD4 <- reorder_long_df_factors(input_CD4, order,datastructureinput_stats )

#by cohort
stats <- kruskal_testing(input_CD4, groupingCount = 1, posthoc = T, datastructureinput_stats) 



stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.05, x = "outcome",  use_posthoc = T,  datastructureinput_stats)



#CD19

input_CD19 <- create.long.Df(projectfile$CD19$Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")

input_CD19 <- reorder_long_df_factors(input_CD19, order,datastructureinput_stats )

#by cohort
stats <- kruskal_testing(input_CD19, groupingCount = 1, posthoc = T, datastructureinput_stats) 
stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.05, x = "outcome",  use_posthoc = T,  datastructureinput_stats)

#CD8#####

input_CD8 <- create.long.Df(projectfile$CD8$Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")

input_CD8 <- reorder_long_df_factors(input_CD8, order,datastructureinput_stats )

#by cohort
stats <- kruskal_testing(input_CD8, groupingCount = 1, posthoc = T, datastructureinput_stats) 
stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.06, x = "outcome",  use_posthoc = T,  datastructureinput_stats)


#CD66B#####

input_CD66B <- create.long.Df(projectfile$CD66B$Normalized_df,metadata_columns)
order <-  c("control", "survivor", "nonsurvivor")

input_CD66B <- reorder_long_df_factors(input_CD66B, order,datastructureinput_stats )
input_CD66B <- RemoveZeroVariance(input_CD66B, datastructureinput_stats)

#by cohort
stats <- kruskal_testing(input_CD66B, groupingCount = 1, posthoc = T, datastructureinput_stats) 
stat.test <-Sig.For.Plotting (stats , scaling_of_plot = "free", filter_of_sig = 0.06, x = "outcome",  use_posthoc = T,  datastructureinput_stats)




###########################
datastructureinput_stats$cohort <- "cohort"
datastructureinput_stats$firstcohort <- "control"
datastructureinput_stats$secondcohort <- "patient"


students_t<- Students_T_Test(MS_experiment_lipid$, groupingCount = 1, datastructureinput_stats)

log2foldchanges <- CalculateFoldChange (MS_experiment_lipid$plotting$untouched, groupingCount = 1,datastructureinput_stats)



#stat.test <-stat.test %>% dplyr:: rename ( Variable = MC)
full_stats <- log2foldchanges%>% left_join( students_t,  by = c("Variable" )) # by = c(Variable,Variable2 ))

#removes infinite values and NAs, joins with other stats
is.na(full_stats) <- do.call(cbind,lapply(full_stats, is.infinite))
full_stats <- drop_na(full_stats)
#full_stats <- full_stats %>% left_join( highandlow[,c("Variable", "Celltype","Firstcohort_mean","Secondcohort_mean")],  by = c(Variable))



#CD14
Variable = "Variable"
Variable2 = "Celltype"
measured = "Measurement"  #name of the measured variable column
cohort = "outcome"
datastructure$sampleID <- "sampleID"




#######################
datastructureinput_stats$cohort <- "cohort"

input <- create.long.Df (projectfile$CD4$Normalized_df, metadata_columns)
projectfile$CD4$stats <- Students_T_Test(input, groupingCount = 1, datastructureinput_stats) %>%
  select(c(Variable, .y., n1, n2, p, p.adj, p.adj.signif, conf.low,conf.high, df, method, rowID ))

firstcohort <- "control"
secondcohort <- "patient"

full_stats <- CalculateFoldChange (projectfile$CD4$replaced_long_df, groupingCount = 1,datastructureinput_stats) %>%
  left_join( projectfile$CD4$stats,  by = c("Variable" ))
#removes infinite values and NAs, joins with other stats
is.na(full_stats) <- do.call(cbind,lapply(full_stats, is.infinite))
full_stats <- drop_na(full_stats)


full_stats %>%
   select (c(Variable, log2foldchanges, Firstcohort_mean,	Secondcohort_mean, n1,	n2, p)) %>% 
   dplyr::rename( "Control:mean"  = Firstcohort_mean,  "Cov19:mean" = Secondcohort_mean) %>% 
  gt() %>%
  #gt::data_color(
   # columns = p ,
    #rows = currency < 50,
  #  method = "numeric",
    #palette = c("red", "white","white", "transparent"),
  #  color=scales::col_bin(
   #   bins=c(0, 0.01, 0.057, 1),
   #   palette = c("red", "orange", "white")),
    #domain = c(0,0.05,0.06,1)
 # )%>%
  tab_footnote( #adds ref number and footnote
    footnote = "Students T Test,  2-tailed",
    locations = cells_body(columns = p, rows = (1: length(full_stats$p)))
)%>%
 tab_header(
    title = md("**T Helper Cells**"), #markdown or #html("Wind,<br>mph")
subtitle =md("healthy controls, severe COV19 cases"))%>%
gtsave(filename ="Stats_CD4B.pdf", expand = 10)





#data in network
#stats: from scaled data, mean and log2 fold: unscaled data
library(RCy3)
input_network <- full_stats %>%
  #filter ( Celltype == "CD14")%>%
  select(c(log2foldchanges,Variable, p, p.adj,Firstcohort_mean,Secondcohort_mean))

loadTableData(input_network,data.key.column="Variable", table.key.column = "label")






#VIP, logfoldchange needs updated datastructure
datastructure <- list ()
datastructure$firstcohort <- "control" 
 datastructure$secondcohort <- "patient" 
 datastructure$sampleID <- "sampleID"
 datastructure$cohort <- "cohort"
 datastructure$Variable <-"Variable"
 datastructure$Variable2 <- "Celltype"
 datastructure$measured <- "Measurement"

#
pp <- VIP_Log2fold(wide_DF_normalised = projectfile$CD14$Normalized_df ,wide_df_not_normalized = projectfile$CD14$replaced_df ,  metadata_columns, datastructure, VIPtitle = "B Cells", VIPmin = 1)
pp$plot

ggsave( "VIP_CD8.pdf" ,height = 10, width = 10, limitsize = FALSE)

pp$logfold
pp$VIP$features

cohort <- "cohort"

plot_LPS_DA(input =projectfile$combined_transformed_Normalized_df , datastructure, imputation = TRUE, ind.names = F, titlePCA = "Comparison of Celltypes", components = c(1:2))


#heatmap, set annot colors correctly by datastructure######
my_colour = list(
cohort = c(survivor = "grey" , nonsurvivor = "black", control = "green"),
cohort2= c(CD14 = "red", CD16 = "orange", CD19 = "blue" , CD4 = "lightblue", CD8 = "pink", CD66B = "yellow")
)


#automatic color setting 
a <- create_color_mapping(projectfile$combined_transformed_Normalized_df,datastructure, "Variable2", palette_name = "ggsci::nrc_npg", palette_style = "discrete") 
b <- create_color_mapping(projectfile$combined_transformed_Normalized_df,datastructure, "cohort", palette_name = "scico::berlin", palette_style = "continuous") 
my_colour <- list(a,b)

#names of color vectors must equal to annotation col names
names(my_colour)[2] <- datastructure$cohort
names(my_colour)[1] <- datastructure$Variable2



# cohort and Variable2 are the annotation rows 
c <- Heatmap.simple( projectfile$CD66B$Normalized_df %>% select(c(pp$VIP$features,metadata_columns)), metadata_columns, breaks = seq(-2, 2, length.out = 100), 
                           sampleID = "sampleID", cohort = "outcome", Variable2 = "Celltype", my_colour, title_heatmap = "title")

ggsave(plot = c, "Heatmap_CD8.pdf" ,height = 20, width = 20, limitsize = FALSE)



datastructure <- list()
datastructure$Variable <-"Variable"
datastructure$Variable2 <- "cohort"
datastructure$measured <- "Measurement"
datastructure$cohort <- "outcome"
datastructure$sampleID <-"ADM_ADMISSIONNR"
datastructure$firstcohort <- "CD14"
datastructure$secondcohort <- "CD16"


#https://github.com/EmilHvitfeldt/r-color-palettes
#create color mapping
a <- create_color_mapping(projectfile$combined_transformed_Normalized_df,datastructure, "Variable2", palette_name = "ggsci::nrc_npg", palette_style = "discrete") 
b <- create_color_mapping(projectfile$combined_transformed_Normalized_df,datastructure, "cohort", palette_name = "scico::berlin", palette_style, palette_style = "continuous") 
group.colors <- c(a,b)


group.colors <- c( "control" ="blue","patient" = "violet", "CD14" = "red", "CD16" = "orange", "CD19" = "blue" , "CD4" = "lightblue", "CD8" = "pink", "CD66B" = "yellow")


box <- Boxplots_simple( create.long.Df(projectfile$combined_transformed_Normalized_df),datastructure, stat.test = NULL, group.colors, order = F,  dotplot = F, P_values = F, title = "Test", headersize = 10, legendsize = 30, alpha_dots = 0.5, alpha_value = 1)
ggsave(plot = box, "test3.pdf" ,height = 60, width = 20, limitsize = FALSE)


datastructure_stats <- datastructure
datastructure_stats$cohort <- "Celltype"
datastructure_stats$Variable2 <- "cohort"


stats <- RemoveZeroVariance(create.long.Df(projectfile$combined_transformed_Normalized_df), datastructure = datastructure_stats)

model  <- Students_T_Test(stats, groupingCount = 2, datastructure_stats) 
model <- kruskal_testing(stats, groupingCount = 2, posthoc = T, datastructure_stats) 

stat.test <-Sig.For.Plotting (model , scaling_of_plot = "free", filter_of_sig = 0.05 , x = "cohort", use_posthoc = F, datastructure_stats)
