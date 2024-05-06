
library("approxmapR")
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)
library(readxl)
#library(xlsx)
library("igraph")





# Subsetting data
dfb <- subset(df, sta3n == 631) %>% select("id", "date", "event")
dfb <- process_varying_aggregation(dfb, scheme = 1) %>% select(id, date, period, event)
dfb_agg <- dfb %>% pre_aggregated(include_date = TRUE, summary_stats = TRUE, output_directory = output_path)

dfa <- subset(df, sta3n == 512) %>% select("id", "date", "event")
dfa <- process_varying_aggregation(dfa, scheme = 1) %>% select(id, date, period, event)
dfa_agg <- dfa %>% pre_aggregated(include_date = TRUE, summary_stats = TRUE, output_directory = output_path)



# Finding optimal k value for clustering, for each data frame

## dataframe a
dfa_ktable <- dfa_agg %>% find_optimal_k(clustering = 'k-medoids', min_k = 2, max_k = 8,
                                         save_table = TRUE, use_cache = TRUE,
                                         file_name = "Dataset A Optimal K.csv" ,
                                         output_directory = output_path)
plot_ktable(dfa_ktable)

dfa_opt_k <- dfa_ktable$k[[which.max(dfa_ktable$average_silhouette_width)]]
dfa_opt_k <- 3
dfa_siloi <- dfa_ktable$silhouette_object[[2]]
plot_silhouette(dfa_siloi)


dfa_clustered_kmed <- dfa_agg %>% cluster_kmedoids(k = dfa_opt_k, use_cache = TRUE)
dfa_patterns <- dfa_clustered_kmed %>% filter_pattern(threshold = .2)
dfa_patterns %>% generate_reports(sil_table = dfa_siloi, output_directory = output_path, end_filename_with = "_sta3n512")


## dataframe b
dfb_ktable <- dfb_agg %>% find_optimal_k(clustering = 'k-medoids', min_k = 2, max_k = 8,
                                         save_table = TRUE, use_cache = TRUE,
                                         file_name = "Dataset B Optimal K.csv" ,
                                         output_directory = output_path)
plot_ktable(dfb_ktable)

dfb_opt_k <- dfb_ktable$k[[which.max(dfb_ktable$average_silhouette_width)]]
dfb_opt_k <- 3
dfb_siloi <- dfb_ktable$silhouette_object[[2]]
plot_silhouette(dfb_siloi)


dfb_clustered_kmed <- dfb_agg %>% cluster_kmedoids(k = dfb_opt_k, use_cache = TRUE)
dfb_patterns <- dfb_clustered_kmed %>% filter_pattern(threshold = .2)
dfb_patterns %>% generate_reports(sil_table = dfb_siloi, output_directory = output_path, end_filename_with = "_sta3n631")





##########################################
## Consensus Distance Begins Below Here ##
##########################################

## NOTE: In the consensus data frames, 'id' = the cluster
dfa_consensus <- dfa_patterns %>% convert_consensus_to_events()
dfb_consensus <- dfb_patterns %>% convert_consensus_to_events()

dfb_id_consensus_distance <- consensus_sequence_distance(dfb_agg, dfa_consensus)
dfa_id_consensus_distance <- consensus_sequence_distance(dfa_agg, dfb_consensus)



