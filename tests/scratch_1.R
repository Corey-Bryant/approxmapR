## TO DO:
##  (1) Add tie breaker for when an -id- has the same distance to > 1 consensus patterns
##  (2) Add ability to retreive -period- information from the consensus pattern


library(devtools)
#devtools::install_github("Corey-Bryant/approxmapR")
install_local("C:/Users/Corey/Documents/GitHub/approxmapR", force = TRUE)
#devtools::install_github("pinformatics/approxmapR")



library(usethis)
library(fpc)
library(cluster)



library(approxmapR)
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)
library(readxl)
#library(xlsx)
library("igraph")








output_path = "C:\\Users\\Corey\\Documents\\Approxmap\\PROCESS\\Development 2023"



### Loading Data and Sub-setting Data to Sta3n = 506 ###
df <- read.csv("C:\\Users\\Corey\\Documents\\Approxmap\\PROCESS\\Data\\Analysis_FY18_level1_CPT_events.csv")
df <- subset(df, IPflag != 1)
names(df) <- c("sta3n", "id", "date", "event", "IPflag")


### Subsetting data and selecting random (reproducible) sample
dfb <- subset(df, sta3n == 631) %>% select("id", "date", "event")
dfb <- process_varying_aggregation(dfb, scheme = 1) %>% select(id, date, period, event)
dfb_agg <- dfb %>% pre_aggregated(include_date = TRUE, summary_stats = TRUE, output_directory = output_path)



dfa <- subset(df, sta3n == 512) %>% select("id", "date", "event")
dfa <- process_varying_aggregation(dfa, scheme = 1) %>% select(id, date, period, event)
dfa_agg <- dfa %>% pre_aggregated(include_date = TRUE, summary_stats = TRUE, output_directory = output_path)

dfa_ktable <- dfa_agg %>% find_optimal_k(clustering = 'k-medoids', min_k = 2, max_k = 8,
                                         save_table = TRUE, use_cache = TRUE,
                                         file_name = "Dataset A Optimal K.csv" ,
                                         output_directory = output_path)
plot_ktable(dfa_ktable)

dfa_siloi <- dfa_ktable$silhouette_object[[2]]
plot_silhouette(dfa_siloi)
dfa_opt_k <- dfa_ktable$k[[which.max(dfa_ktable$average_silhouette_width)]]
dfa_opt_k <- 3

dfa_clustered_kmed <- dfa_agg %>% cluster_kmedoids(k = dfa_opt_k, use_cache = TRUE)
dfa_patterns <- dfa_clustered_kmed %>% filter_pattern(threshold = .2)
dfa_patterns %>% generate_reports(sil_table = dfa_siloi, output_directory = output_path)




##########################################
## Consensus Distance Begins Below Here ##
##########################################
dfa_patterns_formatted <- dfa_patterns %>% format_sequence()

dfa_output <- dfa_patterns_formatted %>%
              mutate(id = paste0("c", cluster)) %>%
              select(id, consensus_pattern)

dfa_consensus <- dfa_output %>% convert_to_events("id", "consensus_pattern")

dfa_consensus_clustered <- dfa_consensus %>% convert_to_sequence() %>% ungroup()





## Works - original before expanding
consensus_sequence_distance <- function(sequence_list, consensus_sequence_list) {

  total_id_sequences = sum(sequence_list$n)

  list_cluster <- numeric(total_id_sequences)
  list_id <- character(total_id_sequences)
  list_consensusid <- character(total_id_sequences)
  list_distance <- character(total_id_sequences)


  # For each individual sequence, calculate the distance from each consensus pattern
  #   and assign
  count <- 0
  for (c in sequence_list$cluster) {

    current_cluster_ids <- subset(sequence_list, cluster == c)$df_sequences[[1]]

    for (i in current_cluster_ids$id) {

      current_cluster_ids_sequence <- current_cluster_ids %>% filter(id == i) %>% pull(sequence)

      for (j in consensus_sequence_list$id) {
        count <- count + 1
        current_consensus_sequence <- dfa_consensus_clustered %>% filter(id == j) %>% pull(sequence)

        current_sequences <- rbind(current_consensus_sequence,
                                   current_cluster_ids_sequence)
        class(current_sequences) <- "Sequence_List"

        print(paste0("Calculating distance for id = ", i, " & consensus pattern id = ", j))
        current_id_consensus_distance <- inter_sequence_distance(current_sequences)

        list_cluster[count] <- c
        list_id[count] <- i
        list_consensusid[count] <- j
        list_distance[count] <- as.numeric(current_id_consensus_distance[1,2])

      }
    }
  }


  tibble("cluster" = list_cluster,
         "id" = list_id,
         "consensusid" = list_consensusid,
         "distance" = list_distance)
}

# dfa_clustered_kmed is the R object that is produced after running the clustering function and is what is passed
# to the filter_pattern function which pulls the consensus patter
id_consensus_distance <- consensus_sequence_distance(dfa_clustered_kmed, dfa_consensus_clustered)







## WORKS! Expanding the original function from above
consensus_sequence_distance <- function(df_aggregated, consensus_sequence_list) {

  total_id_sequences = n_distinct(df_aggregated$id)

  list_cluster <- numeric(total_id_sequences)
  list_id <- character(total_id_sequences)
  list_consensusid <- character(total_id_sequences)
  list_distance <- character(total_id_sequences)


  # For each individual sequence, calculate the distance from each consensus pattern
  #   and assign
  df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()

  count <- 0
  for (i in df_cluster$id) {

      current_cluster_ids_sequence <- df_cluster %>% filter(id == i) %>% pull(sequence)

      for (j in consensus_sequence_list$id) {
        count <- count + 1
        current_consensus_sequence <- consensus_sequence_list %>% filter(id == j) %>% pull(sequence)

        current_sequences <- rbind(current_consensus_sequence,
                                   current_cluster_ids_sequence)
        class(current_sequences) <- "Sequence_List"

        print(paste0("Calculating distance for id = ", i, " & consensus pattern id = ", j))
        current_id_consensus_distance <- inter_sequence_distance(current_sequences)

        list_id[count] <- i
        list_consensusid[count] <- j
        list_distance[count] <- as.numeric(current_id_consensus_distance[1,2])

      }
    }


  # Creates a Tibble frame which contains:
  #   the original cluster (if any) the individual sequence came from
  #
  id_consensus_distance <- tibble("id" = list_id, "consensusid" = list_consensusid, "distance" = list_distance)

  id_consensus_distance <- id_consensus_distance %>% group_by(id) %>%
                            mutate(minDistance = min(distance),
                                   distanceMatch = case_when(distance == minDistance ~ 1, TRUE ~ 0),
                                   nDup = sum(distanceMatch),
                                   assignedConsensusCluster = case_when(distanceMatch == 1 & nDup == 1 ~ consensusid,
                                                                        TRUE ~ "NEED SENSITIVITY ANALYSIS")) %>%
                            arrange(desc(nDup), id, consensusid)

  id_consensus_distance <- id_consensus_distance %>% filter(distance == minDistance)
  id_consensus_distance_dups <- id_consensus_distance %>% filter (nDup > 1)


  return(id_consensus_distance)

}


id_consensus_distance <- consensus_sequence_distance(dfb_agg, dfa_consensus_clustered)







# Artificially creating a duplicate
id_consensus_distance[2, 4] <- id_consensus_distance[1, 4]
id_consensus_distance[5, 4] <- id_consensus_distance[4, 4]

# Filters data down to cases where Distance == minDistance
#   Note: Duplicates are possible at this point. If there are duplicates present, need to
#         conduct sensitivity analysis of assigning the ID to the consensus sequences with the
#         same minDistance --> assign based on assignment criterion
id_consensus_distance <- id_consensus_distance %>% group_by(id) %>%
                          mutate(minDistance = min(distance),
                                 distanceMatch = case_when(distance == minDistance ~ 1, TRUE ~ 0),
                                 nDup = sum(distanceMatch),
                                 assignedConsensusCluster = case_when(distanceMatch == 1 & nDup == 1 ~ consensusid,
                                                                      TRUE ~ "NEED SENSITIVITY ANALYSIS")) %>%
                          arrange(desc(nDup), id, consensusid)

id_consensus_distance <- id_consensus_distance %>% filter(distance == minDistance)
id_consensus_distance_dups <- id_consensus_distance %>% filter (nDup > 1)








id_consensus_distance_dups
id_consensus_distance






# Use the aggregated data -> convert to a sequence -> ungroup
dfa_cluster <- dfa_agg %>% convert_to_sequence() %>% ungroup()
dfa_dm <- inter_sequence_distance(dfa_cluster %>% pull(sequence))


## Creating visualization for ID distance from Consensus ##
test2 <- test
names(test2) <- c("Source (cluster)", "Source (id)", "Target", "Distance")

g <- graph_from_data_frame(test2, directed = F)
g_cluster_edge <- cluster_edge_betweenness(g)


plot(g)

plot(g, mark.groups = test2["Source (cluster)"], layout=layout_as_tree)

plot(g_cluster_edge, g,  mark.groups = test2["Source (cluster)"], edge.label=test2$Distance)




g <- graph_from_data_frame(select(test2, -"Source (cluster)"), directed = F)
plot(g)

make_clusters(g,
              membership = select(test2, "Source (cluster)"))


