
library("approxmapR")
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
library(tidyverse)
library(readxl)
#library(xlsx)
library("igraph")




output_path = "C:\\Users\\Corey\\Documents\\Approxmap\\PROCESS\\Development 2023"


data("demo1")
demo1 <- data.frame(do.call("rbind", strsplit(as.character(demo1$id.date.item), ",")))
names(demo1) <- c("id", "period", "event")


agg_df <-  demo1 %>% aggregate_sequences(format = "%m/%d/%Y",
                                         unit = "month",
                                         n_units = 1,
                                         summary_stats = FALSE, include_date = TRUE)


clustered_kmed <- agg_df %>% cluster_kmedoids(k = 6, use_cache = TRUE)

patterns <- clustered_kmed %>% filter_pattern(threshold = .2)

patterns




##########################################
## Consensus Distance Begins Below Here ##
##########################################
convert_consensus_to_events <- function(df_patterns) {

list_cluster <- numeric(1)
list_period <- numeric(1)
list_elements <- character(1)


ix <- 0
for (c in df_patterns$cluster) {

    current_concensus <- subset(df_patterns, cluster == c)$consensus_pattern[[1]]
    current_concensus_length <- length(current_concensus)

    for (i in 1:current_concensus_length) {

      for (e in current_concensus[[i]]$elements) {

        ix <- ix + 1

        list_cluster[ix] <- c
        list_elements[ix] <- e
        list_period[ix] <- current_concensus[[i]]$period

      }
    }
}


  tibble("id" = list_cluster,
         "date" = list_period,
         "period" = list_period,
         "event" = list_elements)
}

onsensus_sequence_distance <- function(df_aggregated, consensus_sequence_list) {

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

        current_sequences <- c(current_consensus_sequence,
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


consensus_events <- patterns %>% convert_consensus_to_events()





dfa_consensus <- dfa_patterns %>% convert_consensus_to_events()
dfa_consensus_clustered <- dfa_consensus %>% convert_to_sequence() %>% ungroup()

id_consensus_distance <- consensus_sequence_distance(dfb_agg, dfa_consensus_clustered)








