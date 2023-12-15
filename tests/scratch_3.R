


dfa_cluster <- dfa_agg %>% convert_to_sequence() %>% ungroup()
dfb_cluster <- dfb_agg %>% convert_to_sequence() %>% ungroup()


dfa_cluster <- dfa_agg %>% convert_to_sequence() %>% ungroup()
dm <- inter_sequence_distance(dfa_cluster %>% pull(sequence))
dm[is.na(dm)] = 0

res <- pam(dm, k = 3, diss = TRUE)
dfa_cluster$cluster_id <- res$clustering



dfa_cluster <- dfa_cluster %>%
                  group_by(cluster_id) %>%
                  select(-sequence_formatted, , cluster = cluster_id) %>%
                  nest(df_sequences=c("id", "nested_id", "sequence")) %>%
                  mutate(n = map_int(df_sequences, nrow),
                         df_sequences = map(df_sequences, function(df_sequence) {
                           class(df_sequence$sequence) <- c("Sequence_List", class(df_sequence$sequence))
                           names(df_sequence$sequence) <- df_sequence$id
                           df_sequence
                  })) %>%
                  arrange(desc(n)) %>%
                  ungroup()

class(dfa_cluster) <- c("Clustered_Dataframe", class(dfa_cluster))





## Unnests the sequences within each cluster, returns the -df_aggregated- structure EXCEPT FOR "date" and "period" while adding "cluster"
(dfa_clustered_kmed %>% select(-n) %>% unnest(cols = c(df_sequences)) %>% unnest(cols = c(sequence))) %>% unnest(cols = c(sequence))








attributes((dfa_clustered_kmed %>% filter_pattern(threshold = .2))$weighted_sequence[1][[1]])
attr((dfa_clustered_kmed %>% filter_pattern(threshold = .2))$weighted_sequence[1][[1]], "alignments")[1]



wseq <- dfa_clustered_kmed %>% get_weighted_sequence() %>% filter(cluster == 1) %>% pull(weighted_sequence)
wseq2 <- attr((dfa_patterns %>% filter(cluster == 1) %>% pull(weighted_sequence))[[1]], "alignments")





  elements <- unlist(map(wseq[[1]][1][[1]],
                         function(x) {
                           if (length(x$elements) == 0)
                             return(NULL)
                           x$elements
                         }))




sequence <- ((dfa_clustered_kmed %>% filter(cluster == 1) %>% pull(df_sequences))[[1]] %>% pull(sequence))[[1]]
sequence2 <- ((dfa_clustered_kmed %>% filter(cluster == 1) %>% pull(df_sequences))[[1]] %>% pull(sequence))[[2]]
w_sequence <- align_sequences(sequence,
                              sequence,
                              sorenson_distance)

for (i in 1:length(sequence)) {
  w_sequence[[i]]$itemset_weight <- 1
  w_sequence[[i]]$element_weights <- (w_sequence[[i]]$element_weights) / 2
}
attr(w_sequence, "n") <- 1
alignments <- attr(w_sequence, "alignments")[1]




((dfa_clustered_kmed %>% filter(cluster == 1) %>% pull(df_sequences))[[1]] %>% pull(sequence))[[1]][[1]]

(dfa_cluster %>% pull(sequence))[[1]][[2]]
(dfa_cluster %>% pull(sequence))[[1]][2]




## Trying some new with-cluster sum of squares calculation
dfa_clustered_kmed %>% filter(cluster == 1) %>% pull(df_sequences)

dfb %>% group_by(id) %>% nest() %>% mutate(dist = map_dbl(data, ~sum(as.matrix(dist(.x)^2)) / (2 * nrow(.x)))) %>% pull(dist)




(dfa_clustered_kmed %>% select(-n) %>% unnest(cols = c(df_sequences)) %>% unnest(cols = c(sequence))) %>% filter(id == 159) %>% pull(sequence)




(dfa_clustered_kmed %>% select(-n) %>% unnest(cols = c(df_sequences)) %>% unnest(cols = c(sequence))) %>% pull(sequence)













