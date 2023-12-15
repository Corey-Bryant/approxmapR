

> dfa_consensus_clustered
# A tibble: 2 x 3
  id    sequence       sequence_formatted
  <chr> <Sqnc_Lst>     <chr>
1 c1    <Sequence [1]> <(PTMH)>
2 c2    <Sequence [2]> <(MHXP) (OVPC)>






for (c in 1:dim(dfa_clustered_kmed)[1]) {

  current_cluster <- dfa_clustered_kmed$cluster[c]
  current_cluster_ids <- subset(dfa_clustered_kmed, cluster == current_cluster)$df_sequences[[1]]

  print(paste0("i = ", c, " cluster = ", current_cluster))
  for (i in current_cluster_ids$id) {

    print(paste0("Current cluster = ", current_cluster, " Current id = ", i))
    print(current_cluster_ids %>% filter(id == i) %>% pull(sequence))

  }

}




current_cluster_ids <- subset(dfa_clustered_kmed, cluster == 2)$df_sequences[[1]]
current_sequences <- rbind(dfa_consensus_clustered %>% filter(id == "c1") %>% pull(sequence),
                           current_cluster_ids %>% filter(id == 404) %>% pull(sequence))


current_id_consensus_distance <- inter_sequence_distance(current_sequences)



data("demo1")
demo1 <- data.frame(do.call("rbind", strsplit(as.character(demo1$id.date.item), ",")))
names(demo1) <- c("id", "period", "event")


agg <- demo1 %>% aggregate_sequences(format = "%m/%d/%Y", unit = "day", n_units = 1)

ktable <- agg %>% find_optimal_k(clustering = 'k-medoids', min_k = 2, max_k = 8,
                                         save_table = TRUE, use_cache = TRUE,
                                         file_name = "Demo1 Optimal K.csv" ,
                                         output_directory = output_path)
plot_ktable(ktable)

siloi <- ktable$silhouette_object[[2]]
plot_silhouette(siloi)
opt_k <- ktable$k[[which.max(ktable$average_silhouette_width)]]


clustered_kmed <- agg %>% cluster_kmedoids(k = opt_k, use_cache = TRUE)
patterns <- clustered_kmed %>% filter_pattern(threshold = .3)
patterns %>% generate_reports(sil_table =siloi, output_directory = output_path)




##################################
#  CREATING DATA FOR TEST CASES  #
##################################

#### Adding sequence
## Copying ID = 8 to switch the last 2 events
to_add <- demo1 %>% filter(id == 8)
to_add$id <- 808
to_add[6, 3] <- "M"
to_add[7, 3] <- "L"


demo1 <- rbind(demo1, to_add)

(demo1 %>% filter(id == 808))[, 1:3]




#### Adding sequence
## Copying ID 8 again and removing some event sets
to_add <- demo1 %>% filter(id == 8)
to_add$id <- 8082
to_add <- to_add[-c(5), ]

demo1 <- rbind(demo1, to_add)




#### Adding sequence
## Copying ID = 7 to test specificity of regex; only keeping the first event set
to_add <- (demo1 %>% filter(id == 7))[1:3, 1:3]
to_add$id <- 707

demo1 <- rbind(demo1, to_add)




#### Adding sequence
to_add <- demo1 %>% filter(id == 8)
to_add$id <- "900"


to_add[1:4, 3] <- c("I", "K", "Z", "M")
to_add[3:4, 2] <- c("12/20/2015", "12/21/2015")

to_add[6, 2:3] <- c("7/23/2016", "B")
to_add[7, 2:3] <- c("7/24/2016", "C")

to_add <- rbind(to_add, c("900", "2/03/2017", "Z"))
to_add <- rbind(to_add, c("900", "2/03/2017", "I")) # If multiple events occur on same day then it choose the first one (alphabetically sorted)
to_add <- rbind(to_add, c("900", "2/04/2017", "Z")) # If multiple same events occur in event set, it keeps the first one (time sorted)
to_add <- rbind(to_add, c("900", "2/05/2017", "Z"))
to_add <- rbind(to_add, c("900", "2/05/2017", "J"))
to_add <- rbind(to_add, c("900", "2/06/2017", "Z"))
to_add <- rbind(to_add, c("900", "2/07/2017", "Z"))


demo1 <- rbind(demo1, to_add)




#### Adding multiple sequences
to_add <- data.frame(id = c("901", "901", "901", "901", "901", "901", "901", "901", "901", "901"),
                     period = c("01/01/2017", "01/02/2017", "01/03/2017", "01/04/2017", "01/05/2017", "02/04/2017", "02/05/2017", "02/06/2017", "02/07/2017", "02/08/2017"),
                     event = c("I", "K", "M", "Q", "Z", "B", "C", "D", "G", "M"))


to_add <- rbind(to_add, c("500", "01/01/2017", "A"))
to_add <- rbind(to_add, c("500", "01/02/2017", "B"))
to_add <- rbind(to_add, c("500", "01/03/2017", "C"))
to_add <- rbind(to_add, c("500", "03/04/2017", "B"))
to_add <- rbind(to_add, c("500", "03/05/2017", "I"))
to_add <- rbind(to_add, c("500", "03/06/2017", "J"))
to_add <- rbind(to_add, c("500", "04/07/2017", "X"))
to_add <- rbind(to_add, c("500", "06/01/2017", "B"))
to_add <- rbind(to_add, c("500", "06/02/2017", "K"))
to_add <- rbind(to_add, c("500", "06/03/2017", "L"))
to_add <- rbind(to_add, c("500", "06/03/2017", "M"))

to_add <- rbind(to_add, c("500", "02/01/2017", "B"))
to_add <- rbind(to_add, c("500", "02/02/2017", "B"))
to_add <- rbind(to_add, c("500", "02/03/2017", "B"))
to_add <- rbind(to_add, c("500", "04/04/2017", "B"))
to_add <- rbind(to_add, c("500", "04/05/2017", "B"))
to_add <- rbind(to_add, c("500", "04/06/2017", "B"))
to_add <- rbind(to_add, c("500", "04/07/2017", "C"))
to_add <- rbind(to_add, c("500", "07/01/2017", "C"))
to_add <- rbind(to_add, c("500", "07/02/2017", "C"))
to_add <- rbind(to_add, c("500", "07/03/2017", "C"))
to_add <- rbind(to_add, c("500", "07/03/2017", "C"))


to_add <- rbind(to_add, c("501", "01/01/2017", "A2"))
to_add <- rbind(to_add, c("501", "01/02/2017", "B4"))
to_add <- rbind(to_add, c("501", "01/03/2017", "C5"))
to_add <- rbind(to_add, c("501", "03/04/2017", "B5"))
to_add <- rbind(to_add, c("501", "03/05/2017", "I6"))
to_add <- rbind(to_add, c("501", "03/06/2017", "J6"))
to_add <- rbind(to_add, c("501", "04/07/2017", "X4"))
to_add <- rbind(to_add, c("501", "06/01/2017", "B4"))
to_add <- rbind(to_add, c("501", "06/02/2017", "K1"))
to_add <- rbind(to_add, c("501", "06/03/2017", "L2"))
to_add <- rbind(to_add, c("501", "06/03/2017", "M5"))


# Adding cases that are all numeric
to_add <- rbind(to_add, c("502", "01/01/2017", "2"))
to_add <- rbind(to_add, c("502", "01/02/2017", "4"))
to_add <- rbind(to_add, c("502", "01/03/2017", "5"))
to_add <- rbind(to_add, c("502", "03/04/2017", "5"))
to_add <- rbind(to_add, c("502", "03/05/2017", "6"))
to_add <- rbind(to_add, c("502", "03/06/2017", "6"))
to_add <- rbind(to_add, c("502", "04/07/2017", "4"))
to_add <- rbind(to_add, c("502", "06/01/2017", "4"))
to_add <- rbind(to_add, c("502", "06/02/2017", "1"))
to_add <- rbind(to_add, c("502", "06/03/2017", "2"))
to_add <- rbind(to_add, c("502", "06/03/2017", "5"))


# Adding cases that are a combination of numeric_alpha characters
to_add <- rbind(to_add, c("503", "01/01/2017", "2A"))
to_add <- rbind(to_add, c("503", "01/02/2017", "4B"))
to_add <- rbind(to_add, c("503", "01/03/2017", "5C"))
to_add <- rbind(to_add, c("503", "03/04/2017", "5B"))
to_add <- rbind(to_add, c("503", "03/05/2017", "6I"))
to_add <- rbind(to_add, c("503", "03/06/2017", "6J"))
to_add <- rbind(to_add, c("503", "04/07/2017", "4X"))
to_add <- rbind(to_add, c("503", "06/01/2017", "4B"))
to_add <- rbind(to_add, c("503", "06/02/2017", "1K"))
to_add <- rbind(to_add, c("503", "06/03/2017", "2L"))
to_add <- rbind(to_add, c("503", "06/03/2017", "5M"))


# Adding cases to make a second event set within consensus pattern
to_add <- rbind(to_add, c("340", "01/01/2017", "C"))
to_add <- rbind(to_add, c("340", "01/02/2017", "A"))
to_add <- rbind(to_add, c("340", "01/03/2017", "B"))
to_add <- rbind(to_add, c("340", "03/04/2017", "C"))
to_add <- rbind(to_add, c("340", "03/05/2017", "A"))
to_add <- rbind(to_add, c("340", "03/06/2017", "C"))
to_add <- rbind(to_add, c("340", "04/07/2017", "B"))
to_add <- rbind(to_add, c("340", "06/01/2017", "A"))
to_add <- rbind(to_add, c("340", "06/02/2017", "B"))
to_add <- rbind(to_add, c("340", "06/03/2017", "C"))
to_add <- rbind(to_add, c("340", "06/03/2017", "D"))

to_add <- rbind(to_add, c("340", "02/01/2017", "C"))
to_add <- rbind(to_add, c("340", "02/02/2017", "A"))
to_add <- rbind(to_add, c("340", "02/03/2017", "B"))
to_add <- rbind(to_add, c("340", "04/04/2017", "C"))
to_add <- rbind(to_add, c("340", "04/05/2017", "A"))
to_add <- rbind(to_add, c("340", "04/06/2017", "C"))
to_add <- rbind(to_add, c("340", "04/07/2017", "B"))
to_add <- rbind(to_add, c("340", "07/01/2017", "A"))
to_add <- rbind(to_add, c("340", "07/02/2017", "B"))
to_add <- rbind(to_add, c("340", "07/03/2017", "C"))
to_add <- rbind(to_add, c("340", "07/03/2017", "D"))

to_add <- rbind(to_add, c("440", "01/01/2017", "C"))
to_add <- rbind(to_add, c("440", "01/02/2017", "A"))
to_add <- rbind(to_add, c("440", "01/03/2017", "B"))
to_add <- rbind(to_add, c("440", "03/04/2017", "C"))
to_add <- rbind(to_add, c("440", "03/05/2017", "A"))
to_add <- rbind(to_add, c("440", "03/06/2017", "C"))
to_add <- rbind(to_add, c("440", "04/07/2017", "B"))
to_add <- rbind(to_add, c("440", "06/01/2017", "A"))
to_add <- rbind(to_add, c("440", "06/02/2017", "B"))
to_add <- rbind(to_add, c("440", "06/03/2017", "C"))
to_add <- rbind(to_add, c("440", "06/03/2017", "D"))

to_add <- rbind(to_add, c("440", "02/01/2017", "C"))
to_add <- rbind(to_add, c("440", "02/02/2017", "C"))
to_add <- rbind(to_add, c("440", "02/03/2017", "C"))
to_add <- rbind(to_add, c("440", "04/04/2017", "C"))
to_add <- rbind(to_add, c("440", "04/05/2017", "C"))
to_add <- rbind(to_add, c("440", "04/06/2017", "C"))
to_add <- rbind(to_add, c("440", "04/07/2017", "B"))
to_add <- rbind(to_add, c("440", "07/01/2017", "B"))
to_add <- rbind(to_add, c("440", "07/02/2017", "B"))
to_add <- rbind(to_add, c("440", "07/03/2017", "B"))
to_add <- rbind(to_add, c("440", "07/03/2017", "B"))

to_add <- rbind(to_add, c("540", "01/01/2017", "C"))
to_add <- rbind(to_add, c("540", "01/02/2017", "A"))
to_add <- rbind(to_add, c("540", "01/03/2017", "B"))
to_add <- rbind(to_add, c("540", "03/04/2017", "C"))
to_add <- rbind(to_add, c("540", "03/05/2017", "A"))
to_add <- rbind(to_add, c("540", "03/06/2017", "C"))
to_add <- rbind(to_add, c("540", "04/07/2017", "B"))
to_add <- rbind(to_add, c("540", "06/01/2017", "A"))
to_add <- rbind(to_add, c("540", "06/02/2017", "B"))
to_add <- rbind(to_add, c("540", "06/03/2017", "C"))
to_add <- rbind(to_add, c("540", "06/03/2017", "D"))

to_add <- rbind(to_add, c("540", "02/01/2017", "C"))
to_add <- rbind(to_add, c("540", "02/02/2017", "A"))
to_add <- rbind(to_add, c("540", "02/03/2017", "B"))
to_add <- rbind(to_add, c("540", "04/04/2017", "C"))
to_add <- rbind(to_add, c("540", "04/05/2017", "A"))
to_add <- rbind(to_add, c("540", "04/06/2017", "C"))
to_add <- rbind(to_add, c("540", "04/07/2017", "B"))
to_add <- rbind(to_add, c("540", "07/01/2017", "A"))
to_add <- rbind(to_add, c("540", "07/02/2017", "B"))
to_add <- rbind(to_add, c("540", "07/03/2017", "C"))
to_add <- rbind(to_add, c("540", "07/03/2017", "D"))


## Adding to test variable time frame grouping
to_add <- rbind(to_add, c("11151", "01/01/2017", "A"))
to_add <- rbind(to_add, c("11151", "01/02/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/03/2017", "C"))
to_add <- rbind(to_add, c("11151", "01/04/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/05/2017", "I"))
to_add <- rbind(to_add, c("11151", "01/06/2017", "J"))
to_add <- rbind(to_add, c("11151", "01/07/2017", "X"))
to_add <- rbind(to_add, c("11151", "01/08/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/09/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/10/2017", "C"))
to_add <- rbind(to_add, c("11151", "01/11/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/12/2017", "I"))
to_add <- rbind(to_add, c("11151", "01/13/2017", "J"))
to_add <- rbind(to_add, c("11151", "01/14/2017", "X"))
to_add <- rbind(to_add, c("11151", "01/15/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/16/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/17/2017", "C"))
to_add <- rbind(to_add, c("11151", "01/18/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/19/2017", "I"))
to_add <- rbind(to_add, c("11151", "01/20/2017", "J"))
to_add <- rbind(to_add, c("11151", "01/21/2017", "X"))
to_add <- rbind(to_add, c("11151", "01/22/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/23/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/24/2017", "C"))
to_add <- rbind(to_add, c("11151", "01/25/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/26/2017", "I"))
to_add <- rbind(to_add, c("11151", "01/27/2017", "J"))
to_add <- rbind(to_add, c("11151", "01/28/2017", "X"))
to_add <- rbind(to_add, c("11151", "01/29/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/30/2017", "B"))
to_add <- rbind(to_add, c("11151", "01/31/2017", "C"))
to_add <- rbind(to_add, c("11151", "02/02/2017", "B"))

to_add <- rbind(to_add, c("11151", "02/03/2017", "I"))
to_add <- rbind(to_add, c("11151", "02/10/2017", "J"))
to_add <- rbind(to_add, c("11151", "02/15/2017", "I"))
to_add <- rbind(to_add, c("11151", "02/20/2017", "J"))

to_add <- rbind(to_add, c("11151", "03/02/2017", "X"))
to_add <- rbind(to_add, c("11151", "03/07/2017", "B"))
to_add <- rbind(to_add, c("11151", "03/28/2017", "X"))
to_add <- rbind(to_add, c("11151", "03/29/2017", "B"))

to_add <- rbind(to_add, c("11151", "04/01/2017", "X"))
to_add <- rbind(to_add, c("11151", "04/07/2017", "B"))

to_add <- rbind(to_add, c("11151", "04/20/2017", "X"))
to_add <- rbind(to_add, c("11151", "04/21/2017", "B"))



to_add <- rbind(to_add, c("11152", "01/01/2017", "A"))
to_add <- rbind(to_add, c("11152", "01/02/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/03/2017", "C"))
to_add <- rbind(to_add, c("11152", "01/04/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/05/2017", "I"))
to_add <- rbind(to_add, c("11152", "01/06/2017", "J"))
to_add <- rbind(to_add, c("11152", "01/07/2017", "X"))
to_add <- rbind(to_add, c("11152", "01/08/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/09/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/10/2017", "C"))
to_add <- rbind(to_add, c("11152", "01/11/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/12/2017", "I"))
to_add <- rbind(to_add, c("11152", "01/13/2017", "J"))
to_add <- rbind(to_add, c("11152", "01/14/2017", "X"))
to_add <- rbind(to_add, c("11152", "01/15/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/16/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/17/2017", "C"))
to_add <- rbind(to_add, c("11152", "01/18/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/19/2017", "I"))
to_add <- rbind(to_add, c("11152", "01/20/2017", "J"))
to_add <- rbind(to_add, c("11152", "01/21/2017", "X"))
to_add <- rbind(to_add, c("11152", "01/22/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/23/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/24/2017", "C"))
to_add <- rbind(to_add, c("11152", "01/25/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/26/2017", "I"))
to_add <- rbind(to_add, c("11152", "01/27/2017", "J"))
to_add <- rbind(to_add, c("11152", "01/28/2017", "X"))
to_add <- rbind(to_add, c("11152", "01/29/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/30/2017", "B"))
to_add <- rbind(to_add, c("11152", "01/31/2017", "C"))
to_add <- rbind(to_add, c("11152", "02/02/2017", "B"))
to_add <- rbind(to_add, c("11152", "02/15/2017", "I"))
to_add <- rbind(to_add, c("11152", "02/20/2017", "J"))
to_add <- rbind(to_add, c("11152", "03/28/2017", "X"))
to_add <- rbind(to_add, c("11152", "03/29/2017", "B"))
to_add <- rbind(to_add, c("11152", "04/20/2017", "X"))
to_add <- rbind(to_add, c("11152", "04/21/2017", "B"))


to_add <- rbind(to_add, c("11153", "01/01/2017", "A"))
to_add <- rbind(to_add, c("11153", "01/02/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/03/2017", "C"))
to_add <- rbind(to_add, c("11153", "01/04/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/05/2017", "I"))
to_add <- rbind(to_add, c("11153", "01/06/2017", "J"))
to_add <- rbind(to_add, c("11153", "01/07/2017", "X"))
to_add <- rbind(to_add, c("11153", "01/08/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/09/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/10/2017", "C"))
to_add <- rbind(to_add, c("11153", "01/11/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/12/2017", "I"))
to_add <- rbind(to_add, c("11153", "01/13/2017", "J"))
to_add <- rbind(to_add, c("11153", "01/14/2017", "X"))
to_add <- rbind(to_add, c("11153", "01/15/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/16/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/17/2017", "C"))
to_add <- rbind(to_add, c("11153", "01/18/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/19/2017", "I"))
to_add <- rbind(to_add, c("11153", "01/20/2017", "J"))
to_add <- rbind(to_add, c("11153", "01/21/2017", "X"))
to_add <- rbind(to_add, c("11153", "01/22/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/23/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/24/2017", "C"))
to_add <- rbind(to_add, c("11153", "01/25/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/26/2017", "I"))
to_add <- rbind(to_add, c("11153", "01/27/2017", "J"))
to_add <- rbind(to_add, c("11153", "01/28/2017", "X"))
to_add <- rbind(to_add, c("11153", "01/29/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/30/2017", "B"))
to_add <- rbind(to_add, c("11153", "01/31/2017", "C"))
to_add <- rbind(to_add, c("11153", "02/02/2017", "B"))
to_add <- rbind(to_add, c("11153", "02/15/2017", "I"))
to_add <- rbind(to_add, c("11153", "02/20/2017", "J"))
to_add <- rbind(to_add, c("11153", "03/28/2017", "X"))
to_add <- rbind(to_add, c("11153", "03/29/2017", "B"))
to_add <- rbind(to_add, c("11153", "04/20/2017", "X"))
to_add <- rbind(to_add, c("11153", "04/21/2017", "B"))


to_add <- rbind(to_add, c("11154", "01/01/2017", "A"))
to_add <- rbind(to_add, c("11154", "01/02/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/03/2017", "C"))
to_add <- rbind(to_add, c("11154", "01/04/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/05/2017", "I"))
to_add <- rbind(to_add, c("11154", "01/06/2017", "J"))
to_add <- rbind(to_add, c("11154", "01/07/2017", "X"))
to_add <- rbind(to_add, c("11154", "01/08/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/09/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/10/2017", "C"))
to_add <- rbind(to_add, c("11154", "01/11/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/12/2017", "I"))
to_add <- rbind(to_add, c("11154", "01/13/2017", "J"))
to_add <- rbind(to_add, c("11154", "01/14/2017", "X"))
to_add <- rbind(to_add, c("11154", "01/15/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/16/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/17/2017", "C"))
to_add <- rbind(to_add, c("11154", "01/18/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/19/2017", "I"))
to_add <- rbind(to_add, c("11154", "01/20/2017", "J"))
to_add <- rbind(to_add, c("11154", "01/21/2017", "X"))
to_add <- rbind(to_add, c("11154", "01/22/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/23/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/24/2017", "C"))
to_add <- rbind(to_add, c("11154", "01/25/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/26/2017", "I"))
to_add <- rbind(to_add, c("11154", "01/27/2017", "J"))
to_add <- rbind(to_add, c("11154", "01/28/2017", "X"))
to_add <- rbind(to_add, c("11154", "01/29/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/30/2017", "B"))
to_add <- rbind(to_add, c("11154", "01/31/2017", "C"))
to_add <- rbind(to_add, c("11154", "02/02/2017", "B"))
to_add <- rbind(to_add, c("11154", "02/15/2017", "I"))
to_add <- rbind(to_add, c("11154", "02/20/2017", "J"))
to_add <- rbind(to_add, c("11154", "03/28/2017", "X"))
to_add <- rbind(to_add, c("11154", "03/29/2017", "B"))
to_add <- rbind(to_add, c("11154", "04/20/2017", "X"))
to_add <- rbind(to_add, c("11154", "04/21/2017", "B"))


















