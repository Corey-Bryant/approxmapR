#' @export
format_sequence.W_Sequence_Dataframe <-
    function(df,
             compare = FALSE,
             truncate_patterns = FALSE,
             html_format = FALSE) {
        column_patterns <- names(df)[str_detect(names(df), "_pattern")]
        if ("weighted_sequence" %in% names(df)) {
            columns <- c(column_patterns, "weighted_sequence")
        } else {
            columns <- column_patterns
        }

        if (truncate_patterns) {
            df <- df %>% mutate_at(column_patterns, truncate_pattern)
        }



        df <-
            df %>% select(one_of("cluster", "n", columns)) %>% mutate_at(columns, function(x)
                format_sequence(x,
                                html_format = html_format)) %>% mutate(n = as.double(n), n_percent = str_c(round(n /
                                                                                                                     sum(n) * 100,
                                                                                                                 digits = 2), "%")) %>% select(one_of("cluster", "n", "n_percent", columns))

        if (compare) {
            compare_sequences(df)
        } else {
            df
        }

    }

#' @export
view_formatted_sequence <- function(seq) {
    format_sequence(seq, html = TRUE) %>% stringr:::str_view_widget()
}

#' @export
compare_sequences <- function(df) {
    df %>% gather(-cluster,-n,-n_percent, key = "pattern", value = "sequence") %>% arrange(cluster) %>%
        mutate(pattern = stringr::str_replace(pattern, "_pattern", ""))
}




#' @export
class_it <- function(obj, class_name) {
    class(obj) <- c(class_name, class(obj)) %>% unique()
    obj
}

#' @export
truncate_pattern <- function(x, ...) {
    UseMethod("truncate_pattern")
}

#' @export
truncate_pattern.W_Sequence_Pattern_List <-
    function(w_sequence_list) {

        class_it(map(w_sequence_list, truncate_pattern),
                 "W_Sequence_List")
    }

## [ISSUE HERE]
#' @export
truncate_pattern.W_Sequence_Pattern <- function(w_sequence) {

    #browser()


    truncate_index <- rep(FALSE, length(w_sequence))
    for (i in seq_along(w_sequence)) {
        if (i == length(w_sequence))
            (break)()
        e_1 <- sort(w_sequence[[i]]$elements)
        e_2 <- sort(w_sequence[[i + 1]]$elements)

        # Commented out on 08/06/2021 - seems to fix the removal of duplicated event sets
        #   witin the _truncated view. Good thing, we want all to be shown.

        #if (identical(e_1, e_2)) {
        #    truncate_index[i] <- TRUE
        #}
    }
    w_sequence[truncate_index] <- NULL

    ## May need to uncomment this out
    #compressed_n <-
    #    (truncate_index %>%
    #         as.integer() %>%
    #         as.character() %>%
    #         str_c(collapse = "") %>%
    #         str_split("0") %>%
    #         pluck(1) %>%
    #         str_subset(".") %>%
    #         str_count("1")) + 1

    #for(i in seq(1,length(w_sequence))){
    #    w_sequence[[i]]$itemset_weight <- compressed_n[i]
    #}


    w_sequence %>% class_it("W_Sequence_Pattern_Compressed")
}

format_sequence.W_Sequence_Pattern_Compressed <-
    function(w_sequence, html_format = FALSE) {
        n <- attr(w_sequence, "n")
        if (html_format) {
            if(n > 1){
                colors <-
                    rev(colormap::colormap(colormap = "bluered", nshades = n) %>%
                            stringr::str_sub(1, -3))
            } else {
                colors <- colormap::colormap(nshades = 2)[1]
            }


            # cuts <- floor(n*seq(0,1,0.2))[2:5]
            w_sequence %>%
                map_chr(function(w_itemset) {
                    tibble(
                        element = as.character(w_itemset$elements),
                        weight = as.integer(w_itemset$element_weights)
                    ) %>%
                        mutate(
                            ratio = weight / n,
                            # i = ceiling(ratio),
                            # color = map_chr(i, ~colors[.]),
                            color = colors[weight],
                            font_size = paste0(floor((1 + ratio * .6) * 100), "%"),
                            font_weight = signif(460 + ratio * 340, 1),
                            otag = str_c(
                                '<span style="',
                                "color: ",
                                color,
                                "; ",
                                "font-size: ",
                                font_size,
                                "; ",
                                "font-weight: ",
                                font_weight,
                                ";",
                                '">'
                            ),
                            ctag = "</span>",
                            element_html = str_c(otag, element, ":", weight, ctag)
                        ) %>%
                        pull(element_html) %>%
                        str_c(collapse = ", ") %>%
                        paste0("(", ., ")", ":", w_itemset$itemset_weight, "<br>")
                }) %>%
                str_c(collapse = " ") %>%
                paste0("<", ., ">", " : ", n) %>%
                stringr::str_replace("<\\(", " < ( ")

        } else{
            w_sequence %>%
                map_chr(function(w_itemset) {
                    if (length(w_itemset$elements) > 0) {
                        str_c(w_itemset$elements, ":", w_itemset$element_weights) %>%
                            str_c(collapse = ", ") %>%
                            paste0("(", ., ")", ":", w_itemset$itemset_weight)
                    }
                    else{
                        NA
                    }

                }) %>%
                .[!is.na(.)] %>%
                str_c(collapse = " ") %>%
                paste0("<", ., "<br>", ">", " : ", n)
        }

    }


#' @export
w_sequence_to_tibble <- function(w_sequence) {
    tibble(
        element = map(w_sequence, "elements") %>% unlist(),
        element_weight = map(w_sequence, "element_weights") %>%
            unlist(),
        itemset = map2(1:length(w_sequence), w_sequence, ~ rep(.x, length(.y$elements))) %>% unlist()
    ) %>%
        mutate(element_no = row_number())
}

#' @export
plot_weighted_sequence <- function(w_sequence) {
    df_sequence <- w_sequence %>% w_sequence_to_tibble()

    df_itemset <-
        df_sequence %>% group_by(itemset) %>% filter(element_no == max(element_no))

    df_sequence %>% ggplot(aes(element_no, element_weight)) + geom_point() + geom_label(aes(y = element_weight +
                                                                                                0.02 * element_weight, label = element)) + geom_vline(data = df_itemset, aes(xintercept = element_no))
}





#' @export
convert_to_events <- function(data, id_column, sequence_column) {

  data %>%
    mutate(event_set = str_split(data[[sequence_column]], "[\\(\\)]")) %>%
    unnest(cols = c(event_set)) %>%
    filter(event_set != "") %>% filter(event_set != " ") %>%
    group_by({{ sequence_column }}) %>%
    mutate(period = row_number()) %>%
    mutate(event = str_split(event_set, "[, ]")) %>%
    unnest(cols = c(event)) %>%
    filter(event != "") %>% filter(event != " ") %>% ungroup() %>%
    select(id_column, period, event)

}

#' @export
convert_consensus_to_events <- function(df_patterns, consensus_pattern_column = "consensus_pattern") {

  list_cluster <- numeric(1)
  list_period <- numeric(1)
  list_elements <- character(1)

list_cluster <- numeric(1)
list_period <- numeric(1)
list_elements <- character(1)


ix <- 0
for (c in df_patterns$cluster) {

    current_concensus <- subset(df_patterns, cluster == c)[[consensus_pattern_column]][[1]]
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



#' @export
consensus_sequence_distance <- function(df_aggregated, consensus_sequence_list) {

  total_id_sequences = n_distinct(df_aggregated$id)

  list_cluster <- numeric(total_id_sequences)
  list_id <- character(total_id_sequences)
  list_consensusid <- character(total_id_sequences)
  list_distance <- character(total_id_sequences)


  # For each individual sequence, calculate the distance from each consensus pattern
  #   and assign
  df_cluster <- df_aggregated %>% convert_to_sequence() %>% ungroup()
  consensus_sequence_list <- consensus_sequence_list %>% convert_to_sequence() %>% ungroup()

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
  id_consensus_distance$distance <- as.numeric(id_consensus_distance$distance)

  id_consensus_distance <- id_consensus_distance %>% filter(distance == minDistance)

  # Calculate the current average distance for each group
  group_averages <- id_consensus_distance %>% filter(assignedConsensusCluster != "NEED SENSITIVITY ANALYSIS") %>%
                    group_by(consensusid) %>% summarise(average_distance = mean(as.numeric(distance), na.rm = TRUE),
                                                        nID = n_distinct(id)) %>% ungroup()

  id_consensus_distance <- left_join(id_consensus_distance, group_averages, by = "consensusid")

  id_consensus_distance_dups <- id_consensus_distance %>% filter (nDup > 1)
  id_consensus_distance_dups <- id_consensus_distance_dups %>% mutate(average_distance_with = (distance + average_distance) / 2,
                                                                      average_distance_with_delta = average_distance_with - average_distance,
                                                                      nID_with = nID + 1) %>%
                                group_by(id) %>% mutate(min_delta = min(abs(average_distance_with_delta)),
                                                        assignedConsensusCluster_new = case_when(
                                                          abs(average_distance_with_delta) == min_delta ~ consensusid,
                                                          TRUE ~ "NEED SENSITIVITY ANALYSIS")) %>%
                                ungroup() %>% filter(assignedConsensusCluster_new != "NEED SENSITIVITY ANALYSIS") %>%
                                select("id", "assignedConsensusCluster_new")

  id_consensus_distance <- left_join(id_consensus_distance, id_consensus_distance_dups, by = "id")
  id_consensus_distance <- id_consensus_distance %>% mutate(assignedConsensusCluster = case_when(
                                                                                          assignedConsensusCluster == "NEED SENSITIVITY ANALYSIS" ~ assignedConsensusCluster_new,
                                                                                          TRUE ~ assignedConsensusCluster)
                                                            ) %>% filter(consensusid == assignedConsensusCluster) %>%
                                                                  select("id", "assignedConsensusCluster", "distance")

  return(id_consensus_distance)

}





# -sequencer- is a slightly modified version of -format_sequence- in that it adds
#   a comma between event sets in a sequence for the id
#' @export
sequencer <- function(sequence) {

  sequence <- sequence %>% map_chr(function(itemset) {
                              itemset <- str_c(itemset, collapse = ", ")
                              paste0("(", itemset, ")")
                           }) %>%
                           str_c(collapse = ", ")

  sequence <- paste0("<", sequence, ">")

  as.character(sequence)

}



#' @export
pattern_search <- function(Clustered_Dataframe, find_pattern = NULL, event_set = FALSE, exact = FALSE) {

  ## Checking parameters and criteria - checks verified ##
  if (is.null(find_pattern)) {
    stop("Error: find_pattern parameter is NULL.")
  }

  if (event_set & exact){
    stop("Error: The event_set and exact parameters both cannot be TRUE")
  }

  if ("Clustered_Dataframe" %in% class(Clustered_Dataframe)) {
    # This is code to find the pattern for the clustered dataframe. This is
    #   class is produced during the clustering step and/or after filter_pattern
    #   which finds the consensus patterns.
    df_seq <- Clustered_Dataframe %>%
      select(cluster, n, df_sequences) %>%
      unnest(cols = c(df_sequences))
    df_seq <- df_seq %>% mutate(sequences = map_chr(sequence, format_sequence))
  }

  if ("Aggregated_Dataframe" %in% class(Clustered_Dataframe)) {
    # This is code to find the pattern for the clustered dataframe. This is
    #   class is produced during the clustering step and/or after filter_pattern
    #   which finds the consensus patterns.
    df_seq <- Clustered_Dataframe %>% convert_to_sequence()
    names(df_seq) <- c("id", "sequence", "sequences")
  }




  if (event_set) {

    find_pattern <- str_replace_all(find_pattern, fixed("("), "\\(")
    find_pattern <- str_replace_all(find_pattern, fixed(")"), "\\)")


    # Match an event 0 or more times
    find_pattern <- str_replace_all(find_pattern, fixed("event*, "), "([:alnum:]*, )*")
    find_pattern <- str_replace_all(find_pattern, fixed(", event*"), "(, [:alnum:]*)*")


    # Match an event 1 or more times
    find_pattern <- str_replace_all(find_pattern, fixed("event+, "), "([:alnum:]*, )+")
    find_pattern <- str_replace_all(find_pattern, fixed(", event+"), "(, [:alnum:]*)+")


    # Wild card - any alphanumeric ([:alnum:]), punction ([:punct:]), and space characters
    find_pattern <- str_replace_all(find_pattern, fixed("**"), "[[:print:]]*")


    # Match an event set structure 0 or more times
    find_pattern <- str_replace_all(find_pattern, fixed("eventset* "), "(\\([[:alnum:], ]*[[:alnum:]*]+\\) )*")
    find_pattern <- str_replace_all(find_pattern, fixed(" eventset*"), "( \\([[:alnum:], ]*[[:alnum:]*]+\\))*")

    # Match an event set structure 1 or more times
    find_pattern <- str_replace_all(find_pattern, fixed("eventset+ "), "(\\([[:alnum:]*, ]*[[:alnum:]*]+\\) )+")
    find_pattern <- str_replace_all(find_pattern, fixed(" eventset+"), "( \\([[:alnum:]*, ]*[[:alnum:]*]+\\))+")


  } else if (exact) {

    find_pattern <- fixed(find_pattern)

  } else {

    pieces <- (str_extract_all(find_pattern, "\\(|(([:alnum:]*)[:alnum:](?=,|\\)))|\\)| "))[[1]]

    pieces_conv <- str_replace_all(pieces, "\\(", "(?:") %>% str_replace_all(., "\\)", ")")

    pieces <- str_subset(pieces, "[^ ]")


    sets <- str_c(pieces_conv, collapse= "")
    sets <- str_replace_all(sets, "(?<!\\)) ", "|")
    sets <- str_split(sets, " ")[[1]]






    # Building pattern structure
    pattern <- ""
    previous_item <- ""
    item_index <- 1
    end <- length(pieces)

    sets_counter <- 1

    for (item in pieces) {

      if (item == "(" & item_index == 1) {

        pattern <- str_c(pattern, "[\\(([:alnum:], )*([:alnum:])+\\) ]*", "\\(")

      } else if (item == "(" & item_index > 1) {

        pattern <- str_c(pattern, " \\(")

      } else if (item == ")" & item_index != end) {

        pattern <- str_c(pattern, "(, [:alnum:]*)*", "\\)", "[ \\(([:alnum:], )*([:alnum:])+\\)]*")

        sets_counter <- sets_counter + 1

      } else if (item == ")" & item_index == end) {

        pattern <- str_c(pattern, "(, [:alnum:]*)*", "\\)", "[ \\(([:alnum:], )*([:alnum:])+\\)]*")

        sets_counter <- sets_counter + 1

      } else {

        if (pieces[item_index + 1] == ")") {

          pattern <- str_c(pattern,  "([:alnum:]*, )*", sets[sets_counter])

        } else {

          pattern <- str_c(pattern,  "([:alnum:]*, )*", sets[sets_counter], ", ")

        }

      }

      item_index <- item_index + 1

    }

    find_pattern <- pattern

  }

  #print(find_pattern)


  # Now to pull the IDs with the pattern(s)
  # print(find_pattern)
  if (length(find_pattern) > 1) {
    count <- 1
    for (pattern in find_pattern) {
      #print(pattern)
      if (count == 1) {
        to_pull <- str_detect(df_seq$sequences, pattern)
        count = count + 1
      } else {
        to_pull_n <- str_detect(df_seq$sequences, pattern)
        count = count + 1
        to_pull <- replace(to_pull, to_pull_n, TRUE)
      }
    }
  } else {
    pattern <- find_pattern
    to_pull <- str_detect(df_seq$sequences, pattern)
  }

  df_seq <- subset.data.frame(df_seq, subset = to_pull)

  if ("Clustered_Dataframe" %in% class(Clustered_Dataframe)) {
    df_seq <- df_seq %>% select(cluster, id, sequence, sequences) %>% arrange(cluster)
  }

  if ("Aggregated_Dataframe" %in% class(Clustered_Dataframe)) {
    df_seq <- df_seq %>% select(id, sequence, sequences) %>% arrange(id)
  }

  return(df_seq)

}





plot_ktable <- function(ktable,
                        validation_measure = 'silhouette',
                        save_graph = TRUE,
                        graph_file_name = NULL,
                        size_width = 855, size_height = 317,
                        output_directory = "~") {


  # Parameter Checks
  stopifnot("ktable" %in% class(ktable))


  if (save_graph) {

    if (is.null(graph_file_name)) {

      graph_file_name = paste0(attr(ktable, "algorithm"), " Optimal K Plot_", validation_measure, ".png")

    }

    if (!endsWith(graph_file_name, ".png")) {
      stop("Error: The graph file name must end with '.png'. Only PNG images are supported at this time.")
    }


  }



  # Plotting Information
  if (validation_measure == 'silhouette') {

    measure = 'Average Silhouette Width'
    measure_values =ktable$average_silhouette_width

  } else if (validation_measure == 'dunn') {

    measure = 'Dunn Index'
    measure_values = ktable$dunn

  } else if (validation_measure == 'wb_ratio') {

    measure = 'Average Distance Within Cluster / Average Distance Between Clusters'
    measure_values = ktable$wb_ratio

  } else if (validation_measure == 'average_between') {

    measure = 'Average Distance Between Clusters'
    measure_values = ktable$average_between

  } else if (validation_measure == 'average_within') {

    measure = 'Average Distance Within Cluster'
    measure_values = ktable$average_within

  } else if (validation_measure == 'within_cluster_ss') {

    measure = 'Sum of Within Cluster / Cluster Size'
    measure_values = ktable$within_cluster_ss

  } else {

    stop("Only validation measures of silhouette, dunn, wb_ratio, average_between, average_within, and within_cluster_ss are supported.")

  }



  if (validation_measure == 'silhouette') {

    k_plot <- ggplot(ktable, aes(k, average_silhouette_width)) +

      geom_line(color = "#20B2AA") +

      geom_errorbar(aes(ymax = average_silhouette_width_upper_ci,
                        ymin = average_silhouette_width_lower_ci),
                    width = .25,
                    color = "#20B2AA") +

      geom_vline(xintercept = ktable$k[which.max(ktable$average_silhouette_width)],
                 color = "#20B2AA", linetype = 'dashed') +

      labs(title = paste0(attr(ktable, "algorithm"), " Optimal K Plot"),
           subtitle = paste0("k =", ktable$k[which.max(ktable$average_silhouette_width)], "; Max average silhouette width = ", round(ktable$average_silhouette_width[which.max(ktable$average_silhouette_width)], digits = 3)),
           x = "K Value",
           y = measure)  +


      coord_cartesian(xlim = c(min(ktable$k), max(ktable$k))) +

      scale_x_continuous(labels = as.character(ktable$k), breaks = ktable$k)



  } else if (validation_measure == 'dunn'){

    k_plot <- ggplot(ktable, aes(k, dunn)) +

      geom_line(color = "#20B2AA") +

      geom_vline(xintercept = ktable$k[which.max(ktable$dunn)],
                 color = "#20B2AA", linetype = 'dashed') +

      labs(title = paste0(attr(ktable, "algorithm"), " Optimal K Plot"),
           subtitle = paste0("k =", ktable$k[which.max(ktable$dunn)], "; Max Dunn Index = ", round(ktable$dunn[which.max(ktable$dunn)], digits = 3)),
           x = "K Value",
           y = measure)  +


      coord_cartesian(xlim = c(min(ktable$k), max(ktable$k))) +

      scale_x_continuous(labels = as.character(ktable$k), breaks = ktable$k)



  } else {


    k_plot <- ggplot(ktable, aes(k, measure_values)) +

      geom_line(color = "#20B2AA") +

      labs(title = paste0(attr(ktable, "algorithm"), " Optimal K Plot"),
           x = "K Value",
           y = measure) +

      coord_cartesian(xlim = c(min(ktable$k), max(ktable$k))) +

      scale_x_continuous(labels = as.character(ktable$k), breaks = ktable$k)


  }





  # This portion saves the graph if the option is selected
  if (save_graph) {
    output_directory <- create_folder(output_directory, "approxmap_results")
    output_directory_graphs <- create_folder(output_directory, "graphs")

    png(file = paste0(output_directory_graphs, "/", file_check(output_directory_graphs,
                                                               graph_file_name)), width = size_width, height = size_height)
    print(k_plot)
    dev.off()


  }

  print(k_plot)


}





## plot_silhouette
plot_silhouette <- function(ktable,
                            save_graph = TRUE,
                            graph_file_name_individual = NULL,
                            graph_file_name_cluster = NULL,
                            size_width = 855, size_height = 317,
                            save_table = TRUE,
                            table_file_name = NULL,
                            output_directory = "~") {


  # Calculating average silhouette width per cluster and merging into org. data
  ktable$cluster <- as.factor(ktable$cluster)

  clustertable <- ktable %>% group_by(cluster) %>% summarize(cluster_sil_width = mean(sil_width))

  ktable <- merge(ktable, clustertable, on = c("cluster")) %>% select(id, cluster, neighbor, sil_width, cluster_sil_width)



  # Creating plot at the individual level
  id_plot <- ggplot(ktable, aes(x = sil_width, y = reorder(id, cluster, sort),
                                group = cluster, label = round(sil_width, digits = 2))) +

    geom_col(aes(fill = cluster), color = "white", position = "dodge") +
    geom_text(hjust = -0.2) +

    scale_fill_hue(l = 40) +

    labs(title = "Individual Silhoutte Plot",
         subtitle = paste0("n = ", length(ktable$id), "; Average silhouette width = ", round(mean(ktable$sil_width), digits = 3)),
         x = "Silhoutte Width S_i",
         y = "ID") +

    coord_cartesian(xlim = c(min(ktable$sil_width), 1)) +
    guides(fill = guide_legend(title = "Cluster"))



  # Creating plot at the cluster level
  cluster_plot <- ggplot(clustertable, aes(x = cluster_sil_width, y = reorder(cluster, cluster, sort),
                                           label = round(cluster_sil_width, digits = 2))) +

    geom_col(aes(fill = cluster), color = "white", position = "dodge") +
    geom_text(hjust = -0.2) +

    scale_fill_hue(l = 40) +

    labs(title = "Cluster Silhoutte Plot",
         subtitle = paste0("n = ", length(ktable$id), "; Average silhouette width = ", round(mean(ktable$sil_width), digits = 3)),
         x = "Average Silhoutte Width per Cluster",
         y = "Cluster") +

    coord_cartesian(xlim = c(min(ktable$sil_width), 1)) +
    guides(fill = guide_legend(title = "Cluster"))





  # Writing the table and graphs if desired
  if (save_graph) {
    output_directory <- create_folder(output_directory, "approxmap_results")
    output_directory_graphs <- create_folder(output_directory, "graphs")

    # Checking the individual graph information
    if (is.null(graph_file_name_individual)) {

      graph_file_name_individual = paste0("Individual Silhouette Plot.png")

    }

    if (!endsWith(graph_file_name_individual, ".png")) {
      stop("Error: The graph file name at the individual level must end with '.png'. Only PNG images are supported at this time.")
    }


    # Checking the cluster graph information
    if (is.null(graph_file_name_cluster)) {

      graph_file_name_cluster = paste0("Cluster Silhouette Plot.png")

    }

    if (!endsWith(graph_file_name_cluster, ".png")) {
      stop("Error: The graph file name at the cluster level must end with '.png'. Only PNG images are supported at this time.")
    }


    # Writing the graphs
    png(file = paste0(output_directory_graphs, "/", file_check(output_directory_graphs,
                                                               graph_file_name_individual)), width = size_width, height = size_height)
    print(id_plot)
    dev.off()


    png(file = paste0(output_directory_graphs, "/", file_check(output_directory_graphs,
                                                               graph_file_name_cluster)), width = size_width, height = size_height)
    print(cluster_plot)
    dev.off()
  }


  if (save_table) {
    #output_directory <- create_folder(output_directory, "approxmap_results")
    output_directory_table <- create_folder(output_directory, "private")


    if (is.null(table_file_name)) {

      table_file_name = "Silhouette Clustering Information.csv"

    }

    write.csv(ktable,
              file = paste0(output_directory_table, "/", file_check(output_directory_table, table_file_name)),
              row.names = FALSE)

  }

  list(print(id_plot), print(cluster_plot))


}
