#' Search genes from gff files
#'
#' @description The goal of search_genes is to append genes to combined_gwas_result.
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @param output_path An output path.
#' @return combined_gwas_result or NULL if something missing.
#' @keywords Genes
#' @export
#'
search_genes <- function(combined_gwas_result = NULL,
                            output_path = file.path("~"),
                            GFF_file_path = NULL, 
                            GFF_file_name = NULL, 
                            GFF_file_extension = NULL, 
                            GFF_file_named_sequentially_from = NULL, 
                            GFF_file_named_sequentially_to = NULL) {

    #######################################################################
    ## Initialize variables
    #######################################################################

    Count = 0

    #######################################################################
    ## Create folders to store outputs
    #######################################################################

    auto_save_path <- file.path(output_path, "GAPIT_auto_output")
    GAPIT_manhattan_plot_save_path <- file.path(output_path, "GAPIT_Manhattan_Plot")
    GAPIT_qq_plot_save_path <- file.path(output_path, "GAPIT_QQ_Plot")
    GAPIT_significant_save_path <- file.path(output_path, "GAPIT_significant")
    ped_and_info_save_path <- file.path(output_path, "Haploview_PEDandINFO")
    ld_data_save_path <- file.path(output_path, "Haploview_LD_data")
    ld_plot_save_path <- file.path(output_path, "Haploview_LD_plot")
    haplotypes_gabriel_blocks_save_path <- file.path(output_path, "Haploview_Haplotypes_gabriel_blocks")
    gff_save_path <- file.path(output_path, "GFF")

    temp <- c(auto_save_path, GAPIT_manhattan_plot_save_path, GAPIT_qq_plot_save_path, GAPIT_significant_save_path, 
    ped_and_info_save_path, ld_data_save_path, ld_plot_save_path, haplotypes_gabriel_blocks_save_path, gff_save_path)

    for (i in 1:length(temp)) {
        if (!dir.exists(temp[i])){
            try(dir.create(temp[i]))
        }
    }

    #######################################################################
    ## Create GFF reference file path and put everything in a table
    #######################################################################

    GFF_file_table <- tryCatch({
        data.frame(
            "Chromosome" = GFF_file_named_sequentially_from:GFF_file_named_sequentially_to, 
            "File_path" = file.path(GFF_file_path, 
                                    paste0(GFF_file_name, 
                                            GFF_file_named_sequentially_from:GFF_file_named_sequentially_to, 
                                            ".", 
                                            GFF_file_extension
                                            )
                                    )
        )
    }, error = function(e) {
        print("GFF_file_table cannot be created!!!")
        return(NULL)
    })

    if(nrow(GFF_file_table) > 0){
        for(i in 1:nrow(GFF_file_table)){
            if(!file.exists(file.path(GFF_file_table[i, 2]))){
                print("File path for GFF does not exists.")
                return(NULL)
            }
        }
    } else{
        print("GFF_file_table has zero row.")
        return(NULL)
    }

    #######################################################################
    ## Search genes using GFF files
    #######################################################################

    combined_gwas_result <- as.data.frame(combined_gwas_result)
    combined_gwas_result$Gene_name <- NA
    combined_gwas_result$Gene_start <- NA
    combined_gwas_result$Gene_stop <- NA
    combined_gwas_result$Gene_description <- NA

    index <- match(as.integer(combined_gwas_result[1,2]), as.integer(GFF_file_table[,1]))
    index <- ifelse(is.na(index), 1, index)

    gff <- read_file(file_path = file.path(GFF_file_table[index, 2]), header = FALSE)
    if(is.null(gff)){
        print("Unable to read GFF reference file or combined GWAS results are all NA.")

    } else{
        # For each row in each table
        i <- 1
        while(i <= nrow(combined_gwas_result)){

            if(nrow(gff) > 0 & !is.na(as.integer(combined_gwas_result[i,2])) & !is.na(as.numeric(combined_gwas_result[i,3])) &
                !is.na(combined_gwas_result$LD_start[i]) & !is.na(combined_gwas_result$LD_end[i])){

                # If chromosome number is different, read new gff that matches the chromosome number
                if(as.integer(gff[1,1]) != as.integer(combined_gwas_result[i,2])){

                    index <- match(as.integer(combined_gwas_result[i,2]), as.integer(GFF_file_table[,1]))
                    index <- ifelse(is.na(index), 1, index)
                    gff <- read_file(file_path = file.path(GFF_file_table[index, 2]), header = FALSE)
                    if(is.null(gff)){
                        print("Unable to read GFF reference file.")
                        return(NULL)
                    }
                }

                if(all(c("LD_start", "LD_end", "Haploblock_start", "Haploblock_stop") %in% colnames(combined_gwas_result)) & as.integer(gff[1,1]) == as.integer(combined_gwas_result[i,2])){
                    if(nrow(gff) > 0 & !is.na(combined_gwas_result$Haploblock_start[i]) & !is.na(combined_gwas_result$Haploblock_stop[i] & as.integer(gff[1,1]) == as.integer(combined_gwas_result[i,2]))){

                        temp_gff <- gff[(
                        gff[,4] <= combined_gwas_result$Haploblock_start[i] &
                        gff[,4] <= combined_gwas_result$Haploblock_stop[i] &
                        gff[,5] >= combined_gwas_result$Haploblock_start[i] &
                        gff[,5] <= combined_gwas_result$Haploblock_stop[i]
                        ) | (
                        gff[,4] >= combined_gwas_result$Haploblock_start[i] &
                        gff[,4] <= combined_gwas_result$Haploblock_stop[i] &
                        gff[,5] >= combined_gwas_result$Haploblock_start[i] &
                        gff[,5] >= combined_gwas_result$Haploblock_stop[i]
                        ) | (
                        gff[,4] >= combined_gwas_result$Haploblock_start[i] &
                        gff[,4] <= combined_gwas_result$Haploblock_stop[i] &
                        gff[,5] >= combined_gwas_result$Haploblock_start[i] &
                        gff[,5] <= combined_gwas_result$Haploblock_stop[i]
                        ) | (
                        gff[,4] <= combined_gwas_result$Haploblock_start[i] &
                        gff[,4] <= combined_gwas_result$Haploblock_stop[i] &
                        gff[,5] >= combined_gwas_result$Haploblock_start[i] &
                        gff[,5] >= combined_gwas_result$Haploblock_stop[i]
                        ),]

                    } else if(nrow(gff) > 0 & !is.na(combined_gwas_result$LD_start[i]) & !is.na(combined_gwas_result$LD_end[i] & as.integer(gff[1,1]) == as.integer(combined_gwas_result[i,2]))){

                        temp_gff <- gff[(
                        gff[,4] <= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] <= combined_gwas_result$LD_end[i]
                        ) | (
                        gff[,4] >= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] >= combined_gwas_result$LD_end[i]
                        ) | (
                        gff[,4] >= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] <= combined_gwas_result$LD_end[i]
                        ) | (
                        gff[,4] <= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] >= combined_gwas_result$LD_end[i]
                        ),]

                    }
                } else if(all(c("LD_start", "LD_end") %in% colnames(combined_gwas_result)) & as.integer(gff[1,1]) == as.integer(combined_gwas_result[i,2])){
                    if(nrow(gff) > 0 & !is.na(combined_gwas_result$LD_start[i]) & !is.na(combined_gwas_result$LD_end[i]) & as.integer(gff[1,1]) == as.integer(combined_gwas_result[i,2])){

                        temp_gff <- gff[(
                        gff[,4] <= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] <= combined_gwas_result$LD_end[i]
                        ) | (
                        gff[,4] >= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] >= combined_gwas_result$LD_end[i]
                        ) | (
                        gff[,4] >= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] <= combined_gwas_result$LD_end[i]
                        ) | (
                        gff[,4] <= combined_gwas_result$LD_start[i] &
                        gff[,4] <= combined_gwas_result$LD_end[i] &
                        gff[,5] >= combined_gwas_result$LD_start[i] &
                        gff[,5] >= combined_gwas_result$LD_end[i]
                        ),]

                    }
                }



                # If the results after matching have at least 1 row, write all the results to combined_gwas_result
                if(nrow(temp_gff) > 0 & as.integer(temp_gff[1,1]) == as.integer(combined_gwas_result[i,2])){
                    for(m in 1:nrow(temp_gff)){

                        if(m == 1){
                            combined_gwas_result$Gene_name[i] <- temp_gff[m,9]
                            combined_gwas_result$Gene_start[i] <- temp_gff[m,4]
                            combined_gwas_result$Gene_stop[i] <- temp_gff[m,5]
                            combined_gwas_result$Gene_description[i] <- temp_gff[m,10]
                        } else if(m > 1){
                            combined_gwas_result <- DataCombine::InsertRow(combined_gwas_result, NewRow = combined_gwas_result[i,], RowNum = i+1)
                            i <- i + 1
                            combined_gwas_result$Gene_name[i] <- temp_gff[m,9]
                            combined_gwas_result$Gene_start[i] <- temp_gff[m,4]
                            combined_gwas_result$Gene_stop[i] <- temp_gff[m,5]
                            combined_gwas_result$Gene_description[i] <- temp_gff[m,10]
                        }

                        # Remove any row that contains all NA
                        combined_gwas_result <- combined_gwas_result[rowSums(is.na(combined_gwas_result)) != ncol(combined_gwas_result),]
                    }
                }
            }

            i <- i + 1
        }
    }

    #######################################################################
    ## Save all GWAS Results
    #######################################################################

    # Order the table base on position number which is on column 3
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,3]), ]
    # Order the table base on chromosome number which is on column 2 
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,2]), ]
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$Trait), ]
    
    gwas_result_filename <- paste("GAPIT.combined.GWAS.Results.csv", sep = "")
    utils::write.csv(combined_gwas_result, file.path(output_path, gwas_result_filename), row.names = FALSE)

    gwas_result_list <- split(x = combined_gwas_result, f = combined_gwas_result$Trait)
    for (i in 1:length(gwas_result_list)) {
        gwas_result_filename <- paste("GAPIT.", names(gwas_result_list)[i], ".GWAS.Results.csv", sep = "")
        utils::write.csv(gwas_result_list[[i]], file.path(GAPIT_significant_save_path, gwas_result_filename), row.names = FALSE)
    }

    # output records with 5 lowest P.value SNPs
    combined_gwas_result_with_lowest_p_value = combined_gwas_result
    combined_gwas_result_with_lowest_p_value = combined_gwas_result_with_lowest_p_value[order(combined_gwas_result_with_lowest_p_value[,4]),]
    if(length(unique(combined_gwas_result_with_lowest_p_value[,4]))>5 & length(unique(combined_gwas_result_with_lowest_p_value[,1]))>5){
        combined_gwas_result_with_lowest_p_value = combined_gwas_result_with_lowest_p_value[combined_gwas_result_with_lowest_p_value[,1] %in% unique(combined_gwas_result_with_lowest_p_value[,1])[1:5],]
    }
    combined_gwas_result_with_lowest_p_value = combined_gwas_result_with_lowest_p_value[order(combined_gwas_result_with_lowest_p_value[,1]), ]
    gwas_result_filename <- paste("GAPIT.combined.GWAS.five.lowest.p_value.Results.csv", sep = "")
    utils::write.csv(combined_gwas_result_with_lowest_p_value, file.path(output_path, gwas_result_filename), row.names = FALSE)

    # output records with 5 most frequent SNPs
    combined_gwas_result_with_most_occur = combined_gwas_result
    temp_frequent_snps = combined_gwas_result %>%
        dplyr::select(1,2,3,4) %>%
        tidyr::drop_na() %>%
        dplyr::distinct() %>%
        dplyr::group_by_at(1) %>%
        dplyr::summarize(Count = dplyr::n()) %>%
        dplyr::arrange(dplyr::desc(Count)) %>%
        as.data.frame(stringsAsFactors = FALSE)
    combined_gwas_result_with_most_occur = combined_gwas_result_with_most_occur %>% dplyr::left_join(temp_frequent_snps, by = "SNP") %>% as.data.frame(stringsAsFactors = FALSE)
    if(nrow(temp_frequent_snps) > 5 & length(unique(combined_gwas_result_with_most_occur[,1]))>5){
        combined_gwas_result_with_most_occur = combined_gwas_result_with_most_occur[combined_gwas_result_with_most_occur[,1] %in% temp_frequent_snps[1:5,1],]
    }
    combined_gwas_result_with_most_occur = combined_gwas_result_with_most_occur[order(match(combined_gwas_result_with_most_occur[,1], temp_frequent_snps[,1])),]
    gwas_result_filename <- paste("GAPIT.combined.GWAS.five.most.occur.Results.csv", sep = "")
    utils::write.csv(combined_gwas_result_with_most_occur, file.path(output_path, gwas_result_filename), row.names = FALSE)

    if(all(c("LD_start", "LD_end") %in% colnames(combined_gwas_result))){
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$LD_end), ]
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$LD_start), ]
    }

    if(all(c("Haploblock_start", "Haploblock_stop") %in% colnames(combined_gwas_result))){
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$Haploblock_stop), ]
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$Haploblock_start), ]
    }

    # Order the table base on chromosome number which is on column 2
    combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,2]), ]

    #######################################################################
    ## Everything is done, return combined_gwas_result
    #######################################################################

    return(combined_gwas_result)

}