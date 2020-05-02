#' Extract haplotype using Haploview
#'
#' @description The goal of search_genes is to extract haplotype blocks and append the results to combined_gwas_result.
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @param output_path An output path.
#' @return combined_gwas_result or NULL if something missing.
#' @keywords Haplotype_blocks
#' @export
#'
extract_haplotype <- function(combined_gwas_result = NULL,
                                output_path = file.path("~"),
                                Haploview_file_path = NULL, 
                                Haploview_file_name = NULL, 
                                Haploview_file_extension = NULL, 
                                Haploview_file_named_sequentially_from = NULL, 
                                Haploview_file_named_sequentially_to = NULL) {

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
    ## Create Haploview reference file path and put everything in a table
    #######################################################################

    Haploview_file_table <- NULL
    Haploview_file_table <- tryCatch({
        data.frame(
            "Chromosome" = Haploview_file_named_sequentially_from:Haploview_file_named_sequentially_to, 
            "File_path" = file.path(Haploview_file_path, 
                                    paste0(Haploview_file_name, 
                                            Haploview_file_named_sequentially_from:Haploview_file_named_sequentially_to, 
                                            ".", 
                                            Haploview_file_extension
                                            )
                                    )
        )
    }, error = function(e) {
        print("Haploview_file_table cannot be created!!!")
        return(NULL)
    })

    if(!is.null(Haploview_file_table)){
        if(nrow(Haploview_file_table) > 0){
            for(i in 1:nrow(Haploview_file_table)){
                if(!file.exists(file.path(Haploview_file_table[i, 2]))){
                    print("File path for Haploview does not exists.")
                    return(NULL)
                }
            }
        } else{
            print("Haploview_file_table has zero row.")
            return(NULL)
        }
    }

    #######################################################################
    ## Run Haploview here
    #######################################################################

    if(!is.null(combined_gwas_result)){
        combined_gwas_result <- as.data.frame(combined_gwas_result)
        combined_gwas_result$Haploblock_number <- NA
        combined_gwas_result$Haploblock_start <- NA
        combined_gwas_result$Haploblock_stop <- NA

        index <- match(as.integer(combined_gwas_result[1,2]), Haploview_file_table[,1])
        index <- ifelse(is.na(index), 1, index)

        haploview_ref <- read_file(file_path = file.path(Haploview_file_table[index, 2]))
        if(is.null(haploview_ref)){
            print("Unable to read haploview reference file or combined GWAS results are all NA.")

        } else{
            # Change colnames of haploview_ref
            colnames(haploview_ref)[1] <- "Chromosome"
            colnames(haploview_ref)[2] <- "Positions"

            # For each row in each table
            i <- 1
            while(i <= nrow(combined_gwas_result)){

                if(nrow(haploview_ref) > 0 & !is.na(combined_gwas_result[i, 2]) & !is.na(combined_gwas_result[i, 3]) & !is.na(combined_gwas_result$LD_start[i]) & !is.na(combined_gwas_result$LD_end[i])){

                    # If chromosome number is different, read new haploview_ref that matches the chromosome number
                    if(as.numeric(haploview_ref$Chromosome[1]) != as.numeric(combined_gwas_result[i, 2])){

                        index <- match(as.integer(combined_gwas_result[i, 2]), Haploview_file_table[,1])
                        index <- ifelse(is.na(index), 1, index)

                        haploview_ref <- read_file(file_path = file.path(Haploview_file_table[index, 2]))
                        if(is.null(haploview_ref)){
                            print("Unable to read haploview reference file.")
                            return(NULL)
                        }
                        # Change colnames of haploview_ref
                        colnames(haploview_ref)[1] <- "Chromosome"
                        colnames(haploview_ref)[2] <- "Positions"
                    }

                    if(nrow(haploview_ref) > 0 & haploview_ref$Chromosome[1] == as.numeric(combined_gwas_result[i, 2])){
                        # Create a temporary haploview_ref from chromosome, LD start, and LD stop
                        temp_haploview_ref <- haploview_ref[(haploview_ref$Positions >= combined_gwas_result$LD_start[i] &
                                                                haploview_ref$Positions <= combined_gwas_result$LD_end[i] &
                                                                haploview_ref$Chromosome == combined_gwas_result[i, 2]), ]

                        if(nrow(temp_haploview_ref) > 0 & temp_haploview_ref$Chromosome[1] == as.numeric(combined_gwas_result[i, 2])){
                            # Put positions to info file
                            info <- data.frame(temp_haploview_ref$Positions, temp_haploview_ref$Positions)

                            # Remove Chromosome and Positions columns
                            temp_haploview_ref <- temp_haploview_ref[, c(-1, -2)]

                            # Prepare ped file
                            ibd <- paste("IBD",1:ncol(temp_haploview_ref), sep = "")
                            ped <- rbind(temp_haploview_ref[0,], ibd, 0, 0, 7, 1, temp_haploview_ref[c(1:nrow(temp_haploview_ref)),])
                            ped <- t(ped)

                            filename <- paste0(combined_gwas_result$Trait[i], "_", combined_gwas_result[i,2], "_",
                                                combined_gwas_result[i,3], "_", combined_gwas_result$LD_start[i], "-",
                                                combined_gwas_result$LD_end[i])

                            # Name ped and info file
                            ped_file_name <- paste0(filename, ".ped")
                            info_file_name <- paste0(filename, ".info")

                            # Save info and ped file
                            write.table(info, file.path(ped_and_info_save_path, info_file_name), sep = '\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
                            write.table(ped, file.path(ped_and_info_save_path, ped_file_name), sep = '\t', row.names = TRUE, col.names = FALSE, quote=FALSE)

                            # Prepare the Haploview command
                            ld_data_command <- paste("java -jar ", file.path(getwd(), "Haploview.jar"),
                                                        " -n -out ", file.path(ld_data_save_path, filename),
                                                        " -pedfile ", file.path(ped_and_info_save_path, ped_file_name),
                                                        " -info ", file.path(ped_and_info_save_path, info_file_name),
                                                        " -skipcheck -dprime -png -ldcolorscheme DEFAULT -ldvalues DPRIME -blockoutput GAB -minMAF 0.05",
                                                        sep = "")

                            # Run Haploview
                            system(ld_data_command)

                            # Move all the plots into their specific folders
                            if(file.exists(file.path(ld_data_save_path, paste(filename, ".LD.PNG", sep = "")))){
                                system(paste("mv", file.path(ld_data_save_path, paste(filename, ".LD.PNG", sep = "")), ld_plot_save_path, sep = " "))
                            }
                            if(file.exists(file.path(ld_data_save_path, paste(filename, ".GABRIELblocks", sep = "")))){
                                system(paste("mv", file.path(ld_data_save_path, paste(filename, ".GABRIELblocks", sep = "")), haplotypes_gabriel_blocks_save_path, sep = " "))
                            }

                            # Read in LD_data and gabriel_block_string
                            if(file.exists(file.path(ld_data_save_path, paste(filename, ".LD", sep = "")))){
                                LD_data <- try(read.table(file.path(ld_data_save_path, paste(filename, ".LD", sep = "")), check.names = FALSE, header = TRUE))
                            }
                            if(file.exists(file.path(haplotypes_gabriel_blocks_save_path, paste(filename, ".GABRIELblocks", sep = "")))){
                                gabriel_block_string <- try(readLines(file.path(haplotypes_gabriel_blocks_save_path, paste(filename, ".GABRIELblocks", sep = ""))))
                            }

                            # Make sure all the Haploview output exists in the directory and exists in the workspace as well
                            if(file.exists(file.path(ld_data_save_path, paste(filename, ".LD", sep = ""))) &
                                file.exists(file.path(haplotypes_gabriel_blocks_save_path, paste(filename, ".GABRIELblocks", sep = ""))) &
                                exists("LD_data") & exists("gabriel_block_string")){

                                # Make sure the LD_data and gabriel_block_string are not totally empty
                                if(!identical(LD_data, character(0)) & !identical(gabriel_block_string, character(0))){
                                    # Get Haploblock start and stop
                                    ld <- sort(as.double(unique(c(LD_data[,1], LD_data[,2]))))

                                    # Parse gabriel blocks data
                                    gbb <- list()
                                    for (k in 1:length(gabriel_block_string)) {
                                        if(grepl("MARKERS: ", gabriel_block_string[k], ignore.case = TRUE)){
                                            gbb <- append(gbb, strsplit(x = gsub(".*MARKERS: ", "" , gabriel_block_string[k]), split = " "))
                                        }
                                    }

                                    # Write the number of haploblocks to the corresponding row of Haploblock_number column
                                    combined_gwas_result$Haploblock_number[i] <- as.numeric(length(gbb))

                                    m <- length(gbb)
                                    n <- length(gbb[[m]])

                                    if(!is.na(as.integer(gbb[[m]][n])) & length(ld) >= as.integer(gbb[[m]][n])){

                                        # Get all the start and stop of markers from ld and gbb
                                        for(m in 1:length(gbb)){
                                            for (n in 1:length(gbb[[m]])) {
                                                gbb[[m]][n] <- as.double(ld[as.integer(gbb[[m]][n])])
                                            }
                                        }

                                        # Write the Haploblock_start and Haploblock_stop that enclose the position
                                        for(m in 1:length(gbb)){
                                            if(as.double(combined_gwas_result[i,3]) >= as.double(gbb[[m]][1]) &
                                                as.double(combined_gwas_result[i,3]) <= as.double(gbb[[m]][length(gbb[[m]])]) ){
                                                    combined_gwas_result$Haploblock_start[i] <- as.double(gbb[[m]][1])
                                                    combined_gwas_result$Haploblock_stop[i] <- as.double(gbb[[m]][length(gbb[[m]])])
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                i <- i + 1
            }
        }

        # Remove any row that contains all NA
        combined_gwas_result <- combined_gwas_result[rowSums(is.na(combined_gwas_result)) != ncol(combined_gwas_result),]

        #######################################################################
        ## Save all GWAS Results
        #######################################################################

        # Order the table base on position number which is on column 3
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,3]), ]
        # Order the table base on chromosome number which is on column 2
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,2]), ]
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$Trait), ]

        gwas_result_filename <- paste("GAPIT.combined.GWAS.Results.csv", sep = "")
        write.csv(combined_gwas_result, file.path(output_path, gwas_result_filename), row.names = FALSE)

        gwas_result_list <- split(x = combined_gwas_result, f = combined_gwas_result$Trait)
        for (i in 1:length(gwas_result_list)) {
            gwas_result_filename <- paste("GAPIT.", names(gwas_result_list)[i], ".GWAS.Results.csv", sep = "")
            write.csv(gwas_result_list[[i]], file.path(GAPIT_significant_save_path, gwas_result_filename), row.names = FALSE)
        }

        # output records with 5 lowest P.value SNPs
        combined_gwas_result_with_lowest_p_value = combined_gwas_result
        combined_gwas_result_with_lowest_p_value = combined_gwas_result_with_lowest_p_value[order(combined_gwas_result_with_lowest_p_value[,4]),]
        if(length(unique(combined_gwas_result_with_lowest_p_value[,4]))>5 & length(unique(combined_gwas_result_with_lowest_p_value[,1]))>5){
            combined_gwas_result_with_lowest_p_value = combined_gwas_result_with_lowest_p_value[combined_gwas_result_with_lowest_p_value[,1] %in% unique(combined_gwas_result_with_lowest_p_value[,1])[1:5],]
        }
        combined_gwas_result_with_lowest_p_value = combined_gwas_result_with_lowest_p_value[order(combined_gwas_result_with_lowest_p_value[,1]), ]
        gwas_result_filename <- paste("GAPIT.combined.GWAS.five.lowest.p_value.Results.csv", sep = "")
        write.csv(combined_gwas_result_with_lowest_p_value, file.path(output_path, gwas_result_filename), row.names = FALSE)

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
        write.csv(combined_gwas_result_with_most_occur, file.path(output_path, gwas_result_filename), row.names = FALSE)

        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$LD_end), ]
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$LD_start), ]
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$Haploblock_stop), ]
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result$Haploblock_start), ]
        # Order the table base on chromosome number which is on column 2
        combined_gwas_result <- combined_gwas_result[order(combined_gwas_result[,2]), ]

        #######################################################################
        ## Everything is done, return combined_gwas_result
        #######################################################################

        return(combined_gwas_result)

    } else{
        return(NULL)
    }

    return(NULL)
}