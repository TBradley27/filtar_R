#' Generates a miRNA family file in specific targetscan format
#' @param mirna_seed information - a table with miRNA and species identifiers as well as seed sequences
#' @param species three letter species code e.g. hsa,mmu etc - used as the reference species
#' @return A data frame in which miRNA famiilies have been enumerated, with seed sequences and NCBI species taxnonomic IDs. All families contain a seed sequence from the reference species
#' @export

get_mirna_family = function(mirna_seeds, species) {

        mirna_seeds = readr::read_tsv(mirna_seeds, col_names=c("identifier", "seq"))

        # group miRNA families together
        mirna_seeds = mirna_seeds[order(mirna_seeds$seq),]

        # return the three letter identifier for each species
        mirna_seeds$species = stringr::str_sub(mirna_seeds$identifier, 1,3)

        TaxID = list(hsa="9606", ptr="9598", ggo="9595",   pab="9601",   nle="9581",    rma="9544",   mfa="9541",
          pan="9557",   csa="60711",   cja="9483",   sbo="27679",   oga="30611",   tch="246437",   str="43179",   jja="51337",
          moc="79684",   cgr="10029",   mau="10036",   mmu="10090",   rno="10116",   hgl="10181",   cpo="10141",   cla="34839",
          ode="10160",   ocu="9986",   opr="9978",   ssc="9823",   vpa="30538",   cfe="419612",   ttr="9739",   oor="9733",
          pod="59538",   bta="9913",   oar="9940",   chi="9925",   eca="9796",   csi="9807",   fca="9685",   cfa="9615",
          mfu="9669",   ame="9646",   oro="9708",   lwe="9713",   pal="9402",   pva="132908",   mda="225400",   mlu="59463",
          efu="29078",   eeu="9365",   sar="42254",   ccr="143302",   laf="9785",   eed="28737",   tma="127582",   cas="185453",
          ete="9371",   oaf="1230840",   dno="9361",   mdo="13616",   sha="9305",   meu="9315",   oan="9258",   fch="345164",
          fpe="9854",   fal="59894",   zal="44394",   gfo="48883",   tgu="59729",   phu="181119",   mun="13146",   avi="241585",
          ama="176014",   cli="8932",  apl="8839",   gga="9031",  mga="9103",   ami="8496",   cmy="8469",  cpi="8478",
          psi="13735",   asp="55534",   aca="28377",   xtr="8364",   lch="7897")

	map_ids = function(string) {
  		tax_id = TaxID[[string]]

  		return (tax_id)
	}

        mirna_seeds = mirna_seeds[mirna_seeds$species %in% names(TaxID),]

        mirna_seeds$tax_id = purrr::map(mirna_seeds$species, map_ids)
        mirna_seeds$species = NULL

        mirna_seeds$tax_id = as.character(mirna_seeds$tax_id)

        # delete families which do not include a human orhtologue

        unique_seeds = mirna_seeds$seq %>% unique()

	delete_mirs_without_reference = function(string) {

	x = dplyr::filter(mirna_seeds, seq == string)

	if (TaxID[[species]] %in% x$tax_id) {
		return (x)
	}

	else {
		return ()
	}
	}

        mirna_seeds = purrr::map(unique_seeds, delete_mirs_without_reference)

        mirna_seeds = plyr::ldply(mirna_seeds, data.frame) %>% tibble::as.tibble()

#        print(mirna_seeds)

        mirna_seeds$seq = as.factor(mirna_seeds$seq)
        mirna_seeds$identifier = as.integer(mirna_seeds$seq)

        # remove duplicate lines

        mirna_seeds = unique( mirna_seeds[,] )

        return (mirna_seeds)
}

#' get miRNA file needed for computing contextpp scores
#' @param mirna_seed_file A file containing miRNA accessions with corresponding seed sequences
#' @param mature_mirnas_file A file containing miRNA accessions with corresponding mature miRNA sequences
#' @param species - three-letter species code
#' @return A file containing the name of the miRNA family, the NCBI species taxonomic ID, the name of the mature miRNA and also the mature miRNA sequence
#' @export

get_mirna_context = function (mirna_seed_file, mature_mirnas_file, species) {

	TaxID = list(hsa="9606", mmu="10090")

	specific_tax_id = TaxID[[species]]

	mirna_seeds = readr::read_tsv(mirna_seed_file, col_names=c("identifier", "seq", "tax_id"), col_types='dci')
#	print(mirna_seeds)

	mirnas = readr::read_tsv(mature_mirnas_file, col_names = c('miRNA_family','tax_id','miRNA','miRNA_sequence'), col_types='cicc')
#	print(mirnas)
	mirnas$seed = stringr::str_sub(mirnas$miRNA_sequence, 2, 8)

	test = merge(mirnas, mirna_seeds, by.x = 'seed', by.y = 'seq')
	test = subset(test, test$tax_id.y == specific_tax_id)

	test = test[!duplicated(test), ]
	test = test[,c('identifier','tax_id.x', 'miRNA', 'miRNA_sequence')]

return(test)
}

#' Get a full AIR file specific to a given coding transcriptome

#' @param APA_file - The output from APAtrap's APA prediction perl script
#' @param UTR_lengths_file - A dataframe containing transcript accessions with corresponding 3'UTR lengths
#' @export

get_AIR_file = function(APA_file,UTR_lengths_file) {
	APAtrap_output = read.table(APA_file, sep="\t",
				    header=TRUE)
	APAtrap_output = tibble::as.tibble(APAtrap_output)

	#APAtrap_output = APAtrap_output[1:1268,] # remove incomplete rows - I have no idea what this line is for

	reposition_last_APA_site = function (APA_df) {
	  
	 APA_df$strand = stringr::str_sub(APA_df$Gene, -1, -1) # retrieve the last character in the column which is the transcript strand
	 APA_df$last_APA = ''                         # create space for the last APA position
	 APA_df$start_pos = ''                        # create space for what I assume is the transcript start position
	  
	 # Account for strandedness
	 for (i in 1:dim(APA_df)[1]) {                # iterate over rows
	   if (APA_df$strand[i] == '+') {
	      
	     APA_df$last_APA[i] = sub('.*-', '', APA_df$Loci[i]) # remove everything before the hyphen
	     APA_df$start_pos[i] = sub('[0-9A-Z]*:', '', APA_df$Loci[i]) %>% 
		sub(pattern='-.*',replacement='') # capture everything after the colon and before the hyphen
	      
	    } else if (APA_df$strand[i] == '-')
	    {
	      
	      APA_df$last_APA[i] = sub('[0-9A-Z]*:', '', APA_df$Loci[i]) %>% 
		sub(pattern='-.*',replacement='')
	      APA_df$start_pos[i] = sub('.*-', '', APA_df$Loci[i])
	      
	    }
	  }
	  APA_df$Predicted_APA = stringr::str_glue_data(APA_df,"{Predicted_APA},{last_APA}") # collect all APA sites together
	  
	  # clean up
	  APA_df$last_APA = NULL
	  return (APA_df)
	}

	APAtrap_output = reposition_last_APA_site(APAtrap_output)

#	print(APAtrap_output)
	#convert_to_percentages = function (exp_string, total) {
	#  expression_values = str_split(exp_string, ',') %>% unlist %>% as.numeric
	  #percentage_values = ( expression_values / total ) * 100
	#  return (percentage_values)
	#}

	# I think this function is to remove APA sites with a quoted abundance of 0.00 which in effect would mean that they are redundant - probably initially identified because they have the correct motif

	remove_unused_APAsites = function (APA_df) {
	  converted_exp_column = purrr::map(APA_df$Group_1_1_Separate_Exp, function (x)
	    stringr::str_split(x, ',') %>% unlist %>% as.numeric) # convert to numerical vector
	  
	  converted_loci_column = purrr::map(APA_df$Predicted_APA, function (x)
	    stringr::str_split(x, ',') %>% unlist) # convert to numerical vector
	  
	  for (i in 1:length(converted_exp_column)) { # loops over each row of the data frame
	    index = converted_exp_column[[i]] != 0 # A Boolean vector
	    converted_exp_column[[i]] = converted_exp_column[[i]][index]
	    converted_loci_column[[i]] = converted_loci_column[[i]][index]
	  }
	  
	  APA_df$a = converted_exp_column
	  APA_df$b = converted_loci_column
	  
	  return (APA_df)
	  
	}

	APAtrap_output = remove_unused_APAsites(APAtrap_output)
#	print(APAtrap_output)

	get_relative_abundances = function (exp_string, total, tx_id) {
	  
	  # function applied recursively
	  get_cumulative_depreciation = function (percent_vec) {
	    cumulative_depreciation = c() #initialise vector 
	    cumulative_depreciation = c(cumulative_depreciation, 100) # initialise vector 
	    for (i in 1:length(percent_vec) - 1) { 
	      percent_fall = percent_vec[i]
	      cumulative_depreciation = c(
		cumulative_depreciation, 
		tail(cumulative_depreciation, n=1) - percent_fall
	      )
	    }
	    return (cumulative_depreciation)
	  }
	    
	  #expression_values = str_split(exp_string, ',') %>% unlist %>% as.numeric
	  percentage_values = ( exp_string / total ) * 100

	  return (get_cumulative_depreciation(percentage_values) )
	  
	}

	# get relative abundances

	y = purrr::map2(APAtrap_output$a,
		    APAtrap_output$Group_1_1_Total_Exp,
		    get_relative_abundances)

	names(y) = APAtrap_output$Gene
	y = unlist(y) %>% as.data.frame
	y$id = rownames(y)
	rownames(y) = NULL
	y$id = gsub('\\|.*','',y$id)
	y = y[,c(2,1)] # reorder columns

	new = APAtrap_output
	new$Gene = gsub('\\|.*','',APAtrap_output$Gene) 
	new = new[,c('Gene','b','Loci','strand','start_pos')]

	y = merge(y, new, by.x='id', by.y='Gene')

	y$Loci = NULL
	#y$strand = NULL

	#print(y) # end position missing

	tx_ids = y$id %>% as.factor %>% levels

	get_rel_APA_position = function(tx_id) {
	  records = subset(y, y$id == tx_id)
	  
	  absolute_start = as.numeric(records$start_pos[1])
	  
	  for (i in 1:dim(records)[1] ) {
	    
	    if (i == 1) {
		if (records$strand[i] == '+') {
	      		records$rel_start_pos[i] = 1 
	      		records$rel_end_pos[i] = as.numeric(records$b[[1]][i]) - absolute_start + 1 # initial relative end equals first APA site - absolute start
		} else if (records$strand[i] == '-') {
			records$rel_start_pos[i] = 1
			records$rel_end_pos[i] = absolute_start - as.numeric(records$b[[1]][i]) + 1
		}
	    }   else {
			if (records$strand[i] == '+') {
	      			records$rel_start_pos[i] = records$rel_end_pos[i-1] + 1
	      			records$rel_end_pos[i] = as.numeric(records$b[[1]][i]) - absolute_start + 1
		} else if (records$strand[i] == '-') {
				records$rel_start_pos[i] = records$rel_end_pos[i-1] + 1
                                records$rel_end_pos[i] = absolute_start - as.numeric(records$b[[1]][i]) + 1
			}
	    }
	  }
	  #print(absolute_start) 
	  #print(records) 
	  records$b = NULL
	  records$start = NULL
	  records = records[,c(1,5,6,2)]
	  colnames(records) = c('id','rel_start_pos','rel_end_pos','AIR')
	  
	  return (records)
	}

	APA_records = purrr::map(tx_ids, get_rel_APA_position) %>% plyr::ldply(data.frame)
	APA_records = tibble::as.tibble(APA_records)

	all_transcripts = readr::read_tsv(UTR_lengths_file, col_names=TRUE)

	# add version numbers to all APA records

	all_transcripts_sep = tidyr::separate(all_transcripts, tx_id, into=c('id','version'))
#	print(all_transcripts)
	#APA_records = merge(APA_records,all_transcripts_sep,by='id')
        APA_records = plyr::join(APA_records, all_transcripts_sep, by=c('id'))
	APA_records$utr_length = NULL
#	print(APA_records)
	APA_records = tidyr::unite(APA_records, col='id', c('id','version'), sep=".")	
#	print(APA_records)
	non_updated_tx = all_transcripts[!all_transcripts$tx_id %in% APA_records$id,]
#	print(non_updated_tx)
	non_updated_tx$start = 1
	non_updated_tx$AIR = 100
	non_updated_tx = non_updated_tx[,c('tx_id','start','utr_length','AIR')]
	colnames(non_updated_tx) = c('id','rel_start_pos','rel_end_pos','AIR')
#	print(non_updated_tx)

	APA_records = rbind(APA_records,non_updated_tx) %>% tibble::as.tibble()
	
	return(APA_records)
}

#' fix output from a modified version of the targetscan script which finds seed matches

#' @param ts_sites_output - Output from the modified version of the first targetscan script
#' @return A corrected targetscan seed match output file
#' @export

fix_ts_output = function (ts_sites_output) {

        ts_sites = readr::read_tsv(ts_sites_output)

        sort_species = function (string) {
          if (grepl(' ', string)) {
                 sorted_string = string %>%
                 strsplit(split=' ', fixed=TRUE) %>%
                 unlist %>%
                 as.numeric %>%
                 unique %>%
                 sort(method='quick') %>%
                 paste(collapse=" ")
                 return(sorted_string)
          } else {
                return (string)
         }
        }

        # reference in order to sort site-types into correct ordering
        y = c('7mer-1a','7mer-m8','8mer-1a','6mer')

        remove_redundant_site_types = function (record) {
        #  if (grepl('+',record)) {
                new_record = record %>% strsplit('+', fixed=TRUE) %>% unlist %>% unique
                new_record =  new_record[order(match(new_record,y))] %>%
                paste(collapse="+")
                return (new_record)
        #  }
        #  else {
        #       return (record)
        #  }
        }

        dummy = vector(length=dim(ts_sites)[1])
        dummy2 = vector(length=dim(ts_sites)[1])
        dummy3 = vector(length=dim(ts_sites)[1])

        for (i in 1:dim(ts_sites)[1] ) {
 #               print(  (i/dim(ts_sites)[1]) * 100 ) # progress bar
                dummy[i] = sort_species(ts_sites$Species_in_this_group[i])
                dummy2[i] = remove_redundant_site_types(ts_sites$Group_type[i])

                if ( grepl('\\+', dummy2[i]) ) {
                        dummy3[i] =
                        ts_sites$Species_in_this_group_with_this_site_type[i] %>%
                        strsplit(' ', fixed=TRUE) %>%
                        unlist %>%
                        as.numeric %>%
                        unique %>%
                        sort() %>%
                        paste(collapse=" ")
                } else {
                        dummy3[i]=''
                }

        }

        ts_sites$Species_in_this_group = dummy
        ts_sites$Group_type = dummy2
        ts_sites$Species_in_this_group_with_this_site_type = dummy3

        return(ts_sites)
}

#' @export

filter_contextpp_scores = function (contextpp_scores_filename, expression_values_filename) {

	contextpp_scores = readr::read_tsv(contextpp_scores_filename, col_names=TRUE)
	expression_values = readr::read_tsv(expression_values_filename, col_names=TRUE)

	merged_dataset = merge(contextpp_scores, expression_values, by.x='Gene ID', by.y='Name')
	filtered_merged_dataset = dplyr::filter(merged_dataset, TPM >= snakemake@params['tpm_expression_threshold'])

	return(filtered_merged_dataset)
}








