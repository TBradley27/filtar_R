#' Filter full set of mature miRNAs by a list of user selected miRNAs
#' @param input_file A tab delimited file of mature miRNA sequences and associated identifiers
#' @param list_of_mirs A vector of miRNA names
#' @param output_file Name of file to write to
#' @return NULL
#' @export

filter_mature_mirs = function (input_file, list_of_mirs, output_file) {
	list_of_mirs_suffix = substring(list_of_mirs,5)
	if (length(list_of_mirs) == 0 | list_of_mirs_suffix == 'all') { # if the user does not specify a miRNA - then process all species miRNA as default
		file.copy(from=input_file, to=output_file)
		return(NULL)
	}
	else {		
		df = readr::read_tsv(input_file, col_names=c('miR_family_id','species_id','mature_mirna_name','sequence'), col_types='cccc')

		for (mir in list_of_mirs) {
			if (!mir %in% df$mature_mirna_name) {
				write(stringr::str_interp("miRNA identifier '${mir}' not found in data/mature.fa"), stderr())
			}
		}

		df = df[df$mature_mirna_name %in% list_of_mirs,]
		
		if (dim(df)[1] == 0) {
			write("No valid miRNA identifers in input. Halting execution.", stderr())
			quit(save="no", status=1)			
		}

		utils::write.table(df, output_file, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
		return(NULL)
	}
}

#' Filter miRNA family data according to a vector of miRNA name
#' @param mature_mirna_file file name of mature miRNA data
#' @param mirna_families_file file name of miRNA family data
#' @param list_of_mirs A vector of miRNA names
#' @param output_file A filtered set of miRNA family data
#' @return NULL
#' @export

filter_mir_families = function (mature_mirna_file, mirna_families_file, list_of_mirs, output_file) {
        if (length(list_of_mirs) == 0 | list_of_mirs == 'all') { # if the user does not specify a miRNA - then process all species miRNA as default
                file.copy(from=mirna_families_file, to=output_file)
                return(NULL)
        }
        else {
                mature_mirnas = readr::read_tsv(mature_mirna_file, col_names=c('miR_family_id','species_id','mature_mirna_name','sequence'), col_types='cccc')
		mirna_families = readr::read_tsv(mirna_families_file, col_names=c('miR_family_id','seed_sequence','species_id'), col_types='ccc')

		mature_mirnas$seed_sequence = stringr::str_sub(mature_mirnas$sequence, 2, 8)
		mirna_families_wanted = mirna_families[mirna_families$seed_sequence %in% mature_mirnas$seed_sequence,]
                
		utils::write.table(mirna_families_wanted, output_file, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                return(NULL)
        }
}
