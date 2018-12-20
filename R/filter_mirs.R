#' @export

filter_mature_mirs = function (input_file, list_of_mirs, output_file) {
	if (length(list_of_mirs) == 0) { # if the user does not specify a miRNA - then process all species miRNA as default
		file.copy(from=input_file, to=output_file)
		return(NULL)
	}
	else {
		df = readr::read_tsv(input_file, col_names=c('miR_family_id','species_id','mature_mirna_name','sequence'), col_types='cccc')
		df = df[df$mature_mirna_name %in% list_of_mirs,]
		write.table(df, output_file, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
		return(NULL)
	}
}

#' @export

filter_mir_families = function (mature_mirna_file, mirna_families_file, list_of_mirs, output_file) {
        if (length(list_of_mirs) == 0) { # if the user does not specify a miRNA - then process all species miRNA as default
                file.copy(from=input_file, to=output_file)
                return(NULL)
        }
        else {
                mature_mirnas = readr::read_tsv(mature_mirna_file, col_names=c('miR_family_id','species_id','mature_mirna_name','sequence'), col_types='cccc')
		mirna_families = readr::read_tsv(mirna_families_file, col_names=c('miR_family_id','seed_sequence','species_id'), col_types='ccc')

		mature_mirnas$seed_sequence = stringr::str_sub(mature_mirnas$sequence, 2, 8)
		mirna_families_wanted = mirna_families[mirna_families$seed_sequence %in% mature_mirnas$seed_sequence,]
                
		write.table(mirna_families_wanted, output_file, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                return(NULL)
        }
}
