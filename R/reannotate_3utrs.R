# reannotates transcript records which are found in both a normal bed file and an extended bed file in order to match the latter.
# multi-exon 3'UTRs are not reannotated

#  tx_id A single transcript accession
#  An update or unchanged single transcript bed record(s)
#




# reorder bed records within a given transcript ID to ensure correct standard ordering
#  tx_ID A single transcript accession
#  All bed records correspodning to tx_ID with correct ordering 
# 


#' Generates a full comprehensive 3'UTR bed 6 file by integrating information from a canonical bed 6 file and APAtrap reannotations

#' @param normal_bed_file A 3'UTR bed 6 file in standard format
#' @param extended_bed_file APAtrap output. Some records will correspond to accessions in a normal_bed_file, whilst some records will reference accessions not found in a normal_bed_file
#' @param all_transcripts_file A single column list of all transcript files used to re-add transcript version information to all recrods
#' @return A complete bed6 file in standard format in which information from all parameters have been integrated, and in which transcript accession verison information is included
#' @export

get_full_bed = function (normal_bed_file, extended_bed_file, all_transcripts_file) {

	# read in the data
	normal_bed = readr::read_tsv(normal_bed_file, col_names=c('chromosome','start','stop','strand','id'), col_types=list('c','i','i','c','c'))
	extended_utrs = readr::read_tsv(extended_bed_file, col_names=c('chromosome','start','stop','id','dummy','strand'))

	print('bar')

	# tidy the data
	normal_bed = tidyr::separate(normal_bed, id, into=c('id','version'))
	print(extended_utrs)
	extended_utrs = tidyr::separate(extended_utrs, id, into=c('id','dummy2','chrom_dup','strand_dup'))
	print('bar')
	normal_bed$version = NULL
	print('bar')
	extended_utrs = extended_utrs[,c('chromosome','start','stop','strand','id')]
	extended_utrs$chromosome = stringr::str_replace_all(extended_utrs$chromosome, 'chr','')
	normal_bed$chromosome = stringr::str_replace_all(normal_bed$chromosome, 'chr','')

	print('bar')

	##

	print(normal_bed)
	print(extended_utrs)

	tx_ids = normal_bed[normal_bed$id %in% extended_utrs$id,]
	tx_ids = tx_ids$id %>% unique()

	reannotate = function(tx_id) {

		canonical_tx_record = dplyr::filter(normal_bed, id==tx_id)
		apatrap_tx_record = dplyr::filter(extended_utrs, id==tx_id)

		if (dim(apatrap_tx_record)[1] == 0) { # no record
			return (canonical_tx_record)
		} else {

			if (dim(canonical_tx_record)[1] == 1) { # single exon

				new_tx_record = canonical_tx_record
				new_tx_record$start = apatrap_tx_record$start
				new_tx_record$stop = apatrap_tx_record$stop

				return(new_tx_record)

			} else {   # I think it is best to do this as APAtrap does not seem to acknowledge multi-exon 3UTRs 

				return(canonical_tx_record)

		#	if (apatrap_tx_record$strand[1] == '+') {

		#		new_tx_record = canonical_tx_record
		#        	new_tx_record$start[1] = apatrap_tx_record$start
		#        	new_tx_record$stop[dim(canonical_tx_record)[1]] = apatrap_tx_record$stop

		#		return (new_tx_record)

	#}	#	else if (apatrap_tx_record$strand[1] == '-') {

		#		new_tx_record = canonical_tx_record
		#        	new_tx_record$stop[1] = apatrap_tx_record$stop[1]
		#        	new_tx_record$start[dim(canonical_tx_record)[1]] = apatrap_tx_record$start

		#		return (new_tx_record)

			}
		}
	}

	old_records_changed = purrr::map(tx_ids, reannotate)

	old_records_changed = plyr::ldply(old_records_changed, data.frame) %>% as.tibble()

	print('old records changed')
	old_records_changed %>% select(id) %>% table() %>% print()

	new_records = extended_utrs[!extended_utrs$id %in% normal_bed$id,]
	#new_records$chromosome = paste('chr',new_records$chromosome,sep='')
	new_records$strand = gsub('\\+','1',new_records$strand)
	new_records$strand = gsub('-','-1',new_records$strand)

	print('new records')
	print(new_records)

	#new_records %>% select(id) %>% table() %>% print()

	old_records_unchanged = normal_bed[!normal_bed$id %in% extended_utrs$id,]

	print('old records unchanged')
	print(old_records_unchanged)

	full_set = rbind(old_records_changed, new_records, old_records_unchanged)

	# reformat bed file

	full_set = full_set[order(full_set$id, decreasing=FALSE),]
	#full_set = full_set[order(full_set$start, decreasing=FALSE),]

	#full_set = full_set[order(full_set$chromosome, decreasing=FALSE),] # ordering of chromosomes within the BED file
	full_set$chromosome = as.character(full_set$chromosome)

	# add version number back in

	all_transcripts = readr::read_tsv(file=all_transcripts_file, col_names=c('id'))
	all_transcripts = tidyr::separate(all_transcripts, id, into=c('id','version'))
	all_transcripts = all_transcripts %>% distinct # remove duplicate entries

	print(all_transcripts)

	x = merge(full_set, all_transcripts, by.x='id', by.y='id')
	x = tidyr::unite(x, id, c('id','version'), sep = ".", remove = TRUE)

	x = x[,c('chromosome','start','stop','strand','id')]
	print(x)

	full_set = x

	full_set = full_set[order(full_set$start, decreasing=FALSE),]

	print('just after tx start position ordering')

	tx_IDs = full_set$id %>% unique()

	reorder_bed_files = function (tx_ID) { # ordering of exons within a transcript
	transcript_specific_bed_records = subset(full_set, full_set$id == tx_ID)
	#print(transcript_specific_bed_records)
	if (transcript_specific_bed_records$strand[1] == '1') {
	   txs = transcript_specific_bed_records[order(transcript_specific_bed_records$start, decreasing=FALSE),]
	}
	else if (transcript_specific_bed_records$strand[1] == '-1' ) {
	   txs = transcript_specific_bed_records[order(transcript_specific_bed_records$start, decreasing=TRUE),]
	}
	#print(transcript_specific_bed_records)
	return (txs)
	}

	full_set_sorted = purrr::map(tx_IDs, reorder_bed_files)
	print(full_set_sorted[1:30])
	full_set_sorted = plyr::ldply(full_set_sorted, data.frame) %>% as.tibble()
	print(full_set_sorted[1:200,], n=Inf)

	return(full_set_sorted)
}
