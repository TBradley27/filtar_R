#' @export
#' @importFrom magrittr %>%

get_canonical_AIRs = function(utr_lengths) {
	AIRs = readr::read_tsv(utr_lengths, col_names=TRUE)
	AIRs$start = 1
	AIRs$AIR = 100
	AIRs = AIRs[,c('tx_id','start','utr_length','AIR')]

	return(AIRs)
}
