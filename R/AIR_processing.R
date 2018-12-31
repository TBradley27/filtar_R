
#' Generate a file of affected isoform ratio (AIR) data for use in computing context++ scores when executing targetscan7 (10.7554/eLife.05005).
#' @param utr_lengths A two-column tab delimited file: tx_id, utr_length
#' @return AIRs A four-column tibble - tx_id, start, utr_length, AIR
#' @export
#' @importFrom magrittr %>%

get_canonical_AIRs = function(utr_lengths) {
	AIRs = readr::read_tsv(utr_lengths, col_names=TRUE)
	AIRs$start = 1
	AIRs$AIR = 100
	AIRs = AIRs[,c('tx_id','start','utr_length','AIR')]

	return(AIRs)
}
