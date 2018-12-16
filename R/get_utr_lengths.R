#' Returns the length of 3'UTRs from bed data

#' @param utr_bed_file A bed 6 file of 3'UTRs
#' @return A two-column dataframe of transcript IDs and corresponding 3'UTR lengths
#' @export

get_utr_lengths = function(utr_bed_file) {

        bed_records = readr::read_tsv(
                file = utr_bed_file,
                col_names = c('chrom','start','stop','strand','tx_id'),
                col_types = 'ciiic'
                )

        bed_records$exon_length = bed_records$stop - bed_records$start

        get_utr_lengths = function (id) {
                bed_subset = dplyr::filter(bed_records, tx_id == id)
                utr_length = sum(bed_subset$exon_length)

                x = list(tx_id=id, utr_length=utr_length)

                return (x)
        }

        tx_ids = bed_records$tx_id %>% unique

        utr_lengths = purrr::map(tx_ids, get_utr_lengths)
        utr_lengths = plyr::ldply(utr_lengths, data.frame) %>% tibble::as.tibble()

        return (utr_lengths)
}
