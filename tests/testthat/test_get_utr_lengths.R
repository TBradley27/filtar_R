context('test get_utr_lengths')

output = get_utr_lengths('ensembl_example2.bed')
output$tx_id = as.character(output$tx_id)

expected_output = readr::read_tsv('expected_utr_lengths.tsv', col_names=c('tx_id','utr_length'), col_types='ci')

test_that('actual output and expected output is equal',{
	expect_equal(output,expected_output)
})
