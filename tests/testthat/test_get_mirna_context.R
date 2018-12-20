library(plyr)
library(tidyverse)
library(filtar)

context('test get_mirna_context')

run_tests = function(mirna_family_filename,mature_mirna_filename,species) {
	output = get_mirna_context(mirna_family_filename,mature_mirna_filename,species)

	test_that ('The identifier column only contains 1-4 repetitions of 0-9',{
		expect_match( output$identifier %>% as.character, '^[0-9]{1,4}$' )
	})

	test_that('The first column is monotonic', {
		expect_true(all(output$identifier == cummax(output$identifier)))
	})

	if (species == 'hsa') {

		test_that('The tax id column is correct', {
			expect_match(as.character(output$tax_id.x), '9606')
		})

		test_that('Seed sequences do not contain invalid characters', {
			expect_match( output$miRNA, "^hsa-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
		})

		input_seeds = read_tsv(mirna_family_filename, col_names=c('identifier','seq','tax_id')) %>% filter(tax_id=='9606')


	} else if (species == 'mmu') {

		test_that('The tax id column is correct', {
			expect_match(as.character(output$tax_id.x), '10090')
		})

		test_that('Seed sequences do not contain invalid characters', {
			expect_match( output$miRNA, "^mmu-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
		})

		input_seeds = read_tsv(mirna_family_filename, col_names=c('identifier','seq','tax_id')) %>% filter(tax_id=='10090')

	}

	test_that('mature seqeunce has the right length and characters', {
		expect_match( output$miRNA_sequence, '^[ACTGU]{18,24}$' )
	})

	test_that('A single seed can lead to multiple mature sequences', {
		expect_true(output$identifier %>% table %>% max > 1)
	})

	test_that('every reference seed sequence has a corresponding mature sequence', {
		expect_match(input_seeds$seq %in% substring(output$miRNA_sequence, 2, 8) %>% as.character, as.character(TRUE))
	})

	if (species == 'hsa') {

		test_that('identifiers map to the correct sequence name', {
			expect_equal(as.character(output$identifier[output$miRNA == 'hsa-miR-1a-3p']),'1')
		})

		test_that('identifiers map to the correct sequence', {
			expect_equal(as.character(output$identifier[output$miRNA_sequence == 'AAGCGCGCCCTCTCGCGAGA']),'1')
		})

		test_that('sequence maps to the correct sequence name', {
			expect_equal(as.character(output$miRNA[output$miRNA_sequence == 'AAGCGCGCCCTCTCGCGAGA']),'hsa-miR-1a-3p')
		})

		test_that('output does not contain non-reference families',{
			expect_match('CGAAAAG' %in% output$seq %>% as.character, 'FALSE')
		})

		expected_output = tibble(identifier=c(1,1,2,3),tax_id.x=c(as.integer(9606),as.integer(9606),as.integer(9606),as.integer(9606)),miRNA=c('hsa-miR-1a-3p','hsa-miR-1b-3p','hsa-miR-2-5p',
					'hsa-miR-3-5p'),
					miRNA_sequence=c('AAGCGCGCCCTCTCGCGAGA','AAGCGCGCCCTCTCGCGTGA','CCCTCTCGAAAAAAAGCGCG','CGCGCGAGAAAAAAAGCGCG'))
		test_that('expected output matches actual output', {
			expect_equal(expected_output, output)
		})


	} else if (species == 'mmu' ) {

		test_that('identifiers map to the correct sequence name', {
			expect_equal(as.character(output$identifier[output$miRNA == 'mmu-miR-4-5p']),'1')
		})

		test_that('identifiers map to the correct sequence', {
			expect_equal(as.character(output$identifier[output$miRNA_sequence == 'CAAAAAAAGCGCGCGCGCGT']),'1')
		})

		test_that('sequence maps to the correct sequence name', {
			expect_equal(as.character(output$miRNA[output$miRNA_sequence == 'CAAAAAAAGCGCGCGCGCGT']),'mmu-miR-4-5p')
		})

		test_that('output does not contain non-reference families',{
			expect_match('GATTGCG' %in% output$seq %>% as.character, 'FALSE')
		})

		expected_output = tibble(identifier=c(1,1,2),tax_id.x=c(as.integer(10090),as.integer(10090),as.integer(10090)),miRNA=c('mmu-miR-4-5p','mmu-miR-4b-5p','mmu-miR-1-3p'),
					miRNA_sequence=c('CAAAAAAAGCGCGCGCGCGT','CAAAAAAAGCGCGCGTTCGT','ACCTCTCGCGAGAGAGGCGC'))
		test_that('expected output matches actual output', {
			expect_equal(expected_output, output)
		})

	}
}

run_tests('hsa_mirna_families.tsv', 'hsa_mature_mirnas.tsv', 'hsa')
run_tests('mmu_mirna_families.tsv','mmu_mature_mirnas.tsv','mmu')
