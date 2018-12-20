library(filtar)

context('test get_mirna_family')

run_test = function(input, species ) {

        output = get_mirna_family(input,species)

        test_that('id column contains 1-4 numbers', {
                expect_match(as.character(output$identifier), '^[0-9]{1,4}$')
        })

        test_that('Seed sequences contain 7 repetitions of allowed characters', {
                expect_match( as.character(output$seq), '^[ACGTU]{7}$' )
        })

        test_that('tax_id contain 4-6 repetitions of allowed characters', {
                expect_match( output$tax_id, '^[0-9]{4,6}$' )
        })

        test_that('The first column is monotonic', {
                expect_true(all(output$identifier == cummax(output$identifier)))
        })

        test_that('The combined seq and tax_id columns contain no duplicates', {
                expect_match(output[,c('seq','tax_id')] %>% duplicated %>% as.character, 'FALSE')
        })

        tmp =  output[,c('seq','tax_id')] %>% dplyr::filter(seq=='CCTCTCG')
        tmp2 = c('9606','10090') %in% tmp$tax_id

        test_that('The output factorises orthologous seeds as belonging to the same family', {
                expect_match(tmp2 %>% as.character, '^TRUE$')
        })

        if (species == 'hsa') {

                test_that('The output does not contain seed sequences without a reference orthologue', {
                        expect_equal(dim(output %>% dplyr::filter(seq=='AAAAAAA'))[1], 0)
                })

        } else if (species == 'mmu') {

                test_that('The output does not contain seed sequences without a reference orthologue', {
                        expect_equal(dim(output %>% dplyr::filter(seq=='AGCGCGC'))[1], 0)
                })
        }
}

run_test('mirna_seeds.tsv','hsa')
run_test('mirna_seeds.tsv', 'mmu')
