context('test fix_ts_output')

output_canonical=readr::read_tsv('example_targetscan_sites.tsv') # this is output from the unpatched script
output = fix_ts_output('example_targetscan_sites_intermediate.tsv')            #('tmp.tsv') # this is output from patched script

output_canonical[is.na(output_canonical)] <- ""

test_that('Canonical and patched script output is identical',{
        expect_equal(output, output_canonical)
})
