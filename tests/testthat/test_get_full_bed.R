library(plyr)
library(tidyverse)
library(filtar)

context('test get_full_bed')

x = get_full_bed('example_ensembl.bed','example_APAtrap.bed','hsa_all_transcripts.txt')

y = x %>% filter(id=='ENST00000224237.9')

test_that("ENST00000224237.9 has the correct start codon annotation", {
                expect_equal(17237230, y$start[1], )
})

test_that("ENST00000224237.9 has the correct end codon annotation", {
                expect_equal(17237588, y$stop[1], )
})

y = x %>% filter(id=='ENST00000464969.6')

test_that("ENST00000464969.6 has the correct start codon annotation", {
                expect_equal(102157494, y$start[1], )
})

test_that("ENST00000464969.6 has the correct end codon annotation", {
                expect_equal(102162599, y$stop[10], )
})


y = x %>% filter(id=='ENST00000298510.3')

test_that("ENST00000298510.3 has the correct start codon annotation", {
                expect_equal(119167702, y$start[1], )
})

test_that("ENST00000298510.3 has the correct end codon annotation", {
                expect_equal(119168533, y$stop[1], )
})

y = x %>% filter(id=='ENST00000463743.5')

test_that("ENST00000463743.5 has the correct start codon annotation", {
                expect_equal(93323077, y$start[1], )
})

test_that("ENST00000463743.5 has the correct end codon annotation", {
                expect_equal(93307001, y$stop[7], )
})

y = x %>% filter(id=='ENST00000363306.1')

test_that("ENST00000363306.1 has the correct start codon annotation", {
                expect_equal(86889569, y$start[1], )
})

test_that("ENST00000363306.1 has the correct end codon annotation", {
                expect_equal(86889682, y$stop[1], )
})

y = x %>% filter(id=='ENST00000378952.7')

test_that("ENST00000378952.7 has the correct start codon annotation", {
                expect_equal(12167673, y$start[1], )
})

test_that("ENST00000378952.7 has the correct end codon annotation", {
                expect_equal(12167811, y$stop[1], )
})

y = x %>% filter(id=='ENST00000381604.8')

test_that("ENST00000378952.8 has the correct start codon annotation", {
                expect_equal(252470, y$start[1], )
})

test_that("ENST00000378952.8 has the correct end codon annotation", {
                expect_equal(254626, y$stop[1], )
})

y = x %>% filter(id=='ENST00000474119.5')

test_that("ENST00000474119.5 has the correct start codon annotation", {
                expect_equal(4847234, y$start[1], )
})

test_that("ENST00000474119.5 has the correct end codon annotation", {
                expect_equal(4848062, y$stop[2], )
})

y = x %>% filter(id=='ENST00000564130.2')

test_that("ENST00000564130.2 has the correct start codon annotation", {
                expect_equal(46891, y$start[1], )
})

test_that("ENST00000564130.2 has the correct end codon annotation", {
                expect_equal(47056, y$stop[1], )
})

y = x %>% filter(id=='ENST00000567466.1')

test_that("ENST00000567466.1 has the correct start codon annotation", {
                expect_equal(48424, y$start[1], )
})

test_that("ENST00000567466.1 has the correct end codon annotation", {
                expect_equal(48114, y$stop[2], )
})

test_that('bed file has the correct number of columns', {
	expect_equal(dim(x)[2], 5)
})

normal_bed = read_tsv('example_ensembl.bed', col_names=c('chromosome','start','stop','strand','id'))
extended_utrs = read_tsv('example_APAtrap.bed', col_names=c('chromosome','start','stop','id','dummy','strand'))

normal_bed = separate(normal_bed, id, into=c('id','version'))
extended_utrs = separate(extended_utrs, id, into=c('id','dummy2','chrom_dup','strand_dup'))

normal_bed$version = NULL

extended_utrs = extended_utrs[,c('chromosome','start','stop','strand','id')]
extended_utrs$chromosome = stringr::str_replace_all(extended_utrs$chromosome, 'chr','')
normal_bed$chromosome = stringr::str_replace_all(normal_bed$chromosome, 'chr','')

old_records = extended_utrs[extended_utrs$id %in% normal_bed$id,]

test_that('bed file has the correct number of rows', {
        expect_equal(dim(x)[1], dim(normal_bed)[1] + dim(extended_utrs)[1] - dim(old_records)[1])
})
