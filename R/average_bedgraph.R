
#' Take a merged bedgraph dataframe, and return a dataframe with the mean average bedgraph coverage for each record
#'
#' @param df A merged bedgraph dataframe in which each record represent one locus on a give chromosome, and columns represent coverage values from different read mapping analyses
#' @return df A bedgraph dataframe with the artihmetic mean average coverage for each records locus
#' @export

AvgBedgraph = function (df) {

avg = vector(length=dim(df)[1])

# get the average of col_4 ... col_n

cov_values = df[,4:dim(df)[2]]
cov_values = t(cov_values)

for (i in 1:dim(cov_values)[2]) {

        # only retrieve the bed record columns
        avg[i] = mean(cov_values[,i])
}

df$avg = avg

return(df)
}
