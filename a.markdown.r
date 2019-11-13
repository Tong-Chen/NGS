data <- read.table("a", sep="\t", check.names=F, quote="", comment="", header=TRUE, row.names=NULL, nrows=3)

if (1 > 0){
	data <- data[, 1:1, drop=F]
}

a = knitr::kable(data, format="markdown")

write(a, file="a.markdown.table")

