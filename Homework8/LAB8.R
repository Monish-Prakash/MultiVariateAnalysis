library(readr)
library("FactoMineR")
library("arules")


wd = getwd()
if(grepl("charlyo", wd)) {
    setwd("~/Desktop/MVA/Practiques")
} else {
    setwd("~/Documents/UPC-MASTER/2nd Semester/MVA/Homework8")
}
Sys.setlocale("LC_TIME", "en_US")  # OS X, in UTF-8

# 1
tic_tt <- read.csv("tic_tt.csv", header = TRUE, sep =";", stringsAsFactors = TRUE, colClasses=c("NULL", rep('factor', 33)))

# 2
tic_tt.catdes = catdes(as.data.frame(tic_tt), num.var = 28)

# 3
ttr <- as(tic_tt, "transactions")
summary(ttr)

# 4
min_support <- 0.1
min_confidence <- 0.7
max_size <- 15
fis <- apriori (ttr, parameter = list(support=min_support, confidence=min_confidence, maxlen = max_size))

# 5
fsets <- unique(generatingItemsets(fis))
sor.fsets <- sort(fsets, by = "support")
inspect(sor.fsets[1:10, ])

# 6
sor.fis <- sort(fis, by = "lift")
inspect(sor.fis[1:10, ])

# 7
fis.payments <- subset(fis, subset = rhs %in% c("Pagament.a.trav.s.d.Internet.=TRUE","Pagament.a.trav.s.d.Internet.=FALSE") )
sor.fis.payments <- sort(fis.payments, by = "lift")
inspect(sor.fis.payments[1:10, ])
