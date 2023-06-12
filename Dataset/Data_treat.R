
library(xlsx)

data<-read.xlsx("LTE_data.xlsx", sheetName = "SOC")


library(reshape2)
melt(data, variable.name = "site", value.name = "SOC", year ="year")
colnames(data)[3:17]=seq(1:15)

melted_data<-melt(data, id.vars = c('Site', 'year'))

write.xlsx(melted_data, file = "LTE_data.xlsx", sheetName = "SOC3", append=T)
