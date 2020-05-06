library(qs)
popdem=readRDS("../../covid-uk/covidm/data/wpp2019_pop2020.rds");
write.table(popdem[popdem$name=="Portugal",],file="portugal.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="Switzerland",],file="switzerland.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="Poland",],file="poland.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="United Kingdom",],file="uk.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="Slovenia",],file="slovenia.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="Sweden",],file="sweden.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="Spain",],file="spain.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="Germany",],file="germany.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)
write.table(popdem[popdem$name=="France",],file="france.agegroups",sep="\t",col.names=TRUE,row.names=FALSE)

