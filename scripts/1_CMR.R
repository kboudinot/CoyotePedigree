#packages
library(tidyverse)
library(RMark) #also need to have MARK installed on computer to interface with R

#load data
genotypes_22 <- read.csv("genotypes_22.csv")
genotypes_23 <- read.csv("genotypes_23.csv")

#2022 
#capture data
df <- import.chdata("RMARK_22.txt")
#model
null.model <- mark(df, model="Huggins")
summary(null.model)
null.model$results$derived # derived population size estimate 2022

#2023 
#load capture data
df2 <- import.chdata("RMARK_23.txt")
#model
null.model2 <- mark(df2, model="Huggins")
summary(null.model2)
null.model2$results$derived # derived population size estimate 2023
