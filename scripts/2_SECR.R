#packages
library(tidyverse)
library(secr)
library(sf)

#load data
genotypes_22 <- read.csv("genotypes_22.csv")
genotypes_23 <- read.csv("genotypes_23.csv")

#coordinates of center point of each site (60)
coords <- data.frame(x = c(-69.9695280359307,
                           -70.0712271060513,
                           -70.018555871455,
                           -70.0890886801604,
                           -69.9848551835505,
                           -70.0602073590689,
                           -69.9902521543908,
                           -70.2067263516777,
                           -69.953934,
                           -70.050230547121,
                           -70.2031351816076,
                           -69.9984550447941,
                           -69.979132903413,
                           -70.0732003651383,
                           -70.0236950796058,
                           -70.0874575429285,
                           -69.9759421452131,
                           -70.0489714038199,
                           -70.2147520682651,
                           -69.9634048600708,
                           -70.0419550850144,
                           -70.1914043661973,
                           -70.068452807628,
                           -70.0017246128436,
                           -70.1096602389818,
                           -70.1680181863967,
                           -69.971377124273,
                           -70.0569320341383,
                           -70.0134001902317,
                           -69.9766299931947,
                           -70.058942763424,
                           -69.9609186186687,
                           -70.0215523108671,
                           -70.2097108642592,
                           -70.0746584831957,
                           -70.0050057633588,
                           -70.0984358840267,
                           -69.9746595358165,
                           -70.0634818777151,
                           -70.0068574760567,
                           -70.2240452775397,
                           -69.9767333314636,
                           -70.0536559026333,
                           -70.0386871611552,
                           -70.2290203595527,
                           -69.9662595109533,
                           -69.9977123403749,
                           -70.1821616483663,
                           -70.1947268318705,
                           -70.0667477015964,
                           -70.0152749506544,
                           -70.1122277308586,
                           -70.1442268123013,
                           -69.9825958261874,
                           -69.9933129614499,
                           -69.966967208441,
                           -70.0386768188774,
                           -70.2214055048765,
                           -70.070053,
                           -70.0269754809172),
                     y = c(41.9013294094054,
                           41.9597275361084,
                           41.967762253221,
                           42.0473940483435,
                           41.9236964822779,
                           41.9577819412466,
                           41.9432188302937,
                           42.0418194446764,
                           41.841784,
                           41.9832284330378,
                           42.0788696816995,
                           41.954321319523,
                           41.9041711809366,
                           41.9514337507092,
                           41.987390387768,
                           42.0537544002077,
                           41.9208459901205,
                           42.0204507123857,
                           42.0461749648428,
                           41.819638826142,
                           41.9725873270349,
                           42.0773980909623,
                           41.9088124490749,
                           41.945719835656,
                           42.0605256084237,
                           42.0727439204978,
                           41.8930248299205,
                           41.9663901113588,
                           41.9844605185844,
                           41.9126037199558,
                           42.0232208649736,
                           41.8538609371961,
                           41.9955470259118,
                           42.0627141741754,
                           41.9436516488233,
                           41.9371132726037,
                           42.0572533279191,
                           41.8844195567129,
                           41.9491737739248,
                           41.9653481370603,
                           42.0566554015911,
                           41.8765216151937,
                           41.974998284461,
                           42.0173496278006,
                           42.0749089735494,
                           41.8733967447361,
                           41.9621773203485,
                           42.0243901683629,
                           42.0689298022074,
                           41.9405687853829,
                           41.9763694659824,
                           42.0532888413598,
                           42.0682266790484,
                           41.9319812675149,
                           41.9346980150766,
                           41.8367297966292,
                           41.9811951822504,
                           42.0651220661978,
                           41.929345,
                           41.9787828506225)
)

#construct sf object
latlon <- st_as_sf(coords, coords = 1:2)

# specify initial CRS: WGS84 lat-lon
st_crs(latlon) <- 4326
st_coordinates(latlon)

# project to Cartesian coordinate system units meters
# EPSG:32619 is the WGS 84/UTM zone 19N code
trps <- st_transform(latlon, crs = 32619)
st_coordinates(trps)

## Read in capture file and trap layout file to get capthist object
capthist <- secr::read.capthist(captfile="secrcaps.txt",
                                trapfile = "secrtraps.txt",
                                detector="count",
                                fmt="trapID")
spacing <-  937.1729
buffer <- suggest.buffer(capthist)

# loading shapefile for model mask
library(sf)
ccns.geo <- st_read("Boundary/CapeCodNationalSeashoreBoundary.shp")
st_set_crs(ccns.geo, value = 4326)# setting crs to WGS84
ccns.utm <- st_transform (ccns.geo, crs = 32619) 

# make mask 
mask2 <- make.mask(traps(capthist), 
                   buffer = buffer, # 3978 from suggest.buffer
                   spacing = spacing, # avg spacing of centroids from capthist summary
                   type = "polygon",
                   poly = ccns.utm)
summary(mask2)
plot(mask2); plot(traps(capthist), add=T)

# model 2022
m1 <- secr.fit(capthist,
               model = list(D ~ 1, g0 ~ 1, sigma ~ 1),
               mask = mask2)
summary(m1)
predict(m1)

# 2023
capthist2 <- read.capthist(captfile="secrcaps23.txt",
                           trapfile = "secrtraps.txt",
                           detector="count",
                           fmt="trapID")

buffer2 <- suggest.buffer(capthist2)

mask23 <- make.mask(traps(capthist2), 
                    buffer = buffer2, 
                    spacing = spacing, 
                    type = "polygon",
                    poly = ccns.utm)

m2 <- secr.fit(capthist2,
               model = list(D ~ 1),
               mask = mask23)
summary(m2)
predict(m2)


# Rarefaction----
library(iNEXT) #load library
inext_data <- read.csv("iNEXT_data.csv") #read in data 
inext_data <- inext_data[,-1] #wrangle
test <- iNEXT(inext_data, #run iNEXT
              datatype="abundance")
plot(test, #plot with standard errors 
     type=1,
     se=T,
     show.legend=F,
     show.main = F,
     col = "black")
