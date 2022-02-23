## A draft app for visualization of climate change in defined areas
## author: Colin Mahony colin.mahony@gov.bc.ca

# Copyright 2021 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(shiny)
library(RColorBrewer)
library(stinepack) # for interpolation splines
library(scales)
library(plotrix)
library(car)
library(DT)
library(rgdal)
library(mapview)
library(leaflet)
library(leafem)
library(htmlwidgets)
library(leaflet.opacity)
library(sf)
library(raster)

studyarea <- "RKB"

edatopes <- c("B2", "C4", "D6")
edatope.names <- c("Poor-subxeric", "Medium-mesic", "Rich-hygric")

scenarios <- c("ssp126", "ssp245", "ssp370")
scenario.names=c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")

proj.year.names=c("2001-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")

bdy <- readOGR(dsn = paste("bdy/bdy", studyarea, "shp", sep="."))
P4S.epsg <- CRS ("+init=epsg:4326") # web mercator
bdy <- spTransform(bdy, P4S.epsg)

bbox <- as.vector(unlist(extent(bdy)))

options(scipen=999)
## Load the input data

modelMetadata <- read.csv("data/ModelList.csv")
# str(ecoprov.climate)

## CLIMATE DATA

clim.meanChange <- read.csv(paste("data/clim.meanChange", studyarea, "csv", sep="."))[,-c(1:3)]
clim.refMean <- read.csv(paste("data/clim.refMean", studyarea, "csv", sep="."))

variables <- names(clim.meanChange)
variable.names <- read.csv("data/Variables_ClimateBC.csv")
variables.select <- variables[c(grep("_wt|_sp|_sm|_at", variables), 225:272)]
variables.select <- variables.select[-grep("RH|Rad|MAR", variables.select)]

variable.types <- rep(NA, length(variables))
variable.types[grep("PPT|DD|PAS|NFFD|Eref|FFP|CMD|MAP|MSP|AHM|SHM|Rad|MAR", variables)] <- "ratio"
variable.types[grep("Tmax|Tmin|Tave|MAT|MWMT|MCMT|TD|EMT|EXT|bFFP|eFFP|CMI", variables)] <- "interval"
variable.types[grep("RH", variables)] <- "pct"

# Convert absolute change to relative change for zero-limited variables
clim.meanChange.ratio <- clim.meanChange
for(i in which(variable.types=="ratio")){clim.meanChange.ratio[,i] <- clim.meanChange.ratio[,i]/clim.refMean[,i]}  #set zero values to one, to facilitate log-transformation

period <- read.csv(paste("data/clim.meanChange", studyarea, "csv", sep="."), stringsAsFactors = F)[,c(1:3)]
period[,3] <- substr(period[,3],1,4)

scenarios <- sort(unique(period[,2]))
proj.years <- unique(period[,3])
gcms <- unique(period[,1])[-1]
mods <- c("ACC", "BCC", "CNR", "EC", "GFD", "GIS", "MIR", "MPI", "MRI", "ENS")

ColScheme.gcms <- c(brewer.pal(n=length(gcms), "Paired"))

## BGC PROJECTIONS

bgc.names <- read.csv(paste("data/All_BGCs_v11_21.csv", sep="."), stringsAsFactors = F)

bgcs.native <- as.vector(unlist(read.csv(paste("data/bgcs.native", studyarea, "csv", sep="."))))
zones.native <- as.vector(unlist(read.csv(paste("data/zones.native", studyarea, "csv", sep="."))))

## projected area of biogeoclimatic units
bgc.area <- read.csv(paste("data/PredSum.BGC", studyarea, "csv", sep="."))
temp <- bgc.area[-which(bgc.area$GCM=="ensemble"),] #remove ensemble vote
temp.ens <- bgc.area[which(bgc.area$GCM=="ensemble"),] #ensemble vote
bgc.area <- rbind(temp, temp.ens)
bgc.area <- bgc.area[,-c(1:3)] # matches the "period" table
bgc.area[is.na(bgc.area)] <- 0
bgcs.all <- names(bgc.area)
r <- raster(paste("data/BGC.pred", studyarea, "ref.tif", sep="."))
cellarea <- (res(r)[2]*111)*(res(r)[1]*111*cos(mean(extent(r)[3:4])))
bgc.area <- bgc.area*cellarea

#BGC subzone color scheme
bgccolors <- read.csv("data/WNAv11_Subzone_Colours.csv", stringsAsFactors = F)
zonecolors <- read.csv("data/WNAv11_Zone_Colours.csv", stringsAsFactors = F)
zonecolors.BC <- read.csv("data/BGCzone_Colorscheme.csv", stringsAsFactors = F)
zonecolors$colour[match(zonecolors.BC$zone, zonecolors$classification)] <- as.character(zonecolors.BC$HEX)

## projected area of biogeoclimatic zones
bgc.zones <- bgcs.all
for(i in zonecolors$classification){ bgc.zones[grep(i,bgcs.all)] <- i }
zones <- unique(bgc.zones)
zone.area <- data.frame(matrix(NA, length(bgc.area[,1]), length(zones)))
names(zone.area) <- zones
for(zone in zones){ 
  s <- which(bgc.zones==zone)
  zone.area[, which(zones==zone)] <- if(length(s)>1) apply(bgc.area[,s], 1, sum, na.rm=T) else bgc.area[,s]
}
zones.all <- names(zone.area)

## add colors for units not in color scheme. 
needcolor <- bgcs.all[-which(bgcs.all%in%bgccolors$classification)]
needcolor.zone <- needcolor
for(i in zonecolors$classification){ needcolor.zone[grep(i,needcolor)] <- i }
bgccolors <- rbind(bgccolors, data.frame(classification=needcolor, colour=zonecolors$colour[match(needcolor.zone, zonecolors$classification)]))


## simplify and structure the biogeoclimatic area tables
bgc.area <- bgc.area[,rev(order(apply(bgc.area, 2, sum)))] #sort by total projection area
totalarea <- sum(bgc.area[1,]) #historical distribution 
bgc.area <- bgc.area[,-which(apply(bgc.area, 2, sum)/totalarea < 0.25)] #remove small units
bgc.area <- bgc.area[,order(names(bgc.area))] #sort by total projection area
bgcs <- names(bgc.area)

zone.area <- zone.area[,rev(order(apply(zone.area, 2, sum)))] #sort by total projection area
totalarea <- sum(zone.area[1,]) #historical distribution 
zone.area <- zone.area[,-which(apply(zone.area, 2, sum)/totalarea < 0.5)] #remove small units
zones <- names(zone.area)

# biogeoclimatic spatial data
bgc.pred.ref <- raster(paste("data/BGC.pred", studyarea, "ref.tif", sep="."))
bgc.pred.2001 <- raster(paste("data/BGC.pred", studyarea, "hist.2001.tif", sep="."))
gcm="ensemble"
scenario="ssp245"
for(gcm in gcms){
  for(proj.year in proj.years[-1]){
    assign(paste("bgc.pred", gcm, scenario, proj.year, sep="."), raster(paste("data/BGC.pred", studyarea, gcm, scenario, proj.year, "tif", sep=".")))
  }
}
levels.bgc <- read.csv("data/levels.bgc.csv")[,1]

## persistence and expansion of bgc units

    #bgc total table
    bgc.count <- read.csv(paste("data/PredSum.BGC", studyarea, "csv", sep="."))
    bgc.count[is.na(bgc.count)] <- 0
    temp <- bgc.count[-which(bgc.count$GCM=="ensemble"),] #remove ensemble vote
    temp.ens <- bgc.count[which(bgc.count$GCM=="ensemble"),] #ensemble vote
    bgc.count <- rbind(temp, temp.ens)
    bgc.count <- bgc.count[,-c(1:3)] # matches the "period" table

    #zone total table
    zone.count <- read.csv(paste("data/PredSum.zone", studyarea, "csv", sep="."))
    zone.count[is.na(zone.count)] <- 0
    temp <- zone.count[-which(zone.count$GCM=="ensemble"),] #remove ensemble vote
    temp.ens <- zone.count[which(zone.count$GCM=="ensemble"),] #ensemble vote
    zone.count <- rbind(temp, temp.ens)
    zone.count <- zone.count[,-c(1:3)] # matches the "period" table

    #bgc table for home range
    bgc.count.home <- read.csv(paste("data/PredSum.BGC.home", studyarea, "csv", sep="."))
    bgc.count.home[is.na(bgc.count.home)] <- 0
    temp <- bgc.count.home[-which(bgc.count.home$GCM=="ensemble"),] #remove ensemble vote
    temp.ens <- bgc.count.home[which(bgc.count.home$GCM=="ensemble"),] #ensemble vote
    bgc.count.home <- rbind(temp, temp.ens)
    bgc.count.home <- bgc.count.home[,-c(1:3)] # matches the "period" table

    #zone table for home range
    zone.count.home <- read.csv(paste("data/PredSum.zone.home", studyarea, "csv", sep="."))
    zone.count.home[is.na(zone.count.home)] <- 0
    temp <- zone.count.home[-which(zone.count.home$GCM=="ensemble"),] #remove ensemble vote
    temp.ens <- zone.count.home[which(zone.count.home$GCM=="ensemble"),] #ensemble vote
    zone.count.home <- rbind(temp, temp.ens)
    zone.count.home <- zone.count.home[,-c(1:3)] # matches the "period" table

    ## calculate persistence and expansion tables
    bgc.persistence <- sweep(bgc.count.home, MARGIN=2,unlist(bgc.count[1,]), '/' )
    zone.persistence <- sweep(zone.count.home, MARGIN=2,unlist(zone.count[1,match(names(zone.count.home), names(zone.count))]), '/' )
    bgc.expansion <- sweep(bgc.count-bgc.count.home, MARGIN=2,unlist(bgc.count[1,]), '/' )
    zone.expansion <- sweep(zone.count[,match(names(zone.count.home), names(zone.count))]-zone.count.home, MARGIN=2,unlist(zone.count[1,match(names(zone.count.home), names(zone.count))]), '/' )

    ## simplify and structure the persistence and expansion tables
    totalarea <- sum(bgc.count[1,], na.rm = T) #historical distribution
    small <- which(bgc.count[1,]/totalarea < 0.01) # establish insignificant species for removal
    small <- c(small, which(is.na(bgc.count[1,]))) # establish insignificant species for removal
    bgc.persistence <- bgc.persistence[,-small] #remove small units and assign to permanent table
    bgc.expansion <- bgc.expansion[,-small] #remove small units and assign to permanent table


    
# ## SPECIES FEASIBILITIES
# SiteLookup <- read.csv("data/SiteLookup.csv", stringsAsFactors = F)
# SuitLookup <- read.csv("data/SuitLookup.csv", stringsAsFactors = F)
# 
# spps.all <- vector()
# spps.native <- vector()
# edatope="C4"
# for(edatope in edatopes){
#   
#   #fractional feasibility table
#   suit.area <- read.csv(paste("data/PredSum.suit", studyarea, edatope, "csv", sep="."))
#   temp <- suit.area[-which(suit.area$GCM=="ensemble"),] #remove ensemble vote
#   temp.ens <- suit.area[which(suit.area$GCM=="ensemble"),] #ensemble vote
#   suit.area <- rbind(temp, temp.ens)
#   suit.area <- suit.area[,-c(1:3)] # matches the "period" table
#   
#   #species table
#   spp.area <- read.csv(paste("data/PredSum.spp", studyarea, edatope, "csv", sep="."))
#   temp <- spp.area[-which(spp.area$GCM=="ensemble"),] #remove ensemble vote
#   temp.ens <- spp.area[which(spp.area$GCM=="ensemble"),] #ensemble vote
#   spp.area <- rbind(temp, temp.ens)
#   spp.area <- spp.area[,-c(1:3)] # matches the "period" table
#   
#   #fractional feasibility table for home range
#   suit.area.home <- read.csv(paste("data/PredSum.suit.home", studyarea, edatope, "csv", sep="."))
#   temp <- suit.area.home[-which(suit.area.home$GCM=="ensemble"),] #remove ensemble vote
#   temp.ens <- suit.area.home[which(suit.area.home$GCM=="ensemble"),] #ensemble vote
#   suit.area.home <- rbind(temp, temp.ens)
#   suit.area.home <- suit.area.home[,-c(1:3)] # matches the "period" table
#   
#   #species tablefor home range
#   spp.area.home <- read.csv(paste("data/PredSum.spp.home", studyarea, edatope, "csv", sep="."))
#   temp <- spp.area.home[-which(spp.area.home$GCM=="ensemble"),] #remove ensemble vote
#   temp.ens <- spp.area.home[which(spp.area.home$GCM=="ensemble"),] #ensemble vote
#   spp.area.home <- rbind(temp, temp.ens)
#   spp.area.home <- spp.area.home[,-c(1:3)] # matches the "period" table
#   
#   ## spp color scheme
#   colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
#   colors = colors[-grep("yellow", colors)]
#   set.seed(5)
#   sppcolors <- c(brewer.pal(n=12, "Paired")[-11],sample(colors,dim(spp.area)[1]-11)) # removal of "11" is light yellow, doesn't show up well. 
#   
#   ## calculate persistence and expansion tables
#   suit.persistence <- sweep(suit.area.home, MARGIN=2,unlist(suit.area[1,]), '/' )
#   spp.persistence <- sweep(spp.area.home, MARGIN=2,unlist(spp.area[1,]), '/' )
#   suit.expansion <- sweep(suit.area-suit.area.home, MARGIN=2,unlist(suit.area[1,]), '/' )
#   spp.expansion <- sweep(spp.area-spp.area.home, MARGIN=2,unlist(spp.area[1,]), '/' )
#   
#   ## simplify and structure the area tables
#   totalarea <- sum(suit.area[1,]) #historical distribution 
#   small <- which(apply(suit.area, 2, sum, na.rm=T)/totalarea < 0.5) # establish insignificant species for removal
#   assign(paste("suit.area", edatope, sep="."), suit.area[,-small]) #remove small units and assign to permanent table
#   assign(paste("spp.area", edatope, sep="."), spp.area[,-small]) #remove small units and assign to permanent table
#   assign(paste("spps.all", edatope, sep="."), names(spp.area)[-small])
#   
#   ## simplify and structure the persistence and expansion tables
#   totalarea <- sum(spp.area[1,]) #historical distribution 
#   exotic <- which(spp.area[1,]/totalarea < 0.01) # establish insignificant species for removal
#   assign(paste("suit.persistence", edatope, sep="."), suit.persistence[,-exotic]) #remove exotic units and assign to permanent table
#   assign(paste("spp.persistence", edatope, sep="."), spp.persistence[,-exotic]) #remove exotic units and assign to permanent table
#   assign(paste("suit.expansion", edatope, sep="."), suit.expansion[,-exotic]) #remove exotic units and assign to permanent table
#   assign(paste("spp.expansion", edatope, sep="."), spp.expansion[,-exotic]) #remove exotic units and assign to permanent table
#   assign(paste("spps.native", edatope, sep="."), names(spp.area)[-exotic])
#   
#   spps.all <- c(spps.all, names(spp.area)[-small])
#   spps.native <- c(spps.native, names(spp.area)[-exotic])
#   
# }
# 
# spps.all <- unique(spps.all)[order(unique(spps.all))]
# spps.native <- unique(spps.native)[order(unique(spps.native))]
# 
# ## Color Schemes for species change maps
# breakpoints.suit <-   breakseq <- c(0.5,1.5,2.5,3.5)
# palette.suit <-   c("#006400", "#1E90FF", "#EEC900")
# ColScheme.suit <- colorBin(palette.suit, bins=breakpoints.suit, na.color = NA)
# breakpoints.change <- seq(-3,3,0.5)
# palette.change <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,5,6)], brewer.pal(11,"RdBu")[c(6,7,8,9,10,11)])
# ColScheme.change <- colorBin(palette.change, bins=breakpoints.change, na.color = NA)
# labels.change <- breakpoints.change[-median(1:length(breakpoints.change))]
# breakpoints.binary <- seq(-1,1,0.2)
# palette.binary <- c(brewer.pal(11,"RdBu")[c(1:4,6)], brewer.pal(11,"RdBu")[c(6,8:11)])
# ColScheme.binary <- colorBin(palette.binary, bins=breakpoints.binary, na.color = NA)
# labels.binary <- paste(abs(seq(-.9,.9,0.2))*100, "%", sep="")

## SPATIAL DATA
bgc.simple <- st_read("data/bgc.simple.shp")
bgc.maprecord <- as.character(bgc.simple$BGC)
zone.maprecord <- bgc.maprecord
for(i in zonecolors$classification){ zone.maprecord[grep(i,bgc.maprecord)] <- i }
bgc.list <- sort(unique(bgc.maprecord))
zone.list <- sort(unique(zone.maprecord))

# Define UI ----
ui <- fluidPage(
  navbarPage(title = "Climate change summary", theme = "bcgov.css", 
             tabPanel("App", 
                      fluidRow(
                        column(2,
                               helpText("Use this app to explore projected changes in climate climate conditions, expressed as shifts in biogeoclimatic units. All data are from ClimateBC. See the 'Model Info' tab for model names"),
                               
                               tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            });
                            $(window).resize(function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            });
                            ')),
                               
                               
                               radioButtons("type", inline = FALSE, 
                                            label = "Choose the type of map",
                                            choices = list("Climate variables" = 1, "Biogeoclimatic units" = 2),
                                            # choices = list("Climate variables" = 1, "Biogeoclimatic units" = 2, "Species feasibility" = 3),
                                            selected = 2),
                               
                               sliderInput("transparency", label = "Layer transparency", min = 0, 
                                           max = 1, value = 0.7),
                               
                               radioButtons("maptype",
                                            label = "Choose a time period",
                                            choices = list("Reference (1961-1990)" = 1, "Recent (2001-2020)" = 2, "Future" = 3),
                                            selected = 1),
                               
                               conditionalPanel(
                                 condition = "input.zonelevel == true",
                                 h4("Reference biogeoclimatic zone map"),
                                 img(src = paste("refmap",studyarea, "zones.png", sep="."), height = 669*0.45, width = 661*0.45)
                               ),
                               
                               conditionalPanel(
                                 condition = "input.zonelevel == false",
                                 h4("Reference biogeoclimatic variant map"),
                                 img(src = paste("refmap",studyarea, "variants.png", sep="."), height = 669*0.45, width = 661*0.45)
                               )
                        ),    
                        
                        column(5, 
                               
                               leafletOutput(outputId = "map", height="86vh"),
                               downloadButton(outputId = "downloadMap", label = "Download map")
                               
                        ),
                        
                        column(5, 
                               column(6, 
                                      
                                      selectInput("var1", 
                                                  label = "Choose the x-axis variable",
                                                  choices = as.list(variables.select),
                                                  selected = "MAT"),
                                      
                                      selectInput("gcm.focal", 
                                                  label = "Choose a climate model",
                                                  choices = as.list(gcms),
                                                  selected = "ensemble"),
                                      
                                      radioButtons("proj.year", inline = TRUE,
                                                   label = "Choose a future period",
                                                   choices = list("2001-2020" = 1, "2021-2040" = 2, "2041-2060" = 3, "2061-2080" = 4, "2081-2100" = 5),
                                                   # choices = list(proj.year.names[1] = 1, proj.year.names[2] = 2, proj.year.names[3] = 3, proj.year.names[4] = 4, proj.year.names[5] = 5),
                                                   selected = 3),
                                      
                                      checkboxInput("recent", label = "Show recent observed climate (2001-2020)", value = T),
                                      
                               ),
                               column(6, 
                                      
                                      conditionalPanel(
                                        condition = "input.type == 1",
                                        
                                        selectInput("var2", 
                                                    label = "Choose the y-axis variable",
                                                    choices = as.list(variables.select),
                                                    selected = "MAP"),
                                        
                                        checkboxInput("ratioscale", label = "Relative (%) scale for ratio variables", value = T),
                                        
                                      ),
                                      
                                      conditionalPanel(
                                        condition = "input.type == 2",
                                        
                                        checkboxInput("zonelevel", label = "Generalize to BGC zone level", value = T),
                                        
                                        radioButtons("plotbgc", inline = TRUE,
                                                     label = "Choose a plot type",
                                                     choices = list("Area" = 1, "Persistence" = 2),
                                                     selected = 1),
                                        
                                        conditionalPanel(
                                          condition = "input.plotbgc == 1",
                                          
                                          conditionalPanel(
                                            condition = "input.zonelevel == true",
                                            
                                            radioButtons("zone.area.focal", inline = TRUE, 
                                                         label = "Select BGC zone for ensemble detail",
                                                         choices = as.list(c("none", zones)),
                                                         selected = "none"),
                                            
                                          ),
                                          
                                          conditionalPanel(
                                            condition = "input.zonelevel == false",
                                            
                                            selectInput("bgc.area.focal", 
                                                        label = "Select BGC subzone-variant for ensemble detail",
                                                        choices = as.list(c("none", bgcs)),
                                                        selected = "none"),
                                            
                                          ),
                                          
                                        ),
                                        
                                        conditionalPanel(
                                          condition = "input.plotbgc == 2",
                                          
                                          conditionalPanel(
                                            condition = "input.zonelevel == true",
                                            
                                            radioButtons("zone.persistence.focal", inline = TRUE, 
                                                         label = "Select BGC zone for ensemble detail",
                                                         choices = as.list(c("none", zones.native)),
                                                         selected = "none"),
                                            
                                          ),
                                          
                                          conditionalPanel(
                                            condition = "input.zonelevel == false",
                                            
                                            selectInput("bgc.persistence.focal", 
                                                        label = "Select BGC subzone-variant for ensemble detail",
                                                        choices = as.list(c("none", bgcs.native)),
                                                        selected = "none"),
                                            
                                          ),
                                          
                                        ),
                                        
                                      ),
                                      
                                      # conditionalPanel(
                                      #   condition = "input.type == 3",
                                      #   
                                      #   radioButtons("mapspp", inline = TRUE,
                                      #                label = "Choose a map type",
                                      #                choices = list("Feasibility" = 1, "Change" = 2, "Loss/gain" = 3),
                                      #                selected = 1),
                                      #   
                                      #   radioButtons("plotspp", inline = TRUE,
                                      #                label = "Choose a plot type",
                                      #                choices = list("Area" = 1, "Persistence" = 2),
                                      #                selected = 1),
                                      #   
                                      #   checkboxInput("fractional", label = "Use fractional (partial) feasibilities", value = F),
                                      #   
                                      #   radioButtons("edatope", inline = TRUE, 
                                      #                label = "Select an edatope (site type)",
                                      #                choices = as.list(edatopes),
                                      #                selected = edatopes[2]),
                                      #   
                                      #   conditionalPanel(
                                      #     condition = "input.plotspp == 1",
                                      #     
                                      #     conditionalPanel(
                                      #       condition = "input.edatope == 'B2'",
                                      #       
                                      #       radioButtons("spp.focal.1.B2", inline = TRUE,
                                      #                    label = "Select a species for ensemble detail",
                                      #                    choices = as.list(c("none", spps.all.B2)),
                                      #                    selected = "none")
                                      #     ),
                                      #     
                                      #     conditionalPanel(
                                      #       condition = "input.edatope == 'C4'",
                                      #       
                                      #       radioButtons("spp.focal.1.C4", inline = TRUE,
                                      #                    label = "Select a species for ensemble detail",
                                      #                    choices = as.list(c("none", spps.all.C4)),
                                      #                    selected = "none")
                                      #     ),
                                      #     
                                      #     conditionalPanel(
                                      #       condition = "input.edatope == 'D6'",
                                      #       
                                      #       radioButtons("spp.focal.1.D6", inline = TRUE,
                                      #                    label = "Select a species for ensemble detail",
                                      #                    choices = as.list(c("none", spps.all.D6)),
                                      #                    selected = "none")
                                      #     ),
                                      #     
                                      #   ),
                                      #   
                                      #   conditionalPanel(
                                      #     condition = "input.plotspp == 2",
                                      #     
                                      #     conditionalPanel(
                                      #       condition = "input.edatope == 'B2'",
                                      #       
                                      #       radioButtons("spp.focal.2.B2", inline = TRUE,
                                      #                    label = "Select a species for ensemble detail",
                                      #                    choices = as.list(c("none", spps.native.B2)),
                                      #                    selected = "none")
                                      #     ),
                                      #     
                                      #     conditionalPanel(
                                      #       condition = "input.edatope == 'C4'",
                                      #       
                                      #       radioButtons("spp.focal.2.C4", inline = TRUE,
                                      #                    label = "Select a species for ensemble detail",
                                      #                    choices = as.list(c("none", spps.native.C4)),
                                      #                    selected = "none")
                                      #     ),
                                      #     
                                      #     conditionalPanel(
                                      #       condition = "input.edatope == 'D6'",
                                      #       
                                      #       radioButtons("spp.focal.2.D6", inline = TRUE,
                                      #                    label = "Select a species for ensemble detail",
                                      #                    choices = as.list(c("none", spps.native.D6)),
                                      #                    selected = "none")
                                      #     ),
                                      #   ),
                                      #   
                                      # ),
                                      
                                      
                                      
                               ),
                               
                               column(12, 
                                      plotOutput(outputId = "scatterPlot", height="auto"),
                                      downloadButton(outputId = "downloadPlot", label = "Download plot")
                               ),
                        )
                      ),
                      
                      
                      column(width = 12,
                             style = "background-color:#003366; border-top:2px solid #fcba19;position: absolute; bottom:0%; ",
                             
                             tags$footer(class="footer", 
                                         tags$div(class="container", style="display:flex; justify-content:center; flex-direction:column; text-align:center; height:46px;",
                                                  tags$ul(style="display:flex; flex-direction:row; flex-wrap:wrap; margin:0; list-style:none; align-items:center; height:100%;",
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home", "Home", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/disclaimer", "Disclaimer", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/privacy", "Privacy", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/accessibility", "Accessibility", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/copyright", "Copyright", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/StaticWebResources/static/gov3/html/contact-us.html", "Contact", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;"))
                                                  )
                                         )
                             )
                      )
             ),
             
             tabPanel("About",
                      
                      includeMarkdown("about.Rmd"),
                      
                      column(width = 12,
                             style = "background-color:#003366; border-top:2px solid #fcba19;",
                             
                             tags$footer(class="footer",
                                         tags$div(class="container", style="display:flex; justify-content:center; flex-direction:column; text-align:center; height:46px;",
                                                  tags$ul(style="display:flex; flex-direction:row; flex-wrap:wrap; margin:0; list-style:none; align-items:center; height:100%;",
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home", "Home", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/disclaimer", "Disclaimer", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/privacy", "Privacy", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/accessibility", "Accessibility", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/copyright", "Copyright", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/StaticWebResources/static/gov3/html/contact-us.html", "Contact", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;"))
                                                  )
                                         )
                             )
                      )
             ),
             
             tabPanel("Find-a-BEC",
                      sidebarLayout(
                        sidebarPanel(
                          helpText("Choose a BGC zone or subzone-variant to show it on the map"),
                          
                          tags$head(tags$script('$(document).on("shiny:connected", function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            });
                            $(window).resize(function(e) {
                            Shiny.onInputChange("innerWidth", window.innerWidth);
                            });
                            ')),
                          
                          selectInput("showbgc", 
                                      label = "Choose a BGC subzone-variant",
                                      choices = as.list(c("none", bgc.list)),
                                      selected = "none"),
                          
                          selectInput("showzone", 
                                      label = "Choose a BGC zone",
                                      choices = as.list(c("none", zone.list)),
                                      selected = "none"),
                          
                        ),    
                        
                        mainPanel(
                          
                          leafletOutput(outputId = "becmap", height="86vh")
                          
                        )
                      ),
                      column(width = 12,
                             style = "background-color:#003366; border-top:2px solid #fcba19;",
                             
                             tags$footer(class="footer",
                                         tags$div(class="container", style="display:flex; justify-content:center; flex-direction:column; text-align:center; height:46px;",
                                                  tags$ul(style="display:flex; flex-direction:row; flex-wrap:wrap; margin:0; list-style:none; align-items:center; height:100%;",
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home", "Home", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/disclaimer", "Disclaimer", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/privacy", "Privacy", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/accessibility", "Accessibility", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/copyright", "Copyright", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/StaticWebResources/static/gov3/html/contact-us.html", "Contact", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;"))
                                                  )
                                         )
                             )
                      )
             ),
             
             tabPanel("Model Info",
                      DT::dataTableOutput("table"),
                      column(width = 12,
                             style = "background-color:#003366; border-top:2px solid #fcba19;",
                             
                             tags$footer(class="footer",
                                         tags$div(class="container", style="display:flex; justify-content:center; flex-direction:column; text-align:center; height:46px;",
                                                  tags$ul(style="display:flex; flex-direction:row; flex-wrap:wrap; margin:0; list-style:none; align-items:center; height:100%;",
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home", "Home", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/disclaimer", "Disclaimer", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/privacy", "Privacy", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/accessibility", "Accessibility", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/gov/content/home/copyright", "Copyright", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;")),
                                                          tags$li(a(href="https://www2.gov.bc.ca/StaticWebResources/static/gov3/html/contact-us.html", "Contact", style="font-size:1em; font-weight:normal; color:white; padding-left:5px; padding-right:5px; border-right:1px solid #4b5e7e;"))
                                                  )
                                         )
                             )
                      )
             )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  
  
  
  output$map <- renderLeaflet({
    
    
    leaflet() %>% 
      addTiles() %>% 
      # addProviderTiles("Esri.WorldImagery", group = "Satellite view") %>%
      # addProviderTiles("Esri.WorldTerrain", group = "Terrain only") %>%
      fitBounds(lng1 = bbox[1], lat1 = bbox[3], lng2 = bbox[2], lat2 = bbox[4]) %>%
      # addLayersControl(
      #   baseGroups = c("Base map", "Terrain only", "Satellite view"),
      #   options = layersControlOptions(collapsed = FALSE),
      # ) %>%
      addPolygons(data=bdy, fillColor = NA, color="black", smoothFactor = 0.2, fillOpacity = 0, weight=2)
    
    
    
  }
  )
  
  observe({
    
    # proj.year <- 2041
    # zonelevel=F
    # transparency <- 0.7
    # gcm.focal <- gcms[2]
    
    zonelevel <- if(input$zonelevel==T) T else F
    gcm.focal <- input$gcm.focal
    scenario <- "ssp245"
    proj.year <-  proj.years[as.numeric(input$proj.year)+1]
    transparency <- input$transparency
    
    if(input$maptype==1) X <- bgc.pred.ref
    if(input$maptype==2) X <- bgc.pred.2001
    if(input$maptype==3) X <- get(paste("bgc.pred", gcm.focal, scenario, proj.year, sep="."))
    BGC.pred <- levels.bgc[values(X)]
    
    zone.pred <- rep(NA, length(BGC.pred))
    for(i in zones.all){ zone.pred[grep(i,BGC.pred)] <- i }
    
    ColScheme <- if(zonelevel==T) zonecolors$colour[match(zones.all, zonecolors$classification)] else bgccolors$colour[match(bgcs.all, bgccolors$classification)]
    units <- if(zonelevel==T) zones.all else bgcs.all
    pred <- if(zonelevel==T) zone.pred else BGC.pred
    
    values(X) <- factor(pred, levels=units)
    values(X)[1:length(units)] <- 1:length(units) # this is a patch that is necessary to get the color scheme right.
    
    if(input$type==2){
      
      leafletProxy("map") %>%
        addProviderTiles("Esri.WorldTopoMap", group = "Base map") %>%
        addRasterImage(X, colors = ColScheme, method="ngb", opacity = transparency, maxBytes = 6 * 1024 * 1024)%>%
        addPolygons(data=bdy, fillColor = NA, color="black", smoothFactor = 0.2, fillOpacity = 0, weight=2)
      
    } 
    if(input$type==3){
      
      X <- bgc.pred.ref
      values(X) <- NA
      
      edatope <- input$edatope
      if(input$plotspp==1){
        if(edatope=="B2") spp.focal <- input$spp.focal.1.B2
        if(edatope=="C4") spp.focal <- input$spp.focal.1.C4
        if(edatope=="D6") spp.focal <- input$spp.focal.1.D6
      } else if(input$plotspp==2){
        if(edatope=="B2") spp.focal <- input$spp.focal.2.B2
        if(edatope=="C4") spp.focal <- input$spp.focal.2.C4
        if(edatope=="D6") spp.focal <- input$spp.focal.2.D6
      }  
      
      # spp.focal <- get(paste("input$spp.focal", input$plotspp, edatope, sep="."))
      
      
      if(input$mapspp==1){
        values(X) <- NA
        if(spp.focal!="none") {
          suit <- SuitLookup$ESuit[which(SuitLookup$Spp==spp.focal)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), SuitLookup$SS_NoSpace[which(SuitLookup$Spp==spp.focal)])]
          if(input$maptype==3){
            if(gcm.focal=="ensemble"){
              suit.ref <- suit[match(levels.bgc[values(bgc.pred.ref)], SiteLookup$BGC)]
              suit.ref[is.na(suit.ref)] <- 4 #set non-suitable to 4
              change.proj <- values(raster(paste("data/Spp.ChangeSuit", studyarea, spp.focal, edatope, scenario, proj.year, "tif", sep=".")))
              temp <- suit.ref-change.proj
              temp[which(temp>3.5)] <- NA #set non-suitable to NA
              temp[which(temp<1)] <- 1 #set non-suitable to NA
            } else {
              temp <- suit[match(BGC.pred, SiteLookup$BGC)]
              temp[which(temp>3.5)] <- NA #set non-suitable to NA
            }
          } else {
            temp <- suit[match(BGC.pred, SiteLookup$BGC)]
            temp[which(temp>3.5)] <- NA #set non-suitable to NA
          }
          values(X) <- temp
        }
        leafletProxy("map") %>%
          addProviderTiles("Esri.WorldTopoMap", group = "Base map") %>%
          addRasterImage(X, colors =  ColScheme.suit, method="ngb", opacity = transparency, maxBytes = 6 * 1024 * 1024)%>%
          addPolygons(data=bdy, fillColor = NA, color="black", smoothFactor = 0.2, fillOpacity = 0, weight=2)
      } 
      if(input$mapspp==2){
        values(X) <- NA
        if(input$maptype==3){
          if(spp.focal!="none") X <- raster(paste("data/Spp.ChangeSuit", studyarea, spp.focal, edatope, scenario, proj.year, "tif", sep="."))
          leafletProxy("map") %>%
            addProviderTiles("Esri.WorldTopoMap", group = "Base map") %>%
            addRasterImage(X, colors =  ColScheme.change, method="ngb", opacity = transparency, maxBytes = 6 * 1024 * 1024)%>%
            addPolygons(data=bdy, fillColor = NA, color="black", smoothFactor = 0.2, fillOpacity = 0, weight=2)
        }
      } 
      if(input$mapspp==3){
        values(X) <- NA
        if(input$maptype==3){
          if(spp.focal!="none") X <- raster(paste("data/Spp.binary", studyarea, spp.focal, edatope, scenario, proj.year, "tif", sep="."))
          leafletProxy("map") %>%
            addProviderTiles("Esri.WorldTopoMap", group = "Base map") %>%
            addRasterImage(X, colors =  ColScheme.binary, method="ngb", opacity = transparency, maxBytes = 6 * 1024 * 1024)%>%
            addPolygons(data=bdy, fillColor = NA, color="black", smoothFactor = 0.2, fillOpacity = 0, weight=2)
        } # end if(input$maptype==3)
      }
    }
    
  })
  
  # Use a separate observer to recreate the legend as needed.
  observe({
    proxy <- leafletProxy("map")
    
    # Remove any existing legend, and only if the legend is
    # enabled, create a new one.
    proxy %>% clearControls()
    if (input$type==3) {
      if(input$mapspp==1) proxy %>% addLegend(colors =  palette.suit, labels=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), title = "Climatic feasibility")
      if(input$mapspp==2) proxy %>% addLegend(colors = palette.change, labels = labels.change, title = "Ensemble-mean change</br>in climatic feasibility")
      if(input$mapspp==3) proxy %>% addLegend(colors = palette.binary, labels = labels.binary, title = "% of ensemble predicting</br>loss (red) or gain (blue)</br>of climatic feasibility")
    }
  })
  
  #Show popup on click
  observeEvent(input$map_click, {
    gcm.focal <- input$gcm.focal
    scenario <- "ssp245"
    proj.year <-  proj.years[as.numeric(input$proj.year)+1]
    
    if(input$maptype==1) X <- bgc.pred.ref
    if(input$maptype==2) X <- bgc.pred.2001
    if(input$maptype==3) X <- get(paste("bgc.pred", gcm.focal, scenario, proj.year, sep="."))
    BGC.pred <- levels.bgc[values(X)]
    
    click <- input$map_click
    bgc.popup <- BGC.pred[cellFromXY(X, matrix(c(click$lng, click$lat), 1))]
    text<-paste0("<strong>", bgc.popup, "</strong>", "<br/>Zone: ", bgc.names$ZoneName[which(bgc.names$Map_Label==bgc.popup)], "<br/>Subzone/Variant: ",  bgc.names$SubzoneName[which(bgc.names$Map_Label==bgc.popup)])
    proxy <- leafletProxy("map")
    proxy %>% clearPopups() %>%
      addPopups(click$lng, click$lat, text)
  })
  
  # Map download
  output$downloadMap <- downloadHandler(
    filename =  function(){
      paste(studyarea, 
            "map",
            if(input$zonelevel==T) "zones" else "bgcs", 
            if(input$maptype==1) "ref" else if(input$maptype==2) "obs" else paste("proj",input$gcm.focal,proj.years[as.numeric(input$proj.year)+1],sep="."),
            "png", sep=".")
    },
    
    content = function(file) {
      
      zonelevel <- if(input$zonelevel==T) T else F
      gcm.focal <- input$gcm.focal
      scenario <- "ssp245"
      proj.year <-  proj.years[as.numeric(input$proj.year)+1]
      transparency <- input$transparency
      
      if(input$maptype==1) X <- bgc.pred.ref
      if(input$maptype==2) X <- bgc.pred.2001
      if(input$maptype==3) X <- get(paste("bgc.pred", gcm.focal, scenario, proj.year, sep="."))
      BGC.pred <- levels.bgc[values(X)]
      
      zone.pred <- rep(NA, length(BGC.pred))
      for(i in zones.all){ zone.pred[grep(i,BGC.pred)] <- i }
      
      ColScheme <- if(zonelevel==T) zonecolors$colour[match(zones.all, zonecolors$classification)] else bgccolors$colour[match(bgcs.all, bgccolors$classification)]
      units <- if(zonelevel==T) zones.all else bgcs.all
      pred <- if(zonelevel==T) zone.pred else BGC.pred
      
      values(X) <- factor(pred, levels=units)
      values(X)[1:length(units)] <- 1:length(units) # this is a patch that is necessary to get the color scheme right.
      
      png(file, width = 7*300, height = 7.5*300, res = 300)
      par(mar=c(0,0,0,0))
      
      plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
      values(X)[-(1:length(units))] <- NA # cover up the color bar
      image(X, add=T, col="white") # cover up the color bar
      plot(bdy, add=T, lwd=1)  
      
      temp <- table(pred)
      temp <- rev(sort(temp))
      legendunits <- names(temp)
      legendunits <- legendunits[which(temp > (totalarea*0.005))]
      
      if(zonelevel==T) legend("topright", legend=legendunits, fill=zonecolors$colour[match(legendunits, zonecolors$classification)], ncol= if(length(legendunits)>10) 2 else 1, bty="n", cex=1.4)
      
      # box()
      dev.off()
    }
  )
  
  
  
  
  
  ## Plot window (done as a function so that the user can export)
  scatterPlot <- function() {
    
    # proj.year <- 2041
    # scenario <- "ssp245"
    # var1 <- "MAT"
    # var2 <- "MAP"
    # ratioscale <- T
    # zonelevel=T
    # unit.area.focal <- if(zonelevel==T) zones[1] else bgcs[1]
    # unit.persistence.focal <- if(zonelevel==T) zones.native[1] else bgcs.native[1]
    # units.area <- if(zonelevel==T) zones else bgcs
    # fractional=F
    # edatope <- edatopes[1]
    # spp.focal <- "Fd"
    # recent <- T
    # gcm.focal <- "ensemble"
    
    proj.year <-  proj.years[as.numeric(input$proj.year)+1]
    # scenario <- scenarios[as.numeric(input$scenario)]
    var1 <- input$var1
    var2 <- input$var2
    zonelevel <- if(input$zonelevel==T) T else F
    unit.area.focal <- if(zonelevel==T) input$zone.area.focal else input$bgc.area.focal
    unit.persistence.focal <- if(zonelevel==T) input$zone.persistence.focal else input$bgc.persistence.focal
    units.area <- if(zonelevel==T) zones else bgcs
    # fractional <- if(input$fractional==T) T else F
    # edatope <- input$edatope
    recent <- input$recent
    gcm.focal <- input$gcm.focal
    
    variable.type1 <- variable.types[which(variables==var1)]
    variable.type2 <- variable.types[which(variables==var2)]
    
    if(input$type==1) {
      
      #-------------------------
      # climate scatterplot
      #-------------------------
      
      ratioscale <- if(input$ratioscale==T) T else F
      
      data <- if(ratioscale==T & variable.type1=="ratio") clim.meanChange.ratio else if(ratioscale==T & variable.type2=="ratio") clim.meanChange.ratio else clim.meanChange
      
      x <- data[, which(variables==var1)]
      y <- data[, which(variables==var2)]
      
      xlim=range(x)*c(if(min(x)<0) 1.1 else 0.9, if(max(x)>0) 1.1 else 0.9)
      ylim=range(y)*c(if(min(y)<0) 1.1 else 0.9, if(max(y)>0) 1.1 else 0.9)
      
      par(mar=c(3,4,0,1), mgp=c(1.25, 0.25,0), cex=1.5)
      plot(x,y,col="white", tck=0, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, ylab="",
           xlab=paste("Change in", if(var1%in%variable.names$Code) variable.names$Variable[which(variable.names$Code==var1)] else var1))
      par(mgp=c(2.5,0.25, 0))
      title(ylab=paste("Change in", if(var2%in%variable.names$Code) variable.names$Variable[which(variable.names$Code==var2)] else var2))
      lines(c(0,0), c(-99,99), lty=2, col="gray")
      lines(c(-99,99), c(0,0), lty=2, col="gray")
      
      if(recent==T){
        x1 <- data[2, which(variables==var1)]
        y1 <- data[2, which(variables==var2)]
        points(x1,y1, pch=16, col="gray", cex=2.5)
        text(x1,y1, "2001-2020 (observed)", cex=1, font=2, pos=4, col="gray", offset=0.9)  
      }
      
      for(gcm in gcms){
        i=which(gcms==gcm)
        x2 <- data[c(1, which(period[,1]==gcm)), which(variables==var1)]
        y2 <- data[c(1, which(period[,1]==gcm)), which(variables==var2)]
        if(length(unique(sign(diff(x2))))==1){
          x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
          y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
          s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
          lines(s, col=ColScheme.gcms[i], lwd=2)
        } else lines(x2, y2, col=ColScheme.gcms[i], lwd=2)
        j=which(proj.years==proj.year)
        points(x2,y2, pch=21, bg=ColScheme.gcms[i], cex=1)
        points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm==gcm.focal) 3.5 else 3)
        text(x2[j],y2[j], mods[i], cex=if(gcm==gcm.focal) 0.7 else 0.5, font=2)
      }
      
      axis(1, at=pretty(x), labels=if(ratioscale==T & variable.type1=="ratio") paste(pretty(x)*100, "%", sep="") else pretty(x), tck=0)
      axis(2, at=pretty(y), labels=if(ratioscale==T & variable.type2=="ratio") paste(pretty(y)*100, "%", sep="") else pretty(y), las=2, tck=0)
      
    } else if(input$type==2){ 
      
      if(input$plotbgc==1){
        #-------------------------
        # BGC scatterplot
        #-------------------------
        
        
        data <- if(zonelevel==T) zone.area else bgc.area
        clim.data <- clim.meanChange
        ColScheme <- if(zonelevel==T) zonecolors else bgccolors
        units <- if(zonelevel==T) zones else bgcs
        
        x <- clim.data[, which(variables==var1)]
        variable.type1 <- variable.types[which(variables==var1)]
        
        xlim=range(x)*c(if(min(x)<0) 1.1 else 0.9, if(max(x)>0) 1.1 else 0.9)-c(diff(range(x))/4, 0)
        ylim=c(0, max(data, na.rm=T)*1.05)
        
        par(mar=c(3,4,0,1), mgp=c(1.25, 0.25,0), cex=1.5)
        plot(0,col="white", tck=0, xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim, ylab="",
             xlab=paste("Change in", variable.names$Variable[which(variable.names$Code==var1)]), 
        )
        par(mgp=c(2.5,0.25, 0))
        title(ylab=paste("Area of biogeoclimatic unit ('000 sq.km)"))
        
        axis(1, at=pretty(x), labels=pretty(x), tck=0)
        axis(2, at=pretty(ylim), labels=pretty(ylim)/1000, tck=0, las=2)
        
        lines(rep(x[which(period[,1]==gcm.focal & period[,3]==proj.year)], 2), ylim, col="gray90", lwd=2)
        
        increasing <- data[which(period$GCM==gcm.focal & period$proj.year==2081),]>data[1,]
        order.increasing <- rev(order(data[which(period$GCM==gcm.focal & period$proj.year==2081),increasing==T]))
        order.decreasing <- rev(order(data[1,increasing==F]))
        units.increasing <- units.area[increasing==T][order.increasing]
        units.decreasing <- units.area[increasing==F][order.decreasing]
        data.increasing <-  if(length(which(increasing==T))>1) data[,increasing==T][,order.increasing] else data[,increasing==T]
        data.decreasing <- if(length(which(increasing==F))>1) data[,increasing==F][,order.decreasing] else data[,increasing==F]
        units.sort <- c(units.increasing, units.decreasing)
        data.sort <- cbind(data.increasing, data.decreasing)
        increasing.sort <- increasing[match(units.sort, units)]
        for(unit in units.sort){
          i <- which(units.sort==unit)
          col.focal <- if(unit.area.focal=="none") ColScheme$colour[which(ColScheme$classification==unit)] else "lightgray"
          col.focal2 <- if(unit.area.focal=="none") "black" else "darkgray"
          x1 <- x[c(1, which(period[,1]==gcm.focal))]
          y1 <- data.sort[c(1, which(period[,1]==gcm.focal)), which(units.sort==unit)]
          y1[is.na(y1)] <- 0
          if(length(unique(sign(diff(x1))))==1){
            x3 <- if(unique(sign(diff(x1)))==-1) rev(x1) else x1
            y3 <- if(unique(sign(diff(x1)))==-1) rev(y1) else y1
            s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
            lines(s, col=col.focal, lwd=2)
          } else lines(x1, y1, col=col.focal, lwd=3)
          
          # labels
          position <- rep(0:2, times=100)
          side <- if(increasing.sort[i]==T) 4 else 2
          space <- 12
          lines(if(increasing.sort[i]==T) c(x1[length(x1)], x1[length(x1)]+position[i]*diff(range(x))/space) else c(x1[1], x1[1]-position[i]*diff(range(x))/space),
                if(increasing.sort[i]==T) rep(y1[length(y1)],2) else rep(y1[1],2), col=col.focal, lty=2)
          text(if(increasing.sort[i]==T) x1[length(x1)]+position[i]*diff(range(x))/space else x1[1]-position[i]*diff(range(x))/space,
               if(increasing.sort[i]==T) y1[length(y1)] else y1[1],
               unit, col=col.focal, pos=side, font=2, cex=0.7, offset=0.1)
          
          if(recent==T){
            x1 <- x[1:2]
            y1 <- data.sort[1:2, which(units.sort==unit)]
            lines(x1, y1, col=col.focal, lwd=1.25, lty=1)
            points(x1[2],y1[2], pch=21, bg=col.focal, col=col.focal2, cex=1.2)
          }
          
        }
        
        if(unit.area.focal!="none"){
          for(gcm in gcms){
            i=which(gcms==gcm)
            x2 <- x[c(1, which(period[,1]==gcm))]
            y2 <- data[c(1, which(period[,1]==gcm)), which(units==unit.area.focal)]
            if(length(unique(sign(diff(x2))))==1){
              x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
              y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
              s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
              lines(s, col=ColScheme.gcms[i], lwd=if(gcm==gcm.focal) 4 else 2)
            } else lines(x2, y2, col=ColScheme.gcms[i], lwd=if(gcm==gcm.focal) 4 else 2)
            j=which(proj.years==proj.year)
            points(x2,y2, pch=21, bg=ColScheme.gcms[i], cex=1)
            points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm==gcm.focal) 3.5 else 3)
            text(x2[j],y2[j], mods[i], cex=if(gcm==gcm.focal) 0.7 else 0.5, font=2)
          }

          if(recent==T){
            x1 <- x[1:2]
            y1 <- data.sort[1:2, which(units.sort==unit.area.focal)]
            lines(x1, y1, col="gray", lwd=1.25, lty=1)
            points(x1[2],y1[2], pch=21, bg="gray", col=1, cex=2)
            text(x1[2],y1[2], "2001-2020 (observed)", cex=1, font=2, pos=4, col="gray30", offset=0.7)  
          }
          
        }
        box()
        
      } else if(input$plotbgc==2){
        
        #-------------------------
        # bgc bubbleplot
        #-------------------------
        
        persistence <- if(zonelevel==T) zone.persistence else bgc.persistence
        expansion <- if(zonelevel==T) zone.expansion else bgc.expansion
        ColScheme <- if(zonelevel==T) zonecolors else bgccolors
        # units <- names(persistence)
        units.native <- if(zonelevel==T) zones.native else bgcs.native
        units <- units.native
        
        par(mar=c(3,4,0.1,0.1), mgp=c(1.25, 0.25, 0), cex=1.5) 
        
        xlim <- c(0, 1.1)
        ylim <- c(-5,3)
        plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab="Persistence within historically feasible range", ylab="")
        axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
        axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
        par(mgp=c(2.75, 0.25, 0))
        title(ylab="Expansion beyond historically feasible range", cex.lab=1)
        iso <- seq(0,1.2, 0.001)
        lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
        # arctext(x = "Growing feasible range", center = c(-1, -28.7), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
        # arctext(x = "Shrinking feasible range", center = c(-1, -29.3), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
        # mtext(paste(edatope.name[which(edatopes==edatope)], " sites", " (", edatope, ")", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)
        
        unit=units[2]
        for(unit in units){
          col.focal <- if(unit.persistence.focal=="none") ColScheme$colour[which(ColScheme$classification==unit)] else "lightgray"
          col.focal2 <- if(unit.persistence.focal=="none") "black" else "darkgray"
          x <- persistence[which(period[,3]==proj.year), which(names(persistence)==unit)]
          y <- expansion[which(period[,3]==proj.year), which(names(persistence)==unit)]
          y[y<2^(ylim[1])] <- 2^(ylim[1])
          y <- log2(y)
          
          # points(x,y)
          if(length(x)>1){
            if(var(y)==0) lines(range(x, na.rm=T), range(y), col=col.focal) else dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
          }
          points(mean(x),mean(y), pch=21, bg=col.focal, cex=if(unit==unit.persistence.focal) 4.5 else 3, col=col.focal2)
          text(mean(x),mean(y), unit, cex=if(unit==unit.persistence.focal) 1 else 0.7, font=2, col=col.focal2)
        }
        
        if(unit.persistence.focal!="none"){
          
          for(gcm in gcms){
            i=which(gcms==gcm)
            x2 <- persistence[c(1, which(period[,1]==gcm)), which(names(persistence)==unit.persistence.focal)]
            y2 <- expansion[c(1, which(period[,1]==gcm)), which(names(expansion)==unit.persistence.focal)]
            y2[y2<2^(ylim[1])] <- 2^(ylim[1])
            y2 <- log2(y2)
            # if(all(diff(x2) > 0)){
            if(length(unique(sign(diff(x2))))==1 & sum(diff(x2))!=0){
              x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
              y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
              s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
              lines(s, col=ColScheme.gcms[i], lwd=2)
            } else lines(x2, y2, col=ColScheme.gcms[i], lwd=2)
            j=which(proj.years==proj.year)
            points(x2,y2, pch=21, bg=ColScheme.gcms[i], cex=1)
            points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm==gcm.focal) 3.5 else 3)
            text(x2[j],y2[j], mods[i], cex=if(gcm==gcm.focal) 0.7 else 0.5, font=2)
          }          
          
          if(recent==T){
            x1 <- persistence[2, which(names(persistence)==unit.persistence.focal)]
            y1 <- expansion[2, which(names(expansion)==unit.persistence.focal)]
            y1[y1<2^(ylim[1])] <- 2^(ylim[1])
            y1 <- log2(y1)
            points(x1,y1, pch=21, bg="gray", col=1, cex=2)
            text(x1,y1, "2001-2020 (observed)", cex=1, font=2, pos=4, col="gray30", offset=0.7)  
          }
          
        }
        box()
        
      }
    } else if(input$type==3){ 
      
      if(input$plotspp==1){
        
        if(edatope=="B2") spp.focal <- input$spp.focal.1.B2
        if(edatope=="C4") spp.focal <- input$spp.focal.1.C4
        if(edatope=="D6") spp.focal <- input$spp.focal.1.D6
        
        #-------------------------
        # species scatterplot
        #-------------------------
        
        data <- if(fractional==T) get(paste("suit.area", edatope, sep=".")) else get(paste("spp.area", edatope, sep="."))
        clim.data <- clim.meanChange
        ColScheme <- sppcolors
        spps <- names(data)
        
        x <- clim.data[, which(variables==var1)]
        variable.type1 <- variable.types[which(variables==var1)]
        
        xlim=range(x)*c(if(min(x)<0) 1.1 else 0.9, if(max(x)>0) 1.1 else 0.9)-c(diff(range(x))/4, 0)
        ylim=c(0, max(data, na.rm=T)*1.05)
        
        par(mar=c(3,4,0,1), mgp=c(1.25, 0.25,0), cex=1.5)
        plot(0,col="white", tck=0, xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim, ylab="",
             xlab=paste("Change in", variable.names$Variable[which(variable.names$Code==var1)]), 
        )
        par(mgp=c(2.5,0.25, 0))
        title(ylab=paste("Tree species' feasible area ('000 sq.km)"))
        
        axis(1, at=pretty(x), labels=pretty(x), tck=0)
        axis(2, at=pretty(ylim), labels=pretty(ylim)/1000, tck=0, las=2)
        
        lines(rep(x[which(period[,1]==gcm.focal & period[,3]==proj.year)], 2), ylim, col="gray90", lwd=2)
        
        increasing <- data[which(period$GCM==gcm.focal & period$proj.year==2081),]>data[1,]
        order.increasing <- rev(order(data[which(period$GCM==gcm.focal & period$proj.year==2081),increasing==T]))
        order.decreasing <- rev(order(data[1,increasing==F]))
        spps.increasing <- spps[increasing==T][order.increasing]
        spps.decreasing <- spps[increasing==F][order.decreasing]
        data.increasing <- data[,increasing==T][,order.increasing]
        data.decreasing <- data[,increasing==F][,order.decreasing]
        spps.sort <- c(spps.increasing, spps.decreasing)
        data.sort <- cbind(data.increasing, data.decreasing)
        increasing.sort <- increasing[match(spps.sort, spps)]
        for(spp in spps.sort){
          i <- which(spps.sort==spp)
          col.focal <- if(spp.focal=="none") sppcolors[i] else "lightgray"  
          col.focal2 <- if(spp.focal=="none") "black" else "darkgray"  
          x1 <- x[c(1, which(period[,1]==gcm.focal))]
          y1 <- data.sort[c(1, which(period[,1]==gcm.focal)), which(spps.sort==spp)]
          y1[is.na(y1)] <- 0
          if(length(unique(sign(diff(x1))))==1){
            x3 <- if(unique(sign(diff(x1)))==-1) rev(x1) else x1
            y3 <- if(unique(sign(diff(x1)))==-1) rev(y1) else y1
            s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
            lines(s, col=col.focal, lwd=3)
          } else lines(x1, y1, col=col.focal, lwd=3)
          
          # labels
          position <- rep(0:2, times=100)
          side <- if(increasing.sort[i]==T) 4 else 2
          space <- 12
          lines(if(increasing.sort[i]==T) c(x1[length(x1)], x1[length(x1)]+position[i]*diff(range(x))/space) else c(x1[1], x1[1]-position[i]*diff(range(x))/space), 
                if(increasing.sort[i]==T) rep(y1[length(y1)],2) else rep(y1[1],2), col=col.focal, lty=2)
          text(if(increasing.sort[i]==T) x1[length(x1)]+position[i]*diff(range(x))/space else x1[1]-position[i]*diff(range(x))/space, 
               if(increasing.sort[i]==T) y1[length(y1)] else y1[1], 
               spp, col=col.focal, pos=side, font=2, cex=0.7, offset=0.1)
          
          if(recent==T){
            x1 <- x[1:2]
            y1 <- data.sort[1:2, which(spps.sort==spp)]
            lines(x1, y1, col=col.focal, lwd=1.25, lty=1)  
            points(x1[2],y1[2], pch=21, bg=col.focal, col=col.focal2, cex=1.2)
          }
          
        }
        
        if(spp.focal!="none"){
          for(gcm in gcms){
            i=which(gcms==gcm)
            x2 <- x[c(1, which(period[,1]==gcm))]
            y2 <- data[c(1, which(period[,1]==gcm)), which(spps==spp.focal)]
            if(length(unique(sign(diff(x2))))==1){
              x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
              y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
              s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
              lines(s, col=ColScheme.gcms[i], lwd=if(gcm==gcm.focal) 4 else 2)
            } else lines(x2, y2, col=ColScheme.gcms[i], lwd=if(gcm==gcm.focal) 4 else 2)
            j=which(proj.years==proj.year)-1
            points(x2,y2, pch=21, bg=ColScheme.gcms[i], cex=1)
            points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm==gcm.focal) 3.5 else 3)
            text(x2[j],y2[j], mods[i], cex=if(gcm==gcm.focal) 0.7 else 0.5, font=2)
          }
        }
        
        
      } else if(input$plotspp==2){
        
        #-------------------------
        # species bubbleplot
        #-------------------------
        
        if(edatope=="B2") spp.focal <- input$spp.focal.2.B2
        if(edatope=="C4") spp.focal <- input$spp.focal.2.C4
        if(edatope=="D6") spp.focal <- input$spp.focal.2.D6
        
        persistence <- if(fractional==T) get(paste("suit.persistence", edatope, sep=".")) else get(paste("spp.persistence", edatope, sep=".")) 
        expansion <- if(fractional==T) get(paste("suit.expansion", edatope, sep=".")) else get(paste("spp.expansion", edatope, sep=".")) 
        spps <- names(persistence)
        
        par(mar=c(3,4,0.1,0.1), mgp=c(1.25, 0.25, 0), cex=1.5) 
        
        xlim <- c(0, 1.5)
        ylim <- c(-5,3)
        plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab="Persistence within historically feasible range", ylab="")
        axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
        axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
        par(mgp=c(2.75, 0.25, 0))
        title(ylab="Expansion beyond historically feasible range", cex.lab=1)
        iso <- seq(0,1.2, 0.001)
        lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
        # arctext(x = "Growing feasible range", center = c(-1, -28.7), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
        # arctext(x = "Shrinking feasible range", center = c(-1, -29.3), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
        # mtext(paste(edatope.names[which(edatopes==edatope)], " sites", " (", edatope, ")", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)
        
        spp=spps[2]
        for(spp in spps){
          i <- which(spps==spp)
          col.focal <- if(spp.focal=="none") sppcolors[i] else "lightgray"  
          col.focal2 <- if(spp.focal=="none") "black" else "darkgray"  
          x <- persistence[which(period[,3]==proj.year), i]
          y <- expansion[which(period[,3]==proj.year), i]
          y[y<2^(ylim[1])] <- 2^(ylim[1])
          y <- log2(y)
          
          # points(x,y)
          if(length(x)>1){
            if(var(y)==0) lines(range(x), range(y), col=col.focal) else dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
          }
          points(mean(x),mean(y), pch=21, bg=col.focal, cex=if(spp==spp.focal) 4.5 else 3, col=col.focal2)
          text(mean(x),mean(y), spp, cex=if(spp==spp.focal) 1 else 0.7, font=2, col=col.focal2)
        }
        
        if(spp.focal!="none"){
          for(gcm in gcms){
            i=which(gcms==gcm)
            x2 <- persistence[c(1, which(period[,1]==gcm)), which(spps==spp.focal)]
            y2 <- expansion[c(1, which(period[,1]==gcm)), which(spps==spp.focal)]
            y2[y2<2^(ylim[1])] <- 2^(ylim[1])
            y2 <- log2(y2)
            # if(all(diff(x2) > 0)){
            if(length(unique(sign(diff(x2))))==1 & sum(diff(x2))!=0){
              x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
              y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
              s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
              lines(s, col=ColScheme.gcms[i], lwd=2)
            } else lines(x2, y2, col=ColScheme.gcms[i], lwd=2)
            j=which(proj.years==proj.year)-1
            points(x2,y2, pch=21, bg=ColScheme.gcms[i], cex=1)
            points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm==gcm.focal) 3.5 else 3)
            text(x2[j],y2[j], mods[i], cex=if(gcm==gcm.focal) 0.7 else 0.5, font=2)
          }
        }
        
        
      }
    }
  }
  output$scatterPlot <- renderPlot({ scatterPlot() },
                                   height=reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*0.25,0))
  )
  
  # Plot download
  output$downloadPlot <- downloadHandler(
    filename =  function(){
      paste(studyarea, 
            if(input$type==1){
              paste(
                "climate",
                input$var1, 
                input$var2, 
                sep="."
              )
            } else if(input$type==2){
              paste(
                "bgc", 
                if(input$plotbgc==1){
                  paste("area", if(input$zonelevel==T) input$zone.area.focal else input$bgc.area.focal, sep=".")
                } else {
                  paste("persistence", if(input$zonelevel==T) input$zone.persistence.focal else input$bgc.persistence.focal, sep=".")
                },
                sep="."
              )
            } else if(input$type==3){
              paste(
                "spp",
                sep="."
              )
            },
            proj.years[as.numeric(input$proj.year)+1], 
            "png", sep=".")
    },
    
    content = function(file) {
      
      pixelratio <- session$clientData$pixelratio
      width  <- session$clientData$output_scatterPlot_width
      height <- session$clientData$output_scatterPlot_height
      
      png(file, width = width*pixelratio*1.5, height = height*pixelratio*2, res = 120*pixelratio)
      scatterPlot()
      dev.off()
    }
  )
  
  #-------------------------
  # Find-a-BEC
  #-------------------------
  
  # showbgc <- "BGxh1"
  # showzone <- "BG"
  
  output$becmap <- renderLeaflet({
    
    showbgc <- input$showbgc
    showzone <- input$showzone
    
    leaflet() %>% 
      addTiles() %>% 
      addProviderTiles("Esri.WorldImagery", group = "Satellite view") %>%
      addProviderTiles("Esri.WorldTerrain", group = "Terrain only") %>%
      addProviderTiles("Esri.WorldTopoMap", group = "Base map") %>%
      fitBounds(lng1 = extent(bgc.simple)[1], lat1 = extent(bgc.simple)[3], lng2 = extent(bgc.simple)[2], lat2 = extent(bgc.simple)[4]) %>%
      addLayersControl(
        baseGroups = c("Base map", "Terrain only", "Satellite view"),
        options = layersControlOptions(collapsed = FALSE),
      ) %>%
      addPolygons(data=bgc.simple[zone.maprecord == showzone,], fillColor = "red", color="red", smoothFactor = 0.2, fillOpacity = 0.4, weight=2, opacity=1)%>%
      addPolygons(data=bgc.simple[bgc.maprecord == showbgc,], fillColor = "black", color="black", smoothFactor = 0.2, fillOpacity = 0.4, weight=2, opacity=1) 
    
  },
  )
  
  #-------------------------
  # Model Metadata Table
  #-------------------------
  
  output$table <- DT::renderDataTable({
    DT::datatable(modelMetadata, 
                  options = list(pageLength = dim(modelMetadata)[1]), 
                  rownames= FALSE, 
                  caption = HTML("<p><h4>Global climate models featured in this app; see the <a href='https://bcgov-env.shinyapps.io/cmip6-BC/' target='_blank'>cmip-BC app</a> for more info. 
                                 ECS is equilibrium climate sensitivity (long-term temperature change in response to an instant doubling of CO2), and values are quoted from <a href='https://advances.sciencemag.org/content/6/26/eaba1981.abstract' target='_blank'>Meehl et al. (2020)</a>. 
                                 The last two columns are the number of model runs for each scenario that are included in ClimateBC</p></h4>")
    )
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)


