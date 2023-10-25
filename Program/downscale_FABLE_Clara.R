

require(tidyr)
require(dplyr)
library(progress)
#devtools::install_github("tkrisztin/downscalr",ref="HEAD")
library(downscalr)
library(readr)


## load prior input (should be there already), betas and targets
load("./Outputs_Rdata/prior_input_FABLE.Rdata")
load("./Outputs_Rdata/betas_FABLE.Rdata")
load("./Outputs_Rdata/targets_FABLE.Rdata")


targets <- targets %>% filter(lu.from != lu.to)
##want to rename the targets at  some point ! 

#targets$value[which(targets$year==2020)] <- 0
## changed times to years as Warning message: Unknown or uninitialised column: `times`. 


lu_levels <- start_map
##### define Ns based on targets...
lu_levels$Ts <- 2020


#lu_levels_tot <- lu_levels %>% pivot_wider(id_cols = c(ns,Ts), names_from = lu.class, values_from = value)

lu_levels_tot <- lu_levels # %>% pivot_wider(id_cols = c(ns,Ts), names_from = lu.class, values_from = value)

lu_levels_tot[is.na(lu_levels_tot)] <- 0
lu_levels <- lu_levels_tot %>% pivot_longer(cols = !c(ns,Ts),names_to = "lu.from", values_to = "value")


lu_levels <- lu_levels %>% left_join(.,mapping[, c("uniqueID","country")], by=c("ns"="uniqueID")) %>% rename(Ns=country)

## rename year and country to match following code 
colnames(targets)[2] ="Ns"
colnames(targets)[1] ="times"

Ns <- unique(targets$Ns)
#Ns <- Ns[-which(is.na(Ns))]



## edited for UK? or removed? have removed because I dont need to filter as the targets are just UK
### targets <- targets %>% filter(!Ns=="UK", times>=2020)

## I dont have priors (Created in compile_priors) so dont need the code
##priors$ns <- as.character(priors$ns)
# start.areas = start.areas
lu_levels = lu_levels %>% dplyr::filter(Ts == 2020) %>% 
  dplyr::select(-c(Ts))



# areas.sum_to = function(res, curr.areas, priors, xmat, xmat.proj) {
#   curr.areas = res$out.res %>%
#     group_by(.data$ns,.data$lu.to) %>%
#     summarize(value = sum(.data$value),.groups = "keep") %>%
#     rename("lu.from" = "lu.to")
#   # correct for small numerical mistakes
#   if (min(curr.areas$value) < 0 && min(curr.areas$value) > -10^-10) {
#     curr.areas$value[curr.areas$value < 0] = 0
#   }
#   return(curr.areas)
# }
# 
# 
# 
# areas.identity = function(res, curr.areas, priors, xmat, xmat.proj) {
#   curr.areas = res$out.res
#   # correct for small numerical mistakes
#   if (min(curr.areas$value) < 0 && min(curr.areas$value) > -10^-10) {
#     curr.areas$value[curr.areas$value < 0] = 0
#   }
#   return(curr.areas)
# }




xmat.dyn.fun = function(res, curr.areas, priors, xmat, xmat.proj) {
  tmp.mat = curr.areas %>% pivot_wider(names_from = lu.from,values_from = value)
  tmp.mat = tmp.mat[match(row.names(xmat),tmp.mat$ns),-1]
  
  tmp.mat = tmp.mat %>% left_join()
  return(tmp.mat)
}

est.grid = expand.grid(Ns = Ns)

pb <- progress_bar$new(
  format = "  [:bar] :percent in :eta",
  total = nrow(est.grid), clear = FALSE, width= 60)



classes <- c()#c("Urban","NotRel")#, "Rives and lakes", "Marine inlets and transitional waters", "Marine", "Wetlands" )

iter <- 1
for (iter in 1:nrow(est.grid)) {
  curr.N = est.grid$Ns[iter]
  
  curr.ns <- as.character(mapping$uniqueID[mapping$country==curr.N])
  ##curr.region <- unique(mapping$REGION_37[mapping$country==curr.N])
  curr.country <- unique(mapping$country[mapping$country==curr.N])
  curr.targets <- targets[targets$Ns==curr.N,] %>% dplyr::select(-Ns)
  
  
  curr.targets <- curr.targets %>% filter(!lu.from%in%classes&!lu.to%in%classes)
  
  # curr.targets$lu.from <- gsub(" ", ".", curr.targets$lu.from)
  # curr.targets$lu.to <- gsub(" ", ".", curr.targets$lu.to)
  
  
  
  
  curr.betas <- betas[betas$Ns==curr.country,] %>% dplyr::select(!c(Ns))
  #curr.betas$lu.from <- gsub(" ", ".", curr.betas$lu.from)
  curr.betas <- curr.betas %>% filter(!lu.from%in%classes,!lu.to%in%classes)
  curr.betas <- data.frame(curr.betas)
  
 # curr.priors <- priors[priors$Ns==curr.N,] %>% dplyr::select(!c(Ns,REGION_37))
  #curr.priors$ns <- as.character(curr.priors$ns)
  # curr.priors <- data.frame(curr.priors)
  
  
  X.temp <- X[[curr.country]]
  curr.X <- X.temp %>% filter(rownames(X.temp)%in%curr.ns)
  curr.X$ns <- as.character(curr.ns)
  curr.X <- curr.X %>% pivot_longer(cols = !c(ns), names_to = "ks", values_to = "value")
  curr.X <- data.frame(curr.X)
  
  curr.lu_levels <- lu_levels %>% filter(Ns==curr.N, lu.from!="NODATA",!lu.from%in%classes) %>% dplyr::select(-Ns)
  #curr.lu_levels <- na.omit(curr.lu_levels)
  # curr.lu_levels <- subset(curr.lu_levels,curr.lu_levels$lu.from!="NODATA")
  # curr.lu_levels <- subset(curr.lu_levels,!curr.lu_levels$lu.from%in%classes)
  # curr.lu_levels$lu.from <- gsub(" ", ".", curr.lu_levels$lu.from)
  curr.lu_levels$ns <- as.character(curr.lu_levels$ns)
  curr.lu_levels <- data.frame(curr.lu_levels)
  
  
  # xmatcols <- data.frame(ks=colnames(X.temp)[-ncol(X.temp)], value=c("static",rep("dynamic",2),"static","static",rep("dynamic",4),rep("static",ncol(X.temp) - 10)))
  
  
  
  ###
  res1 = downscale(targets = curr.targets, 
                   start.areas = curr.lu_levels,
                   xmat = curr.X,
                   betas = curr.betas) # ,priors=curr.priors)

  
  ## out.res is output of solve_biascorr.R - for which you call downscale function
  ## chk targets
  
  
  ttemp.res <- res1$out.res %>% filter(value!=0)
  ttemp.agg <- res1$out.res %>% group_by(lu.to, times, ns) %>% summarise(value = sum(value), .groups = 'drop')
  
  
  chck.targets =   curr.targets %>%
    left_join(
      res1$out.res %>%
        group_by(lu.from,lu.to,times) %>%
        summarize(downscale.value = sum(value),.groups = "keep"),by = c("lu.from", "lu.to","times") ) %>%
    mutate(diff = value - downscale.value)
  
  
  chck.targets <- chck.targets %>% filter(diff>=1)
  fault.years <- unique(chck.targets$times)
  
  if(nrow(ttemp.res)!=0){
    ttemp.res <- data.frame(ttemp.res, downscaling_area_flag=0)
    ttemp.res$downscaling_area_flag[ttemp.res$times%in%fault.years] <- 1
  }
  
  ttemp.agg <- data.frame(ttemp.agg, downscaling_area_flag=0)
  ttemp.agg$downscaling_area_flag[ttemp.agg$times%in%fault.years] <- 1
  
  
  if(iter==1){
    luc.res <- ttemp.res
    lu.res <- ttemp.agg
  } else {
    luc.res <- rbind(luc.res,ttemp.res)
    lu.res <- rbind(lu.res,ttemp.agg)
  }
  
  pb$tick()
}


# mapping$uniqueID <- as.character(mapping$uniqueID)
# luc.res.colrow <- left_join(luc.res, mapping, by=c("ns"="uniqueID")) %>% rename(Ns=country)
luc.res_check <- left_join(luc.res, mapping[,c("uniqueID","country")], by=c("ns"="uniqueID")) %>% dplyr::select(country, ns, times,lu.from,lu.to,value, downscaling_area_flag) %>% rename(Ns="country")


chck.targets.tot.Ns =   targets %>%
  left_join(
    luc.res_check %>%
      group_by(lu.from,lu.to,times, Ns) %>%
      summarize(downscale.value = sum(value, na.rm=T),.groups = "keep"),by = c("lu.from", "lu.to","times", "Ns") ) %>%
  mutate(diff = value - downscale.value) %>% filter(!is.na(diff))


chck.targets.tot.LUC = chck.targets.tot.Ns %>%
  group_by(lu.from,lu.to,times) %>%
  summarise(value=sum(value),
            downscale.value=sum(downscale.value),
            diff=sum(diff))


chck.targets.tot.Ns$percent_mismatch = chck.targets.tot.Ns$downscale.value/chck.targets.tot.Ns$value
chck.targets.tot.LUC$percent_mismatch = chck.targets.tot.LUC$downscale.value/chck.targets.tot.LUC$value
# 
# write.csv(chck.targets.tot, file = "missed_targets_colrow_15feb.csv", row.names = FALSE)
# write.csv(chck.targets.tot.indo, file = "missed_targets_indo_15feb.csv", row.names = FALSE)

write.csv(chck.targets.tot.LUC, file = "./Outputs/check_targets_LUC_291122.csv", row.names = FALSE)
write.csv(chck.targets.tot.Ns, file = "./Outputs/check_targets_Ns_291122.csv", row.names = FALSE)

## Because FABLE does not have a final classification, the results are contained in lu.res (area levels) and luc.res(flows)
## so skip the first set and go to ## add col row below (added in here for ease!)

# add col row mapping 

luc.res_fin <- left_join(luc.res, mapping[,c("uniqueID","country")], by=c("ns"="uniqueID")) %>% dplyr::select(country, ns, times,lu.from,lu.to,value, downscaling_area_flag)
lu.res_fin <- left_join(lu.res, mapping[,c("uniqueID","country")], by=c("ns"="uniqueID")) %>% dplyr::select(country, ns, times,lu.to,value, downscaling_area_flag)


colnames(luc.res_fin) <- c("country", "ns", "year","lu.from","lu.to","value","downscaling.flag")
colnames(lu.res_fin) <- c("country", "ns", "year","lu.class","value","fivearcmnArea")

# luc.res_ETL2 <- left_join(luc.res,share_ETL2_map, by=c("ns","lu.to"="GLOB_LUC_class")) %>% mutate(value=value*init_share) %>% 
#   group_by(times,ns,Ecosystem_types_level_2) %>% summarize(value=sum(value)) %>% drop_na()

## Don't think I need to do this as dont have ETL2 map
##lu.res_ETL2 <- left_join(lu.res,share_ETL2_map, by=c("ns","lu.to"="GLOB_LUC_class")) %>% mutate(value=value*init_share) %>% 
 ##group_by(times,ns,!!rlang::sym(LUC_final_classification)) %>% summarize(value=sum(value)) %>% drop_na()



#### add colrow mapping
#luc.res_fin <- left_join(luc.res, mapping[,c("uniqueID","country")], by=c("ns"="uniqueID")) %>% dplyr::select(country, ns, times,lu.from,lu.to,value, downscaling_area_flag)
##lu.res_fin <- left_join(lu.res_ETL2, mapping[,c("uniqueID","country", "fivearcmnArea")], by=c("ns"="uniqueID")) %>% dplyr::select(country, ns, times, !!rlang::sym(LUC_final_classification), value, fivearcmnArea)

#colnames(luc.res_fin) <- c("country", "ns", "year","lu.from","lu.to","value","downscaling.flag")
#colnames(lu.res_fin) <- c("country", "ns", "year","lu.class","value","fivearcmnArea")


#asdf <- lu.res_fin %>% group_by(ns,times) %>% summarise(area_covered=sum(value)/fivearcmnArea)
# 
# 
# 
# 
# # 
# # 
# # 
# # 
# # 
# # test.restrictions = res1$out.res %>% mutate(ns = as.double(ns)) %>%
# #   left_join(restrictions %>% rename("restr" = "value"),by = c("lu.from", "ns", "lu.to"))  %>%
# #   group_by(times, lu.from, lu.to) %>%
# #   filter(restr > 0) %>%
# #   summarise(value = sum(value),restr = sum(restr),.groups = "keep")
# # cat("Restrictions: ",all(test.restrictions$value == 0),"\n")
# # 
# # 
# require(raster)
# geosims <- raster("../Downscalin+g_all_models/Data/simu_raster/w001001.adf")
#  sp::proj4string(geosims)<-sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #apply projection manually
#  
#  mapping_table <- read_csv("inputs/full_simu_map_biodiv.csv")
#  
#  #check which regions
#  #unique(mapping_table$REGION_37)
#  
#  eu.regions <- c("EU_South", "EU_North", "EU_MidWest", "EU_Baltic", "RCEU", "ROWE")
#  
#  mapping_table <- subset(mapping_table, mapping_table$REGION_37%in%eu.regions)
#  excluding.countries <- c("Greenland")
#  mapping_table <- subset(mapping_table, mapping_table$country!=excluding.countries)

# save_geovals <- data.frame(getValues(geosims)) #just the SimU IDs
# colnames(save_geovals)<-"SimUID"
# 
# targets %>% group_by(lu.from) %>% summarise(value=sum(value))
