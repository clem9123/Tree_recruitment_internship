library(tidyverse)
library(beepr)
library(R2jags)
library(bayesplot)

Run <- function(sp){

data = readRDS("data/full_data.RDS")

# scale all variables in data
data <- data %>%
    filter(essence == sp) %>%
        select("id_pe",
            "year_measured", "latitude", "longitude", "epmatorg","ph_humus",
            "tree_ba", "partial_logging", "logging", "burn", "outbreak",
            "is_partial_logging", "is_logging", "is_species", "all_cl",
            "is_burn", "is_outbreak", "presence_gaule")

data <- data %>% mutate(
    year_measured = scale(year_measured)[,],
    latitude = scale(latitude)[,],
    longitude = scale(longitude)[,],
    epmatorg = scale(epmatorg)[,],
    ph_humus = scale(ph_humus)[,],
    tree_ba = scale(tree_ba)[,],
    partial_logging = scale(partial_logging)[,],
    logging = scale(logging)[,],
    burn = scale(burn)[,],
    outbreak = scale(outbreak)[,]
) %>% mutate(
    logging = ifelse(is.na(logging), 0, logging),
    partial_logging = ifelse(is.na(partial_logging), 0, partial_logging),
    burn = ifelse(is.na(burn), 0, burn),
    outbreak = ifelse(is.na(outbreak), 0, outbreak),
    is_partial_logging = ifelse(is.na(is_partial_logging), 0, is_partial_logging),
    is_logging = ifelse(is.na(is_logging), 0, is_logging),
    is_burn = ifelse(is.na(is_burn), 0, is_burn),
    is_outbreak = ifelse(is.na(is_outbreak), 0, is_outbreak)
) %>% mutate(
    id_pe =  as.numeric(as.factor(data$id_pe))
)

data <- data %>% na.omit()

jags_data <- list(
    N = nrow(data),
    pa = data$presence_gaule,
    nb = data$all_cl,
    year = data$year_measured,
    latitude = data$latitude,
    longitude = data$longitude,
    epmatorg = data$epmatorg,
    ph_humus = data$ph_humus,
    total_ba = data$tree_ba,
    is_species = data$is_species,
    TSpl = data$partial_logging,
    TSl = data$logging,
    TSb = data$burn,
    TSo = data$outbreak,
    is_pl = data$is_partial_logging,
    is_l = data$is_logging,
    is_b = data$is_burn,
    is_o = data$is_outbreak,
    placette = data$id_pe,
    n_placettes = max(data$id_pe)
)

parameters = c(
    # Intercept
    "pa_intercept",
    "nb_intercept",
    # Spatial-temporal variables
    "pa_t", "pa_lat", "pa_lon", "pa_interaction_lat_t",
    "nb_t", "nb_lat", "nb_lon", "nb_interaction_lat_t",
    # Soil variables
    "pa_emo", "pa_ph", "pa_ph2",
    "nb_emo", "nb_ph", "nb_ph2",
    # Biotic variables
    "pa_ba", "pa_sp",
    "nb_ba", "nb_sp",
    # Perturbation
    "pa_beta_pl", "pa_TSD_pl", "pa_TSD2_pl",
    "pa_beta_l", "pa_TSD_l", "pa_TSD2_l",
    "pa_beta_b", "pa_TSD_b", "pa_TSD2_b",
    "pa_beta_o", "pa_TSD_o", "pa_TSD2_o",
    "nb_beta_pl", "nb_TSD_pl", "nb_TSD2_pl",
    "nb_beta_l", "nb_TSD_l", "nb_TSD2_l",
    "nb_beta_b", "nb_TSD_b", "nb_TSD2_b",
    "nb_beta_o", "nb_TSD_o", "nb_TSD2_o",
    # Random effect of placette
    "tau1", "tau2"
)

# temps au début
begin = Sys.time()

out <- jags.parallel(
    model.file = "models/model_pert2.txt",
    data = jags_data,
    parameters.to.save = parameters,
    n.chains = 3
)

# temps à la fin
end = Sys.time()

# temps d'exécution
Tex <- end - begin

out$runtime <- Tex

save(out, file = paste("output/", sp, ".RData", sep = ""))

}

### LES ESPECES ###

#tbl <- data %>% group_by(essence) %>% summarize(nb = sum(presence_gaule))
#
#sp_type <- read.csv("Documentation/Tree_sp.csv")
#
#full <- merge(tbl, sp_type)
#
#ggplot(full %>% filter(nb >0)) +
#geom_bar(aes(x = essence, y = nb, fill = zone), stat = "identity")

# ajouter boréale ou tepérée

# BOREALES
Run("SAB")
Run("ERE")
Run("EPN")
Run("EPB")
# +
Run("AUR")
Run("ERP")
Run("SAL")
Run("PIG")

# TEMPEREES
Run("BOP")
Run("ERR")
Run("ERS")
Run("BOJ")
# +
Run("PET")
Run("PRP")
Run("SOA")
Run("THO")

# All BOREALES

data = readRDS("data/full_data.RDS")

# scale all variables in data
data <- data %>%
    filter(essence %in% c("SAB", "ERE", "EPN", "EPB", "AUR", "ERP", "SAL", "PIG")) %>%
        select("id_pe",
            "year_measured", "latitude", "longitude", "epmatorg","ph_humus",
            "tree_ba", "partial_logging", "logging", "burn", "outbreak",
            "is_partial_logging", "is_logging", "is_species", "all_cl",
            "is_burn", "is_outbreak", "presence_gaule")

data <- data %>% mutate(
    year_measured = scale(year_measured)[,],
    latitude = scale(latitude)[,],
    longitude = scale(longitude)[,],
    epmatorg = scale(epmatorg)[,],
    ph_humus = scale(ph_humus)[,],
    tree_ba = scale(tree_ba)[,],
    partial_logging = scale(partial_logging)[,],
    logging = scale(logging)[,],
    burn = scale(burn)[,],
    outbreak = scale(outbreak)[,]
) %>% mutate(
    logging = ifelse(is.na(logging), 0, logging),
    partial_logging = ifelse(is.na(partial_logging), 0, partial_logging),
    burn = ifelse(is.na(burn), 0, burn),
    outbreak = ifelse(is.na(outbreak), 0, outbreak),
    is_partial_logging = ifelse(is.na(is_partial_logging), 0, is_partial_logging),
    is_logging = ifelse(is.na(is_logging), 0, is_logging),
    is_burn = ifelse(is.na(is_burn), 0, is_burn),
    is_outbreak = ifelse(is.na(is_outbreak), 0, is_outbreak)
) %>% mutate(
    id_pe =  as.numeric(as.factor(data$id_pe))
)

data <- data %>% na.omit()

jags_data <- list(
    N = nrow(data),
    pa = data$presence_gaule,
    nb = data$all_cl,
    year = data$year_measured,
    latitude = data$latitude,
    longitude = data$longitude,
    epmatorg = data$epmatorg,
    ph_humus = data$ph_humus,
    total_ba = data$tree_ba,
    is_species = data$is_species,
    TSpl = data$partial_logging,
    TSl = data$logging,
    TSb = data$burn,
    TSo = data$outbreak,
    is_pl = data$is_partial_logging,
    is_l = data$is_logging,
    is_b = data$is_burn,
    is_o = data$is_outbreak,
    placette = data$id_pe,
    n_placettes = max(data$id_pe)
)

parameters = c(
    # Intercept
    "pa_intercept",
    "nb_intercept",
    # Spatial-temporal variables
    "pa_t", "pa_lat", "pa_lon", "pa_interaction_lat_t",
    "nb_t", "nb_lat", "nb_lon", "nb_interaction_lat_t",
    # Soil variables
    "pa_emo", "pa_ph", "pa_ph2",
    "nb_emo", "nb_ph", "nb_ph2",
    # Biotic variables
    "pa_ba", "pa_sp",
    "nb_ba", "nb_sp",
    # Perturbation
    "pa_beta_pl", "pa_TSD_pl", "pa_TSD2_pl",
    "pa_beta_l", "pa_TSD_l", "pa_TSD2_l",
    "pa_beta_b", "pa_TSD_b", "pa_TSD2_b",
    "pa_beta_o", "pa_TSD_o", "pa_TSD2_o",
    "nb_beta_pl", "nb_TSD_pl", "nb_TSD2_pl",
    "nb_beta_l", "nb_TSD_l", "nb_TSD2_l",
    "nb_beta_b", "nb_TSD_b", "nb_TSD2_b",
    "nb_beta_o", "nb_TSD_o", "nb_TSD2_o",
    # Random effect of placette
    "tau1", "tau2"
)

# temps au début
begin = Sys.time()

out <- jags.parallel(
    model.file = "models/model_pert2.txt",
    data = jags_data,
    parameters.to.save = parameters,
    n.chains = 3
)

# temps à la fin
end = Sys.time()

# temps d'exécution
Tex <- end - begin

out$runtime <- Tex

save(out, file = paste("output/", "BOREALES", ".RData", sep = ""))

# All TEMPEREES

data = readRDS("data/full_data.RDS")

# scale all variables in data
data <- data %>%
    filter(essence %in% c("BOP", "ERR", "ERS", "BOJ", "PET", "PRP", "SOA", "THO")) %>%
        select("id_pe",
            "year_measured", "latitude", "longitude", "epmatorg","ph_humus",
            "tree_ba", "partial_logging", "logging", "burn", "outbreak",
            "is_partial_logging", "is_logging", "is_species", "all_cl",
            "is_burn", "is_outbreak", "presence_gaule")

data <- data %>% mutate(
    year_measured = scale(year_measured)[,],
    latitude = scale(latitude)[,],
    longitude = scale(longitude)[,],
    epmatorg = scale(epmatorg)[,],
    ph_humus = scale(ph_humus)[,],
    tree_ba = scale(tree_ba)[,],
    partial_logging = scale(partial_logging)[,],
    logging = scale(logging)[,],
    burn = scale(burn)[,],
    outbreak = scale(outbreak)[,]
) %>% mutate(
    logging = ifelse(is.na(logging), 0, logging),
    partial_logging = ifelse(is.na(partial_logging), 0, partial_logging),
    burn = ifelse(is.na(burn), 0, burn),
    outbreak = ifelse(is.na(outbreak), 0, outbreak),
    is_partial_logging = ifelse(is.na(is_partial_logging), 0, is_partial_logging),
    is_logging = ifelse(is.na(is_logging), 0, is_logging),
    is_burn = ifelse(is.na(is_burn), 0, is_burn),
    is_outbreak = ifelse(is.na(is_outbreak), 0, is_outbreak)
) %>% mutate(
    id_pe =  as.numeric(as.factor(data$id_pe))
)

data <- data %>% na.omit()

jags_data <- list(
    N = nrow(data),
    pa = data$presence_gaule,
    nb = data$all_cl,
    year = data$year_measured,
    latitude = data$latitude,
    longitude = data$longitude,
    epmatorg = data$epmatorg,
    ph_humus = data$ph_humus,
    total_ba = data$tree_ba,
    is_species = data$is_species,
    TSpl = data$partial_logging,
    TSl = data$logging,
    TSb = data$burn,
    TSo = data$outbreak,
    is_pl = data$is_partial_logging,
    is_l = data$is_logging,
    is_b = data$is_burn,
    is_o = data$is_outbreak,
    placette = data$id_pe,
    n_placettes = max(data$id_pe)
)

parameters = c(
    # Intercept
    "pa_intercept",
    "nb_intercept",
    # Spatial-temporal variables
    "pa_t", "pa_lat", "pa_lon", "pa_interaction_lat_t",
    "nb_t", "nb_lat", "nb_lon", "nb_interaction_lat_t",
    # Soil variables
    "pa_emo", "pa_ph", "pa_ph2",
    "nb_emo", "nb_ph", "nb_ph2",
    # Biotic variables
    "pa_ba", "pa_sp",
    "nb_ba", "nb_sp",
    # Perturbation
    "pa_beta_pl", "pa_TSD_pl", "pa_TSD2_pl",
    "pa_beta_l", "pa_TSD_l", "pa_TSD2_l",
    "pa_beta_b", "pa_TSD_b", "pa_TSD2_b",
    "pa_beta_o", "pa_TSD_o", "pa_TSD2_o",
    "nb_beta_pl", "nb_TSD_pl", "nb_TSD2_pl",
    "nb_beta_l", "nb_TSD_l", "nb_TSD2_l",
    "nb_beta_b", "nb_TSD_b", "nb_TSD2_b",
    "nb_beta_o", "nb_TSD_o", "nb_TSD2_o",
    # Random effect of placette
    "tau1", "tau2"
)

# temps au début
begin = Sys.time()

out <- jags.parallel(
    model.file = "models/model_pert2.txt",
    data = jags_data,
    parameters.to.save = parameters,
    n.chains = 3
)

# temps à la fin
end = Sys.time()

# temps d'exécution
Tex <- end - begin

out$runtime <- Tex

save(out, file = paste("output/", "TEMPEREES", ".RData", sep = ""))