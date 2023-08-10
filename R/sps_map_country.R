#' Creates a world map with target population and selected countries from SPS indicated.
#' @param out Output from function \code{sps()}.
#' @param ... Arguments passed onto \code{fuzzymatch::stringdist_join()}.
#' @import ggplot2
#' @import fuzzyjoin
#' @importFrom dplyr group_by_at slice_min
#' @importFrom countrycode countrycode
#' @return \code{sps_map_country} visualizes the selection of study sites on a world map.
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export

sps_map_country <- function(out, ...){
  # World map data
  world_map <- map_data("world")
  world_map$region <- ifelse(world_map$region %in% c('Trinidad', 'Tobago'), 'Trinidad and Tobago', world_map$region)
  world_map$iso3 <- countrycode::countrycode(world_map$region, "country.name", "iso3c",
                                             custom_match = c('Micronesia'     = 'FSM',
                                                              'Virgin Islands' = 'VIR',
                                                              'Saint Martin'   = 'MAF',
                                                              'Kosovo'         = 'XKX'), warn = FALSE)
  world_map$iso3[world_map$subregion == 'Hong Kong']  <- 'HKG'
  world_map$iso3[world_map$subregion == 'Somaliland'] <- 'SML'
  world_map$iso3[world_map$subregion == 'Zanzibar']   <- 'ZZB'
  world_map$iso3[world_map$subregion == 'British']    <- 'VGB'
  world_map$iso3[world_map$subregion == 'Gaza Strip'] <- 'PSG'

  # ISO3 matching
  if (unique(nchar(row.names(out$X)))==1 & nchar(row.names(out$X))==3){
    world_map$target <- ifelse(world_map$iso3 %in% out$selected_sites, 'Selected Sites',
                               ifelse(world_map$iso3 %in% row.names(out$X), 'Target Population', 'Rest of the World'))
  }
  else{
    # Match with country names
    fuzzy <- stringdist_join(world_map[!duplicated(world_map[,c('iso3', 'region')]) & !is.na(world_map$iso3),c('iso3', 'region')],
                             data.frame(region = row.names(out$X), selected = as.numeric(row.names(out$X) %in% out$selected_sites)),
                             by = "region",
                             distance_col = "distance",
                             mode = "inner",
                             max_dist = 20,
                             ...)
    fuzzy <- dplyr::group_by_at(fuzzy, c('region.y'))
    fuzzy <- dplyr::slice_min(fuzzy, order_by = .data[['distance']], n = 1, with_ties = FALSE)
    if (max(fuzzy$distance)>0){
      warning(paste('\nFollowing countries were matched using fuzzy-match:\n',
                    paste(paste(fuzzy$region.y[fuzzy$distance>0], '=>', fuzzy$region.x[fuzzy$distance>0]), collapse = '\n '),
                    '\nPlease review and if necessary, recode the country names in your dataset or supply the ISO3 code instead.'))
    }
    fuzzy <- fuzzy[, c('iso3', 'selected')]
    world_map <- merge(world_map, fuzzy, by = 'iso3', all.x = TRUE)
    world_map <- world_map[order(world_map$order),]
    world_map$target <- ifelse(is.na(world_map$selected), 'Rest of the World',
                               ifelse(world_map$selected == 1, 'Selected Sites', 'Target Population'))
  }

  world_map$target <- factor(world_map$target, levels = c('Selected Sites', 'Target Population', 'Rest of the World'))

  gmap <- ggplot(data = world_map, aes(map_id = world_map$region)) +
    geom_map(data = world_map,
             map  = world_map,
             aes(fill = world_map$target)) +
    scale_fill_manual(na.value = 'white', values =  c('red4', 'chartreuse4', 'gray')) +
    xlab('') + ylab('') +
    expand_limits(x = world_map$long, y = world_map$lat) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = 'bottom',
          axis.text = element_blank(),
          axis.ticks = element_blank())
  return(gmap)
}
