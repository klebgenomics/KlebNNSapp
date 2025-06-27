library(readr)
library(fs)
library(dplyr)
library(sf)
library(stringr)
library(leaflet)
library(shinyWidgets)
library(shinycssloaders)

source('src/R/kleborate.R')

# Define regexes ---------------------------------------------------------------
LOCUS_TYPE_REGEX <- stringr::regex('[KO][L]?[0-9/vOKabcfg⍺β𝛾,]{1,6}') # Misses -D1
YEAR_REGEX <- stringr::regex('\\b(19|20)\\d{2}\\b') # Realistic year regex

# Data globals -----------------------------------------------------------------
DATA <- fs::path('data/TableS3_IsolatesIncluded_NEEDSACCESSIONS.tsv') |> 
  readr::read_tsv(show_col_types=FALSE) |> 
  dplyr::mutate(O_type=O_genotype)
# Manually changed ‘GBS’ to ‘GBS-COP’

convert_o_type <- function(.x) {
  return(  # From Kat
    dplyr::case_match(
      .x=.x, .default=.x,
      "O1.v1" ~ "O1⍺β,2⍺",
      "O1.v2" ~ "O1⍺β,2β",
      "O1.v3" ~ "O1⍺β,2⍺𝛾",
      "O2a.v1" ~ "O2⍺",
      "O2a.v3" ~ "O2⍺𝛾",
      "O2afg.v2" ~ "O2β",
      "O3/O3a" ~ "O3⍺/O3β",
      "O3b" ~ "O3𝛾",
      "OL13" ~ "O13",
      "OL102" ~ "O14",
      "OL103" ~ "O10",
      "OL104" ~ "O15"
    )
  )
}

MODELLED <- fs::dir_ls('data/modelled', glob='*global_estimates.csv') |>
  readr::read_csv(id='id', show_col_types=FALSE) |>
  dplyr::mutate(id=fs::path_ext_remove(fs::path_file(id)), subgroup=NULL, mean=NULL) |>
  tidyr::separate_wider_delim(
    id, '_', names=c('Antigen', 'Infection', 'Category', 'raw_adj'),
    too_many='drop'
  ) |>
  # dplyr::filter(Category=='min10') |>
  dplyr::select(!Category, prop=median) |>
  tidyr::pivot_wider(
    names_from=raw_adj, values_from=c(prop, lower, upper), names_sep='.') |> 
  dplyr::mutate(.by=c(1, 2), denominator.adj=sum(prop.adj), denominator.raw=sum(prop.raw)) |>
  dplyr::arrange(Antigen, Infection, -denominator.adj, -prop.adj) |>
  dplyr::mutate(
    .by=c(1, 2),
    rank.raw=dplyr::row_number(dplyr::desc(denominator.raw)),
    rank.adj=dplyr::row_number(dplyr::desc(denominator.adj))
  ) |> 
  dplyr::mutate(
    Antigen=dplyr::if_else(Antigen=='K', 'K_locus', 'O_type'),
    locus=convert_o_type(locus),
    Infection=dplyr::case_match(
      Infection, .default=Infection,
      'Carba' ~ 'Carbapenemase+', 
      'ESBL' ~ 'ESBL+', 
      'Full' ~ 'All samples',
      'Fatal' ~ 'Fatal cases'
    )
  )

MODELLED_STRATIFIED <- fs::dir_ls('data/modelled/', glob='*subgroup_estimates.csv') |>
  readr::read_csv(id='id', show_col_types=FALSE) |>
  dplyr::mutate(id=fs::path_ext_remove(fs::path_file(id)), mean=NULL) |>
  tidyr::separate_wider_delim(
    id, '_', names=c('Antigen', 'Infection', 'Category', 'raw_adj'),
    too_many='drop'
  ) |>
  # dplyr::filter(Category=='min10') |>
  dplyr::select(Antigen, Infection, Region=subgroup, raw_adj, locus, prop=median, lower, upper) |>
  tidyr::pivot_wider(
    names_from=raw_adj, values_from=c(prop, lower, upper), names_sep='.') |> 
  dplyr::mutate(.by=c(1, 2, 3), denominator.adj=sum(prop.adj), denominator.raw=sum(prop.raw)) |>
  dplyr::arrange(Antigen, Infection, -denominator.adj, -prop.adj) |>
  dplyr::mutate(
    .by=c(1, 2),
    rank.raw=dplyr::row_number(dplyr::desc(denominator.raw)),
    rank.adj=dplyr::row_number(dplyr::desc(denominator.adj))
  ) |>
  dplyr::mutate(
    denominator='Region', 
    Antigen=dplyr::if_else(Antigen=='K', 'K_locus', 'O_type'),
    locus=convert_o_type(locus),
    Infection=dplyr::case_match(
      Infection, .default=Infection,
      'Carba' ~ 'Carbapenemase+', 
      'ESBL' ~ 'ESBL+', 
      'Full' ~ 'All samples',
      'Fatal' ~ 'Fatal cases'
    )
  )

# Variable globals -------------------------------------------------------------
KLEBORATE_COLS <- names(kleborate_column_spec$cols)
KLEBORATE_VARIABLES <- c("ST", "K_locus", "K_type", "O_locus", "O_type",
                         "resistance_score")
X_AXIS_VARIABLES <- c(
  'Country', 'Year', 'Region', 'Study', KLEBORATE_VARIABLES)
COUNTRIES <- unique(DATA$Country)
YEARS <- unique(DATA$Year)
MIN_YEARS <- min(YEARS)
MAX_YEARS <- max(YEARS)
REGIONS <- unique(DATA$Region)
STUDIES <- unique(DATA$Study)
SITES <- unique(DATA$Site)
# SPECIES <- unique(DATA$species)
ST <- unique(DATA$ST)
KL <- unique(DATA$K_locus)
OL <- unique(DATA$O_locus)
AMR <- c('All samples', 'ESBL+', 'Carbapenemase+')
MODELLED_INFECTION_TYPES <- unique(MODELLED$Infection)

# Map globals ------------------------------------------------------------------
WORLD_SF <- sf::read_sf('data/world-administrative-boundaries') |> 
  dplyr::rename('Region'='region', 'Country'='name')
SF_COUNTRIES <- unique(WORLD_SF$Country)
SF_REGIONS <- unique(WORLD_SF$Region)
MAP_PALETTE <- leaflet::colorBin(
  "viridis", seq(0, 1), na.color="transparent", bins=seq(0, 1, length.out=9),
  reverse=TRUE
)

# Cosmetic globals -------------------------------------------------------------
HI_FILL <- '#1E88E5' #D81B60
LO_FILL <- '#FFC107'
BAR_FILL <- '#D1C1E1FF'
## CSS Loader Spinner Formatting: https://projects.lukehaas.me/css-loaders/
SPINNER_TYPE <- 8
SPINNER_COLOR <- BAR_FILL
SPINNER_SIZE <- .9
SIDEBAR_WIDTH <- 320
DEFAULT_VALENCY <- 20

## Wrapper function for the hover dropdown menu to save lines of repetitive code
DROPDOWN_FUNCTION <- function(...) {
  shinyWidgets::dropdown(
    ...,
    label=NULL,
    style="jelly", 
    icon=shiny::icon("wrench"),
    status="danger",
    size="sm",
    tooltip=shinyWidgets::tooltipOptions(title='Configure plot'),
    animate=shinyWidgets::animateOptions(
      enter=shinyWidgets::animations$fading_entrances$fadeInUp,
      exit=shinyWidgets::animations$fading_exits$fadeOutDown,
      duration=0.2
    )
  )
}
