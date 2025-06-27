library(dplyr)
library(tidyselect)

# Adjusted prevalence function -------------------------------------------------
prevalence <- function(
    data, stratify_by, adjust_for=NULL, n_distinct=NULL, denominator=NULL) {
  
  stopifnot(  # Argument checks ------------------------------------------------
    is.data.frame(data),
    is.character(stratify_by) & stratify_by %in% names(data),
    is.character(adjust_for) & adjust_for %in% names(data) & !adjust_for %in% stratify_by,
    is.null(denominator) | (is.character(denominator) & denominator %in% stratify_by),
    is.null(n_distinct) | (is.character(n_distinct) & adjust_for %in% names(data) & !n_distinct %in% stratify_by)
  )
  
  # Set denominator as first stratification variable if it is not set ----------
  if(length(stratify_by) > 1){
    if(is.null(denominator)){denominator <- stratify_by[1]}
  } else {
    if(!is.null(denominator)){
      warning("Only stratifying by 1 variable, setting denominator to NULL")
      denominator <- NULL
    }
  }
  return(  # Main algorithm ----------------------------------------------------
   data |>
     dplyr::mutate(
       .by=tidyselect::all_of(c(denominator)),  # Group by denominator ---------
       
       # Calculate denominator values ------------------------------------------
       denominator.raw = dplyr::n(),
       denominator.adj = if(!is.null(adjust_for)) {
         dplyr::n_distinct(dplyr::across(tidyselect::all_of(c(adjust_for))))
       } else {
         NULL
       }
     ) |>
     {\(.) if (!is.null(adjust_for)) dplyr::mutate(
       .,
       .by=tidyselect::all_of(c(stratify_by)),  # Group by strata --------------
       
       # Calculate adjusted conditionally --------------------------------------
       count.adj=dplyr::n_distinct(dplyr::across(tidyselect::all_of(c(adjust_for)))),
       prop.adj=count.adj/denominator.adj,
       se.adj=sqrt(prop.adj*(1-prop.adj)/denominator.adj),
       lower.adj=prop.adj-(1.96*se.adj),
       upper.adj=prop.adj+(1.96*se.adj)
       
       ) else . }() |> 
     dplyr::mutate(
       .by=tidyselect::all_of(c(stratify_by)),  # Group by strata --------------
       
       # Calculate raw ---------------------------------------------------------
       count.raw=dplyr::n(), 
       prop.raw=count.raw/denominator.raw,
       se.raw=sqrt(prop.raw*(1-prop.raw)/denominator.raw),
       lower.raw=prop.raw-(1.96*se.raw),
       upper.raw=prop.raw+(1.96*se.raw),
       
       # Restrict intervals to between 0-1, note, min()/max() cause errors here
       dplyr::across(
         tidyselect::matches('^(prop|upper|lower)'),
         ~dplyr::case_when(.x<0~0, .x>1~1, .default=.x)
       ),
       
       # Calculate n distinct --------------------------------------------------
       dplyr::across(tidyselect::any_of(c(n_distinct)), dplyr::n_distinct,
                     .names='# {.col}')
     ) |> 
     ## Here we are 'artificially' summarising by selecting the first
     ## line per stratification group. This is more efficient than
     ## dplyr::distinct, and we know all values will be the same for the group.
     dplyr::slice_head(n=1, by=tidyselect::all_of(c(stratify_by))) |> 
    
     # Sort data ---------------------------------------------------------------
     {\(.) if (!is.null(adjust_for))
       dplyr::arrange(., -denominator.adj, -count.adj) else
         dplyr::arrange(., -denominator.raw, -count.raw)}() |>

     # Calculate ranks ---------------------------------------------------------
     ## Here the ranks refer to the top level, i.e. denominator
     dplyr::mutate(
       .by=tidyselect::all_of(c(denominator)),
       dplyr::across(
         tidyselect::matches("prop.(raw|adj)"),
         ~dplyr::row_number(dplyr::desc(.x)),
         .names="{stringr::str_replace(.col, 'prop', 'rank')}"
       )
     ) |>
     # Add a denominator column, if multiple groups ----------------------------
     dplyr::mutate(denominator=denominator) |> 
     
     # Sort columns ------------------------------------------------------------
     dplyr::select(
       tidyselect::all_of(c(stratify_by)),
       tidyselect::matches("\\.(raw|adj)"),
       denominator,
       tidyselect::any_of(paste('#', n_distinct))
      )
  )
}
