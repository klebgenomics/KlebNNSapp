library(dplyr)
library(tidyr)
library(tidyselect)
library(ggplot2)

# Pyramid plotting function ----------------------------------------------------
pyramid_plot <- function(
    data, y_order=NULL, low_col="#FFC107", high_col="#1E88E5", fill='count'){
  
  stopifnot(  # Argument checks
    is.data.frame(data),
    'prop.raw' %in% names(data),
    !'denominator' %in% names(data),  # Data is stratified by >1 group!
    is.null(y_order) | is.character(y_order),
    is.character(low_col), is.character(high_col),
    fill %in% c('count', 'prop', 'rank')
  )
  
  y_index <- 1
  y_column <- names(data)[y_index]
  y_column_plural <- plural_name(y_column, nice=TRUE)
  data[[y_index]] <- as.character(data[[y_index]])  # Force character vector
  if(is.null(y_order)) {y_order <- unique(data[[y_index]])}
  
  return(
    data |>
      tidyr::pivot_longer(
        tidyselect::matches('.(raw|adj)$'),
        names_to=c(".value", "type"),  # type is either raw or adj
        names_sep="\\."  # replace full-stop with space
      ) |> 
      dplyr::mutate(
        type=factor(
          dplyr::if_else(type=='raw', 'Raw proportion', 'Adjusted proportion'), 
          c('Raw proportion', 'Adjusted proportion') # Ensure adjusted is on the right
        ),  # Make raw values negative so they appear on the left 
        dplyr::across(c(prop, lower, upper), ~dplyr::if_else(type=='Raw proportion', -.x, .x)),
        hover_text=paste0(
          nice_name(y_column), ": ", 
          .data[[y_column]], ' (rank = ', rank, ")", "\n",
          type, ": ", sprintf("%.3f", abs(prop)), "\n",
          "Upper: ", sprintf("%.3f", abs(upper)), "\n",
          "Lower: ", sprintf("%.3f", abs(lower)), "\n",
          'Fill = ', nice_name(fill), ' (', sprintf("%.3f", abs(.data[[fill]])), ')'
        )
      ) |> 
      ggplot2::ggplot(ggplot2::aes(
        y=.data[[y_column]], fill=abs(.data[[fill]]), 
        x=prop, label=rank, xmin=lower, xmax=upper, text=hover_text
      )) +
      ggplot2::geom_col(col='black', linewidth=.6) +
      ggplot2::scale_fill_gradient(fill, low=low_col, high=high_col) +
      ggplot2::scale_x_continuous(labels=~sprintf("%.2f", abs(.x))) +
      ggplot2::scale_y_discrete(limits=rev(y_order)) +
      ggplot2::geom_errorbar(col="black", lwd=0.8, width=0.4) +
      ggplot2::facet_grid(~type, scales="free_x", axes='all_y') +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position='top', legend.direction='horizontal',
        legend.title=ggplot2::element_text(colour='black', face="bold"),
        legend.text=ggplot2::element_text(colour='black', face="bold"),
        legend.background=ggplot2::element_rect(fill='transparent', color=NA),
        legend.box.background=ggplot2::element_rect(fill='transparent', color=NA),
        axis.text.y=ggplot2::element_text(colour='black', face="bold"),
        axis.text.x=ggplot2::element_text(colour='black', face="bold"),
        axis.title=ggplot2::element_text(colour='black', face="bold"),
        panel.background=ggplot2::element_rect(fill='transparent', color=NA),
        panel.border=ggplot2::element_blank(),
        panel.spacing=ggplot2::unit(0, "lines"),
        strip.background=ggplot2::element_blank(),
        strip.text=ggplot2::element_text(colour='black', face="bold"),
        plot.background=ggplot2::element_rect(fill='transparent', color=NA)
    )
  )
} 

# Heatmap plotting function ----------------------------------------------------
heatmap_plot <- function(
    data, y_order=NULL, low_col='#FFC107', high_col='#1E88E5', 
    fill_column='prop.raw', max_x=Inf){
  
  stopifnot(  # Argument checks
    is.data.frame(data),
    is.null(y_order) | is.character(y_order),
    max_x == Inf | (is.numeric(max_x) & max_x > 1),
    'denominator' %in% names(data),  # Data is not stratified
    is.character(low_col), is.character(high_col),
    c(fill_column) %in% names(data)
  )
  
  denominator_column <- data[['denominator']][1]
  denominator_index <-  which(colnames(data) == denominator_column)
  
  subgroup_column <- names(data)[ifelse(denominator_index == 1, 2, 1)]
  
  y_column <- names(data)[1]
  data[[y_column]] <- as.character(data[[y_column]])  # Force character vector
  y_column_plural <- plural_name(y_column)
  if(is.null(y_order)) {y_order <- unique(data[[y_column]])}
  
  x_column <- names(data)[2]
  data[[x_column]] <- as.character(data[[x_column]])  # Force character vector
  x_column_plural <- plural_name(x_column)
  
  x_order <- data |>
    # dplyr::summarise(.by=tidyselect::all_of(c(x_column)), 
    #                  total_n=sum(count.raw, na.rm=TRUE)) |>
    dplyr::summarise(.by=tidyselect::all_of(c(x_column)), total_n=dplyr::n()) |>
    dplyr::arrange(-total_n) |>
    dplyr::pull(1)
  
  if(length(x_order) > max_x){  # Logic for "Other" column on heatmap X-axis
    x_order <- x_order[1:max_x]  #  Trim the order to max_x
    other_column <- paste("Other", x_column_plural)  # Get other column name
    data <- data |> 
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(c(x_column)),
          ~dplyr::if_else(.x %in% x_order, .x, other_column)
        )
      )
    x_order <- c(x_order, other_column)  # Add other column name to the order
  }
  return(
    data |>
      dplyr::summarise(
        .by=tidyselect::all_of(c(y_column, x_column)),
        dplyr::across(tidyselect::all_of(c(fill_column)), ~sum(.x, na.rm=TRUE))
      ) |> 
      dplyr::mutate(
        hover_text=paste0(
          .data[[subgroup_column]], " prevalence in ", .data[[denominator_column]],
          ": ", sprintf("%.3f", .data[[fill_column]])
        )
      ) |> 
      ggplot2::ggplot(ggplot2::aes(
        x=.data[[x_column]], y=.data[[y_column]], fill=.data[[fill_column]], 
        text=hover_text)
      ) +
      ggplot2::geom_tile(colour='black') +
      ggplot2::scale_fill_gradient(low=low_col, high=high_col) +
      ggplot2::scale_x_discrete(limits=x_order) +
      ggplot2::scale_y_discrete(limits=rev(y_order)) +
      ggplot2::labs(x=x_column_plural, y=y_column_plural) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position='top', legend.direction='horizontal',
        legend.text=ggplot2::element_text(colour='black', face="bold"),
        legend.title=ggplot2::element_text(colour='black', face="bold"),
        legend.background=ggplot2::element_rect(fill='transparent', color=NA),
        legend.box.background=ggplot2::element_rect(fill='transparent', color=NA),
        axis.text.y=ggplot2::element_text(colour='black', face="bold"),
        axis.text.x=ggplot2::element_text(colour='black', face="bold", angle=45, hjust=1),
        axis.title=ggplot2::element_text(colour='black', face="bold"),
        panel.background=ggplot2::element_rect(fill='transparent', color=NA),
        panel.border=ggplot2::element_blank(),
        panel.spacing=ggplot2::unit(0, "lines"),
        strip.background=ggplot2::element_blank(),
        strip.text=ggplot2::element_text(colour='black', face="bold"),
        plot.background=ggplot2::element_rect(fill='transparent', color=NA)
    )
  )
}

# Bar plotting function --------------------------------------------------------
bar_plot <- function(data, x_column, y_order=NULL, fill='black'){
  
  stopifnot(  # Argument checks
    is.data.frame(data),
    is.character(x_column) & x_column %in% names(data) & is.numeric(data[[x_column]]),
    is.null(y_order) | is.character(y_order),
    is.character(fill)
  )
  x_column_plural <- plural_name(x_column)
  
  y_column <- names(data)[1]
  data[[y_column]] <- as.character(data[[y_column]])  # Force character vector
  y_column_plural <- plural_name(y_column)
  if(is.null(y_order)) {y_order <- unique(data[[y_column]])}
  
  return(
    data |>
      dplyr::mutate(
        hover_text=paste0(
          nice_name(y_column), ": ", .data[[y_column]], "\n",
          x_column_plural, ": ", .data[[x_column]]
        )
      ) |> 
      ggplot2::ggplot(ggplot2::aes(x=.data[[x_column]], y=.data[[y_column]], 
                                   text=hover_text)) +
      ggplot2::geom_col(colour='black', fill=fill, linewidth=.5) +
      ggplot2::scale_y_discrete(limits=rev(y_order)) +
      ggplot2::labs(x=x_column_plural, y=y_column_plural) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y=ggplot2::element_text(colour='black', face="bold"),
        axis.text.x=ggplot2::element_text(colour='black', face="bold"),
        axis.title=ggplot2::element_text(colour='black', face="bold"),
        panel.background=ggplot2::element_rect(fill='transparent', color=NA),
        panel.border=ggplot2::element_blank(),
        panel.spacing=ggplot2::unit(0, "lines"),
        strip.background=ggplot2::element_blank(),
        strip.text=ggplot2::element_text(colour='black', face="bold"),
        plot.background=ggplot2::element_rect(fill='transparent', color=NA)
    ) 
  )
}

# Coverage plotting function ---------------------------------------------------
coverage_plot <-  function(data, x_order=NULL){
  
  stopifnot(  # Argument checks
    is.data.frame(data),
    c('prop.raw', 'upper.raw', 'lower.raw') %in% names(data),
    is.null(x_order) | is.character(x_order),
    'denominator' %in% names(data)
  )
  
  denominator_index <- which(colnames(data) == data[['denominator']][1])
  data[[denominator_index]] <- as.character(data[[denominator_index]])  # Force character vector
  denominator_column <- names(data)[denominator_index]
  denominator_column_plural <- plural_name(denominator_column)
  
  x_index <- ifelse(denominator_index==1, 2, 1)
  x_column <- names(data)[x_index]
  x_column_plural <- plural_name(x_column)
  data[[x_index]] <- as.character(data[[x_index]])  # Force character vector
  
  if(is.null(x_order)){  # Get the top x-axis variables
    x_order <- data |>
      dplyr::summarise(.by=x_index, total_n=sum(count.raw, na.rm=TRUE)) |>
      dplyr::arrange(-total_n) |>
      dplyr::pull(1)
  }
  
  return(
    data |>
      dplyr::select(
        denominator=denominator_index, x=x_index, prop=prop.raw, 
        upper=upper.raw, lower=lower.raw
      ) |> 
      tidyr::complete(denominator, x, fill=list(prop=0, upper=0, lower=0)) |> 
      dplyr::filter(x %in% x_order) |> 
      dplyr::mutate(x=factor(x, x_order)) |>
      dplyr::arrange(denominator, x) |> 
      dplyr::mutate(
        .by=denominator, 
        upper_diff=abs(prop-upper), lower_diff=abs(prop-lower), cum_prop=cumsum(prop), 
        cum_upper=cum_prop+upper_diff, cum_lower=cum_prop-lower_diff
      ) |> 
      dplyr::mutate(
        dplyr::across(
          c(cum_lower, cum_upper),
          ~dplyr::case_when(.x<0~0, .x>1~1, .default=.x)
        ),
        hover_text=paste0(
          x, " (rank = ", as.integer(x), ") coverage of ", denominator, ": ", 
          sprintf("%.3f", prop), "\n",
          "Current coverage of ", denominator, ": ", sprintf("%.3f", cum_prop), "\n",
          "Upper: ", sprintf("%.3f", abs(upper)), "\n",
          "Lower: ", sprintf("%.3f", abs(lower)), "\n",
          "Current upper: ", sprintf("%.3f", abs(cum_upper)), "\n",
          "Current lower: ", sprintf("%.3f", abs(cum_lower))
        )
      ) |> 
      ggplot2::ggplot(ggplot2::aes(
        x=x, y=cum_prop, group=denominator, fill=denominator, ymin=cum_lower, 
        ymax=cum_upper, text=hover_text
      )) +
      ggplot2::geom_point(
        alpha=.8, size=3, colour='black', shape=21, show.legend=TRUE) +
      ggplot2::geom_ribbon(alpha=0.1, show.legend=FALSE) +
      ggplot2::labs(fill=denominator_column_plural, x=x_column_plural, 
                    y='Cumulative infection coverage') +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.text=ggplot2::element_text(colour='black', face="bold"),
        legend.title=ggplot2::element_text(colour='black', face="bold"),
        legend.background=ggplot2::element_rect(fill='transparent', color=NA),
        legend.box.background=ggplot2::element_rect(fill='transparent', color=NA),
        axis.text.y=ggplot2::element_text(colour='black', face="bold"),
        axis.text.x=ggplot2::element_text(colour='black', face="bold", angle=45, hjust=1),
        axis.title=ggplot2::element_text(colour='black', face="bold"),
        panel.background=ggplot2::element_rect(fill='transparent', color=NA),
        panel.border=ggplot2::element_blank(),
        strip.background=ggplot2::element_blank(),
        strip.text=ggplot2::element_text(colour='black', face="bold"),
        plot.background=ggplot2::element_rect(fill='transparent', color=NA)
    ) 
  )
}

# Helper functions -------------------------------------------------------------
nice_name <- function(x, replace_chars='_.') {
  replace_regex <- paste0("[", replace_chars, "]")
  return(stringr::str_to_title(stringr::str_replace_all(x, replace_regex, ' ')))
}

plural_name <- function(x, .default='s', nice=FALSE) {
  if(isTRUE(nice)){x <- nice_name(x)}
  return(
    dplyr::case_when(
      .default=paste0(x, .default),
      endsWith(x, 'us') ~ paste0(stringr::str_sub(x, end=-3), 'i'),
      endsWith(x, 's') ~ x,
      endsWith(x, 'y') ~ paste0(stringr::str_sub(x, end=-2), 'ies'),
    ) 
  )
}
