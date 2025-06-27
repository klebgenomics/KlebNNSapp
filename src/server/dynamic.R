library(shiny)
library(dplyr)
library(stringr)
library(plotly)
library(shinyWidgets)
library(rclipboard)
library(leaflet)
library(htmltools)
library(shinyjqui)

source('src/R/prevalence.R')
source('src/R/plotting.R')

# Reactive data ----------------------------------------------------------------
reactive_dynamic_data <- shiny::reactive({
  d <- DATA |>
    dplyr::filter(
      dplyr::between(Year, input$dynamic_year_selector[1], input$dynamic_year_selector[2]),
      Study %in% input$dynamic_study_selector,
      Region %in% input$dynamic_region_selector,
      resistance_score >= dplyr::case_match(
        input$dynamic_amr_selector, .default=0, 'ESBL+' ~ 1, 'Carbapenemase+' ~ 2
      ),
      if(input$dynamic_outcome_selector != 'All samples') Mortality == 1 else TRUE,
    )
  if(nrow(d) == 0){
    shiny::showNotification('No strains matching filter', type='warning')
  }
  return(d)
})

# Data selector reset ----------------------------------------------------------
shiny::observeEvent(input$dynamic_data_reset, {
  shiny::showNotification('Resetting data filters', type='message')
  shinyWidgets::updatePickerInput(session, 'dynamic_amr_selector', selected='All samples')
  shiny::updateSliderInput(session, 'dynamic_year_selector', value=c(MIN_YEARS, MAX_YEARS))
  shinyWidgets::updatePickerInput(session, 'dynamic_region_selector', selected=REGIONS)
  shinyWidgets::updatePickerInput(session, 'dynamic_study_selector', selected=STUDIES)
  shinyWidgets::updatePickerInput(session, 'dynamic_outcome_selector', selected='All samples')
})

# Reactive prevalences ---------------------------------------------------------
reactive_dynamic_prevalence <- shiny::reactive({
  return(
    prevalence(
      data=reactive_dynamic_data(),
      stratify_by=input$dynamic_antigen_selector,
      adjust_for=c('Cluster'),
      n_distinct=setdiff(X_AXIS_VARIABLES, input$dynamic_antigen_selector)
    )
  )
})

reactive_dynamic_prevalence_stratified <- shiny::reactive({
  return(
    prevalence(
      data=reactive_dynamic_data(),
      stratify_by=c(input$dynamic_antigen_selector, input$dynamic_heatmap_x),
      adjust_for=c('Cluster'),
      denominator=ifelse(
        input$dynamic_heatmap_swap_denominator,
        input$dynamic_antigen_selector, input$dynamic_heatmap_x
      )
    )
  )
})

reactive_dynamic_prevalence_coverage <- shiny::reactive({
  return(
    prevalence(
      data=reactive_dynamic_data(),
      stratify_by=c(input$dynamic_coverage_group, input$dynamic_antigen_selector),
      adjust_for=c('Cluster')
    ) 
  )
})

# Reactive antigen lists -------------------------------------------------------
shiny::observeEvent(input$dynamic_antigen_selector, {
  new_vars <- setdiff(X_AXIS_VARIABLES, input$dynamic_antigen_selector)
  shinyWidgets::updatePickerInput(
    session=session, 
    inputId="dynamic_heatmap_x", 
    choices=new_vars
  )
  shinyWidgets::updatePickerInput(
    session=session, 
    inputId='dynamic_bars_x', 
    choices=new_vars
  )
})

shiny::observeEvent(list(input$dynamic_antigen_reset, 
                         input$dynamic_antigen_selector,
                         input$dynamic_vaccine_valency), {
  # Reset the text list
  shiny::updateTextInput(
     session=session,
     inputId='dynamic_text_antigens', 
     value=''
  )
                           
  # Reset the click and drag lists                   
  all_antigens <- unique(reactive_dynamic_prevalence()[[1]])  # First column
  # Unique just in case, but should be unique anyway
  valency <- min(input$dynamic_vaccine_valency, length(all_antigens))
  default_antigens <- all_antigens[1:valency]
  
  shinyjqui::updateOrderInput(
    session=session,
    inputId='dynamic_sorted_antigens',
    items=default_antigens
  )
  
  shinyjqui::updateOrderInput(
    session=session,
    inputId='dynamic_available_antigens',
    items=setdiff(all_antigens, default_antigens)
  )
})

dynamic_antigen_list <- shiny::reactive({
  user_list <- unique(base::unlist(stringr::str_extract_all(
    input$dynamic_text_antigens, LOCUS_TYPE_REGEX)))
  
  if(length(user_list) > 0){
    valency <- min(input$dynamic_vaccine_valency, length(user_list))
    return(user_list[1:valency])
  } else {
    valency <- min(input$dynamic_vaccine_valency, length(input$dynamic_sorted_antigens))
    return(input$dynamic_sorted_antigens[1:valency])
  }
})

# Antigen clipboard copy -------------------------------------------------------
output$dynamic_antigen_copy <- shiny::renderUI({
  rclipboard::rclipButton(
    "dynamic_antigen_copy", shiny::icon("clipboard"),
    dynamic_antigen_list(),
    tooltip='Copy antigen order to clipboard'
  )
})

output$dynamic_summary <- shiny::renderText({
  return(
    paste0(
      "Samples: ", nrow(reactive_dynamic_data()), '/', nrow(DATA), "; ",
      "Studies: ", length(input$dynamic_study_selector), '/', length(STUDIES), "; ",
      "Regions: ", length(input$dynamic_region_selector), '/', length(REGIONS), "; ",
      "Years: ", input$dynamic_year_selector[2] - input$dynamic_year_selector[1],
      '/', MAX_YEARS-MIN_YEARS, "; ",
      "Resistance: ", input$dynamic_amr_selector, "; ",
      "Outcome: ", input$dynamic_outcome_selector
    )
  )
})

# Output merged plot -----------------------------------------------------------
output$dynamic_global <- plotly::renderPlotly({
  prevalence <- reactive_dynamic_prevalence()
  if(is.null(prevalence) || nrow(prevalence) == 0){return(NULL)}
  prevalence_stratified <- reactive_dynamic_prevalence_stratified()
  if(is.null(prevalence_stratified) || nrow(prevalence_stratified) == 0){return(NULL)}
  y_order <- dynamic_antigen_list()
  if(is.null(y_order) | length(y_order) == 0){return(NULL)}
  
  pyramid_gg <- pyramid_plot(
      data=prevalence,
      low_col=input$dynamic_pyramid_low_col,
      high_col=input$dynamic_pyramid_hi_col,
      y_order=y_order
    ) |> 
    plotly::ggplotly(tooltip = "text")

  heatmap_gg <- heatmap_plot(
      data=prevalence_stratified,
      low_col=input$dynamic_heatmap_low_col,
      high_col=input$dynamic_heatmap_hi_col,
      y_order=y_order,
      max_x=input$dynamic_heatmap_num_x
    ) |> 
    plotly::ggplotly(tooltip = "text")

  bars_gg <- bar_plot(
      data=prevalence,
      x_column=paste('#', input$dynamic_bars_x),
      y_order=y_order,
      fill=input$dynamic_bars_fill
    ) |> 
    plotly::ggplotly(tooltip = "text")

  return(
    plotly::subplot(
      pyramid_gg, heatmap_gg, bars_gg, widths=c(.4, .4, .2), margin=0.01, 
      shareY=TRUE, nrows=1, titleX=TRUE, titleY=TRUE
    ) |> 
      plotly::hide_guides()
  )
})

# Output cumulative plot -------------------------------------------------------
output$dynamic_coverage_plot <- plotly::renderPlotly({
  prevalence_coverage <- reactive_dynamic_prevalence_coverage()
  if(is.null(prevalence_coverage) || nrow(prevalence_coverage) == 0){return(NULL)}
  return(
    plotly::ggplotly(
      coverage_plot(data=prevalence_coverage, x_order=dynamic_antigen_list()),
      tooltip = "text"
    )
  )
})

# Plot map ---------------------------------------------------------------------
output$dynamic_map <- leaflet::renderLeaflet({
  prevalence_coverage <- reactive_dynamic_prevalence_coverage()
  if(is.null(prevalence_coverage) || nrow(prevalence_coverage) == 0){
    return(
      leaflet::leaflet() |>
        leaflet::addProviderTiles(leaflet::providers$CartoDB.Voyager)
    )
  }
  map_data <- WORLD_SF |>
    dplyr::inner_join(
      prevalence_coverage |>
        dplyr::filter(dplyr::if_all(2, ~.x %in% dynamic_antigen_list())) |> 
        dplyr::summarise(
          .by=1, 
          Infections=sum(count.raw), 
          Coverage=min(sum(prop.raw), 1)
        ),
      by=input$dynamic_coverage_group
    )
  
  map_text <- paste0(
    input$dynamic_coverage_group, ": ", 
    map_data[[input$dynamic_coverage_group]], "<br/>",
    "Infections: ", map_data$Infections, "<br/>",
    "Coverage: ", sprintf("%.3f", map_data$Coverage)
    ) |>
    lapply(htmltools::HTML)
  
  map_data |> 
    leaflet::leaflet() |>
    leaflet::addProviderTiles(leaflet::providers$CartoDB.Voyager) |>
    leaflet::addPolygons(
      fillColor=~MAP_PALETTE(Coverage),
      color="white", weight=2, opacity=.8,
      label=map_text,
      labelOptions=leaflet::labelOptions(
        style=list("font-weight"="normal", padding="3px 8px"),
        textsize="13px", direction="auto"
      ),
      highlightOptions=leaflet::highlightOptions(
        weight=3, color="black", fillOpacity=.8, bringToFront=TRUE
      )
    ) |> 
    leaflet::addLegend(
      pal=MAP_PALETTE, values=~Coverage, opacity=0.8,
      title="Infection Coverage", position="bottomleft"
    )
})
