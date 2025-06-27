library(shiny)
library(dplyr)
library(stringr)
library(glue)
library(plotly)
library(shinyWidgets)
library(rclipboard)
library(leaflet)
library(htmltools)

source('src/R/prevalence.R')
source('src/R/plotting.R')

# Reactive prevalences ---------------------------------------------------------
reactive_modelled_prevalence <- shiny::reactive({
  MODELLED |> 
    dplyr::filter(
      Antigen == input$modelled_antigen_selector,
      Infection == input$modelled_infection_selector
    ) |> 
    dplyr::select(!c(Antigen, Infection))
})

reactive_modelled_prevalence_stratified <- shiny::reactive({
  MODELLED_STRATIFIED |> 
    dplyr::filter(
      Antigen == input$modelled_antigen_selector,
      Infection == input$modelled_infection_selector
    ) |> 
    dplyr::select(!c(Antigen, Infection))
})

# Reactive loci lists ----------------------------------------------------------
shiny::observeEvent(list(input$modelled_antigen_reset, 
                         input$modelled_antigen_selector,
                         input$modelled_vaccine_valency), {
   # Reset the text list
   shiny::updateTextInput(
     session=session,
     inputId='modelled_text_antigens', 
     value=''
   )
   
   # Reset the click and drag lists                   
   all_antigens <- unique(reactive_modelled_prevalence()[[1]])  # First column
   # Unique just in case, but should be unique anyway
   valency <- min(input$modelled_vaccine_valency, length(all_antigens))
   default_antigens <- all_antigens[1:valency]
   
   shinyjqui::updateOrderInput(
     session=session,
     inputId='modelled_sorted_antigens',
     items=default_antigens
   )
   
   shinyjqui::updateOrderInput(
     session=session,
     inputId='modelled_available_antigens',
     items=setdiff(all_antigens, default_antigens)
   )
 })

modelled_antigen_list <- shiny::reactive({
  user_list <- unique(base::unlist(stringr::str_extract_all(
    input$modelled_text_antigens, LOCUS_TYPE_REGEX)))
  
  if(length(user_list) > 0){
    valency <- min(input$modelled_vaccine_valency, length(user_list))
    return(user_list[1:valency])
  } else {
    valency <- min(input$modelled_vaccine_valency, length(input$modelled_sorted_antigens))
    return(input$modelled_sorted_antigens[1:valency])
  }
})

# Antigen clipboard copy -------------------------------------------------------
output$modelled_antigen_copy <- shiny::renderUI({
  rclipboard::rclipButton(
    "modelled_antigen_copy", shiny::icon("clipboard"),
    modelled_antigen_list(),
    tooltip='Copy antigen order to clipboard'
  )
})

# Output merged plot -----------------------------------------------------------
output$modelled_merged_plot <- plotly::renderPlotly({
  prevalence <- reactive_modelled_prevalence()
  if(is.null(prevalence) || nrow(prevalence) == 0){return(NULL)}
  prevalence_stratified <- reactive_modelled_prevalence_stratified()
  if(is.null(prevalence_stratified) || nrow(prevalence_stratified) == 0){return(NULL)}
  y_order <- modelled_antigen_list()
  if(is.null(y_order) | length(y_order) == 0){return(NULL)}
  
  pyramid_gg <- pyramid_plot(
    data=prevalence,
    low_col=input$modelled_pyramid_low_col,
    high_col=input$modelled_pyramid_hi_col,
    y_order=y_order,
    fill='prop'
  ) |> 
    plotly::ggplotly(tooltip = "text")
  
  heatmap_gg <- heatmap_plot(
    data=dplyr::select(prevalence_stratified, 2, 1, dplyr::everything()),
    low_col=input$modelled_heatmap_low_col,
    high_col=input$modelled_heatmap_hi_col,
    y_order=y_order
  ) |> 
    plotly::ggplotly(tooltip = "text")
  
  return(
    plotly::subplot(
      pyramid_gg, heatmap_gg, widths=c(.5, .5), margin=0.01, 
      shareY=TRUE, nrows=1, titleX=TRUE, titleY=TRUE
    ) |> 
      plotly::hide_guides()
  )
})

# Output cumulative plot -------------------------------------------------------
output$modelled_coverage_plot <- plotly::renderPlotly({
  prevalence_coverage <- reactive_modelled_prevalence_stratified()
  if(is.null(prevalence_coverage) || nrow(prevalence_coverage) == 0){return(NULL)}
  return(
    plotly::ggplotly(
      coverage_plot(data=prevalence_coverage, x_order=modelled_antigen_list()),
      tooltip='text'
    )
  )
})

# Plot map ---------------------------------------------------------------------
output$modelled_map <- leaflet::renderLeaflet({
  map_data <- WORLD_SF |>
    dplyr::inner_join(
      reactive_modelled_prevalence_stratified() |>
        dplyr::filter(dplyr::if_all(2, ~.x %in% modelled_antigen_list())) |> 
        dplyr::summarise(
          .by=1, 
          Coverage=min(sum(prop.raw), 1)
        ),
      by='Region'
    )
  
  map_text <- paste0(
    "Region: ", map_data[['Region']], "<br/>", 
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

# output$modelled_summary <- shiny::renderText({
#   return(
#     paste0(
#       "Samples: ", nrow(reactive_modelled_prevalence()), '/', nrow(MODELLED), "; ",
#       "Infection: ", input$modelled_infection_selector
#     )
#   )
# })
