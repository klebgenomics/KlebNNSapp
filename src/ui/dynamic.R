library(shiny)
library(bslib)
library(shinyWidgets)
library(plotly)
library(leaflet)
library(colourpicker)
library(htmltools)
library(shinyjqui)

# Antigen selector -------------------------------------------------------------
dynamic_antigen_selector <- bslib::accordion_panel(
  title=htmltools::h5('Antigen selector'), value='dynamic_antigen_selector',
  icon=shiny::icon('syringe'),
  shiny::fluidRow(
    shiny::column(
      6,
      shinyWidgets::pickerInput(
        'dynamic_antigen_selector', htmltools::h6("Select antigen:"), 
        choices=c('K_locus', 'O_type'), 
        selected='K_locus'
      ),
    ),
    shiny::column(
      6,
      shiny::sliderInput(
        'dynamic_vaccine_valency', htmltools::h6("# of antigens"), 
        min=1, max=length(KL), value=DEFAULT_VALENCY
      )
    )
  ),
  htmltools::em(
  'Hit update to order antigens by the most prevalent in the current data'
  ),
  htmltools::br(),
  shiny::fluidRow(
    shiny::column(
      10,
      shiny::actionButton(
        'dynamic_antigen_reset', "Update antigen order", class="btn-primary", 
        icon=shiny::icon("spinner"), width='100%'
      ) 
    ),
    shiny::column(1,  shiny::uiOutput("dynamic_antigen_copy"))
  ),
  htmltools::br(),
  htmltools::p('Type or paste in the box below to order antigens manually; 
      this will overide the click-and-drag order unless the box is empty:'),
  shiny::textInput('dynamic_text_antigens', label=NULL),
  shiny::fluidRow(
    htmltools::h6('Drag to order antigens'),
    shiny::column(
      width = 6,
      shiny::wellPanel(
        shinyjqui::orderInput(
          width='100px',
          label="Selected",
          items=list(),
          inputId="dynamic_sorted_antigens",
          connect='dynamic_available_antigens',
          item_class = "warning"
        )
      )
    ),
    shiny::column(
      width = 6,
      shiny::wellPanel(
        shinyjqui::orderInput(
          width='100px',
          label="Available",
          items=list(),
          inputId="dynamic_available_antigens",
          connect='dynamic_sorted_antigens'
        )
      )
    )
  )
)

# Data filters -----------------------------------------------------------------
dynamic_data_selector <- bslib::accordion_panel(
  title=htmltools::h5('Data selector'), value='dynamic_data_selector',
  icon=shiny::icon('filter'),
  htmltools::h6("Select which samples to include"),
  htmltools::em("Remember to 'update antigen order' in the antigen selector to 
                reset plot axis orders to the current data!"),
  htmltools::br(),
  shinyWidgets::pickerInput(
    'dynamic_amr_selector', "Resistance:",
    choices=AMR, 
    selected='All samples'
  ),
  shinyWidgets::pickerInput(
    'dynamic_region_selector', "Regions:", multiple=TRUE,
    choices=REGIONS, 
    selected=REGIONS,
    options=shinyWidgets::pickerOptions(actionsBox=TRUE, liveSearch=TRUE)
  ),
  shinyWidgets::pickerInput(
    'dynamic_study_selector', "Studies:", multiple=TRUE,
    choices=STUDIES, 
    selected=STUDIES, 
    options=shinyWidgets::pickerOptions(actionsBox=TRUE, liveSearch=TRUE)
  ),
  shinyWidgets::pickerInput(
    'dynamic_outcome_selector', "Outcome:",
    choices=c('All samples', 'Fatal cases'), 
    selected='All samples'
  ),
  shiny::sliderInput(
    'dynamic_year_selector', "Sampling years:",
    MIN_YEARS, MAX_YEARS, c(MIN_YEARS, MAX_YEARS), sep=""
  ),
  shiny::actionButton(
    'dynamic_data_reset', "Reset data selector", class="btn-primary", 
    icon=shiny::icon("arrows-rotate"), width='100%'
  )
)

# Sidebar ----------------------------------------------------------------------
dynamic_sidebar <- bslib::sidebar(
  width=SIDEBAR_WIDTH, bg='white',
  bslib::accordion(
    multiple=FALSE, open='dynamic_data_selector',
    dynamic_antigen_selector, dynamic_data_selector
  )
)

# Global panel -----------------------------------------------------------------
dynamic_prevalence_panel <- bslib::accordion_panel(
  title=htmltools::h4('Global prevalence'), value='dynamic_prevalence_panel',
  icon=shiny::icon('earth-africa'),
  shiny::fluidRow(
    shiny::column(
      width=3,
      DROPDOWN_FUNCTION(
        shiny::fluidRow(
          shiny::column(
            width=6, 
            colourpicker::colourInput('dynamic_pyramid_low_col', "Low", LO_FILL)
          ),
          shiny::column(
            width=6, 
            colourpicker::colourInput('dynamic_pyramid_hi_col', "High", HI_FILL)
          )
        )
      ),
      htmltools::h5('Serotype prevalence')
    ),
    shiny::column(width=2),
    shiny::column(
      width=3,
      DROPDOWN_FUNCTION(
        shinyWidgets::pickerInput(
          'dynamic_heatmap_x', "Subgroup", X_AXIS_VARIABLES, 'ST'
        ),
        shiny::sliderInput(
          'dynamic_heatmap_num_x', "Maximum X-axis values", 
          min=4, max=50, value=DEFAULT_VALENCY
        ),
        shinyWidgets::switchInput(
          'dynamic_heatmap_swap_denominator', size='large', inline=FALSE, 
          label='Denominator', onLabel="Serotype", offLabel="Subgroup",
        ),
        shiny::fluidRow(
          shiny::column(
            width=6, 
            colourpicker::colourInput('dynamic_heatmap_low_col', "Low", LO_FILL)
          ),
          shiny::column(
            width=6, 
            colourpicker::colourInput('dynamic_heatmap_hi_col', "High", HI_FILL)
          )
        )
      ),
      htmltools::h5('Prevalence by subgroup')
    ),
    shiny::column(width=2),
    shiny::column(
      width=2,
      DROPDOWN_FUNCTION(
        shinyWidgets::pickerInput(
          'dynamic_bars_x', "Subgroup for barplot", X_AXIS_VARIABLES, 'ST'
        ),
        colourpicker::colourInput('dynamic_bars_fill', "Bar colour", BAR_FILL)
      ),
      htmltools::h5('Unique subgroup values per serotype')
    )
  ),
  plotly::plotlyOutput("dynamic_global", height="600px") |>
    shinycssloaders::withSpinner(color=SPINNER_COLOR, type=SPINNER_TYPE, size=SPINNER_SIZE)
  # bslib::card(
  #   full_screen=TRUE,
  #   bslib::card_body(
  #     class="p-0",
  #    
  #   )
  # )
)

# Coverage panel ---------------------------------------------------------------
dynamic_coverage_panel <- bslib::accordion_panel(
  title=htmltools::h4('Geographical coverage'), value='dynamic_coverage_panel',
  icon=shiny::icon('map'),
  DROPDOWN_FUNCTION(
    shinyWidgets::pickerInput(
      'dynamic_coverage_group', "Calculate coverage of", 
      choices=c('Country', 'Region'),
      selected='Country'
    ),
  ),
  bslib::layout_column_wrap(
    width=1/2,
    height=400,
    bslib::card(
      full_screen=TRUE,
      bslib::card_body(
        class="p-0",
        plotly::plotlyOutput("dynamic_coverage_plot") |> 
          shinycssloaders::withSpinner(color=SPINNER_COLOR, type=SPINNER_TYPE, size=SPINNER_SIZE)
      )
    ),
    bslib::card(
      full_screen=TRUE,
      bslib::card_body(
        class="p-0",
        leaflet::leafletOutput("dynamic_map") |> 
          shinycssloaders::withSpinner(color=SPINNER_COLOR, type=SPINNER_TYPE, size=SPINNER_SIZE)
      )
    )
  )
)

# Main panel -------------------------------------------------------------------
dynamic <- bslib::nav_panel(
  title='Dynamic antigen prevalence', icon=shiny::icon('dna'),
  bslib::layout_sidebar(
    sidebar=dynamic_sidebar,
    htmltools::p(
      "
      These plots show the raw and adjusted prevalences of the chosen antigen 
      amongst the infections passing the filter in the 'Data selector' in the sidebar. 
      The default antigen order is based on adjusted prevalence, however this can be 
      changed using the 'Antigen selector' in the sidebar. Adjustments are performed 
      by only counting one isolate per cluster, thereby reducing the impact of localised
      outbreaks on prevalence estimates.
      "
    ),
    shiny::textOutput('dynamic_summary'),
    bslib::accordion(
      multiple=TRUE, open='dynamic_prevalence_panel', 
      dynamic_prevalence_panel, dynamic_coverage_panel
    )
  )
)
