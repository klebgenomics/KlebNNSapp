library(shiny)
library(bslib)
library(shinyWidgets)
library(plotly)
library(leaflet)
library(colourpicker)
library(htmltools)
library(shinyjqui)


# UI ---------------------------------------------------------------------------
modelled_sidebar <- bslib::sidebar(
  width=SIDEBAR_WIDTH, bg='white', 
  bslib::accordion(
    multiple=FALSE, open='modelled_data_selector',
    bslib::accordion_panel(
      title=htmltools::h5('Antigen selector'), value='modelled_antigen_selector',
      icon=shiny::icon('syringe'),
      shiny::fluidRow(
        shiny::column(
          width=6,
          shinyWidgets::pickerInput(
            'modelled_antigen_selector', htmltools::h6("Select antigen:"), 
            choices=c('K_locus', 'O_type'),
            selected='K_locus'
          )
        ), 
        shiny::column(
          width=6,
          shiny::sliderInput(
            'modelled_vaccine_valency', htmltools::h6("# of antigens"), 
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
          width=10,
          shiny::actionButton(
            'modelled_antigen_reset', "Update antigen order", class="btn-primary", 
            icon=shiny::icon("spinner"), width='100%'
          ) 
        ),
        shiny::column(1,  shiny::uiOutput("modelled_antigen_copy"))
      ),
      htmltools::br(),
      htmltools::p('Type or paste in the box below to order antigens manually; 
      this will overide the click-and-drag order unless the box is empty:'),
      shiny::textInput('modelled_text_antigens', label=NULL),
      shiny::fluidRow(
        htmltools::h6('Drag to order antigens'),
        shiny::column(
          width = 6,
          shiny::wellPanel(
            shinyjqui::orderInput(
              width='100px',
              label="Selected",
              items=list(),
              inputId="modelled_sorted_antigens",
              connect='modelled_available_antigens',
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
              inputId="modelled_available_antigens",
              connect='modelled_sorted_antigens'
            )
          )
        )
      )
    ),
    bslib::accordion_panel(
      title=htmltools::h5('Data selector'), value='modelled_data_selector',
      icon=shiny::icon('filter'),
      htmltools::h6("Select which samples to include"),
      htmltools::em("Remember to 'update antigen order' in the antigen selector 
                    to reset plot axis orders to the current data!"),
      htmltools::br(),
      shinyWidgets::pickerInput(
        'modelled_infection_selector', "Select infection type:", 
        choices=MODELLED_INFECTION_TYPES, 
        selected='All samples'
      ),
      shiny::actionButton(
        'modelled_data_reset', "Reset data selector", class="btn-primary", 
        icon=shiny::icon("arrows-rotate"), width='100%'
      )
    )
  )
)

modelled_global_panel <- bslib::accordion_panel(
  title=htmltools::h4('Global prevalence'), value='modelled_global_panel',
  icon=shiny::icon('earth-africa'),
  shiny::fluidRow(
    shiny::column(
      width=6,
      DROPDOWN_FUNCTION(
        shiny::fluidRow(
          shiny::column(
            width=6, 
            colourpicker::colourInput('modelled_pyramid_low_col', "Low", LO_FILL)
          ),
          shiny::column(
            width=6, 
            colourpicker::colourInput('modelled_pyramid_hi_col', "High", HI_FILL)
          )
        )
      ),
      htmltools::h5('Serotype prevalence')
    ),
    shiny::column(
      width=6, 
      DROPDOWN_FUNCTION(
        shiny::fluidRow(
          shiny::column(
            width=6, 
            colourpicker::colourInput('modelled_heatmap_low_col', "Low", LO_FILL)
          ),
          shiny::column(
            width=6, 
            colourpicker::colourInput('modelled_heatmap_hi_col', "High", HI_FILL)
          )
        )
      ),
      htmltools::h5('Prevalence by subgroup')
    )
  ),
  plotly::plotlyOutput("modelled_merged_plot", height="600px") |> 
    shinycssloaders::withSpinner(color=SPINNER_COLOR, type=SPINNER_TYPE, size=SPINNER_SIZE),
)

modelled_coverage_panel <- bslib::accordion_panel(
  title=htmltools::h4('Geographical coverage'), value='modelled_coverage_panel',
  icon=shiny::icon('map'),
  bslib::layout_column_wrap(
    width=1/2,
    height=400,
    bslib::card(
      full_screen=TRUE,
      bslib::card_body(
        class="p-0",
        plotly::plotlyOutput("modelled_coverage_plot") |> 
          shinycssloaders::withSpinner(color=SPINNER_COLOR, type=SPINNER_TYPE, size=SPINNER_SIZE)
      )
    ),
    bslib::card(
      full_screen=TRUE,
      bslib::card_body(
        class="p-0",
        leaflet::leafletOutput("modelled_map") |> 
          shinycssloaders::withSpinner(color=SPINNER_COLOR, type=SPINNER_TYPE, size=SPINNER_SIZE)
      )
    )
  )
)
 
modelled <- bslib::nav_panel(
  title='Modelled antigen prevalence', icon=shiny::icon('calculator'),
  bslib::layout_sidebar(
    sidebar=modelled_sidebar,
    htmltools::HTML(
      "<p>These plots allow exploration of pre-calculated global and regional 
      prevalence estimates, modelled using Bayesian hierarchical meta-analysis. 
      Subgroup analyses are limited to those modelled and reported in the paper 
      (Stanton <em>et al</em>, 2025), i.e. ESBL+, carbapenemase+ or fatal cases.
      </p>"
    ),
    # shiny::textOutput('modelled_summary'),
    bslib::accordion(
      multiple=TRUE, open='modelled_global_panel', 
      modelled_global_panel, modelled_coverage_panel
    )
  )
)
