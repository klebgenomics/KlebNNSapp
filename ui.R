library(shiny)
library(bslib)
library(rclipboard)
library(shinyWidgets)
library(htmltools)
# Load in the separate UI files
source('src/ui/dynamic.R')
source('src/ui/modelled.R')

# Convenience function for adding images with clickable hyperlinks
create_logo_link <- function(id, src, url, width="100%", tooltip_text=NULL) {
  img_tag <- htmltools::img(
    src=src, width=width, style="vertical-align: middle;"
  )  # Ensure vertical alignment
  link_tag <- shiny::actionLink(inputId=id, label=img_tag, onclick=sprintf("window.open('%s')", url))
  if (!is.null(tooltip_text)) {
    bslib::tooltip(link_tag, tooltip_text, placement='bottom')
  } else {
    link_tag
  }
}

# Define all hyperlinked images here -------------------------------------------
kaptive_logo <- create_logo_link(
  id="kaptive_logo",
  src="kaptive.png",
  url="https://kaptive.readthedocs.io",
  width="100px",
  tooltip_text='Read the docs'
)
kleborate_logo <- create_logo_link(
  id="kleborate_logo",
  src="kleborate.png",
  url="https://kleborate.readthedocs.io",
  width="100px",
  tooltip_text='Read the docs'
)
monash_logo <- create_logo_link(
  id='monash_logo',
  src="monash.svg",
  width="100px",
  url="https://www.monash.edu"
)
lshtm_logo <- create_logo_link(
  id='lshtm_logo',
  src="lshtm.png",
  width="70px",
  url="https://www.lshtm.ac.uk"
)

# Define footer ----------------------------------------------------------------
footer <- htmltools::div(
  style="footer; margin: 10px auto; width: 100%;", # Ensure it takes the full width
  # htmltools::hr(),
  htmltools::div(
    style="display: flex; justify-content: space-between; align-items: center;",
    htmltools::div(
      style="display: flex; align-items: center; flex-wrap: nowrap;",
      lshtm_logo, monash_logo
    ),
    htmltools::div(
      style="display: flex; align-items: center; flex-wrap: nowrap;",
      kaptive_logo, kleborate_logo,
    )
  )
)

# Define home panel ------------------------------------------------------------
home <- bslib::nav_panel(
  title='Home', icon=shiny::icon('home'),
  htmltools::tags$head(
    rclipboard::rclipboardSetup(),
    htmltools::tags$style(  # This is needed for the click-drag antigen order UI
        htmltools::HTML("
        .bttn-jelly.bttn-danger {
          background: #543391;
        }
        .ui-sortable-handle {
          border-radius: 3px;
          display: block;
          padding: 1px 3px;
          background-color: #f8f8f8;
          border: 1px solid #ddd;
          overflow: hidden;
          width: 100%;
        }
      ")
    )
  ),
  shiny::fluidRow(
    style="margin-top: 20px;",
    shiny::column(
      width=10, offset=1,
      bslib::card(
        bslib::card_body(
          htmltools::HTML(
            "
            <p>This app allows users to explore the distribution of predicted K 
            and O serotypes for <i>Klebsiella pneumoniae</i> isolated from neonatal 
            sepsis cases in 13 studies across countries in Africa and Southern 
            Asia, reported in the paper Stanton et al, 2025.</p>
            
            <p>The functionality is geared towards exploring sets of K/O antigens, 
            in terms of their prevalence and distribution across geographical 
            regions and theoretical coverage of infection isolates, to inform 
            vaccine design.</p> 
            
            <p>Prevalence estimates are adjusted for localised nosocomial 
            clustering, to reduce the bias introduced by random outbreaks during 
            surveillance periods. Coverage estimates are based on total isolates, 
            not adjusted for clustering.</p>
            
            <p>The <b>Modelled antigen prevalence</b> tab is populated with 
            pre-calculated global and regional prevalence estimates modelled 
            using Bayesian hierarchical meta-analysis, as described in the 
            paper. Subgroup analyses are limited to those modelled and reported 
            in the paper (geographic regions, fatal cases, ESBL- or carbapenemase- carrying isolates).
            </p>
            
            <p>The <b>Dynamic antigen prevalence</b> tab is populated with simple pooled 
            estimates calculated on the fly, allowing users to interactively explore prevalence 
            and coverage more flexibly by country, study, year, and multi-locus sequence type 
            (ST).</p>
            
            <b>Exploring plots</b>
            
            <p>All plots are interactive and can be downloaded via Plotly,
            so play around, zoom in and hover over the plots to get a 
            closer look at specific variables. The plot colours and
            X-axis variables/size can also be adjusted, just click on the
            cog buttons above each plot!</p>
            
            <b>Citations and feedback</b>
            
            <p>If you use this in your work, please cite the paper Stanton <em>et al</em>, 2025.
            The R shiny code is available in <a href='https://github.com/klebgenomics/KlebNNSapp/'>GitHub</a>,
            so you can also download and run the app locally offline.
            If you encounter issues please post them <a href='https://github.com/klebgenomics/KlebNNSapp/issues'>here</a>.
            </p>
            "
          )
        )
      )
    )
  ),
  htmltools::hr()
) 

# Define the UI ----------------------------------------------------------------
bslib::page_navbar(
  theme=bslib::bs_theme(version=5, bootswatch='pulse'),
  window_title="KlebSeroEpi", fillable=FALSE, title='KlebSeroEpi',
  # Main body ------------------------------------------------------------------
  bslib::navset_tab(home, dynamic, modelled),
  # Footer ---------------------------------------------------------------------
  footer=footer
)
