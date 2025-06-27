library(dplyr)
library(vroom)
library(tidyselect)

source('src/R/utils.R', local=TRUE)

kleborate_column_spec <- vroom::cols(
  `Genome ID` = vroom::col_character(),
  `Genome Name` = vroom::col_character(),
  `Kleborate version` = vroom::col_character(),
  `Wrapper version` = vroom::col_double(),
  species = vroom::col_character(),
  species_match = vroom::col_character(),
  contig_count = vroom::col_double(),
  N50 = vroom::col_double(),
  largest_contig = vroom::col_double(),
  total_size = vroom::col_double(),
  ambiguous_bases = vroom::col_character(),
  QC_warnings = vroom::col_character(),
  ST = vroom::col_character(),
  gapA = vroom::col_character(),
  infB = vroom::col_double(),
  mdh = vroom::col_double(),
  pgi = vroom::col_double(),
  phoE = vroom::col_character(),
  rpoB = vroom::col_double(),
  tonB = vroom::col_character(),
  YbST = vroom::col_character(),
  Yersiniabactin = vroom::col_character(),
  ybtS = vroom::col_character(),
  ybtX = vroom::col_character(),
  ybtQ = vroom::col_character(),
  ybtP = vroom::col_character(),
  ybtA = vroom::col_character(),
  irp2 = vroom::col_character(),
  irp1 = vroom::col_character(),
  ybtU = vroom::col_character(),
  ybtT = vroom::col_character(),
  ybtE = vroom::col_character(),
  fyuA = vroom::col_character(),
  spurious_ybt_hits = vroom::col_character(),
  CbST = vroom::col_double(),
  Colibactin = vroom::col_character(),
  clbA = vroom::col_character(),
  clbB = vroom::col_character(),
  clbC = vroom::col_character(),
  clbD = vroom::col_character(),
  clbE = vroom::col_character(),
  clbF = vroom::col_character(),
  clbG = vroom::col_character(),
  clbH = vroom::col_character(),
  clbI = vroom::col_character(),
  clbL = vroom::col_character(),
  clbM = vroom::col_character(),
  clbN = vroom::col_character(),
  clbO = vroom::col_character(),
  clbP = vroom::col_character(),
  clbQ = vroom::col_character(),
  spurious_clb_hits = vroom::col_character(),
  AbST = vroom::col_character(),
  Aerobactin = vroom::col_character(),
  iucA = vroom::col_character(),
  iucB = vroom::col_character(),
  iucC = vroom::col_character(),
  iucD = vroom::col_character(),
  iutA = vroom::col_character(),
  spurious_abst_hits = vroom::col_character(),
  SmST = vroom::col_character(),
  Salmochelin = vroom::col_character(),
  iroB = vroom::col_character(),
  iroC = vroom::col_character(),
  iroD = vroom::col_character(),
  iroN = vroom::col_character(),
  spurious_smst_hits = vroom::col_character(),
  RmST = vroom::col_character(),
  RmpADC = vroom::col_character(),
  rmpA = vroom::col_character(),
  rmpD = vroom::col_character(),
  rmpC = vroom::col_character(),
  spurious_rmst_hits = vroom::col_character(),
  virulence_score = vroom::col_double(),
  spurious_virulence_hits = vroom::col_character(),
  rmpA2 = vroom::col_character(),
  AGly_acquired = vroom::col_character(),
  Col_acquired = vroom::col_character(),
  Fcyn_acquired = vroom::col_character(),
  Flq_acquired = vroom::col_character(),
  Gly_acquired = vroom::col_character(),
  MLS_acquired = vroom::col_character(),
  Phe_acquired = vroom::col_character(),
  Rif_acquired = vroom::col_character(),
  Sul_acquired = vroom::col_character(),
  Tet_acquired = vroom::col_character(),
  Tgc_acquired = vroom::col_character(),
  Tmt_acquired = vroom::col_character(),
  Bla_acquired = vroom::col_character(),
  Bla_inhR_acquired = vroom::col_character(),
  Bla_ESBL_acquired = vroom::col_character(),
  Bla_ESBL_inhR_acquired = vroom::col_character(),
  Bla_Carb_acquired = vroom::col_character(),
  Bla_chr = vroom::col_character(),
  SHV_mutations = vroom::col_character(),
  Omp_mutations = vroom::col_character(),
  Col_mutations = vroom::col_character(),
  Flq_mutations = vroom::col_character(),
  truncated_resistance_hits = vroom::col_character(),
  spurious_resistance_hits = vroom::col_character(),
  resistance_score = vroom::col_double(),
  resistance_class_count = vroom::col_logical(),
  resistance_gene_count = vroom::col_logical(),
  wzi = vroom::col_character(),
  K_locus = vroom::col_character(),
  K_type = vroom::col_character(),
  K_locus_confidence = vroom::col_character(),
  K_locus_problems = vroom::col_character(),
  K_locus_identity = vroom::col_character(),
  K_Missing_expected_genes = vroom::col_character(),
  O_locus = vroom::col_character(),
  O_type = vroom::col_character(),
  O_locus_confidence = vroom::col_character(),
  O_locus_problems = vroom::col_character(),
  O_locus_identity = vroom::col_character(),
  O_Missing_expected_genes = vroom::col_character()
)

read_kleborate <- function(file, ..., qc=TRUE) {
  data <- exec2(vroom::vroom, ..., args=list(file=file, col_types=kleborate_column_spec))
  stopifnot(length(intersect(names(data), names(kleborate_column_spec$cols))) >=
              length(kleborate_column_spec$cols) / 10)  # 10% of expected cols
  return(
    data |>
      # dplyr::mutate(
      #   dplyr::across(where(is.character), ~dplyr::na_if(., "-")),
      # ) |>
      (\(.) if(isTRUE(qc)) klebnet_qc(.) else .)()
  )
}

summarise_kleborate <- function(data, cols=c('ST', 'K_locus', 'O_locus')) {
  data |>
    dplyr::summarise(
      dplyr::across(tidyselect::all_of(cols), dplyr::n_distinct)
    )
}

klebnet_qc <- function(data, max_contigs=500, min_size=4969898, max_size=6132846) {
  # https://bigsdb.pasteur.fr/klebsiella/genome-quality-check/
  # stopifnot(all(c() %in% names(data)))
  return(
    data |>
      dplyr::filter(
        dplyr::if_any(tidyselect::matches('locus_confidence'), ~.x == 'Typeable'),
        dplyr::if_any(tidyselect::matches('species_match'), ~.x != 'weak'),
        dplyr::if_any(tidyselect::matches('contig_count'), ~.x <= max_contigs),
        dplyr::if_any(tidyselect::matches('total_size'), ~dplyr::between(.x, min_size, max_size))
      )
  )
}


update_kaptive_kpsc_o_type <- function(data) {
  base::stopifnot(
    is.data.frame(data),
    "O_locus" %in% names(data),
    "O_type" %in% names(data)
  )
  
  data |> 
    dplyr::mutate(  # Change loci first
      O_locus = dplyr::case_match(
        .x=O_locus, .default=O_locus,
        "O1/O2v1" ~ "OL2α.1",
        "O1/O2v2" ~ "OL2α.2",
        "O1/O2v3" ~ "OL2α.3",
        "O3/O3a" ~ "OL3α/β",
        "O3b" ~ "OL3𝛾",
        "OL102" ~ "OL14",
        "OL103" ~ "OL10",
        "OL104" ~ "OL15"
      ),
      
      O_type = dplyr::case_when(
        .default=O_type,
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c('wbbY', 'wbbZ'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O1αβ,2α", 
        stringr::str_detect(O_locus, "OL2α") & stringr::str_detect(O_locus_extra_genes, "wbbY") ~ "O1α,2α",
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c('wbbY', 'wbbZ', 'gml2β'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O1αβ,2β", 
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c('wbbY', 'gml2β'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O1α,2β",
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c('wbbY', 'wbbZ', 'orf8'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O1αβ,2γ",
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c(), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O11αβ,2α",
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c(), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O11α,2α", 
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c(), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O11αβ,2β", 
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c(), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O11α,2β", 
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c(, 'orf8'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O11αβ,2γ ",
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c('a', 'b'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O2β",
        stringr::str_detect(O_locus, "OL2α") & all(purrr::map_lgl(c('orf8'), ~stringr::str_detect(O_locus_extra_genes, .x))) ~ "O2αγ ",
        stringr::str_detect(O_locus, "OL2α") ~ "O2α",
        O_locus == "OL3" & stringr::str_detect(O_locus_extra_genes, "") ~ "O3α + O3β", 
        O_locus == "OL3" & stringr::str_detect(O_locus_extra_genes, "") ~ "O3γ",
        O_locus == "OL4" ~ "O4",
        O_locus == "OL5" ~ "O5", 
        O_locus == "OL10" ~ "O10", 
        O_locus == "OL12" ~ "O12", 
        O_locus == "OL13" ~ "O13",  
        O_locus == "OL14" ~ "O14",
        O_locus == "OL15" ~ "O15",
      )
    )
}

# ========================== ===================================================== =============================================== ==============================================
# New serotype designation   Required genes/loci (implemented in Kaptive v.3.1+)   Prior Kaptive designation (v.2.0.8–v.3.0.0b6)   Prior Kaptive genes/loci (v.2.0.8–v.3.0.0b6)
# ========================== ===================================================== =============================================== ==============================================
# O1αβ,2α                    OL2α.(1/2/3), wbbYZ                                   O1ab                                            O1/O2v1, wbbYZ
# O1α,2α                     OL2α.(1/2/3), wbbY                                    O1a                                             O1/O2v1, wbbY
# O1αβ,2β                    OL2α.(1/2/3), gml2β, wbbYZ                            O1ab                                            O1/O2v2, wbbYZ
# O1α,2β                     OL2α.(1/2/3), gml2β, wbbY                             O1a                                             O1/O2v2, wbbY
# O1αβ,2γ                    OL2α.(1/2/3), orf8, wbbYZ                             O1ab                                            O1/O2v3, wbbYZ
# O2α                        OL2α.(1/2/3)                                          O2a                                             O1/O2v1
# O2β                        OL2α.(1/2/3), gml2β                                   O2afg                                           O1/O2v2
# O2αγ                       OL2α.(1/2/3), orf8                                    O2a                                             O1/O2v3
# O3α + O3β                  OL3α/β                                                O3/O3a                                          O3/O3a
# O3γ                        OL3γ                                                  O3b                                             O3b
# O4                         OL4                                                   O4                                              O4
# O5                         OL5                                                   O5                                              O5
# O10                        OL10                                                  OL103                                           OL103
# O11αβ,2α                   OL2α.(1/2/3), wbmVWX                                  O2ac                                            O1/O2v1, wbmVWX
# O11α,2α                    OL2α.(1/2/3), wbmVW                                   O2ac                                            O1/O2v1, wbmVW
# O11αβ,2β                   OL2α.(1/2/3), gml2β, wbmVWX                           O2ac                                            O1/O2v2, wbmVWX
# O11α,2β                    OL2α.(1/2/3), gml2β, wbmVW                            O2ac                                            O1/O2v2, wbmVW
# O11αβ,2γ                   OL2α.(1/2/3), orf8, wbmVWX                            O2ac                                            O1/O2v3, wbmVW
# O12                        OL12                                                  O12                                             O12
# O13                        OL13                                                  O13                                             OL13
# O14                        OL14                                                  OL102                                           OL102
# O15                        OL15                                                  OL104                                           OL104