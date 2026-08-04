# The data will still change a bit. This scripts cleans it up, imports the latest version and removes the old files

# Set working directory
setwd("C:/Users/vivie/OneDrive - Université de Fribourg/Article_Screening/Collection_antagonistic_activity_against_6_pathogens")

# Load packages
library(readxl)
library(tidyverse)


# Import full-collection phylogenetic tree
## Check previous version
if(file.exists("data/0_full_tree.nwk")){
  file.rename("data/0_full_tree.nwk", 
              "data/99_old_full_tree.nwk")
}

## Copy full-collection phylogenetic tree
file.copy("../phylogenetic tree creation/bactbank_tree.nwk",
          "data/0_full_tree.nwk")

# Import filamentous strains phylogenetic tree
## Check previous version
if(file.exists("data/0_filamentous_tree.nwk")){
  file.rename("data/0_filamentous_tree.nwk", 
              "data/99_old_filamentous_tree.nwk")
}

## Copy full-collection phylogenetic tree
file.copy("../phylogenetic tree creation/actino_tree/actino_tree.nwk",
          "data/0_filamentous_tree.nwk")


# Clean and import data from filamentous strains
## Checks if a previous version of the data exists. If so, renames it.
if(file.exists("data/0_actino_score_database_cleaned.csv")){
  file.rename("data/0_actino_score_database_cleaned.csv", 
              "data/99_old_actino_score_database.csv")
}

if(file.exists("data/0_strain_metadata.csv")){
  file.rename("data/0_strain_metadata.csv", 
              "data/99_strain_metadata.csv")
}


## List strains to be removed
strains_remove <- c("AR012b", "AR031a", "AR043a", "AR043b", "AR113a", "AS041a", 
                    "AS076a", "AS076b", "AS097a1", "AS097b", "AS104a", "BL038b",
                    "BL038d", "BL038a", "BL038c", "BL047b2", "BL047b3", "BL047b4",
                    "BR084a", "BR084b", "BR091a", "BR091b", "BR091c", "BR092c",
                    "BR092d", "BR164b", "BS093a2", "BS107", "BS126a", "BS137", "BS152a", 
                    "BS195", "BS216", "BS217a", "BS022b", "BS233b", "BS219a",
                    "AR010"
)


# List strains to rename
strains_to_rename <- c("AR012a", "AR031b", "AR113b", "AS041b", "AS097a2", "AS104b",
                       "BL047b1", "BR164a", "BS093a1", "BS126b", "BS152b", "BS217b",
                       "BS233a")


# List the new names of the strains
strains_new_names <- c("AR012", "AR031", "AR113", "AS041", "AS097", "AS104",
                       "BL047b", "BR164", "BS093a", "BS126", "BS152", "BS217",
                       "BS233")



# Import filamentous data
tbl_score_old <- readxl::read_excel( "../actino_results_db.xlsx", 
                                     sheet = "score" )

# import strain characeristics
tbl_info <- readxl::read_excel( "../results_db.xlsx", 
                                sheet = "strain_characteristics" ) |> 
  mutate(compartment = replace_values(compartment, "Rhizospheric soil" ~ "Rhizosphere soil"),
         Medium = replace_values(Medium,
                                 "actino" ~ "ActIA"),
         Medium = factor(Medium,
                         levels = c("ARE", "TSA", "ActIA", "PDA"))
         )

# Remove unwanted strains, and rename what should be renamed.
tbl_score_new <- tbl_score_old |> 
  # Remove contaminated strains AR010, and the strains that we want to remove because they are not the original
  filter(!(strain_name %in% strains_remove)) |> 
  # Rename strains
  mutate(strain_name = replace_values(strain_name,
                                      from = strains_to_rename,
                                      to = strains_new_names)) |>
  # Add isolation and phylogenetic information
  left_join(tbl_info, by = join_by("strain_name" == "strain")) |>
  ## Format phylogeny and factors
  mutate( Genus = if_else( 
    # Fill missing phylogeny
    is.na(Genus),
    paste0("NA_", Family),
    Genus),
    # Limit taxa name length to 15
    Genus = if_else(
      nchar(Genus) > 15,
      str_c(str_sub(Genus, 1, 14), "."),
      Genus),
    patho_factor = recode_values( pathogen,
                                  "Pythium" ~ "Pult",
                                  "Dickeya" ~ "Dsol",
                                  "Alternaria" ~ "Asol",
                                  "Rhizoctonia" ~ "Rsol",
                                  "Pectobacterium" ~ "Pcar",
                                  "Phytophthora" ~ "Pinf",
                                  default = "error"
    ) |> ## Order the pathogens
      factor(levels = c("Dsol", "Pcar", "Asol",
                        "Rsol", "Pinf", "Pult" )),
    strain = strain_name 
    ) |> 
  ## Remove columns that won't be used
  select( -RDP_NCBI_closest_hit,
          -Antibiotic,
          -Sample,
          -strain_name )
  

# Export actino results table
write.csv(tbl_score_new, 
          "data/0_actino_score_database_cleaned.csv",
          row.names = FALSE)

# Export info table
write.csv(tbl_info, 
          "data/0_strain_metadata.csv",
          row.names = FALSE)

# Import data from non-filamentous strains

## Checks if a previous version of the data exists. If so, renames it.
if(file.exists("data/0_nonfilamentous_database_cleaned.csv")){
  file.rename("data/0_nonfilamentous_database_cleaned.csv",
              "data/99_old_nonfilamentous_database_cleaned.csv")
}

## Import all the tables before export as csv

### Score
tbl_score_nf <- readxl::read_excel( "../results_db.xlsx", 
                                    sheet = "score" )

# Plate-strain correspondance
tbl_plate_strain_nf <- readxl::read_excel(
  "../results_db.xlsx", 
  sheet = "plate_strain_corr")

# import scores of the controls
tbl_ctrl_score_nf <- readxl::read_excel( "../results_db.xlsx", 
                                         sheet = "score_ctrl" )


## Add strain name to score table
### Create join keys
tbl_score_nf = tbl_score_nf |> 
  mutate(key_well = paste( plate, well ),
         key_ctrl = paste( plate, well, start_date ))

tbl_plate_strain_nf = tbl_plate_strain_nf |> 
  mutate(key_well = paste(plate, well))

tbl_ctrl_score_nf = tbl_ctrl_score_nf |> 
  mutate(key_ctrl = paste(plate, well, date))


## Add strain name to the score table
tbl_score_nf <- tbl_score_nf |> left_join(tbl_plate_strain_nf[, c("strain",
                                                                  "key_well")],
                                          by = "key_well")

## Add control score to the score table
tbl_score_nf <- tbl_score_nf |> left_join(tbl_ctrl_score_nf[, c("ctrl_growth",
                                                                "key_ctrl")],
                                          by = "key_ctrl")

## List contaminated strains
contaminations <- c("BR014a", "BS201", "BS222", "AL056", "BR046", "BR047")

## Filter and transform variables
tbl_score_nf <- tbl_score_nf |> 
  ## Remove rows corresponding to empty wells
  filter( strain != "empty", 
          ## Remove invalid strains
          !(strain %in% c("AR017", "AR201", "BL078")),
          !(strain %in% contaminations)) |>
  ## Replace short names by full genus
  mutate(
    patho_factor = recode_values( pathogen,
                                  "Pythu" ~ "Pult",
                                  "Dickeya" ~ "Dsol",
                                  "Altsol" ~ "Asol",
                                  "Rhisol" ~ "Rsol",
                                  "Pecto" ~ "Pcar",
                                  "Pinf" ~ "Pinf",
                                  default = "error"
    ) |> ## Order the pathogens
      factor(levels = c("Dsol", "Pcar", "Asol",
                        "Rsol", "Pinf", "Pult" )),
    ## Adapt old scoring system to the new one
    inhibition_score = case_when(
      (patho_factor == "Asol" | patho_factor == "Rsol") & inhibition_score == 0 ~ 0,
      (patho_factor == "Asol" | patho_factor == "Rsol") & inhibition_score %in% 1:3 ~ inhibition_score - 1,
      .default = inhibition_score
    )
  )

## Add metadata
tbl_score_nf <- tbl_score_nf |>
  # Add sequencing and phylogenetic information
  left_join(tbl_info, by = "strain") |>
  ## Format phylogeny and factors
  mutate( Genus = if_else( 
    # Fill missing phylogeny
    is.na(Genus),
    paste0("NA_", Family),
    Genus),
    # Limit taxa name length to 15
    Genus = if_else(
      nchar(Genus) > 15,
      str_c(str_sub(Genus, 1, 14), "."),
      Genus)) |> 
  ## Remove columns that won't be used
  select( -key_well, 
          -key_ctrl,
          -RDP_NCBI_closest_hit,
          -Antibiotic,
          -Sample)

# Export updated and cleaned data
write.csv(tbl_score_nf, 
          "data/0_nonfilamentous_database_cleaned.csv",
          row.names = FALSE)
