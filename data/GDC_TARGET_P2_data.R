
# Cleaning Data -----------------------------------------------------------

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(progressr)
library(tibble)

handlers(global = TRUE)
handlers("txtprogressbar")

#clinical
clinical = fread("../Target-ALL-P2/clinical.cart.2025-01-03/clinical.tsv")

#sample info
sample_sheet = fread("../Target-ALL-P2/gdc_sample_sheet.2025-01-03.tsv")

#biospecimen info ince some samples ahve multiple files
biospecimen = fread("../Target-ALL-P2/biospecimen.cart.2025-01-03/aliquot.tsv")

#manifest file
manifest = fread("../Target-ALL-P2/gdc_manifest.2025-01-03.095525.txt")

#since there are multiple aliquots for the same sample, will sort lexicographically and pick highest
#https://www.biostars.org/p/311017/
sample_sheet = sample_sheet[match(stringr::str_sort(sample_sheet$`File ID`), sample_sheet$`File ID`),] %>%
  group_by(`Case ID`) %>%
  slice(1) %>%
  filter(`Case ID` %in% clinical$case_submitter_id)

# Key Genes for TARGET-ALL-P2
# 
# IKZF1 - Frequently altered in ALL and associated with poor prognosis.
# PAX5 - A key regulator of B-cell development; often mutated in ALL.
# CD19 - B-cell marker, relevant for CAR T-cell therapy.
# CD10 (MME) - Common marker for pre-B ALL.
# TP53 - Tumor suppressor gene, mutations often associated with poor outcomes.
# NOTCH1 - Frequently mutated in T-cell ALL.
# CRLF2 - Overexpression linked to high-risk B-cell ALL.
# ABL1 - Associated with the BCR-ABL fusion, common in Ph-positive ALL.
# BCR - Partner in the BCR-ABL translocation.
# ETV6 - Involved in translocations leading to ETV6-RUNX1 fusions, a common driver in ALL.
# RUNX1 - Partner in the ETV6-RUNX1 fusion.
# JAK2 - Mutations lead to activation of the JAK-STAT pathway, relevant in some ALL subtypes.
# IL7R - Important for T-cell development, mutated in some ALL cases.
# CD3E - A marker for T cells, helpful in differentiating T-ALL.
# CEBPA - Regulates hematopoiesis, alterations seen in some leukemia cases.
# FLT3 - Often mutated in various leukemias.
# NRAS - Oncogene, mutations drive proliferation in some ALL cases.
# KRAS - Similar to NRAS, involved in proliferation signaling.
# CXCR4 - Chemokine receptor influencing cell migration and homing.
# MYC - Oncogene driving cell proliferation, frequently overexpressed in ALL.
# 
# Why These Genes?
#   
# These genes cover:
# Markers for Diagnosis: CD19, CD10, CD3E.
# Prognostic Factors: IKZF1, TP53, CRLF2.
# Therapeutic Targets: BCR-ABL, JAK2, IL7R.
# Pathway Insights: NOTCH1, MYC, NRAS, KRAS.

#important genes in TARGET-ALL
genes <- c("IKZF1", "PAX5", "CD19", "MME", "TP53", "NOTCH1", "CRLF2", "ABL1", 
           "BCR", "ETV6", "RUNX1", "JAK2", "IL7R", "CD3E", "CEBPA", "FLT3", 
           "NRAS", "KRAS", "CXCR4", "MYC")

#process all samples
with_progress({
  p = progressor(along = 1:nrow(sample_sheet))
  exp_data = lapply(seq(nrow(sample_sheet)), function(i){
    x = sample_sheet$`File Name`[i]
    key = basename(x)
    case_id = sample_sheet$`Case ID`[sample_sheet$`File Name` == key]
    p(sprintf("Processing iteration %s", case_id))
    tmp = fread(file.path("../Target-ALL-P2/gdc_download_20250103_133609.419604", sample_sheet$`File ID`[i], x)) %>%
      filter(!is.na(tpm_unstranded)) %>%
      select(gene_id, gene_name, gene_type, tpm_unstranded) %>%
      rename(!!case_id := tpm_unstranded) %>%
      filter(gene_name %in% genes)
    return(tmp)
  })
})

#merge all together
tpm_dat = exp_data %>%
  purrr::reduce(inner_join, by = join_by(gene_id, gene_name, gene_type)) %>%
  filter(gene_id != "ENSG00000205755.12_PAR_Y") %>% #remove the gene that looks to be on a contig
  column_to_rownames("gene_name") %>%
  select(-gene_id, -gene_type) %>%
  t() %>% 
  data.table(keep.rownames = TRUE) %>%
  rename("case_submitter_id" = 'rn')

#save data
fwrite(tpm_dat, "data/TARGET_P2_tpm.csv")



#remove clinical columns with all missing value
clinical = clinical[,!apply(clinical, 2, function(x) all(x == "'--")), with = FALSE]
fwrite(clinical, "data/TARGET_P2_clinical.tsv")

#followup
followup = fread("../Target-ALL-P2/clinical.cart.2025-01-03/follow_up.tsv")
followup = followup[,!apply(followup, 2, function(x) all(x == "'--")), with = FALSE]

#fileID
#d4ddbad4-3ec4-48d2-9607-db615731ff6d