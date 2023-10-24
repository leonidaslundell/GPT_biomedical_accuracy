wd <- Sys.getenv("wd")
data <- paste0(wd, "/data/")
api_key <- Sys.getenv("api_key")

source(paste0(wd, "/src/functions.R"))
library(data.table)
library(progress)

####################################################
#Downloading data from a textmining resource
# from https://diseases.jensenlab.org/Downloads
# https://www.sciencedirect.com/science/article/pii/S1046202314003831

#using it in the broadest sense without any prefiltering in confidence

diseases <- fread(paste0(data, "human_disease_textmining_filtered.tsv"))[,1:4]
colnames(diseases) <- c("identifier", "symbol", "DO_ID", "disease")
#drop noncoding
diseases <- diseases[!grepl("rRNA|hsa_circ|sno|hsa|AS|orf|ENS|\\.|RNA", symbol)]
#Obsolete DOID as as as Type 2 diabebets labeled
diseases[disease == "DOID:0081062", "disease"] <- "Type 2 diabetes mellitus"

####################################################
#donwload data from chatGPT
n_genes <- 100

#The openai api calls stall for some weird reason. Manual work around using the function below.
#stupid function to restart the loops...
remaining <- function(x){
  sapply(x, \(x) class(x) == 'try-error' | is.null(x)) |> which() |> names()
}

#############
#well known genes (top 5%)

well_known_subset <- generate_lists(c(0,0.05), diseases)
well_known_subset$genes <- sample(well_known_subset$genes, n_genes, replace = F)

well_known_genes_3.5 <- well_known_genes_4 <- well_known_freeform <- vector(mode = "list", length = n_genes)
names(well_known_genes_3.5) <- names(well_known_genes_4) <- names(well_known_freeform) <- well_known_subset$genes

pb <- progress_bar$new(total = 100)
for(gene in remaining(well_known_genes_3.5)){
  pb$tick()
  Sys.sleep(2)
  well_known_genes_3.5[[gene]] <- try(GPT_query_list(gene, disease_list = well_known_subset$disease, api_key = api_key, model = "gpt-3.5-turbo"))
}

pb <- progress_bar$new(total = 100)
for(gene in remaining(well_known_genes_4)){
  pb$tick()
  Sys.sleep(2)
  well_known_genes_4[[gene]] <- try(GPT_query_list(gene, disease_list = well_known_subset$disease, api_key = api_key, model = "gpt-4"))
}

pb <- progress_bar$new(total = 100)
for(gene in remaining(well_known_freeform)){
  pb$tick()
  Sys.sleep(2)
  well_known_freeform[[gene]] <- try(GPT_query_freeform(gene, api_key = api_key, model = "gpt-4"))
}

#############
#intermediate known genes, 20%-25%

intermediate_known_subset <- generate_lists(c(.20, 0.25), diseases)
remove <- intersect(intermediate_known_subset$genes, well_known_subset$genes)
intermediate_known_subset$genes <- intermediate_known_subset$genes[!intermediate_known_subset$genes %in% remove]
intermediate_known_subset$genes <- sample(intermediate_known_subset$genes, n_genes, replace = F)

intermediate_known_genes_3.5 <- intermediate_known_genes_4 <- vector(mode = "list", length = n_genes)
names(intermediate_known_genes_3.5) <- names(intermediate_known_genes_4) <- intermediate_known_subset$genes

pb <- progress_bar$new(total = 100)
for(gene in remaining(intermediate_known_genes_3.5)){
  pb$tick()
  Sys.sleep(2)
  intermediate_known_genes_3.5[[gene]] <- try(GPT_query_list(gene, disease_list = intermediate_known_subset$disease, api_key = api_key, model = "gpt-3.5-turbo"))
}

pb <- progress_bar$new(total = 100)
for(gene in remaining(intermediate_known_genes_4)){
  pb$tick()
  Sys.sleep(2)
  intermediate_known_genes_4[[gene]] <- try(GPT_query_list(gene, disease_list = intermediate_known_subset$disease, api_key = api_key, model = "gpt-4"))
}

# save(list = c("well_known_genes_3.5", "well_known_genes_4", "intermediate_known_genes_3.5", "intermediate_known_genes_4", "well_known_freeform"),
#      file = paste0(data, "GPT_raw_response.Rdata"))

####################################################
#Manually review all responses to format them so that they are comparable

intermediate_known_genes_3.5 <- manual_review(intermediate_known_genes_3.5)
intermediate_known_genes_4 <- manual_review(intermediate_known_genes_4)
well_known_genes_3.5 <- manual_review(well_known_genes_3.5)
well_known_genes_4 <- manual_review(well_known_genes_4)#when several diseases predicted pick last

# save(list = c("well_known_genes_3.5", "well_known_genes_4", "intermediate_known_genes_3.5", "intermediate_known_genes_4"),
#      file = paste0(data, "GPT_manual_review_response.Rdata"))

####################################################
#evaluate response
load(paste0(data, "GPT_manual_review_response.Rdata"))

pb <- progress_bar$new(total = 100)
for(gene in names(intermediate_known_genes_3.5)){
  pb$tick()
  intermediate_known_genes_3.5[[gene]] <- evaluate_responses(intermediate_known_genes_3.5[[gene]], gene)
}

pb <- progress_bar$new(total = 100)
for(gene in names(intermediate_known_genes_4)){
  pb$tick()
  intermediate_known_genes_4[[gene]] <- evaluate_responses(intermediate_known_genes_4[[gene]], gene)
}

pb <- progress_bar$new(total = 100)
for(gene in names(well_known_genes_3.5)){
  pb$tick()
  well_known_genes_3.5[[gene]] <- evaluate_responses(well_known_genes_3.5[[gene]], gene)
}

pb <- progress_bar$new(total = 100)
for(gene in names(well_known_genes_4)){
  pb$tick()
  well_known_genes_4[[gene]] <- evaluate_responses(well_known_genes_4[[gene]], gene)
}

well_known_freeform <- evaluate_freeform(well_known_freeform)

#############
#comments on how I associated. This is not everything...

# T2D == metabolic syndrome
# Alzehimer == dementia
# Congenital Secretory Sodium Diarrhea
# Major Depressive Disorder == Mood disorder
# Tuberous sclerosis == Tuberous Sclerosis Complex
# Rheumatic disease == Rheumatic arthitis
# Epilepsy == Action myoclonus-renal failure syndrome
# Hemolytic anemia == anemia

#intersting how close: Arginase deficiency, Argininosuccinic aciduria
# intersting variation... Diabetes Mellitus Type 2

# save(list = c("well_known_genes_3.5", "well_known_genes_4", "intermediate_known_genes_3.5", "intermediate_known_genes_4", "well_known_freeform"),
#      file = paste0(data, "GPT_response_full.Rdata"))

####################################################
# Reformat into tidy
#load(paste0(data, "GPT_response_full.Rdata"))

#typo in T/F
well_known_freeform$CSF3$accuracy <- well_known_freeform$IL18$accuracy <- T

well_known <- reformat_results(well_known_genes_3.5, well_known_genes_4)
intermediate_known <- reformat_results(intermediate_known_genes_3.5, intermediate_known_genes_4)
#note that there is a major warning here, rownames have been deduplicated by adding 1 at the end...
# all(stringr::str_sub(rownames(well_known)[101:200], end = -2) == rownames(well_known)[1:100])

sapply(well_known_freeform, \(x) as.logical(x$accuracy)) |> table()

well_known_freeform <- sapply(well_known_freeform, \(x) as.logical(x$accuracy))

# save(list = c("well_known","intermediate_known","well_known_freeform"), 
#      file = paste0(data, "GPT_response.Rdata"))