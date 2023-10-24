library(rjson)
library(openai)
# library(svMisc)
library(R.utils)

pubmed_query <- function(title = "Time-restricted feeding alters lipid and amino acid metabolite rhythmicity without perturbing clock gene expression",
                         gene = NULL){
  if(is.null(gene)){
    pubmed_api_query <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term=”ARTICLETITLE”[Title:~0]"
    pubmed_api_query <- gsub('ARTICLETITLE', title, pubmed_api_query, fixed = TRUE)
  }else{
    pubmed_api_query <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term='GENE'"#
    pubmed_api_query <- gsub('GENE', gene, pubmed_api_query, fixed = TRUE)
  }
  
  pubmed_api_query <- gsub(" ", "%20", pubmed_api_query)
  reference <- withTimeout({
    try(fromJSON(file=pubmed_api_query, simplify=FALSE))
    }, timeout = 60, onTimeout = "silent")
  
  if(class(reference) == "try-error"){
    reference <- NA
  }else{
    reference <- as.numeric(reference$esearchresult$count) 
  }
  
  return(reference)
}

reformat <- function(x){
  x <- gsub('"', "", x, fixed = T)
  x <- gsub('\n', "---", x, fixed = T)
  x <- strsplit(x, "---", fixed = T)
  x <- sapply(x, trimws)
  x <- sapply(x, tolower)
  x <- sapply(x, \(xx){
    xx <- gsub("[[:digit:]]\\.", "", xx)
    xx <- gsub("\\.", "", xx)
    xx <- gsub("- ", "", xx)
    xx
  })
  x <- sapply(x, trimws)
  names(x) <- NULL
  return(x)
}

GPT_query_list <- function(symbol = "GLUT4", disease_list, api_key, model = "gpt-3.5-turbo", full = F){
  
  #adding a time-out because it seems to crash some times
  withTimeout({
    prompt_prime <- "Answer as an expert medical doctor, molecular biologist, and geneticist."
    
    prompt_disease <- "Which disease from the diseaselist below is gene GENESYMBOL involved in? Reply only with the name of the disease. Diseaselist:"
    prompt_reference <- "Give me 5 or fewer scientific literature references to further my understanding of gene GENESYMBOL. Answer only with article title. Do not elaborate. Do not include authors. Separate by new line. Text:\n"
    
    prompt_disease <- gsub("GENESYMBOL", symbol, prompt_disease)
    prompt_reference <- gsub("GENESYMBOL", symbol, prompt_reference)
    
    prompt_gene_disease <- paste0(prompt_disease, paste0(" ", disease_list, collapse = ", "))
    
    # print(prompt_gene_disease)
    
    gpt_answer_disease <- create_chat_completion(openai_api_key = api_key,
                                                 model = model,
                                                 # temperature = 0.5,
                                                 messages = list(list("role" = "system",
                                                                      "content" = prompt_prime),
                                                                 list("role" = "user",
                                                                      "content" = prompt_gene_disease)
                                                 )
    )
    
    gpt_answer_reference <- create_chat_completion(openai_api_key = api_key,
                                                   model = "gpt-3.5-turbo",
                                                   # temperature = 0.5,
                                                   messages = list(list("role" = "system",
                                                                        "content" = prompt_prime),
                                                                   list("role" = "user",
                                                                        "content" = prompt_reference)
                                                   )
    )
    gpt_answer_disease_reformat <- gsub(".*: |\\.","", gpt_answer_disease$choices$message.content)
    gpt_answer_reference <- reformat(gpt_answer_reference$choices$message.content)
    
    if(full){
      gpt_answer <- list(reference = gpt_answer_reference,
                         disease = gpt_answer_disease_reformat,
                         full_response = gpt_answer_disease)
    }else{
      gpt_answer <- list(reference = gpt_answer_reference,
                         disease = gpt_answer_disease_reformat) 
    }
    
    return(gpt_answer)
  }, timeout = 60, onTimeout = "silent")
}

GPT_query_freeform <- function(symbol = "GLUT4", disease_list, api_key, model = "gpt-3.5-turbo", full = F){
  
  #adding a time-out because it seems to crash some times
  withTimeout({
    prompt_prime <- "Answer as an expert medical doctor, molecular biologist, and geneticist."
    prompt_gene_disease <- 'Which relevant conditions, specific diseases, or syndromes does the gene GENESYMBOL play a role in. Format answer only with a list of specific names separated by new line without any elaboration. Limit your answer only 1 disease, the most common.'
    prompt_gene_disease <- gsub("GENESYMBOL", symbol, prompt_gene_disease)
    
    gpt_answer_disease <- create_chat_completion(openai_api_key = api_key,
                                                 model = model,
                                                 # temperature = 0.5,
                                                 messages = list(list("role" = "system",
                                                                      "content" = prompt_prime),
                                                                 list("role" = "user",
                                                                      "content" = prompt_gene_disease))
    )
    gpt_answer <- gpt_answer_disease$choices$message.content
    
    return(gpt_answer)
  }, timeout = 60, onTimeout = "silent")
}

generate_lists <- function(interval, disease_table){
  disease_counts <- disease_table$disease |> table() |> sort(decreasing = T)
  
  disease_interval <- quantile(1:length(disease_counts), na.rm = T, probs = c(interval[1],interval[2]))
  disease_interval <- disease_interval[1]:disease_interval[2] |> ceiling()
  disease_subset <- names(disease_counts)[disease_interval]
  
  gene_counts <- disease_table[disease %in% disease_subset, symbol] |> table() |> sort(decreasing = T)
  gene_interval <- quantile(1:length(gene_counts), na.rm = T, probs = c(interval[1],interval[2]))
  gene_interval <- gene_interval[1]:gene_interval[2] |> ceiling()
  gene_subset <- names(gene_counts)[disease_interval]
  
  subset <- list(
    disease = disease_subset,
    genes = gene_subset
  )
  
  return(subset)
}

manual_review <- function(GTP_response){
  for(i in names(GTP_response)){
    response <- GTP_response[[i]]$disease
    cat(response)
    disease <- readline()
    
    if(nchar(disease)==1 & !is.na(suppressWarnings(as.numeric("1")))){
      disease <- as.numeric(disease)
      diseases_words <- strsplit(GTP_response[[i]]$disease, " ")[[1]]
      start <- length(diseases_words) - (disease - 1)
      end <- length(diseases_words)
      
      GTP_response[[i]]$disease <- paste0(diseases_words[start:end], collapse = " ")
      } else if (nchar(disease)>0){
      GTP_response[[i]]$disease <- disease
    }
    
    cat(paste0(which(i == names(GTP_response)), " ------------\n"))
  }
  
  return(GTP_response)
}

evaluate_responses <- function(GPT_response, gene){
  
  GPT_response$gene_hits <- pubmed_query(gene = gene)
  GPT_response$disease_correct <- any(GPT_response$disease %in% diseases[symbol == gene, disease])
  GPT_response$reference_correct <- sapply(GPT_response$reference, pubmed_query, USE.NAMES = F)
  
  return(GPT_response)
}

evaluate_freeform <- function(GPT_response){
  
  for(gene in names(GPT_response)){
    print("\n")
    print(diseases[symbol == gene, sort(unique(disease))])
    print("--------------")
    print(GPT_response[[gene]])
    i <- which(names(GPT_response) == gene)
    
    disease <- GPT_response[[gene]]
    GPT_response[[gene]] <- NULL#reformating for nested list
    GPT_response[[gene]]$disease <- disease
    
    print(paste0(i, "--------------"))
    response_accuracy <- readline("True/False\n")
    GPT_response[[gene]]$accuracy <- response_accuracy
    }
  
  
  return(GPT_response)
}

reformat_results <- function(res_3.5, res_4){
  res_3.5 <- sapply(res_3.5, \(x){
    c(correct = x$disease_correct, 
      gene_hits = x$gene_hits, 
      ref_correct = sum(x$reference_correct>0)/length(x$reference_correct),
      model = "GPT-3.5")
  }) |> t() |> data.frame() 
  res_4 <- sapply(res_4, \(x){
    c(correct = x$disease_correct, 
      gene_hits = x$gene_hits, 
      ref_correct = sum(x$reference_correct>0)/length(x$reference_correct),
      model = "GPT-4")
  }) |> t() |> data.frame() 
  res <- rbind(
    res_3.5, res_4
  )
  return(res)
}
