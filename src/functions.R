library(rjson)
library(openai)
library(svMisc)

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
  reference <- fromJSON(file=pubmed_api_query, simplify=FALSE)
  reference <- as.numeric(reference$esearchresult$count)
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

GPT_query <- function(symbol = "GLUT4", api_key){
  
  prompt_prime <- "Answer as an expert medical doctor, molecular biologist, and geneticist."
  prompt_disease <- "Which clinically relevant conditions, specific diseases, or syndromes does the gene GENESYMBOL play a role in."
  prompt_disease_followup <- "From the text below, give me only the conditions, specific diseases, or syndromes. Format answer only with a list of specific names separated by new line without any elaboration. Text:\n"
  
  prompt_reference <- "For the text below, give me up to 5 scientific literature references to further my understanding about GENESYMBOL. Answer only with article title. Do not elaborate. Do not include authors. Separate by new line. Text:\n"
  
  prompt_disease <- gsub('GENESYMBOL', symbol, prompt_disease)
  prompt_reference <- gsub('GENESYMBOL', symbol, prompt_reference)
  
  gpt_answer_disease <- create_chat_completion(openai_api_key = api_key,
                                               model = "gpt-3.5-turbo",
                                               temperature = 0.5,
                                               messages = list(list("role" = "system",
                                                                    "content" = prompt_prime),
                                                               list("role" = "user",
                                                                    "content" = prompt_disease)
                                               )
  )
  
  prompt_disease_followup_gene <- paste0(prompt_disease_followup, gpt_answer_disease$choices$message.content)
  gpt_answer_disease_followup <- create_chat_completion(openai_api_key = api_key,
                                                        model = "gpt-3.5-turbo",
                                                        temperature = 0.5,
                                                        messages = list(list("role" = "system",
                                                                             "content" = prompt_prime),
                                                                        list("role" = "user",
                                                                             "content" = prompt_disease_followup_gene)
                                                        )
  )
  
  prompt_reference_gene <- paste0(prompt_reference, gpt_answer_disease$choices$message.content)
  gpt_answer_reference <- create_chat_completion(openai_api_key = api_key,
                                                 model = "gpt-3.5-turbo",
                                                 temperature = 0.5,
                                                 messages = list(list("role" = "system",
                                                                      "content" = prompt_prime),
                                                                 list("role" = "user",
                                                                      "content" = prompt_reference_gene)
                                                 )
  )
  gpt_answer_disease_list <- reformat(gpt_answer_disease_followup$choices$message.content)
  gpt_answer_reference <- reformat(gpt_answer_reference$choices$message.content)
  
  gpt_answer <- list(reference = gpt_answer_reference,
                     disease = gpt_answer_disease_list,
                     full_response = gpt_answer_disease$choices$message.content)
  
  return(gpt_answer)
}

gene_to_uniprot_disease <- function(symbol = "GLUT4"){
  symbol <- toupper(symbol)
  query <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+", 
    symbol,
    "&format=json"
  )
  uniprot <- fromJSON(file=query, simplify=FALSE)
  
  disease <- sapply(uniprot[["results"]], \(res){
    sapply(res$comments, \(comment){
      c(comment$disease$diseaseId, comment$disease$acronym)
    })
  }) |> unlist(use.names = FALSE)
  
  return(disease)
}

word_by_word_lcs <- function(GPT, GO, mean_calc = T){
  GPT_word <- strsplit(tolower(GPT), "\\s+")[[1]]
  
  GO <- paste0(gsub(',', '', tolower(GO)), collapse = ' ')
  
  GPT_GO_sim <- sapply(GPT_word, \(x){
    GPT_similarity <- CEGO::distanceStringLCStr(x, GO)
    GPT_similarity <- nchar(GO) - GPT_similarity
    GPT_similarity <- GPT_similarity/nchar(x)
    GPT_similarity
  })
  if(mean_calc){
    GPT_GO_sim <- mean(GPT_GO_sim)
  }
  
  return(GPT_GO_sim)
}