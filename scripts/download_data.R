wd <- Sys.getenv("wd")
key <- Sys.getenv("api_key")

source(paste0(wd, "/src/functions.R"))
library(R.utils)

# Genes to download

library(org.Hs.eg.db)
library(reactome.db)

reactome_genes <- reactomeEXTID2PATHID |> as.list()
human <- sapply(reactome_genes, \(x){
  grepl('HSA', x) |> any()  
})
entrez <- reactome_genes[human] |> names() |> unique()
human_genes <- select(org.Hs.eg.db, keys=entrez, columns=c("SYMBOL"), keytype="ENTREZID")

# Query chatGPT and evaluate results

GPT_response <- NULL

while(length(GPT_response)<=100){
  
  #some APIs seem to crash sometimes. Restart after 1 muntes.
  withTimeout({
    symbol <- sample(na.omit(human_genes$SYMBOL), 1)
    print(symbol)
    
    disease <- gene_to_uniprot_disease(symbol) 
    
    #some genes have no diseases connected to them.
    if(!is.null(disease)){
      GPT_response[[symbol]] <- NULL
      GPT_response[[symbol]]$disease <- disease
      
      #to not pass the per minute query limits
      gpt <- try(GPT_query(symbol, key))
      if(class(gpt) == "try-error"){
        #Sometimes limit is reached in anycase. 
        gpt <- GPT_query(symbol, key)
      }
      
      #check if the response is as requested in terms of structure
      length_diseases <- sapply(gpt$disease, \(xx){
          stringr::str_count(xx, '\\w+')<30
        }) |> all()
      
      if(length_diseases){
        accuracy <- sapply(gpt$disease, GO = disease, mean_calc = T, word_by_word_lcs)
        
        gpt_refernce_predictions <- try(sapply(gpt$reference, pubmed_query))
        
        if(class(gpt_refernce_predictions) == "try-error"){
          gpt_refernce_predictions <- "pubmed query failed"
        }
        
        gene_popularity <- try(pubmed_query(gene = symbol))
        if(class(gene_popularity) == "try-error"){
          gene_popularity <- "pubmed query failed"
        }
        
        GPT_response[[symbol]]$accuracy <- accuracy
        GPT_response[[symbol]]$gene_references <- gene_popularity
        GPT_response[[symbol]]$gpt_out <- gpt
        GPT_response[[symbol]]$gpt_refernce_predictions <- gpt_refernce_predictions
      }else{
        GPT_response[[symbol]] <- NULL
        print(gpt$disease)
      }
      
    }else{
      print("No gene context")
    }
    
  }, timeout = 60, onTimeout = "silent")
}

# Cleaning up the data

# sometimes the response has padding in the form of GENE plays a role \n \n. Below cleans that up.
plays_a_role <- sapply(GPT_response, \(x){
  nchar(x$gpt_out$disease)[1]>60 & nchar(x$gpt_out$disease)[2]==0
})

for(gene in which(plays_a_role)){
  GPT_response[[gene]]$gpt_out$disease <- GPT_response[[gene]]$gpt_out$disease[-1:-2]
  GPT_response[[gene]]$accuracy <- GPT_response[[gene]]$accuracy[-1:-2]
}

#droping any edge cases
GPT_response <- GPT_response[sapply(GPT_response, length)>1]

save(GPT_response, file = paste0(wd, '/data/GPT_response.Rdata'))
