install.packages("taxize")
install.packages("myTAI")
library("taxize")
library("myTAI")
library(dplyr)


lst <- read.csv("/home/sbi6dap/Projects/NON-ATHAL/BI31007/REFS/genera.txt")
lst <- read.csv("/home/sbi6dap/species.txt", header=FALSE)

shortlist <- head(lst,50)

taxonomy( organism = "Utamphorophora" , db = "ncbi", output = "classification" )

outputlst <- apply(lst, 1, function(x) taxonomy( organism = x , db = "ncbi", output = "classification" ))

taxdata = data.frame()
for(x in 1:length(outputlst))
  {
  tryCatch({
      phylum=filter(outputlst[[x]], rank =="phylum")$name
      class=filter(outputlst[[x]], rank =="class")$name
      order=filter(outputlst[[x]], rank =="order")$name
      family=filter(outputlst[[x]], rank =="family")$name
      genus=filter(outputlst[[x]], rank =="genus")$name
      species=filter(outputlst[[x]], rank =="species")$name

    row <- data.frame(cbind(phylum=phylum,class=class,order=order,family=family,genus=genus))
    #row <- bind_cols(list(phylum=phylum,
     #                     if(length(class) == 0){ "NA" }else{ class=class},
      #                    if(length(order) == 0){ "NA" }else{ order=order},
       #                   if(length(family) == 0){ "NA" }else{ family=family},
        #                  if(length(genus) == 0){ "NA" }else{ genus=genus}))
    #print(row)
    
    taxdata <- bind_rows(taxdata, row)    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  
}


write.csv(taxdata, file="taxdata-full.txt")


