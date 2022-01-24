idir <- '/home/pvermees/Documents/geotopes/extra/'
odir <- '/home/pvermees/Desktop/GEOL0017/'

solution <- FALSE
results <- FALSE
multidoc <- FALSE
randomise <- FALSE

preamble <- paste0("---\n",
ifelse(solution,ifelse(multidoc,"",
       "title: Solutions to additional exercises\n"),
       "title: Additional exercises\n"),
"output:
  pdf_document:
    includes:
      keep_tex: yes
    number_sections: yes
---

```{r settings, include = FALSE}
knitr::opts_chunk$set(echo=",solution,", include=",results,")
set.seed(1)
```
")

fread <- function(fname){
    readChar(fname, file.info(fname)$size)
}
set.seed(1)
FNAME <- paste0(odir,'extra.Rmd')
cat(preamble,file=FNAME)
if (solution){
    cat("\\newif\\ifsol\n\\soltrue\n",file=FNAME,append=TRUE)
} else {
    cat("\\newif\\ifsol\n\\solfalse\n",file=FNAME,append=TRUE)
}
if (solution & multidoc){
    cat("\\pagestyle{empty}\n",file=FNAME,append=TRUE)
} else {
    cat(fread(paste0(idir,'header.txt')),file=FNAME,append=TRUE)
    cat("# Theory (5 questions per week)\n\n",file=FNAME,append=TRUE)
}
nq <- 21
NN <- round(runif(n=nq,min=0,max=100))
for (q in 1:nq){
    if (solution & multidoc){
        qname <- paste0(odir,'q',q,ifelse(randomise,NN[q],''),'.Rmd')
        print(qname)
        qfile <- cat(fread(FNAME),file=qname)
        cat(q,file=qname,append=TRUE)
        cat("\\. ",file=qname,append=TRUE)
        cat(fread(paste0(idir,"q",q,".txt")),file=qname,append=TRUE)
        rmarkdown::render(qname, 'pdf_document')
    } else {
        cat("1. ",file=FNAME,append=TRUE)
        cat(fread(paste0(idir,"q",q,".txt")),file=FNAME,append=TRUE)
    }
}
if (!solution | !multidoc){
    cat("\n# Practicals (1 problem per week)\n\n",file=FNAME,append=TRUE)
}
for (p in 1:5){
    if (solution & multidoc){
        pname <- paste0(odir,'p',p,ifelse(randomise,NN[p],''),'.Rmd')
        pfile <- cat(fread(FNAME),file=pname)
        cat(p,file=pname,append=TRUE)
        cat("\\. ",file=pname,append=TRUE)
        cat(fread(paste0(idir,"p",p,".txt")),file=pname,append=TRUE)
        rmarkdown::render(pname, 'pdf_document')        
    } else {
        cat("1. ",file=FNAME,append=TRUE)
        cat(fread(paste0(idir,"p",p,".txt")),file=FNAME,append=TRUE)
    }
}
if (!solution | !multidoc){
    rmarkdown::render(FNAME, 'pdf_document')
}
