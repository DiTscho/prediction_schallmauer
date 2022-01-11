dimn = function(...){map(list(...), dim)}

check_save = function(data, datafile, datapath = data_path){
  if(!file.exists(paste(datapath, datafile, sep = ""))){
    print(paste("Saving", datafile, "in", datapath))
    save(data, file=paste(datapath, datafile, sep = ""))}
  else{print(paste(datafile, "already exists"))}
  }

gl = compose(get, load) # syntactic sugar

illumina2symbol = function(listOfIlluminaCodes){
  if(!"illuminaHumanv3.db" %in% (.packages())){
    library("illuminaHumanv3.db")
  }
  listOfgeneSymbols = as.character(listOfIlluminaCodes)
  HugoSymbol= data.frame(Gene=unlist(mget(x = listOfIlluminaCodes,envir = illuminaHumanv3SYMBOL)))
  HugoSymbol$Gene
  }

filter_and_annotate_mb = function(df){
  #Get gene names, normalize with limma, and throw away columns without gene names
  df = t(df)
  rnames = illumina2symbol(rownames(df))
  rownames(df) = rnames
  df = limma::avereps(df, rownames(df))
  df = t(df)
  return(as.data.frame(df[,!is.na(colnames(df))]))
  }