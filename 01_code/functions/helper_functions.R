sim_data_path <- function(...){
  
  # Laptop Jens
  if(base::Sys.info()["nodename"] == "OEK-NB-JK01" & 
     base::Sys.info()["effective_user"] == "Jens Klenke"){
    path <- 'C:\\Users\\Jens Klenke\\sciebo - Klenke, Jens (snjeklen@uni-duisburg-essen.de)@uni-duisburg-essen.sciebo.de\\BART_IV-data\\00_sim_data'
  }
  
  # IBES Server
  if(base::Sys.info()["nodename"] == 'IBES-CALC01' & 
     base::Sys.info()["effective_user"] == 'jens.klenke'){
    path <- 'E:\\jens.klenke\\00_sim_data'
  }
  
  cat('\n', 'Input path set to: \n' , 
      '\t -', crayon::yellow(path), '\n' 
      )

  return(path)
}
