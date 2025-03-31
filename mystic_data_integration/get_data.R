
get_dataset <- function(dataset=c("clinical", "clinical_longitudinal", 
                                  "response_longitudinal",
                                  "measurements_longitudinal", "dictionary"), 
                        data_source=c("sdtm","adam")){
  
  # It takes the object path, retrieves it, gets the ID, and uses that to
  # retrieve the actual data. Probably not needed and just getting the ID as in
  # quartzbio would be enough, but IDs seem frail and its not taking much more
  # time anyway
  
  if(dataset[1]=="clinical"){
    
    clinical_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/Clinical/clinical.D419AC00001.20241216_150504.csv")
    clinical_csv <- solvebio::Object.query(id=clinical_db$id, paginate=TRUE)
    clinical_csv
  }else if(dataset[1]=="clinical_longitudinal"){
    if(data_source[1]=="sdtm"){
      clinical_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/epu_lb.csv")
    }else if(data_source[1]=="adam"){
      clinical_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/epu_adlb.csv")
    }else{
      warning("Wrong data_source arg!")
      break()
    }
    
    clinical_long_csv <- Object.query(id=clinical_long_db$id, paginate=TRUE)
    clinical_long_csv
  }else if(dataset[1]=="response_longitudinal"){
    if(data_source[1]=="sdtm"){
      response_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/epu_rs.csv")
    }else if(data_source[1]=="adam"){
      response_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/epu_rsrs.csv")
    }else{
      warning("Wrong data_source arg!")
      break()
    }
    
    response_long_csv <- Object.query(id=response_long_db$id, paginate=TRUE)
    response_long_csv
  }else if(dataset[1]=="measurements_longitudinal"){
    if(data_source[1]=="sdtm"){
      measu_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/epu_tr.csv")
    }else if(data_source[1]=="adam"){
      measu_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/epu_rstr.csv")
    }else{
      warning("Wrong data_source arg!")
      break()
    }
    
    measu_long_csv <- Object.query(id=measu_long_db$id, paginate=TRUE)
    measu_long_csv
  }else if(dataset[1]=="dictionary"){
    
  }else{
    print("Error! use help_get_dataset() to see valid arguments!")
    break()
  }
}


help_get_dataset <- function(information=c("general", "dataset_type")){
  if(information[1]=="dataset_type"){
    cat(paste("'clinical': Information about clinical data, 1 row per patient, baseline measurements, outcome information, etc.",
              "'clinical_longitudinal': Information about clinical data at a longitudinal level. Multiple rows per patient (LDH, albumin, ...). TAKES A WHILE TO RUN",
              "'response_longitudinal': Information about response in the CTscan by RECIST 1.1. Multiple rows per patient, summ of lesions NOT included.",
              "'measurements_longitudinal': Information about lesions for each ctSCAN and pacient. Multiple rows per patient. TAKES A LONG TIME!",
              "'dictionary': A list of human names for some databases code column names", sep="\n"))
  }else if(information[1]=="general"){
    cat("Use 'download_dataset_API()' to get the datasets (use download=TRUE if its the first time you use this function on Domino).",
        "'get_dataset()' might fail due to reading the csv directly from quartzbio.", sep="\n")
    
  }else{
    warning("Information requirement not known. Please check the argument!")
    break()
  }
  
}

download_dataset_API <- function(dataset=c("clinical", "clinical_longitudinal", 
                                           "response_longitudinal",
                                           "measurements_longitudinal", "dictionary"),
                                 download=FALSE,
                                 data_source=c("sdtm", "adam")){
  if(dataset[1]=="clinical"){
    if(download==TRUE){
      warning("'clinical' might give unexpected errors when 'download=TRUE', load it locally if that is the case!")
      clinical_db <- Dataset.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/Clinical/clinical.D419AC00001.20241216_150504")
      clinical_csv <- DatasetExport.create(
        clinical_db$id,
        format = 'csv',
        params = NULL,
        send_email_on_completion = FALSE
      )
      url <- DatasetExport.get_download_url(clinical_db$id)
      url <- DatasetExport.get_download_url(clinical_csv$id)
      download.file(url, '/mnt/data/PIONEER_2025_Historical_data/clinical_data.csv')
    }
    read.csv("/mnt/data/PIONEER_2025_Historical_data/clinical_data.csv")
    
  }else if(dataset[1]=="clinical_longitudinal"){
    if(download==TRUE){
      if(data_source[1]=="sdtm"){
        clinical_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/epu_lb.csv")
      }else if(data_source[1]=="adam"){
        clinical_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/epu_adlb.csv")
      }else{
        warning("Error, data_source arg not valid!")
        break()
      }
      
      url <- solvebio::Object.get_download_url(clinical_long_db$id)
      download.file(url, paste0('/mnt/data/PIONEER_2025_Historical_data/clinical_long_data_',data_source[1],'.csv'))
    }
    read.csv(paste0("/mnt/data/PIONEER_2025_Historical_data/clinical_long_data_",data_source[1],".csv"))
  }else if(dataset[1]=="response_longitudinal"){
    if(download==TRUE){
      if(data_source[1]=="sdtm"){
        response_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/epu_rs.csv")
      }else if(data_source[1]=="adam"){
        response_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/epu_rsrs.csv")
      }else{
        warning("Error, data_source arg not valid!")
        break()
      }
      
      url <- Object.get_download_url(response_long_db$id)
      download.file(url, paste0('/mnt/data/PIONEER_2025_Historical_data/response_long_data_',data_source[1],'.csv'))
    }
    read.csv(paste0("/mnt/data/PIONEER_2025_Historical_data/response_long_data_",data_source[1],".csv"))
  }else if(dataset[1]=="measurements_longitudinal"){
    if(download==TRUE){
      if(data_source[1]=="sdtm"){
        measu_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/epu_tr.csv")
      }else if(data_source[1]=="adam"){
        measu_long_db <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/epu_rstr.csv")
      }else{
        warning("Error, data_source arg not valid!")
        break()
      }
      url <- Object.get_download_url(measu_long_db$id)
      download.file(url, paste0('/mnt/data/PIONEER_2025_Historical_data/measurements_long_data_',data_source[1],'.csv'))
    }
    read.csv(paste0("/mnt/data/PIONEER_2025_Historical_data/measurements_long_data_",data_source[1],".csv"))
  }else if(dataset[1]=="dictionary"){
    if(download==TRUE){
      if(data_source[1]=="sdtm"){
        dict <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/sdtm/cdap_20220908/docs/deid_contents_sdtm.xlsx")
      }else if(data_source[1]=="adam"){
        dict <- Object.get_by_full_path("astrazeneca:_INT-20230324-1259:/HISTORICAL_D419AC00001/entimice_source_files/adam/cdap_20220908/docs/deid_contents_adam.xlsx")
      }else{
        warning("Error, data_source arg not valid!")
        break()
      }
      
      url <- Object.get_download_url(dict$id)
      download.file(url, paste0('/mnt/data/PIONEER_2025_Historical_data/deid_contents_',data_source[1],'.xlsx'))
    }
    read_excel(paste0("/mnt/data/PIONEER_2025_Historical_data/deid_contents_",data_source[1],".xlsx"))
  }else{
    print("Error! use help_get_dataset() to see valid arguments!")
    break()
  }
}
