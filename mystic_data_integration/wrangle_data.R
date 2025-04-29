transform_numbers_debug <- function(database){
  warning("Slow function (for large datasets). Use 'transform_numbers_speed()' instead for faster performance!")
  c_db <- database
  t_c_db <- apply(c_db, 2, as.numeric)
  t_c_db <- as.data.frame(t_c_db)
  for(i in 1:ncol(t_c_db)){
    if(sum(is.na(t_c_db[,i]))==nrow(t_c_db)){
      t_c_db[,i] <- c_db[,i]
    }
  }
  t_c_db
  
}

transform_numbers_speed <- function(database) {
  warning("This function is stable. However if you are wary of unintended behaviour use 'transform_numbers_debug()' which is somewhat slower but safer!")
  if (!is.data.frame(database)) {
    database <- as.data.frame(database, stringsAsFactors = FALSE)
  }
  is_numeric <- sapply(database, function(x) {
    tryCatch({
      as.numeric(x)
      TRUE
    }, warning = function(w) FALSE)
  })
  database[is_numeric] <- lapply(database[is_numeric], as.numeric)
  return(database)
}

transform_dates <- function(database){
  clinical_db <- database
  c_db_d <- apply(clinical_db, 2, function(x)as.Date(x, '%Y-%m-%d'))
  c_db_d <- as.data.frame(c_db_d)
  for(i in 1:ncol(c_db_d)){
    if(sum(is.na(c_db_d[,i]))!=nrow(c_db_d)){
      c_db_d[,i] <- as.numeric(c_db_d[,i])
      c_db_d[,i] <- as.Date(c_db_d[,i])
    }else{
      c_db_d[,i] <- clinical_db[,i]
    }
  }
  c_db_d
}

transform_dates_speed <- function(database) {
  warning("Use 'transformation_dates()' for a safer option!")
  dt <- as.data.table(database)
  convert_to_date <- function(x) {
    if (is.character(x)) {
      parsed_date <- parse_date_time(x, orders = c("ymd", "dmy", "mdy"), quiet = TRUE)
      if (!all(is.na(parsed_date))) {
        return(as.Date(parsed_date))
      }
    } else if (is.numeric(x)) {
      date_from_numeric <- as.Date(x, origin = "1970-01-01")
      if (!all(is.na(date_from_numeric))) {
        return(date_from_numeric)
      }
    }
    warning("No trasformation was perfomed!")
    return(x)
  }
  dt <- dt[, lapply(.SD, convert_to_date)]
  return(dt)
}

transform_names <- function(database, database_name=c("clinical",
                                                      "clinical_longitudinal",
                                                      "response_longitudinal", 
                                                      "measurements_longitudinal"),
                            data_source=c("sdtm", "adam")){
  tryCatch({
    dictio <-  read_xlsx(paste0("/mnt/data/PIONEER_2025_Historical_data/deid_contents_",data_source,".xlsx"))
  },error=function(e){
    stop(paste0("Make sure the file 'deid_contents_",data_source,".xlsx' exists. Run 'download_dataset_API(dictionary, download=TRUE)' if not!"))
  })
  
  if(!hasArg(database_name)){
    database_name <- deparse(substitute(database))
  }
  if(database_name[1]=="clinical"){
    stop("'clinical' database already has human names. Please, refer to the
            documentation if further information is needed.")
    
  }else if(database_name[1]=="response_longitudinal" | 
           database_name[1]=="measurements_longitudinal" | 
           database_name[1]=="clinical_longitudinal"){
    if(database_name[1]=="response_longitudinal"){
      if(data_source[1]=="sdtm"){
        dictio_rl <- subset(dictio, dictio$memname=="EPU_RS")
      }else if(data_source[1]=="adam"){
        dictio_rl <- subset(dictio, dictio$memname=="EPU_RSRS")
      }else{
        stop("Not a valid data_source. Please review!")
      }
      
    }else if(database_name[1]=="measurements_longitudinal"){
      if(data_source[1]=="sdtm"){
        dictio_rl <- subset(dictio, dictio$memname=="EPU_TR")
      }else if(data_source[1]=="adam"){
        dictio_rl <- subset(dictio, dictio$memname=="EPU_RSTR")
      }else{
        stop("Not a valid data_source. Please review!")
      }
      
    }else if(database_name[1]=="clinical_longitudinal"){
      if(data_source[1]=="sdtm"){
        dictio_rl <- subset(dictio, dictio$memname=="EPU_LB")
      }else if(data_source[1]=="adam"){
        dictio_rl <- subset(dictio, dictio$memname=="EPU_ADLB")
      }else{
        stop("Not a valid data_source. Please review!")
      }
      
    }else{
      stop("Check arguments and database!")
    }
    if(data_source[1]=="adam"){
      colnames(database) <- toupper(colnames(database))
    }
    for(i in 1:ncol(database)){
      tryCatch({
        colnames(database)[i] <- dictio_rl$label[dictio_rl$name==colnames(database)[i]]
      },error=function(e){
        
      })
    }
    colnames(database) <- gsub(" ", "_", colnames(database))
    colnames(database) <- gsub("/", "_", colnames(database))
    colnames(database) <- tolower(colnames(database))
    colnames(database) <- gsub("\n", "_", colnames(database))
    database
  }else{
    stop("Database name not valid. Specify it in the 'database_name' argument or make sure it's a valid database type")
  }
}

generate_assessment_visit_date_dataset <- function(measurements_longitudinal_dataset = measurements_longitudinal,
                                                   response_longitudinal_dataset = response_longitudinal){
  warning("ULTRA SLOW FUNCTION! Intended only for debugging as it uses for loops.
          use 'claude_opt_generate_assessment_visit_date_dplyr' instead!")
  # We get the user ids as unique is study/user, we split by /
  str_splitted <- NULL
  for(i in 1:nrow(measurements_longitudinal)){
    str_splitted[i] <- strsplit(measurements_longitudinal$unique_subject_identifier[i], "/")[[1]][2] 
  }
  
  # we get a patient based table, however not the desired one as we want the total sums, not individual lesions
  pre_avd <- data.frame(studyid=measurements_longitudinal$study_identifier,
                        usubjid = str_splitted,
                        visitnum = measurements_longitudinal$visit_number,
                        ady = measurements_longitudinal$study_day_of_tumor_measurement,
                        day = rep(NA, nrow(measurements_longitudinal)),
                        week = measurements_longitudinal$study_day_of_tumor_measurement/7,
                        mmsumdiam = measurements_longitudinal$numeric_result_finding_in_standard_units,
                        response = rep(NA, nrow(measurements_longitudinal)),
                        tatn = measurements_longitudinal$tumor_assessment_test_name,
                        arf = measurements_longitudinal$accepted_record_flag)
  
  # We get rid of change from baseline / nadir / other and retain only diameters, actual ctScan
  pre_avd2 <- subset(pre_avd, pre_avd$tatn=="Sum of Diameter" & pre_avd$arf=="Y") # Sum of Diameter?
  
  visits <- sort(as.numeric(unique(pre_avd2$ady))) # we'll get the chronological visits
  df <- NULL
  fdf <- NULL
  sum_diam <- NULL
  visit <- NULL
  id <- NULL
  print("GENERATING INTERMEDIATE DATASET")
  for(i in unique(pre_avd2$usubjid)){
    print_progress(i=which(unique(pre_avd2$usubjid==i)), n=len(unique(pre_avd2$usubjid)), initial_iter = 1)
    for(j in 1:len(visits)){
      tdf <- subset(pre_avd2, pre_avd2$usubjid==i)
      sum_diam[j] <- sum(tdf$mmsumdiam[tdf$ady==visits[j]]) # we get the sum of the visit we are interested on 
      visit[j] <- visits[j] # j is an order indicator, we need visits[j] to get the real visit number
      id[j] <- i
    }
    fdf <- data.frame(sum_diam=sum_diam,
                      visit= visit,
                      id=id)
    if(is.null(df)){
      df <- fdf
    }else{
      df <- rbind(df,fdf )
    }
  }
  
  # We need response per visit (aka, for every ctSCAN which response we got)
  df2 <- NULL
  for(i in unique(response_longitudinal$unique_subject_identifier)){ # Patient level
    srl <- subset(response_longitudinal, response_longitudinal$unique_subject_identifier == i)
    ssrl <- srl[, c("response_assessment_result_in_std_format", "study_day_of_response_assessment", "evaluator_identifier")]
    sssrl <- subset(ssrl, ssrl$evaluator_identifier == "") # We get rid of duplicates. Usually Radiologist 1 and 2 don't use the nomenclature were after
    
    if(length(unique(sssrl$study_day_of_response_assessment)) == 0){ # When there are no visits
      fsr <- data.frame(response = NA, visit = NA, id = i)
      df2 <- rbind(df2, fsr)
    } else {
      for(j in unique(sssrl$study_day_of_response_assessment)){ # Visit level
        response_value <- NA
        
        responses <- sssrl$response_assessment_result_in_std_format[sssrl$study_day_of_response_assessment == j]
        
        if(length(responses) > 1) {
          non_Y_responses <- responses[responses != "Y"] # If we have more than one response per visit, means we have 3 nomenclatures, we get rid of the flag "Y"
          
          if(length(non_Y_responses) > 0) { 
            # Pick the first matching response
            matched_responses <- non_Y_responses[non_Y_responses %in% c("PD", "SD", "PR", "CR", "NE")]
            if(length(matched_responses) > 0) {
              response_value <- matched_responses[1]
            } else {
              response_value <- non_Y_responses[1]
            }
          }
        }
        
        # Ensure only single row with consistent lengths
        fsr <- data.frame(response = response_value, visit = j, id = i, stringsAsFactors = FALSE)
        df2 <- rbind(df2, fsr)
      }
    }
  }
  
  
  df$usubjid <- paste0("D419AC00001/",df$id)
  colnames(df) <- c("sumdiam", "ady", "id", "usubjid")
  colnames(df2) <- c("response", "ady", "usubjid")
  mpdf <- dplyr::full_join(df2, df, by=c("usubjid", "ady"))
  
  
  assessment_visit_date <- data.frame(
    studyid = rep("D419AC00001", nrow(mpdf)),
    usubjid = mpdf$id,
    visitnum = rep(NA, nrow(mpdf)),
    ady = mpdf$ady,
    day = rep(NA, nrow(mpdf)),
    week = rep(NA, nrow(mpdf)),
    mmsumdiam = mpdf$sumdiam,
    response = mpdf$response
  )
  print("GENERATING FINAL DATASET")
  for(i in 1:nrow(assessment_visit_date)){
    print_progress(i=i, n=nrow(assessment_visit_date), initial_iter=1)
    assessment_visit_date$visitnum[i] <- pre_avd2$visitnum[pre_avd2$usubjid==assessment_visit_date$usubjid[i] & pre_avd2$ady==assessment_visit_date$ady[i]][1]
    assessment_visit_date$week[i] <- pre_avd2$week[pre_avd2$usubjid==assessment_visit_date$usubjid[i] & pre_avd2$ady==assessment_visit_date$ady[i]][1]
  }
  assessment_visit_date
}

claude_opt_generate_assessment_visit_date_dplyr <- function(measurements_longitudinal_dataset = measurements_longitudinal,
                                                            response_longitudinal_dataset = response_longitudinal) {
  baseline_ids <- measurements_longitudinal_dataset %>%
    group_by(subject_identifier_for_the_study) %>%
    filter(baseline_record_flag=="Y") %>%
    select(date_time_of_tumor_measurement)
  baseline_visits_by_id <- baseline_ids %>%
    group_by(subject_identifier_for_the_study) %>%
    summarise(
      baseline_visit = first(sort(date_time_of_tumor_measurement))
    )

  mldj <- full_join(baseline_visits_by_id, measurements_longitudinal_dataset,
                    by="subject_identifier_for_the_study")
  # 1. Pre-processing measurements_longitudinal
  pre_avd <- mldj %>%
    mutate(
      usubjid = str_split_fixed(unique_subject_identifier, "/", 2)[,2],
      ady2 = as.numeric(date_time_of_tumor_measurement - baseline_visit), 
      week = as.numeric(date_time_of_tumor_measurement - baseline_visit-1) %/% 7 +1
    ) %>%
    select(
      studyid = study_identifier,
      usubjid,
      visitnum = visit_number,
      ady = study_day_of_tumor_measurement,
      ady2 = ady2,
        #study_day_of_tumor_measurement,
      week,
      mmsumdiam = numeric_result_finding_in_standard_units,
      tatn = tumor_assessment_test_name,
      accepted_record_flag = accepted_record_flag
    ) %>%
    filter(tatn == "Sum of Diameter", accepted_record_flag == "Y")
  
  # 2. Aggregate measurements
  df <- pre_avd %>%
    group_by(usubjid, ady) %>%
    summarise(
      sumdiam = first(mmsumdiam),
      visitnum = as.character(first(visitnum)),
      week = first(week),
      ady2 = ady2,
      .groups = "drop"
    )
  
  # 3. Process response_longitudinal
  df2 <- response_longitudinal_dataset %>%
    filter(accepted_record_flag=="Y") %>%
    filter(parameter=="Overall Response") %>%
    mutate(usubjid = str_split_fixed(unique_subject_identifier, "/", 2)[,2]) %>%
    group_by(usubjid, study_day_of_response_assessment) %>%
    summarise(
      response = case_when(
        any(analysis_value_c %in% c("PD", "SD", "PR", "CR", "NE")) ~
          analysis_value_c[analysis_value_c %in% c("PD", "SD", "PR", "CR", "NE")][1],
        any(analysis_value_c != "Y") ~
          analysis_value_c[analysis_value_c != "Y"][1],
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    rename(ady = study_day_of_response_assessment)
  
  # 4. Join datasets
  assessment_visit_date <- df %>%
    full_join(df2, by = c("usubjid", "ady")) %>%
    mutate(
      studyid = "D419AC00001",
      day = NA,
      ady = ady2,
      week = (ady2-1)%/%7+1
    ) %>%
    select(studyid, usubjid, visitnum, ady, day, week, mmsumdiam = sumdiam, response)
  assessment_visit_date$usubjid <- as.factor(assessment_visit_date$usubjid)
  return(assessment_visit_date)
}

generate_patient_data_dataset <- function(clinical_dataset = clinical_db_no_SF,
                                          assessment_visit_data = assessment_visit_data,
                                          measurements_longitudinal_data = measurements_longitudinal_no_SF,
                                          clinical_longitudinal_data=clinical_longitudinal_no_SF){
  first_patient_sd <- sort(clinical_dataset$treatment_start_date)[1] # we get the first visit of the first patient
  week_filtered <- assessment_visit_data %>% ## NOT USED BUT MERGED LATER FIX.
    group_by(usubjid) %>%
    filter(week>1) %>% # week 1 in this dataset is alway baseline measurement for the ctscan
    summarise(
      week = first(week), # we are interested in teh first one after the bl ctscan
      .groups = "drop"
    )
  
  week_filtered2 <- assessment_visit_data %>% 
    group_by(usubjid) %>%
    filter(week>1) %>%
    summarise(
      lweek = last(week), # idem with the last one
      .groups = "drop"
    )
  
  colnames(week_filtered) <- c("subjid", "week")
  colnames(week_filtered2) <- c("subjid", "lweek")
  
  visit_dates <- measurements_longitudinal_data %>% 
    mutate(usubjid = str_split_fixed(unique_subject_identifier, "/", 2)[,2]) %>% 
    group_by(usubjid) %>%
    filter(!is.na(date_of_visit)) %>% 
    arrange(usubjid, date_of_visit) %>%
    summarise(first_date_of_visit=first(date_of_visit),
              last_date_of_visit = last(date_of_visit))

  
  colnames(visit_dates) <- c("subjid", "first_date_of_visit",
                             "last_date_of_visit")      
  
  ldh <- clinical_longitudinal_data %>% 
    group_by(subject_identifier_for_the_study) %>%
    filter(lab_test_or_examination_name=="Lactate Dehydrogenase") %>%
    filter(visit_name=="Screening") %>%
    select(result_or_finding_in_original_units)
  colnames(ldh) <- c("subjid", "ldh")
  # probably don't need to join, however, this way we ensure the order to be 
  # correct when merging both datasets later for the final one
  # fjd = full join dataset
  
  fjd <- full_join(week_filtered, clinical_dataset, by="subjid")  
  fjd <- full_join(fjd, week_filtered2, by="subjid")
  fjd <- full_join(fjd, visit_dates, by="subjid")
  fjd <- full_join(fjd, ldh, by="subjid")
  studyids <- fjd %>% mutate(
    USUBJID = str_split_fixed(USUBJID, "/", 2)[,1]) %>%
    select(USUBJID)
  
  fjdc <- fjd
  colnames(fjdc) <- c("usubjid", colnames(fjdc)[2:100])
  avd2 <- full_join(assessment_visit_data,
                    fjdc[,c("usubjid","progression_free_survival_time")],
                    by="usubjid")
  
  avd2$progression_free_survival_time <- (avd2$progression_free_survival_time-1)%/%7+1
  
  week_filtered3 <- avd2 %>%
    group_by(usubjid) %>%
    filter(week>1 & week<progression_free_survival_time) %>% 
    summarise(
      week = last(week),
      .groups = "drop"
    )
  colnames(week_filtered3) <- c("subjid", "lbpfs")
  fjd <- full_join(fjd,week_filtered3,  by="subjid")
  
  df <- data.frame(studyid= studyids$USUBJID,
                   usubjid=fjd$subjid,
                   trtsdt=fjd$treatment_start_date,
                   trtedt=fjd$treatment_end_date,
                   treatment_end_week=(as.numeric(fjd$treatment_end_date - fjd$treatment_start_date)-1)%/%7+1,
                   treatment_end_day=as.numeric(fjd$treatment_end_date - fjd$treatment_start_date),
                   calendar_week=floor(as.numeric(fjd$treatment_start_date-first_patient_sd)/7)+1, # This one needs the floor format
                   calendar_day=as.numeric(fjd$treatment_start_date-first_patient_sd),
                   patient_min_t=(as.numeric(fjd$first_date_of_visit-fjd$treatment_start_date)-1)%/%7+1,
                   patient_max_t=fjd$lweek, # I don't think we need the floor
                   patient_first_visit=fjd$first_date_of_visit,
                   patient_last_visit=fjd$last_date_of_visit,
                   patient_t_width=(as.numeric(fjd$last_date_of_visit-fjd$first_date_of_visit)-1)%/%7+1,
                   death = ifelse(fjd$overall_survival_censor==1, FALSE, TRUE),
                   death_week = (as.numeric(fjd$death_date-fjd$treatment_start_date)-1)%/%7+1,
                   progression_before_death = ifelse(fjd$progression_free_survival_time==fjd$overall_survival_time, 0, 1),
                   right_censored= ifelse( # Two things to check (regarding death)
                     (((as.numeric(fjd$death_date-fjd$treatment_start_date)-1)%/%7+1 > fjd$lweek) # If death happened AFTER the last treatment week
                      | is.na((as.numeric(fjd$death_date-fjd$treatment_start_date)-1)%/%7+1)) # Or death does not exist (NA)
                     & ((ifelse( # And if two other things to check (regarding PD)
                       ifelse(fjd$progression_free_survival_time==fjd$overall_survival_time, 0, 1)==1, # If PD happens 
                       (fjd$progression_free_survival_time-1)%/%7, NA)> floor(fjd$lweek))| # AFTER the last treatment week
                         is.na(ifelse(ifelse(fjd$progression_free_survival_time==fjd$overall_survival_time, 0, 1)==1, # Or if PD never happens (NA)
                                      (fjd$progression_free_survival_time-1)%/%7+1, NA))),TRUE,FALSE), # Yeah, sorry about this. Basically if the death is after the last treatment week (or there's no death), and if the PD is after the last week of treatment (or no PD) we deem it Right cens.
                   progress_week=ifelse(
                     ifelse(fjd$progression_free_survival_time==fjd$overall_survival_time, 0, 1)==1, # PD happens
                     (fjd$progression_free_survival_time-1)%/%7+1, NA), # We turn the PFS to weeks (PFS as Upper bound)
                   pfs=fjd$lbpfs,
                   interval_censored=ifelse(
                     ifelse(fjd$progression_free_survival_time==fjd$overall_survival_time, 0, 1)==1, # If Patient had a PD
                     (fjd$progression_free_survival_time-1)%/%7+1, NA)-fjd$lbpfs, # PFS (UB) - PFS (LB), otherwise NA (as NA-number is NA)
                   age=fjd$age,
                   age_group=ifelse(fjd$age<18,"<18",ifelse(fjd$age<40, "18-40", ifelse(fjd$age<65,"40-65",ifelse(fjd$age<75,"65-75",">75")))),
                   sex=fjd$sex,
                   country=fjd$country,
                   ecogbl=ifelse(fjd$ecog_score=="ecog performance status 0",0, ifelse(fjd$ecog_score=="ecog performance status 1 or higher", 1, fjd$ecog_score)),
                   bmibl=fjd$bmi,
                   baseline_albumin=fjd$lab_results_albumin_value,
                   baseline_creatinine=fjd$lab_results_creatinine_clearance_value,
                   baseline_hemoglobin=fjd$lab_results_hemoglobin_value,
                   baseline_ldh=fjd$ldh)
  df
}

remove_patient_data_dups <- function(patient_data=patient_data){
  dups <- which(duplicated(patient_data$usubjid))
  leave <- NULL
  j <- 1
  for(i in dups){
    fst <- sum(colSums(is.na(patient_data[(i-1),])))
    dp <- sum(colSums(is.na(patient_data[i,])))
    if(fst>dp){
      leave[j] <- i-1
    }else if(dp>fst){
      leave[j] <- i
    }else{
      leave[j] <- i
    }
   j <- j+1
  }
  patient_data <- patient_data[-leave,]
  #patient_data <- patient_data %>% arrange(usubjid)
  rownames(patient_data) <- 1:nrow(patient_data)
  patient_data
}



calc_visit_date <- function(data) {
  data |> 
    mutate(
      # day = if_else(day > 0, # Is post-treatment day? 
      #               day - 1, 
      #               day),
      week = ((ady-1) %/% 7) + 1, # Last pre-screening week is 0. 
      #treated_week = week > 0, # Was this a post-treatment week?
    ) 
}



# ======================
#        SPARE
# ======================


# resp <- NULL
# df2 <- NULL
# fsr <- NULL

# for(i in unique(response_longitudinal$unique_subject_identifier)){
#   #print(i)
#   # print("entro i")
#   srl <- subset(response_longitudinal, response_longitudinal$unique_subject_identifier==i) #i
#   ssrl <- srl[,c("response_assessment_result_in_std_format","study_day_of_response_assessment","evaluator_identifier")]
#   sssrl <- subset(ssrl, ssrl$evaluator_identifier=="")
#   if(len(unique(sssrl$study_day_of_response_assessment))==0){
#     #print(paste0("holis soc ", i))
#     fsr <- data.frame(response=NA,
#                       visit=NA,
#                       id=i)
#   }else{
#     for(j in unique(sssrl$study_day_of_response_assessment)){
#       #print(j)
#       # print("entro j")
#       tryCatch({
#         if(len(sssrl$"response_assessment_result_in_std_format"[sssrl$study_day_of_response_assessment==j])>1){ # j
#           if(len(which((sssrl$"response_assessment_result_in_std_format"[sssrl$study_day_of_response_assessment==j]!="Y")==TRUE))>1){ # j
#             resp[j] <- sssrl$"response_assessment_result_in_std_format"[sssrl$study_day_of_response_assessment==j][which(sssrl$"response_assessment_result_in_std_format"[sssrl$study_day_of_response_assessment==j] %in% c("PD","SD","PR","CR","NE"))] #j
#           }else{
#             resp[j] <- sssrl$"response_assessment_result_in_std_format"[sssrl$study_day_of_response_assessment==j][which(sssrl$"response_assessment_result_in_std_format"[sssrl$study_day_of_response_assessment==j]!="Y")] #j
#           }
#         }
#       },error=function(e){
#         resp[j] <- sssrl$response_assessment_result_in_std_format[sssrl$study_day_of_response_assessment==j][1]
#       })
#       fsr <- data.frame(response=resp[j],
#                         visit=j,
#                         id=i)
#     }
#   }
#   if(is.null(df2)){
#     df2 <- fsr
#   }else{
#     df2 <- rbind(df2, fsr)
#   }
# }
# df2



# 
# 
# claude_opt_generate_assessment_visit_date_dplyr <- function(measurements_longitudinal_dataset = measurements_longitudinal,
#                                                             response_longitudinal_dataset = response_longitudinal) {
#   
#   # 1. Pre-processing measurements_longitudinal
#   pre_avd <- measurements_longitudinal_dataset %>%
#     mutate(
#       usubjid = str_split_fixed(unique_subject_identifier, "/", 2)[,2],
#       week = study_day_of_tumor_measurement / 7
#     ) %>%
#     select(
#       studyid = study_identifier,
#       usubjid,
#       visitnum = visit_number,
#       ady = study_day_of_tumor_measurement,
#       week,
#       mmsumdiam = numeric_result_finding_in_standard_units,
#       tatn = tumor_assessment_test_name,
#       accepted_record_flag = accepted_record_flag
#     ) %>%
#     filter(tatn == "Sum of Diameter") %>%
#     filter(accepted_record_flag=="Y")
#   
#   # 2. Aggregate measurements
#   df <- pre_avd %>%
#     group_by(usubjid, ady) %>%
#     summarise(
#       sumdiam = mmsumdiam,
#       visitnum = first(visitnum),
#       week = first(week)
#     ) %>%
#     ungroup()
#   
#   # 3. Process response_longitudinal
#   df2 <- response_longitudinal_dataset %>%
#     filter(evaluator_identifier == "") %>%
#     group_by(unique_subject_identifier, study_day_of_response_assessment) %>%
#     summarise(
#       response = case_when(
#         any(response_assessment_result_in_std_format %in% c("PD", "SD", "PR", "CR", "NE")) ~ 
#           response_assessment_result_in_std_format[response_assessment_result_in_std_format %in% c("PD", "SD", "PR", "CR", "NE")][1],
#         any(response_assessment_result_in_std_format != "Y") ~ 
#           response_assessment_result_in_std_format[response_assessment_result_in_std_format != "Y"][1],
#         TRUE ~ NA_character_
#       )
#     ) %>%
#     ungroup() %>%
#     rename(usubjid = unique_subject_identifier, ady = study_day_of_response_assessment)
#   
#   # 4. Join datasets
#   assessment_visit_date <- df %>%
#     full_join(df2, by = c("usubjid", "ady")) %>%
#     mutate(
#       studyid = "D419AC00001",
#       day = NA
#     ) %>%
#     select(studyid, usubjid, visitnum, ady, day, week, mmsumdiam = sumdiam, response)
#   
#   return(assessment_visit_date)
# }
