# Load libraries
options(warn=-1) 

library(data.table)
library(tidyverse)
library(lubridate)
library(reticulate)
library(dplyr)
library(tidyr)
library(reshape2)

use_python('/Users/quicl/Anaconda3')
np   <- import('numpy', convert = FALSE)
pd   <- import('pandas', convert = FALSE)
js   <- import('json', convert = FALSE)
re   <- import('re', convert = FALSE)

# Read Json data 

load_json_event <- function(df = "berlin_ards.json", window = 15, back = 0, use_diff = FALSE ) {
  
  py <- import_builtins()
  with(py$open("berlin_ards.json") %as% json_file, {
    df = js$load(json_file) 
  })
  
  # Convert py to r
  df = pd$DataFrame(df)$transpose()
  df = py_to_r(df)

  # Load oraganised json data
  source('load_json_data.R')
  info_timex <- load_json_data()
  
  print(paste0('Same feature names across each subject : ',(length(unique(df[[2]]))==1)))
  matrix_names   <- unlist(df[[2]][1])
  
  
  # find a way to detect length in each vital signs
  a <- b <- mearsured_points  <- d  <- e              <- c()
  ok_id  <- info_timex$subject_id
  id <- info_timex$subject_id
  strt <- Sys.time()
  
  for (i in 1:length(id)){ #length(id)
    if (id[i] %in% ok_id){
      i_idx <- info_timex$idx[info_timex$subject_id==id[i]]#info$idx[i]
      #print(i_idx)
      for (k in 1:length(matrix_names)){
        pos_k <- length(df[[4]][[i_idx]][[k]])
        a     <- append(a,pos_k)
      }
      if(sum(a)%%length(matrix_names) != 0 | a[3] == 0){
        #print(paste0('a sum ', sum(a)))
        b <- append(b, id[i])
        e <- append(e, i_idx)
      }  # get the time point length of each subject  
      if (a[5] != 0 & (!i_idx %in% d)) { mearsured_points <- append(mearsured_points, a[1])
      d <- append(d, i_idx)}  
      # get the row id
      
      a <- c()
    }
  }
  
  remove_incomplete_idx <- setdiff(d,info_timex$idx)
  
  end <- Sys.time()
  (end-strt) 
  
  print(paste0('No need to remove extra ids ; ', length(remove_incomplete_idx)<=0))
  
  ###########################  ###########################  ###########################
  
  # select feature data from time_to_event window
  
  compare_tp <- data.frame(cbind(mearsured_points,d)) %>%
    mutate(idx = d) %>%
    right_join(info_timex[,c('idx', 'total_timpoints','prior_to_event','is_Berlin_ARDS')])%>%
    mutate(reject = total_timpoints - mearsured_points )%>%
    mutate(death = ifelse(is_Berlin_ARDS==TRUE, 1, 0) )%>%
    select(-c(idx,is_Berlin_ARDS))%>%
    filter(prior_to_event > (window+1))

  
  # Decide window for prediction, add column names 
  time_to_events <- window
  time_b4_events <- ifelse(window > back, back, 0)
  recon_feature <- data.frame(matrix(NA,nrow(compare_tp)*(time_to_events-time_b4_events),9)) #sum(compare_tp$prior_to_event)
  colnames(recon_feature) <- unlist(df[[2]][1])
  
  row_nx          <- c()
  b               <- a <- c <- 0
  diff_change     <- FALSE
  for (x in compare_tp$d){
    
    lst1            <- lapply(df[[4]][[x]], unlist)
    c               <- c+1
    feature_mx      <- as.data.frame(lst1)[(compare_tp$prior_to_event[c]-time_to_events+1):(compare_tp$prior_to_event[c]-time_b4_events),]
    b               <- a+1
    a               <- a+nrow(feature_mx)
                  
    recon_feature[b:a,]   <- feature_mx
    row_nx          <- append(row_nx,rep(x,nrow(feature_mx)))
    
  }
  recon_feature$row_idx <- row_nx
  end <- Sys.time()
  (end-strt)
  
  ##########################  ###########################  ###########################
  # feature engineer
  Develop <- TRUE
  
  if(Develop){
    na_abp  <- which(is.na(recon_feature$SysABP))
    na_WBC  <- which(is.na(recon_feature$WBC))
    na_SpO2 <- which(is.na(recon_feature$SpO2))
    na_HR   <- which(is.na(recon_feature$HR))
    na_Creatinine   <- which(is.na(recon_feature$Creatinine))
    na_Platelets    <- which(is.na(recon_feature$Platelets))
    na_Temp         <- which(is.na(recon_feature$Temp))
    na_RespRate     <- which(is.na(recon_feature$RespRate))
    
  }else{
    print(matrix_names)
   
    hist(a$SysABP[!is.na(a$SysABP)]) 
    
    which(is.na(a$DiasABP))[1:5]
    hist(a$DiasABP[!is.na(a$DiasABP)])
    
    setdiff(which(is.na(a$SysABP)),which(is.na(a$DiasABP)))
    
    # correlation between 2 apbs
    cor_dis <- a$DiasABP[!is.na(a$DiasABP)][-2246]
    cor_sys <- a$SysABP[!is.na(a$SysABP)]#[-2246]
    cor(cor_dis,cor_sys) 
    
    length(which(is.na(a$SpO2)))#[1:5]
    hist(a$WBC[!is.na(a$WBC)])
    
  }

 dvar <- function(value1,idx = row_idx){
     df <- data.frame(cbind(idx, value1))
     colnames(df) <- c("row_idx","value1")
     
     df %>%
     group_by(row_idx) %>%
     mutate(valuex = as.numeric(sub("NaN", NA, value1)))%>%
     fill(valuex, .direction = "updown")%>%
     mutate(diff_value =  c(0,diff(valuex, lag = 1))) -> dfx
     
     dfx$diff_value
       }

  a <- recon_feature %>% 
       mutate(correct_creatinine = ifelse(Creatinine < 10, Creatinine*88.4017, Creatinine))%>%
       mutate(Creatinine = correct_creatinine)%>%
       mutate(class_Creatinine = case_when(correct_creatinine < 120 & correct_creatinine > 52 ~ 0,
                                           correct_creatinine >= 150 ~ 2,
                                           correct_creatinine >= 120 ~ 1,
                                           correct_creatinine <= 32  ~ -2,
                                           correct_creatinine <= 52  ~ -1,
                                           
                                            T ~ -99))%>% 

       mutate(class_ABP = case_when(SysABP < 135 & SysABP > 85 ~ 0,
                                    SysABP >= 155 ~ 2,
                                    SysABP >= 135 ~ 1,
                                    SysABP <= 65  ~ -2,
                                    SysABP <= 85  ~ -1,
                                    T ~ -99))%>%
    
      mutate(class_Temp = case_when(Temp < 36.8 & Temp > 36.5 ~ 0,
                                 Temp >= 38 ~ 2,
                                 Temp >= 36.8 ~ 1,
                                 Temp <= 36  ~ -2,
                                 Temp <= 36.5  ~ -1,
                                 T ~ -99))%>%
    
      mutate(class_RespRate= case_when(RespRate < 17 & RespRate > 20 ~ 0,
                                 RespRate >= 25 ~ 2,
                                 RespRate >= 20 ~ 1,
                                 RespRate <= 14  ~ -2,
                                 RespRate <= 17  ~ -1,
                                 T ~ -99))%>%

      mutate(class_HR = case_when(HR < 100 & HR > 60 ~ 0,
                                   HR >= 120 ~ 2,
                                   HR >= 100 ~ 1,
                                   HR <= 45 ~ -2,
                                   HR <= 60  ~ -1,
                                   T ~ -99))%>% 

       mutate(class_SpO2 = case_when(SpO2 < 100 & SpO2 > 95 ~ 0,
                                     SpO2 <= 80  ~ -3,
                                     SpO2 <= 85  ~ -2,
                                     SpO2 <= 95  ~ -1,
                                     T ~ -99))%>% 
    
       mutate(class_WBC = case_when(WBC < 10 & WBC > 4.5 ~ 0,
                                    WBC <= 3  ~ -2,
                                    WBC <= 4.5  ~ -1,
                                    WBC >= 15 ~ 2,
                                    WBC >= 10 ~ 1,
                                      T ~ -99))%>%

      mutate(class_Platelets =  case_when(Platelets < 400 & Platelets > 150 ~ 0,
                                          Platelets >= 620 ~ 2,
                                          Platelets >= 400 ~ 1,
                                          Platelets <= 100 ~ -2,
                                          Platelets <= 150  ~ -1,
                                          T ~ -99))%>% 
       mutate(dCreatinine = dvar(Creatinine,idx = row_idx))%>%
       mutate(dSysABP = dvar(SysABP,idx = row_idx))%>%
       mutate(dHR = dvar(HR,idx = row_idx))%>%
       mutate(dSpO2 = dvar(SpO2,idx = row_idx))%>%
       mutate(dWBC = dvar(WBC,idx = row_idx))%>%
       mutate(dPlatelets = dvar(Platelets,idx = row_idx))%>%
       mutate(dRespRate = dvar(RespRate,idx = row_idx))%>%
       mutate(dTemp = dvar(Temp,idx = row_idx))%>%
       select(-c(correct_creatinine, DiasABP ))
  
  a$class_Creatinine[na_Creatinine] <- -99
  a$class_ABP[na_abp] <- -99
  a$class_SpO2[na_SpO2] <- -99
  a$class_HR[na_HR] <- -99
  a$class_WBC[na_WBC] <- -99
  a$class_Platelets [na_Platelets] <- -99
  a$class_RespRate [na_RespRate] <- -99
  a$class_Temp [na_Temp] <- -99

  
  ##########################  ###########################  ###########################
  # Apply carry-forward imputation
  
  data <- function(df, diff_change = FALSE){
    if(diff_change){
      
      # fill NA with next value
      df %>%
        group_by(row_idx, variable)%>%
        mutate(value1 = as.numeric(sub("NaN", NA, value)))%>%
        fill(value1, .direction = "updown")%>%
        mutate(diff_value =  c(0,diff(value1, lag = 1)))
      #%>%        dplyr::ungroup()

         
    }
    else{
        df %>%
        group_by(row_idx, variable)%>%
        mutate(value1 = sub("NaN", NA, value))%>%
        fill(value1, .direction = "updown")
        }
  }
  
  # convert long to wide format, add a time point number to each feature
  
  strt <- Sys.time()
  
  diff_changex = use_diff 

  
  xxx_ori_ck <- gather(recon_feature, variable, value, -row_idx) 
  xxx_ori    <- gather(a, variable, value, -row_idx) 
  
  xxx_ori$number <-  rep(1:time_to_events,(nrow(xxx_ori)/(time_to_events-time_b4_events))) 
  xxx_ori_ck$number <- rep(1:time_to_events,(nrow(xxx_ori_ck)/(time_to_events-time_b4_events)))
  
  ##########################  ###########################  ###########################
  
  # check data
  rz <- sample(xxx_ori$row_idx ,1)
  ry <- sample(1:9,1)
  t_2_event <- info_timex$prior_to_event[info_timex$idx==rz]
  ori_a   <- df[[4]][[rz]][[ry]][t_2_event-time_b4_events]
  recon_a <- xxx_ori_ck$value[xxx_ori_ck$row_idx == rz & xxx_ori_ck$number == time_to_events & xxx_ori_ck$variable == matrix_names[ry]]
  
  print(paste0('PASS Sanity CHECK,  processed data is equal to the original json : ', ori_a == recon_a))
  print(paste0("Total subject id used in the analysis: ",length(unique(xxx_ori$row_idx))))
  print(paste0("Record trace back to : ",0.5*window, "hours, and before ", time_b4_events*0.5, "hours to the event"))

  dt <- data(xxx_ori,diff_change = use_diff)

  dt
  }
