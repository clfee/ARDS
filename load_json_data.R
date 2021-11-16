# Load libraries

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

load_json_data <- function(df = "berlin_ards.json" ) {
  
  py <- import_builtins()
  with(py$open("berlin_ards.json") %as% json_file, {
    df = js$load(json_file) 
  })
    
    # Convert py to r
    df = pd$DataFrame(df)$transpose()
    df = py_to_r(df)
    
    # Data loading 
    lst <- lapply(df[[1]], unlist)
    info     <- as.data.frame(lst)
    
    id = as.numeric(sub("(?i).*subject_id.*?(\\d+).*", "\\1", colnames(info))) # make sure ids are the same
    
    info$key <- rownames(info)
    
    info <- reshape2::melt(info, id = "key")
    info <- reshape2::dcast(info, formula = variable ~ key)%>%
      select(-c(variable))%>%
      mutate(first_measurement_time = as_datetime(sub("T", " ", first_measurement_time)))%>%
      mutate(is_Berlin_time = as_datetime(sub("T", " ", is_Berlin_time)))%>%
      mutate(last_measurement_time = as_datetime(sub("T", " ", last_measurement_time)))%>%
      mutate(prior_to_event = floor(as.numeric(difftime(is_Berlin_time,first_measurement_time, units = "mins"))/30))%>%
      mutate(after_event = ceiling(as.numeric(difftime(last_measurement_time,is_Berlin_time, units = "mins"))/30))%>%
      mutate(T.prior_to_event = as.numeric(difftime(is_Berlin_time,first_measurement_time, units = "mins")))%>%
      mutate(T.after_event =    as.numeric(difftime(last_measurement_time,is_Berlin_time, units = "mins")))%>%
      mutate(total_timpoints = ceiling(as.numeric(difftime(last_measurement_time,first_measurement_time, units = "mins"))/30))%>%
      mutate(check_prior = as_datetime(prior_to_event*30 + as.numeric(first_measurement_time)))%>%
      mutate(prior_ok = ifelse(check_prior <= is_Berlin_time, 0, 1))%>%
      mutate(check_timpoints = difftime(last_measurement_time,first_measurement_time, units = "mins"))%>%
      mutate(age_group = cut_interval(as.numeric(age), n = 5))%>%
      select(-c(first_measurement_time,last_measurement_time,is_in_hospital_death,T.after_event,check_timpoints))
    
    info$idx <- c(1:nrow(info))
    
    info_timex <- info %>%
      filter(T.prior_to_event > 0) %>%
      filter(!duplicated(subject_id))%>%
      filter(!subject_id %in% c(88146,31290))
    
    print(paste0(' PASS Sanity CHECK, No reduandant id : ',all.equal(unique(id),as.numeric(unique((info$subject_id))))))
    print(paste0(' PASS Sanity CHECK, All Piror Event are measured : ',sum(info_timex$prior)==0))
    
    if (all.equal(unique(id),as.numeric(unique((info$subject_id))))) {
              info_timex } else {print('Check Data Quality 2')}
}

