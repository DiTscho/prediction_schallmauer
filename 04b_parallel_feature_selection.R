source("exp/nb_04.R")

data_path = "../data/"
dfs = readRDS(paste(data_path,"dfs_filtered.rds"))
dfs = map(dfs, impute_and_clean)


methods =  c("Univariable_cox_filter",
             "RSF_variable_importance", "RSF_minimal_depth", "RSF_variable_hunting",
             "MSRS_RF_variable_importance", "mRMR", "Sequential_selection", "SIS")


feature_selections = readRDS("../data/feature_selections.rds")
saveRDS(feature_selections, file = "../data/feature_selections_backup.rds")


### reduction for gse96058 only
df = dfs[[3]]
df_var = sapply(df[3:ncol(df)], var)
cols = names(df_var[df_var > 0.25])
dfs[[3]] = df[,c("survival_time", "event",cols)]
dim(dfs[[3]])
###

for (i in i:7){
  
  colnames(dfs[[i]]) = gsub("\\s*\\([^\\)]+\\)","", colnames(dfs[[i]]))
  df = dfs[[i]]
  
  cat(paste0('Processing ', names(dfs)[i],'\n'))
  feature_selections[[1,i]] = list(cox_fi(df)                    )
  saveRDS(feature_selections, file = "../data/feature_selections.rds")

  cat(paste0('Processing ', names(dfs)[i], ' 1 of 8 ','\n'))
  feature_selections[[2,i]] = list(rsf_fi(df, method = "vh.vimp"))
  saveRDS(feature_selections, file = "../data/feature_selections.rds")

  cat(paste0('Processing ', names(dfs)[i], ' 2 of 8 ','\n'))
  feature_selections[[3,i]] = list(rsf_fi(df, method = "md")     ) 
  saveRDS(feature_selections, file = "../data/feature_selections.rds")
  
  cat(paste0('Processing ', names(dfs)[i], ' 3 of 8 ','\n'))
  feature_selections[[4,i]] = list(rsf_fi(df, method = "vh")     ) 
  saveRDS(feature_selections, file = "../data/feature_selections.rds")
  
  cat(paste0('Processing ', names(dfs)[i], ' 4 of 8 ','\n'))
  feature_selections[[5,i]] = list(msrs_rf_fi(df)                ) 
  saveRDS(feature_selections, file = "../data/feature_selections.rds")
  
  cat(paste0('Processing ', names(dfs)[i], ' 5 of 8 ','\n'))
  feature_selections[[6,7]] = list(mrmr_fi(df)                   )
  saveRDS(feature_selections, file = "../data/feature_selections.rds")
  
}


