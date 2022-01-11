suppressMessages(source("exp/nb_04.R"))

### LOAD ###########


dfs = readRDS("../data/dfs_filtered_for_feature_selection.rds")

#export
learning_algs =   c("Cox_PH_model",
                    "Ridge",
                    "Elastic_Net",
                    "Lasso",
                    "Gradient_Boosting_tree_based",
                    "Gradient_Boosting_linear_model_based",
                    "Random_Survival_Forests",
                    "Maximally_selected_rank_statistics_Random_Forests",
                    "Survival_Tree"
)

#export
learners = list(makeLearner("surv.coxph",           id = learning_algs[[1]]),
                makeLearner("surv.cvglmnet",        id = learning_algs[[2]], alpha = 0,   nfolds=20),
                makeLearner("surv.cvglmnet",        id = learning_algs[[3]], alpha = 0.5, nfolds=20, s="lambda.min"),
                makeLearner("surv.cvglmnet",        id = learning_algs[[4]], alpha = 1,   nfolds=20, s="lambda.min"),
                makeLearner("surv.gamboost",        id = learning_algs[[5]], baselearner = "bols" ),
                makeLearner("surv.gamboost",        id = learning_algs[[6]], baselearner = "btree"),
                makeLearner("surv.randomForestSRC", id = learning_algs[[7]]),
                makeLearner("surv.ranger",          id = learning_algs[[8]]),
                makeLearner("surv.rpart",           id = learning_algs[[9]])
)
names(learners) = learning_algs

#export
metalearner = function(df, learner, feature_selector){
  task  = makeSurvTask(data = df, target = c("survival_time", "event"))
  inner = makeResampleDesc("CV", iters=5)
  n = 20  # number of features
  feature_selectors = c("univariate_model_score", "mrmr", "randomForestSRC_importance",
                        "randomForestSRC_var_select_md", "randomForestSRC_var_select_vh", 
                        "party_cforest_importance")
  
  if(!(feature_selector %in% feature_selectors)){
    stop("feature_selector must be one of ", feature_selectors)   
  }
  
  if (feature_selector == "univariate_model_score"){
    lrn = makeFilterWrapper(learner = learner, fw.method="univariate.model.score",    fw.abs = n, 
                            perf.learner=learner)
  }
  else if (feature_selector == "mrmr"){
    lrn = makeFilterWrapper(learner = learner, fw.method="mrmr",                        fw.abs = n)
  }
  else if (feature_selector == "randomForestSRC_importance"){
    lrn  = makeFilterWrapper(learner = learner, fw.method="randomForestSRC_importance", fw.abs = n)
  }
  else if (feature_selector == "randomForestSRC_var_select_md"){
    lrn  = makeFilterWrapper(learner = learner, fw.method="randomForestSRC_var.select", fw.abs = n, 
                             more.args = list("randomForestSRC_var.select"=list(method="md")))
  }
  else if (feature_selector == "randomForestSRC_var_select_vh"){
    lrn  = makeFilterWrapper(learner = learner, fw.method="randomForestSRC_var.select", fw.abs = n, 
                             more.args = list("randomForestSRC_var.select"=list(method="vh")))
  }
  else if (feature_selector == "party_cforest_importance"){
    lrn  = makeFilterWrapper(learner = learner, fw.method="party_cforest.importance",   fw.abs = n)
  }
  
  res   = resample(learner = lrn, task = task, resampling=inner, models=TRUE, show.info  = FALSE)
  return(res)
}


feature_selectors = c("univariate_model_score", 
                      #"mrmr", 
                      "randomForestSRC_importance",
                      #"randomForestSRC_var_select_md", 
                      "randomForestSRC_var_select_vh", 
                      "party_cforest_importance")

#"randomForestSRC_var_select_md"

data_path = "../data/metalearners/metalearners_from_feature_selections/"

### summary ############

#names(dfs)
#"METABRIC" "GSE11121" "GSE96058" "GSE7390"  "GSE9893"  "NKI"      "TCGA"     "GSE4922" 
#> names(learners)
#[1] "Cox_PH_model"                                     
#[2] "Ridge"                                            
#[3] "Elastic_Net"                                      
#[4] "Lasso"                                            
#[5] "Gradient_Boosting_tree_based"                     
#[6] "Gradient_Boosting_linear_model_based"             
#[7] "Random_Survival_Forests"                          
#[8] "Maximally_selected_rank_statistics_Random_Forests"
#[9] "Survival_Tree"             

#> feature_selectors
#[1] "univariate_model_score"        "mrmr"                          "randomForestSRC_importance"   
#[4] "randomForestSRC_var_select_md" "randomForestSRC_var_select_vh" "party_cforest_importance" 


### RUN ################

not_working = c(paste(data_path,"ml_", names(dfs)[[1]],"_",names(learners)[[6]], "_", "party_cforest_importance.rds", sep=""))



for (d in 2:length(dfs)){                          # dataframes
  print(names(dfs)[[d]])
  for (l in 1:length(learners)){                 # learners
    print(names(learners)[[l]])
    for (f in 2:length(feature_selectors)){    # feature selections except univariate_model_score
      
      filename = paste(data_path,"ml_",
                       names(dfs)[[d]],"_",
                       names(learners)[[l]], "_",
                       feature_selectors[[f]],
                       ".rds", sep="")
      
      if (file.exists(filename) | filename %in% not_working ){next}
      else {
        print(feature_selectors[[f]])
        df               = dfs[[d]]
        learner          = learners[[l]]
        feature_selector = feature_selectors[[f]]
        
        res = metalearner(df, learner, feature_selector)
        
        saveRDS(res, filename)
        
      }
    }
  }
}

### tmp ###


for (d in 1:length(dfs)){                          # dataframes
  for (l in 1:length(learners)){                 # learners
    for (f in 2:length(feature_selectors)){    # feature selections except univariate_model_score
      
      filename = paste(data_path,"ml_",
                       names(dfs)[[d]],"_",
                       names(learners)[[l]], "_",
                       feature_selectors[[f]],
                       ".rds", sep="")
      print(filename)

}
}}










