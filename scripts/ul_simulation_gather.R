ul_simulation_gather <- function(result_path, scale) {
    checkmate::check_choice(
        scale, choices = c("1km", "2km", "4km", "6km", "8km", "10km"))
    
    # Read data
    load(file.path(result_path, sprintf("runs_%s.rda", scale)))
    
    # Get the mean of prediction, marginal_responses, 
    # independent_responses, shap_dependences, and variable_analysis
    ## Prediction
    message("Collect prediction raster.")
    prediction <- merge(do.call(c, lapply(runs, function(run) {
        merge(do.call(c, lapply(run, function (each) {
            each$prediction
        }))) %>% 
            st_apply(c("x", "y"), mean, na.rm = TRUE)
    }))) %>% st_apply(c("x", "y"), mean, na.rm = TRUE)
    
    write_stars(prediction, 
                file.path(result_path, sprintf("ul_suit_space_%s.tif", scale)))
    
    ## Variable analysis
    message("Collect variable analysis.")
    var_analysis <- mean_var_analysis(runs)
    
    ## marginal_responses
    message("Collect marginal responses.")
    marginal_responses <- collect_response(runs, "marginal_responses")
    
    ## independent_responses
    message("Collect independent responses.")
    independent_responses <- collect_response(runs, "independent_responses")
    
    ## shap_dependences
    message("Collect shap dependences.")
    shap_dependences <- collect_response(runs, "shap_dependences")
    
    # return
    list(prediction = prediction,
         variable_analysis = var_analysis,
         marginal_responses = marginal_responses,
         independent_responses = independent_responses,
         shap_dependences = shap_dependences)
}

mean_var_analysis <- function(runs) {
    variables <- runs[[1]][[1]]$variable_analysis$variables
    
    pearson_correlation <- reduce(lapply(runs, function(run) {
        reduce(lapply(run, function(each) {
        each$variable_analysis$pearson_correlation
    }), full_join, by = c("variable", "method", "usage")) %>% 
        mutate(value = rowMeans(.[sapply(., is.numeric)])) %>% 
        select(variable, method, usage, value)
    }), full_join, by = c("variable", "method", "usage")) %>% 
        mutate(value = rowMeans(.[sapply(., is.numeric)])) %>% 
        select(variable, method, usage, value)
    
    full_AUC_ratio <- do.call(rbind, lapply(runs, function(run) {
        do.call(rbind, lapply(run, function(each) {
        each$variable_analysis$full_AUC_ratio
    })) %>% summarise(across(all_of(names(.)), mean))
    })) %>% summarise(across(all_of(names(.)), mean))
    
    AUC_ratio <- reduce(lapply(runs, function(run) {
        reduce(lapply(run, function(each) {
        each$variable_analysis$AUC_ratio
    }), full_join, by = c("variable", "method", "usage")) %>% 
        mutate(value = rowMeans(.[sapply(., is.numeric)])) %>% 
        select(variable, method, usage, value)
    }), full_join, by = c("variable", "method", "usage")) %>% 
        mutate(value = rowMeans(.[sapply(., is.numeric)])) %>% 
        select(variable, method, usage, value)
    
    shap_train <- do.call(rbind, lapply(runs, function(run) {
        do.call(rbind, lapply(run, function(each) {
        each$variable_analysis$SHAP$train %>%
            summarise(across(all_of(names(.)), function(v) mean(abs(v))))
    })) %>% summarise(across(all_of(names(.)), function(v) mean(abs(v))))
    })) %>% summarise(across(all_of(names(.)), function(v) mean(abs(v))))
    
    shap_test <- do.call(rbind, lapply(runs, function(run) {
        do.call(rbind, lapply(run, function(each) {
        each$variable_analysis$SHAP$test %>%
            summarise(across(all_of(names(.)), function(v) mean(abs(v))))
    })) %>% summarise(across(all_of(names(.)), function(v) mean(abs(v))))
    })) %>% summarise(across(all_of(names(.)), function(v) mean(abs(v))))
    
    SHAP <- list(train = shap_train, test = shap_test)
    
    aly <- list(variables = variables,
                pearson_correlation = pearson_correlation,
                full_AUC_ratio = full_AUC_ratio,
                AUC_ratio = AUC_ratio,
                SHAP = SHAP)
    class(aly) <- append("VariableAnalysis", class(aly))
    aly
}

collect_response <- function(runs, type) {
    checkmate::check_choice(
        type, choices = c("marginal_responses", "independent_responses", 
                          "shap_dependences"))
    
    runs <- do.call(c, runs)
    if (type %in% c("marginal_responses", "independent_responses")) {
        nms_cont <- names(runs[[1]][[type]]$responses_cont)
        nms_cat <- names(runs[[1]][[type]]$responses_cat)
        
        if (!is.null(nms_cont)) {
            responses_cont <- lapply(nms_cont, function(var) {
                do.call(rbind, lapply(1:length(runs), function(n) {
                    runs[[n]][[type]]$responses_cont[[var]] %>% 
                        mutate(run = n)
                }))
            })
            names(responses_cont) <- nms_cont
        } else responses_cont <- NULL
        
        if (!is.null(nms_cat)) {
            responses_cat <- lapply(nms_cat, function(var) {
                do.call(rbind, lapply(1:length(runs), function(n) {
                    runs[[n]][[type]]$responses_cat[[var]] %>% 
                        mutate(run = n)
                }))
            })
            names(responses_cat) <- nms_cat
        } else responses_cat <- NULL
        
        list(responses_cont = responses_cont,
             responses_cat = responses_cat)
    } else {
        nms_cont <- names(runs[[1]][[type]]$dependences_cont)
        nms_cat <- names(runs[[1]][[type]]$dependences_cat)
        
        if (!is.null(nms_cont)) {
            dependences_cont <- lapply(nms_cont, function(var) {
                do.call(rbind, lapply(1:length(runs), function(n) {
                    runs[[n]][[type]]$dependences_cont[[var]] %>% 
                        mutate(run = n)
                }))
            })
            names(dependences_cont) <- nms_cont
        } else dependences_cont <- NULL
        
        if (!is.null(nms_cat)) {
            dependences_cat <- lapply(nms_cat, function(var) {
                do.call(rbind, lapply(1:length(runs), function(n) {
                    runs[[n]][[type]]$dependences_cat[[var]] %>% 
                        mutate(run = n)
                }))
            })
            names(dependences_cat) <- nms_cat
        } else dependences_cat <- NULL
        
        feature_values <- do.call(rbind, 
                                   lapply(1:length(runs), 
                                          function(n) {
                runs[[n]][[type]]$feature_values %>% 
                    mutate(run = n)
            }))
        
        list(dependences_cont = dependences_cont,
             dependences_cat = dependences_cat,
             feature_values = feature_values)
    }
}
