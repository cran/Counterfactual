
counterfactual <- function (formula, data, weights, na.action = na.exclude, group, treatment = FALSE, decomposition = FALSE, counterfactual_var, transformation = FALSE, quantiles = c(1:9)/10, method = "qr", trimming = 0.005, nreg = 100, scale_variable, counterfactual_scale_variable, censoring = 0, right = FALSE, nsteps = 3, firstc = 0.1, secondc = 0.05, noboot = FALSE, weightedboot = FALSE, seed = 8, robust = FALSE, reps = 100, alpha = 0.05, first = 0.1, last = 0.9, cons_test = 0, printdeco = TRUE, sepcore = FALSE, ncore = 1){

    taus = quantiles
    nreg = round(nreg)
    if (nreg < 1) { 
        stop("The option nreg must be a strictly positive integer.")
    }
    if (missing(group) && missing(counterfactual_var)) {
        stop("One of the options group and counterfactual are required")
    }

    if (method != "qr" && method != "loc" && method != "locsca" && method != "cqr" && method != "cox" && method != "logit" && method != "probit" && method != "lpm"){
        stop("The selected method has not been implemented")
    }
    if (method == "locsca" && !missing(scale_variable) && missing(group) && missing(counterfactual_scale_variable)){
        stop("If the location scale estimator is used with the option scale, then either the group option or the counterfactual_scale_variable option must be provided.")
    }
    if (method == "cqr"){
        if (nsteps < 3)         stop("The options nsteps must be at least 3")
        if (missing(censoring)) stop("The option censoring must be provided to use the censored quantile regression estimator.")
    }
    if (trimming > 0.1){
        stop("Trimming must be less than 10% of the observations.") 
    }     

    mf = match.call()
    if (!missing(group)){
        m                     = match(c("formula", "data", "weights", "na.action", "group"), names(mf), 0)
        mf                    = mf[c(1, m)]
        mf$drop.unused.levels = TRUE
        mf[[1]]               = as.name("model.frame")
        mf                    = eval.parent(mf)
        mt                    = attr(mf, "terms")
        dep                   = model.response(mf, "numeric")
        reg                   = model.matrix(mt, mf)
        reg                   = as.matrix(reg[, -1])
        weight                = model.weights(mf)

        if(!is.null(weight) && !is.numeric(weight)){ stop("'weights' must be a numeric vector")}
        if(!is.null(weight)){weight = weight} else {weight = rep(1, length(dep))}
        weight                = weight/sum(weight)

        groupN                = mf[, ncol(mf)]
        group0                = (groupN == 0)
        group1                = (groupN == 1)

        dep0                  = subset(dep, group0)
        reg0                  = subset(reg, group0)
        dep1                  = subset(dep, group1)
        counterfactual_var    = subset(reg, group1)
        weight0               = subset(weight, group0)
        weight1               = subset(weight, group1)
        obs0                  = length(dep0)
        obs1                  = length(dep) - obs0

    } else {
        m                     = match(c("formula", "data", "weights", "na.action","counterfactual_var"), names(mf), 0)
        mf                    = mf[c(1, m)]
        mf$drop.unused.levels = TRUE
        mf[[1]]               = as.name("model.frame")
        mf                    = eval.parent(mf)
        mt                    = attr(mf, "terms")
        dep                   = model.response(mf, "numeric")
        reg                   = model.matrix(mt, mf)
        reg                   = as.matrix(reg[, -1])
        weight                = model.weights(mf)

        if (!is.null(weight) && !is.numeric(weight)) { stop("'weights' must be a numeric vector")}
        if (!is.null(weight)) { weight = weight} else { weight = rep(1, length(dep))}
        weight                = weight/sum(weight)

        dep0                  = dep
        reg0                  = reg
        dep1                  = dep0
        counterfactual_var    = mf[, ncol(mf)]
        weight0               = weight
        weight1               = weight
        obs0                  = length(dep)
        obs1                  = obs0
    }    

    depeval = sort(unique(dep0))
    if (!missing(nreg) && (nreg != -1) && (nreg < length(depeval))) {
        wq         = (1:(nreg - 2) - 0.5)/(nreg - 2)
        seqU       = floor(wq * length(depeval))
        thredIndex = c(1, seqU, length(depeval))
        depeval    = depeval[thredIndex]
    }

    if (treatment) {
        depeval_1 = sort(unique(dep1))
        if (!missing(nreg) && (nreg != -1) && (nreg < length(depeval_1))) {
            wq           = (1:(nreg - 2) - 0.5)/(nreg - 2)
            seqU_1       = floor(wq * length(depeval_1))
            thredIndex_1 = c(1, seqU_1, length(depeval_1))
            depeval_1    = depeval_1[thredIndex_1]
        }
    } else {
        depeval_1 = depeval
    }

    #################################################### Estimation Begins ##########################################   
    if (missing(scale_variable)) {scale_variable = reg0}
    if (missing(counterfactual_scale_variable)) {counterfactual_scale_variable = counterfactual_var}
    QE_Results = QteDistEst(dep0, depeval, reg0, weight0, counterfactual_var, dep1, depeval_1, treatment, decomposition, weight1, taus, method, trimming, nreg, scale_variable, counterfactual_scale_variable, censoring, right, nsteps, firstc, secondc)
    #################################################################################################################           
    quantile_effect               = QE_Results$quantile_effect
    marginal_obs                  = QE_Results$marginal_obs
    marginal_fitted               = QE_Results$marginal_fitted
    marginal_counterfactual       = QE_Results$marginal_counterfactual 
    if (treatment && decomposition) {
        composition_effect  = QE_Results$composition_effect
        total_effect        = QE_Results$total_effect

        qte_cov_obs         = c(quantile_effect, marginal_obs, marginal_fitted, marginal_counterfactual, composition_effect, total_effect, QE_Results$marginal_obs_1, QE_Results$marginal_fitted_1, QE_Results$F_ms)
        if(method=="logit" | method=="probit" | method=="lpm"){
            qte_cov_obs     <- c(qte_cov_obs, QE_Results$marginal_fitted_cdf, QE_Results$cdf_obs, QE_Results$marginal_counterfactual_cdf, QE_Results$marginal_fitted_1_cdf, QE_Results$cdf_obs1)
        }
    } else {
        qte_cov_obs         = c(quantile_effect, marginal_obs, marginal_fitted, marginal_counterfactual)
        if(method=="logit" | method=="probit" | method=="lpm"){
            qte_cov_obs     <- c(qte_cov_obs, QE_Results$marginal_fitted_cdf, QE_Results$cdf_obs, QE_Results$marginal_counterfactual_cdf)
        }
    }

    # Inference & Testing 
    if (!noboot) {

        data0             = list(dep = dep0, reg = reg0, weight = weight0, scale_variable = scale_variable)
        data1             = list(counterfactual_var = counterfactual_var, weight = weight1, scale_variable = counterfactual_scale_variable)
        
        qte_cov_booti     = BootstrapProcedure(transformation, reps, data0, depeval, data1, weightedboot, obs0, obs1, dep1, depeval_1, treatment, decomposition, taus, method, trimming, nreg, censoring, right, nsteps, firstc, secondc, seed, sepcore, ncore)
        Bootstrap_Results = InferenceTestingEval(method, qte_cov_booti, reps, taus, robust, qte_cov_obs, last, first, alpha, cons_test, treatment, decomposition, depeval, depeval_1)
      
    } else {

        qte           = quantile_effect
        distributions = cbind(marginal_obs, marginal_fitted, marginal_counterfactual)
        if (treatment && decomposition) {
            composition_effect = composition_effect
            total_effect       = total_effect
        }
    }

    # Display Confidence Interval
    if (printdeco) {
        cat("\n")
        cat(format("Conditional Model: ", width = 40))
        if (method == "qr") { cat("linear quantile regression\n") }
        if (method == "loc") { cat("location model\n")}
        if (method == "locsca") { cat("location scale model\n") }
        if (method == "cqr") { cat("linear censored quantile regression\n") }
        
        if (method == "logit") { cat("logit\n") }
        if (method == "probit") { cat("probit\n") }
        if (method == "lpm") { cat("linear probability model\n") }
        if (method == "cox") { cat("cox duration model\n") }
        if (method != "cox") { cat(format("Number of regressions estimated:", width = 40), QE_Results$nreg, "\n\n")}

        if (!noboot) {

            cat("The variance has been estimated by bootstraping the results", reps, "times.\n\n")
            cat(format("No. of obs. in the reference group:", width = 40), obs0, "\n")
            cat(format("No. of obs. in the counterfactual group:", width = 40), obs1, "\n\n\n")

            if (!treatment && !decomposition) {

                ## CE
                bootstrap_CI_result_print_ce = Bootstrap_Results$resCE

                cat(format("Quantile Effects -- Composition", width = 70, justify = "centre"), "\n")
                cat(paste(rep("-", 70), collapse = ""), "\n")
                cat(format("", width = 20), format("Pointwise", width = 8, justify = "centre"), format("Pointwise", width = 20, justify = "centre"), format("Functional", width = 18, justify = "centre"), "\n")
                if (robust) {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Robust S.E.", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                } else {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Std.Err", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                }
                for (i in 1:length(taus)) {
                    cat(format(taus[i], width = 8, justify = "left"), formatC(bootstrap_CI_result_print_ce[i, ], digits = 3, width = 10), "\n", sep = "")
                }
                if (printdeco) {

                    Counter_Test_CE = Bootstrap_Results$testCE
                    row_of_test     = nrow(Counter_Test_CE)
                    constest1       = sort(unique(cons_test))
                    constestNo0     = setdiff(constest1, 0)
                    nc              = length(constestNo0)

                    cat("\n\n")
                    cat(format("Bootstrap inference on the counterfactual quantile process", width = 70, justify = "centre"), "\n")
                    cat(paste(rep("-", 70), collapse = ""), "\n")
                    cat("\t\t\t\t\t\t    P-values\t\n")
                    cat("\t\t\t\t\t       ", paste(rep("=", 22), collapse = ""), "\n")
                    cat(format("NULL-Hypthoesis", width = 38), "\t", format("KS-statistic", width = 12), format("CMS-statistic", width = 12), "\n")
                    cat(paste(rep("=", 70), collapse = ""), "\n")
                    cat(format("Correct specification of the parametric model", width = 40), format(Counter_Test_CE[1, 1], digits = 2, nsmall = 2, width = 10, justify = "centre"), format(Counter_Test_CE[1, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    cat(format("No effect: QE(tau)=0 for all taus            ", width = 40), format(Counter_Test_CE[2, 1], digits = 2, nsmall = 2, width = 10, justify = "centre"), format(Counter_Test_CE[2, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    if (nc > 0) {
                        for (i in 1:nc) {
                            cat(format(paste("Constant effect: QE(tau)=", constestNo0[i], "for all taus"), width = 40), format(Counter_Test_CE[i + 2, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_CE[i + 2, 2], width = 10, digits = 2, nsmall = 2), "\n")
                        }
                    }
                    cat(format("Constant effect: QE(tau)=QE(0.5) for all taus", width = 40), format(Counter_Test_CE[row_of_test - 2, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_CE[row_of_test - 2, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)>0 for all taus ", width = 40), format(Counter_Test_CE[row_of_test - 1, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_CE[row_of_test - 1, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)<0 for all taus ", width = 40), format(Counter_Test_CE[row_of_test, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_CE[row_of_test, 2], width = 10, digits = 2, nsmall = 2), "\n")
                }

            } else if (treatment && !decomposition) {

                bootstrap_CI_result_print_se = Bootstrap_Results$resSE

                cat(format("Quantile Effects -- Structure", width = 70, justify = "centre"), "\n")
                cat(paste(rep("-", 70), collapse = ""), "\n")
                cat(format("", width = 20), format("Pointwise", width = 8, justify = "centre"), format("Pointwise", width = 20, justify = "centre"), format("Functional", width = 18, justify = "centre"), "\n")
                if (robust) {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Robust S.E.", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                } else {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Std.Err", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                }
                for (i in 1:length(taus)) {
                    cat(format(taus[i], width = 8, justify = "left"), formatC(bootstrap_CI_result_print_se[i, ], digits = 3, width = 10), "\n", sep = "")
                }
                if (printdeco) {

                    Counter_Test_SE = Bootstrap_Results$testSE
                    row_of_test     = nrow(Counter_Test_SE)
                    constest1       = sort(unique(cons_test))
                    constestNo0     = setdiff(constest1, 0)
                    nc              = length(constestNo0)

                    cat("\n\n")
                    cat(format("Bootstrap inference on the counterfactual quantile process", width = 70, justify = "centre"), "\n")
                    cat(paste(rep("-", 70), collapse = ""), "\n")
                    cat("\t\t\t\t\t\t    P-values\t\n")
                    cat("\t\t\t\t\t       ", paste(rep("=", 22), collapse = ""), "\n")
                    cat(format("NULL-Hypthoesis", width = 38), "\t", format("KS-statistic", width = 12), format("CMS-statistic", width = 12), "\n")
                    cat(paste(rep("=", 70), collapse = ""), "\n")
                    cat(format("Correct specification of the parametric model", width = 40), format(Counter_Test_SE[1, 1], digits = 2, nsmall = 2, width = 10, justify = "centre"), format(Counter_Test_SE[1, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    cat(format("No effect: QE(tau)=0 for all taus            ", width = 40), format(Counter_Test_SE[2, 1], digits = 2, nsmall = 2, width = 10, justify = "centre"), format(Counter_Test_SE[2, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    if (nc > 0) {
                        for (i in 1:nc) {
                        cat(format(paste("Constant effect: QE(tau)=", constestNo0[i], "for all taus"), width = 40), format(Counter_Test_SE[i + 2, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_SE[i + 2, 2], width = 10, digits = 2, nsmall = 2), "\n")
                        }
                    }
                    cat(format("Constant effect: QE(tau)=QE(0.5) for all taus", width = 40), format(Counter_Test_SE[row_of_test - 2, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_SE[row_of_test - 2, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)>0 for all taus ", width = 40), format(Counter_Test_SE[row_of_test - 1, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_SE[row_of_test - 1, 2], width = 10, digits = 2, nsmall = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)<0 for all taus ", width = 40), format(Counter_Test_SE[row_of_test, 1], digits = 2, nsmall = 2, width = 10), format(Counter_Test_SE[row_of_test, 2], width = 10, digits = 2, nsmall = 2), "\n")
                }

            } else if (treatment && decomposition) {

                ##SE
                bootstrap_CI_result_print_se = Bootstrap_Results$resSE

                cat(format("Quantile Effects -- Structure", width = 70, justify = "centre"), "\n")
                cat(paste(rep("-", 70), collapse = ""), "\n")
                cat(format("", width = 20), format("Pointwise", width = 8, justify = "centre"), format("Pointwise", width = 20, justify = "centre"), format("Functional", width = 18, justify = "centre"), "\n")
                if (robust) {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Robust S.E.", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                } else {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Std.Err", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                }
                for (i in 1:length(taus)) {
                    cat(format(taus[i], width = 8, justify = "left"), formatC(bootstrap_CI_result_print_se[i, ], digits = 3, width = 10), "\n", sep = "")
                }

                if (printdeco) {

                    Counter_Test_SE = Bootstrap_Results$testSE
                    row_of_test     = nrow(Counter_Test_SE)
                    constest1       = sort(unique(cons_test))
                    constestNo0     = setdiff(constest1, 0)
                    nc              = length(constestNo0)

                    cat("\n\n")
                    cat(format("Bootstrap inference on the counterfactual quantile process", width = 70, justify = "centre"), "\n")
                    cat(paste(rep("-", 70), collapse = ""), "\n")
                    cat("\t\t\t\t\t\t    P-values\t\n")
                    cat("\t\t\t\t\t       ", paste(rep("=", 22), collapse = ""), "\n")
                    cat(format("NULL-Hypthoesis", width = 38),  "\t", format("KS-statistic", width = 12), format("CMS-statistic", width = 12), "\n")
                    cat(paste(rep("=", 70), collapse = ""), "\n")
                    cat(format("Correct specification of the parametric model", width = 40), format(Counter_Test_SE[1, 1], digits = 2, width = 10, justify = "centre"), format(Counter_Test_SE[1, 2], width = 10, digits = 2), "\n")
                    cat(format("No effect: QE(tau)=0 for all taus            ", width = 40), format(Counter_Test_SE[2, 1], digits = 2, width = 10, justify = "centre"), format(Counter_Test_SE[2, 2], width = 10, digits = 2), "\n")
                    if (nc > 0) {
                        for (i in 1:nc) {
                            cat(format(paste("Constant effect: QE(tau)=", constestNo0[i], "for all taus"), width = 40), format(Counter_Test_SE[i + 2, 1], digits = 2, width = 10), format(Counter_Test_SE[i + 2, 2], width = 10, digits = 2), "\n")
                        }
                    }
                    cat(format("Constant effect: QE(tau)=QE(0.5) for all taus", width = 40), format(Counter_Test_SE[row_of_test - 2, 1], digits = 2, width = 10), format(Counter_Test_SE[row_of_test - 2, 2], width = 10, digits = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)>0 for all taus ", width = 40), format(Counter_Test_SE[row_of_test - 1, 1], digits = 2, width = 10), format(Counter_Test_SE[row_of_test - 1, 2], width = 10, digits = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)<0 for all taus ", width = 40), format(Counter_Test_SE[row_of_test, 1], digits = 2, width = 10), format(Counter_Test_SE[row_of_test, 2], width = 10, digits = 2), "\n")
                }

                ## CE
                bootstrap_CI_result_print_ce = Bootstrap_Results$resCE
                cat("\n\n")
                cat(format("Quantile Effects -- Composition", width = 70, justify = "centre"), "\n")
                cat(paste(rep("-", 70), collapse = ""), "\n")
                cat(format("", width = 20), format("Pointwise", width = 8, justify = "centre"), format("Pointwise", width = 20, justify = "centre"), format("Functional", width = 18, justify = "centre"), "\n")
                if (robust) {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Robust S.E.", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                } else {
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Std.Err", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                }
                for (i in 1:length(taus)) {
                    cat(format(taus[i], width = 8, justify = "left"), formatC(bootstrap_CI_result_print_ce[i, ], digits = 3, width = 10), "\n", sep = "")
                }
                if (printdeco) {

                    Counter_Test_CE = Bootstrap_Results$testCE
                    cat("\n\n")
                    cat(format("Bootstrap inference on the counterfactual quantile process", width = 70, justify = "centre"), "\n")
                    cat(paste(rep("-", 70), collapse = ""), "\n")
                    cat("\t\t\t\t\t\t    P-values\t\n")
                    cat("\t\t\t\t\t       ", paste(rep("=", 22), collapse = ""), "\n")
                    cat(format("NULL-Hypthoesis", width = 38), "\t", format("KS-statistic", width = 12), format("CMS-statistic", width = 12), "\n")
                    cat(paste(rep("=", 70), collapse = ""), "\n")
                    cat(format("Correct specification of the parametric model", width = 40), format(Counter_Test_CE[1, 1], digits = 2, width = 10, justify = "centre"), format(Counter_Test_CE[1, 2], width = 10, digits = 2), "\n")
                    cat(format("No effect: QE(tau)=0 for all taus            ", width = 40), format(Counter_Test_CE[2, 1], digits = 2, width = 10, justify = "centre"), format(Counter_Test_CE[2, 2], width = 10, digits = 2), "\n")
                    if (nc > 0) {
                        for (i in 1:nc) {
                            cat(format(paste("Constant effect: QE(tau)=", constestNo0[i], "for all taus"), width = 40), format(Counter_Test_CE[i + 2, 1], digits = 2, width = 10), format(Counter_Test_CE[i + 2, 2], width = 10, digits = 2), "\n")
                        }
                    }
                    cat(format("Constant effect: QE(tau)=QE(0.5) for all taus", width = 40), format(Counter_Test_CE[row_of_test - 2, 1], digits = 2, width = 10), format(Counter_Test_CE[row_of_test - 2, 2], width = 10, digits = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)>0 for all taus ", width = 40), format(Counter_Test_CE[row_of_test - 1, 1], digits = 2, width = 10), format(Counter_Test_CE[row_of_test - 1, 2], width = 10, digits = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)<0 for all taus ", width = 40), format(Counter_Test_CE[row_of_test, 1], digits = 2, width = 10), format(Counter_Test_CE[row_of_test, 2], width = 10, digits = 2), "\n")
                }

                ## TE
                bootstrap_CI_result_print_te = Bootstrap_Results$resTE
                cat("\n\n")
                cat(format("Quantile Effects -- Total", width = 70, justify = "centre"), "\n")
                cat(paste(rep("-", 70), collapse = ""), "\n")
                cat(format("", width = 20), format("Pointwise", width = 8, justify = "centre"), format("Pointwise", width = 20, justify = "centre"), format("Functional", width = 18, justify = "centre"), "\n")
                if(robust){
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Robust S.E.", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                }else{
                    cat(format("Quantile", width = 10, justify = "centre"), format("Est.", width = 10, justify = "centre"), format("Std.Err", width = 8, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 20, justify = "centre"), format(paste((1 - alpha) * 100, "% Conf.Interval", sep = ""), width = 18, justify = "centre"), "\n")
                }
                for (i in 1:length(taus)) {
                    cat(format(taus[i], width = 8, justify = "left"), formatC(bootstrap_CI_result_print_te[i, ], digits = 3, width = 10), "\n", sep = "")
                }
                if (printdeco) {

                    Counter_Test_TE = Bootstrap_Results$testTE
                    cat("\n\n")
                    cat(format("Bootstrap inference on the counterfactual quantile process", width = 70, justify = "centre"), "\n")
                    cat(paste(rep("-", 70), collapse = ""), "\n")
                    cat("\t\t\t\t\t\t    P-values\t\n")
                    cat("\t\t\t\t\t       ", paste(rep("=", 22), collapse = ""), "\n")
                    cat(format("NULL-Hypthoesis", width = 38), "\t", format("KS-statistic", width = 12), format("CMS-statistic", width = 12), "\n")
                    cat(paste(rep("=", 70), collapse = ""), "\n")
                    cat(format("Correct specification of the parametric model", width = 40), format(Counter_Test_TE[1, 1], digits = 2, width = 10, justify = "centre"), format(Counter_Test_TE[1, 2], width = 10, digits = 2), "\n")
                    cat(format("No effect: QE(tau)=0 for all taus            ", width = 40), format(Counter_Test_TE[2, 1], digits = 2, width = 10, justify = "centre"), format(Counter_Test_TE[2, 2], width = 10, digits = 2), "\n")
                    if (nc > 0) {
                        for (i in 1:nc) {
                            cat(format(paste("Constant effect: QE(tau)=", constestNo0[i], "for all taus"), width = 40), format(Counter_Test_TE[i + 2, 1], digits = 2, width = 10), format(Counter_Test_TE[i + 2, 2], width = 10, digits = 2), "\n")
                        }
                    }
                    cat(format("Constant effect: QE(tau)=QE(0.5) for all taus", width = 40), format(Counter_Test_TE[row_of_test - 2, 1], digits = 2, width = 10), format(Counter_Test_TE[row_of_test - 2, 2], width = 10, digits = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)>0 for all taus ", width = 40), format(Counter_Test_TE[row_of_test - 1, 1], digits = 2, width = 10), format(Counter_Test_TE[row_of_test - 1, 2], width = 10, digits = 2), "\n")
                    cat(format("Stochastic dominance: QE(tau)<0 for all taus ", width = 40), format(Counter_Test_TE[row_of_test, 1], digits = 2, width = 10), format(Counter_Test_TE[row_of_test, 2], width = 10, digits = 2), "\n")
                }
            }

        }else {
            cat("The variance has not been computed.\n")
            cat("Do not turn the option noboot on if you want to compute it.\n\n")
            cat(format("No. of obs. in the reference group:", width = 40), obs0, "\n")
            cat(format("No. of obs. in the counterfactual group:", width = 40), obs1, "\n\n\n")

            if (treatment) {
                cat(format("Quantile Effects -- Structure", width = 50, justify = "centre"), "\n")
            } else {
                cat(format("Quantile Effects -- Composition", width = 50, justify = "centre"), "\n")
            }
            cat(paste(rep("-", 50), collapse = ""), "\n")
            cat(format("Quantile", width = 18, justify = "right"), format("Est.", width = 20, justify = "centre"),  "\n")
            for (i in 1:length(quantile_effect)) {
                cat(format(taus[i], width = 15, justify = "centre"), "\t", format(quantile_effect[i], justify = "centre", width = 15, digits = 2), "\n", sep = "")
            }

            if (treatment && decomposition) {
                cat("\n\n")
                cat(format("Quantile Effects -- Composition", width = 50, justify = "centre"), "\n")
                cat(paste(rep("-", 50), collapse = ""), "\n")
                cat(format("Quantile", width = 18, justify = "right"), format("Est.", width = 20, justify = "centre"), "\n")
                for (i in 1:length(composition_effect)) {
                    cat(format(taus[i], width = 15, justify = "centre"), "\t", format(composition_effect[i], justify = "centre", width = 15, digits = 2), "\n", sep = "")
                }
                cat("\n\n")
                cat(format("Quantile Effects -- Total", width = 50, justify = "centre"), "\n")
                cat(paste(rep("-", 50), collapse = ""), "\n")
                cat(format("Quantile", width = 18, justify = "right"), format("Est.", width = 20, justify = "centre"), "\n")
                for (i in 1:length(total_effect)) {
                    cat(format(taus[i], width = 15, justify = "centre"), "\t", format(total_effect[i], justify = "centre", width = 15, digits = 2), "\n", sep = "")
                }
            }
        }
    }

    ## returns
    results = NULL
    if (treatment && decomposition) {
        results = list(quantiles = taus, structral_effect = quantile_effect, composition_effect = composition_effect, total_effect = total_effect, nreg = nreg, marginal_counterfactual = marginal_counterfactual)
    }else if (treatment && !decomposition) {
        results = list(quantiles = taus, structral_effect = quantile_effect, nreg = nreg, marginal_counterfactual = marginal_counterfactual)
    }else {
        results = list(quantiles = taus, composition_effect = quantile_effect, nreg = nreg, marginal_counterfactual = marginal_counterfactual)
    }

    if(method=="logit" | method=="probit" | method=="lpm"){
        results$depeval <- depeval
        if(treatment) {
            results$depeval_1 <- depeval_1
        }
    }
    if (!noboot) {
        results        = c(results, Bootstrap_Results)
        class(results) <- c("counterfactual", class(results))
    }
    return(results)
}

QteDistEst <- function (dep0, depeval, reg0, weight0, counterfactual_var, dep1, depeval_1, treatment, decomposition, weight1, taus, method, trimming, nreg, scale_variable, counterfactual_scale_variable, censoring, right, nsteps, firstc, secondc) {
    
    reg1 = counterfactual_var

    if (method == "qr") {

        wq                          = (1:nreg - 0.5)/nreg
        wqNew                       = wq[(wq >= trimming) && wq <= (1 - trimming)]

        wqW                         = rep(1, length(wqNew))
        wei0                        = as.matrix(weight0) %*% wqW
        wei1                        = as.matrix(weight1) %*% wqW

        conditional_coef            = coef(rq(dep0 ~ reg0, tau = wqNew, weights = weight0))
        marginal_fitted             = wtd.quantile(as.vector(trimming + cbind(1, reg0) %*% conditional_coef), as.vector(wei0), probs = taus, na.rm = TRUE)
        marginal_counterfactual     = wtd.quantile(as.vector(trimming + cbind(1, counterfactual_var) %*% conditional_coef), as.vector(wei1), probs = taus, na.rm = TRUE)
        
        if (treatment) {
            conditional_coef_1      = coef(rq(dep1 ~ reg1, tau = wqNew, weights = weight1))
            marginal_fitted_1       = wtd.quantile(as.vector(trimming + cbind(1, reg1) %*% conditional_coef_1), as.vector(wei1), probs = taus, na.rm = TRUE)
        }

        if (treatment && decomposition) {
            y_thred                 = c(depeval, depeval_1)
            y_te                    = c(dep0, dep1)

            conditionTaus_0         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse((trimming + cbind(1, reg0) %*% conditional_coef) <= xxx, 1, 0))}))
            conditionTaus_1         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse((trimming + cbind(1, reg1) %*% conditional_coef_1) <= xxx, 1, 0))}))

            F_0hat                  = length(dep0)/(length(dep0) + length(dep1))
            F_est                   = conditionTaus_0 * F_0hat + conditionTaus_1 * (1 - F_0hat)
            F_obs                   = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))
            F_ms                    = F_obs - F_est
        }
    }

    if (method == "loc") {

        wq                          = (1:nreg - 0.5)/nreg

        wqW                         = rep(1, length(wq))
        wei0                        = as.matrix(weight0) %*% wqW
        wei1                        = as.matrix(weight1) %*% wqW

        loc                         = lm(dep0 ~ reg0, weights = weight0)
        conditional_coef            = kronecker(as.matrix(coef(loc)), matrix(1, 1, nreg))
        conditional_coef[1, ]       = conditional_coef[1, ] + wtd.quantile(as.matrix(resid(loc)), weight0, wq, na.rm = TRUE)
        marginal_fitted             = wtd.quantile(as.vector(cbind(1, reg0) %*% conditional_coef), as.vector(wei0), taus, na.rm = TRUE)
        marginal_counterfactual     = wtd.quantile(as.vector(cbind(1, counterfactual_var) %*% conditional_coef), as.vector(wei1), taus, na.rm = TRUE)

        if (treatment) {
            loc_1                   = lm(dep1 ~ reg1, weights = weight1)
            conditional_coef_1      = kronecker(as.matrix(coef(loc_1)), matrix(1, 1, nreg))
            conditional_coef_1[1, ] = conditional_coef_1[1, ] + wtd.quantile(as.matrix(resid(loc_1)), weight1, wq, na.rm = TRUE)
            marginal_fitted_1       = wtd.quantile(as.vector(cbind(1, reg1) %*% conditional_coef_1), as.vector(wei1), taus, na.rm = TRUE)
        }

        if (treatment && decomposition) {

            y_thred                 = c(depeval, depeval_1)
            y_te                    = c(dep0, dep1)

            conditionTaus_0         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse((cbind(1, reg0) %*% conditional_coef) <= xxx, 1, 0))}))
            conditionTaus_1         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse((cbind(1, reg1) %*% conditional_coef_1) <= xxx, 1, 0))}))

            F_0hat                  = length(dep0)/(length(dep0) + length(dep1))
            F_est                   = conditionTaus_0 * F_0hat + conditionTaus_1 * (1 - F_0hat)
            F_obs                   = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))
            F_ms                    = F_obs - F_est
        }
    }

    if (method == "locsca") {

        wq                          = (1:nreg - 0.5)/nreg

        wqW                         = rep(1, length(wq))
        wei0                        = as.matrix(weight0) %*% wqW
        wei1                        = as.matrix(weight1) %*% wqW

        loc                         = lm(dep0 ~ reg0, weights = weight0)
        conditional_resid           = resid(loc)
        lmR                         = lm(log(conditional_resid^2) ~ scale_variable, weights = weight0)
        predsca                     = sqrt(exp(predict(lmR)))
        residN                      = wtd.quantile(conditional_resid/predsca, weight0, wq, na.rm = TRUE)
        predall                     = kronecker(predict(loc), matrix(1, 1, nreg)) + kronecker(predsca, t(residN))
        marginal_fitted             = wtd.quantile(as.vector(predall), wei0, probs = taus, na.rm = TRUE)

        predCounter                 = cbind(1, counterfactual_var) %*% coef(loc)
        predscaCounter              = sqrt(exp(cbind(1, counterfactual_scale_variable) %*% coef(lmR)))
        predallCounter              = kronecker(predCounter, matrix(1, 1, nreg)) + kronecker(predscaCounter, t(residN))
        marginal_counterfactual     = wtd.quantile(as.vector(predallCounter), wei1, probs = taus, na.rm = TRUE)
        
        if (treatment) {
            loc_1                   = lm(dep1 ~ reg1, weights = weight1)
            conditional_resid_1     = resid(loc_1)
            lmR_1                   = lm(log(conditional_resid_1^2) ~ counterfactual_scale_variable, weights = weight1)
            predsca_1               = sqrt(exp(predict(lmR_1)))
            residN_1                = wtd.quantile(conditional_resid_1/predsca_1, weight1, wq, na.rm = TRUE)
            predall_1               = kronecker(predict(loc_1), matrix(1, 1, nreg)) + kronecker(predsca_1, t(residN_1))
            marginal_fitted_1       = wtd.quantile(as.vector(predall_1), wei1, probs = taus, na.rm = TRUE)
        }

        if (treatment && decomposition) {

            y_thred                 = c(depeval, depeval_1)
            y_te                    = c(dep0, dep1)

            conditionTaus_0         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse(predall <= xxx, 1, 0))}))
            conditionTaus_1         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse(predall_1 <= xxx, 1, 0))}))

            F_0hat                  = length(dep0)/(length(dep0) + length(dep1))
            F_est                   = conditionTaus_0 * F_0hat + conditionTaus_1 * (1 - F_0hat)
            F_obs                   = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))
            F_ms                    = F_obs - F_est
        }
    }

    if (method == "cqr") {
        
        wq                          = (1:nreg - 0.5)/nreg

        est_cqrl <- function(dep, reg, weight, wq, censoring, nsteps, c1, c2) {

            coef_logit              = coef(glm((dep > censoring) ~ reg, weights = weight, family = binomial(link = logit)))
            pred                    = plogis(cbind(1, reg) %*% coef_logit)
            
            coef_wq                 = NULL
            for (i in wq) {
                delta1              = quantile(pred[pred > 1 - i], c1)
                select1             = (pred >= delta1)
                temp                = coef(rq(dep[select1] ~ as.matrix(reg)[select1, ], i, weights = weight[select1]))

                steps   = 3
                while (steps <= nsteps) {
                    pred1           = cbind(1, reg) %*% temp
                    delta2          = quantile((pred1 - censoring)[pred1 > censoring], c2)
                    select2         = (pred1 >= delta2)
                    temp            = coef(rq(dep[select2] ~ as.matrix(reg)[select2, ], i, weights = weight[select2]))
                    steps           = steps + 1
                    select1         = select2
                }
                coef_wq             = cbind(coef_wq, temp)
            }
            return(coef_wq)
        }
        if (!right) {
            conditional_coef        = est_cqrl(dep0, reg0, weight0, wq, censoring, nsteps, firstc, secondc)
        } else {
            conditional_coef        = -est_cqrl(-dep0, reg0, weight0, 1 - wq, -censoring, nsteps, firstc, secondc)
        }
        marginal_fitted             = quantile(as.vector(cbind(1, reg0) %*% conditional_coef), probs = taus, na.rm = TRUE)
        marginal_counterfactual     = quantile(as.vector(cbind(1, counterfactual_var) %*% conditional_coef), probs = taus, na.rm = TRUE)

        if (treatment) {
            if (!right) {
                conditional_coef_1  = est_cqrl(dep1, reg1, weight1, wq, censoring, nsteps, firstc, secondc)
            } else {
                depR_1              = -dep1
                censoringR_1        = -censoring
                conditional_coef_1  = -est_cqrl(depR_1, reg1, weight1, 1 - wq, censoringR_1, nsteps, firstc, secondc)
            }
            marginal_fitted_1       = quantile(as.vector(cbind(1, reg1) %*% conditional_coef_1), probs = taus, na.rm = TRUE)
        }

        if (treatment && decomposition) {

            y_thred                 = c(depeval, depeval_1)
            y_te                    = c(dep0, dep1)

            conditionTaus_0         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse((cbind(1, reg0) %*% conditional_coef) <= xxx, 1, 0)) }))
            conditionTaus_1         = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {rowMeans(ifelse((cbind(1, reg1) %*% conditional_coef_1) <= xxx, 1, 0))}))

            F_0hat                  = length(dep0)/(length(dep0) + length(dep1))
            F_est                   = conditionTaus_0 * F_0hat + conditionTaus_1 * (1 - F_0hat)
            F_obs                   = colMeans(apply(as.matrix(y_thred), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))
            F_ms                    = F_obs - F_est
        }
    }

    if (method == "cox" | method == "logit" | method == "probit" | method == "lpm") {

        if (method == "cox") {
            if (min(dep0) < 0) {
                stop("Only positive dependent variable allowed")
            }

            fit0                    = coxph(Surv(dep0) ~ reg0, weights = weight0, ties = "breslow")
            depeval                 = basehaz(fit0)$time
            coef0                   = coef(fit0)
            S0                      = (survfit(fit0)$surv)^c(exp(-colMeans(reg0) %*% coef0))
            NS0                     = length(S0)

            reg0                    = as.matrix(reg0)
            counterfactual_var      = as.matrix(counterfactual_var)

            index0                  = exp(reg0 %*% coef0)
            index1                  = exp(counterfactual_var %*% coef0)

            NReg0                   = nrow(reg0)
            Sc0                     = kronecker(matrix(S0, 1, NS0), matrix(1, NReg0, 1))^kronecker(matrix(index0, NReg0, 1), matrix(1, 1, NS0))
            F0                      = sort(apply(as.matrix(1-Sc0) , 2, function(a) wtd.mean(a, weights = weight0)))
            marginal_fitted         = getquantile(depeval, F0, taus)

            NReg1                   = nrow(counterfactual_var)
            Sc1                     = kronecker(matrix(S0, 1, NS0), matrix(1, NReg1, 1))^kronecker(matrix(index1, NReg1, 1), matrix(1, 1, NS0))
            F1                      = sort(apply(as.matrix(1-Sc1), 2, function(a) wtd.mean(a, weights = weight1)))
            marginal_counterfactual = getquantile(depeval, F1, taus)
            
            if (treatment) {
                fit_1               = coxph(Surv(dep1) ~ reg1, weights = weight1, ties = "breslow")
                depeval_1           = basehaz(fit_1)$time
                coef_1              = coef(fit_1)
                reg1                = as.matrix(reg1)
                S_1                 = (survfit(fit_1)$surv)^c(exp(-colMeans(reg1) %*% coef_1))
                NS_1                = length(S_1)

                index_1             = exp(reg1 %*% coef(fit_1))

                NReg_1              = nrow(reg1)
                Sc0_1               = kronecker(matrix(S_1, 1, NS_1), matrix(1, NReg_1, 1))^kronecker(matrix(index_1, NReg_1, 1), matrix(1, 1, NS_1))
                F0_1                = sort(apply(as.matrix(1-Sc0_1), 2, function(a) wtd.mean(a, weights = weight1)))
                marginal_fitted_1   = getquantile(depeval_1, F0_1, taus)
            }

            if (treatment && decomposition) {
                y_thred             = c(depeval, depeval_1)
                thredIndex          = c(1, floor((1:(nreg - 2) - 0.5)/(nreg - 2) * length(y_thred)), length(y_thred))
                y_thred_approx      = y_thred[thredIndex]
                y_te                = c(dep0, dep1)
                F_obs               = colMeans(apply(as.matrix(y_thred_approx), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))

                F_approx_fun0       = approxfun(depeval, F0, yleft = 0, yright = 1)
                F_approx_fun1       = approxfun(depeval_1, F0_1, yleft = 0, yright = 1)
                F_0hat              = length(dep0)/(length(dep0) + length(dep1))
                F_est               = F_approx_fun0(y_thred_approx) * F_0hat + F_approx_fun1(y_thred_approx) * (1 - F_0hat)
                F_ms                = F_obs - F_est
            }
        }

        if (method == "logit" | method == "probit" | method ==  "lpm") {
            
            dep2_0                      = depeval
            nreg                        = length(dep2_0)
            cdf_obs                     <- apply(as.matrix(dep2_0), 1, function(xxx) wtd.mean(dep0 <= xxx, weights = weight0))

            if (treatment) {
                dep2_1                  = depeval_1
                cdf_obs1                <- apply(as.matrix(dep2_1), 1, function(xxx) wtd.mean(dep1 <= xxx, weights = weight1))
            }

            if (method == "logit") {
                conditional_coef                = apply(as.matrix(dep2_0), 1, function(xxx) coef(glm((dep0 <= xxx) ~ reg0, weights = weight0, family = binomial(link = "logit"))))
                conditional_coef                = conditional_coef[, colSums(!is.nan(conditional_coef)) > 0]
                pred                            = plogis(cbind(1, reg0) %*% conditional_coef)
                conditional_fitted              = sort(apply(pred, 2, function(a) wtd.mean(a, weights = weight0)))
                marginal_fitted_cdf             <- conditional_fitted
   
                if (max(conditional_fitted) == 1) {
                    conditional_fitted_s        = conditional_fitted
                } else { conditional_fitted_s        = c(conditional_fitted, 1)}
                marginal_fitted                 = getquantile(dep2_0, conditional_fitted_s, taus)

                predCounter                     = plogis(cbind(1, counterfactual_var) %*% conditional_coef)
                conditional_fittedCounter       = sort(apply(predCounter, 2, function(a) wtd.mean(a, weights = weight1)))
                marginal_counterfactual_cdf     <- conditional_fittedCounter

                if (max(conditional_fittedCounter) == 1) {
                    predCounter_s               = conditional_fittedCounter
                } else { predCounter_s          = c(conditional_fittedCounter, 1)}
                marginal_counterfactual         = getquantile(dep2_0, predCounter_s, taus)

                if (treatment) {
                    conditional_coef_1          = apply(as.matrix(dep2_1), 1, function(xxx) coef(glm((dep1 <= xxx) ~ reg1, weights = weight1, family = binomial(link = "logit"))))
                    conditional_coef_1          = conditional_coef_1[, colSums(!is.nan(conditional_coef_1)) > 0]
                    pred_1                      = plogis(cbind(1, reg1) %*% conditional_coef_1)
                    conditional_fitted_1        = sort(apply(pred_1, 2, function(a) wtd.mean(a, weights = weight1)))
                    marginal_fitted_1_cdf       <- conditional_fitted_1

                    if (max(conditional_fitted_1) == 1) {
                        conditional_fitted_1_s  = conditional_fitted_1
                    } else {
                        conditional_fitted_1_s  = c(conditional_fitted_1, 1)
                    }
                    marginal_fitted_1           = getquantile(dep2_1, conditional_fitted_1_s, taus)
                }

                if (treatment && decomposition) {
                    y_thred                     = c(depeval, depeval_1)
                    thredIndex                  = c(1, floor((1:(nreg - 2) - 0.5)/(nreg - 2) * length(y_thred)), length(y_thred))
                    y_thred_approx              = y_thred[thredIndex]
                    y_te                        = c(dep0, dep1)
                    F_obs                       = colMeans(apply(as.matrix(y_thred_approx), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))

                    F_approx_fun0               = approxfun(dep2_0, conditional_fitted, yleft = 0, yright = 1)
                    F_approx_fun1               = approxfun(dep2_1, conditional_fitted_1, yleft = 0, yright = 1)
                    F_0hat                      = length(dep0)/(length(dep0) + length(dep1))
                    F_est                       = F_approx_fun0(y_thred_approx) * F_0hat + F_approx_fun1(y_thred_approx) * (1 - F_0hat)
                    F_ms                        = F_obs - F_est
                }
            }

            if (method == "probit") {
                conditional_coef                = apply(as.matrix(dep2_0), 1, function(xxx) coef(glm((dep0 <= xxx) ~ reg0, weights = weight0, family = binomial(link = "probit"))))
                conditional_coef                = conditional_coef[, colSums(!is.nan(conditional_coef)) > 0]
                pred                            = pnorm(cbind(1, reg0) %*% conditional_coef)
                conditional_fitted              = sort(apply(pred, 2, function(a) wtd.mean(a, weights = weight0)))
                marginal_fitted_cdf             <- conditional_fitted
           
                if (max(conditional_fitted) == 1) {
                    conditional_fitted_s        = conditional_fitted
                } else { conditional_fitted_s   = c(conditional_fitted, 1)}
                marginal_fitted                 = getquantile(dep2_0, conditional_fitted_s, taus)

                predCounter                     = pnorm(cbind(1, counterfactual_var) %*% conditional_coef)
                conditional_fittedCounter       = sort(apply(predCounter, 2, function(a) wtd.mean(a, weights = weight1)))
                marginal_counterfactual_cdf     <- conditional_fittedCounter

                if (max(conditional_fittedCounter) == 1) {
                    predCounter_s               = conditional_fittedCounter
                } else { predCounter_s          = c(conditional_fittedCounter, 1)}
                marginal_counterfactual         = getquantile(dep2_0, predCounter_s, taus)

                if (treatment) {
                    conditional_coef_1          = apply(as.matrix(dep2_1), 1, function(xxx) coef(glm((dep1 <= xxx) ~ reg1, weights = weight1, family = binomial(link = "probit"))))
                    conditional_coef_1          = conditional_coef_1[, colSums(!is.nan(conditional_coef_1)) > 0]
                    pred_1                      = pnorm(cbind(1, reg1) %*% conditional_coef_1)
                    conditional_fitted_1        = sort(apply(pred_1, 2, function(a) wtd.mean(a, weights = weight1)))
                    marginal_fitted_1_cdf       <- conditional_fitted_1

                    if (max(conditional_fitted_1) == 1) {
                        conditional_fitted_1_s  = conditional_fitted_1
                    } else { 
                        conditional_fitted_1_s  = c(conditional_fitted_1, 1)
                    }
                    marginal_fitted_1           = getquantile(dep2_1, conditional_fitted_1_s, taus)
                }

                if (treatment && decomposition) {
                    y_thred                     = c(depeval, depeval_1)
                    thredIndex                  = c(1, floor((1:(nreg - 2) - 0.5)/(nreg - 2) * length(y_thred)), length(y_thred))
                    y_thred_approx              = y_thred[thredIndex]
                    y_te                        = c(dep0, dep1)
                    F_obs                       = colMeans(apply(as.matrix(y_thred_approx), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))

                    F_approx_fun0               = approxfun(dep2_0, conditional_fitted, yleft = 0, yright = 1)
                    F_approx_fun1               = approxfun(dep2_1, conditional_fitted_1, yleft = 0, yright = 1)
                    F_0hat                      = length(dep0)/(length(dep0) + length(dep1))
                    F_est                       = F_approx_fun0(y_thred_approx) * F_0hat + F_approx_fun1(y_thred_approx) * (1 - F_0hat)
                    F_ms                        = F_obs - F_est
                }
            }

            if (method == "lpm") {
                conditional_coef                = apply(as.matrix(dep2_0), 1, function(xxx) coef(glm((dep0 <= xxx) ~ reg0, weights = weight0)))
                conditional_coef                = conditional_coef[, colSums(!is.nan(conditional_coef)) > 0]
                pred                            = cbind(1, reg0) %*% conditional_coef
                conditional_fitted              = sort(apply(pred, 2, function(a) wtd.mean(a, weights = weight0)))
                marginal_fitted_cdf             <- conditional_fitted

                if (max(conditional_fitted) == 1) {
                    conditional_fitted_s = conditional_fitted
                }else {conditional_fitted_s     = c(conditional_fitted, 1)}
                marginal_fitted                 = getquantile(dep2_0, conditional_fitted_s, taus)

                predCounter                     = cbind(1, counterfactual_var) %*% conditional_coef
                conditional_fittedCounter       = sort(apply(predCounter, 2, function(a) wtd.mean(a, weights = weight1)))
                marginal_counterfactual_cdf     <- conditional_fittedCounter

                if (max(conditional_fittedCounter) == 1) {
                    predCounter_s               = conditional_fittedCounter
                } else { predCounter_s          = c(conditional_fittedCounter, 1)}
                marginal_counterfactual         = getquantile(dep2_0, predCounter_s, taus)
                
                if (treatment) {
                    conditional_coef_1          = apply(as.matrix(dep2_1), 1, function(xxx) coef(glm((dep1 <= xxx) ~ reg1, weights = weight1)))
                    conditional_coef_1          = conditional_coef_1[, colSums(!is.nan(conditional_coef_1)) > 0]
                    pred_1                      = cbind(1, reg1) %*% conditional_coef_1
                    conditional_fitted_1        = sort(apply(pred_1, 2, function(a) wtd.mean(a, weights = weight1)))
                    marginal_fitted_1_cdf       <- conditional_fitted_1

                    if (max(conditional_fitted_1) == 1) { 
                        conditional_fitted_1_s  = conditional_fitted_1
                    } else { 
                        conditional_fitted_1_s  = c(conditional_fitted_1, 1) 
                    }
                    marginal_fitted_1           = getquantile(dep2_1, conditional_fitted_1_s, taus)
                }

                if (treatment && decomposition) {
                    y_thred                     = c(depeval, depeval_1)
                    thredIndex                  = c(1, floor((1:(nreg - 2) - 0.5)/(nreg - 2) * length(y_thred)), length(y_thred))
                    y_thred_approx              = y_thred[thredIndex]
                    y_te                        = c(dep0, dep1)
                    F_obs                       = colMeans(apply(as.matrix(y_thred_approx), 1, function(xxx) {ifelse(y_te <= xxx, 1, 0)}))

                    F_approx_fun0               = approxfun(dep2_0, conditional_fitted, yleft = 0, yright = 1)
                    F_approx_fun1               = approxfun(dep2_1, conditional_fitted_1, yleft = 0, yright = 1)
                    F_0hat                      = length(dep0)/(length(dep0) + length(dep1))
                    F_est                       = F_approx_fun0(y_thred_approx) * F_0hat + F_approx_fun1(y_thred_approx) * (1 - F_0hat)
                    F_ms                        = F_obs - F_est
                }
            }
        }
    }

    if (treatment && decomposition) {
        structral_effect            = marginal_fitted_1 - marginal_counterfactual
        composition_effect          = marginal_counterfactual - marginal_fitted
        total_effect                = marginal_fitted_1 - marginal_fitted
        quantile_effect             = structral_effect
    } else if (treatment && !decomposition) {
        structral_effect            = marginal_fitted_1 - marginal_counterfactual
        quantile_effect             = structral_effect
    } else {
        composition_effect          = marginal_counterfactual - marginal_fitted
        quantile_effect             = composition_effect
    }

    marginal_obs                    = wtd.quantile(dep0, weights = weight0 * length(weight0), probs = taus, na.rm = TRUE)
    names(marginal_fitted)          = NULL
    names(marginal_counterfactual)  = NULL
    names(marginal_obs)             = NULL
    
    if (treatment) { 
        marginal_obs_1              = wtd.quantile(dep1, weights = weight1 * length(weight1), probs = taus, na.rm = TRUE)
    }

    if (treatment && decomposition) {
        if(method=="logit" | method=="probit" | method=="lpm"){
            return(list(taus = taus, quantile_effect = quantile_effect, composition_effect = composition_effect, total_effect = total_effect, marginal_obs = marginal_obs, marginal_fitted = marginal_fitted, marginal_counterfactual = marginal_counterfactual, nreg = nreg, F_ms = F_ms, marginal_obs_1 = marginal_obs_1, marginal_fitted_1 = marginal_fitted_1, marginal_fitted_cdf= marginal_fitted_cdf, cdf_obs = cdf_obs, marginal_counterfactual_cdf = marginal_counterfactual_cdf, marginal_fitted_1_cdf = marginal_fitted_1_cdf, cdf_obs1=cdf_obs1))
        } else{
            return(list(taus = taus, quantile_effect = quantile_effect, composition_effect = composition_effect, total_effect = total_effect, marginal_obs = marginal_obs, marginal_fitted = marginal_fitted, marginal_counterfactual = marginal_counterfactual, nreg = nreg, F_ms = F_ms, marginal_obs_1 = marginal_obs_1, marginal_fitted_1 = marginal_fitted_1))
        }
    } else if (treatment && !decomposition) {
        if(method=="logit" | method=="probit" | method=="lpm"){
            return(list(taus = taus, quantile_effect = quantile_effect, marginal_obs = marginal_obs, marginal_fitted = marginal_fitted, marginal_counterfactual = marginal_counterfactual, nreg = nreg, F_ms = F_ms, marginal_obs_1 = marginal_obs_1, marginal_fitted_1 = marginal_fitted_1, marginal_fitted_cdf= marginal_fitted_cdf, cdf_obs = cdf_obs, marginal_counterfactual_cdf = marginal_counterfactual_cdf, marginal_fitted_1_cdf = marginal_fitted_1_cdf, cdf_obs1 = cdf_obs1))
        } else{
            return(list(taus = taus, quantile_effect = quantile_effect, marginal_obs = marginal_obs, marginal_fitted = marginal_fitted, marginal_counterfactual = marginal_counterfactual, nreg = nreg, F_ms = F_ms, marginal_obs_1 = marginal_obs_1, marginal_fitted_1 = marginal_fitted_1))
        }
    } else {
        if(method=="logit" | method=="probit" | method=="lpm"){
            return(list(taus = taus, quantile_effect = quantile_effect, marginal_obs = marginal_obs, marginal_fitted = marginal_fitted, marginal_counterfactual = marginal_counterfactual, nreg = nreg, marginal_fitted_cdf = marginal_fitted_cdf, cdf_obs = cdf_obs, marginal_counterfactual_cdf = marginal_counterfactual_cdf))        
        } else{
            return(list(taus = taus, quantile_effect = quantile_effect, marginal_obs = marginal_obs, marginal_fitted = marginal_fitted, marginal_counterfactual = marginal_counterfactual, nreg = nreg))
        }
    }
}

BootstrapProcedure <- function (transformation, reps, data0, depeval, data1, weightedboot, obs0, obs1, dep1, depeval_1, treatment, decomposition, taus, method, trimming, nreg, censoring, right, nsteps, firstc, secondc, seed, sepcore, ncore){

    subsamplek <- function(transformation, data0, data1, weightedboot, obs0, obs1, dep1, depeval_1, treatment, decomposition, taus, method, trimming, nreg, censoring, right, nsteps, firstc, secondc) {
        
        if (weightedboot) {

            #bootstrap_weight                    = data0$weight*rexp(obs0, rate = 1)
            bootstrap_weight                    = rexp(obs0, rate = 1)
            if(transformation){
                bootstrap_counterfactual_weight = bootstrap_weight
            }else{
                #bootstrap_counterfactual_weight = data1$weight*rexp(obs1, rate = 1)
                bootstrap_counterfactual_weight =  rexp(obs1, rate = 1)
            }
            sub_PE_res                          = QteDistEst(data0$dep, depeval, data0$reg, bootstrap_weight, data1$counterfactual, dep1, depeval_1, treatment, decomposition, bootstrap_counterfactual_weight, taus, method, trimming, nreg, data0$scale_variable, data1$counterfactual_scale_variable, censoring, right, nsteps, firstc, secondc)       
       
        } else{
            subsample_Index0                    = sample(obs0, obs0, replace = TRUE)

            subsampled_dep0                     = (data0$dep)[subsample_Index0]
            subsampled_reg0                     = as.matrix(data0$reg)[subsample_Index0, ]
            subsampled_weight0                  = (data0$weight)[subsample_Index0]
            subsampled_scale_variable           = as.matrix(data0$scale_variable)[subsample_Index0, ]

            if (transformation) {
                subsample_Index1                = subsample_Index0
            } else { 
                subsample_Index1                = sample(obs1, obs1, replace = TRUE)
            }
            subsampled_counterfactual_var       = as.matrix(data1$counterfactual)[subsample_Index1, ]
            subsampled_weight1                  = (data1$weight)[subsample_Index1]
            subsampled_counterfactual_scale_var = as.matrix(data1$scale_variable)[subsample_Index1, ]
            if (treatment) { 
                subsampled_dep1                 = dep1[subsample_Index1] 
            }
            sub_PE_res                          = QteDistEst(subsampled_dep0, depeval, subsampled_reg0, subsampled_weight0, subsampled_counterfactual_var, subsampled_dep1, depeval_1, treatment, decomposition, subsampled_weight1, taus, method, trimming, nreg, subsampled_scale_variable, subsampled_counterfactual_scale_var, censoring, right, nsteps, firstc, secondc)
        }

        if (treatment && decomposition) {
            if(method=="logit" | method=="probit" | method=="lpm"){
                qte_cov_boot = c(sub_PE_res$quantile_effect, sub_PE_res$marginal_obs, sub_PE_res$marginal_fitted, sub_PE_res$marginal_counterfactual, sub_PE_res$composition_effect, sub_PE_res$total_effect, sub_PE_res$marginal_obs_1, sub_PE_res$marginal_fitted_1, sub_PE_res$F_ms, sub_PE_res$marginal_fitted_cdf, sub_PE_res$cdf_obs,sub_PE_res$marginal_counterfactual_cdf, sub_PE_res$marginal_fitted_1_cdf, sub_PE_res$cdf_obs1)
            } else{
                qte_cov_boot = c(sub_PE_res$quantile_effect, sub_PE_res$marginal_obs, sub_PE_res$marginal_fitted, sub_PE_res$marginal_counterfactual, sub_PE_res$composition_effect, sub_PE_res$total_effect, sub_PE_res$marginal_obs_1, sub_PE_res$marginal_fitted_1, sub_PE_res$F_ms)
            }
        } else if (treatment && !decomposition) {
            if(method=="logit" | method=="probit" | method=="lpm"){
                qte_cov_boot = c(sub_PE_res$quantile_effect, sub_PE_res$marginal_obs, sub_PE_res$marginal_fitted, sub_PE_res$marginal_counterfactual, sub_PE_res$marginal_obs_1, sub_PE_res$marginal_fitted_1, sub_PE_res$F_ms, sub_PE_res$marginal_fitted_cdf, sub_PE_res$cdf_obs,sub_PE_res$marginal_counterfactual_cdf, sub_PE_res$marginal_fitted_1_cdf,sub_PE_res$cdf_obs1)
            } else {
                qte_cov_boot = c(sub_PE_res$quantile_effect, sub_PE_res$marginal_obs, sub_PE_res$marginal_fitted, sub_PE_res$marginal_counterfactual, sub_PE_res$marginal_obs_1, sub_PE_res$marginal_fitted_1, sub_PE_res$F_ms)
            }
        } else {
            if(method=="logit" | method=="probit" | method=="lpm"){
                qte_cov_boot = c(sub_PE_res$quantile_effect, sub_PE_res$marginal_obs, sub_PE_res$marginal_fitted, sub_PE_res$marginal_counterfactual, sub_PE_res$marginal_fitted_cdf, sub_PE_res$cdf_obs,sub_PE_res$marginal_counterfactual_cdf)
            } else{
                qte_cov_boot = c(sub_PE_res$quantile_effect, sub_PE_res$marginal_obs, sub_PE_res$marginal_fitted, sub_PE_res$marginal_counterfactual)
            }
        }
        return(qte_cov_boot=qte_cov_boot)
    }

    if (sepcore) {
        if (ncore >= 1) {
            cores <- detectCores()
            if (ncore >= cores) {
                stop("The number of cores specified is bigger than or equal to the computer has.")
            } else {
                if (ncore > 1) {cat(format(paste("cores used=", ncore, "\n"), width = 50))}
                cl            <- makeCluster(ncore)
                registerDoParallel(cl)
                qte_cov_booti <- foreach(i = 1:reps, .combine = cbind, .packages = c("quantreg", "Hmisc", "survival","doRNG"), .export = c("QteDistEst","getquantile"), .options.RNG = seed) %dorng% {
                    tryCatch({subsamplek(transformation, data0, data1, weightedboot, obs0, obs1, dep1, depeval_1, treatment, decomposition, taus, method, trimming, nreg, censoring, right, nsteps, firstc, secondc)}, error = function(e) NULL)
                }
                stopCluster(cl)
            }
        } else { 
            stop("The number of cores specified needs to be an integer.")
        }
    } else {
        qte_cov_booti <- foreach(i = 1:reps, .combine = cbind) %do% {
            tryCatch({subsamplek(transformation, data0, data1, weightedboot, obs0, obs1, dep1, depeval_1, treatment, decomposition, taus, method, trimming, nreg, censoring, right, nsteps, firstc, secondc)}, error = function(e) NULL)
        }
    }
    qte_cov_booti = t(qte_cov_booti)
    return(qte_cov_booti = qte_cov_booti)
}

InferenceTestingEval <- function (method, qte_cov_booti, reps, taus, robust, qte_cov_obs, last, first, alpha, cons_test, treatment, decomposition, depeval, depeval_1) {
    
    nqs     = length(taus)
    reps    = nrow(qte_cov_booti)
    sel0    = ((sum(taus < first)) + 1):(sum(taus <= last))
    neval   <- length(depeval)
    neval1  <- length(depeval_1)
    ntotal  <- length(qte_cov_obs)

    if (is.element(0.5, taus)) {
        median_sel                      = which(taus == 0.5)
    } else {
        print("Q_0.5 is not evaluated. This is the median of the specified rather than the Q_0.5")
        median_sel                      = ceiling(nqs/2)
    }
    sel_median                          = c((sum(taus < first) + 1):(median_sel - 1), (median_sel + 1):sum(taus <= last))
    selall                              = 1:length(sel_median)

    quantile_effect_bootstrap           = qte_cov_booti[, 1:nqs]
    marginal_obs_bootstrap              = qte_cov_booti[, (nqs + 1):(2 * nqs)]
    marginal_fitted_bootstrap           = qte_cov_booti[, (2 * nqs + 1):(3 * nqs)]
    marginal_counterfactual_bootstrap   = qte_cov_booti[, (3 * nqs + 1):(4 * nqs)]

    quantile_effect                     = qte_cov_obs[1:nqs]
    marginal_obs                        = qte_cov_obs[(nqs + 1):(2 * nqs)]
    marginal_fitted                     = qte_cov_obs[(2 * nqs + 1):(3 * nqs)]
    marginal_counterfactual             = qte_cov_obs[(3 * nqs + 1):(4 * nqs)]

    if (treatment && decomposition) {
        composition_effect_bootstrap    = qte_cov_booti[, (4 * nqs + 1):(5 * nqs)]
        total_effect_bootstrap          = qte_cov_booti[, (5 * nqs + 1):(6 * nqs)]
        marginal_obs_1_bootstrap        = qte_cov_booti[, (6 * nqs + 1):(7 * nqs)]
        marginal_fitted_1_bootstrap     = qte_cov_booti[, (7 * nqs + 1):(8 * nqs)]

        composition_effect              = qte_cov_obs[(4 * nqs + 1):(5 * nqs)]
        total_effect                    = qte_cov_obs[(5 * nqs + 1):(6 * nqs)]
        marginal_obs_1                  = qte_cov_obs[(6 * nqs + 1):(7 * nqs)]
        marginal_fitted_1               = qte_cov_obs[(7 * nqs + 1):(8 * nqs)]

        if(method=="logit" | method=="probit" | method=="lpm"){
            cdf_fitted_bootstrap        <- qte_cov_booti[, (ntotal-3*neval-2*neval1+1):(ntotal-2*neval-2*neval1)]
            cdf_obs_bootstrap           <- qte_cov_booti[, (ntotal-2*neval-2*neval1+1):(ntotal-neval-2*neval1)]
            cdf_counter_bootstrap       <- qte_cov_booti[, (ntotal-neval-2*neval1+1):(ntotal-2*neval1)]
            cdf_fitted1_bootstrap       <- qte_cov_booti[, (ntotal-2*neval1+1):(ntotal-neval1)]
            cdf_obs1_bootstrap          <- qte_cov_booti[, (ntotal-neval1+1):ntotal]

            cdf_fitted                  <- qte_cov_obs[(ntotal-3*neval-2*neval1+1):(ntotal-2*neval-2*neval1)]
            cdf_obs                     <- qte_cov_obs[(ntotal-2*neval-2*neval1+1):(ntotal-neval-2*neval1)]
            cdf_counter                 <- qte_cov_obs[(ntotal-neval-2*neval1+1):(ntotal-2*neval1)]
            cdf_fitted1                 <- qte_cov_obs[(ntotal-2*neval1+1):(ntotal-neval1)]
            cdf_obs1                    <- qte_cov_obs[(ntotal-neval1+1):ntotal]

            F_ms_bootstrap              = qte_cov_booti[, (8 * nqs+1):(ntotal-3*neval-2*neval1)]
            F_ms_obs                    = qte_cov_obs[(8 * nqs+1):(ntotal-3*neval-2*neval1)]
            sel_ms                      = 1:length(F_ms_obs)   
        } else{
            F_ms_bootstrap              = qte_cov_booti[,-(1:8 * nqs)]
            F_ms_obs                    = qte_cov_obs[-(1:8 * nqs)]
            sel_ms                      = 1:length(F_ms_obs)
        }

    } else if (treatment && !decomposition) {
        marginal_obs_1_bootstrap        = qte_cov_booti[, (4 * nqs + 1):(5 * nqs)]
        marginal_fitted_1_bootstrap     = qte_cov_booti[, (5 * nqs + 1):(6 * nqs)]

        marginal_obs_1                  = qte_cov_obs[(4 * nqs + 1):(5 * nqs)]
        marginal_fitted_1               = qte_cov_obs[(5 * nqs + 1):(6 * nqs)]

        if(method=="logit" | method=="probit" | method=="lpm"){
            cdf_fitted_bootstrap        <- qte_cov_booti[, (ntotal-3*neval-2*neval1+1):(ntotal-2*neval-2*neval1)]
            cdf_obs_bootstrap           <- qte_cov_booti[, (ntotal-2*neval-2*neval1+1):(ntotal-neval-2*neval1)]
            cdf_counter_bootstrap       <- qte_cov_booti[, (ntotal-neval-2*neval1+1):(ntotal-2*neval1)]
            cdf_fitted1_bootstrap       <- qte_cov_booti[, (ntotal-2*neval1+1):(ntotal-neval1)]
            cdf_obs1_bootstrap          <- qte_cov_booti[, (ntotal-neval1+1):ntotal]

            cdf_fitted                  <- qte_cov_obs[(ntotal-3*neval-2*neval1+1):(ntotal-2*neval-2*neval1)]
            cdf_obs                     <- qte_cov_obs[(ntotal-2*neval-2*neval1+1):(ntotal-neval-2*neval1)]
            cdf_counter                 <- qte_cov_obs[(ntotal-neval-2*neval1+1):(ntotal-2*neval1)]
            cdf_fitted1                 <- qte_cov_obs[(ntotal-2*neval1+1):(ntotal-neval1)]
            cdf_obs1                    <- qte_cov_obs[(ntotal-neval1+1):ntotal]

            F_ms_bootstrap              =  qte_cov_booti[,(6 * nqs+1):(ntotal-3*neval-2*neval1)]
            F_ms_obs                    =  qte_cov_obs[(6 * nqs+1):(ntotal-3*neval-2*neval1)]
            sel_ms                      =  1:length(F_ms_obs)            
        } else{
            F_ms_bootstrap              = qte_cov_booti[, -(1:6 * nqs)]
            F_ms_obs                    = qte_cov_obs[-(1:6 * nqs)]
            sel_ms                      = 1:length(F_ms_obs)
        }

    } else if(!treatment && !decomposition) {
        marginal_obs_1_bootstrap        = qte_cov_booti[, (nqs + 1):(2 * nqs)]
        marginal_fitted_1_bootstrap     = qte_cov_booti[, (2 * nqs + 1):(3 * nqs)]

        marginal_obs_1                  = qte_cov_obs[(nqs + 1):(2 * nqs)]
        marginal_fitted_1               = qte_cov_obs[(2 * nqs + 1):(3 * nqs)]

        if(method=="logit" | method=="probit" | method=="lpm"){
            cdf_fitted_bootstrap        <- qte_cov_booti[, (4 * nqs + 1):(4* nqs+neval)]
            cdf_obs_bootstrap           <- qte_cov_booti[, (4 * nqs + neval + 1):(4 * nqs + 2*neval)]
            cdf_counter_bootstrap       <- qte_cov_booti[, (4 * nqs + 2*neval + 1):(4 * nqs + 3*neval)]

            cdf_fitted                  <- qte_cov_obs[(4 * nqs + 1):(4* nqs+neval)]
            cdf_obs                     <- qte_cov_obs[(4 * nqs + neval + 1):(4 * nqs + 2*neval)]
            cdf_counter                 <- qte_cov_obs[(4 * nqs + 2*neval + 1):(4 * nqs + 3*neval)]
        }
    }

    ## variance for quantile_effect
    pe_boot                             = VarianceEval(quantile_effect_bootstrap, quantile_effect, reps, nqs, alpha, robust, sel0)
    pe_lb_point                         = pe_boot$qte_cov - pe_boot$se * qnorm(1 - alpha/2)
    pe_ub_point                         = pe_boot$qte_cov + pe_boot$se * qnorm(1 - alpha/2)

    if (!treatment && !decomposition) {        
        
        resCE                           = cbind(pe_boot$qte_cov, pe_boot$se, pe_lb_point, pe_ub_point, pe_boot$lb, pe_boot$ub)
        colnames(resCE)                 = c("composition.effect", "se.ce", "lb.ce.point", "ub.ce.point", "lb.ce", "ub.ce")

    } else if (treatment && !decomposition) {        
        
        resSE                           = cbind(pe_boot$qte_cov, pe_boot$se, pe_lb_point, pe_ub_point, pe_boot$lb, pe_boot$ub)
        colnames(resSE)                 = c("structure.effect", "se.se", "lb.se.point", "ub.se.point", "lb.se", "ub.se")

    } else if (treatment && decomposition) {

        ## variance for structure effect
        resSE                           = cbind(pe_boot$qte_cov, pe_boot$se, pe_lb_point, pe_ub_point, pe_boot$lb, pe_boot$ub)
        colnames(resSE)                 = c("structure.effect", "se.se", "lb.se.point", "ub.se.point", "lb.se", "ub.se")
        
        ## variance for composition effect
        ce_boot                         = VarianceEval(composition_effect_bootstrap, composition_effect, reps, nqs, alpha, robust, sel0)
        ce_lb_point                     = ce_boot$qte_cov - ce_boot$se * qnorm(1 - alpha/2)
        ce_ub_point                     = ce_boot$qte_cov + ce_boot$se * qnorm(1 - alpha/2)
        resCE                           = cbind(ce_boot$qte_cov, ce_boot$se, ce_lb_point, ce_ub_point, ce_boot$lb, ce_boot$ub)
        colnames(resCE)                 = c("composition.effect", "se.ce", "lb.ce.point", "ub.ce.point", "lb.ce", "ub.ce")

        ## variance for total effect
        te_boot                         = VarianceEval(total_effect_bootstrap, total_effect, reps, nqs, alpha, robust, sel0)
        te_lb_point                     = te_boot$qte_cov - te_boot$se * qnorm(1 - alpha/2)
        te_ub_point                     = te_boot$qte_cov + te_boot$se * qnorm(1 - alpha/2)
        resTE                           = cbind(te_boot$qte_cov, te_boot$se, te_lb_point, te_ub_point, te_boot$lb, te_boot$ub)
        colnames(resTE)                 = c("total.effect", "se.te", "lb.te.point", "ub.te.point", "lb.te", "ub.te")
    }

    ## Distribution 
    sample_quantile_ref0                = VarianceEval(marginal_obs_bootstrap, marginal_obs, reps, nqs, alpha, robust, sel0)
    colnames(sample_quantile_ref0)      = c("ME.ref0", "se.ref0", "lb.ref0", "ub.ref0")

    model_quantile_ref0                 = VarianceEval(marginal_fitted_bootstrap, marginal_fitted, reps, nqs, alpha, robust, sel0)
    colnames(model_quantile_ref0)       = c("ME.fitted0", "se.fitted0", "lb.fitted0", "ub.fitted0")

    model_quantile_counter              = VarianceEval(marginal_counterfactual_bootstrap, marginal_counterfactual, reps, nqs, alpha, robust, sel0)
    colnames(model_quantile_counter)    = c("ME.counter", "se.counter", "lb.counter", "ub.counter")

    if (treatment) {
        sample_quantile_ref1            = VarianceEval(marginal_obs_1_bootstrap, marginal_obs_1, reps, nqs, alpha, robust, sel0)
        colnames(sample_quantile_ref1)  = c("ME.ref1", "se.ref1", "lb.ref1", "ub.ref1")

        model_quantile_ref1             = VarianceEval(marginal_fitted_1_bootstrap, marginal_fitted_1, reps, nqs, alpha, robust, sel0)
        colnames(model_quantile_ref1)   = c("ME.fitted1", "se.fitted1", "lb.fitted1", "ub.fitted1")
    }

    ##********************************************************* Testing *************************************
    ## test of no misspecification, obs-fitted
    if (method == "logit" | method == "lpm"){ 
        test_MS                         = c(NA, NA)
    }else{
        qte_cov_boot_ms                 = marginal_obs_bootstrap - marginal_fitted_bootstrap
        qte_cov_ms                      = marginal_obs - marginal_fitted
        test_MS                         = TestingEval((qte_cov_boot_ms - kronecker(matrix(qte_cov_ms, 1, nqs), matrix(1, reps, 1)))[, sel0], qte_cov_ms[sel0], qte_cov_boot_ms, sel0, reps, robust)
    }
    ## test of no effect, test_const == 0
    test_0                              = TestingEval((quantile_effect_bootstrap - kronecker(matrix(quantile_effect, 1, nqs), matrix(1, reps, 1)))[, sel0], quantile_effect[sel0], quantile_effect_bootstrap, sel0, reps, robust)
    ## test of const effect
    constest1                           = sort(unique(cons_test))
    constestNo0                         = setdiff(constest1, 0)
    nc                                  = length(constestNo0)
    if (nc > 0) {
        test_const = matrix(0, nc, 2)
        for (i in 1:nc) {
            test_const[i, ]             = TestingEval((quantile_effect_bootstrap - kronecker(matrix(quantile_effect, 1, nqs), matrix(1, reps, 1)))[, sel0], quantile_effect[sel0] - constestNo0[i], quantile_effect_bootstrap, sel0, reps, robust)
        }
    } else {
        test_const                      = NULL
    }
    ## test of median
    qte_cov_boot_median                 = quantile_effect_bootstrap[, sel_median] - quantile_effect_bootstrap[, median_sel]
    qte_cov_def_median                  = quantile_effect[sel_median] - quantile_effect[median_sel]
    test_median                         = TestingEval(qte_cov_boot_median - kronecker(matrix(qte_cov_def_median, 1, ncol(qte_cov_boot_median)), matrix(1, reps, 1)), qte_cov_def_median, qte_cov_boot_median, selall, reps, robust) 
    ## test of stochastic dominance
    qte_cov_boot_SD_ex                  = (quantile_effect_bootstrap - kronecker(matrix(quantile_effect, 1, nqs), matrix(1, reps, 1)))[, sel0]
    test_SD                             = TestingEval((qte_cov_boot_SD_ex * (qte_cov_boot_SD_ex <= 0)), (quantile_effect * (quantile_effect <= 0))[sel0], quantile_effect_bootstrap, sel0, reps, robust)
    ## test of being stochastically dominated 
    test_SDD                            = TestingEval((qte_cov_boot_SD_ex * (qte_cov_boot_SD_ex >= 0)), (quantile_effect * (quantile_effect >= 0))[sel0], quantile_effect_bootstrap, sel0, reps, robust)
    
    ## Tests
    if (!treatment && !decomposition) {
        testCE = rbind(test_MS, test_0, test_const, test_median, test_SD, test_SDD)
    } else if (treatment && !decomposition) {
        testSE = rbind(test_MS, test_0, test_const, test_median, test_SD, test_SDD)
    } else if (treatment && decomposition) {

        ## *** SE ***
        testSE = rbind(test_MS, test_0, test_const, test_median, test_SD, test_SDD)

        ## *** CE ***
        ## test of no misspecification, obs-fitted
        if (method == "logit" | method == "lpm") {
            test_MS_CE                  = c(NA, NA)
        } else {
            qte_cov_boot_ms_CE          = marginal_obs_1_bootstrap - marginal_fitted_1_bootstrap
            qte_cov_ms_CE               = marginal_obs_1 - marginal_fitted_1
            test_MS_CE                  = TestingEval((qte_cov_boot_ms_CE - kronecker(matrix(qte_cov_ms_CE, 1, length(qte_cov_ms_CE)), matrix(1, reps, 1)))[, sel0], qte_cov_ms_CE[sel0], qte_cov_boot_ms_CE, sel0, reps, robust)
        }
        ## test of no effect, test_const = 0
        test_0_CE                       = TestingEval((composition_effect_bootstrap - kronecker(matrix(composition_effect, 1, nqs), matrix(1, reps, 1)))[, sel0], composition_effect[sel0], composition_effect_bootstrap, sel0, reps, robust)
        ## test of const effect
        if (nc > 0) {
            test_const_CE = matrix(0, nc, 2)
            for (i in 1:nc) {
                test_const_CE[i, ]      = TestingEval((composition_effect_bootstrap - kronecker(matrix(composition_effect, 1, nqs), matrix(1, reps, 1)))[, sel0], composition_effect[sel0] - constestNo0[i], composition_effect_bootstrap, sel0, reps, robust)
            }
        } else {
            test_const_CE               = NULL
        }
        ## test of median
        ce_cov_boot_median              = composition_effect_bootstrap[, sel_median] - composition_effect_bootstrap[, median_sel]
        ce_cov_def_median               = composition_effect[sel_median] - composition_effect[median_sel]
        test_median_CE                  = TestingEval(ce_cov_boot_median - kronecker(matrix(ce_cov_def_median, 1, ncol(ce_cov_boot_median)), matrix(1, reps, 1)), ce_cov_def_median, ce_cov_boot_median, selall, reps, robust)
        ## test of stochastic dominance
        ce_cov_boot_SD_ex               = (composition_effect_bootstrap - kronecker(matrix(composition_effect, 1, nqs), matrix(1, reps, 1)))[, sel0]
        test_SD_CE                      = TestingEval((ce_cov_boot_SD_ex * (ce_cov_boot_SD_ex <= 0)), (composition_effect * (composition_effect <= 0))[sel0], composition_effect_bootstrap, sel0, reps, robust)
        ## test of being stochastically dominated
        test_SDD_CE                     = TestingEval((ce_cov_boot_SD_ex * (ce_cov_boot_SD_ex >= 0)), (composition_effect * (composition_effect >= 0))[sel0], composition_effect_bootstrap, sel0, reps, robust)
        testCE                          = rbind(test_MS_CE, test_0_CE, test_const_CE, test_median_CE, test_SD_CE, test_SDD_CE)
        
        ### ***** TE ****
        ## test of no misspecification, obs-fitted
        if (method == "logit" | method == "lpm") {
            test_MS_TE                  = c(NA, NA)
        } else {
            qte_cov_boot_ms_TE          = F_ms_bootstrap
            qte_cov_ms_TE               = F_ms_obs
            test_MS_TE                  = TestingEval((qte_cov_boot_ms_TE - kronecker(matrix(qte_cov_ms_TE, 1, length(qte_cov_ms_TE)), matrix(1, reps, 1))), qte_cov_ms_TE, qte_cov_boot_ms_TE, sel_ms, reps, robust)
        }
        ## test of no effect, test_const = 0 
        test_0_TE                       = TestingEval((total_effect_bootstrap - kronecker(matrix(total_effect, 1, nqs), matrix(1, reps, 1)))[, sel0], total_effect[sel0], total_effect_bootstrap, sel0, reps, robust)
        ## test of constant effect
        if (nc > 0) {
            test_const_TE = matrix(0, nc, 2)
            for (i in 1:nc) {
                test_const_TE[i, ]      = TestingEval((total_effect_bootstrap - kronecker(matrix(total_effect, 1, nqs), matrix(1, reps, 1)))[, sel0], total_effect[sel0] - constestNo0[i], total_effect_bootstrap, sel0, reps, robust)
            }
        } else {
            test_const_TE               = NULL
        }
        ## test of median
        te_cov_boot_median              = total_effect_bootstrap[, sel_median] - total_effect_bootstrap[, median_sel]
        te_cov_def_median               = total_effect[sel_median] - total_effect[median_sel]
        test_median_TE                  = TestingEval(te_cov_boot_median - kronecker(matrix(te_cov_def_median, 1, ncol(te_cov_boot_median)), matrix(1, reps, 1)), te_cov_def_median, te_cov_boot_median, selall, reps, robust)
        ## test of stochastic dominance
        te_cov_boot_SD_ex               = (total_effect_bootstrap - kronecker(matrix(total_effect, 1, nqs), matrix(1, reps, 1)))[, sel0]
        test_SD_TE                      = TestingEval((te_cov_boot_SD_ex * (te_cov_boot_SD_ex <= 0)), (total_effect * (total_effect <= 0))[sel0], total_effect_bootstrap, sel0, reps, robust)
        ## test of being stochastically dominated
        test_SDD_TE                     = TestingEval((te_cov_boot_SD_ex * (te_cov_boot_SD_ex >= 0)), (total_effect * (total_effect >= 0))[sel0], total_effect_bootstrap, sel0, reps, robust)
        testTE                          = rbind(test_MS_TE, test_0_TE, test_const_TE, test_median_TE, test_SD_TE, test_SDD_TE)
    }

    res_boot                            = NULL
    res_boot$sample_quantile_ref0       = sample_quantile_ref0
    res_boot$model_quantile_ref0        = model_quantile_ref0
    res_boot$model_quantile_counter     = model_quantile_counter
    if (!treatment && !decomposition) {
        res_boot$resCE                  = resCE
        res_boot$testCE                 = testCE

    } else if (treatment && !decomposition) {
        res_boot$sample_quantile_ref1   = sample_quantile_ref1
        res_boot$model_quantile_ref1    = model_quantile_ref1

        res_boot$resSE                  = resSE
        res_boot$testSE                 = testSE

    } else if (treatment && decomposition) {
        res_boot$sample_quantile_ref1   = sample_quantile_ref1
        res_boot$model_quantile_ref1    = model_quantile_ref1

        res_boot$resSE                  = resSE
        res_boot$testSE                 = testSE

        res_boot$resCE                  = resCE
        res_boot$testCE                 = testCE

        res_boot$resTE                  = resTE
        res_boot$testTE                 = testTE
    }

    if(method=="logit" || method=="probit" || method=="lpm"){

        model_cdf_ref0                  <- VarianceEval(cdf_fitted_bootstrap, cdf_fitted, reps, neval, alpha, robust, cdf_fitted>=first & cdf_fitted<=last)
        colnames(model_cdf_ref0)        <- c("ME.cdf0", "se.cdf0", "lb.cdf0", "ub.cdf0")
        res_boot$model_cdf_ref0         <- model_cdf_ref0

        model_cdf_counter               <- VarianceEval(cdf_counter_bootstrap, cdf_counter, reps, neval, alpha, robust, cdf_counter>=first & cdf_counter<=last)
        colnames(model_cdf_counter)     <- c("ME.cdfcounter", "se.cdfcounter", "lb.cdfcounter", "ub.cdfcounter")
        res_boot$model_cdf_counter      <- model_cdf_counter

        cdf_obs0                        <- VarianceEval(cdf_obs_bootstrap, cdf_obs, reps, neval, alpha, robust, cdf_obs>=first & cdf_obs<=last)
        colnames(cdf_obs0)              <- c("ME.cdf0", "se.cdf0", "lb.cdf0", "ub.cdf0")
        res_boot$cdf_obs0               <- cdf_obs0

        if(treatment){
            model_cdf_ref1              <- VarianceEval(cdf_fitted1_bootstrap, cdf_fitted1, reps, neval1, alpha, robust, cdf_fitted1>=first & cdf_fitted1<=last)
            colnames(model_cdf_ref1)    <- c("ME.cdf1", "se.cdf1", "lb.cdf1", "ub.cdf1")
            res_boot$model_cdf_ref1     <- model_cdf_ref1

            cdf_obs1                    <- VarianceEval(cdf_obs1_bootstrap, cdf_obs1, reps, neval, alpha, robust, cdf_obs1>=first & cdf_obs1<=last)
            colnames(cdf_obs1)          <- c("ME.cdf1", "se.cdf1", "lb.cdf1", "ub.cdf1")
            res_boot$cdf_obs1           <- cdf_obs1
        }
    }
    return(res_boot)
}

VarianceEval <- function (qte_cov_boot, qte_cov, reps, nqs, alpha, robust, sel){
    if (robust) {
        seuqf   = apply(qte_cov_boot, 2, function(x) {(quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))/1.34 })
    } else {
        seuqf   = sqrt(diag(var(qte_cov_boot)))
    }
    Vuqf        = seuqf^2

    Kuqf        = sqrt((qte_cov_boot - kronecker(matrix(qte_cov, 1, nqs), matrix(1, reps, 1)))^2/kronecker(matrix(Vuqf, 1, nqs), matrix(1, reps, 1)))[, sel]
    Kuqfsel     = Kuqf[, apply(Kuqf, 2, function(x) all(is.finite(x)))]
    Kmaxuqf     = apply(Kuqfsel, 1, max)

    Kalpha      = quantile(Kmaxuqf, 1 - alpha, na.rm = TRUE, names = FALSE)
    lb          = qte_cov - seuqf * Kalpha
    ub          = qte_cov + seuqf * Kalpha
    simbootres  = data.frame(qte_cov = qte_cov, se = seuqf, lb = lb, ub = ub)
    return(simbootres)
}

TestingEval <- function (boot_test_numerator, obs_test_numerator, variable_for_variance, sel, reps, robust) {
    if (robust) {
        Vuqf    = (apply(variable_for_variance, 2, function(x) { (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))/1.34}))^2
    } else {
        Vuqf    = diag(var(variable_for_variance))
    }
    Vuqf        = Vuqf[sel] + 1e-09
    nqs         = ncol(boot_test_numerator)
    Kuqf        = sqrt(boot_test_numerator^2/kronecker(matrix(Vuqf, 1, nqs), matrix(1, reps, 1)))

    Kmaxuqf     = apply(Kuqf, 1, max)
    KSstat      = max(sqrt(obs_test_numerator^2/Vuqf))

    Kuqf2       = Kuqf^2
    Kmeanuqf    = rowMeans(Kuqf2)
    CMSstat     = mean(obs_test_numerator^2/Vuqf)

    testboot    = c(mean(Kmaxuqf > KSstat), mean(Kmeanuqf > CMSstat))
    return(testboot)
}

getquantile <- function(depevalg, pred, tausg) {
    Q = NULL
    for (i in 1:length(tausg)) {
        Q = rbind(Q, depevalg[min(sum(pred <= tausg[i]) + 1, length(depevalg))])
    }
    return(Q)
}

