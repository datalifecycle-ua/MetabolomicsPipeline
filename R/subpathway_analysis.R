#' Subpathway Analysis
#'
#' Subpathway analysis for metabolite data.
#'
#' @details For each metabolite, we test three models using using ANOVA.
#' \enumerate{
#'
#'  \item Interaction:  \eqn{log Peak = Treatment + block + Treatment*block}
#'
#'  \item Parallel: \eqn{log Peak = Treatment + block}
#'
#'  \item Single: \eqn{log Peak =  Treatment}
#'
#' }
#'
#' For the interaction model, we are focusing only on the interaction term
#'  "Treatment*block" to test if there is a significant interaction between our
#'  treatment and the block variable. The parallel model tests if the block
#' variable explains a significant amount of the metabolite variance, and the
#' treatment model tests if the treatment explains a significant proportion of
#'  the variance for each metabolite. Then, we use the Combined Fisher
#'   probability to test each model at the subpathway level.
#'
#' \deqn{\tilde{X} = -2\sum_{i=1}^k ln(p_i)}
#'
#' where \eqn{k} is the number of metabolites in the subpathway. We can
#' get a p-value from \eqn{P(X \geq\tilde{X})}, knowing that
#' \eqn{\tilde{X}\sim \chi^2_{2k}}. You will notice that smaller p-values will
#' lead to a larger \eqn{\tilde{X}}.
#'
#' @seealso [Loughin, Thomas M. "A systematic comparison of methods for combining p-values from independent tests." Computational statistics & data analysis 47.3 (2004): 467-485.](https://www.sciencedirect.com/science/article/pii/S0167947303002950)

#'
#' @param data SummarizedExperiment with metabolomics experiment data.
#'
#'
#' @param treat_var This is the name of the variable in the analysis data that
#'  is the main variable of interest.
#'
#' @param block_var This is the name of the blocking variable in the dataset.
#' If the the experimental design does not include a blocking variable, then the
#' value of block_var=NULL.
#'
#' @param strat_var Variable to stratify the subpathway analysis by. This is set
#' to NULL by default and will not stratify the analysis unless specified.
#'
#' @param Assay Name of the assay to be used for the pairwise analysis
#'  (default='normalized')
#'
#' @param subPathwayName Column name for subpathway variable as defined in the
#' chemical annotation worksheet.
#'
#' @param chemName Column name for chemical name variable as defined in the
#' chemical annotation worksheet.
#'
#' @param superPathwayName Column name for superpathway variable as defined in
#' the chemical annotation worksheet.
#'
#'
#'
#'
#' @returns A data frame with "CHEM_ID","sub_pathway","chem_name",
#' "interaction_pval","interaction_fisher","parallel_pval","parallel_fisher",
#' "single_pval","single_fisher",and "model" for each metabolite.
#'
#'
#' @examples
#' # Load data
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#' dat <- demoDataSmall
#'
#' # Runsubpathay analysis
#' sub_analysis <- subpathway_analysis(dat,
#'     treat_var = "GROUP_NAME",
#'     block_var = "TIME1",
#'     strat_var = NULL,
#'     Assay = "normalized"
#' )
#'
#' #############################################################################
#' ### Results Plots ###########################################################
#' #############################################################################
#'
#' # significant subpathways by model type
#' subpath_by_model(sub_analysis)
#'
#' # Percentage of signficant subpathways within superpathways
#' subpath_within_superpath(sub_analysis)
#'
#' met_within_sub(sub_analysis, subpathway = "Aminosugar Metabolism")
#'
#' # All signifiicant subpathways
#' all_sig_subpath(sub_analysis)
#'
#' @importFrom  dplyr rename
#' @importFrom stats lm
#' @importFrom stats pchisq
#'
#'
#' @export
#'
#'




subpathway_analysis <- function(data, treat_var, block_var = NULL,
                                strat_var = NULL, Assay = "normalized",
                                subPathwayName = "SUB_PATHWAY",
                                chemName = "CHEMICAL_NAME",
                                superPathwayName = "SUPER_PATHWAY") {
    # Create analysis data
    analysis <- SummarizedExperiment::colData(data) %>%
        merge(t(SummarizedExperiment::assay(data, Assay)), by = "row.names") 
    
    #%>%
     #   dplyr::rename(PARENT_SAMPLE_NAME = Row.names)


    # chemdata
    chem <- SummarizedExperiment::rowData(data)

    if (is.null(strat_var)) {
        if (!is.null(block_var)) {
            # Create results dataframe
            path_data <- data.frame(
                CHEM_ID = rownames(chem),
                sub_pathway = chem[, subPathwayName],
                chem_name = chem[, chemName],
                super_pathway = chem[, superPathwayName],
                interaction_pval = NA,
                interaction_fisher = NA,
                parallel_pval = NA,
                parallel_fisher = NA,
                single_pval = NA,
                single_fisher = NA,
                model = NA
            ) %>%
                dplyr::mutate(sub_pathway = ifelse(is.na(sub_pathway),
                    "Unknown Metabolite", sub_pathway
                ))


            # For each pathways run the models for each metabolite within the
            #  pathway.
            for (pathway in unique(path_data$sub_pathway)) {
                var_names <- path_data$CHEM_ID[which(
                    path_data$sub_pathway == pathway
                )]


                # Calculate p values for each linear model
                for (var in var_names) {
                    outcome <- analysis[, as.character(var)]

                    if (is.null(outcome)) next

                    # This Generalizes Model formulation
                    inter_formula <- as.formula(
                        paste("outcome",
                            paste(
                                c(
                                    treat_var, block_var,
                                    paste0(treat_var, "*", block_var)
                                ),
                                collapse = " + "
                            ),
                            sep = "~"
                        )
                    )

                    int_mod <- lm(inter_formula, data = analysis)

                    # TryCatch for Anova Modeling
                    tryCatch(
                        {
                            path_data$interaction_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(int_mod)[
                                    paste0(treat_var, ":", block_var), 5
                                ]
                        },
                        warning = function(w) {
                            warning(w, " for CHEM ID ", var)
                        },
                        finally = {
                            path_data$interaction_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(int_mod)[
                                    paste0(treat_var, ":", block_var), 5
                                ]
                        }
                    )

                    # Parallel model, this generalizes the model formulation
                    par_formula <- as.formula(
                        paste("outcome",
                            paste(c(treat_var, block_var), collapse = " + "),
                            sep = "~"
                        )
                    )
                    par_mod <- lm(par_formula, data = analysis)

                    # Anova modling for the parallel model.
                    tryCatch(
                        {
                            path_data$parallel_pval[which(
                                path_data$CHEM_ID == var
                            )] <-
                                anova(par_mod)[block_var, 5]
                        },
                        warning = function(w) {
                            warning(w, " for CHEM ID ", var)
                        },
                        finally = {
                            path_data$parallel_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(par_mod)[block_var, 5]
                        }
                    )

                    # Treat model
                    treat_formula <- as.formula(
                        paste("outcome",
                            paste(c(treat_var), collapse = " + "),
                            sep = "~"
                        )
                    )
                    treat_mod <- lm(treat_formula, data = analysis)

                    # Anova model for treatment.
                    tryCatch(
                        {
                            path_data$single_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(treat_mod)[treat_var, 5]
                        },
                        warning = function(w) {
                            warning(w, " for CHEM ID ", var)
                        },
                        finally = {
                            path_data$single_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(treat_mod)[treat_var, 5]
                        }
                    )
                }

                # Compute combined fisher probabilities for each model.
                if (length(var_names) == 1) {
                    path_data$interaction_fisher[
                        which(path_data$CHEM_ID == var)
                    ] <-
                        anova(int_mod)[paste0(treat_var, ":", block_var), 5]

                    path_data$parallel_fisher[
                        which(path_data$CHEM_ID == var)
                    ] <-
                        anova(par_mod)[block_var, 5]

                    path_data$single_fisher[which(path_data$CHEM_ID == var)] <-
                        anova(treat_mod)[treat_var, 5]

                    path_data$model[which(path_data$sub_pathway == pathway)] <-
                        case_when(
                            anova(int_mod)[
                                paste0(treat_var, ":", block_var), 5
                            ] < 0.05 ~
                                "Interaction",
                            anova(par_mod)[2, 5] < 0.05 ~ "Parallel",
                            anova(treat_mod)[1, 5] < 0.05 ~ "Single",
                            TRUE ~ "None"
                        )
                } else {
                    interaction_fisher <- -2 * sum(
                        log(path_data$interaction_pval[
                            which(path_data$sub_pathway == pathway)
                        ]),
                        na.rm = TRUE
                    )

                    interaction_p <- round(
                        pchisq(interaction_fisher, 2 * length(var_names),
                            lower.tail = FALSE
                        ), 2
                    )

                    path_data$interaction_fisher[
                        which(path_data$sub_pathway == pathway)
                    ] <- interaction_p

                    parallel_fisher <- -2 * sum(log(path_data$parallel_pval[
                        which(path_data$sub_pathway == pathway)
                    ]), na.rm = TRUE)

                    parallel_p <- round(
                        pchisq(parallel_fisher, 2 * length(var_names),
                            lower.tail = FALSE
                        ), 2
                    )

                    path_data$parallel_fisher[
                        which(path_data$sub_pathway == pathway)
                    ] <- parallel_p

                    single_fisher <- -2 * sum(
                        log(path_data$single_pval[
                            which(path_data$sub_pathway == pathway)
                        ]),
                        na.rm = TRUE
                    )


                    single_p <- round(pchisq(single_fisher,
                        2 * length(var_names),
                        lower.tail = FALSE
                    ), 2)


                    path_data$single_fisher[
                        which(path_data$sub_pathway == pathway)
                    ] <- single_p

                    path_data$model[which(path_data$sub_pathway == pathway)] <-
                        case_when(
                            interaction_p < 0.05 ~ "Interaction",
                            parallel_p < 0.05 ~ "Parallel",
                            single_p < 0.05 ~ "Single",
                            TRUE ~ "None"
                        )
                }
            }
        }

        # Case when treatment only = True
        if (is.null(block_var) == TRUE) {
            # Create path data dataframe.
            path_data <- data.frame(
                CHEM_ID = rownames(chem),
                sub_pathway = chem[, subPathwayName],
                chem_name = chem[, chemName],
                super_pathway = chem[, superPathwayName],
                single_pval = NA,
                single_fisher = NA,
                model = NA
            ) %>%
                dplyr::mutate(sub_pathway = ifelse(is.na(sub_pathway),
                    "Unknown Metabolite", sub_pathway
                ))

            for (pathway in unique(path_data$sub_pathway)) {
                var_names <- path_data$CHEM_ID[
                    which(path_data$sub_pathway == pathway)
                ]


                # Calculate p values for each linear model
                for (var in var_names) {
                    outcome <- analysis[, as.character(var)]

                    if (is.null(outcome)) next


                    # Treat model
                    treat_formula <- as.formula(
                        paste("outcome",
                            paste(c(treat_var), collapse = " + "),
                            sep = "~"
                        )
                    )
                    treat_mod <- lm(treat_formula, data = analysis)

                    tryCatch(
                        {
                            path_data$single_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(treat_mod)[treat_var, 5]
                        },
                        warning = function(w) {
                            warning(w, " for CHEM ID ", var)
                        },
                        finally = {
                            path_data$single_pval[
                                which(path_data$CHEM_ID == var)
                            ] <-
                                anova(treat_mod)[treat_var, 5]
                        }
                    )
                }

                if (length(var_names) == 1) {
                    path_data$single_fisher[which(path_data$CHEM_ID == var)] <-
                        anova(treat_mod)[treat_var, 5]

                    path_data$model[which(path_data$sub_pathway == pathway)] <-
                        dplyr::case_when(
                            anova(treat_mod)[1, 5] < 0.05 ~ "Single",
                            TRUE ~ "None"
                        )
                } else {
                    single_fisher <- -2 * sum(
                        log(path_data$single_pval[
                            which(path_data$sub_pathway == pathway)
                        ]),
                        na.rm = TRUE
                    )

                    single_p <- round(pchisq(single_fisher,
                        2 * length(var_names),
                        lower.tail = FALSE
                    ), 2)

                    path_data$single_fisher[
                        which(path_data$sub_pathway == pathway)
                    ] <- single_p

                    path_data$model[which(path_data$sub_pathway == pathway)] <-
                        dplyr::case_when(
                            single_p < 0.05 ~ "Single",
                            TRUE ~ "None"
                        )
                }
            }
        }

        # Return path data
        return(path_results = path_data)
    }

    if (!is.null(strat_var)) {
        analysisFull <- analysis

        # get stratas
        strats <- unique(analysis[, strat_var])


        # Subpathresults
        sub_results <- list()

        for (i in seq_len(length(strats))) {
            # Subsetted analysis data
            analysis <- analysisFull[analysisFull[, strat_var] == strats[i], ]

            if (!is.null(block_var)) {
                # Create results dataframe
                path_data <- data.frame(
                    CHEM_ID = rownames(chem),
                    sub_pathway = chem[, subPathwayName],
                    chem_name = chem[, chemName],
                    super_pathway = chem[, superPathwayName],
                    interaction_pval = NA,
                    interaction_fisher = NA,
                    parallel_pval = NA,
                    parallel_fisher = NA,
                    single_pval = NA,
                    single_fisher = NA,
                    model = NA
                ) %>%
                    dplyr::mutate(sub_pathway = ifelse(is.na(sub_pathway),
                        "Unknown Metabolite", sub_pathway
                    ))


                # For each pathways run the models for each metabolite
                # within the pathway.
                for (pathway in unique(path_data$sub_pathway)) {
                    var_names <- path_data$CHEM_ID[
                        which(path_data$sub_pathway == pathway)
                    ]


                    # Calculate p values for each linear model
                    for (var in var_names) {
                        outcome <- analysis[, as.character(var)]

                        if (is.null(outcome)) next

                        # This Generalizes Model formulation
                        inter_formula <- as.formula(
                            paste("outcome",
                                paste(
                                    c(treat_var, block_var, paste0(
                                        treat_var,
                                        "*", block_var
                                    )),
                                    collapse = " + "
                                ),
                                sep = "~"
                            )
                        )

                        int_mod <- lm(inter_formula, data = analysis)

                        # TryCatch for Anova Modeling
                        tryCatch(
                            {
                                path_data$interaction_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(int_mod)[
                                        paste0(treat_var, ":", block_var), 5
                                    ]
                            },
                            warning = function(w) {
                                warning(w, " for CHEM ID ", var)
                            },
                            finally = {
                                path_data$interaction_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(int_mod)[
                                        paste0(treat_var, ":", block_var), 5
                                    ]
                            }
                        )

                        # Parallel model, this generalizes the model formulation
                        par_formula <- as.formula(
                            paste("outcome",
                                paste(c(treat_var, block_var),
                                    collapse = " + "
                                ),
                                sep = "~"
                            )
                        )
                        par_mod <- lm(par_formula, data = analysis)

                        # Anova modling for the parallel model.
                        tryCatch(
                            {
                                path_data$parallel_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(par_mod)[block_var, 5]
                            },
                            warning = function(w) {
                                warning(w, " for CHEM ID ", var)
                            },
                            finally = {
                                path_data$parallel_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(par_mod)[block_var, 5]
                            }
                        )

                        # Treat model
                        treat_formula <- as.formula(
                            paste("outcome",
                                paste(c(treat_var), collapse = " + "),
                                sep = "~"
                            )
                        )
                        treat_mod <- lm(treat_formula, data = analysis)

                        # Anova model for treatment.
                        tryCatch(
                            {
                                path_data$single_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(treat_mod)[treat_var, 5]
                            },
                            warning = function(w) {
                                warning(w, " for CHEM ID ", var)
                            },
                            finally = {
                                path_data$single_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(treat_mod)[treat_var, 5]
                            }
                        )
                    }

                    # Compute combined fisher probabilities for each model.
                    if (length(var_names) == 1) {
                        path_data$interaction_fisher[
                            which(path_data$CHEM_ID == var)
                        ] <-
                            anova(int_mod)[paste0(treat_var, ":", block_var), 5]

                        path_data$parallel_fisher[
                            which(path_data$CHEM_ID == var)
                        ] <-
                            anova(par_mod)[block_var, 5]

                        path_data$single_fisher[
                            which(path_data$CHEM_ID == var)
                        ] <-
                            anova(treat_mod)[treat_var, 5]

                        path_data$model[
                            which(path_data$sub_pathway == pathway)
                        ] <-
                            case_when(
                                anova(int_mod)[
                                    paste0(
                                        treat_var, ":", block_var
                                    ), 5
                                ] < 0.05 ~ "Interaction",
                                anova(par_mod)[2, 5] < 0.05 ~ "Parallel",
                                anova(treat_mod)[1, 5] < 0.05 ~ "Single",
                                TRUE ~ "None"
                            )
                    } else {
                        interaction_fisher <- -2 * sum(
                            log(path_data$interaction_pval[
                                which(path_data$sub_pathway == pathway)
                            ]),
                            na.rm = TRUE
                        )

                        interaction_p <- round(pchisq(
                            interaction_fisher, 2 * length(var_names),
                            lower.tail = FALSE
                        ), 2)

                        path_data$interaction_fisher[
                            which(path_data$sub_pathway == pathway)
                        ] <-
                            interaction_p

                        parallel_fisher <- -2 * sum(
                            log(path_data$parallel_pval[
                                which(path_data$sub_pathway == pathway)
                            ]),
                            na.rm = TRUE
                        )

                        parallel_p <- round(
                            pchisq(parallel_fisher,
                                2 * length(var_names),
                                lower.tail = FALSE
                            ), 2
                        )

                        path_data$parallel_fisher[
                            which(path_data$sub_pathway == pathway)
                        ] <- parallel_p

                        single_fisher <- -2 * sum(
                            log(path_data$single_pval[
                                which(path_data$sub_pathway == pathway)
                            ]),
                            na.rm = TRUE
                        )

                        single_p <- round(pchisq(single_fisher,
                            2 * length(var_names),
                            lower.tail = FALSE
                        ), 2)

                        path_data$single_fisher[which(
                            path_data$sub_pathway == pathway
                        )] <- single_p

                        path_data$model[which(
                            path_data$sub_pathway == pathway
                        )] <- case_when(
                            interaction_p < 0.05 ~ "Interaction",
                            parallel_p < 0.05 ~ "Parallel",
                            single_p < 0.05 ~ "Single",
                            TRUE ~ "None"
                        )
                    }
                }
            }


            if (is.null(block_var) == TRUE) {
                # Create path data dataframe.
                path_data <- data.frame(
                    CHEM_ID = rownames(chem),
                    sub_pathway = chem[, subPathwayName],
                    chem_name = chem[, chemName],
                    super_pathway = chem[, superPathwayName],
                    single_pval = NA,
                    single_fisher = NA,
                    model = NA
                ) %>%
                    dplyr::mutate(sub_pathway = ifelse(is.na(sub_pathway),
                        "Unknown Metabolite", sub_pathway
                    ))

                for (pathway in unique(path_data$sub_pathway)) {
                    var_names <- path_data$CHEM_ID[
                        which(path_data$sub_pathway == pathway)
                    ]


                    # Calculate p values for each linear model
                    for (var in var_names) {
                        outcome <- analysis[, as.character(var)]

                        if (is.null(outcome)) next


                        # Treat model
                        treat_formula <- as.formula(
                            paste("outcome",
                                paste(c(treat_var), collapse = " + "),
                                sep = "~"
                            )
                        )
                        treat_mod <- lm(treat_formula, data = analysis)

                        tryCatch(
                            {
                                path_data$single_pval[which(
                                    path_data$CHEM_ID == var
                                )] <-
                                    anova(treat_mod)[treat_var, 5]
                            },
                            warning = function(w) {
                                warning(w, " for CHEM ID ", var)
                            },
                            finally = {
                                path_data$single_pval[
                                    which(path_data$CHEM_ID == var)
                                ] <-
                                    anova(treat_mod)[treat_var, 5]
                            }
                        )
                    }

                    if (length(var_names) == 1) {
                        path_data$single_fisher[
                            which(path_data$CHEM_ID == var)
                        ] <-
                            anova(treat_mod)[treat_var, 5]

                        path_data$model[
                            which(path_data$sub_pathway == pathway)
                        ] <-
                            dplyr::case_when(
                                anova(treat_mod)[1, 5] < 0.05 ~ "Single",
                                TRUE ~ "None"
                            )
                    } else {
                        single_fisher <- -2 * sum(
                            log(path_data$single_pval[
                                which(path_data$sub_pathway == pathway)
                            ]),
                            na.rm = TRUE
                        )

                        single_p <- round(
                            pchisq(single_fisher, 2 * length(var_names),
                                lower.tail = FALSE
                            ), 2
                        )

                        path_data$single_fisher[
                            which(path_data$sub_pathway == pathway)
                        ] <- single_p

                        path_data$model[
                            which(path_data$sub_pathway == pathway)
                        ] <-
                            dplyr::case_when(
                                single_p < 0.05 ~ "Single",
                                TRUE ~ "None"
                            )
                    }
                }
            }

            sub_results[[strats[i]]] <- path_data
        }

        return(path_results = sub_results)
    }
}
