#' Pairwise function
#'
#' This is the main function for metabolite_pairwise
#'
#' @param out Outcome used as reponse
#'
#' @param form form of the model
#'
#' @param data data used for modeling
#' 
#' @param adjust TRUE/FALSE for whether or not to adjust p-values
#' @returns Pairwise comparisons for a single metabolite.
#'
#'
#' @importFrom emmeans emmeans
#' @importFrom stats anova

pairwise <- function(out, form, data, adjust = NULL) {
    # Define model
    model <- as.formula(
        paste("outcome",
            paste(form, collapse = " + "),
            sep = "~"
        )
    )

    # define outcome
    outcome <- out

    # Run model
    mod <- lm(model, data = data)

    mod2 <- lm(outcome ~ 1, data = data)

    anov <- stats::anova(mod, mod2)

    # Get the F statistic pvalue
    overall <- anov$`Pr(>F)`[2]


    # Create model for the pairwise comparisons
    means_mod <- as.formula(
        paste("pairwise",
            paste(form, collapse = " + "),
            sep = "~"
        )
    )


    # Find pairwise comparisons
    pairs <- emmeans::emmeans(mod, means_mod, adjust = adjust)$contrasts



    cols <- c(
        "Overall_pval", paste0(summary(pairs)[, 1], "_ESTS"),
        paste0(summary(pairs)[, 1], "_PVALS")
    )

    results <- data.frame(matrix(nrow = 1, ncol = length(cols)))

    colnames(results) <- cols

    results[1, "Overall_pval"] <- overall

    results[1, paste0(summary(pairs)[, 1], "_ESTS")] <-
       summary(pairs)[, "estimate"]

    results[1, paste0(summary(pairs)[, 1], "_PVALS")] <-
       summary(pairs)[, "p.value"]

    return(results)
}
