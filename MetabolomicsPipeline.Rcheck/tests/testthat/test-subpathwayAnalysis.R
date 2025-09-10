data("demoDat")


testDat <- demoDat[rowData(demoDat)$SUB_PATHWAY == "Phospholipid Metabolism" & !is.na(rowData(demoDat)$SUB_PATHWAY), ]

# Subpathway analysis without stratification or block
singleMod <- subpathwayAnalysis(testDat,
    treat_var = "GROUP_NAME",
    Assay = "normalized"
)



# create data for test
metaboliteDat <- t(assay(testDat, "normalized"))

metaData <- colData(testDat)

modData <- cbind(metaData, metaboliteDat)


trueModel <- anova(lm(modData$`267` ~ GROUP_NAME, data = modData), lm(modData$`267` ~ 1, modData))


test_that("Testing pvalues for single mod", {
    expect_equal(trueModel$`Pr(>F)`[2], singleMod[singleMod$CHEM_ID == "267", "single_pval"])
})


# test subpathway pvals
xtilde <- -2 * sum(log(singleMod$single_pval))




test_that("Testing subpathway for single mod", {
    expect_equal(round(pchisq(xtilde, 14, lower.tail = F), 2), singleMod[1, "single_fisher"])
})


###############################################################################
## Test with block variable
##############################################################################
blockMod <- subpathwayAnalysis(testDat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    Assay = "normalized"
)




# Test pvalues for chem_ID 267
TrueInter <- anova(lm(modData$`267` ~ GROUP_NAME * TIME1, data = modData), lm(modData$`267` ~ GROUP_NAME + TIME1, modData))
TrueParallel_block <- anova(lm(modData$`267` ~ GROUP_NAME + TIME1, data = modData), lm(modData$`267` ~ GROUP_NAME, modData))
TrueParallel_treatment <- anova(lm(modData$`267` ~ TIME1+GROUP_NAME, data = modData), lm(modData$`267` ~ TIME1, modData))
TrueSingle <- anova(lm(modData$`267` ~ GROUP_NAME, data = modData), lm(modData$`267` ~ 1, modData))


test_that("Test pvalues for block mod on chem ID 267", {
    expect_equal(TrueInter$`Pr(>F)`[2], blockMod[blockMod$CHEM_ID == "267", "interaction_pval"])
    expect_equal(TrueParallel_block$`Pr(>F)`[2], blockMod[blockMod$CHEM_ID == "267", "parallel_block_pval"])
    expect_equal(TrueParallel_treatment$`Pr(>F)`[2], blockMod[blockMod$CHEM_ID == "267", "parallel_treatment_pval"])
    expect_equal(TrueSingle$`Pr(>F)`[2], blockMod[blockMod$CHEM_ID == "267", "single_pval"])
})



InterFisher <- round(pchisq(-2 * sum(log(blockMod$interaction_pval)), 14, lower.tail = F), 2)
ParaFisher_block <- round(pchisq(-2 * sum(log(blockMod$parallel_block_pval)), 14, lower.tail = F), 2)
ParaFisher_treatment <- round(pchisq(-2 * sum(log(blockMod$parallel_treatment_pval)), 14, lower.tail = F), 2)
SingleFisher <- round(pchisq(-2 * sum(log(blockMod$single_pval)), 14, lower.tail = F), 2)

test_that("Test Fisher pvalues for block mod on chem ID 267", {
    expect_equal(InterFisher, blockMod[blockMod$CHEM_ID == "267", "interaction_fisher"])
    expect_equal(ParaFisher_block, blockMod[blockMod$CHEM_ID == "267", "parallel_block_fisher"])
    expect_equal(ParaFisher_treatment, blockMod[blockMod$CHEM_ID == "267", "parallel_treatment_fisher"])
    expect_equal(SingleFisher, blockMod[blockMod$CHEM_ID == "267", "single_fisher"])
})

###############################################################################
## Test with block variable and stratified valieable
##############################################################################
blockMod <- subpathwayAnalysis(testDat,
    treat_var = "GROUP_NAME",
    strat_var = "Gender",
    block_var = "TIME1",
    Assay = "normalized"
)




# Test pvalues for chem_ID 267 and gender ==Female
TrueInter <- anova(lm(modData$`267` ~ GROUP_NAME * TIME1, data = modData, subset = Gender == "Female"), lm(modData$`267` ~ GROUP_NAME + TIME1, modData, subset = Gender == "Female"))
TrueParallel_block <- anova(lm(modData$`267` ~ GROUP_NAME + TIME1, data = modData, subset = Gender == "Female"), lm(modData$`267` ~ GROUP_NAME, modData, subset = Gender == "Female"))
TrueParallel_treatment <- anova(lm(modData$`267` ~ TIME1 +GROUP_NAME, data = modData, subset = Gender == "Female"), lm(modData$`267` ~ TIME1, modData, subset = Gender == "Female"))
TrueSingle <- anova(lm(modData$`267` ~ GROUP_NAME, data = modData, subset = Gender == "Female"), lm(modData$`267` ~ 1, modData, subset = Gender == "Female"))


test_that("Test pvalues for block mod on chem ID 267", {
    expect_equal(TrueInter$`Pr(>F)`[2], blockMod$Female[blockMod$Female$CHEM_ID == "267", "interaction_pval"])
    expect_equal(TrueParallel_block$`Pr(>F)`[2], blockMod$Female[blockMod$Female$CHEM_ID == "267", "parallel_block_pval"])
    expect_equal(TrueParallel_treatment$`Pr(>F)`[2], blockMod$Female[blockMod$Female$CHEM_ID == "267", "parallel_treatment_pval"])
    expect_equal(TrueSingle$`Pr(>F)`[2], blockMod$Female[blockMod$Female$CHEM_ID == "267", "single_pval"])
})



InterFisher <- round(pchisq(-2 * sum(log(blockMod$Female$interaction_pval)), 14, lower.tail = F), 2)
ParaFisher_block <- round(pchisq(-2 * sum(log(blockMod$Female$parallel_block_pval)), 14, lower.tail = F), 2)
ParaFisher_treatment <- round(pchisq(-2 * sum(log(blockMod$Female$parallel_treatment_pval)), 14, lower.tail = F), 2)
SingleFisher <- round(pchisq(-2 * sum(log(blockMod$Female$single_pval)), 14, lower.tail = F), 2)

test_that("Test Fisher pvalues for block mod on chem ID 267", {
    expect_equal(InterFisher, blockMod$Female[blockMod$Female$CHEM_ID == "267", "interaction_fisher"])
    expect_equal(ParaFisher_block, blockMod$Female[blockMod$Female$CHEM_ID == "267", "parallel_block_fisher"])
    expect_equal(ParaFisher_treatment, blockMod$Female[blockMod$Female$CHEM_ID == "267", "parallel_treatment_fisher"])
    expect_equal(SingleFisher, blockMod$Female[blockMod$Female$CHEM_ID == "267", "single_fisher"])
})


###############################################################################
## Test with no block variable and stratified valieable
##############################################################################
singleMod <- subpathwayAnalysis(testDat,
                               treat_var = "GROUP_NAME",
                               strat_var = "Gender",
                               block_var = NULL,
                               Assay = "normalized"
)

trueModel <- anova(lm(modData$`267` ~ GROUP_NAME, data = modData, subset = Gender=="Female"), lm(modData$`267` ~ 1, modData,subset = Gender=="Female"))

test_that("Testing pvalues for single mod with stratified variable", {
    expect_equal(trueModel$`Pr(>F)`[2], singleMod$Female[singleMod$Female$CHEM_ID == "267", "single_pval"])
})

# test subpathway pvals
xtilde <- -2 * sum(log(singleMod$Female$single_pval))

test_that("Testing subpathway for single mod with stratified variable", {
    expect_equal(round(pchisq(xtilde, 14, lower.tail = F), 2), singleMod$Female[singleMod$Female$CHEM_ID == "267", "single_fisher"])
})

