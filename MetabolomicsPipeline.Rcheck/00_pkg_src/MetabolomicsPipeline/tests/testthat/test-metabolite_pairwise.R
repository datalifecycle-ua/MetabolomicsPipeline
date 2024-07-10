data("demoDat")

# Subset for a single metabolite
testDat <- demoDat[c("30", "35"), ]


# run metaboloite pairwise with strat variable
strat <- metabolite_pairwise(testDat, form = "GROUP_NAME*TIME1", strat_var = "Gender")


# test with metabolite 35
met <- assay(testDat, "normalized")["35", ]

meta <- colData(demoDat)

meta$metabolite <- met

modFemale <- lm(metabolite ~ GROUP_NAME * TIME1, data = meta, subset = Gender == "Female")
modMale <- lm(metabolite ~ GROUP_NAME * TIME1, data = meta, subset = Gender == "Male")

pairsFemales <- emmeans::emmeans(modFemale, "pairwise" ~ GROUP_NAME * TIME1, adjust = "none")$contrasts
pairsMales <- emmeans::emmeans(modMale, "pairwise" ~ GROUP_NAME * TIME1, adjust = "none")$contrasts


FemalesSummary <- data.frame(summary(pairsFemales))
MalesSummary <- data.frame(summary(pairsMales))



# Test stratas for estimates
test_that("Test Strata one has the same estimates", {
    expect_equal(round(FemalesSummary$estimate, 4), as.numeric(round(strat$Female["35", paste0(FemalesSummary$contrast, "_ESTS")], 4)))
})

test_that("Test Strata two has the same estimates", {
    expect_equal(round(MalesSummary$estimate, 4), as.numeric(round(strat$Male["35", paste0(MalesSummary$contrast, "_ESTS")], 4)))
})


# Test strata pvalues
test_that("Test Strata one has the same pvalues", {
    expect_equal(round(FemalesSummary$p.value, 4), as.numeric(round(strat$Female["35", paste0(FemalesSummary$contrast, "_PVALS")], 4)))
})

test_that("Test Strata two has the same estimates", {
    expect_equal(round(MalesSummary$p.value, 4), as.numeric(round(strat$Male["35", paste0(MalesSummary$contrast, "_PVALS")], 4)))
})


##############################################################################
## Test non stratified comparison
##############################################################################
metPairwise <- metabolite_pairwise(testDat, form = "GROUP_NAME*TIME1")

mod <- lm(metabolite ~ GROUP_NAME * TIME1, data = meta)
pairs <- emmeans::emmeans(mod, "pairwise" ~ GROUP_NAME * TIME1, adjust = "none")$contrasts

pairsSummary <- data.frame(summary(pairs))

trueEsts <- pairsSummary$estimate
testEsts <- metPairwise["35", paste0(pairsSummary$contrast, "_ESTS")]

truePvals <- pairsSummary$p.value
testPvals <- metPairwise["35", paste0(pairsSummary$contrast, "_PVALS")]

test_that("Non stratified estimates are correct", {
    expect_equal(as.numeric(truePvals), as.numeric(testPvals))
})
