# =============================================================================
# MORE PERMISSIVE, MORE POPULIST?
# Replication script
#
# Study 1: cross-national multilevel logit, CSES Module 5
# Study 2: within-respondent comparison of two simultaneous Swiss ballots
#
# Input:  data/cses5.rdata (CSES Module 5 Full Release, 25 July 2023)
# Output: figures and tables in output/
# Run with project/ as the working directory.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------
library(lme4)       # multilevel logits
library(ggplot2)    # figures
library(dplyr)      # data prep
library(scales)     # percent axes
library(survival)   # clogit(), Study 2 fixed-effects logit

Sys.setlocale("LC_ALL", "C.UTF-8")

CSES_PATH <- file.path("data", "cses5.rdata")
stopifnot(file.exists(CSES_PATH))
OUT <- "output"
dir.create(OUT, showWarnings = FALSE)

# Caches the slow glmer fits to disk; delete output/cache/ to force a refit.
CACHE <- file.path(OUT, "cache")
dir.create(CACHE, showWarnings = FALSE)
cache <- function(name, expr) {
  f <- file.path(CACHE, paste0(name, ".rds"))
  if (file.exists(f)) { cat("[cache] reusing", name, "\n"); return(readRDS(f)) }
  value <- expr
  saveRDS(value, f)
  value
}

# Fixes the draws behind the predicted probabilities reported in the paper.
set.seed(20191020)         # date of the 2019 Swiss federal election

# -----------------------------------------------------------------------------
# 2. IMPORT CSES MODULE 5
# -----------------------------------------------------------------------------
load(CSES_PATH)
d <- cses5
cat("CSES Module 5 loaded:", nrow(d), "respondents,", ncol(d), "variables\n")

# -----------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# -----------------------------------------------------------------------------
# CSES codes refusals/don't-knows as out-of-range values (7/8/9 on 1-5 items,
# 97/98/99 on 0-10 items, 999999 on party codes); every recode below clears
# these to NA first.

# Values above the ceiling are missing codes.
valid_up_to <- function(x, ceiling)
  ifelse(!is.na(x) & x <= ceiling, x, NA)

# Flips a 5-point agree/disagree item so high = anti-elite.
reverse_5 <- function(x)
  ifelse(!is.na(x) & x <= 5, 6 - x, NA)

z_score <- function(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Cronbach's alpha on complete cases.
cronbach_alpha <- function(m) {
  m <- m[complete.cases(m), ]
  k <- ncol(m)
  (k / (k - 1)) * (1 - sum(apply(m, 2, var)) / var(rowSums(m)))
}

stars <- function(p)
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*", ifelse(p < 0.1, "+", ""))))

# Adds a p-value column when summary() omits one.
coef_matrix <- function(model) {
  s <- summary(model)$coefficients
  if (ncol(s) < 4) s <- cbind(s, "Pr(>|z|)" = 2 * pnorm(-abs(s[, 3])))
  s
}

# Coefficient table with 95% Wald intervals, ready for export.
fmt_model <- function(model) {
  s <- coef_matrix(model)
  data.frame(term  = rownames(s),
             b     = round(s[, 1], 3),
             se    = round(s[, 2], 3),
             p     = s[, 4],
             star  = stars(s[, 4]),
             ci_lo = round(s[, 1] - 1.96 * s[, 2], 3),
             ci_hi = round(s[, 1] + 1.96 * s[, 2], 3),
             row.names = NULL)
}

# One row of the robustness table: a labelled interaction estimate.
interaction_row <- function(label, b, se, p, n) {
  data.frame(specification = label,
             b = round(b, 3), se = round(se, 3),
             ci = sprintf("[%.3f, %.3f]", b - 1.96 * se, b + 1.96 * se),
             p = round(p, 3), n = n)
}

# -----------------------------------------------------------------------------
# 4. STUDY 1: DEPENDENT VARIABLE (POPULIST VOTE)
# -----------------------------------------------------------------------------
# Lower-house vote: list ballot where the system has one, candidate ballot
# otherwise (list-only would drop every single-member-district country).
d$vote_list      <- ifelse(d$E3013_LH_PL < 900000, d$E3013_LH_PL, NA)
d$vote_candidate <- ifelse(d$E3013_LH_DC < 900000, d$E3013_LH_DC, NA)
d$vote_party     <- ifelse(!is.na(d$vote_list), d$vote_list, d$vote_candidate)
d$vote_source    <- ifelse(!is.na(d$vote_list), "list ballot",
                    ifelse(!is.na(d$vote_candidate), "candidate ballot", NA))

# Match the chosen party to its expert populism rating (up to nine parties
# per study: identifiers E5000_A..I, ratings E5020_A..I).
d$party_populism <- NA_real_
for (letter in LETTERS[1:9]) {
  party_code   <- d[[paste0("E5000_", letter)]]
  expert_score <- valid_up_to(d[[paste0("E5020_", letter)]], 10)
  matched <- !is.na(d$vote_party) & !is.na(party_code) &
             party_code < 900000 & d$vote_party == party_code
  d$party_populism[matched] <- expert_score[matched]
}

# -----------------------------------------------------------------------------
# 5. STUDY 1: INDIVIDUAL-LEVEL PREDICTORS
# -----------------------------------------------------------------------------
# Five items from the CSES anti-elite battery, recoded so high = more
# anti-elite. Item 6 ("the people, not politicians, should decide") is
# excluded: in a direct-democracy country it taps support for an existing
# institution rather than anti-elite resentment (see alpha check below).
anti_elite_items <- cbind(
  compromise_is_selling_out = reverse_5(d$E3004_1),
  politicians_dont_care     = reverse_5(d$E3004_2),
  politicians_trustworthy   = valid_up_to(d$E3004_3, 5),
  politicians_main_problem  = reverse_5(d$E3004_4),
  serve_rich_and_powerful   = reverse_5(d$E3004_7)
)

# At least three answered items are required for the average to be meaningful.
items_answered <- rowSums(!is.na(anti_elite_items))
d$anti_elite_raw <- ifelse(items_answered >= 3,
                           rowMeans(anti_elite_items, na.rm = TRUE), NA)

alpha_pooled <- cronbach_alpha(anti_elite_items)
cat("Cronbach's alpha, pooled five-item scale:", round(alpha_pooled, 3), "\n")

# Alpha by election study, as a diagnostic for cross-national comparability.
alpha_by_study <- sapply(split(seq_len(nrow(d)), d$E1004), function(i) {
  m <- anti_elite_items[i, , drop = FALSE]
  m <- m[complete.cases(m), , drop = FALSE]
  if (nrow(m) < 50) return(NA_real_)
  cronbach_alpha(m)
})

# Controls, recoded from their CSES codes.
d$age    <- valid_up_to(d$E2001_A, 120)
d$female <- ifelse(d$E2002 %in% c(1, 2), d$E2002 - 1, NA)
d$edu    <- valid_up_to(d$E2003, 9)
d$lr     <- valid_up_to(d$E3020, 10)

# Satisfaction with democracy is coded 1, 2, 4, 5 (3 unused); mapped to a
# clean 1-4 scale and flipped so that high = more satisfied.
sat_raw  <- ifelse(d$E3023 %in% c(1, 2, 4, 5), d$E3023, NA)
d$satdem <- 5 - c(1, 2, NA, 3, 4)[sat_raw]

# -----------------------------------------------------------------------------
# 6. STUDY 1: STUDY-LEVEL PREDICTORS
# -----------------------------------------------------------------------------
# District magnitude (E4001), falling back to the nationwide figure (E4001_N)
# for single-national-district systems.
d$magnitude_ind <- ifelse(d$E4001 < 990, d$E4001,
                   ifelse(d$E4001_N < 990, d$E4001_N, NA))

# Fix: for the two US studies, E4001 records electoral college votes, not
# House-district magnitude. House districts are single-member, so this is 1.
US_STUDIES <- c("USA_2016", "USA_2020")
d$magnitude_ind[d$E1004 %in% US_STUDIES] <- 1

# Main moderator: log of the study's median district magnitude (median:
# magnitude is right-skewed within countries; log: diminishing marginal
# effect of added seats).
study_magnitude <- d %>%
  filter(!is.na(magnitude_ind)) %>%
  group_by(E1004) %>%
  summarise(magnitude_med = median(magnitude_ind), .groups = "drop")
d <- left_join(d, study_magnitude, by = "E1004")
d$log_magnitude <- log(d$magnitude_med)

# Study-level control: any party rated 7+ on populism held cabinet
# portfolios before the election (populists in office attract votes for
# reasons unrelated to anti-elite protest, cf. Jungkunz et al. 2021).
study_rows <- d[!duplicated(d$E1004), ]
populist_in_gov <- sapply(seq_len(nrow(study_rows)), function(i) {
  scores     <- valid_up_to(unlist(study_rows[i, paste0("E5020_", LETTERS[1:9])]), 10)
  portfolios <- valid_up_to(unlist(study_rows[i, paste0("E5011_", LETTERS[1:9])]), 90)
  as.integer(any(!is.na(scores) & scores >= 7 &
                 !is.na(portfolios) & portfolios > 0))
})
d <- left_join(d, data.frame(E1004 = study_rows$E1004,
                             populist_in_gov = populist_in_gov), by = "E1004")

# Drop unused CSES columns before the sample gets rebuilt repeatedly across
# the robustness refits (keeps memory use manageable).
d <- d[, c("E1004", "vote_source", "party_populism", "anti_elite_raw",
           "age", "female", "edu", "lr", "satdem",
           "magnitude_ind", "magnitude_med", "log_magnitude", "populist_in_gov",
           "E3013_LH_PL", "E3013_UH_DC_1", "E2021")]
rm(cses5, anti_elite_items, study_rows, study_magnitude)
gc(verbose = FALSE)

# -----------------------------------------------------------------------------
# 7. STUDY 1: ANALYTIC SAMPLE
# -----------------------------------------------------------------------------
# Retention rule (post listwise deletion): >=300 respondents and >=5 populist
# votes per study, so each study-specific intercept/slope is estimable.
# Wrapped in a function so the threshold robustness checks (section 10b)
# reuse the same sample definition.
make_study1_sample <- function(pop_threshold = 7, min_n = 300, min_pop = 5) {
  dd <- d
  dd$populist_vote <- ifelse(!is.na(dd$party_populism),
                             as.integer(dd$party_populism >= pop_threshold), NA)
  model_vars <- c("populist_vote", "anti_elite_raw", "age", "female",
                  "edu", "lr", "satdem", "log_magnitude", "populist_in_gov")
  s <- dd[complete.cases(dd[, model_vars]), ]

  keep <- s %>% group_by(E1004) %>%
    summarise(n = n(), n_pop = sum(populist_vote), .groups = "drop") %>%
    filter(n >= min_n, n_pop >= min_pop)
  s <- s %>% filter(E1004 %in% keep$E1004)

  # Standardised/centred within the analytic sample.
  s$anti_elite <- z_score(s$anti_elite_raw)
  s$age_c      <- (s$age - mean(s$age)) / 10   # units of ten years
  s$edu_c      <- s$edu    - mean(s$edu)
  s$lr_c       <- s$lr     - mean(s$lr)
  s$satdem_c   <- s$satdem - mean(s$satdem)
  s$logmag_c   <- s$log_magnitude - mean(s$log_magnitude)
  s
}

s1 <- make_study1_sample()
cat("Study 1 analytic sample:", nrow(s1), "respondents in",
    length(unique(s1$E1004)), "election studies\n")

# Reliability inside the retained studies.
alpha_retained <- alpha_by_study[names(alpha_by_study) %in% unique(s1$E1004)]
cat("Alpha by study: min", round(min(alpha_retained, na.rm = TRUE), 3),
    "median", round(median(alpha_retained, na.rm = TRUE), 3),
    "max", round(max(alpha_retained, na.rm = TRUE), 3), "\n")
write.csv(data.frame(study = names(alpha_retained),
                     alpha = round(as.vector(alpha_retained), 3)),
          file.path(OUT, "table_alpha_by_study.csv"), row.names = FALSE)

# Descriptive table of the retained studies (Table 1 of the article).
study_table <- s1 %>% group_by(E1004) %>% summarise(
  n = n(),
  populist_share = round(100 * mean(populist_vote), 1),
  magnitude      = median(magnitude_med),
  ballot         = names(which.max(table(vote_source))),
  pop_in_gov     = first(populist_in_gov),
  .groups = "drop")
write.csv(study_table, file.path(OUT, "table_studies.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. STUDY 1: MAIN MODEL
# -----------------------------------------------------------------------------
# Random intercept and random anti-elite slope by election study; without
# the random slope, the cross-level interaction would be tested against a
# standard error that assumes one common slope. anti_elite:logmag_c is H1.
m1 <- cache("m1", glmer(
  populist_vote ~ anti_elite * logmag_c +
    age_c + female + edu_c + lr_c + satdem_c + populist_in_gov +
    (1 + anti_elite | E1004),
  data = s1, family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 200000))))
print(summary(m1))
write.csv(fmt_model(m1), file.path(OUT, "table_model1.csv"), row.names = FALSE)

m1_coef  <- coef_matrix(m1)
m1_Sigma <- VarCorr(m1)$E1004
write.csv(data.frame(
  n_obs     = nobs(m1),
  n_studies = ngrps(m1),
  aic       = round(AIC(m1), 1),
  var_int   = round(as.data.frame(VarCorr(m1))$vcov[1], 3),
  var_slope = round(as.data.frame(VarCorr(m1))$vcov[2], 3)),
  file.path(OUT, "table_model1_fit.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. STUDY 1: PREDICTED PROBABILITIES AND FIGURE 1
# -----------------------------------------------------------------------------
# Population-averaged predicted probabilities via parametric simulation:
# coefficient vectors drawn from their sampling distribution (chol()
# preserves the coefficient correlations), random effects integrated out by
# averaging over draws. Describes a randomly drawn voter, not a median study.
predict_marginal <- function(model, newdata, formula_rhs, re_draws,
                             n_beta = 500, n_re = 500) {
  X        <- model.matrix(formula_rhs, newdata)
  beta_hat <- fixef(model)
  X        <- X[, names(beta_hat), drop = FALSE]

  betas <- matrix(rnorm(n_beta * length(beta_hat)), n_beta) %*%
             chol(as.matrix(vcov(model)))
  betas <- sweep(betas, 2, beta_hat, "+")

  eta <- X %*% t(betas)          # grid points x coefficient draws
  U   <- re_draws(newdata, n_re) # grid points x random-effect draws

  out <- t(sapply(seq_len(nrow(X)), function(i) {
    p <- colMeans(plogis(outer(U[i, ], eta[i, ], "+")))
    c(mean(p), unname(quantile(p, 0.025)), unname(quantile(p, 0.975)))
  }))
  data.frame(newdata, fit = out[, 1], lo = out[, 2], hi = out[, 3])
}

# Correlated intercept/slope deviations, drawn jointly; slope deviation
# scaled by the attitude value at each grid point.
re_draws_study1 <- function(nd, n) {
  Z <- matrix(rnorm(n * 2), n, 2) %*% chol(m1_Sigma)
  t(sapply(nd$anti_elite, function(x) Z[, 1] + Z[, 2] * x))
}

magnitudes <- c(1, 5, 15, 40)
grid1 <- expand.grid(anti_elite = seq(-2, 2, by = 0.1),
                     magnitude_shown = magnitudes)
grid1$logmag_c        <- log(grid1$magnitude_shown) - mean(s1$log_magnitude)
grid1$age_c           <- 0
grid1$female          <- mean(s1$female)
grid1$edu_c           <- 0
grid1$lr_c            <- 0
grid1$satdem_c        <- 0
grid1$populist_in_gov <- mean(s1$populist_in_gov)

pred1 <- predict_marginal(m1, grid1,
  ~ anti_elite * logmag_c + age_c + female + edu_c + lr_c +
    satdem_c + populist_in_gov, re_draws_study1)
pred1$magnitude_label <- factor(
  pred1$magnitude_shown, levels = magnitudes,
  labels = paste0(magnitudes, ifelse(magnitudes == 1, " seat", " seats")))

fig1 <- ggplot(pred1, aes(anti_elite, fit,
                          colour = magnitude_label, fill = magnitude_label)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_colour_brewer(palette = "Dark2", name = "District magnitude") +
  scale_fill_brewer(palette = "Dark2", name = "District magnitude") +
  labs(x = "Anti-elite attitudes (standard deviations from the mean)",
       y = "Predicted probability of a populist vote") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig1_crossnational.png"), fig1,
       width = 6.5, height = 4.2, dpi = 300)
write.csv(pred1[, c("anti_elite", "magnitude_shown", "fit", "lo", "hi")],
          file.path(OUT, "pred_study1.csv"), row.names = FALSE)

# Marginal effect of one SD of attitudes, in pp, at each magnitude. Quoted
# in the text.
me1 <- pred1 %>% filter(anti_elite %in% c(-0.5, 0.5)) %>%
  group_by(magnitude_shown) %>%
  summarise(effect_pp = round(100 * (fit[anti_elite == 0.5] -
                                     fit[anti_elite == -0.5]), 1),
            .groups = "drop")
write.csv(me1, file.path(OUT, "table_marginal_effects1.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. STUDY 1: ROBUSTNESS
# -----------------------------------------------------------------------------
# (a) Magnitude at the respondent's own district rather than the study
#     median: uses within- as well as between-country variation, at the
#     cost of no longer being a cross-level interaction.
s1$logmag_ind_c <- log(s1$magnitude_ind) -
                   mean(log(s1$magnitude_ind), na.rm = TRUE)
m1_rob <- cache("m1_rob", glmer(
  populist_vote ~ anti_elite * logmag_ind_c +
    age_c + female + edu_c + lr_c + satdem_c + populist_in_gov +
    (1 + anti_elite | E1004),
  data = s1, family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 200000))))
write.csv(fmt_model(m1_rob), file.path(OUT, "table_model1_robust.csv"),
          row.names = FALSE)
rob <- coef_matrix(m1_rob)["anti_elite:logmag_ind_c", ]
row_district <- interaction_row("H1, respondent's own district magnitude",
                                rob[1], rob[2], rob[4], nobs(m1_rob))
rm(m1_rob); gc(verbose = FALSE)


# (b) Alternative populist-vote thresholds: the 7+ cut is a convention;
#     moving it changes which votes count as populist and, via the
#     retention rule, which studies survive.
fit_h1 <- function(dat) glmer(
  populist_vote ~ anti_elite * logmag_c +
    age_c + female + edu_c + lr_c + satdem_c + populist_in_gov +
    (1 + anti_elite | E1004),
  data = dat, family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))

run_threshold <- function(label, threshold, cache_name) {
  cache(cache_name, {
    dat <- make_study1_sample(pop_threshold = threshold)
    cat(label, ":", length(unique(dat$E1004)), "studies,", nrow(dat),
        "respondents\n")
    fit <- fit_h1(dat)
    cf  <- coef_matrix(fit)["anti_elite:logmag_c", ]
    row <- interaction_row(label, cf[1], cf[2], cf[4], nobs(fit))
    rm(fit, dat); gc(verbose = FALSE)
    row
  })
}
row_t6 <- run_threshold("H1, populist threshold at 6", 6, "row_t6")
row_t8 <- run_threshold("H1, populist threshold at 8", 8, "row_t8")

# -----------------------------------------------------------------------------
# 11. STUDY 2: SWISS SAMPLE AND LONG-FORMAT RESHAPE
# -----------------------------------------------------------------------------
# Switzerland elects both chambers the same day: National Council by PR
# (canton as district), Council of States by majority rule except in Jura
# and Neuchatel. Same voters, attitudes measured once, complete both ballots.
che <- d[d$E1004 == "CHE_2019", ]

che$vote_nc <- ifelse(che$E3013_LH_PL   < 900000, che$E3013_LH_PL, NA)
# Upper-house outcome: the first named Council of States candidate, the
# choice comparable to the single National Council list vote.
che$vote_cs <- ifelse(che$E3013_UH_DC_1 < 900000, che$E3013_UH_DC_1, NA)
che$canton  <- che$E2021   # official FSO numbering, 1 = Zurich .. 26 = Jura

# Cantons with a single National Council seat are excluded: their "PR"
# election is a plurality contest in practice, so the contrast disappears.
SINGLE_SEAT   <- c(4, 6, 7, 8, 15, 16)   # UR, OW, NW, GL, AR, AI
SVP_CODE      <- 756001                  # Swiss People's Party
PR_CANTONS_CS <- c(24, 26)               # Neuchatel and Jura: upper house
                                         # also elected proportionally

che <- che %>% filter(!canton %in% SINGLE_SEAT,
                      !is.na(vote_nc), !is.na(vote_cs))
che <- che[complete.cases(che[, c("anti_elite_raw", "age", "female",
                                  "edu", "lr", "satdem")]), ]
che$respondent_id <- seq_len(nrow(che))

# Standardised within the Swiss sample.
che$anti_elite <- z_score(che$anti_elite_raw)
che$age_c      <- (che$age - mean(che$age)) / 10
che$edu_c      <- che$edu    - mean(che$edu)
che$lr_c       <- che$lr     - mean(che$lr)
che$satdem_c   <- che$satdem - mean(che$satdem)

# Long format: two rows per respondent, one per ballot. The repeated
# respondent identifier lets a random intercept absorb everything stable
# about the person.
base_cols <- c("respondent_id", "canton", "anti_elite", "age_c", "female",
               "edu_c", "lr_c", "satdem_c")
long <- rbind(
  transform(che[, base_cols], chamber = "National Council",
            svp_vote = as.integer(che$vote_nc == SVP_CODE)),
  transform(che[, base_cols], chamber = "Council of States",
            svp_vote = as.integer(che$vote_cs == SVP_CODE)))

# The treatment indicator is the proportional ballot, not the chamber: in
# Neuchatel and Jura the Council of States ballot is proportional too.
long$pr <- as.integer(
  long$chamber == "National Council" |
  (long$chamber == "Council of States" & long$canton %in% PR_CANTONS_CS))

cat("Study 2 analytic sample:", length(unique(long$respondent_id)),
    "respondents,", nrow(long), "vote observations,",
    length(unique(long$canton)), "cantons\n")

sw_desc <- long %>% group_by(chamber) %>%
  summarise(n = n(), svp = sum(svp_vote),
            share = round(100 * mean(svp_vote), 1), .groups = "drop")
write.csv(sw_desc, file.path(OUT, "table_swiss_desc.csv"), row.names = FALSE)

# Only respondents whose SVP choice differs across ballots identify the
# within-respondent comparison; this count, not the sample size, is the
# effective basis of the H2 test.
write.csv(data.frame(
  discordant  = sum((che$vote_nc == SVP_CODE) != (che$vote_cs == SVP_CODE)),
  svp_both    = sum(che$vote_nc == SVP_CODE & che$vote_cs == SVP_CODE),
  svp_nc_only = sum(che$vote_nc == SVP_CODE & che$vote_cs != SVP_CODE),
  svp_cs_only = sum(che$vote_nc != SVP_CODE & che$vote_cs == SVP_CODE)),
  file.path(OUT, "table_discordant.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 12. STUDY 2: MULTILEVEL LOGIT AND CONDITIONAL LOGIT
# -----------------------------------------------------------------------------
# Cross-classified random intercepts: (1|respondent_id) absorbs everything
# fixed about a person, identifying the ballot coefficient within persons;
# (1|canton) absorbs local SVP strength. anti_elite:pr is the H2 quantity.
m2 <- cache("m2", glmer(
  svp_vote ~ anti_elite * pr + age_c + female + edu_c + lr_c + satdem_c +
    (1 | respondent_id) + (1 | canton),
  data = long, family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 200000))))
print(summary(m2))
write.csv(fmt_model(m2), file.path(OUT, "table_model2.csv"), row.names = FALSE)

m2_coef <- coef_matrix(m2)
m2_sds  <- as.data.frame(VarCorr(m2))$sdcor
write.csv(data.frame(
  n_obs      = nobs(m2),
  n_resp     = ngrps(m2)[["respondent_id"]],
  n_canton   = ngrps(m2)[["canton"]],
  aic        = round(AIC(m2), 1),
  var_resp   = round(as.data.frame(VarCorr(m2))$vcov[1], 3),
  var_canton = round(as.data.frame(VarCorr(m2))$vcov[2], 3)),
  file.path(OUT, "table_model2_fit.csv"), row.names = FALSE)

# Conditional logit stratified by respondent: no distributional assumption
# on the person effects; everything constant within a person drops out of
# the likelihood, leaving only pr and its interaction with attitudes.
m2_clogit <- clogit(svp_vote ~ pr + anti_elite:pr + strata(respondent_id),
                    data = long)
print(summary(m2_clogit))
cl <- summary(m2_clogit)$coefficients
write.csv(data.frame(term = rownames(cl),
                     b  = round(cl[, "coef"], 3),
                     se = round(cl[, "se(coef)"], 3),
                     p  = round(cl[, "Pr(>|z|)"], 4)),
          file.path(OUT, "table_model2_clogit.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 13. STUDY 2: PREDICTED PROBABILITIES AND FIGURE 2
# -----------------------------------------------------------------------------
# Study 2 random effects: respondent and canton deviations, independent and
# constant across the grid.
re_draws_study2 <- function(nd, n) {
  u <- rnorm(n, 0, m2_sds[1]) + rnorm(n, 0, m2_sds[2])
  matrix(rep(u, each = nrow(nd)), nrow = nrow(nd))
}

grid2 <- expand.grid(anti_elite = seq(-2, 2, by = 0.1), pr = c(0, 1))
grid2$age_c    <- 0
grid2$female   <- mean(long$female)
grid2$edu_c    <- 0
grid2$lr_c     <- 0
grid2$satdem_c <- 0

pred2 <- predict_marginal(m2, grid2,
  ~ anti_elite * pr + age_c + female + edu_c + lr_c + satdem_c,
  re_draws_study2)
pred2$ballot <- factor(pred2$pr, levels = c(1, 0),
  labels = c("Proportional ballot",
             "Majoritarian ballot"))

fig2 <- ggplot(pred2, aes(anti_elite, fit, colour = ballot, fill = ballot)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("#1B6CA8", "#B03A2E"), name = NULL) +
  scale_fill_manual(values = c("#1B6CA8", "#B03A2E"), name = NULL) +
  labs(x = "Anti-elite attitudes (standard deviations from the Swiss mean)",
       y = "Predicted probability of an SVP vote") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig2_switzerland.png"), fig2,
       width = 6.5, height = 4.2, dpi = 300)
write.csv(pred2[, c("anti_elite", "pr", "fit", "lo", "hi")],
          file.path(OUT, "pred_study2.csv"), row.names = FALSE)
fig2

me2 <- pred2 %>% filter(anti_elite %in% c(-0.5, 0.5)) %>% group_by(pr) %>%
  summarise(effect_pp = round(100 * (fit[anti_elite == 0.5] -
                                     fit[anti_elite == -0.5]), 1),
            .groups = "drop")
write.csv(me2, file.path(OUT, "table_marginal_effects2.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 14. ROBUSTNESS TABLE AND SESSION INFORMATION
# -----------------------------------------------------------------------------
h1_main <- m1_coef["anti_elite:logmag_c", ]
h2_main <- m2_coef["anti_elite:pr", ]
robustness <- rbind(
  interaction_row("H1, main specification (study median magnitude)",
                  h1_main[1], h1_main[2], h1_main[4], nobs(m1)),
  row_district, row_t6, row_t8,
  interaction_row("H2, main specification (multilevel logit)",
                  h2_main[1], h2_main[2], h2_main[4], nrow(long)),
  interaction_row("H2, conditional logit stratified by respondent",
                  cl["pr:anti_elite", "coef"], cl["pr:anti_elite", "se(coef)"],
                  cl["pr:anti_elite", "Pr(>|z|)"], nrow(long)))
write.csv(robustness, file.path(OUT, "table_robustness.csv"), row.names = FALSE)
print(robustness)

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("Done. Outputs written to", OUT, "\n")
