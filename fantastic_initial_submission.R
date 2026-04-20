rm(list=ls())
library(data.table)
library(survival)
library(survminer)
library(openxlsx)
library(car)
library(gtsummary)
library(officer)
library(flextable)
library(broom)       
library(dplyr)      
library(tibble)
library(janitor)
library(ggalluvial)
library(psych)
library(ggplot2)
library(reshape2)
library(scales)
library(vautils)
source('~/gitlab/translationaloncology/tmed_stat.R')


#---- primary endpoint ---- 

load("merge_dataset3.Rdat")

table(merge_dataset$Site.of.recurrence,useNA='always')

dim(merge_dataset)
median(merge_dataset$RFS)
table(merge_dataset$Age_group)
prop.table(table(merge_dataset$Age_group))
summary(merge_dataset$Age)
table(merge_dataset$Site_Mets)
prop.table(table(merge_dataset$Site_Mets))*100
table(merge_dataset$ctDNA_MRD_4wk)
table(merge_dataset$ctDNA_postACT_7mo)
table(merge_dataset[!is.na(ctDNA_MRD_4wk),FOLFOXIRI])



dt.pat.table.source = merge_dataset[, .(FOLFOXIRI,ACT, Sex, `Age group`=Age_group, 
                                        `Tumor location`=Tumor_location, `Pathological T stage`=pT.Stage, `Pathological N stage`=pN.Stage, 
                                        `Microsatellite Instability`=MSI, RAS, BRAF, `ECOG PS`=ECOG_PS, Resectability=Margins, 
                                        `Site of Metastases`=Site_Mets,`Synchronicity`=Mets.type)]
dt.pat.table.source[, FOLFOXIRI := factor(FOLFOXIRI, levels = c("Yes", "No"))]

reset_gtsummary_theme()
theme_gtsummary_compact()

table <- dt.pat.table.source |>
  tbl_summary(
    by = FOLFOXIRI
  ) |>
  add_overall() |>  # adds "All Patients"
  # add_n() |>        # adds patient count to headers
  modify_header(    # custom headers with patient count
    stat_0 ~ "**All Patients**<br>N = {n}",
    # stat_1 ~ "**mFOLFOXIRI: No**<br>N = {n}",
    # stat_2 ~ "**mFOLFOXIRI: Yes**<br>N = {n}"
    stat_1 ~ "**mFOLFOXIRI Group**<br>N = {n}",
    stat_2 ~ "**SOC Group**<br>N = {n}"
  )

table


table_flextable=as_flex_table(table)
doc <- read_docx() %>%
  body_add_flextable(table_flextable)
print(doc,target='Table1.docx')

dt.pat.table.source = merge_dataset[FOLFOXIRI=='No', .(ACT_regimen, Sex, `Age group`=Age_group, 
                                        `Tumor location`=Tumor_location, `Pathological T stage`=pT.Stage, `Pathological N stage`=pN.Stage, 
                                        `Microsatellite Instability`=MSI, RAS, BRAF, `ECOG PS`=ECOG_PS, Resectability=Margins, 
                                        `Site of Metastases`=Site_Mets,`Synchronicity`=Mets.type)]
reset_gtsummary_theme()
theme_gtsummary_compact()

table <- dt.pat.table.source |>
  tbl_summary(
    by = ACT_regimen
  ) |>
  add_overall()

table


# table_flextable=as_flex_table(table)
# doc <- read_docx() %>%
#   body_add_flextable(table_flextable)
# print(doc,target='TableS1.docx')

## other number for basic stats ##
dim(merge_dataset[ctDNA_MRD_4wk=='POSITIVE' & RFS.event==1])
table(merge_dataset[ctDNA_MRD_4wk=='POSITIVE'  & RFS.event==1,Site.of.recurrence],useNA='always')
prop.table(table(merge_dataset[ctDNA_MRD_4wk=='POSITIVE'  & RFS.event==1,Site.of.recurrence],useNA='always'))
dim(merge_dataset[ctDNA_MRD_4wk=='NEGATIVE' & RFS.event==1])
table(merge_dataset[ctDNA_MRD_4wk=='NEGATIVE'  & RFS.event==1,Site.of.recurrence],useNA='always')
prop.table(table(merge_dataset[ctDNA_MRD_4wk=='NEGATIVE'  & RFS.event==1,Site.of.recurrence],useNA='always'))


table(merge_dataset[ctDNA_postACT_7mo=='POSITIVE'  & RFS.event==1,Site.of.recurrence],useNA='always')
prop.table(table(merge_dataset[ctDNA_postACT_7mo=='POSITIVE'  & RFS.event==1,Site.of.recurrence],useNA='always'))
table(merge_dataset[ctDNA_postACT_7mo=='NEGATIVE'  & RFS.event==1,Site.of.recurrence],useNA='always')
prop.table(table(merge_dataset[ctDNA_postACT_7mo=='NEGATIVE'  & RFS.event==1,Site.of.recurrence],useNA='always'))

merge_dataset[FOLFOXIRI=='Yes' & ctDNA_MRD_4wk=='POSITIVE' & ctDNA_postACT_7mo=='NEGATIVE',]
summary(merge_dataset[FOLFOXIRI=='Yes' & ctDNA_MRD_4wk=='POSITIVE',MTM4w])
merge_dataset[FOLFOXIRI=='Yes' & ctDNA_MRD_4wk=='NEGATIVE' & ctDNA_postACT_7mo=='POSITIVE',]

table(merge_dataset[FOLFOXIRI=='No' & ctDNA_MRD_4wk=='POSITIVE',ACT_regimen])
merge_dataset[FOLFOXIRI=='No' & ctDNA_MRD_4wk=='POSITIVE']
merge_dataset[FOLFOXIRI=='No' & ctDNA_MRD_4wk=='POSITIVE' & ctDNA_postACT_7mo=='NEGATIVE',]
summary(merge_dataset[FOLFOXIRI=='No' & ctDNA_MRD_4wk=='POSITIVE' & ctDNA_postACT_7mo %in% c('POSITIVE','NEGATIVE'),MTM4w])

table(merge_dataset[FOLFOXIRI=='No' & ctDNA_MRD_4wk=='NEGATIVE' ,ctDNA_postACT_7mo])

####

table(merge_dataset[,.(ACT_regimen,FOLFOXIRI)]) # confirm that we can use FOLFOXIRI label
table(merge_dataset[,.(Cycle_mFOLFOXIRI,`Completion.of.mFOLFOXIRI`)]) # comlete mean >=8 cycle

table(merge_dataset[,.(Completion.of.mFOLFOXIRI,FOLFOXIRI)]) # confirm that all completion is only for FALFOXSIRI


completion_logical <- as.logical(merge_dataset$Completion.of.mFOLFOXIRI)
completion_logical <- completion_logical[!is.na(completion_logical)]
num_success <- sum(completion_logical)
num_total <- length(completion_logical)
ci_80 <- binom.test(num_success, num_total, conf.level = 0.80)
ci_95 <- binom.test(num_success, num_total, conf.level = 0.95)
cat("Sample size:", num_total, "\n")
cat("Proportion completing mFOLFOXIRI:", round(ci_95$estimate, 3), "\n")
cat("80% CI: (", round(ci_80$conf.int[1], 3), ",", round(ci_80$conf.int[2], 3), ")\n")
cat("95% CI: (", round(ci_95$conf.int[1], 3), ",", round(ci_95$conf.int[2], 3), ")\n")


n <- num_total
x <- num_success  # observed number of responses = 27
p0 <- 0.35            # null threshold

# One-sided binomial test
bino_p_value=binom.test(x, n, p = p0, alternative = "greater")
bino_p_value$p.value

plot_dt <- data.table(
  CompletionRate = ci_95$estimate,
  CI80_lower = ci_80$conf.int[1],
  CI80_upper = ci_80$conf.int[2],
  CI95_lower = ci_95$conf.int[1],
  CI95_upper = ci_95$conf.int[2],
  p_value=bino_p_value$p.value,
  Label = "Overall",
  N = num_total
)

# Create formatted labels
plot_dt[, CI80_label := sprintf("80%% CI: %.2f–%.2f", CI80_lower, CI80_upper)]
plot_dt[, CI95_label := sprintf("95%% CI: %.2f–%.2f", CI95_lower, CI95_upper)]
plot_dt[, PointLabel := sprintf("%.2f (N = %d)", CompletionRate, N)]
plot_dt[, PValue_label := sprintf("p = %.2g", p_value)]

# Plot
ggplot(plot_dt, aes(x = Label, y = CompletionRate)) +
  geom_errorbar(aes(ymin = CI95_lower, ymax = CI95_upper), width = 0.1, size = 0.6, color = "black") +
  geom_errorbar(aes(ymin = CI80_lower, ymax = CI80_upper), width = 0.1, size = 2, color = "darkgray") +
  geom_point(size = 4, color = "black") +
  
  # Add point estimate + N label beside the dot
  geom_text(
    aes(label = PointLabel),
    position = position_nudge(x = 0.26),
    vjust = 0.5,
    size = 5
  ) +
  geom_text(
    aes(y = CompletionRate - 0.08, label = PValue_label),
    position = position_nudge(x = 0.24),
    vjust = 0.5,
    size = 5,
    color = "black"
  ) +
    geom_text(
    aes(y = CI80_upper - 0.01, label = CI80_label),
    position = position_nudge(x = 0.3),
    size = 5, color = "darkgray"
  ) +
  geom_text(
    aes(y = CI95_upper + 0.01, label = CI95_label),
    position = position_nudge(x = 0.3),
    size = 5, color = "black"
  ) +

  
  labs(
    title = "Completion of mFOLFOXIRI with 80% and 95% Confidence Intervals",
    x = NULL,
    y = "Proportion Completed"
  ) +
  ylim(0, 1) +
  theme_minimal()




dt <- merge_dataset[!is.na(Completion.of.mFOLFOXIRI) & !is.na(Site_Mets)]
dt[, Completion := as.logical(Completion.of.mFOLFOXIRI)]

site_dt <- dt[, {
  successes <- sum(Completion)
  total <- .N
  ci <- binom.test(successes, total, conf.level = 0.95)
  list(
    successes=successes,
    total=total,
    CompletionRate = ci$estimate,
    CI_lower = ci$conf.int[1],
    CI_upper = ci$conf.int[2],
    N = total
  )
}, by = Site_Mets]

site_dt[, Label := sprintf("%.2f \n(N = %d)", CompletionRate, N)]

setorder(site_dt, -CompletionRate)
site_dt[, Site_Mets := factor(Site_Mets, levels = Site_Mets)]

# Plot
ggplot(site_dt, aes(x = Site_Mets, y = CompletionRate)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, size = 0.6, color = "black") +
  geom_point(size = 4, color = "black") +
  geom_text(
    aes(label = Label),
    position = position_nudge(x = 0.3),
    vjust = 0.5,
    size = 4.2
  ) +
  labs(
    title = "Completion of mFOLFOXIRI by Metastasis sites with 95% Confidence Intervals",
    x = "Site of Metastasis",
    y = "Proportion Completed"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


dt <- merge_dataset[!is.na(Completion.of.mFOLFOXIRI) & !is.na(Site_Mets)]
dt[, Complete_mFOLFOXIRI := fifelse(Completion.of.mFOLFOXIRI == "TRUE", 1L,
                                    fifelse(Completion.of.mFOLFOXIRI == "FALSE", 0L, NA_integer_))]
dt=dt[,.(Complete_mFOLFOXIRI, Age_group, Sex, ECOG_PS,Liver_met,Site_Mets,Tumor_location)]
dt
dt[, Sex := factor(Sex, levels = c("Female","Male"))]  
dt[, Age_group := factor(Age_group, levels = c("<65",">=65"))]
dt[, ECOG_PS := factor(ECOG_PS,levels=c(0,1))]
dt[, Liver_met := factor(Liver_met, levels = c("others", "liver_met"))]  
dt[,Site_Mets:=factor(Site_Mets)]
dt

#####################################
##---- landmark consideration ---- ##
#####################################
merge_dataset %>% names()
# merge_dataset[,.((Date_recurrence-Date_surgery)/30,RFS)]
avg_month <- 30
merge_dataset[, `:=`(
  p_mrel_ctDNA4wk  = round(as.numeric(`4wk.DNA.collected` - Date_surgery) / avg_month, 2),
  p_mrel_ctDNA7mo  = round(as.numeric(`7M.ctDNA.collected` - Date_surgery) / avg_month, 2),
  p_mrel_ACT_start = round(as.numeric(Start.of.ACT - Date_surgery) / avg_month, 2),
  p_mrel_ACT_end   = round(as.numeric(Last.date.ACT.treatment - Date_surgery) / avg_month, 2)
)]

cat("--- Summary: Relative Months from Surgery ---\n")
merge_dataset[, lapply(.SD, function(x) {
  c(max = round(max(x, na.rm = TRUE), 2), 
    Median = round(median(x, na.rm = TRUE), 2))
}), .SDcols = c("p_mrel_ctDNA4wk", "p_mrel_ctDNA7mo", "p_mrel_ACT_start", "p_mrel_ACT_end")]


# cat("--- Summary: Relative days from Surgery ---\n")
# merge_dataset[, lapply(.SD, function(x) {
#   c(max = round(max(x*30, na.rm = TRUE), 2), 
#     Median = round(median(x*30, na.rm = TRUE), 2))
# }), .SDcols = c("p_mrel_ctDNA4wk", "p_mrel_ctDNA7mo", "p_mrel_ACT_start", "p_mrel_ACT_end")]


landmarkctDNA4wk=1
landmarkctDNA7mo=7
landmark_medianACTend=7.5



#################################
##---- secondary endpoint ---- ##
#################################
#---- secondary endpoint: positivity ---- 


dt <- merge_dataset[!is.na(ctDNA_MRD_4wk) & !is.na(FOLFOXIRI)]
dt[, ctDNA_MRD_4wk := factor(ctDNA_MRD_4wk, levels = c("NEGATIVE", "POSITIVE"))]
dt[, FOLFOXIRI := factor(FOLFOXIRI, levels = c("No", "Yes"))]

tab <- table(dt[, .(ctDNA_MRD_4wk, FOLFOXIRI)])
prop_tab <- prop.table(tab, margin = 2)
col_totals <- colSums(tab)
total_n <- sum(tab)

summary_table <- data.table(
  `MRD Status` = c("POSITIVE", "NEGATIVE")
)

summary_table[, `FOLFOXIRI = Yes` := sprintf("%d (%.1f%%)", 
                                             tab["POSITIVE", "Yes"], 100 * prop_tab["POSITIVE", "Yes"])]
summary_table[2, `FOLFOXIRI = Yes` := sprintf("%d (%.1f%%)", 
                                              tab["NEGATIVE", "Yes"], 100 * prop_tab["NEGATIVE", "Yes"])]

summary_table[, `FOLFOXIRI = No` := sprintf("%d (%.1f%%)", 
                                            tab["POSITIVE", "No"], 100 * prop_tab["POSITIVE", "No"])]
summary_table[2, `FOLFOXIRI = No` := sprintf("%d (%.1f%%)", 
                                             tab["NEGATIVE", "No"], 100 * prop_tab["NEGATIVE", "No"])]

summary_table[, `Total` := sprintf("%d (%.1f%%)", 
                                   sum(tab["POSITIVE", ]), 100 * sum(tab["POSITIVE", ]) / total_n)]
summary_table[2, Total := sprintf("%d (%.1f%%)", 
                                  sum(tab["NEGATIVE", ]), 100 * sum(tab["NEGATIVE", ]) / total_n)]

summary_table <- rbind(summary_table, data.table(
  `MRD Status` = "Total",
  `FOLFOXIRI = Yes` = as.character(col_totals["Yes"]),
  `FOLFOXIRI = No` = as.character(col_totals["No"]),
  `Total` = as.character(total_n)
))

p_val <- fisher.test(tab)$p.value
summary_table <- rbind(summary_table, data.table(
  `MRD Status` = "P-value (Fisher's)",
  `FOLFOXIRI = Yes` = sprintf("P = %.3g", p_val),
  `FOLFOXIRI = No` = "",
  `Total` = ""
))
summary_table



dt <- merge_dataset[!is.na(ctDNA_postACT_7mo) & ctDNA_postACT_7mo %in% c('POSITIVE','NEGATIVE') & !is.na(FOLFOXIRI)]
dt[, ctDNA_postACT_7mo := factor(ctDNA_postACT_7mo, levels = c("NEGATIVE", "POSITIVE"))]
dt[, FOLFOXIRI := factor(FOLFOXIRI, levels = c("No", "Yes"))]
prop.table(table(dt[,.(ctDNA_postACT_7mo)]))


tab <- table(dt[, .(ctDNA_postACT_7mo, FOLFOXIRI)])

count_tab <- addmargins(tab)
prop_tab <- prop.table(tab, margin = 2)
fisher_test <- fisher.test(tab)

cat("\n=== MRD Status by FOLFOXIRI Treatment ===\n")
print(count_tab)
cat("\nProportions by Treatment Arm:\n")
print(round(100 * prop_tab, 1))  # Show as percentages
cat("\nFisher's Exact Test p-value:", signif(fisher_test$p.value, 3), "\n")

###########################################
##---- secondary endpoint: RFS & OS ---- ##
###########################################
#### ctDNA positive/negative
## RFS 4wk

landmark_month <-landmarkctDNA4wk
# landmark_month <- 12*14/30

# df <- as.data.table(merge_dataset)
# df <- as.data.table(merge_dataset[FOLFOXIRI=='Yes'])
df <- as.data.table(merge_dataset[FOLFOXIRI!='Yes'])
dim(df)


df <- df[!is.na(ctDNA_MRD_4wk) & !is.na(RFS) & !is.na(RFS.event)]
df[, ctDNA_MRD_4wk := factor(ctDNA_MRD_4wk, levels = c("NEGATIVE", "POSITIVE"))]
df <- df[RFS > landmark_month]
df[, RFS := RFS - landmark_month]
dim(df)
cox_fit <- coxph(Surv(RFS, RFS.event) ~ ctDNA_MRD_4wk, data = df)
cox_summary <- summary(cox_fit)
hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
lower_ci <- round(cox_summary$conf.int[,"lower .95"], 2)
upper_ci <- round(cox_summary$conf.int[,"upper .95"], 2)
hr_label <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", hr, lower_ci, upper_ci)

km_fit <- survfit(Surv(RFS, RFS.event) ~ ctDNA_MRD_4wk, data = df)

medians <- summary(km_fit)$table
median_dt <- data.table(
  group = rownames(medians),
  median = medians[,"median"]
)
median_dt[, group := gsub("^.*=", "", group)]
median_dt[, median_text := ifelse(is.na(median), "Median: NR", paste0("Median: ", round(median, 1), " mo"))]

summary_fit <- summary(km_fit, times = c(12,24, 36))
surv_df <- data.table(
  group = gsub("^.*=", "", summary_fit$strata),
  time = summary_fit$time,
  surv = summary_fit$surv
)

surv_text <- surv_df[!is.na(surv), .(
  text = paste0(time / 12, "y: ", scales::percent(surv, digits = 1))
), by = group][
  , .(rfs = paste(text, collapse = "; ")), by = group
]

surv_text <- merge(surv_text, median_dt[, .(group, median_text)], by = "group", all.x = TRUE)
surv_text[, label := paste0(group, ": ", median_text, "; RFS: ", rfs)]
surv_line <- paste(surv_text$label, collapse = "\n")
hr_label <- paste0(hr_label, "\n", surv_line)

g=ggsurvplot(
  km_fit,
  data = df,
  conf.int = FALSE,
  risk.table = TRUE,
  break.time.by = 6,
  lwd = 2,
  font.legend = c(14, "bold"),      
  pval = FALSE,
  xlab = "Time (months from landmark)",
  ylab = "Recurrence-Free Survival Probability",
  # title = "RFS by MRD Status at 4 Weeks Post-surgery",
  # subtitle = hr_label,
  legend.title = "ctDNA MRD (4 wk)",
  legend.labs = c("NEGATIVE", "POSITIVE"),
  palette = c("#46AC46", "#004677"),
  surv.median.line = "hv",
  censor = TRUE,
  ggtheme = theme_bw() + theme(plot.subtitle = element_text(size = 16, face="bold")),
  tables.theme = theme(
    axis.title.y = element_text(size = 8)
  )
)

g
# pdf('Figure2A_DFS_ctDNA4wk_allpop.pdf')
# pdf('Figure2B_DFS_ctDNA4wk_mFOLFOXIRI.pdf')
pdf('Figure2C_DFS_ctDNA4wk_SOC.pdf')

print(g, newpage = FALSE)
dev.off()
report_surv_extra(km_fit, cox_fit)

## RFS 7mo

landmark_month <- landmarkctDNA7mo

# df <- as.data.table(merge_dataset)
# df <- as.data.table(merge_dataset[FOLFOXIRI=='Yes'])
df <- as.data.table(merge_dataset[FOLFOXIRI!='Yes'])

df <- df[!is.na(ctDNA_postACT_7mo) & !is.na(RFS) & !is.na(RFS.event)]
df[, ctDNA_postACT_7mo := factor(ctDNA_postACT_7mo, levels = c("NEGATIVE", "POSITIVE"))]
df <- df[RFS > landmark_month]
df[, RFS := RFS - landmark_month]

cox_fit <- coxph(Surv(RFS, RFS.event) ~ ctDNA_postACT_7mo, data = df)
cox_summary <- summary(cox_fit)
hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
lower_ci <- round(cox_summary$conf.int[,"lower .95"], 2)
upper_ci <- round(cox_summary$conf.int[,"upper .95"], 2)
hr_label <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", hr, lower_ci, upper_ci)

km_fit <- survfit(Surv(RFS, RFS.event) ~ ctDNA_postACT_7mo, data = df)

medians <- summary(km_fit)$table
median_dt <- data.table(
  group = rownames(medians),
  median = medians[,"median"]
)
median_dt[, group := gsub("^.*=", "", group)]
median_dt[, median_text := ifelse(is.na(median), "Median: NR", paste0("Median: ", round(median, 1), " mo"))]

summary_fit <- summary(km_fit, times = c(12,24, 36))
surv_df <- data.table(
  group = gsub("^.*=", "", summary_fit$strata),
  time = summary_fit$time,
  surv = summary_fit$surv
)

surv_text <- surv_df[!is.na(surv), .(
  text = paste0(time / 12, "y: ", scales::percent(surv, digits = 1))
), by = group][
  , .(rfs = paste(text, collapse = "; ")), by = group
]

surv_text <- merge(surv_text, median_dt[, .(group, median_text)], by = "group", all.x = TRUE)
surv_text[, label := paste0(group, ": ", median_text, "; RFS: ", rfs)]
surv_line <- paste(surv_text$label, collapse = "\n")
hr_label <- paste0(hr_label, "\n", surv_line)

km_fit <- survfit(Surv(RFS, RFS.event) ~ ctDNA_postACT_7mo, data = df)
g=ggsurvplot(
  km_fit,
  data = df,
  lwd = 2,
  font.legend = c(14, "bold"), 
  conf.int = FALSE,
  risk.table = TRUE,
  break.time.by = 6,
  pval = FALSE,
  xlab = "Time (months from landmark)",
  ylab = "Recurrence-Free Survival Probability",
  # title = sprintf("RFS by post-ACT Status at 7 Months Post surgery"),
  # subtitle = hr_label,
  # title = sprintf("RFS by post-ACT Status at 7 Months Post surgery (Landmark: %g months)", landmark_month),
  legend.title = "ctDNA post-ACT (7 mo)",
  legend.labs = c("NEGATIVE", "POSITIVE"),
  palette = c("#46AC46", "#004677"),
  surv.median.line = "hv",
  censor = TRUE,
  ggtheme = theme_bw() + theme(plot.subtitle = element_text(size = 16, face="bold")),
  tables.theme = theme(
    axis.title.y = element_text(size = 8)
  )
)

g
# pdf('Figure2D_DFS_ctDNA7mo_allpop.pdf')
# pdf('Figure2E_DFS_ctDNA7mo_mFOLFOXIRI.pdf')
pdf('Figure2F_DFS_ctDNA7mo_SOC.pdf')

print(g, newpage = FALSE)
dev.off()
report_surv_extra(km_fit, cox_fit)


# ## OS
# 
# landmark_month <- landmarkctDNA4wk
# merge_dataset[,OS:=Date_LAST_fu-Date_surgery,]
# 
# df <- as.data.table(merge_dataset)
# # df <- as.data.table(merge_dataset[FOLFOXIRI=='Yes'])
# # df <- as.data.table(merge_dataset[FOLFOXIRI!='Yes'])
# dim(df)
# 
# df <- df[!is.na(ctDNA_MRD_4wk) & !is.na(OS) & !is.na(OS.event)]
# df[, ctDNA_MRD_4wk := factor(ctDNA_MRD_4wk, levels = c("NEGATIVE", "POSITIVE"))]
# df <- df[OS > landmark_month]
# df[, OS := OS - landmark_month]
# 
# cox_fit <- coxph(Surv(OS, OS.event) ~ ctDNA_MRD_4wk, data = df)
# cox_summary <- summary(cox_fit)
# hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
# lower_ci <- round(cox_summary$conf.int[,"lower .95"], 2)
# upper_ci <- round(cox_summary$conf.int[,"upper .95"], 2)
# hr_label <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", hr, lower_ci, upper_ci)
# 
# km_fit <- survfit(Surv(OS, OS.event) ~ ctDNA_MRD_4wk, data = df)
# group_names <- names(km_fit$strata)
# group_labels <- gsub("^.*=", "", group_names)
# 
# summary_fit <- summary(km_fit, times = c(12,24, 36))
# 
# surv_df <- data.table(
#   group = gsub("^.*=", "", summary_fit$strata),
#   time = summary_fit$time,
#   surv = summary_fit$surv
# )
# 
# surv_text <- surv_df[!is.na(surv), .(
#   text = paste0(time / 12, "y: ", scales::percent(surv, accuracy = 0.1))
# ), by = group][
#   , .(text = paste(text, collapse = "; ")), by = group
# ][
#   , paste0(group, " OS: ", text)
# ]
# 
# surv_line <- paste(surv_text, collapse = "\n")
# hr_label <- paste0(hr_label, "\n", surv_line)
# 
# 
# km_fit <- survfit(Surv(OS, OS.event) ~ ctDNA_MRD_4wk, data = df)
# ggsurvplot(
#   km_fit,
#   data = df,
#   conf.int = FALSE,
#   risk.table = TRUE,
#   break.time.by = 6,
#   pval = TRUE,
#   xlab = "Time (months from landmark)",
#   ylab = "Recurrence-Free Survival Probability",
#   title = sprintf("OS by MRD Status at 4 Weeks Post-surgery"),
#   subtitle = hr_label,
#   # title = sprintf("OS by MRD Status at 4 Weeks Post-surgery (Landmark: %g months)", landmark_month),
#   legend.title = "ctDNA MRD (4 wk)",
#   legend.labs = c("NEGATIVE", "POSITIVE"),
#   palette = c("#46AC46", "#004677"),
#   surv.median.line = "hv",
#   censor = TRUE,
#   ggtheme = theme_bw() + theme(plot.subtitle = element_text(size = 12))  
# )
# 
# 
# landmark_month <- landmarkctDNA7mo
# df <- as.data.table(merge_dataset)
# # df <- as.data.table(merge_dataset[FOLFOXIRI=='Yes'])
# # df <- as.data.table(merge_dataset[FOLFOXIRI!='Yes'])
# 
# df <- df[!is.na(ctDNA_postACT_7mo) & !is.na(OS) & !is.na(OS.event) & ctDNA_postACT_7mo %in% c('POSITIVE','NEGATIVE')]
# df[, ctDNA_postACT_7mo := factor(ctDNA_postACT_7mo, levels = c("NEGATIVE", "POSITIVE"))]
# df <- df[OS > landmark_month]
# df[, OS := OS - landmark_month]
# 
# cox_fit <- coxph(Surv(OS, OS.event) ~ ctDNA_postACT_7mo, data = df)
# cox_summary <- summary(cox_fit)
# hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
# lower_ci <- round(cox_summary$conf.int[,"lower .95"], 2)
# upper_ci <- round(cox_summary$conf.int[,"upper .95"], 2)
# hr_label <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", hr, lower_ci, upper_ci)
# 
# km_fit <- survfit(Surv(OS, OS.event) ~ ctDNA_postACT_7mo, data = df)
# group_names <- names(km_fit$strata)
# group_labels <- gsub("^.*=", "", group_names)
# 
# summary_fit <- summary(km_fit, times = c(12,24, 36))
# 
# surv_df <- data.table(
#   group = gsub("^.*=", "", summary_fit$strata),
#   time = summary_fit$time,
#   surv = summary_fit$surv
# )
# 
# surv_text <- surv_df[!is.na(surv), .(
#   text = paste0(time / 12, "y: ", scales::percent(surv, accuracy = 0.1))
# ), by = group][
#   , .(text = paste(text, collapse = "; ")), by = group
# ][
#   , paste0(group, " OS: ", text)
# ]
# 
# surv_line <- paste(surv_text, collapse = "\n")
# hr_label <- paste0(hr_label, "\n", surv_line)
# 
# km_fit <- survfit(Surv(OS, OS.event) ~ ctDNA_postACT_7mo, data = df)
# ggsurvplot(
#   km_fit,
#   data = df,
#   conf.int = FALSE,
#   risk.table = TRUE,
#   break.time.by = 6,
#   pval = TRUE,
#   xlab = "Time (months from landmark)",
#   ylab = "Recurrence-Free Survival Probability",
#   title = sprintf("OS by post-ACT Status at 7 Months Post surgery"),
#   subtitle = hr_label,
#   # title = sprintf("OS by post-ACT Status at 7 Months Post surgery (Landmark: %g months)", landmark_month),
#   legend.title = "ctDNA post-ACT (7 mo)",
#   legend.labs = c("NEGATIVE", "POSITIVE"),
#   palette = c("#46AC46", "#004677"),
#   surv.median.line = "hv",
#   censor = TRUE,
#   ggtheme = theme_bw() + theme(plot.subtitle = element_text(size = 12))  
# )

#### observation vs treatment

### RFS
# landmark_month <- 535/30
landmark_month <- landmark_medianACTend

df <- as.data.table(merge_dataset)
# df <- as.data.table(merge_dataset[ctDNA_MRD_4wk=='POSITIVE'])
# df <- as.data.table(merge_dataset[ctDNA_MRD_4wk=='NEGATIVE'])
# df <- as.data.table(merge_dataset[ctDNA_postACT_7mo=='POSITIVE'])
# df <- as.data.table(merge_dataset[ctDNA_postACT_7mo=='NEGATIVE'])
dim(df)

df <- df[!is.na(Group) & !is.na(RFS) & !is.na(RFS.event)]
df[, Group := factor(Group,levels=c("SOC","mFOLFOXIRI"))]
df <- df[RFS > landmark_month]
df[, RFS := RFS - landmark_month]

cox_fit <- coxph(Surv(RFS, RFS.event) ~ Group, data = df)
cox_summary <- summary(cox_fit)
hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
lower_ci <- round(cox_summary$conf.int[,"lower .95"], 2)
upper_ci <- round(cox_summary$conf.int[,"upper .95"], 2)
hr_label <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", hr, lower_ci, upper_ci)

km_fit <- survfit(Surv(RFS, RFS.event) ~ Group, data = df)

medians <- summary(km_fit)$table
median_dt <- data.table(
  group = rownames(medians),
  median = medians[,"median"]
)
median_dt[, group := gsub("^.*=", "", group)]
median_dt[, median_text := ifelse(is.na(median), "Median: NR", paste0("Median: ", round(median, 1), " mo"))]

summary_fit <- summary(km_fit, times = c(12,24, 36))
surv_df <- data.table(
  group = gsub("^.*=", "", summary_fit$strata),
  time = summary_fit$time,
  surv = summary_fit$surv
)

surv_text <- surv_df[!is.na(surv), .(
  text = paste0(time / 12, "y: ", scales::percent(surv, digits = 1))
), by = group][
  , .(rfs = paste(text, collapse = "; ")), by = group
]

surv_text <- merge(surv_text, median_dt[, .(group, median_text)], by = "group", all.x = TRUE)
surv_text[, label := paste0(group, ": ", median_text, "; RFS: ", rfs)]
surv_line <- paste(surv_text$label, collapse = "\n")
hr_label <- paste0(hr_label, "\n", surv_line)


km_fit <- survfit(Surv(RFS, RFS.event) ~ Group, data = df)

g=ggsurvplot(
  km_fit,
  data = df,
  conf.int = FALSE,
  risk.table = TRUE,
  break.time.by = 6,
  lwd = 2,
  font.legend = c(14, "bold"),
  pval = FALSE,
  xlab = "Time (months from landmark)",
  ylab = "Recurrence-Free Survival Probability",
  # title = sprintf("RFS by Treatment Group"),
  # subtitle = hr_label,
  # title = sprintf("RFS by Treatment Group  (Landmark: %g months)", landmark_month),
  legend.title = "Group",
  legend.labs = c("SOC","mFOLFOXIRI"),
  palette = c("orange","purple"),
  surv.median.line = "hv",
  censor = TRUE,
  ggtheme = theme_bw() + theme(plot.subtitle = element_text(size = 16, face="bold")),
  tables.theme = theme(
    axis.title.y = element_text(size = 8)
  )
  )


g
pdf('Figure3A_DFS_all_pop.pdf')
# pdf('Figure3B_DFS_4wk_ctDNApositive.pdf')
# pdf('Figure3D_DFS_4wk_ctDNAnegative.pdf')
# pdf('Figure3C_DFS_7mo_ctDNApositive.pdf')
# pdf('Figure3E_DFS_7mo_ctDNAnegative.pdf')

print(g, newpage = FALSE)
dev.off()
report_surv_extra(km_fit, cox_fit)

# ### OS
# 
# landmark_month <- landmark_medianACTend
# 
# df <- as.data.table(merge_dataset)
# # df <- as.data.table(merge_dataset[ctDNA_MRD_4wk=='POSITIVE'])
# # df <- as.data.table(merge_dataset[ctDNA_MRD_4wk=='NEGATIVE'])
# # df <- as.data.table(merge_dataset[ctDNA_postACT_7mo=='POSITIVE'])
# # df <- as.data.table(merge_dataset[ctDNA_postACT_7mo=='NEGATIVE'])
# dim(df)
# 
# df <- df[!is.na(Group) & !is.na(OS) & !is.na(OS.event)]
# df[, Group := factor(Group,levels=c("SOC","mFOLFOXIRI"))]
# 
# df <- df[OS > landmark_month]
# df[, OS := OS - landmark_month]
# 
# cox_fit <- coxph(Surv(OS, OS.event) ~ Group, data = df)
# cox_summary <- summary(cox_fit)
# hr <- round(cox_summary$coefficients[,"exp(coef)"], 2)
# lower_ci <- round(cox_summary$conf.int[,"lower .95"], 2)
# upper_ci <- round(cox_summary$conf.int[,"upper .95"], 2)
# hr_label <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", hr, lower_ci, upper_ci)
# 
# km_fit <- survfit(Surv(OS, OS.event) ~ Group, data = df)
# group_names <- names(km_fit$strata)
# group_labels <- gsub("^.*=", "", group_names)
# 
# summary_fit <- summary(km_fit, times = c(12,24, 36))
# 
# surv_df <- data.table(
#   group = gsub("^.*=", "", summary_fit$strata),
#   time = summary_fit$time,
#   surv = summary_fit$surv
# )
# 
# surv_text <- surv_df[!is.na(surv), .(
#   text = paste0(time / 12, "y: ", scales::percent(surv, accuracy = 0.1))
# ), by = group][
#   , .(text = paste(text, collapse = "; ")), by = group
# ][
#   , paste0(group, " OS: ", text)
# ]
# 
# surv_line <- paste(surv_text, collapse = "\n")
# hr_label <- paste0(hr_label, "\n", surv_line)
# 
# km_fit <- survfit(Surv(OS, OS.event) ~ Group, data = df)
# ggsurvplot(
#   km_fit,
#   data = df,
#   conf.int = FALSE,
#   risk.table = TRUE,
#   break.time.by = 6,
#   pval = TRUE,
#   xlab = "Time (months from landmark)",
#   ylab = "Recurrence-Free Survival Probability",
#   title = sprintf("OS by Treatment Group"),
#   subtitle = hr_label,
#   # title = sprintf("OS by Treatment Group  (Landmark: %g months)", landmark_month),
#   legend.title = "Group",
#   legend.labs = c("SOC","mFOLFOXIRI"),
#   palette = c("orange","purple"),
#   surv.median.line = "hv",
#   censor = TRUE,
#   ggtheme = theme_bw() + theme(plot.subtitle = element_text(size = 12))  
# )


#---- clearance ---- 


# dt <- merge_dataset[FOLFOXIRI == "Yes"]
dt <- merge_dataset[FOLFOXIRI == "No"]

# Count transitions
transition_counts <- dt[, .N, by = .(ctDNA_MRD_4wk, ctDNA_postACT_7mo)]
transition_counts[is.na(ctDNA_MRD_4wk),ctDNA_MRD_4wk:='Not tested']

transition_counts[is.na(ctDNA_postACT_7mo),ctDNA_postACT_7mo:='Not tested']
transition_counts[ctDNA_postACT_7mo=='RECURRENCE',ctDNA_postACT_7mo:='Not tested']

# Rename columns to match ggalluvial expectations
setnames(transition_counts,
         c("ctDNA_MRD_4wk", "ctDNA_postACT_7mo"),
         c("MRD_4wk", "MRD_end"))

## Add midpoint positions between axis1 and axis2 for label placement
transition_counts[, label_pos := 1]  # dummy placeholder, not used

ggplot(transition_counts,
       aes(axis1 = MRD_4wk, axis2 = MRD_end, y = N)) +
  
  # Flows colored by destination
  geom_alluvium(aes(fill = MRD_end), width = 1/12, alpha = 0.7) +
  
  # Strata bars colored by MRD state
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  
  # Flow count labels
  geom_text(stat = "alluvium",
            aes(label = N),
            position = "identity",
            size = 3,
            color = "black",
            vjust = 0.5,
            fontface = "bold") +
  
  scale_x_discrete(limits = c("ctDNA at MRD", "ctDNA at end of ACT"),
                   expand = c(.1, .1)) +
  
  scale_fill_manual(values = c(
    "NEGATIVE" = "lightgreen",
    "POSITIVE" = "coral",
    "Not tested" = "grey"
  )) +
  
  labs(y = "Number of patients",
       # title = "Clearance Flow: mFOLFOXIRI-treated Patients",
       title = "Clearance Flow: SOC Patients",
       fill = "MRD Status") +
  theme_minimal()

binom.test(1, 16, conf.level = 0.95)
binom.test(3, 9, conf.level = 0.95)

matrix <- matrix(c(3, 6, 1, 15), nrow = 2,
                 dimnames = list(Group = c("Group1", "Group2"),
                                 Outcome = c("Success", "Failure")))

fisher.test(matrix)

#---- Cox ----
landmark_month <- landmark_medianACTend

dt=as.data.table(merge_dataset)
# dt=as.data.table(merge_dataset[Site_Mets=='Liver'])
# dt=as.data.table(merge_dataset[FOLFOXIRI=='Yes'])
dim(dt)


dt[ctDNA_MRD_4wk %in% c('POSITIVE','NEGATIVE','RECURRENCE') & ctDNA_postACT_7mo %in% c('POSITIVE','NEGATIVE','RECURRENCE') ,clearance:=FALSE]

dt[ctDNA_MRD_4wk=='POSITIVE' & ctDNA_postACT_7mo=='NEGATIVE',clearance:=TRUE]
summary(dt$RFS)
dim(dt)
dt <- dt[RFS > landmark_month]
dt[, RFS := RFS - landmark_month]
summary(dt$RFS)
dim(dt)
dt[ctDNA_postACT_7mo=='RECURRENCE',ctDNA_postACT_7mo:= NA_character_]
table(dt$ACT_regimen)
dt[, ACT_regimen := factor(ACT_regimen, levels = c("None", "FU", "CapeOX", "mFOLFOXIRI"))]
dt[, FOLFOXIRI := factor(FOLFOXIRI)]

dt[, `Metastasis type` := factor(`Mets.type`,levels=c('Synchronous','Metachronous'))]
dt[, `Metastasis size (mm)` := factor(`Metastasis size (mm)`,levels=c('<50','>=50'))]
dt[, `Number of metastasis` := factor(`Number of metastasis`,levels=c('1','>1'))]
dt[, `RAS` := factor(`RAS`,levels=c('WT','MUT'))]
dt[, `pT.Stage` := factor(`pT.Stage`,levels=c('T1','T2','T3','T4'))]
dt[, `Margins` := factor(`Margins`,levels=c('R0','R1','R2'))]


cox_model <- coxph(Surv(RFS, RFS.event) ~ ctDNA_MRD_4wk  + Sex  + FOLFOXIRI+ `Metastasis type`+RAS , data = dt,x=TRUE)



summary(cox_model)
g=forest_plot_cox(cox_model)
g
pdf('Figure4.pdf', width = 12, height = 6)
print(g, newpage = FALSE)
dev.off()




#---- table 2 ----
data=flexread('Clinical Data FANTASTIC_20250502updated.xlsx')
table(data[FOLFOXIRI==1,Specify])# grade 2
table(data[FOLFOXIRI==1,specify]) #grade 3
table(data[FOLFOXIRI==1,.(Specify,specify)]) #both grade

dim(data[FOLFOXIRI==1])


################
## propensity ##
################


# Load necessary libraries for IPTW and weighting
if(!require(WeightIt)) install.packages("WeightIt")
if(!require(cobalt)) install.packages("cobalt")
library(WeightIt)
library(cobalt)
library(survival)
library(survminer)

# 1. Prepare the Landmark Dataset (Survival at 7.5 months)
ps_dt <- as.data.table(merge_dataset)
ps_dt <- ps_dt[!is.na(Group) & !is.na(RFS) & !is.na(RFS.event)]
ps_dt <- ps_dt[RFS > landmark_medianACTend] # Entry into Landmark cohort [cite: 53-55]

# 2. Convert Treatment to Binary (0/1) 
ps_dt[, treat := ifelse(Group == "mFOLFOXIRI", 1, 0)]

# 3. Estimate Weights using IPTW
# NOTE: MSI is excluded because it only has 1 level in the landmark cohort
# 'ate' calculates weights for the Average Treatment Effect
weight_obj <- weightit(treat ~ Age + Sex + Tumor_location + pT.Stage + 
                         pN.Stage + RAS + BRAF + ECOG_PS + 
                         Margins + Site_Mets + Mets.type,
                       data = ps_dt,
                       method = "ps",      # Propensity score weighting
                       estimand = "ATE")   # Average Treatment Effect

# 4. Check Balance (The Love Plot)
# This will have a normal x-axis because we removed MSI
love=love.plot(weight_obj, 
          threshold = 0.1, 
          abs = TRUE, 
          main = "Covariate Balance After IPTW Weighting")

love
pdf('FigureS3a.pdf')
print(love, newpage = FALSE)
dev.off()


# 5. Prepare data for Weighted Analysis
ps_dt$weights <- weight_obj$weights

# 6. Analyze Outcomes (Weighted Kaplan-Meier)
# Re-adjust RFS to start from the landmark
ps_dt[, RFS_adj := RFS - landmark_medianACTend]

# Use the 'weights' argument in survfit
km_weighted <- survfit(Surv(RFS_adj, RFS.event) ~ Group, 
                       data = ps_dt, 
                       weights = weights)

# 7. Plot the Weighted Results
g=ggsurvplot(
  km_weighted,
  data = ps_dt,
  risk.table = TRUE,
  pval = TRUE,
  palette = c("orange", "purple"),
  title = "Landmark Analysis (7.5 mo) - IPTW Weighted Cohort",
  xlab = "Months from Landmark",
  legend.labs = c("SOC", "mFOLFOXIRI"),
  ggtheme = theme_minimal()
)
g
pdf('FigureS3b.pdf')
print(g, newpage = FALSE)
dev.off()

ps_dt[, Group := factor(Group, levels = c("SOC", "mFOLFOXIRI"))]

# Run the Weighted Cox Model
cox_weighted <- coxph(
  Surv(RFS_adj, RFS.event) ~ Group, 
  data = ps_dt, 
  weights = weights, 
  robust = TRUE
)

report_surv_extra(km_weighted, cox_weighted)


# 8. Weighted Cox Model (for Hazard Ratio)
# Must use robust standard errors (robust = TRUE) when using weights
weighted_cox <- coxph(Surv(RFS_adj, RFS.event) ~ Group, 
                      data = ps_dt, 
                      weights = weights, 
                      robust = TRUE)

summary(weighted_cox)


