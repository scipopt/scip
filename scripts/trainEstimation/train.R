#! /usr/bin/Rscript

# @author Gregor Hendel
#
# 1. Analyze accuracy of different tree size prediction methods
# 2. Train regression forests and linear regression to approximate search completion
# 3. Save the results as settings and RFCSV files that SCIP can read
#
# NOTE: this script should not be executed by itself, but within run_training.sh
#

#
# fix random seed to make train/test splits and random forest training reproducible
#
Random.Seed <- 123

#
# set this number to the (positive) report frequency used within SCIP
#
Report.Frequency <- 100

#
# necessary functions
#

#
# import library, but suppress the distracting startup messages
#
quibrary <- function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}
quibrary("readr")
quibrary("magrittr")
quibrary("ggplot2")
quibrary("dplyr")
quibrary("knitr")
quibrary("rpart")
quibrary("randomForest")
quibrary("reshape2")


#
# returns logical vector whether or not an input Ratio x is k-accurate.
#
iskaccurrate <- function(x,k) {
  (x >= (1/k)) & (x <= k)
}

#
# compute estimation from a search completion proxy
#
computeEstimation <- function(x,current) {
  x <- pmax(x,1e-6)
  x <- pmin(x, 1.0)

  current / x
}

#
# make strings prettier
#
printMethod <- function(x) {
  x %>% tolower() %>% gsub(".", " ", ., fixed = T)
}

#
# turn 1-based index positions for R into 0-based indices for C
#
zeroIndex <- function(x) {
  (x %>% as.integer()) - 1
}

#
# save regression forest in a file that SCIP can read.
#
turnIntoCSV <- function(rf.model, filename) {
  ntrees <- rf.model$ntree
  tree_as_data <- do.call(rbind, lapply(1:ntrees, function(x) {randomForest::getTree(rf.model, k = x)})) %>% as.data.frame()

  index_cols <- c("left daughter", "right daughter", "split var")
  tree_as_data[index_cols] <- lapply(tree_as_data[index_cols], zeroIndex)
  tree_as_data["node"] <- zeroIndex(rownames(tree_as_data))
  tree_as_data <- tree_as_data %>% mutate(value=ifelse(`split var` == -1, prediction, `split point`))

  "### NTREES=%d FEATURE_DIM=%d LENGTH=%d\n" %>%
    sprintf(ntrees, max(tree_as_data$`split var`) + 1, nrow(tree_as_data)) %>%
    cat(file = filename, append = FALSE)

  write.table(tree_as_data[c("node", "left daughter", "right daughter", "split var", "value")],
            file = filename,
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE,
            sep=",")
}

#
# function that summarizes records by custom grouping
# data must have columns 'ApproxError' and 'Ratio'
#
SummarizeErrors = function(data, ...) {
  data %>% group_by(.dots = lazyeval::lazy_dots(...)) %>%
    summarise(n=n(),
              "MSE"=mean(ApproxError),
              "MeanRatio" = 2**mean(abs(log2(Ratio))),
              "2Accurate" = mean(iskaccurrate(Ratio,2)),
              "3Accurate" = mean(iskaccurrate(Ratio,3)),
              "4Accurate" = mean(iskaccurrate(Ratio,4))
    )
}

#
# process command line arguments
#
args <- commandArgs(trailingOnly = TRUE)
if( length(args) < 1 ) {
    stop("Missing positional argument for output directory")
}

output.dir <- args[1]
actual.filename <- "%s/actual.csv" %>% sprintf(output.dir)
raw.table.filename <- "%s/table1.csv" %>% sprintf(output.dir)
search.completion.mse.plot.filename <- "%s/searchcompletion_mse.pdf" %>% sprintf(output.dir)
rf.model.filename <- "%s/rf_model.rfcsv" %>% sprintf(output.dir)
monotone.filename <- "%s/monotone.set" %>% sprintf(output.dir)

# read actual nodes
actual <- read.csv(actual.filename, header=FALSE, col.names = c("File", "Actual", "TotTime"), stringsAsFactors = FALSE)
row.names(actual) <- actual$File


# read estimations
table1 <- read_csv(raw.table.filename, col_names = c(
  "File",
  "Method",
  "Report",
  "Time",
  "Current",
  "Leaves",
  "Unsolved",
  "TreeWeight",
  "Estim",
  "Value",
  "Trend",
  "Resolution",
  "Smooth"
  ), col_types = cols(
    File=col_character(),
    Method=col_character(),
    .default = col_double())
  )
# combine them
table1$Actual <- actual[table1$File, "Actual"]
table1$TotTime <- actual[table1$File, "TotTime"]
table1$CorrEstim <- ifelse(table1$Estim < table1$Current, table1$Current + table1$Unsolved, table1$Estim)
table1$Ratio <- (table1$CorrEstim / table1$Actual)
table1$Level <- (table1$TreeWeight * Report.Frequency) %>% as.integer() / Report.Frequency
table1$RelTime <- ((table1$Time/table1$TotTime) * 100) %>% as.integer() / 100.0
table1$ApproxError <- NaN
table1$Prob <- basename(table1$File) %>% gsub("bzfhende.(miplib2017|MMMc).\\d+_([^.]*)", "\\2", ., perl = T)

#
# Step 2 Filtering the results. We use a fresh data frame
#

table2 <- table1 %>% dplyr::filter(Actual >= 100) # keep only instances that require at least 100 nodes
# in rare cases, there might be two observations recorded at the same level. Keep the first
table2 <- table2 %>% dplyr::distinct(Method,Prob,Level, .keep_all = TRUE) %>% as.data.frame()

estim.summary <- table2 %>% SummarizeErrors(Method)

#
# create a new data frame with feature data for training/testing
#
dataset <- data.frame(
  (table2 %>% filter(Method == "tree-weight") %>% select(Value)),
  (table2 %>% filter(Method == "tree-weight") %>% select(Trend)),
  (table2 %>% filter(Method == "ssg") %>% select(Value)),
  (table2 %>% filter(Method == "ssg") %>% select(Trend)),
  (table2 %>% filter(Method == "leaf-frequency") %>% select(Value)),
  (table2 %>% filter(Method == "leaf-frequency") %>% select(Trend)),
  (table2 %>% filter(Method == "gap") %>% select(Value)),
  (table2 %>% filter(Method == "gap") %>% select(Trend)),
  (table2 %>% filter(Method == "open-nodes") %>% select(Trend)) %>% { . < 0 }
)

colnames(dataset) <- c("TreeWeight.Value", "TreeWeight.Trend",
                       "Ssg.Value", "Ssg.Trend",
                       "Leaffreq.Value", "Leaffreq.Trend",
                       "Gap.Value", "Gap.Trend",
                       "OpenTrend"
                       )

# create
singletable <- table2 %>% filter(Method == "ssg")
singletable <- singletable %>% dplyr::mutate(SearchCompletion = Current / Actual)
searchCompletion <- singletable$SearchCompletion

set.seed(Random.Seed)
isTrain <- (1:nrow(dataset)) %in% sample(1:nrow(dataset), 0.8 * nrow(dataset))
trainingset <- dataset[isTrain,]
testset <- dataset[!isTrain,]

#
# learn monotone linear regression model
#
linear.monotone <- lm(searchCompletion[isTrain]~TreeWeight.Value+Ssg.Value,data = trainingset)

#
# normalize coefficients such that their sum is equal to 1
#
coeffs <- linear.monotone$coefficients[2:3]
# SSG needs to be reversed to 1 - SSG
coeffs[2] <- -coeffs[2]
coeffs <- pmax(coeffs,0.01)
coeffs <- pmin(coeffs,0.99)
normalized.coeffs <- coeffs + ((1 - sum(coeffs)) / 2)

#
# store monotone regression coefficients in a SCIP settings file
#
cat(file=monotone.filename, sprintf("estimation/coefmonoweight = %.4f\nestimation/coefmonossg = %.4f\n", normalized.coeffs[1], normalized.coeffs[2]))

#
# learn regression forest
#
set.seed(Random.Seed)
rf.model <- randomForest(searchCompletion[isTrain]~.,data = trainingset, ntree = 150, nodesize=50)

#
# store regression forest in CSV
#
turnIntoCSV(rf.model, rf.model.filename)

# combine predicted labels into test results data frame
p <- lapply(list(rf.model, linear.monotone), function(x) {predict(x, dataset)})
testresults <- do.call(cbind, p)
colnames(testresults) <- c("Random.Forest","linear.monotone")
testresults <- cbind(testresults,searchCompletion)
colnames(testresults)[3] <- "SearchCompletion"
testresults <- testresults %>% as.data.frame()
testresults <- testresults %>% dplyr::mutate(
                              "Linear.Easy.Estim"=computeEstimation(linear.monotone, singletable$Current),
                              "Random.Forest.Estim"=computeEstimation(Random.Forest, singletable$Current),
                                )

#
# compute an error table for the simple search completion measures gap,ssg, tree-weight, leaf frequency
#
approx.data <- dataset %>%
  mutate(treeweight=TreeWeight.Value,
         ssg= 1 - Ssg.Value,
         gap=Gap.Value,
         "leaf-frequency"=2 * pmax(Leaffreq.Value, 0.0)) %>%
  select(treeweight, ssg, gap, `leaf-frequency`) %>%
  reshape2::melt(value.name = "Approx")

approx.data <- approx.data %>%
  mutate(ApproxError=(Approx - singletable$SearchCompletion) ** 2,
         Estim=computeEstimation(Approx, singletable$Current),
         Ratio=Estim / singletable$Actual,
         Training=(rep(isTrain, 4) %>% ifelse("Training", "Test"))
  )

#
# summarize the obtained accuracy for the search completion approximations
#
simple.search.completion.errors <- approx.data %>%
  SummarizeErrors(variable) %>%
  as.data.frame() %>%
  mutate(Group="SearchCompletion")
colnames(simple.search.completion.errors)[1] <- "Method"

#
# treat WBE and tree profile estimation separately
#
custom.methods <- c("tree-profile", "wbe", "treeprofile") # we wrote tree-profile w/o hyphen in an earlier version

custom.method.errors <- estim.summary %>%
  filter(Method %in% custom.methods) %>%
  as.data.frame() %>%
  mutate(Group="Custom")

#
# compute errors for the remaining double exponential smoothing methods
#
double.exponential.errors <- estim.summary %>%
  filter(!Method %in% custom.methods) %>%
  as.data.frame() %>%
  mutate(Group="Forecast")

# learned search completion methods (monotone) linear regression and regression forest
learned.data.search.completion.1 <- testresults %>%
  select(linear.monotone,
        Random.Forest) %>%
  reshape2::melt(value.name = "Approx")

learned.data.search.completion.1 <- learned.data.search.completion.1 %>%
  mutate(ApproxError=(Approx - singletable$SearchCompletion) ** 2,
         Estim=computeEstimation(Approx, singletable$Current),
         Ratio=Estim / singletable$Actual,
         Training=(rep(isTrain, 2) %>% ifelse("Training", "Test"))
  )

learned.data.search.completion.1.errors <- learned.data.search.completion.1 %>%
  SummarizeErrors(variable) %>%
  as.data.frame() %>%
  mutate(Group="Learned")
colnames(learned.data.search.completion.1.errors)[1] <- "Method"


#
# create and output combined summary of the different methods
#
errors.summary.table <- rbind(
      simple.search.completion.errors,
      custom.method.errors,
      double.exponential.errors,
      learned.data.search.completion.1.errors
      )

cat("Methods ranked by geometric mean approximation ratio")
knitr::kable(errors.summary.table[order(errors.summary.table$MeanRatio),], digits=3) %>% print()


#
# create a plot to compare performance on training and test data
#
learned.data.search.completion.1.bytraining.errors <- learned.data.search.completion.1 %>%
  SummarizeErrors(variable, Training) %>%
  as.data.frame()
colnames(learned.data.search.completion.1.bytraining.errors)[1] <- "Method"

simple.search.completion.bytraining.errors <- approx.data %>%
  SummarizeErrors(variable, Training) %>%
  as.data.frame()
colnames(simple.search.completion.bytraining.errors)[1] <- "Method"


bytraining.errors <- rbind(simple.search.completion.bytraining.errors,
                           learned.data.search.completion.1.bytraining.errors)

bytraining.errors <- bytraining.errors %>% mutate(Method= Method %>% printMethod())
bytraining.errors$Training <- factor(bytraining.errors$Training)
bytraining.errors$Training <- factor(bytraining.errors$Training, levels=levels(bytraining.errors$Training)[c(2,1)])
bytraining.errors$Method <-factor(bytraining.errors$Method)
bytraining.errors$Method <- factor(bytraining.errors$Method, levels = levels(bytraining.errors$Method)[c(1,9,5,2,4,3,8,7,6)])

bytraining.errors <- bytraining.errors[order(bytraining.errors$Training),]
ssg.pos <- which(bytraining.errors$Method == "ssg") %>% rep(each=6)
bytraining.errors$MSEQ <- bytraining.errors$MSE / bytraining.errors$MSE[ssg.pos]
bytraining.errors <- bytraining.errors %>%
  mutate(MSELabel=sprintf("%.3f (x%.2f)", MSE, MSEQ))
ggplot(bytraining.errors, aes(Method, MSE, label=MSELabel)) +
  geom_col(position = "dodge") +
  geom_text(check_overlap = TRUE, position = position_dodge(1.0), hjust = -0.3) +
  facet_wrap(~Training) + coord_flip() +
  theme_light() + xlab(NULL) + ylim(0.0, 0.8)
ggsave(search.completion.mse.plot.filename, width = 8.5, height=5.5)
