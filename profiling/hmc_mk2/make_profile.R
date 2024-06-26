require(dplyr)
require(data.tree)
require(stringr)
require(rmarkdown)
require(treemap)
require(ggplot2) 
require(knitr) 
require(hadron) 
require(kableExtra) 
require(ggrepel)


# wrapper around base::gregexpr which makes it operate like microseq::gregexpr
extract_gregexpr <- function(pattern, text, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE){
  lst <- base::gregexpr(pattern, text, ignore.case, perl, fixed, useBytes)
  lapply(1:length(lst), function(i){substring(text[i], lst[[i]], lst[[i]]+attr(lst[[i]], "match.length") - 1)})
}

# extract the monomial name from the timing path string
extract_monomial <- function(str){
  l1 <- extract_gregexpr(text = str,
                         pattern = '(?<=/)(.*?)(?=:)',
                         perl = TRUE)

  unlist(extract_gregexpr(text = unlist(l1),
                          pattern = '(?<=/).*',
                          perl = TRUE))
}

# extract which part of the monomial cost the timing with this path string
# belongs to
extract_type <- function(str){
  tmp <- str
  tmp[!( grepl('acc', str) | grepl('heatbath',str) | grepl('derivative', str) )] <- 'other'
  tmp[grepl('acc', str)] <- 'acc'
  tmp[grepl('heatbath', str)] <- 'heatbath'
  tmp[grepl('derivative', str)] <- 'derivative'
  return(tmp)
}

if(!interactive()){
  args <- commandArgs(trailingOnly=TRUE)
  if( length(args) < 2 ){
    stop("usage: Rscript make_profile.R <logfile> <output_basename> [include_all]")
  }
  infile <<- args[1]
  outbase <<- args[2]
  if( length(args) == 3 ){
    other_threshold <<- 0.0
  }
} else {
  if( !exists("infile") | !exists("outbase") ){
    stop("When not running interactively, the variables 'infile' (path to log file) and 'outbase' (basename for ouput files) must be defined!")
  }
}

if( !exists("other_threshold") ){
  other_threshold <<- 0.05
}

stopifnot(file.exists(infile))

# extract the different monomial names from the log file
monomial_names <- system(paste("grep \"Initialised monomial\" ", infile, "| awk '{print $4}'"),
                         intern = TRUE)

# extract the various timings from the log file
raw_data <- system(paste("grep \"Time for\"",
                         infile, 
                         "| grep level",
                         "| awk '{print $2 \" \" $5 \" \" $6 \" \" $9 \" \" $12}'"),
                   intern = TRUE)

# extract QUDA's overall profiling data if present
quda_raw_data <- system(paste("grep -A8 \"QUDA Total time\"",
                              infile,
                              "| tail -n 8", # skip the first line
                              "| awk '{print $3 \" \" $5}'"),
                        intern = TRUE)

if( length(quda_raw_data) > 1 ){
  ## the format of the QUDA output has changed historically
  ## if the above has not worked, try again in the old format
  if( strsplit(quda_raw_data[1], split=' ')[[1]][1] != "download" ){
    quda_raw_data <- system(paste("grep -A8 \"QUDA Total time\"",
                                  infile,
                                  "| tail -n 8", # skip the first line
                                  "| awk '{print $1 \" \" $3}'"),
                            intern = TRUE)
  }
}

quda_data <- NA

stopifnot(length(raw_data) > 0)

data <- read.table(text = raw_data, stringsAsFactors = FALSE)
colnames(data) <- c("prefix", "name", "time", "level", "pathString")

if(length(quda_raw_data) > 0){
  quda_data <- read.table(text = quda_raw_data, stringsAsFactors = FALSE)
  colnames(quda_data) <- c("name", "time")
}

data <- dplyr::mutate(data,
                      ## R is 1-indexed, level starts at 0, hence 'level+1' is the parent level
                      ## because we have a leading forward-slash (word 1 is ""),
                      monomial = extract_monomial(pathString),
                      type = extract_type(pathString), 
                      parent = stringr::word(pathString, start = 1, end = level+1, sep='/'))

sum_data <- dplyr::group_by(data, pathString) %>%
            dplyr::summarise(time = sum(time),
                             invocations = n(),
                             unit_time = sum(time)/n(),
                             prefix = unique(prefix),
                             name = unique(name),
                             parent = unique(parent),
                             monomial = unique(monomial),
                             type = unique(type),
                             level = unique(level)) %>%
            dplyr::ungroup()

t_tree <- data.tree::as.Node(sum_data)

# at this level we want to produce summarised data so we collapse everthing
# which accounts for less than 5% of the total time at each level
# and logical unit
sum_data_per_mon <- dplyr::filter(sum_data, nchar(monomial) > 0) %>%
                    dplyr::group_by(monomial, type, level) %>%
                    dplyr::mutate(name = ifelse(time < other_threshold * sum(time),
                                                'other', name)) %>%
                    dplyr::group_by(monomial, type, level, name) %>%
                    dplyr::summarise(time = sum(time),
                                     invocations = sum(invocations),
                                     unit_time = sum(time)/sum(invocations),
                                     prefix = unique(prefix),
                                     name = unique(name),
                                     level = unique(level),
                                     type = unique(type),
                                     name = unique(name)) %>%
                    dplyr::ungroup() %>%
                    dplyr::group_by(monomial, type, level) %>%
                    dplyr::mutate(prop = 100 * time / sum(time)) %>%
                    dplyr::ungroup()

total_time <- dplyr::filter(sum_data, name == 'HMC')$time

type_per_mon <- dplyr::group_by(dplyr::filter(sum_data, 
                                              level == 1 & 
                                                nchar(monomial) > 0 & 
                                                !grepl('init', name) &
                                                !grepl('trlog', name) ),
                                monomial, name) %>%
                dplyr::summarise(time = sum(time),
                                 invocations = sum(invocations),
                                 unit_time = sum(time)/sum(invocations),
                                 type = unlist(extract_gregexpr(text = unique(name),
                                                                pattern = '(?<=_).*',
                                                                perl = TRUE))) %>%
                dplyr::ungroup()

# if we're analysing an incomplete run or one which is still
# in progress, then we estimate the total time by summing
# all the time spent at level 1
if( length(total_time) == 0 ){
  warning("total time could not be determined, summing up all times at level 1 instead")
  total_time <- sum(dplyr::filter(sum_data,level==1)$time)
}

total_per_mon <- dplyr::group_by(type_per_mon, monomial) %>%
                 dplyr::summarise(time = sum(time),
                                  invocations = sum(invocations),
                                  unit_time = sum(time)/sum(invocations)) %>%
                 dplyr::ungroup() %>%
                 dplyr::mutate(prop = 100 * time / sum(total_time)) %>%
                 dplyr::arrange(desc(prop))

hb_acc_der_data <- dplyr::filter(sum_data_per_mon, level == 1 &
                                  !grepl('init', name) &
                                  !grepl('trlog', name)) %>%
                   dplyr::group_by(monomial) %>%
                   dplyr::mutate(prop = 100 * (time / sum(time))) %>%
                   dplyr::mutate(prop_label = sprintf("%.1f %%", prop)) %>%
                   dplyr::ungroup()

save(monomial_names,
     raw_data,
     data,
     quda_data,
     sum_data,
     sum_data_per_mon,
     hb_acc_der_data,
     t_tree,
     total_per_mon,
     type_per_mon,
     total_time,
     file = "profile.RData")

rmarkdown::render("profile.Rmd")

message(sprintf("Renaming the output data and report to %s.[RData, pdf]", outbase))

# rename the report
system(sprintf("mv profile.RData %s.RData", outbase))
system(sprintf("mv profile.pdf %s.pdf", outbase))


