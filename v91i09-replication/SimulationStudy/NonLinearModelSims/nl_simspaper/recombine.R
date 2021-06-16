files <- list.files("out", full.names = TRUE)
get_indices <- function(files) {
  split <- strsplit(files, "out/out")
  indices <- unlist(lapply(split, function(x) {
    x[2]
  }))
  indices <- as.numeric(indices)
  indices <- indices[is.na(indices) == FALSE]
  indicies <- sort(indices)
}

#The function takes numeric indices arranged in alphabetical
#order and reorders them according the numeric order.
reorder_files <- function(files, n){
  bad_ordering <- sort(as.character(1:n))
  files <- files[order(as.numeric(bad_ordering))]
  files
}
files <- reorder_files(files, 400)

indices <- get_indices(files)
rf250 <- matrix(NA, 1, 10)
ff250 <- matrix(NA, 1, 10)
cif250 <- matrix(NA, 1, 10)
rf500 <- matrix(NA, 1, 10)
ff500 <- matrix(NA, 1, 10)

for (i in 201:400) {
  load(files[i])
  rf500 <- rbind(rf500, out$rf)
  ff500 <- rbind(ff500, out$ff)
}
for (i in 1:200) {
  load(files[i])
  rf250 <- rbind(rf250, out$rf)
  ff250 <- rbind(ff250, out$ff)
  cif250 <- rbind(cif250, out$cif)
}
rf500 <- rf500[-1, ]
ff500 <- ff500[-1, ]
rf250 <- rf250[-1, ]
ff250 <- ff250[-1, ]
cif250 <- cif250[-1, ]

feat_names <- paste("V", c(1:5, 76:80), sep = "")
rf500 <- apply(rf500, 2, sum)/1000
ff500 <- apply(ff500, 2, sum)/1000
ff250 <- apply(ff250, 2, sum)/1000
rf250 <- apply(rf250, 2, sum)/1000
cif250 <- apply(cif250, 2, sum)/1000

n250 <- vector("list", 3)
n250[[1]] <- ff250
n250[[2]] <- rf250
n250[[3]] <- cif250

n500 <- vector("list", 2)
n500[[1]] <- ff500
n500[[2]] <- rf500

nl_sims <- list(n250, n500)
save(nl_sims, file = "nl_sims.RData")
