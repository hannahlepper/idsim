#Negative frequency dependent selection simulation

#key functions:
#   newpop: generate a new population
#   popimmunity: generate an immune response and remove pathogens
#   pathpopreplicate: pathogens replicate within hosts
#   hostsurvival: hosts with too many pathogens die
#   hostpopreplicate: copy survivng hosts
#   poptransmit: pool and redistribute the pathogens in the population
#   round: one round of this process

library(rlist)

unbind <- function(seq) {unlist(strsplit(seq, ""))}
rebind <- function(chars) {paste(chars, collapse = "")}

#==========Create hosts==========

seqgen <- function(len) {sample(letters[1:3], len, replace = T)}
seqs <- function (n, len) {lapply(1:n, function(n) seqgen(len))}
namegen <- function(orgtype, n) {
  paste(rep(orgtype, n), sample(LETTERS, n, replace = F), sample(10:99, n, replace = F), sep = "")
}

newhost <- function(npath, lenseqhost, lenseqpath){
  seqs <- c(seqs(1,lenseqhost), list(seqs(npath, lenseqpath)))
  names(seqs) <- c("host", "path")
  names(seqs$path) <- namegen("p", npath)
  return(seqs)
}

newpop <- function(nhosts, npath, lenseqhost, lenseqpath) { #*** key function ***
  newpop <- lapply(1:nhosts, function(nhosts) newhost(npath, lenseqhost, lenseqpath))
  names(newpop) <- namegen("host", nhosts)
  return(newpop)
}

examplepop <- newpop(5, 5, 20, 10)

#==========Hosts get rid of some pathogens==========

snips <- function (tosnip, sniplen, start) {rebind(tosnip[start:(start+sniplen-1)])}

snipper <- function(chars, len) {
  snipnum <- length(chars) - len + 1
  lapply(1:snipnum, function (snipnum) snips(chars, len, snipnum))
}

tags <- function(seq, immun, path) {
  if(immun %in% path){c(seq, "z")} else {return(seq)}
}

antibody <- function (immun, path) {
  immunsnips <- snipper(immun, 4)
  pathsnips <- snipper(path, 4)
  for (i in immunsnips){
    path <- tags(seq = path, immun = i, path = pathsnips)}
  return(path)
}

Bcell <- function(host){
  n <- length(host$path)
  for (i in 1:n) {
    host$path[[i]]<-antibody(immun = host$host, path = host$path[[i]])
  }
 return(host)
}

taggedhost <- Bcell(examplepop$hostM65)

spotstag <- function(pathseq) {
  if ("z" %in% pathseq) {return(TRUE)} else (return(FALSE))}

CD8 <- function(host) {
  host <- list(host$host, list.clean(host$path, spotstag))
  names(host) <- c("host", "path")
  return(host)
}

decolonisedhost <- CD8(taggedhost)

popimmunity <- function(pop) {# *** key function ***
  lapply(lapply(pop, Bcell), CD8)
}

decolonisedpop <- popimmunity(examplepop)

#==========Replication machinery==========

proofreader <- function() {
  mean <- 0.5
  sd <- 0.1
  error <- rnorm(n = 1, mean = mean, sd = sd)
  error > (mean - sd) & error < (mean + sd)
}

# rbinom(1, T, 0.67)

basebind <- function(char) {
  if (proofreader()){return(char)}
  else (sample(letters[1:3],1))
}

pol <- function(seq) {
  sapply(seq, basebind, USE.NAMES = F)
}

#==========Pathogens replicate within hosts==========

fiss <- function(host){
  newpath <- lapply(host$path, pol)
  host$path <- c(newpath, newpath)
  return(host)
}

pathgenerationhost <- fiss(decolonisedhost)

pathpopreplicate <- function(pop) { # *** key function ***
  lapply(pop, fiss)
}

pathgenerationpop <- pathpopreplicate(decolonisedpop)

#==========Hosts die==========

hostdies <- function(host) {
  length(host$path) > 10
}

hostsurvival <- function(pop) { # *** key function ***
  list.clean(pop, hostdies)
}

survivingpop <- hostsurvival(pathgenerationpop)

#==========Hosts replicate==========

replicater <- function(host){
  list(host = pol(host$host), path = host$path)
}

hostpopreplicate <- function(pop) { # *** keyfunction ***
  newhosts <- lapply(pop, replicater)
  names(newhosts) <- namegen("host", length(newhosts))
  pop <- c(pop, newhosts)
  return(pop)
}

copiedpop <- hostpopreplicate(survivingpop)

#==========Transmission==========

pathextract <- function(host) {
  list(host$path)
}

hoststrip <- function(host) {
  host$path <- NA
  return(host)
}

pathpool <- function(pop) {
  lapply(pop, pathextract)
}

pathpool1 <- pathpool(copiedpop)

popstrip <- function(pop) {lapply(pop, hoststrip)}

strippedpop <- popstrip(copiedpop)

reassignpathsnums <- function(pop, pathpool) {
  sample <- 0
  while (!all(1:length(pop) %in% sample)){
    sample <- sample(x = 1:length(pop), size = length(list.ungroup(pathpool)), replace = T)
  }
  return(sample)
}

assigns <- reassignpathsnums(strippedpop, pathpool1)

assignpathtohost <- function(host, hostnum, pathpool, assignmentnums) {
  paths <- pathpool[which(assignmentnums==hostnum)]
  host$path <- paths
  return(host)
}

assignpathtohost(strippedpop$hostM65, 1, pathpool1, assigns)

popassign <- function(pop, pathpool, assignmentnums) {
  n <- length(pop)
  popnames <- names(pop)
  pop <- lapply(1:n, function(n)
    assignpathtohost(pop[[n]], n, pathpool, assignmentnums))
  names(pop) <- popnames
  return(pop)
}

popassign(strippedpop, list.ungroup(pathpool1), assigns)

poptransmit <- function(pop) { # *** key function ***
  pathpool1 <- pathpool(pop)
  assignnums <- reassignpathsnums(pop, pathpool1)
  pop <- popassign(pop, pathpool1, assignnums)
  return(pop)
}

freshpop <- poptransmit(copiedpop)

#==========One round==========

round <- function(pop) { # *** key function ***
  poptransmit(hostpopreplicate(hostsurvival(pathpopreplicate(popimmunity(pop)))))
}

start <- Sys.time()
round(newpop(5, 1, 10, 20, 10))
Sys.time() - start

#==========Repeat rounds==========

rounds <- function(startpop, nrounds){
  i <- 0
  pop <- startpop
  repeat {
    pop <- poptransmit(hostpopreplicate(hostsurvival(pathpopreplicate(popimmunity(pop)))))
    print(c("round", i, ":", pop))
    i <- i+1
    if (i >= nrounds) {break}
  }
  return(pop)
}

test <- rounds(newpop(5, 10, 20, 10),2)

#BUGS: Extra list level in $path after transmission

pop <- newpop(5, 10, 20, 10)
pop2 <- popimmunity(pop)
pop3 <- pathpopreplicate(pop)
pop4 <- hostsurvival(pop)
pop5 <- hostpopreplicate(pop)
pop6 <- poptransmit(pop)
