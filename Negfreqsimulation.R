#Negative frequency dependent selection simulation

#key functions:
#   newpop: generate a new population
#   popimmunity: generate an immune response and remove pathogens
#   pathpopreplicate: pathogens replicate within hosts
#   hostsurvival: hosts with too many pathogens die
#   hostpopreplicate: copy survivng hosts
#   poptransmit: pool and redistribute the pathogens in the population
#   round: one round of this process

#Naming system: start with host type letter and one letter of aplphabet (e.g. HA)
#Every time a host or pathogen is copied, is gets an addition to its name (e.g. HA -> HA1)
# (HA1 and HA -> HA2 and HA11) (HA, HA1, HA11, HA2 -> HA3, HA12, HA111, HA21)

library(rlist)
library(dplyr)
library(ggplot2)

unbind <- function(seq) {unlist(strsplit(seq, ""))}
rebind <- function(chars) {paste(chars, collapse = "")}

#====Create hosts====

seqgen <- function(len) {sample(letters[1:4], len, replace = T)}
seqs <- function (n, len) {lapply(1:n, function(n) seqgen(len))}

sysnamegen <- function(orgtype, n) {
  paste(rep(orgtype, n), LETTERS[1:n], sep = "")
}

newhost <- function(lenseqhost){
  seqs <- seqs(1,lenseqhost)
  return(seqs)
}

newpaths <- function(npath, lenseqpath, hostname) {
  seqs <- seqs(npath, lenseqpath)
  names(seqs) <- sysnamegen(paste("P", hostname, sep = ""), npath)
  return (seqs)
}

wholehost <- function(listnum, lenseqhost, npath, lenseqpath){
  host <- newhost(lenseqhost)
  paths <- newpaths(npath, lenseqpath, LETTERS[listnum])
  host <- list.append(host, paths)
  names(host) <- c("host", "path")
  return(host)
}

newpop <- function(nhosts, npath, lenseqhost, lenseqpath) { #*** key function ***
  newpop <- lapply(1:nhosts, function(nhosts) wholehost(nhosts, lenseqhost, npath, lenseqpath))
  names(newpop) <- sysnamegen("H", nhosts)
  return(newpop)
}

examplepop <- newpop(5, 5, 20, 10)

#====Hosts get rid of some pathogens====

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

taggedhost <- Bcell(examplepop$HA)

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

#====Replication machinery====

proofreader <- function() {
  mean <- 0.5
  sd <- 0.2
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

#====Names====

nestedlistnames <- function(nestedlist) {
  names <- names(nestedlist)
  for (elem in nestedlist){
    if(is.list(elem)) {names <- c(names, names(elem))}
    for (elem2 in elem) {
      if (is.list(elem2)) {names <- c(names, names(elem2))}}}
  return(names)
}

repnamegen <- function(name, pop, gennum) {
  newnames <- paste(rep(name, gennum), 1:gennum, sep = "")
  newnames[min(which(!(newnames %in% nestedlistnames(pop))))]
}


repnamegen2 <- function(names, pop, gennum) {
  newnames <- NULL
  for (i in names) {
    newnames <- c(newnames, repnamegen(i, pop, gennum))
  }
  return(newnames)
}

newnamegen <- function(hostname, initpop) {
  newnames <- paste(rep(rebind(c("P", tail(unbind(hostname),nchar(hostname)-1))), 26), LETTERS, sep = "")
  newnames[(which(!(newnames %in% nestedlistnames(initpop)))[1:2])]
}

#==========Pathogens replicate within hosts==========

fiss <- function(host, pop, gennum){
  newpath <- lapply(host$path, pol)
  names(newpath) <- repnamegen2(names(host$path), pop, gennum)
  host$path <-  c(host$path, newpath)
  return(host)
}

pathgenerationhost <- fiss(decolonisedhost, examplepop, 3)

pathpopreplicate <- function(pop, gennum) { # *** key function ***
  pop2 <- pop
  lapply(pop, function(pop) fiss(pop, pop2, gennum))
}

pathgenerationpop <- pathpopreplicate(decolonisedpop, 2) 

#==========Hosts die==========

hostdies <- function(host, maxlen) {
  length(host$path) > maxlen
}

hostsurvival <- function(pop, maxlen) { # *** key function ***
  list.clean(pop, function(pop) hostdies(pop, maxlen))
}

survivingpop <- hostsurvival(pathgenerationpop, 5)

#==========Hosts replicate========== 

anypaths <- function(host) {
  !length(host$path)==0
}

rebuildhost <- function(host, paths, hostname, pathnames) {
  if(!(is.na(pathnames[1]))){
    names(paths) <- pathnames
  }
  newhost <- list(list.append(host, paths))
  names(newhost) <- hostname
  names(newhost[[1]]) <- c("host", "path")
  return(newhost)
}

replenishpop <- function(pop, pathlen) {
  empties <- pop[which(!(sapply(pop, anypaths)))]
  for (i in 1:length(empties)){
    current <- empties[i]
    pop[which(names(pop) == names(current))] <- rebuildhost(current[[1]][1], seqs(2, pathlen), names(current), 
                                                                        newnamegen(names(current), pop))
  }
  return(pop)
}


replicater <- function(host, pop, gennum){
  copy <- list(pol(host$host), host$path)
  names(copy) <- c("host", "path")
  names(copy$path) <- repnamegen2(names(copy$path), pop, gennum)
  return(copy)
}

hostpopreplicate <- function(pop, gennum, pathlen) { # *** keyfunction ***
  if(length(which(!(sapply(pop, anypaths))))!=0) {pop <- replenishpop(pop, pathlen)}
  pop2 <- pop
  newhosts <- lapply(pop, function(pop) replicater(pop, pop2, gennum))
  names(newhosts) <- repnamegen2(names(newhosts), pop, gennum)
  pop <- c(pop, newhosts)
  return(pop)
}

copiedpop <- hostpopreplicate(survivingpop, 2, 10)

#==========Transmission==== 

genpairindex <- function(nhosts, type) {
  if (type == "rand" | type == 0) {
    matrix(sample(1:nhosts, nhosts, replace = F), nrow = nhosts/2, byrow = T)}
  else if(type == "partrand" | type == 1) {
    matrix(sample(1:nhosts, nhosts, prob = dpois(1:nhosts, 1)), nrow = nhosts/2, byrow = T)}
  else if(type == "nonrand" | type == 2) {
    matrix(1:nhosts, nrow = nhosts/2, byrow = T)}
}

pathsamp <- function(host) {
  pathlist <- host$path
  list.sample(pathlist, replace = F)
}

rebuilder <- function(host1, host2, pathpool) {
  newhost1 <- rebuildhost(host1[[1]][1], 
                          pathpool[1:(length(pathpool)/2)], 
                          names(host1), pathnames = NA)
  newhost2 <- rebuildhost(host2[[1]][1], 
                          pathpool[((length(pathpool)/2)+1):length(pathpool)], 
                          names(host2), pathnames = NA)                          
  newhosts <- c(newhost1, newhost2)
  return(newhosts)
}

pairtransmit <- function(pop, index1, index2){
  host1 <- pop[index1]
  host2 <- pop[index2]
  pool <-list.sample(c(pathsamp(host1[[1]]), pathsamp(host2[[1]])), replace = F)
  rebuilder(host1, host2, pool)
}

poptransmit <- function(pop, type) { ### Key function ###
  n <- length(pop)/2
  pairs <- genpairindex(length(pop), type)
  newpop <- NULL
  for (i in 1:n){
    newpop <- c(newpop, pairtransmit(pop, pairs[i, 1], pairs[i , 2]))
  }
  return(newpop)
}

poptransmit(copiedpop, 0)

#==========Repeat rounds==========

run <- function(pop, generation, hoststrength, mixing, pathlen){
  pop <- pop %>%
    popimmunity() %>%
    pathpopreplicate(generation + 2) %>%
    hostsurvival(hoststrength)
    if(length(pop)==0 | length(pop)>100) {
      print(c("Hosts all dead or too many hosts."))
    } else {
      pop %>%
      hostpopreplicate(generation + 3, pathlen) %>%
      poptransmit(mixing)
    }
}

main <- function(numgens, initpop, hoststrength, mixing, pathlen) {
  gen <- list(inpoitpop)
  names(gen) <- "gen001"
  for(i in 1:numgens) {
    if(gen[[i]][1]=="Hosts all dead or too many hosts.") {return(gen)}
    else {
    gen <- list.append(gen, run(gen[[i]], i+1, hoststrength, mixing, pathlen))
    names(gen)[i+1] <- rebind(c("gen", as.character(stri_pad_left(i+1, "0", width = 3))))
    }
  }
  return(gen)
}

pop <- newpop(nhosts = 20, npath =  30, lenseqhost = 15, lenseqpath = 10)
test <- main(numgens = 20, initpop = pop, hoststrength = 50, mixing = 0, pathlen = 10)
gensize <- data.frame(generation = names(test), gensize = sapply(test, length))
ggplot(gensize, aes(generation, gensize)) + geom_bar(stat = "identity")






