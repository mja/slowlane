library(plyr)
library(nadiv)

# get the resident members from the pedigree
resident <- function(ped) {

        subset(ped, is.na(death) & is.na(emigrated))
}

# given the current pedigree return a list of living males
# and females of breeding age
living.and.breeding <- function(current, primiparity, t) {
	
	breeding.f <- subset(current, sex == 0 & (t - birth) >= primiparity[1]) # breeding females 
	breeding.m <- subset(current, sex == 1 & (t - birth) >= primiparity[2])
	
	return(list(f=breeding.f$id, m=breeding.m$id))
}

# return birth_rate * n random subset of females who will
# breed this season
will.breed <- function(breeding.f, birth_rate) {
	n <- length(breeding.f)
	breeding.f[rbinom(n, 1, prob=birth_rate) == 1]
}

# pick a mate for a give female
pick.mate <- function(female, males, A, inbreeding_tol) { # A is the relationship matrix
	pair <- NULL
	relatedness <- A[as.character(female),as.character(males)]
	relatedness[relatedness > inbreeding_tol] <- 1 # if relatedness is greater than inbreeding tolerance
	try(pair <- data.frame(female, male=sample(males, 1, p=1 - relatedness)), silent=TRUE)
									   # sample 1 male      exclude related males
	
	return(pair)
}

breed <- function(pairs, t, last.id) {
	
	n <- dim(pairs)[1]
	return(data.frame(id=(last.id+1):(last.id+n),
			   dam=pairs$female, sire=pairs$male,
			   birth=t, death=NA, emigrated=NA, enter=t, exit=NA,
			   sex=rbinom(n, 1, 0.5)))
}

# P[age+1] gives survival
survival.fun <- function(ages, P) {

	i <- ages + 1
	# use P[last] for all ages greater than length P
	i[i > length(P)] <- length(P)
	
	rbinom(length(i), 1, P[i]) == 0
}

#' Generate a random pedigree of given depth and troop size
#' 
#' Given demographic parameters, this function will simulate
#' birth, mating, death, and migration events within a troop
#' over a specified number of seasons, then return the pedigree 
#' of the troop.
#' 
#' The returned pedigree contains columns for individual, dam, and
#' sire IDs as well as information on sex (\code{0 = female, 1 = 
#' male}), and season in which the individual was born, died, and
#' emigrated. If \code{death = NA} and \code{emigrated = NA}, then
#' the individual is still in the troop at the end of the
#' simulation.
#'
#' @param founders number of female, male founders as a vector.
#'   Founders all start at age of sexual maturity.
#' @param capacity carrying capacity of the troop. When the troop
#'   troop size is > 1.5 * capacity group fission occurs and a 
#'   random half of the population emigrates.
#' @param im_rate vector of immigration rate for females, males.
#'   Default simulates a matrilocal society.
#' @param em_rate vector of emigration rates for females, males.
#'   Default simulates a matrilocal society
#' @param primiparity age of sexual maturation for females, males.
#' @param birth_rate probability of a female reproducing each season
#' @param seasons number of seasons (years) to simulate
#' @param inbreeding_tol inbreeding tolerance. Individuals will not
#'   mate with another individuals whose 2 * coefficient of 
#'   coancestry is >= this value.
#' @param morality vector of age based mortalities, where 
#'    \code{mortality[i]} is the mortality rate at age i - 1.
#'    If the vector is shorter than the age of a given individual,
#'    the last value in the vector is recycled. The default 
#'    assumes infant (age = 0) mortality rate of .3 and of .04
#'    for every age after that. 
#' @export
#' @examples
#' library(pedantics)
#' troop <- pedgen(c(50, 50), 140, seasons=20)
#' pedigree <- troop[,1:3]
#' phen <- phensim(pedigree, randomA=1, randomE=1)
#' 
pedgen <- function(founders=c(35, 35),
				   capacity=70,
				   im_rate=c(0.0, 0.05),
				   em_rate=c(0, 0.05),
           primiparity=c(5,5),
				   birth_rate = .4,
				   seasons=100,
				   inbreeding_tol=0.1,
				   mortality=c(rep(.05, 10), seq(.05, 1, by=.1))) {
					

	# create the initial population
	ped <- data.frame(id=1:sum(founders), dam=NA, sire=NA, sex=rep(c(0, 1), times=founders), birth=rep(-primiparity, times=founders), death=NA, emigrated=NA, enter=0, exit=NA)

	# simulate for n seasons
	for(t in 1:seasons) {
		
   	current <- resident(ped)
    n <- dim(current)[1]
		A <- diag(1, nrow=n, ncol=n)
		try(A <- makeA(current[,1:3]), silent=TRUE) # fails if parents are all NA
		rownames(A) <- colnames(A) <- current$id
		
		############
		# breeding #
		############
	
		breeding <- living.and.breeding(current, primiparity, t)

    # ratio of males to females
    sex_ratio <- length(breeding$m) / length(breeding$f)
		
		# select subset of females who will breed
		breeding.f <- will.breed(breeding$f, birth_rate)
		
		# pick a male for each female
		pairs <- adply(breeding.f, .margins=1, .fun=pick.mate, males=breeding$m, A=A, inbreeding_tol=inbreeding_tol)
		
		if(dim(pairs)[1] != 0) { # there are some breeding pairs.
			
			new.generation <- breed(pairs, t, max(ped$id))
			
			ped <- rbind(ped, new.generation)
			
		}
		
		#############
		# migration #
		#############

    # choose emigrants
    # from among breeding individuals
    # adjust for sex ratio (restrict to [0, 1]
    procrust <- function(x, a, b) { if(x < a) return(a); if(x > b) return(b); return(x)}
    em_rate.f <- procrust(em_rate[1] * 1/sex_ratio, 0, 1)
    em_rate.m <- procrust(em_rate[2] * sex_ratio, 0, 1)
    emigrant.f <- breeding$f[rbinom(length(breeding$f), 1, em_rate.f) == 1]
    emigrant.m <- breeding$m[rbinom(length(breeding$m), 1, em_rate.m) == 1]
    
    emigrants <- c(emigrant.f, emigrant.m)
    # individuals only leave if the troop is at capacity
    if(length(emigrants) >= 1 & n > capacity) {
      ped[ped$id %in% emigrants,]$emigrated <- t
      ped[ped$id %in% emigrants,]$exit <- t
    }
		
		# number of new migrants
    # depends on sex ratio
    im_rate.f <- procrust(im_rate[1] * sex_ratio, 0, 1)
    im_rate.m <- procrust(im_rate[2] * 1/sex_ratio, 0, 1)
		immigrant.f <- rbinom(1, round(capacity/2), im_rate.f)
		immigrant.m <- rbinom(1, round(capacity/2), im_rate.m)
		
		last_id = max(ped$id)
		no_immigrants = c(immigrant.f, immigrant.m)
		if(sum(no_immigrants) >= 1) {
			immigrants <- data.frame(id=(last_id + 1):(last_id + sum(no_immigrants)),
									 dam=NA, sire=NA,
									 sex=rep(c(0, 1), times=no_immigrants),
									 birth=rep(t - primiparity, times=no_immigrants), # for now assume new migrants
									 death=NA,				# are exactly at breeding age
                   emigrated=NA, enter=t, exit=NA)
								
			ped <- rbind(ped, immigrants)
		}


		#############
		# mortality #
		#############
		
		survives <- survival.fun(t - current$birth, mortality)
    
		try(ped[ped$id %in% current$id[!survives],]$death <- ped[ped$id %in% current$id[!survives],]$exit <- t, silent=TRUE)
		# because the expression fails if !survives are all FALSE
		
		###########
		# fission #
		###########

    # reassess residence
    current <- resident(ped)
		
		# if the population is 50% over carrying capacity
		pop.size <- dim(current)[1]
		
		if(pop.size > 1.1 * capacity) {
			emigrants <- sample(current$id, pop.size-.9*capacity, replace=FALSE)
			
			ped[ped$id %in% emigrants,]$emigrated <- t
			ped[ped$id %in% emigrants,]$exit <- t
		}
	}	
	
	return(ped)
					
}

resident_in_season <- function(pop, season) subset(pop, enter <= season & (season <= exit | is.na(exit)))


season_ply <- function(pop, .fun) {

        seasons <- sort(unique(c(pop$enter, pop$exit)))

        calc <- sapply(seq_along(seasons), function(s) .fun(resident_in_season(pop, s)))

        names(calc) <- seasons

        return(calc)
}

