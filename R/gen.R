require(pedantics)
require(plyr)

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
pick.suitor <- function(female, males, A=NULL) { # possible future argument: A matrix
	return(data.frame(female, male=sample(males, 1)))
}

breed <- function(pairs, t, last.id) {
	
	n <- dim(pairs)[1]
	return(data.frame(id=(last.id+1):(last.id+n),
			   dam_id=pairs$female, sire_id=pairs$male,
			   birth=t, death=NA, emigrated=NA,
			   sex=rbinom(n, 1, 0.5)))
}

# P[age+1] gives survival
survival.fun <- function(ages, P) {

	i <- ages + 1
	# use P[last] for all ages greater than length P
	i[i > length(P)] <- length(P)
	
	rbinom(length(i), 1, P[i]) == 0
}

# Generate a random pedigree of given depth and population size
genped <- function(founders=c(20, 20),
				   capacity=70,
				   im_rate=c(0, 0.05),
				   em_rate=c(0, 0.1),
           primiparity=c(5,5),
				   birth_rate = .3,
				   seasons=100,
				   inbreeding_tol=0.1,
				   mortality=c(.4, .04)) {
					

	# create the initial population
	ped <- data.frame(id=1:sum(founders), dam_id=NA, sire_id=NA, sex=rep(c(0, 1), times=founders), birth=rep(-primiparity, times=founders), death=NA, emigrated=NA)

	# simulate for n seasons
	for(t in 1:seasons) {
		
    current <- resident(ped)

		############
		# breeding #
		############
	
		breeding <- living.and.breeding(current, primiparity, t)
		
		# select subset of females who will breed
		breeding.f <- will.breed(breeding$f, birth_rate)
		
		# pick a male for each female
		pairs <- adply(breeding.f, .margins=1, .fun=pick.suitor, males=breeding$m)
		
		if(dim(pairs)[1] != 0) { # there are some breeding pairs.
			
			new.generation <- breed(pairs, t, max(ped$id))
			
			ped <- rbind(ped, new.generation)
			
		}
		
		#############
		# migration #
		#############
		
		# number of new migrants
		immigrant.f <- rbinom(1, capacity, im_rate[1])
		immigrant.m <- rbinom(1, capacity, im_rate[2])
		
		last_id = max(ped$id)
		no_immigrants = c(immigrant.f, immigrant.m)
		if(sum(no_immigrants) >= 1) {
			immigrants <- data.frame(id=(last_id + 1):(last_id + sum(no_immigrants)),
									 dam_id=NA, sire_id=NA,
									 sex=rep(c(0, 1), times=no_immigrants),
									 birth=rep(t - primiparity, times=no_immigrants), # for now assume new migrants
									 death=NA,				# are exactly at breeding age
                   emigrated=NA)
								
			ped <- rbind(ped, immigrants)
		}

    # choose emigrants
    # from among breeding individuals
    emigrant.f <- sample(breeding$f, length(breeding$f)*em_rate[1], replace=FALSE)
    emigrant.m <- sample(breeding$m, length(breeding$f)*em_rate[2], replace=FALSE)
    
    emigrants <- c(emigrant.f, emigrant.m)
    if(length(emigrants) >= 1) {
      ped[ped$id %in% emigrants,]$emigrated <- t
    }

		#############
		# mortality #
		#############
		
		survives <- survival.fun(t - current$birth, mortality)
    
		try(ped[ped$id %in% current$id[!survives],]$death <- t, silent=TRUE)
		# because the expression fails if !survives are all FALSE
		
		##########
		# fusion #
		##########
		
		# if the population is 50% over carrying capacity
		pop.size <- dim(current)[1]
		
		if(pop.size > 1.5 * capacity) {
			emigrants <- sample(current$id, pop.size/2, replace=FALSE)
			
			ped[ped$id %in% emigrants,]$emigrated <- t
		}
	}	
	
	return(ped)
					
}
