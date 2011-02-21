require(pedantics)
require(plyr)

# given a pedigree return a list of living males and females of breeding age
living.and.breeding <- function(ped, primiparity, t) {
	
	breeding.f <- subset(ped, is.na(death) & sex == 0 & (t - birth) >= primiparity[1]) # breeding females 
	breeding.m <- subset(ped, is.na(death) & sex == 1 & (t - birth) >= primiparity[2])
	
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
			   birth=t, death=NA,
			   sex=rbinom(n, 1, 0.5)))
}

# Generate a random pedigree of given depth and population size
genped <- function(founders=c(20, 20),
				   capacity=70,
				   im_rate=c(0, 0.05),
				   em_rate=c(0, 0.1),
				   im_age=c(7, 7),
				   em_age=c(7, 7),
				   primiparity=c(5,5),
				   birth_rate = .3,
				   seasons=100,
				   inbreeding_tol=0.1,
				   survival_fun) {
					

	# create the initial population
	ped <- data.frame(id=1:sum(founders), dam_id=NA, sire_id=NA, sex=rep(c(0, 1), times=founders), birth=rep(-primiparity, times=founders), death=NA)

	# simulate for n seasons
	for(t in 1:seasons) {
		
		############
		# breeding #
		############
	
		breeding <- living.and.breeding(ped, primiparity, t)
		
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
									 death=NA)				# are exactly at breeding age
								
			ped <- rbind(ped, immigrants)
		}
		
		#############
		# mortality #
		#############
		
		
		
	}	
	
	return(ped)
					
}
