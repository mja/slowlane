\name{pedgen}
\alias{pedgen}
\title{Generate a random pedigree of given depth and troop size}
\usage{
  pedgen(founders = c(20, 20), capacity = 70,
    im_rate = c(0, 0.1), em_rate = c(0, 0.25),
    primiparity = c(5, 5), birth_rate = 0.3, seasons = 100,
    inbreeding_tol = 0.1, mortality = c(0.3, 0.04))
}
\arguments{
  \item{founders}{number of female, male founders as a
  vector.  Founders all start at age of sexual maturity.}

  \item{capacity}{carrying capacity of the troop. When the
  troop troop size is > 1.5 * capacity group fission occurs
  and a random half of the population emigrates.}

  \item{im_rate}{vector of immigration rate for females,
  males.  Default simulates a matrilocal society.}

  \item{em_rate}{vector of emigration rates for females,
  males.  Default simulates a matrilocal society}

  \item{primiparity}{age of sexual maturation for females,
  males.}

  \item{birth_rate}{probability of a female reproducing
  each season}

  \item{seasons}{number of seasons (years) to simulate}

  \item{inbreeding_tol}{inbreeding tolerance. Individuals
  will not mate with another individuals whose 2 *
  coefficient of coancestry is >= this value.}

  \item{morality}{vector of age based mortalities, where
  \code{mortality[i]} is the mortality rate at age i - 1.
  If the vector is shorter than the age of a given
  individual, the last value in the vector is recycled. The
  default assumes infant (age = 0) mortality rate of .3 and
  of .04 for every age after that.}
}
\description{
  Given demographic parameters, this function will simulate
  birth, mating, death, and migration events within a troop
  over a specified number of seasons, then return the
  pedigree of the troop.
}
\details{
  The returned pedigree contains columns for individual,
  dam, and sire IDs as well as information on sex (\code{0
  = female, 1 = male}), and season in which the individual
  was born, died, and emigrated. If \code{death = NA} and
  \code{emigrated = NA}, then the individual is still in
  the troop at the end of the simulation.
}
\examples{
library(pedantics)
troop <- pedgen(c(50, 50), 140, seasons=20)
pedigree <- troop[,1:3]
phen <- phensim(pedigree, randomA=1, randomE=1)
}

