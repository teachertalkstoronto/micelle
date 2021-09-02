//myrandom.h

#ifndef MYRANDOM_H
#define MYRANDOM_H

#include <ctime> 
#include <vector>

//Uses the boost package.
#include "/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/boost/random/mersenne_twister.hpp"
#include "/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/boost/random/normal_distribution.hpp"
#include "/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/boost/random/uniform_real_distribution.hpp"
#include "/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/boost/random/uniform_on_sphere.hpp"
#include "/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/boost/random/uniform_int_distribution.hpp"
#include "/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/boost/random/variate_generator.hpp"


typedef boost::random::mt19937 ENG;	//Mersenne Twister Algorithm. Good uniform distribution in up to 623 dimensions.

typedef boost::normal_distribution<double> N_DIST;	            //Counts outcomes of (infinitely) repeated Bernoulli experiments
typedef boost::variate_generator<ENG, N_DIST> N_GEN;	        //The random number generator.

typedef boost::random::uniform_real_distribution<double> U_DIST;	    //Continuous uniform distribution on some range [min, max) of real numbers
typedef boost::variate_generator<ENG, U_DIST> U_GEN;	        //The random number generator.

typedef boost::uniform_on_sphere<double> U_SPH_DIST;	        //Uniform distribution on a unit sphere of arbitrary dimension
typedef boost::variate_generator<ENG, U_SPH_DIST> U_SPH_GEN;	//The random number generator.

typedef boost::random::uniform_int_distribution<int> U_INT_DIST;	    //Continuous uniform distribution on some range [min, max) of real numbers
typedef boost::variate_generator<ENG, U_INT_DIST> U_INT_GEN;	        //The random number generator.

#endif
