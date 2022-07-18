#ifndef RandomNumberGenerator_H
#define RandomNumberGenerator_H

#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace RevBayesCore {

    class RandomNumberGenerator {

    public:

                                                    RandomNumberGenerator(void);                            //!< Default constructor using time seed
                                            
        // Regular functions
        unsigned int                                getNewSeed(void) const;                                 //!< Get the new seed values
        unsigned int                                getSeed(void) const;                                    //!< Get the seed values
        void                                        setSeed(unsigned int s);                                //!< Set the seeds of the RNG
        double                                      uniform01(void);                                        //!< Get a random [0,1) var
        double                                      exponential(void);                                      //!< Get a random variate from the standard exponential distribution

    private:
        
        double                                      last_u;
        boost::uniform_01<boost::mt19937>           zeroone;
        unsigned int seed;

    };
    
}

#endif

