#ifndef RandomNumberGenerator_H
#define RandomNumberGenerator_H

#include <cstdint>

namespace RevBayesCore {

    /**
    * @brief
    * ### Pseudo Random Number Generator
    *
    * This class provides a pseudo random number generator using the xoshiro256** algorithm.
    *
    * The code for this class was adapted from the pure C implementation by David Blackman
    * and Sebastiano Vigna in the paper: *Scrambled Linear Pseudorandom Number Generators*.
    *
    * The original C implementation provided by Blackman and Vigna is licenced as public
    * domain. See <http://creativecommons.org/publicdomain/zero/1.0/>.
    *
    * More information can be found at http://prng.di.unimi.it.
    *
    * Note that in the original generator we are able to use uint64_t seeds.
    * To be compliant with the revbayes api, we instead use unsigned int seeds in this
    * implementation. Everything internal still uses uint64_t types.
    *
    */

    class RandomNumberGenerator {

        public:

            RandomNumberGenerator     ( void );            //!< Default constructor using time seed

            unsigned int  getSeed     ( void ) const;      //!< Get the seed values
            void          setSeed     ( unsigned int x );  //!< Set the seeds of the RNG
            double        uniform01   ( void );            //!< Get a random [0,1) var


        private:

            unsigned int  seed;                            //!< Seed for PRNG
            uint64_t      state[4];                        //!< Internal state vector of PRNG

            void          updateState ( void );            //!< Update the internal state based on seed
            uint64_t      next        ( void );            //!< Return next pseudo random number

    };

} //-- End RandomNumberGenerator namespace


#endif //-- RandomNumberGenerator.h



