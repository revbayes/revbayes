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

    private:
        
        double                                      last_u;
        boost::uniform_01<boost::mt19937>           zeroone;
        unsigned int seed;

    };
}

namespace deprecated
{
    template< class RandomIt >
    void random_shuffle( RandomIt first, RandomIt last )
    {
        typename std::iterator_traits<RandomIt>::difference_type i, n;
        n = last - first;
        for (i = n-1; i > 0; --i) {
            using std::swap;
            swap(first[i], first[std::rand() % (i+1)]);
            // rand() % (i+1) isn't actually correct, because the generated number
            // is not uniformly distributed for most values of i. A correct implementation
            // will need to essentially reimplement C++11 std::uniform_int_distribution,
            // which is beyond the scope of this example.
        }
    }
}


#endif

