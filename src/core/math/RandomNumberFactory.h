

#ifndef RandomNumberFactory_H
#define RandomNumberFactory_H

#include <set>

namespace RevBayesCore {

    #define GLOBAL_RNG RandomNumberFactory::randomNumberFactoryInstance().getGlobalRandomNumberGenerator()
//    #define NEW_RNG    RandomNumberFactory::randomNumberFactoryInstance().getRandomNumberGenerator()

    class RandomNumberGenerator;

    /**
     * @brief RandomNumberFactory class declaration
     * The class RandomNumberFactory is
     * used to manage random number generating objects. The class has a pool
     * of random number objects that it can hand off as needed. This singleton
     * class has two seeds it manages: one is a global seed and the other is
     * is a so called local seed.
     *
     */
    class RandomNumberFactory {

	public:
		static RandomNumberFactory&                 randomNumberFactoryInstance(void)                                                      //!< Return a reference to the singleton factory
                                                    {
                                                        static RandomNumberFactory singleRandomNumberFactory;
                                                        return singleRandomNumberFactory;
                                                    }
		void                                        deleteRandomNumberGenerator(RandomNumberGenerator* r);                                 //!< Return a random number object to the pool
		RandomNumberGenerator*                      getGlobalRandomNumberGenerator(void) { return seedGenerator; }                         //!< Return a pointer to the global random number object

	private:
                                                    RandomNumberFactory(void);                                                             //!< Default constructor
                                                    RandomNumberFactory(const RandomNumberFactory&);                                       //!< Copy constructor
                                                    RandomNumberFactory& operator=(const RandomNumberFactory&);                            //!< Assignment operator
                                                   ~RandomNumberFactory(void);                                                             //!< Destructor
		RandomNumberGenerator*                      seedGenerator;                                                                         //!< A random number object that generates seeds
		std::set<RandomNumberGenerator*>            allocatedRandomNumbers;                                                                //!< The pool of random number objects
    };
}

#endif


