/**
 * @file
 * This file contains the declaration of RandomNumberFactory, which is
 * used to manage random number generating objects. The class has a pool
 * of random number objects that it can hand off as needed. This singleton
 * class has two seeds it manages: one is a global seed and the other is
 * is a so called local seed.
 *
 * @brief Declaration of RandomNumberFactory
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#ifndef RandomNumberFactory_H
#define RandomNumberFactory_H

#include <set>
#include <vector>

namespace RevBayesCore {

    #define GLOBAL_RNG RandomNumberFactory::randomNumberFactoryInstance().getGlobalRandomNumberGenerator()
//    #define NEW_RNG    RandomNumberFactory::randomNumberFactoryInstance().getRandomNumberGenerator()

    class RandomNumberGenerator;
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


