/**
 * @file
 * This file contains the declaration of the random variable class for MPEst pseudo-likelihood approach to compute the score of a species tree given gene trees (actually, the score of a rooted triplet distribution given another rooted triplet distribution).
 * This class is derived from the stochastic node and each instance will represent a random variable.
 *
 * @brief Declaration of the constant rate Birth-Death process class.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date:$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-17, version 1.0
 * @interface TypedDagNode
 *
 * $Id:$
 */

#ifndef MPEST_H
#define MPEST_H

#include "RootedTripletDistribution.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class Clade;
    
    class MPEST : public TypedDistribution<RootedTripletDistribution> {
        
    public:
        MPEST(const TypedDagNode<RootedTripletDistribution> *st, /*const TypedDagNode<RootedTripletDistribution> *gt, */bool useSp);
        virtual                                            ~MPEST(void);                                                      //!< Virtual destructor
        
        // public member functions
        MPEST*                                              clone(void) const;                                                                 //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);          //!< Swap a parameter
        
    private:
        
        
        // helper functions
        void                                                computeLnW();                                                       //!< Compute the log of the w factor, used to compute the loglk
        void                                                extractGeneTriplets() ;                                           //!< Extract gene triplets
        double                                              computeMultinomialLogLk(double num1, double num2, double num3, double bl); //!< Compute the multinomial loglk of a triplet
        
        // members
        const TypedDagNode<RootedTripletDistribution>*                                                 speciesTree;
        double                                                                                             logProb;
        bool                                                                                            useSpecies;
        double                                                                                                 lnW;
        std::map < std::pair < std::string, std::pair < std::string, std::string > >, size_t > geneSpeciesTriplets;
        std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, size_t >                    geneTaxonTriplets ;
        //const TypedDagNode<RootedTripletDistribution>*    geneTrees;

    };
    
}

#endif
