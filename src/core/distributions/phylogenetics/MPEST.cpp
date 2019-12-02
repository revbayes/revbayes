//
//  MPEST.cpp
//  RevBayes
//
//  Created by Bastien on 26/09/14.
//  Copyright (c) 2014 Bastien Boussau. All rights reserved.
//

#include "Clade.h"
#include "ConstantNode.h"
#include "DistributionExponential.h"
#include "MPEST.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

MPEST::MPEST(const TypedDagNode<RootedTripletDistribution> *st, /*const TypedDagNode<RootedTripletDistribution> *gt, */bool useSp) : TypedDistribution<RootedTripletDistribution>( NULL ),
    speciesTree( st ),
//geneTrees( gt ),
    logProb (0.0),
    useSpecies(useSp),
    lnW (0.0),
    geneSpeciesTriplets(),
    geneTaxonTriplets()
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( speciesTree );
    
    
}



MPEST::~MPEST()
{
    
}


MPEST* MPEST::clone( void ) const
{
    
    return new MPEST( *this );
}


double MPEST::computeLnProbability( void )
{
    
    // variable declarations and initialization
    double lnProbCoal = 0;
    
    // Here, we need to go through all the rooted triplets in the speciesTree,
    // and for each rooted triplet compute its multinomial likelihood
    // according to the geneTrees rooted triplet distribution.
    // For instance, for triplet (A,(B,C):bl) in the species tree,
    // we need to find the three triplets (A,(B,C)), (B,(A,C)), (C,(B,A)) in the
    // geneTrees rooted triplet distribution, and get their frequencies.
    // Then these frequencies are used to compute the multinomial likelihood,
    // using the branch length bl from the species tree.
    
    const RootedTripletDistribution &sp = speciesTree->getValue();
    //    std::cerr << sp << std::endl;
    
    if (useSpecies) {
        std::map < std::pair < std::string, std::pair < std::string, std::string > >, std::vector<double> > spTriplets = sp.getSpeciesTripletsWithBranchLengths();
        for ( std::map < std::pair < std::string, std::pair < std::string, std::string > >, std::vector<double> >::const_iterator it = spTriplets.begin(); it != spTriplets.end(); ++it ) {
            size_t num1 = geneSpeciesTriplets[it->first];
            std::pair < std::string, std::string > duo1 (it->first.first, it->first.second.first);
            std::pair < std::string, std::pair < std::string, std::string > > alternatePair1 (it->first.second.second, duo1);
            size_t num2 = geneSpeciesTriplets[alternatePair1];
            std::pair < std::string, std::string > duo2 (it->first.first, it->first.second.second);
            std::pair < std::string, std::pair < std::string, std::string > > alternatePair2 (it->first.second.first, duo2);
            size_t num3 = geneSpeciesTriplets[alternatePair2];
            // Now we have the 3 counts, and their total, we can compute the multinomial lk.
            // First, get one branch length out of the vector of branch lengths. We choose the mean.
            double bl = 0.0;
            for (size_t i = 0; i < it->second.size(); ++i) {
                bl+=it->second[i];
            }
            bl = bl/(double) (it->second.size());
            //we have the branch length
            lnProbCoal += computeMultinomialLogLk( (double) num1, (double) num2, (double) num3, bl);
        }
    }
    else {
        std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, std::vector<double> > taxTriplets = sp.getTaxonTripletsWithBranchLengths();
        for ( std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, std::vector<double> >::const_iterator it = taxTriplets.begin(); it != taxTriplets.end(); ++it ) {
            size_t num1 = geneTaxonTriplets[it->first];
            std::pair < Taxon, Taxon > duo1 (it->first.first, it->first.second.first);
            std::pair < Taxon, std::pair < Taxon, Taxon > > alternatePair1 (it->first.second.second, duo1);
            size_t num2 = geneTaxonTriplets[alternatePair1];
            std::pair < Taxon, Taxon > duo2 (it->first.first, it->first.second.second);
            std::pair < Taxon, std::pair < Taxon, Taxon > > alternatePair2 (it->first.second.first, duo2);
            size_t num3 = geneTaxonTriplets[alternatePair2];
            // Now we have the 3 counts, and their total, we can compute the multinomial lk.
            // First, get one branch length out of the vector of branch lengths. We choose the mean.
            double bl = 0.0;
            for (size_t i = 0; i < it->second.size(); ++i) {
                bl+=it->second[i];
            }
            bl = bl/(double) (it->second.size());
            //we have the branch length
            lnProbCoal += computeMultinomialLogLk( (double) num1, (double) num2, (double) num3, bl);
        }
    }
    
    return lnProbCoal;
    
}


void MPEST::computeLnW()
{
    
    if (useSpecies)
    {
        std::map < std::pair < std::string, std::pair < std::string, std::string > >, std::vector<double> > spTriplets = speciesTree->getValue().getSpeciesTripletsWithBranchLengths();
        for ( std::map < std::pair < std::string, std::pair < std::string, std::string > >, std::vector<double> >::const_iterator it = spTriplets.begin(); it != spTriplets.end(); ++it ) {
            size_t num1 = geneSpeciesTriplets[it->first];
            std::pair < std::string, std::string > duo1 (it->first.first, it->first.second.first);
            std::pair < std::string, std::pair < std::string, std::string > > alternatePair1 (it->first.second.second, duo1);
            size_t num2 = geneSpeciesTriplets[alternatePair1];
            std::pair < std::string, std::string > duo2 (it->first.first, it->first.second.second);
            std::pair < std::string, std::pair < std::string, std::string > > alternatePair2 (it->first.second.first, duo2);
            size_t num3 = geneSpeciesTriplets[alternatePair2];
            size_t tot = num1+num2+num3;
            lnW = RbMath::lnFactorial((int)tot) - RbMath::lnFactorial((int)num1) - RbMath::lnFactorial((int)num2) - RbMath::lnFactorial((int)num3);
        }
    }
    else
    {
        std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, std::vector<double> > taxTriplets = speciesTree->getValue().getTaxonTripletsWithBranchLengths();
        for ( std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, std::vector<double> >::const_iterator it = taxTriplets.begin(); it != taxTriplets.end(); ++it ) {
            size_t num1 = geneTaxonTriplets[it->first];
            std::pair < Taxon, Taxon > duo1 (it->first.first, it->first.second.first);
            std::pair < Taxon, std::pair < Taxon, Taxon > > alternatePair1 (it->first.second.second, duo1);
            size_t num2 = geneTaxonTriplets[alternatePair1];
            std::pair < Taxon, Taxon > duo2 (it->first.first, it->first.second.second);
            std::pair < Taxon, std::pair < Taxon, Taxon > > alternatePair2 (it->first.second.first, duo2);
            size_t num3 = geneTaxonTriplets[alternatePair2];
            size_t tot = num1+num2+num3;
            lnW = RbMath::lnFactorial((int)tot) - RbMath::lnFactorial((int)num1) - RbMath::lnFactorial((int)num2) - RbMath::lnFactorial((int)num3);
        }

    }
    return;
}


void MPEST::extractGeneTriplets() {
    if (useSpecies) {
        geneSpeciesTriplets = value->getSpeciesTriplets();
    }
    else {
        geneTaxonTriplets = value->getTaxonTriplets();
    }
    return;
}


double MPEST::computeMultinomialLogLk(double num1, double num2, double num3, double bl) {
    double lnProbCoal = num1*log(1-(2/3)*exp(-bl));
    lnProbCoal += (num2+num3)*log(1-(1/3)*exp(-bl));
    return lnProbCoal;
}


void MPEST::redrawValue( void ) {
    
    std::cout << "Sorry, cannot simulate with the MPEst distribution." <<std::endl;
    
}

/** Swap a parameter of the distribution */
void MPEST::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == speciesTree)
    {
        speciesTree = static_cast<const TypedDagNode< RootedTripletDistribution > * >( newP );
    }
   /* else if ( oldP == geneTrees)
    {
        geneTrees = static_cast<const TypedDagNode< RootedTripletDistribution > *  >( newP );
        extractGeneTriplets();
        computeLnW();
    }*/
    
}
