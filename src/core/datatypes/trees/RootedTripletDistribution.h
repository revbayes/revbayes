/**
 * @file
 * This file contains the declaration of a rooted triplet distribution,
 * that holds a vector of triplets along with their frequencies.
 *
 *
 * @brief Declaration of the RootedTripletDistribution
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: RootedTripletDistribution.h 1651 2014-08-31 14:47:08Z Boussau $
 */

#ifndef RootedTripletDistribution_H
#define RootedTripletDistribution_H

#include <map>
#include <vector>
#include <string>
#include <cstddef>
#include <iosfwd>
#include <utility>

#include "Cloneable.h"
#include "MemberObject.h"
#include "Taxon.h"
#include "Tree.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


namespace RevBayesCore {
class DagNode;
class TopologyNode;
    
    class RootedTripletDistribution : public Cloneable, public MemberObject<std::vector < std::string > > {
        
    public:
    
        RootedTripletDistribution(void);                                                                           //!< Default constructor
        RootedTripletDistribution( const RbVector<Tree>& ts, const std::vector< std::string > spNames, bool recordBls = false ) ;
        RootedTripletDistribution( const RbVector< std::string > spNames ) ;
        RootedTripletDistribution( const RbVector<Tree>& ts, const std::vector< Taxon > tax, bool recordBls = false ) ;
        RootedTripletDistribution( const std::vector< Taxon > tax ) ;
        virtual                                    ~RootedTripletDistribution(void);                                                                                                        //!< Destructor
        
        // Basic utility functions
        RootedTripletDistribution*                                                              clone(void) const;                                                      //!< Clone object
        
        void                                                                                    executeMethod(const std::string &n, const std::vector<const DagNode *> &args, std::vector < std::string > &rv) const ;                                 //!< Execute the member-method
        
        // RootedTripletDistribution functions
        void                                                                                    extractTriplets( const Tree& t );                           //!< Getting the triplets from a tree and adding them to the tripletDistribution
        void                                                                                    extractTriplets(  );                                        //!< Getting the triplets from all trees and adding them to the tripletDistribution
        void                                                                                    populateTripletDistribution ( const TopologyNode* node, std::vector< size_t >& allTips ) ;
        void                                                                                    addAllTriplets(std::vector< size_t >& leftTips, std::vector< size_t >& rightTips) ; //!< Get all rooted triplets given vectors of left and right tips
        void                                                                                    addAllTripletsOneWay( std::vector< size_t >& leftTips, std::vector< size_t >& rightTips, size_t leftSize,size_t rightSize ); //!< Get rooted triplets given vectors of left and right tips, one way only

        
        void                                                                                    populateTripletDistribution ( const Tree& tree, const TopologyNode* node, std::vector< size_t >& allTips, std::vector< double >& distancesToTips, std::map < std::pair < size_t, size_t >, double >& distancesToAncestors ) ;
        void                                                                                    addAllTriplets(std::vector< size_t >& leftTips, std::vector< size_t >& rightTips, std::vector< double >& leftDistancesToTips, std::vector< double >& rightDistancesToTips, std::map < std::pair < size_t, size_t >, double >& distancesToAncestors ) ; //!< Get all rooted triplets given vectors of left and right tips, and keep distances
        void                                                                                    addAllTripletsOneWay( std::vector< size_t >& leftTips, std::vector< size_t >& rightTips, size_t leftSize, size_t rightSize,
                                                                                                                     std::vector<double>& leftDistancesToTips,
                                                                                                                     std::vector<double>& rightDistancesToTips,
                                                                                                                     std::map < std::pair < size_t, size_t >, double >& distancesToAncestors ); //!< Get rooted triplets given vectors of left and right tips, one way only, and keep distances
        size_t                                                                                  getNumberOfTrees() const;                                                                     //!< Get the number of trees that were used to build the object
        size_t                                                                                  getNumberOfTriplets() const;                                                                  //!< Get the number of triplets in the object
        std::map < std::pair < size_t, std::pair < size_t, size_t > >, size_t >                 getTriplets() const;         //!< Get triplets, with their total numbers of occurences
        std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, size_t >                    getTaxonTriplets() const;       //!< Get triplets of taxa, with their total numbers of occurences
        std::map < std::pair < std::string, std::pair < std::string, std::string > >, size_t >  getSpeciesTriplets() const; //!< Get triplets of species, with their total numbers of occurences
        std::map < std::pair < Taxon, std::pair < Taxon, Taxon > >, std::vector<double> >       getTaxonTripletsWithBranchLengths() const;       //!< Get triplets of taxa, with their total numbers of occurences
        std::map < std::pair < std::string, std::pair < std::string, std::string > >, std::vector<double> >     getSpeciesTripletsWithBranchLengths() const; //!< Get triplets of species, with their total numbers of occurences

        void                                                                                    setSpecies ( std::vector< std::string > s );
        void                                                                                    setTaxa ( std::vector< Taxon > t );
        void                                                                                    setTrees ( const RbVector< Tree >& ts );
        void                                                                                    setRecordBranchLengths ( const bool b );
        std::vector< std::string >                                                              getSpecies ( ) const;
        std::vector< Taxon >                                                                    getTaxa ( ) const;
        std::string                                                                             getSpecies ( size_t i ) const;
        Taxon                                                                                   getTaxon ( size_t i ) const;
        void                                                                                    resetDistribution (  ) ;
        const                                                                                   std::string getStringRepresentation() const ;
        
 //   protected:
        // Nothing
        
    private:

        
        //Members:
        RbVector< Tree >                                                                        trees;
        std::vector< Taxon >                                                                    taxa;
        std::vector< std::string >                                                              species;
        std::map < std::pair < size_t, std::pair < size_t, size_t > >, size_t >                 tripletDistribution;
        std::map < std::pair < size_t, std::pair < size_t, size_t > >, std::vector<double> >    tripletDistributionAndBranchLength;
        size_t                                                                                  numberOfTrees;
        std::map< Taxon, size_t >                                                               taxonToIndex;
        std::map< std::string, size_t >                                                         speciesToIndex;
        bool                                                                                    speciesOnly;
        bool                                                                                    recordBranchLengths;
        
        
    };
    
    // Global functions using the class
    std::ostream& operator<<(std::ostream& o, const RootedTripletDistribution& x);       //!< Overloaded output operator
    
}


#endif

