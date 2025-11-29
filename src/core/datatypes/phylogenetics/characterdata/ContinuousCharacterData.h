#ifndef ContinuousCharacterData_H
#define ContinuousCharacterData_H

#include <set>
#include <vector>
#include <cstddef>
#include <iosfwd>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "ContinuousTaxonData.h"
#include "DiscretizedContinuousCharacterData.h"
#include "DistanceMatrix.h"
#include "HomologousCharacterData.h"
#include "MemberObject.h"
#include "Printable.h"

namespace RevBayesCore {
class AbstractCharacterData;
class DagNode;
    
    /**
     * Class for all Continuous character data objects.
     *
     * The continuous character data class is class for all continuous character data objects.
     * Here the data is stored via a set of ContinuousTaxonData objects.
     * Several convenience functions are provided to access, store and delete the data.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-02-16, version 1.0
     */
    class ContinuousCharacterData : public Printable, public HomologousCharacterData, public MemberObject<double> {

    public:
                                                        ContinuousCharacterData(void);                                              //!< Default constructor
        
        // Overloaded operators
        const ContinuousTaxonData&                      operator[](size_t i) const;                                                 //!< Subscript operator (const)
        ContinuousTaxonData&                            operator[](size_t i);                                                       //!< Subscript operator (const)
        
        // implemented methods of the Cloneable interface
        ContinuousCharacterData*                        clone(void) const;
        
        // implemented methods of the Serializable interface
        void                                            initFromFile( const path &dir, const std::string &fn );              //!< Read and resurrect this object from a file in its default format.
        void                                            initFromString( const std::string &s );                                     //!< Serialize (resurrect) the object from a string value
        void                                            writeToFile(const path &dir, const std::string &fn) const;           //!< Write this object into a file in its default format.
	json                                            toJSON() const;
        void                                            printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const;                                  //!< print object for user (in user-formatted way)
        void                                            printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten ) const;           //!< print object with standard rounding
        void                                            printForComplexStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten ) const;          //!< print object with maximum precision

        // CharacterData functions
        void                                            concatenate(const ContinuousCharacterData &d, std::string type = "");                              //!< Concatenate data matrices
        void                                            concatenate(const HomologousCharacterData &d, std::string type = "");                              //!< Concatenate data matrices
        void                                            concatenate(const AbstractCharacterData &d, std::string type = "");                                //!< Concatenate data matrices
        DiscretizedContinuousCharacterData*             discretizeCharacter(double error, std::int64_t num_bins, double cushion_width) const; //!< Discretize the ith character into k bins
        void                                            excludeAllCharacters(void);                                                 //!< Exclude all characters
        void                                            excludeCharacter(size_t i);                                                 //!< Exclude character
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;     //!< Map the member methods to internal function calls
        const double&                                   getCharacter(size_t tn, size_t cn) const;                                   //!< Return a reference to a character element in the character matrix
        std::string                                     getDataType(void) const;
        double                                          getMaxDifference(size_t index) const;                                       //!< Get the maximum difference between two observations
        double                                          getMaxSpeciesDifference(size_t index) const;                                //!< Get the maximum difference between two observations
        double                                          getMeanDifference(size_t index) const;                                      //!< Get the mean of the differences between two observations
        double                                          getMeanSpeciesDifference(size_t index) const;                               //!< Get the mean of the differences between two observations
        double                                          getMeanValue(size_t index) const;                                           //!< Get the mean of the differences between two observations
        double                                          getMinDifference(size_t index) const;                                       //!< Get the minimum difference between two observations
        double                                          getMinSpeciesDifference(size_t index) const;                                //!< Get the minimum difference between two observations
        size_t                                          getNumberOfCharacters(void) const;                                          //!< Number of characters
        size_t                                          getNumberOfIncludedCharacters(void) const;                                  //!< Number of characters
        DistanceMatrix                                  getPairwiseSpeciesDifference( size_t char_index ) const;                    //!< Get the pairwise species distance matrix
        double                                          getSpeciesMean( size_t species, size_t site ) const;                        //!< Get the mean value of the species for column 'i'
        ContinuousTaxonData&                            getTaxonData(size_t tn);                                                    //!< Return a reference to a sequence in the character matrix
        const ContinuousTaxonData&                      getTaxonData(size_t tn) const;                                              //!< Return a reference to a sequence in the character matrix
        ContinuousTaxonData&                            getTaxonData(const std::string &tn);                                        //!< Return a reference to a sequence in the character matrix
        const ContinuousTaxonData&                      getTaxonData(const std::string &tn) const;                                  //!< Return a reference to a sequence in the character matrix
        double                                          getVarDifference(size_t index) const;                                       //!< Get the variance of the differences between two observations
        double                                          getVarSpeciesDifference(size_t index) const;                                //!< Get the variance of the differences between two observations
        double                                          getVarValue( size_t index ) const;                                          //!< Get the variance of the taxon values for column 'i'
        double                                          getWithinSpeciesVariance( size_t species, size_t site ) const;              //!< Get the variance within the species for column 'i'
        void                                            includeCharacter(size_t i);                                                 //!< Include character
        bool                                            isCharacterExcluded(size_t i) const;                                        //!< Is the character excluded
        bool                                            isCharacterResolved(size_t txIdx, size_t chIdx) const;                      //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        bool                                            isCharacterResolved(const std::string &tn, size_t chIdx) const;             //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        void                                            removeExcludedCharacters(void);                                             //!< Remove all the excluded characters
        void                                            restoreCharacter(size_t i);                                                 //!< Restore character

        
    protected:
        std::set<size_t>                                deleted_characters;                                                          //!< Set of deleted characters

    };
    
    
}

#endif

