#ifndef AbstractHomologousDiscreteCharacterData_H
#define AbstractHomologousDiscreteCharacterData_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "DistanceMatrix.h"
#include "HomologousCharacterData.h"
#include "MatrixReal.h"

namespace RevBayesCore {
class AbstractCharacterData;
class AbstractDiscreteTaxonData;
class AbstractTaxonData;
class DiscreteCharacterState;
    
    /**
     * Abstract class for all DISCRETE character data objects.
     *
     * The DISCRETE character data class is derived from the interface character data
     * and simply specifies that all derived classes must contain discrete characters.
     * This simplifies return values because it is known that all are derived from discrete characters
     * and there are several function that only work on discrete characters and not on
     * continuous characters.
     *
     * Note that this class is a pure interface and thus contains only pure virtual functions!
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-02-16, version 1.0
     */
    class AbstractHomologousDiscreteCharacterData : public HomologousCharacterData {
        
    public:
        enum SFS_AMBIGUITY_TREATMENT { ANCESTRAL, DERIVED, RESCALE, SKIP_COLUMN };

        virtual                                                ~AbstractHomologousDiscreteCharacterData(void) {}
        
        // Overloaded operators
        virtual const AbstractTaxonData&                        operator[](size_t i) const = 0;                                                             //!< Subscript operator (const)
        
        // methods of the Cloneable interface
        virtual AbstractHomologousDiscreteCharacterData*        clone(void) const = 0;

        // methods of the Serializable interface
        virtual void                                            initFromString( const std::string &s ) = 0;                                                 //!< Serialize (resurrect) the object from a string value
        virtual void                                            writeToFile(const path &dir, const std::string &fn) const;
        
        // CharacterData functions
        void                                                    applyMissingSitesMask( const std::vector<std::vector<bool> >& mask_gap, const std::vector<std::vector<bool> >& mask_missing );
        virtual void                                            concatenate(const AbstractCharacterData &d, std::string type = "") = 0;                     //!< Concatenate data matrices
        virtual void                                            concatenate(const HomologousCharacterData &d, std::string type = "") = 0;                   //!< Concatenate two sequences
        virtual void                                            concatenate(const AbstractHomologousDiscreteCharacterData &d, std::string type = "") = 0;   //!< Concatenate data matrices
        virtual AbstractHomologousDiscreteCharacterData*        combineCharacters(const AbstractHomologousDiscreteCharacterData &d) const = 0;              //!< Combine/expand data matrices
        virtual double                                          computeMultinomialProfileLikelihood( void ) const = 0;
        virtual std::vector<long>                               computeSiteFrequencySpectrum(bool folded, SFS_AMBIGUITY_TREATMENT ambig_treat) const = 0;
        virtual MatrixReal                                      computeStateFrequencies(void) const = 0;                                                    //!< Compute the state frequencies for this character data object
        virtual void                                            excludeCharacter(size_t i) = 0;                                                             //!< Exclude character
        void                                                    fillMissingSitesMask( std::vector<std::vector<bool> >& mask_gap, std::vector<std::vector<bool> >& mask_missing ) const;
        virtual const DiscreteCharacterState&                   getCharacter(size_t tn, size_t cn) const = 0;                                               //!< Return a reference to a character element in the character matrix
        virtual std::string                                     getDataType(void) const = 0;                                                                //!< Return the data type of this character data matrix
        virtual std::vector<double>                             getEmpiricalBaseFrequencies(void) const = 0;                                                //!< Compute the empirical base frequencies
        virtual std::vector<size_t>                             getIncludedSiteIndices(void) const = 0;
        virtual std::vector<size_t>                             getInvariantSiteIndices(bool excl) const = 0;                                               //!< Get the indices of invariant characters in this matrix
        virtual size_t                                          getNumberOfCharacters(void) const = 0;                                                      //!< Number of characters
        virtual size_t                                          getMaxObservedStateIndex(void) const = 0;                                                   //!< Get the number of observed states for the characters in this matrix
        virtual size_t                                          getNumberOfSegregatingSites(bool excl) const = 0;                                           //!< Compute the number of segregating sites
        virtual size_t                                          getNumberOfStates(void) const = 0;                                                          //!< Get the number of states for the characters in this matrix
        virtual size_t                                          getNumberOfInvariantSites(bool excl) const = 0;                                             //!< Number of invariant sites
        virtual double                                          getAveragePairwiseSequenceDifference(bool excl) const = 0;                                   //!< Get the average pairwise sequence distance.
        virtual size_t                                          getMaxPairwiseSequenceDifference(bool excl) const = 0;                                       //!< Get the average pairwise sequence distance.
        virtual size_t                                          getMinPairwiseSequenceDifference(bool excl) const = 0;                                       //!< Get the average pairwise sequence distance.
        virtual DistanceMatrix                                  getPairwiseSequenceDifference(bool excl) const = 0;                                       //!< Get the average pairwise sequence distance.
        virtual AbstractDiscreteTaxonData&                      getTaxonData(size_t tn) = 0;                                                                //!< Return a reference to a sequence in the character matrix
        virtual const AbstractDiscreteTaxonData&                getTaxonData(size_t tn) const = 0;                                                          //!< Return a reference to a sequence in the character matrix
        virtual AbstractDiscreteTaxonData&                      getTaxonData(const std::string &tn) = 0;                                                    //!< Return a reference to a sequence in the character matrix
        virtual const AbstractDiscreteTaxonData&                getTaxonData(const std::string &tn) const = 0;                                              //!< Return a reference to a sequence in the character matrix
        virtual bool                                            isCharacterExcluded(size_t i) const = 0;                                                    //!< Is the character excluded
        virtual bool                                            isCharacterResolved(size_t txIdx, size_t chIdx) const = 0;                                  //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        virtual bool                                            isCharacterResolved(const std::string &tn, size_t chIdx) const = 0;                         //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
        
        virtual double                                          maxGcContent(bool excl) const = 0;                                                          //!< Maximum GC-content of a sequence
        virtual size_t                                          maxInvariableBlockLength(bool excl) const = 0;                                              //!< Maximum length of a block of invariant sites
        virtual size_t                                          maxVariableBlockLength(bool excl) const = 0;                                                //!< Maximum length of a block of variant sites
        virtual double                                          meanGcContent(bool excl) const = 0;                                                         //!< Mean GC-content of all sequence
        virtual double                                          meanGcContentByCodon(size_t n, bool excl) const = 0;                                        //!< Mean GC-content of all sequences by codon position
        virtual double                                          minGcContent(bool excl) const = 0;                                                          //!< Number of invariant sites
        virtual size_t                                          numInvariableSiteBlocks(bool excl) const = 0;                                               //!< Number of invariant sites
        virtual size_t                                          numberTaxaMissingSequence(double p) const = 0;                                              //!< Number of taxa missing x percent of the sequence
        virtual double                                          varGcContent(bool excl) const = 0;                                                          //!< Mean GC-content of all sequence
        virtual double                                          varGcContentByCodon(size_t n, bool excl) const = 0;                                         //!< Mean GC-content of all sequences by codon position
        
        void                                                    excludeMissingSites( void );
        void                                                    replaceRandomSitesByMissingData( double p );
        virtual void                                            removeExcludedCharacters(void) = 0;                                                         //!< Remove all the excluded characters
        virtual void                                            restoreCharacter(size_t i) = 0;                                                             //!< Restore character
        
        virtual AbstractHomologousDiscreteCharacterData*        expandCharacters(size_t n) const = 0;
        virtual AbstractHomologousDiscreteCharacterData*        translateCharacters(const std::string &type) const = 0;

    protected:
                                                                AbstractHomologousDiscreteCharacterData(void) {}                                            //!< Constructor requires character type
    };
    
    // Global functions using the class
    std::ostream&                                               operator<<(std::ostream& o, const AbstractHomologousDiscreteCharacterData& x);              //!< Overloaded output operator
    
}

#endif
