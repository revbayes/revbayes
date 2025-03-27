#ifndef HomologousDiscreteCharacterData_H
#define HomologousDiscreteCharacterData_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "DiscreteCharacterState.h"
#include "MatrixReal.h"
#include "DiscreteTaxonData.h"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

namespace RevBayesCore
{

  /**
   * Template class for all discrete character data objects.
   *
   * The discrete character data class is template class for all discrete character data objects.
   * Here the data is stored via a set of DiscreteTaxonData objects.
   * Several convenience functions are provided to access, store and delete the data.
   *
   * @copyright Copyright 2009-
   * @author The RevBayes Development Core Team (Sebastian Hoehna)
   * @since 2013-04-15, version 1.0
   */
  template <class charType>
  class HomologousDiscreteCharacterData : public AbstractHomologousDiscreteCharacterData
  {

  public:
    HomologousDiscreteCharacterData(); //!< Default constructor

    // Overloaded operators
    const DiscreteTaxonData<charType> &operator[](size_t i) const;  //!< Subscript operator (const)
    bool operator<(const HomologousDiscreteCharacterData &x) const; //!< Less than operator

    // implemented methods of the Cloneable interface
    HomologousDiscreteCharacterData<charType> *clone(void) const;

    // implemented methods of the Serializable interface
    void initFromFile(const path &dir, const std::string &fn); //!< Read and resurrect this object from a file in its default format.
    void initFromString(const std::string &s);                 //!< Serialize (resurrect) the object from a string value

    // CharacterData functions
    double computeMultinomialProfileLikelihood(void) const;
    std::vector<long> computeSiteFrequencySpectrum(bool folded, SFS_AMBIGUITY_TREATMENT ambig_treat) const;
    MatrixReal computeStateFrequencies(void) const;
    void concatenate(const HomologousDiscreteCharacterData &d, std::string type = "");         //!< Concatenate data matrices
    void concatenate(const AbstractCharacterData &d, std::string type = "");                   //!< Concatenate data matrices
    void concatenate(const HomologousCharacterData &d, std::string type = "");                 //!< Concatenate data matrices
    void concatenate(const AbstractHomologousDiscreteCharacterData &d, std::string type = ""); //!< Concatenate data matrices
    void excludeAllCharacters(void);                                                           //!< Exclude all characters
    void excludeCharacter(size_t i);                                                           //!< Exclude character
    const charType &getCharacter(size_t tn, size_t cn) const;                                  //!< Return a reference to a character element in the character matrix
    std::string getDataType(void) const;
    std::vector<double> getEmpiricalBaseFrequencies(void) const; //!< Compute the empirical base frequencies
    const std::set<size_t> &getExcludedCharacters(void) const;   //!< Returns the name of the file the data came from
    std::vector<size_t> getIncludedSiteIndices(void) const;
    std::vector<size_t> getInvariantSiteIndices(bool excl) const;                 //!< Get the indices of invariant characters in this matrix
    size_t getMaxObservedStateIndex(void) const;                                  //!< Get the number of observed states for the characters in this matrix
    size_t getNumberOfCharacters(void) const;                                     //!< Number of characters
    size_t getNumberOfIncludedCharacters(void) const;                             //!< Number of characters
    size_t getNumberOfInvariantSites(bool excl) const;                            //!< Number of invariant sites
    size_t getNumberOfSegregatingSites(bool excl) const;                          //!< Compute the number of segregating sites
    size_t getNumberOfStates(void) const;                                         //!< Get the number of states for the characters in this matrix
    double getAveragePairwiseSequenceDifference(bool excl) const;                 //!< Get the average pairwise sequence distance.
    size_t getMaxPairwiseSequenceDifference(bool excl) const;                     //!< Get the average pairwise sequence distance.
    size_t getMinPairwiseSequenceDifference(bool excl) const;                     //!< Get the average pairwise sequence distance.
    DistanceMatrix getPairwiseSequenceDifference(bool excl) const;                //!< Get the average pairwise sequence distance.
    DiscreteTaxonData<charType> &getTaxonData(size_t tn);                         //!< Return a reference to a sequence in the character matrix
    const DiscreteTaxonData<charType> &getTaxonData(size_t tn) const;             //!< Return a reference to a sequence in the character matrix
    DiscreteTaxonData<charType> &getTaxonData(const std::string &tn);             //!< Return a reference to a sequence in the character matrix
    const DiscreteTaxonData<charType> &getTaxonData(const std::string &tn) const; //!< Return a reference to a sequence in the character matrix

    double maxGcContent(bool excl) const;                   //!< Maximum GC-content of a sequence
    size_t maxInvariableBlockLength(bool excl) const;       //!< Maximum length of a block of invariant sites
    size_t maxVariableBlockLength(bool excl) const;         //!< Maximum length of a block of variant sites
    double meanGcContent(bool excl) const;                  //!< Mean GC-content of all sequence
    double meanGcContentByCodon(size_t n, bool excl) const; //!< Mean GC-content of all sequences by codon position
    double minGcContent(bool excl) const;                   //!< Number of invariant sites
    size_t numInvariableSiteBlocks(bool excl) const;        //!< Number of invariant blocks
    size_t numberTaxaMissingSequence(double p) const;       //!< Number of taxa missing x percent of the sequence
    double varGcContent(bool excl) const;                   //!< Mean GC-content of all sequence
    double varGcContentByCodon(size_t n, bool excl) const;  //!< Mean GC-content of all sequences by codon position

    HomologousDiscreteCharacterData<NaturalNumbersState> *combineCharacters(const HomologousDiscreteCharacterData &d) const; //!< Combine/expand data matrices
    AbstractHomologousDiscreteCharacterData *combineCharacters(const AbstractHomologousDiscreteCharacterData &d) const;      //!< Combine/expand data matrices
    AbstractHomologousDiscreteCharacterData *expandCharacters(size_t n) const;
    AbstractHomologousDiscreteCharacterData *expandStandard(size_t n) const;

    void includeCharacter(size_t i);                                     //!< Include character
    bool isCharacterExcluded(size_t i) const;                            //!< Is the character excluded
    bool isCharacterResolved(size_t txIdx, size_t chIdx) const;          //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
    bool isCharacterResolved(const std::string &tn, size_t chIdx) const; //!< Returns whether the character is fully resolved (e.g., "A" or "1.32") or not (e.g., "AC" or "?")
    void removeExcludedCharacters(void);                                 //!< Remove all the excluded characters
    void restoreCharacter(size_t i);                                     //!< Restore character
    AbstractHomologousDiscreteCharacterData *translateCharacters(const std::string &type) const;

  protected:
    // Utility functions
    bool isCharacterMissingOrAmbiguous(size_t idx) const; //!< Does the character have missing or ambiguous data?

    // Member variables
    std::set<size_t> deletedCharacters; //!< Set of deleted characters
  };

}

#include "CharacterTranslator.h"
#include "DnaState.h"
#include "NaturalNumbersState.h"
#include "ExpandedStandardState.h"
#include "NclReader.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"

#include <cmath>
#include <fstream>
#include <sstream>

/**
 * Default constructor,
 * Does nothing except instanciating the object.
 */
template <class charType>
RevBayesCore::HomologousDiscreteCharacterData<charType>::HomologousDiscreteCharacterData()
{
}

/**
 * Index (const) operator to access a TaxonData object at position i.
 *
 * \param[in]    i    The position of the TaxonData object.
 *
 * \return            The TaxonData object at position i.
 */
template <class charType>
const RevBayesCore::DiscreteTaxonData<charType> &RevBayesCore::HomologousDiscreteCharacterData<charType>::operator[](const size_t i) const
{

  return getTaxonData(i);
}

/**
 * Less than operator.
 */
template <class charType>
bool RevBayesCore::HomologousDiscreteCharacterData<charType>::operator<(const HomologousDiscreteCharacterData<charType> &x) const
{

  return taxa.size() < x.taxa.size();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the object.
 */
template <class charType>
RevBayesCore::HomologousDiscreteCharacterData<charType> *RevBayesCore::HomologousDiscreteCharacterData<charType>::clone(void) const
{

  return new HomologousDiscreteCharacterData<charType>(*this);
}

/**
 * Compute the state frequencies per site.
 *
 * \return       A matrix of character frequencies where each column is a character and each row a taxon.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::computeMultinomialProfileLikelihood(void) const
{

  std::vector<double> pattern_counts;

  // resize the matrices
  size_t num_sequences = this->taxa.size();

  // create a vector with the correct site indices
  // some of the sites may have been excluded
  std::vector<size_t> site_indices = getIncludedSiteIndices();

  size_t num_sites = getNumberOfIncludedCharacters();

  // find the unique site patterns and compute their respective frequencies
  std::map<std::string, size_t> patterns;
  for (size_t site = 0; site < num_sites; ++site)
  {
    // create the site pattern
    std::string pattern = "";
    for (size_t i = 0; i < num_sequences; ++i)
    {

      const DiscreteTaxonData<charType> &seq = this->getTaxonData(i);
      const DiscreteCharacterState &c = seq.getCharacter(site_indices[site]);
      pattern += c.getStringValue();
    }

    // check if we have already seen this site pattern
    std::map<std::string, size_t>::const_iterator index = patterns.find(pattern);
    if (index != patterns.end())
    {

      // we have already seen this pattern
      // increase the frequency counter
      pattern_counts[index->second]++;
    }
    else
    {

      // insert this pattern with the corresponding index in the map
      patterns.insert(std::pair<std::string, size_t>(pattern, pattern_counts.size()));

      // create a new pattern frequency counter for this pattern
      pattern_counts.push_back(1);
    }
  }

  double lnl = 0.0;
  for (size_t i = 0; i < pattern_counts.size(); ++i)
  {
    double c = pattern_counts[i];
    lnl += c * log(c);
  }
  lnl -= double(num_sites) * log(num_sites);

  return lnl;
}

/**
 * Compute the state frequencies per site.
 *
 * \return       A matrix of character frequencies where each column is a character and each row a taxon.
 */
template <class charType>
RevBayesCore::MatrixReal RevBayesCore::HomologousDiscreteCharacterData<charType>::computeStateFrequencies(void) const
{

  charType tmp = this->getTaxonData(0)[0];
  size_t num_states = tmp.getNumberOfStates();
  size_t numSequences = this->taxa.size();
  MatrixReal m(numSequences, num_states);

  double MIN_THRESHOLD = 1E-3;
  for (size_t i = 0; i < numSequences; ++i)
  {
    const DiscreteTaxonData<charType> &seq = this->getTaxonData(i);
    size_t l = seq.getNumberOfCharacters();
    double nonGapSeqLength = 0.0;
    std::vector<double> stateCounts(num_states, MIN_THRESHOLD);
    for (size_t j = 0; j < l; ++j)
    {
      const charType &c = seq[j];

      if (c.isMissingState() == true)
      {
        nonGapSeqLength++;

        for (size_t index = 0; index < num_states; ++index)
        {
          stateCounts[index] += 1.0 / double(num_states);
        }
      }
      else if (c.isGapState() == false)
      {
        ++nonGapSeqLength;

        double numObservedStates = c.getNumberObservedStates();

        for (size_t k = 0; k < num_states; ++k)
        {

          if (c.isStateSet(k) == true)
          {
            // add a uniform probability of having observed each of the ambiguous characters
            stateCounts[k] += 1.0 / numObservedStates;
          }
        }
      }

    } // finished loop over sequence

    // set the observed state frequencies for this sequence into the matrix
    std::vector<double> &observedFreqs = m[i];
    for (size_t j = 0; j < num_states; ++j)
    {
      observedFreqs[j] = stateCounts[j] / (nonGapSeqLength + 20 * MIN_THRESHOLD);
    }
  }

  return m;
}

/**
 * Compute the site frequency spectrum.
 *
 * \return       A vector of occurrences for the site frequency spectrum.
 */
template <class charType>
std::vector<long> RevBayesCore::HomologousDiscreteCharacterData<charType>::computeSiteFrequencySpectrum(bool folded, SFS_AMBIGUITY_TREATMENT ambig_treat) const
{

  charType tmp = this->getTaxonData(0)[0];
  size_t num_states = tmp.getNumberOfStates();
  size_t num_sequences = this->taxa.size();
  std::vector<long> sfs;

  if (folded == true)
  {
    sfs = std::vector<long>(int(num_sequences / 2) + 1, 0);
  }
  else
  {
    sfs = std::vector<long>(num_sequences + 1, 0);
  }

  const DiscreteTaxonData<charType> &tmp_seq = this->getTaxonData(0);
  size_t seq_len = tmp_seq.getNumberOfCharacters();
  for (size_t i = 0; i < seq_len; ++i)
  {

    size_t derived_count = 0;
    size_t this_SNP_count = num_sequences;
    bool skip_this_SNP = false;
    for (size_t j = 0; j < num_sequences; ++j)
    {
      const charType &c = this->getTaxonData(j)[i];

      if (c.isAmbiguous() == true)
      {
        if (ambig_treat == DERIVED)
        {
          ++derived_count;
        }
        else if (ambig_treat == RESCALE)
        {
          // decrease the counter of the SNP;
          --this_SNP_count;
        }
        else if (ambig_treat == SKIP_COLUMN)
        {
          skip_this_SNP = true;
          break;
        }
      }
      else if (c.isGapState() == false)
      {

        for (size_t k = 1; k < num_states; ++k)
        {

          if (c.isStateSet(k) == true)
          {
            ++derived_count;
            break;
          }
        }
      }
      else
      {
        throw RbException("Unexpected character in SFS calculation.");
      }

    } // finished loop over taxa

    if (skip_this_SNP == false)
    {

      //            if ( folded == true )
      //            {
      //                derived_count = RbMath::min(int(derived_count), int(this_SNP_count-derived_count));
      //            }

      if (ambig_treat == RESCALE)
      {
        double rescaled_derived_count = double(derived_count) / double(num_sequences) * double(this_SNP_count);
        derived_count = int(round(rescaled_derived_count));
      }

      if (folded == true)
      {
        derived_count = RbMath::min(int(derived_count), int(num_sequences - derived_count));
      }

      sfs[derived_count]++;
    }

  } // finished loop over sequence

  return sfs;
}

/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template <class charType>
RevBayesCore::AbstractHomologousDiscreteCharacterData *RevBayesCore::HomologousDiscreteCharacterData<charType>::combineCharacters(const AbstractHomologousDiscreteCharacterData &obsd) const
{

  const HomologousDiscreteCharacterData<charType> *rhs = dynamic_cast<const HomologousDiscreteCharacterData<charType> *>(&obsd);
  if (rhs == NULL)
  {
    throw RbException("Cannot combine CharacterData of two different types!!!");
  }

  return combineCharacters(*rhs);
}

/**
 * Add another character data object to this character data object.
 *
 * \param[in]    second    The CharacterData object that should be added.
 */
template <class charType>
RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::NaturalNumbersState> *RevBayesCore::HomologousDiscreteCharacterData<charType>::combineCharacters(const HomologousDiscreteCharacterData<charType> &second) const
{

  HomologousDiscreteCharacterData<charType> *tmp_first = this->clone();
  HomologousDiscreteCharacterData<charType> *tmp_second = second.clone();
  HomologousDiscreteCharacterData<NaturalNumbersState> *combined_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();

  std::string type = "union";
  //    size_t sequence_length = getNumberOfCharacters();

  // check if both have the same number of taxa
  //    if ( taxa.size() != obsd.getNumberOfTaxa() )
  //    {
  //        throw RbException("Cannot combine two character data objects with different number of taxa!");
  //    }

  std::vector<string> to_delete;
  //    std::vector<bool> used = std::vector<bool>(obsd.getNumberOfTaxa(),false);

  for (size_t i = 0; i < second.getNumberOfTaxa(); ++i)
  {
    const std::string &n = second.getTaxonNameWithIndex(i);
    try
    {
      getIndexOfTaxon(n);
    }
    catch (RbException &e)
    {
      if (type == "union")
      {
        tmp_first->addMissingTaxon(n);
      }
      else if (type != "intersection")
      {
        to_delete.push_back(n);
        //                throw RbException() << "Cannot concatenate two character data objects because first character data object has no taxon with name '" << obsd.getTaxonNameWithIndex(i) << "n'!";
      }
    }
  }

  for (size_t i = 0; i < taxa.size(); ++i)
  {
    const std::string &n = taxa[i].getName();
    const DiscreteTaxonData<charType> &taxon = tmp_first->getTaxonData(n);

    try
    {
      // size_t idx = obsd.getIndexOfTaxon( n );
      //            tmp_second->getIndexOfTaxon( n );
      DiscreteTaxonData<NaturalNumbersState> *new_taxon_data = taxon.combineCharacters(tmp_second->getTaxonData(n));
      combined_data->addTaxonData(*new_taxon_data);

      delete new_taxon_data;
    }
    catch (RbException &e)
    {
      if (type == "intersection")
      {
        to_delete.push_back(n);
      }
      else if (type == "union")
      {
        tmp_second->addMissingTaxon(n);
        DiscreteTaxonData<NaturalNumbersState> *new_taxon_data = taxon.combineCharacters(tmp_second->getTaxonData(n));
        combined_data->addTaxonData(*new_taxon_data);

        delete new_taxon_data;
      }
      else
      {
        throw RbException() << "Cannot concatenate two character data objects because second character data object has no taxon with name '" << n << "n'!";
      }
    }
  }

  for (size_t i = 0; i < to_delete.size(); i++)
  {
    combined_data->deleteTaxon(to_delete[i]);
  }

  //    const std::set<size_t> &exclChars = obsd.getExcludedCharacters();
  //    for (std::set<size_t>::const_iterator it = exclChars.begin(); it != exclChars.end(); ++it)
  //    {
  //        combined_data->deletedCharacters.insert( *it + sequence_length );
  //    }

  return combined_data;
}

/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::concatenate(const AbstractCharacterData &obsd, std::string type)
{

  const HomologousDiscreteCharacterData<charType> *rhs = dynamic_cast<const HomologousDiscreteCharacterData<charType> *>(&obsd);
  if (rhs == NULL)
  {
    throw RbException("Adding wrong character data type into CharacterData!!!");
  }

  concatenate(*rhs, type);
}

/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::concatenate(const HomologousCharacterData &obsd, std::string type)
{

  const HomologousDiscreteCharacterData<charType> *rhs = dynamic_cast<const HomologousDiscreteCharacterData<charType> *>(&obsd);
  if (rhs == NULL)
  {
    throw RbException("Adding wrong character data type into CharacterData!!!");
  }

  concatenate(*rhs, type);
}

/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::concatenate(const AbstractHomologousDiscreteCharacterData &obsd, std::string type)
{

  const HomologousDiscreteCharacterData<charType> *rhs = dynamic_cast<const HomologousDiscreteCharacterData<charType> *>(&obsd);
  if (rhs == NULL)
  {
    throw RbException("Adding wrong character data type into CharacterData!!!");
  }

  concatenate(*rhs, type);
}

/**
 * Add another character data object to this character data object.
 *
 * \param[in]    obsd    The CharacterData object that should be added.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::concatenate(const HomologousDiscreteCharacterData<charType> &obsd, std::string type)
{

  size_t sequence_length = getNumberOfCharacters();

  // check if both have the same number of taxa
  if (taxa.size() != obsd.getNumberOfTaxa() && type != "union" && type != "intersection")
  {
    throw RbException("Cannot concatenate two character data objects with different number of taxa!");
  }

  std::vector<string> to_delete;
  //    std::vector<bool> used = std::vector<bool>(obsd.getNumberOfTaxa(),false);

  for (size_t i = 0; i < obsd.getNumberOfTaxa(); ++i)
  {
    const std::string &n = obsd.getTaxonNameWithIndex(i);
    try
    {
      getIndexOfTaxon(n);
    }
    catch (RbException &e)
    {
      if (type == "union")
      {
        addMissingTaxon(n);
      }
      //            else if ( type != "intersection" )
      //            {
      //                toDelete.push_back( n );
      //                throw RbException() << "Cannot concatenate two character data objects because first character data object has no taxon with name '" << obsd.getTaxonNameWithIndex(i) << "n'!";
      //            }
    }
  }

  for (size_t i = 0; i < taxa.size(); ++i)
  {
    const std::string &n = taxa[i].getName();
    DiscreteTaxonData<charType> &taxon = getTaxonData(n);

    try
    {
      // size_t idx = obsd.getIndexOfTaxon( n );
      obsd.getIndexOfTaxon(n);
      taxon.concatenate(obsd.getTaxonData(n));
    }
    catch (RbException &e)
    {
      if (type == "intersection")
      {
        to_delete.push_back(n);
      }
      else if (type == "union")
      {
        AbstractTaxonData *taxon_data = obsd.getTaxonData(0).clone();
        taxon_data->setAllCharactersMissing();
        taxon.concatenate(*taxon_data);
        delete taxon_data;
      }
      else
      {
        throw RbException() << "Cannot concatenate two character data objects because second character data object has no taxon with name '" << n << "n'!";
      }
    }
  }

  for (size_t i = 0; i < to_delete.size(); i++)
  {
    deleteTaxon(to_delete[i]);
  }

  const std::set<size_t> &exclChars = obsd.getExcludedCharacters();
  for (std::set<size_t>::const_iterator it = exclChars.begin(); it != exclChars.end(); ++it)
  {
    deletedCharacters.insert(*it + sequence_length);
  }
}

/**
 * Exclude all characters.
 * We don't actually delete the characters but mark them for exclusion.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::excludeAllCharacters(void)
{

  for (size_t i = 0; i < getTaxonData(0).getNumberOfCharacters(); ++i)
  {
    deletedCharacters.insert(i);
  }
}

/**
 * Exclude a character.
 * We don't actually delete the character but mark it for exclusion.
 *
 * \param[in]    i    The position of the character that will be excluded.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::excludeCharacter(size_t i)
{

  if (i >= getTaxonData(0).getNumberOfCharacters())
  {
    std::stringstream o;
    o << "Only " << getNumberOfCharacters() << " characters in matrix";
    throw RbException(o.str());
  }

  deletedCharacters.insert(i);
}

/**
 * Compute the state frequencies per site.
 *
 * \return       A matrix of character frequencies where each column is a character and each row a taxon.
 */
template <class charType>
RevBayesCore::AbstractHomologousDiscreteCharacterData *RevBayesCore::HomologousDiscreteCharacterData<charType>::expandCharacters(size_t n) const
{

  size_t num_sequences = this->taxa.size();
  HomologousDiscreteCharacterData<NaturalNumbersState> *trans_char_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();

  for (size_t i = 0; i < num_sequences; ++i)
  {
    const DiscreteTaxonData<charType> &seq = this->getTaxonData(i);
    DiscreteTaxonData<NaturalNumbersState> expanded_seq = DiscreteTaxonData<NaturalNumbersState>(seq.getTaxon());

    size_t seq_length = seq.getNumberOfCharacters();
    for (size_t j = 0; j < seq_length; ++j)
    {
      const charType &org_char = seq[j];

      size_t num_states = org_char.getNumberOfStates();
      NaturalNumbersState number_state = NaturalNumbersState(0, int(num_states * n));

      if (org_char.isMissingState() == true)
      {
        number_state.setMissingState(true);
      }
      else if (org_char.isGapState() == true)
      {
        number_state.setGapState(true);
      }
      else
      {

        bool first = true;
        for (size_t k = 0; k < num_states; ++k)
        {

          if (org_char.isStateSet(k) == true)
          {

            // set the initial state
            if (first == true)
            {
              first = false;
              number_state.setStateByIndex(k);
            }
            else
            {
              number_state.addState(int(k));
            }

            // now we set also the expanded states
            for (size_t l = 1; l < n; ++l)
            {
              number_state.addState(int(k) + int(num_states * l));
            }
          }
        }
      }

      expanded_seq.addCharacter(number_state);
    }

    trans_char_data->addTaxonData(expanded_seq);
  }

  return trans_char_data;
}

/**
 * Expand Character states
 *
 * \return       expanded character data matrix
 */
template <class charType>
RevBayesCore::AbstractHomologousDiscreteCharacterData *RevBayesCore::HomologousDiscreteCharacterData<charType>::expandStandard(size_t n) const
{

  size_t num_sequences = this->taxa.size();
  HomologousDiscreteCharacterData<ExpandedStandardState> *trans_char_data = new HomologousDiscreteCharacterData<ExpandedStandardState>();

  for (size_t i = 0; i < num_sequences; ++i)
  {
    const DiscreteTaxonData<charType> &seq = this->getTaxonData(i);
    DiscreteTaxonData<ExpandedStandardState> expanded_seq = DiscreteTaxonData<ExpandedStandardState>(seq.getTaxon());

    size_t seq_length = seq.getNumberOfCharacters();
    for (size_t j = 0; j < seq_length; ++j)
    {
      const charType &org_char = seq[j];

      size_t num_states = org_char.getNumberOfStates();
      ExpandedStandardState number_state = ExpandedStandardState(0, int(num_states * n));

      if (org_char.isMissingState() == true)
      {
        number_state.setMissingState(true);
      }
      else if (org_char.isGapState() == true)
      {
        number_state.setGapState(true);
      }
      else
      {

        bool first = true;
        for (size_t k = 0; k < num_states; ++k)
        {

          if (org_char.isStateSet(k) == true)
          {

            // set the initial state
            if (first == true)
            {
              first = false;
              number_state.setStateByIndex(k);
            }
            else
            {
              number_state.addState(int(k));
            }

            // now we set also the expanded states
            for (size_t l = 1; l < n; ++l)
            {
              number_state.addState(int(k) + int(num_states * l));
            }
          }
        }
      }

      expanded_seq.addCharacter(number_state);
    }

    trans_char_data->addTaxonData(expanded_seq);
  }

  return trans_char_data;
}

/**
 * Get the cn-th character of the tn-th taxon.
 *
 * \param[in]    tn     The index/position of the taxon.
 * \param[in]    cn     The position of the character.
 *
 * \return              The cn-th character of the tn-th taxon.
 */
template <class charType>
const charType &RevBayesCore::HomologousDiscreteCharacterData<charType>::getCharacter(size_t tn, size_t cn) const
{

  if (cn >= getNumberOfCharacters())
  {
    throw RbException("Character index out of range");
  }

  return getTaxonData(tn)[cn];
}

/**
 * Get the data type of the character stored in this object.
 *
 * \return      The data type (e.g. DNA, RNA or Standard).
 */
template <class charType>
std::string RevBayesCore::HomologousDiscreteCharacterData<charType>::getDataType(void) const
{

  std::string dt = "";
  if (taxa.size() > 0)
  {
    const DiscreteTaxonData<charType> &t = getTaxonData(taxa[0].getName());
    if (t.getNumberOfCharacters() > 0)
    {
      dt = t[0].getDataType();
    }
  }

  return dt;
}

/**
 * Get the set of excluded character indices.
 *
 * \return    The excluded character indices.
 */
template <class charType>
const std::set<size_t> &RevBayesCore::HomologousDiscreteCharacterData<charType>::getExcludedCharacters(void) const
{

  return deletedCharacters;
}

/**
 * Get the set of excluded character indices.
 *
 * \return    The excluded character indices.
 */
template <class charType>
std::vector<double> RevBayesCore::HomologousDiscreteCharacterData<charType>::getEmpiricalBaseFrequencies(void) const
{
  size_t n_states = this->getTaxonData(0)[0].getNumberOfStates();
  std::vector<double> ebf = std::vector<double>(n_states, 0.0);
  double total = 0.0;
  size_t nt = this->getNumberOfTaxa();
  for (size_t i = 0; i < nt; ++i)
  {

    const AbstractDiscreteTaxonData &taxon_data = this->getTaxonData(i);
    size_t nc = taxon_data.getNumberOfCharacters();
    for (size_t j = 0; j < nc; ++j)
    {
      const DiscreteCharacterState &o = taxon_data[j];
      if (o.isAmbiguous() == false)
      {
        ++total;
        ++ebf[o.getStateIndex()];
      }
    }
  }

  for (size_t i = 0; i < n_states; ++i)
  {
    ebf[i] /= total;
  }

  return ebf;
}

template <class charType>
std::vector<size_t> RevBayesCore::HomologousDiscreteCharacterData<charType>::getIncludedSiteIndices(void) const
{
  // create a vector with the correct site indices
  // some of the sites may have been excluded
  std::vector<size_t> site_indices;
  size_t site_index = 0;

  size_t num_sites = getNumberOfIncludedCharacters();
  for (size_t i = 0; i < num_sites; ++i)
  {
    while (this->isCharacterExcluded(site_index))
    {
      site_index++;
      if (site_index >= this->getNumberOfCharacters())
      {
        throw RbException("The character matrix cannot set to this variable because it does not have enough included characters.");
      }
    }

    site_indices.push_back(site_index);
    site_index++;
  }

  // test if there were additional sites that we did not use
  while (site_index < this->getNumberOfCharacters())
  {
    if (!this->isCharacterExcluded(site_index))
    {
      throw RbException("The character matrix cannot set to this variable because it has too many included characters.");
    }
    site_index++;
  }

  return site_indices;
}

/**
 * Get the indices of invariant characters in taxon data object.
 * This is regardless of whether the character are included or excluded.
 *
 * \return    The total number of characters
 */
template <class charType>
std::vector<size_t> RevBayesCore::HomologousDiscreteCharacterData<charType>::getInvariantSiteIndices(bool exclude_missing) const
{
  // create a vector with the correct site indices
  // some of the sites may have been excluded
  std::vector<size_t> inv_site_indices;
  size_t nt = this->getNumberOfTaxa();

  const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(0);
  size_t nc = firstTaxonData.getNumberOfCharacters();
  for (size_t j = 0; j < nc; j++)
  {
    const DiscreteCharacterState *a = &firstTaxonData[j];
    size_t k = 1;
    while (exclude_missing == true && a->isAmbiguous() && k < nt)
    {
      const AbstractDiscreteTaxonData &td = this->getTaxonData(k);
      a = &td[j];
      ++k;
    }

    bool invariant = true;
    for (size_t i = 1; i < nt; i++)
    {
      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(i);
      const DiscreteCharacterState &b = secondTaxonData[j];

      if (exclude_missing == false || (b.isAmbiguous() == false && a->isAmbiguous() == false))
      {
        invariant &= (*a == b);
      }
    }

    if (invariant == true)
    {
      inv_site_indices.push_back(j + 1);
    }
  }

  return inv_site_indices;
}

/**
 * Get the number of characters in taxon data object.
 * This i regardless of whether the character are included or excluded.
 * For simplicity we assume that all taxon data objects contain the same number
 * of character and thus we simply return the number from the first taxon data object.
 *
 * \return    The total number of characters
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getNumberOfCharacters(void) const
{

  if (getNumberOfTaxa() > 0)
  {
    return getTaxonData(0).getNumberOfCharacters();
  }

  return 0;
}

/**
 * Get the number of characters in taxon data object.
 * This i regardless of whether the character are included or excluded.
 * For simplicity we assume that all taxon data objects contain the same number
 * of character and thus we simply return the number from the first taxon data object.
 *
 * \return    The total number of characters
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getNumberOfIncludedCharacters(void) const
{

  if (getNumberOfTaxa() > 0)
  {
    return getTaxonData(0).getNumberOfCharacters() - deletedCharacters.size();
  }
  return 0;
}

/**
 * Get the number of observed states for the characters in this object.
 *
 * \return      The number of observed states for the characters.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getMaxObservedStateIndex(void) const
{

  // Get the first character in the matrix
  if (getNumberOfTaxa() == 0)
  {
    return 0;
  }

  RbBitSet observed(getNumberOfStates());

  for (size_t j = 0; j < getNumberOfCharacters(); j++)
  {
    if (isCharacterExcluded(j))
      continue;

    for (size_t i = 0; i < getNumberOfTaxa(); i++)
    {
      const DiscreteTaxonData<charType> &sequence = getTaxonData(i);
      if (sequence[j].isMissingState() == false)
      {
        observed |= sequence[j].getState();
      }
    }
  }

  int max;
  for (max = int(observed.size()) - 1; max >= 0; max--)
  {
    if (observed.test(max))
    {
      break;
    }
  }

  return size_t(max);
}

/**
 * Get the number of states for the characters in this object.
 * We assume that all of the characters in the matrix are of the same
 * type and have the same number of potential states.
 *
 * \return      The number of states for the characters.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getNumberOfStates(void) const
{

  // Get the first character in the matrix
  if (getNumberOfTaxa() == 0)
  {
    return 0;
  }

  const DiscreteTaxonData<charType> &sequence = getTaxonData(0);
  if (sequence.getNumberOfCharacters() == 0)
  {
    return 0;
  }

  return sequence[0].getNumberOfStates();
}

/**
 * Get the number of invariant characters.
 *
 * \return    The number of invariant sites.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getNumberOfInvariantSites(bool exclude_missing) const
{
  std::vector<size_t> inv_site_indices = this->getInvariantSiteIndices(exclude_missing);
  size_t invSites = inv_site_indices.size();
  return invSites;
}

/**
 * Get the set of excluded character indices.
 * We use the fact that the number of segregating sites is num_sites - numInvariantSites.
 *
 * \return    The excluded character indices.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getNumberOfSegregatingSites(bool exclude_missing) const
{
  const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(0);
  size_t nc = firstTaxonData.getNumberOfCharacters();

  return nc - getNumberOfInvariantSites(exclude_missing);
}

/**
 * Get the average pairwise distance between the sequences.
 *
 * \return    The average pairwise distance.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::getAveragePairwiseSequenceDifference(bool exclude_missing) const
{
  double pairwiseDistance = 0.0;
  size_t nt = this->getNumberOfIncludedTaxa();

  for (size_t i = 0; i < (nt - 1); i++)
  {

    const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(i);
    size_t nc = firstTaxonData.getNumberOfCharacters();

    for (size_t j = i + 1; j < nt; j++)
    {

      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(j);
      double pd = 0.0;

      for (size_t k = 0; k < nc; k++)
      {
        const DiscreteCharacterState &a = firstTaxonData[k];
        const DiscreteCharacterState &b = secondTaxonData[k];
        if (exclude_missing == true || (a.isAmbiguous() == false && b.isAmbiguous() == false))
        {
          if (a != b)
          {
            ++pd;
          }
        }
      }

      pairwiseDistance += pd;

    } // end loop over all second taxa

  } // end loop over all first taxa

  pairwiseDistance *= 2.0 / (nt * (nt - 1.0));

  return pairwiseDistance;
}

/**
 * Get the maximum pairwise distance between the sequences.
 *
 * \return    The max pairwise distance.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getMaxPairwiseSequenceDifference(bool exclude_missing) const
{
  size_t max_pd = 0.0;
  size_t nt = this->getNumberOfIncludedTaxa();

  for (size_t i = 0; i < (nt - 1); ++i)
  {

    const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(i);
    size_t nc = firstTaxonData.getNumberOfCharacters();

    for (size_t j = i + 1; j < nt; ++j)
    {

      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(j);
      size_t pd = 0.0;

      for (size_t k = 0; k < nc; ++k)
      {

        const DiscreteCharacterState &a = firstTaxonData[k];
        const DiscreteCharacterState &b = secondTaxonData[k];
        if (exclude_missing == true || (a.isAmbiguous() == false && b.isAmbiguous() == false))
        {
          if (a != b)
          {
            ++pd;
          }
        }
      }

      if (max_pd < pd)
      {
        max_pd = pd;
      }

    } // end loop over all second taxa

  } // end loop over all first taxa

  return max_pd;
}

/**
 * Get the minimum pairwise distance between the sequences.
 *
 * \return    The min pairwise distance.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::getMinPairwiseSequenceDifference(bool include_missing) const
{
  size_t min_pd = RbConstants::Size_t::max;
  size_t nt = this->getNumberOfIncludedTaxa();

  for (size_t i = 0; i < (nt - 1); i++)
  {

    const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(i);
    size_t nc = firstTaxonData.getNumberOfCharacters();

    for (size_t j = i + 1; j < nt; j++)
    {

      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(j);
      size_t pd = 0.0;

      for (size_t k = 0; k < nc; k++)
      {
        const DiscreteCharacterState &a = firstTaxonData[k];
        const DiscreteCharacterState &b = secondTaxonData[k];
        if (include_missing == true || (a.isAmbiguous() == false && b.isAmbiguous() == false))
        {
          if (a != b)
          {
            ++pd;
          }
        }
      }

      if (min_pd > pd)
      {
        min_pd = pd;
      }

    } // end loop over all second taxa

  } // end loop over all first taxa

  return min_pd;
}

/**
 * Get the minimum pairwise distance between the sequences.
 *
 * \return    The min pairwise distance.
 */
template <class charType>
RevBayesCore::DistanceMatrix RevBayesCore::HomologousDiscreteCharacterData<charType>::getPairwiseSequenceDifference(bool include_missing) const
{
  size_t n_taxa_included = this->getNumberOfIncludedTaxa();
  size_t n_taxa_total = this->getNumberOfTaxa();
  MatrixReal distances = MatrixReal(n_taxa_included);

  size_t index_included_i = 0;
  size_t index_total_i = 0;
  for (; index_total_i < (n_taxa_total - 1); ++index_total_i)
  {
    // skip over excluded taxa
    if (isTaxonExcluded(index_total_i) == true)
      continue;

    const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(index_total_i);
    size_t nc = firstTaxonData.getNumberOfCharacters();

    size_t index_total_j = index_total_i + 1;
    size_t index_included_j = index_included_i + 1;
    for (; index_total_j < n_taxa_total; ++index_total_j)
    {

      if (isTaxonExcluded(index_total_j) == true)
        continue;

      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(index_total_j);
      size_t pd = 0.0;

      for (size_t k = 0; k < nc; k++)
      {
        const DiscreteCharacterState &a = firstTaxonData[k];
        const DiscreteCharacterState &b = secondTaxonData[k];
        if (include_missing == true || (a.isAmbiguous() == false && b.isAmbiguous() == false))
        {
          if (a != b)
          {
            ++pd;
          }
        }
      }

      distances[index_included_i][index_included_j] = pd;
      distances[index_included_j][index_included_i] = pd;

      ++index_included_j;

    } // end loop over all second taxa

    distances[index_included_i][index_included_i] = 0;

    ++index_included_i;

  } // end loop over all first taxa

  return DistanceMatrix(distances, getIncludedTaxa());
}

/**
 * Get the taxon data object with index tn.
 *
 * \return     A const reference to the taxon data object at position tn.
 */
template <class charType>
const RevBayesCore::DiscreteTaxonData<charType> &RevBayesCore::HomologousDiscreteCharacterData<charType>::getTaxonData(size_t tn) const
{

  if (tn >= getNumberOfTaxa())
  {
    throw RbException("Taxon index out of range");
  }

  const std::string &name = taxa[tn].getName();
  const typename std::map<std::string, AbstractTaxonData *>::const_iterator &i = taxonMap.find(name);

  if (i != taxonMap.end())
  {
    return static_cast<const DiscreteTaxonData<charType> &>(*i->second);
  }
  else
  {
    throw RbException() << "Cannot find taxon '" << name << "' in the HomologousDiscreteCharacterData matrix.";
  }
}

/**
 * Get the taxon data object at position tn.
 *
 * \return     A non-const reference to the taxon data object at position tn.
 */
template <class charType>
RevBayesCore::DiscreteTaxonData<charType> &RevBayesCore::HomologousDiscreteCharacterData<charType>::getTaxonData(size_t tn)
{

  if (tn >= getNumberOfTaxa())
  {
    throw RbException("Taxon index out of range");
  }

  const std::string &name = taxa[tn].getName();
  const typename std::map<std::string, AbstractTaxonData *>::iterator &i = taxonMap.find(name);

  if (i != taxonMap.end())
  {
    return static_cast<DiscreteTaxonData<charType> &>(*i->second);
  }
  else
  {
    throw RbException() << "Cannot find taxon '" << name << "' in the HomologousDiscreteCharacterData matrix.";
  }
}

/**
 * Get the taxon data object with name tn.
 *
 * \return     A non-const reference to the taxon data object with name tn.
 */
template <class charType>
const RevBayesCore::DiscreteTaxonData<charType> &RevBayesCore::HomologousDiscreteCharacterData<charType>::getTaxonData(const std::string &tn) const
{

  if (tn == "")
  {
    throw RbException("Ambiguous taxon name.");
  }

  const typename std::map<std::string, AbstractTaxonData *>::const_iterator &i = taxonMap.find(tn);

  if (i != taxonMap.end())
  {
    return static_cast<const DiscreteTaxonData<charType> &>(*i->second);
  }
  else
  {
    throw RbException() << "Cannot find taxon '" << tn << "' in the HomologousDiscreteCharacterData matrix.";
  }
}

/**
 * Get the taxon data object with name tn.
 *
 * \return     A const reference to the taxon data object with name tn.
 */
template <class charType>
RevBayesCore::DiscreteTaxonData<charType> &RevBayesCore::HomologousDiscreteCharacterData<charType>::getTaxonData(const std::string &tn)
{

  if (tn == "")
  {
    throw RbException("Ambiguous taxon name.");
  }

  const typename std::map<std::string, AbstractTaxonData *>::iterator &i = taxonMap.find(tn);

  if (i != taxonMap.end())
  {
    return static_cast<DiscreteTaxonData<charType> &>(*i->second);
  }
  else
  {
    throw RbException() << "Cannot find taxon '" << tn << "' in the HomologousDiscreteCharacterData matrix.";
  }
}

/**
 * Include a character.
 * Since we didn't actually deleted the character but marked it for exclusion
 * we can now simply remove the flag.
 *
 * \param[in]    i    The position of the character that will be included.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::includeCharacter(size_t i)
{

  if (i >= getTaxonData(0).getNumberOfCharacters())
  {
    std::stringstream o;
    o << "Only " << getNumberOfCharacters() << " characters in matrix";
    throw RbException(o.str());
  }

  deletedCharacters.erase(i);
}

/**
 * Initialize this object from a file
 *
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::initFromFile(const RevBayesCore::path &dir, const std::string &fn)
{
  auto filepath = dir / (fn + ".nex");

  // get an instance of the NCL reader
  NclReader reader = NclReader();

  std::string myFileType = "nexus";
  std::string dType = this->getDataType();

  std::string suffix = "|" + dType;
  suffix += "|noninterleaved";
  myFileType += suffix;

  // read the content of the file now
  std::vector<AbstractCharacterData *> m_i = reader.readMatrices(filepath, myFileType);

  if (m_i.size() < 1)
  {
    const std::set<std::string> &warnings = reader.getWarnings();
    for (std::set<std::string>::iterator it = warnings.begin(); it != warnings.end(); ++it)
    {
      std::cerr << "NCL-Warning:\t\t" << *it << std::endl;
    }
    throw RbException() << "Could not read character data matrix from file " << filepath << ".";
  }

  HomologousDiscreteCharacterData<charType> *coreM = static_cast<HomologousDiscreteCharacterData<charType> *>(m_i[0]);

  *this = *coreM;

  delete coreM;
}

/**
 * Initialize this object from a file
 *
 * \param[in]   s    The string value.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::initFromString(const std::string &s)
{
  throw RbException("Cannot initialize a discrete character data matrix from a string.");
}

///**
// * Is this character pattern constant at site idx?
// *
// * \param[in]   idx    The site at which we want to know if it is constant?
// */
// template<class charType>
// bool RevBayesCore::HomologousDiscreteCharacterData<charType>::isCharacterConstant(size_t idx) const
//{
//
//    const CharacterState* f = NULL;
//    for ( size_t i=0; i<getNumberOfTaxa(); ++i )
//    {
//        if ( isTaxonExcluded(i) == false )
//        {
//            if ( f == NULL )
//            {
//                f = &getCharacter( i, idx );
//            }
//            else
//            {
//                const CharacterState* s = &getCharacter( i , idx );
//                if ( (*f) != (*s) )
//                {
//                    return false;
//                }
//
//            }
//
//        }
//
//    }
//
//    return true;
//}

/**
 * Is the character excluded?
 *
 * \param[in]    i   The position of the character.
 */
template <class charType>
bool RevBayesCore::HomologousDiscreteCharacterData<charType>::isCharacterExcluded(size_t i) const
{

  std::set<size_t>::const_iterator it = deletedCharacters.find(i);
  if (it != deletedCharacters.end())
  {
    return true;
  }

  return false;
}

template <class charType>
bool RevBayesCore::HomologousDiscreteCharacterData<charType>::isCharacterResolved(size_t txIdx, size_t chIdx) const
{

  const RevBayesCore::DiscreteTaxonData<charType> &td = getTaxonData(txIdx);
  return td.isCharacterResolved(chIdx);
}

template <class charType>
bool RevBayesCore::HomologousDiscreteCharacterData<charType>::isCharacterResolved(const std::string &tn, size_t chIdx) const
{

  const RevBayesCore::DiscreteTaxonData<charType> &td = getTaxonData(tn);
  return td.isCharacterResolved(chIdx);
}

/**
 * Does the character have missing or ambiguous characters?
 *
 * \param[in]    idx    The position of the character in question.
 */
template <class charType>
bool RevBayesCore::HomologousDiscreteCharacterData<charType>::isCharacterMissingOrAmbiguous(size_t idx) const
{

  for (size_t i = 0; i < getNumberOfTaxa(); i++)
  {
    if (isTaxonExcluded(i) == false)
    {
      const CharacterState &c = getCharacter(i, idx);
      if (c.isGapState() == true || c.isAmbiguous() == true)
        return true;
    }
  }

  return false;
}

/**
 * Compute the maximum gc content for any taxon
 *
 * \return    The max GC content.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::maxGcContent(bool exclude_missing) const
{
  double max_gc = 0;
  size_t nt = this->getNumberOfTaxa();

  DnaState G = DnaState("G");
  DnaState C = DnaState("C");

  for (size_t i = 0; i < nt; ++i)
  {
    int num_gc = 0;
    const AbstractDiscreteTaxonData &taxonData = this->getTaxonData(i);
    size_t nc = taxonData.getNumberOfCharacters();
    size_t n_char_this_seq = 0;
    for (size_t j = 0; j < nc; ++j)
    {
      const DiscreteCharacterState &b = taxonData[j];
      if (exclude_missing == false || b.isAmbiguous() == false)
      {
        ++n_char_this_seq;

        if (b == G || b == C)
        {
          ++num_gc;
        }
      }
    }

    double this_gc_proportion = double(num_gc) / double(n_char_this_seq);
    if (max_gc < this_gc_proportion)
    {
      max_gc = this_gc_proportion;
    }
  }

  return max_gc;
}

/**
 * Compute the maximum length of an invariable block.
 *
 * \return    The max block length.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::maxInvariableBlockLength(bool exclude_missing) const
{
  size_t max_length = 0;
  size_t nt = this->getNumberOfTaxa();

  const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(0);
  size_t nc = firstTaxonData.getNumberOfCharacters();

  size_t this_block_length = 0;
  for (size_t j = 0; j < nc; j++)
  {
    const DiscreteCharacterState *a = &firstTaxonData[j];
    size_t k = 1;
    while (exclude_missing == true && a->isAmbiguous() == true && k < nt)
    {
      const AbstractDiscreteTaxonData &td = this->getTaxonData(k);
      a = &td[j];
      ++k;
    }
    bool invariant = true;
    for (size_t i = 1; i < nt; i++)
    {
      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(i);
      const DiscreteCharacterState &b = secondTaxonData[j];

      if (exclude_missing == false || (a->isAmbiguous() == false && b.isAmbiguous() == false))
      {
        invariant &= (*a == b);
      }
    }

    if (invariant == true)
    {
      ++this_block_length;
      if (this_block_length > max_length)
      {
        max_length = this_block_length;
      }
    }
    else
    {
      if (this_block_length > max_length)
      {
        max_length = this_block_length;
      }
      this_block_length = 0;
    }
  }

  return max_length;
}

/**
 * Compute the maximum length of a variable block.
 *
 * \return    The max block length.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::maxVariableBlockLength(bool exclude_missing) const
{
  size_t max_length = 0;
  size_t nt = this->getNumberOfTaxa();

  const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(0);
  size_t nc = firstTaxonData.getNumberOfCharacters();

  size_t this_block_length = 0;
  for (size_t j = 0; j < nc; j++)
  {
    const DiscreteCharacterState *a = &firstTaxonData[j];
    size_t k = 1;
    while (exclude_missing == true && a->isAmbiguous() && k < nt)
    {
      const AbstractDiscreteTaxonData &td = this->getTaxonData(k);
      a = &td[j];
      ++k;
    }
    bool invariant = true;
    for (size_t i = 1; i < nt; i++)
    {
      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(i);
      const DiscreteCharacterState &b = secondTaxonData[j];

      if (exclude_missing == false || (a->isAmbiguous() == false && b.isAmbiguous() == false))
      {
        invariant &= (*a == b);
      }
    }

    if (invariant == false)
    {
      ++this_block_length;
      if (this_block_length > max_length)
      {
        max_length = this_block_length;
      }
    }
    else
    {
      if (this_block_length > max_length)
      {
        max_length = this_block_length;
      }
      this_block_length = 0;
    }
  }

  return max_length;
}

/**
 * Compute the maximum gc content for any taxon
 *
 * \return    The max GC content.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::meanGcContent(bool exclude_ambiguous) const
{
  double mean_gc = 0;
  size_t nt = this->getNumberOfTaxa();

  DnaState G = DnaState("G");
  DnaState C = DnaState("C");

  for (size_t i = 0; i < nt; ++i)
  {
    int num_gc = 0;
    const AbstractDiscreteTaxonData &taxonData = this->getTaxonData(i);
    size_t nc = taxonData.getNumberOfCharacters();
    size_t n_char_this_seq = 0;
    for (size_t j = 0; j < nc; ++j)
    {
      const DiscreteCharacterState &b = taxonData[j];
      if (exclude_ambiguous == false || b.isAmbiguous() == false)
      {
        ++n_char_this_seq;

        if (b == G || b == C)
        {
          ++num_gc;
        }
      }
    }

    double this_gc_proportion = double(num_gc) / double(n_char_this_seq);
    mean_gc += this_gc_proportion;
  }

  return mean_gc / double(nt);
}

/**
 * Compute the maximum gc content for any taxon
 *
 * \return    The max GC content.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::meanGcContentByCodon(size_t n, bool exclude_ambiguous) const
{
  assert(n >= 1 && n <= 3);

  double mean_gc = 0;
  size_t nt = this->getNumberOfTaxa();

  DnaState G = DnaState("G");
  DnaState C = DnaState("C");

  for (size_t i = 0; i < nt; ++i)
  {
    int num_gc = 0;
    const AbstractDiscreteTaxonData &taxonData = this->getTaxonData(i);
    size_t nc = taxonData.getNumberOfCharacters();
    size_t n_char_this_seq = 0;
    for (size_t j = n - 1; j < nc; j += 3)
    {
      const DiscreteCharacterState &b = taxonData[j];
      if (exclude_ambiguous == false || b.isAmbiguous() == false)
      {
        ++n_char_this_seq;

        if (b == G || b == C)
        {
          ++num_gc;
        }
      }
    }

    double this_gc_proportion = double(num_gc) / double(n_char_this_seq);
    mean_gc += this_gc_proportion;
  }

  return mean_gc / double(nt);
}

/**
 * Compute the minimum gc content for any taxon
 *
 * \return    The min GC content.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::minGcContent(bool exclude_ambiguous) const
{
  double min_gc = 1;
  size_t nt = this->getNumberOfTaxa();

  DnaState G = DnaState("G");
  DnaState C = DnaState("C");

  for (size_t i = 0; i < nt; i++)
  {
    int num_gc = 0;
    const AbstractDiscreteTaxonData &taxonData = this->getTaxonData(i);
    size_t nc = taxonData.getNumberOfCharacters();
    size_t n_char_this_seq = 0;
    for (size_t j = 0; j < nc; ++j)
    {
      const DiscreteCharacterState &b = taxonData[j];
      if (exclude_ambiguous == false || b.isAmbiguous() == false)
      {
        ++n_char_this_seq;

        if (b == G || b == C)
        {
          ++num_gc;
        }
      }
    }

    double this_gc_proportion = double(num_gc) / double(n_char_this_seq);
    if (min_gc > this_gc_proportion)
    {
      min_gc = this_gc_proportion;
    }
  }

  return min_gc;
}

/**
 * Compute the maximum length of an invariable block.
 *
 * \return    The max block length.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::numInvariableSiteBlocks(bool exclude_missing) const
{
  size_t num_blocks = 0;
  size_t nt = this->getNumberOfTaxa();

  const AbstractDiscreteTaxonData &firstTaxonData = this->getTaxonData(0);
  size_t nc = firstTaxonData.getNumberOfCharacters();
  bool previous_site_was_invariable = false;
  for (size_t j = 0; j < nc; j++)
  {
    const DiscreteCharacterState *a = &firstTaxonData[j];
    size_t k = 1;
    while (exclude_missing == true && a->isAmbiguous() && k < nt)
    {
      const AbstractDiscreteTaxonData &td = this->getTaxonData(k);
      a = &td[j];
      ++k;
    }
    bool invariant = true;
    for (size_t i = 1; i < nt; i++)
    {
      const AbstractDiscreteTaxonData &secondTaxonData = this->getTaxonData(i);
      const DiscreteCharacterState &b = secondTaxonData[j];

      if (exclude_missing == false || (a->isAmbiguous() == false && b.isAmbiguous() == false))
      {
        invariant &= (*a == b);
      }
    }

    if (previous_site_was_invariable == false && invariant == true)
    {
      ++num_blocks;
    }

    previous_site_was_invariable = invariant;
  }

  return num_blocks;
}

/**
 * Compute the maximum length of an invariable block.
 *
 * \return    The max block length.
 */
template <class charType>
size_t RevBayesCore::HomologousDiscreteCharacterData<charType>::numberTaxaMissingSequence(double x) const
{
  size_t num_taxa_missing = 0;
  size_t num_taxa_total = this->getNumberOfTaxa();

  for (size_t i = 0; i < num_taxa_total; ++i)
  {
    const AbstractDiscreteTaxonData &taxon = this->getTaxonData(i);
    const std::string &n = taxon.getTaxonName();
    double p = getPercentageMissing(n);

    if (p < x)
    {
      ++num_taxa_missing;
    }
  }

  return num_taxa_missing;
}

/**
 * Remove all the excluded character.
 *
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::removeExcludedCharacters(void)
{

  for (typename std::map<std::string, AbstractTaxonData *>::iterator it = taxonMap.begin(); it != taxonMap.end(); ++it)
  {
    it->second->removeCharacters(deletedCharacters);
  }

  deletedCharacters.clear();
}

/**
 * Restore a character. We simply do not mark the character as excluded anymore.
 *
 * \param[in]    i    The position of the character to restore.
 */
template <class charType>
void RevBayesCore::HomologousDiscreteCharacterData<charType>::restoreCharacter(size_t i)
{

  if (i >= getNumberOfCharacters())
    throw RbException("Character index out of range");

  deletedCharacters.erase(i);
}

/**
 * Compute the state frequencies per site.
 *
 * \return       A matrix of character frequencies where each column is a character and each row a taxon.
 */
template <class charType>
RevBayesCore::AbstractHomologousDiscreteCharacterData *RevBayesCore::HomologousDiscreteCharacterData<charType>::translateCharacters(const std::string &type) const
{

  size_t num_sequences = this->taxa.size();
  HomologousDiscreteCharacterData<DnaState> *trans_char_data = new HomologousDiscreteCharacterData<DnaState>();

  for (size_t i = 0; i < num_sequences; ++i)
  {
    const DiscreteTaxonData<charType> &seq = this->getTaxonData(i);
    //        DiscreteTaxonData<DnaState>& trans_seq = seq.convertToDnaCharacters();
    //
    AbstractDiscreteTaxonData *trans_seq = CharacterTranslator::translateCharacters(seq, type);
    trans_char_data->addTaxonData(*trans_seq);

    delete trans_seq;
  }

  return trans_char_data;
}

/**
 * Compute the maximum gc content for any taxon
 *
 * \return    The max GC content.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::varGcContent(bool exclude_ambiguous) const
{

  size_t nt = this->getNumberOfTaxa();

  double mean_gc = 0;
  std::vector<double> gc_content = std::vector<double>(nt, 0);

  DnaState G = DnaState("G");
  DnaState C = DnaState("C");

  for (size_t i = 0; i < nt; ++i)
  {
    int num_gc = 0;
    const AbstractDiscreteTaxonData &taxonData = this->getTaxonData(i);
    size_t nc = taxonData.getNumberOfCharacters();
    size_t n_char_this_seq = 0;
    for (size_t j = 0; j < nc; ++j)
    {
      const DiscreteCharacterState &b = taxonData[j];
      if (exclude_ambiguous == false || b.isAmbiguous() == false)
      {
        ++n_char_this_seq;

        if (b == G || b == C)
        {
          ++num_gc;
        }
      }
    }

    double this_gc_proportion = double(num_gc) / double(n_char_this_seq);
    gc_content[i] = this_gc_proportion;
    mean_gc += this_gc_proportion;
  }

  // normalize the mean;
  mean_gc /= double(nt);

  // now compute the variance
  double var_gc = 0;
  for (size_t i = 0; i < nt; ++i)
  {
    double diff = gc_content[i] - mean_gc;
    var_gc += diff * diff;
  }
  return var_gc / double(nt - 1);
}

/**
 * Compute the maximum gc content for any taxon
 *
 * \return    The max GC content.
 */
template <class charType>
double RevBayesCore::HomologousDiscreteCharacterData<charType>::varGcContentByCodon(size_t n, bool exclude_ambiguous) const
{
  assert(n >= 1 && n <= 3);

  size_t nt = this->getNumberOfTaxa();

  double mean_gc = 0;
  std::vector<double> gc_content = std::vector<double>(nt, 0);

  DnaState G = DnaState("G");
  DnaState C = DnaState("C");

  for (size_t i = 0; i < nt; ++i)
  {
    int num_gc = 0;
    const AbstractDiscreteTaxonData &taxonData = this->getTaxonData(i);
    size_t nc = taxonData.getNumberOfCharacters();
    size_t n_char_this_seq = 0;
    for (size_t j = n - 1; j < nc; j += 3)
    {
      const DiscreteCharacterState &b = taxonData[j];
      if (exclude_ambiguous == false || b.isAmbiguous() == false)
      {
        ++n_char_this_seq;

        if (b == G || b == C)
        {
          ++num_gc;
        }
      }
    }

    double this_gc_proportion = double(num_gc) / double(n_char_this_seq);
    gc_content[i] = this_gc_proportion;
    mean_gc += this_gc_proportion;
  }

  // normalize the mean;
  mean_gc /= double(nt);

  // now compute the variance
  double var_gc = 0;
  for (size_t i = 0; i < nt; ++i)
  {
    double diff = gc_content[i] - mean_gc;
    var_gc += diff * diff;
  }
  return var_gc / double(nt - 1);
}

#endif
