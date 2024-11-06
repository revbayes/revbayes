#ifndef JointConditionalSiteMixtureAllocationMonitor_H
#define JointConditionalSiteMixtureAllocationMonitor_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "VariableMonitor.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"
#include "DnaState.h"

#include <fstream>
#include <cstdint>
#include <string>
#include <vector>
#include <typeinfo>

namespace RevBayesCore {
    
    /**
     * @brief Declaration and implementation of the SiteMixtureAllocationMonitor class.
     *
     * @file
	 * Declaration and implementation of the SiteMixtureAllocationMonitor class which
	 * monitors the allocation of sites to mixture categories and prints their value into a file.
     *
     * @copyright Copyright 2015-
     * @author The RevBayes Development Core Team (Michael R. May)
     * @version 1.0
     *
     */
    template<class characterType>
    class SiteMixtureAllocationMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
		SiteMixtureAllocationMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, const std::string &del);                                  //!< Constructor
		SiteMixtureAllocationMonitor(StochasticNode<Tree>* ch, std::uint64_t g, const std::string &fname, const std::string &del, bool wt, bool wss);
        virtual ~SiteMixtureAllocationMonitor(void);
        
        SiteMixtureAllocationMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // functions you may want to overwrite
        virtual void                                    monitorVariables(std::uint64_t gen);                                //!< Monitor at generation gen
        virtual void                                    printFileHeader(void);                                              //!< Print header
        virtual void                                    swapNode(DagNode *oldN, DagNode *newN);

    private:
        
        // members
		StochasticNode<AbstractHomologousDiscreteCharacterData>* ctmc;

		// flags
		bool has_site_rate_mixture;
		bool has_site_matrix_mixture;

    };
    
}


#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "RbFileManager.h"
#include "PhyloCTMCSiteHomogeneous.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"

using namespace RevBayesCore;


/* Constructor for CTMC */
template<class characterType>
SiteMixtureAllocationMonitor<characterType>::SiteMixtureAllocationMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, const std::string &del) :
    VariableMonitor(ch, g, fname, del, false, false, false),
    ctmc( ch ),
	has_site_rate_mixture(false),
	has_site_matrix_mixture(false)
{
    addVariable( ctmc );

    AbstractPhyloCTMCSiteHomogeneous<characterType> *dist_ctmc = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );

    // determine if the CTMC is using a site-rate mixture
    has_site_rate_mixture = dist_ctmc->hasSiteRateMixture();

    // determine if the CTMC is using a site-matrix mixture
    has_site_matrix_mixture = dist_ctmc->hasSiteMatrixMixture();

}


/**
 * Destructor.
 */
template<class characterType>
SiteMixtureAllocationMonitor<characterType>::~SiteMixtureAllocationMonitor()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
template<class characterType>
SiteMixtureAllocationMonitor<characterType>* SiteMixtureAllocationMonitor<characterType>::clone(void) const
{
    
    return new SiteMixtureAllocationMonitor<characterType>(*this);
}



/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
template<class characterType>
void SiteMixtureAllocationMonitor<characterType>::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
	// make sure the CTMC has been already sampled mixture components
	AbstractPhyloCTMCSiteHomogeneous<characterType> *dist_ctmc = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
	size_t num_sites = dist_ctmc->getValue().getNumberOfIncludedCharacters();
	dist_ctmc->drawSiteMixtureAllocations();

	// now get the sampled mixture components
	std::vector<size_t> rate_components(num_sites);
	std::vector<size_t> matrix_components(num_sites);

	size_t rate_component, matrix_component;
	for(size_t i = 0; i < num_sites; ++i)
	{
		dist_ctmc->getSampledMixtureComponents(i, rate_component, matrix_component);
		rate_components[i]   = rate_component;
		matrix_components[i] = matrix_component;
	}

    // determine if the CTMC is using a site-rate mixture
    has_site_rate_mixture = dist_ctmc->hasSiteRateMixture();

    // determine if the CTMC is using a site-matrix mixture
    has_site_matrix_mixture = dist_ctmc->hasSiteMatrixMixture();

	// print the rate components to file
	if ( has_site_rate_mixture )
	{
		for(size_t i = 0; i < num_sites; ++i)
		{
			// add a separator before every new element
			out_stream << separator;
			out_stream << rate_components[i] + 1;
		}
	}

	// print the matrix components to file
	if ( has_site_matrix_mixture )
	{
		for(size_t i = 0; i < num_sites; ++i)
		{
			// add a separator before every new element
			out_stream << separator;
			out_stream << matrix_components[i] + 1;
		}
	}

}


/**
 * Print header for monitored values
 */
template<class characterType>
void SiteMixtureAllocationMonitor<characterType>::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;

        // get the number of sites
	AbstractPhyloCTMCSiteHomogeneous<characterType> *dist_ctmc = static_cast<AbstractPhyloCTMCSiteHomogeneous<characterType>* >( &ctmc->getDistribution() );
	const TypedDagNode<Tree> *tree = dist_ctmc->getTree();
	size_t num_sites = dist_ctmc->getValue().getNumberOfCharacters();

    // determine if the CTMC is using a site-rate mixture
    has_site_rate_mixture = dist_ctmc->hasSiteRateMixture();

    // determine if the CTMC is using a site-matrix mixture
    has_site_matrix_mixture = dist_ctmc->hasSiteMatrixMixture();

	// print header for site-rates
	if ( has_site_rate_mixture )
	{
		for(size_t i = 0; i < num_sites; ++i)
		{
			out_stream << separator;
            out_stream << "site_rate[" << i + 1 << "]";
		}
	}

	// print header for site-matrices
	if ( has_site_matrix_mixture )
	{
		for(size_t i = 0; i < num_sites; ++i)
		{
			out_stream << separator;
            out_stream << "site_matrix[" << i + 1 << "]";
		}
	}

}


template<class characterType>
void SiteMixtureAllocationMonitor<characterType>::swapNode(DagNode *oldN, DagNode* newN)
{
	if ( oldN == ctmc )
	{
		ctmc = static_cast< StochasticNode<AbstractHomologousDiscreteCharacterData> *>( newN );
	}
    
    Monitor::swapNode( oldN, newN );
}

#endif 
