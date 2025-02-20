#ifndef CTMCProcess_H
#define CTMCProcess_H

#include "Simplex.h"
#include "RateGenerator.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     */
    template<class charType>
    class CTMCProcess : public TypedDistribution< AbstractDiscreteTaxonData > {

    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        CTMCProcess(size_t nc, size_t ns );
        virtual                                                            ~CTMCProcess(void);                                                     //!< Virtual destructor

        // public member functions
        // pure virtual
        virtual CTMCProcess*                                                clone(void) const;                                                                      //!< Create an independent clone

        // non-virtual
        virtual double                                                      computeLnProbability(void);
        void                                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixReal &rv) const;     //!< Map the member methods to internal function calls
        virtual void                                                        redrawValue(void);
        void                                                                reInitialized(void);

        void                                                                setProcessTime(const TypedDagNode< double > *t);
        void                                                                setRateMatrix(const TypedDagNode< RateGenerator > *rm);
        void                                                                setRateMatrix(const TypedDagNode< RbVector< RateGenerator > > *rm);
        void                                                                setRootFrequencies(const TypedDagNode< Simplex > *f);
        void                                                                setSiteMatricesProbs(const TypedDagNode< Simplex > *rp);
        void                                                                setSiteRates(const TypedDagNode< RbVector< double > > *r);
        void                                                                setSiteRatesProbs(const TypedDagNode< Simplex > *rp);
        void                                                                setValue(charType *v, bool f=false);                         //!< Set the current value, e.g. attach an observation (clamp)


    protected:

        // helper method for this and derived classes
        virtual void                                                        compress(void);
        virtual void                                                        computeSiteLikelihoods( std::vector< double > &rv ) const;
        virtual void                                                        computeSiteLikelihoodsPerSiteMatrix( std::vector< std::vector< double > > &rv ) const;
        virtual void                                                        computeSiteLikelihoodsPerSiteRate( std::vector< std::vector< double > > &rv ) const;
        virtual void                                                        computeSiteLikelihoodsPerSiteRateAndMatrix( std::vector< std::vector< std::vector< double > > > &rv ) const;
        virtual std::vector<size_t>                                         getIncludedSiteIndices();
        virtual std::vector<double>                                         getMixtureProbs( void ) const;
        virtual std::vector<double>                                         getRootFrequencies( size_t mixture = 0 ) const;
        virtual void                                                        getRootFrequencies( std::vector<std::vector<double> >& ) const;
        virtual void                                                        setActivePIDSpecialized(size_t i, size_t n);                                                 //!< Set the number of processes for this distribution.
        virtual void                                                        updateTransitionProbabilities( void ) const;


        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                             //!< Swap a parameter


        // virtual methods that may be overwritten, but then the derived class should call this methods


        // members
//        double                                                              lnProb;
//        double                                                              storedLnProb;
        size_t                                                              num_sites;
        const size_t                                                        num_chars;
        size_t                                                              num_site_rates;
        size_t                                                              num_site_matrices;
        mutable std::vector< std::vector<TransitionProbabilityMatrix> >     transition_prob_matrices;


        // the data
        std::vector<RbBitSet>                                               ambiguous_char_vector;
        std::vector<unsigned long>                                          char_vector;
        std::vector<bool>                                                   gap_vector;
        std::vector<size_t>                                                 pattern_counts;
        size_t                                                              num_patterns;
        bool                                                                compressed;
        std::vector<size_t>                                                 site_pattern;    // an array that keeps track of which pattern is used for each site


        // flags
        bool                                                                using_ambiguous_characters;
        bool                                                                using_weighted_characters;

        // members
        const TypedDagNode< RateGenerator >*                                rate_matrix;
        const TypedDagNode< RbVector< RateGenerator > >*                    site_rate_matrices;
        const TypedDagNode< Simplex >*                                      root_frequencies;
        const TypedDagNode< RbVector< double > >*                           site_rates;
        const TypedDagNode< Simplex >*                                      site_matrix_probs;
        const TypedDagNode< Simplex >*                                      site_rates_probs;
        const TypedDagNode< double >*                                       process_time;


        // flags specifying which model variants we use
        bool                                                                rate_variation_across_sites;

        // MPI variables
        size_t                                                              pattern_block_start;
        size_t                                                              pattern_block_end;
        size_t                                                              pattern_block_size;

        charType                                                            template_state;                                 //!< Template state used for ancestral state estimation. This makes sure that the state labels are preserved.


    private:


    };

}


#include "DiscreteCharacterState.h"
#include "DistributionExponential.h"
#include "HomologousDiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMatrix_JC.h"
#include "StochasticNode.h"

#include <cmath>

#ifdef RB_MPI
#include <mpi.h>
#endif


template<class charType>
RevBayesCore::CTMCProcess<charType>::CTMCProcess(size_t nc, size_t ns) :
    TypedDistribution< AbstractDiscreteTaxonData >(  NULL ),
//    lnProb( 0.0 ),
//    storedLnProb( 0.0 ),
    num_sites( ns ),
    num_chars( nc ),
    num_site_rates( 1 ),
    num_site_matrices( 1 ),
    transition_prob_matrices( std::vector< std::vector<TransitionProbabilityMatrix> >(num_site_matrices, std::vector<TransitionProbabilityMatrix>(num_site_rates,TransitionProbabilityMatrix(num_chars) ) ) ),
    ambiguous_char_vector(),
    char_vector(),
    gap_vector(),
    pattern_counts(),
    num_patterns( num_sites ),
    compressed( false ),
    site_pattern( std::vector<size_t>(num_sites, 0) ),
    using_ambiguous_characters( false ),
    using_weighted_characters( false ),
    pattern_block_start( 0 ),
    pattern_block_end( num_patterns ),
    pattern_block_size( num_patterns ),
    template_state()
{

    // initialize with default parameters
    rate_matrix                 = NULL;
    site_rate_matrices          = NULL;
    root_frequencies            = NULL;
    site_rates                  = NULL;
    site_matrix_probs           = NULL;
    site_rates_probs            = NULL;
    process_time                = NULL;

    // flags specifying which model variants we use
    rate_variation_across_sites = false;


    // initially we use only a single processor until someone else tells us otherwise
    this->setActivePID( this->pid, 1 );

}



/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
template<class charType>
RevBayesCore::CTMCProcess<charType>::~CTMCProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

}


template<class charType>
RevBayesCore::CTMCProcess<charType>* RevBayesCore::CTMCProcess<charType>::clone( void ) const
{
    
    return new CTMCProcess( *this );
}


template<class charType>
double RevBayesCore::CTMCProcess<charType>::computeLnProbability( void )
{
    // @todo: #thread
    // This part should be done on several threads if possible
    // That means we should probabily call this function as a job,
    // where a job is defined as computing the lnProbability for a subset of the data (block)
    // Sebastian: this call is very slow; a lot of work happens in nextCycle()
    

    std::vector<double> likelihoods = std::vector<double>(num_patterns, 0.0);
    computeSiteLikelihoods( likelihoods );
            
    // sum the partials up
    double ln_prob = 0.0;
    for (size_t i=0; i<num_patterns; ++i)
    {
        ln_prob += pattern_counts[i] * log( likelihoods[i] );
    }

    return ln_prob;
}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::computeSiteLikelihoods( std::vector< double > &rv ) const
{
    std::vector< std::vector< std::vector< double > > > likelihoods = std::vector< std::vector< std::vector< double > > >( num_patterns, std::vector< std::vector< double > >( num_site_matrices, std::vector<double>(num_site_rates, 0.0) ) );
    computeSiteLikelihoodsPerSiteRateAndMatrix( likelihoods );
    
    
    std::vector<double> rates_probs(num_site_rates, 1.0/num_site_rates);
    std::vector<double> matrix_probs(num_site_matrices, 1.0/num_site_matrices);

    if ( site_rates_probs != NULL )
    {
        rates_probs = site_rates_probs->getValue();
    }

    if ( site_matrix_probs != NULL )
    {
        matrix_probs = site_matrix_probs->getValue();
    }
    
    // iterate over the number of sites
    for (size_t site_index=0; site_index<num_patterns; ++site_index)
    {
        double site_likelihood = 0;
        
        // iterate over the site matrices
        for (size_t matrix_index=0; matrix_index<num_site_matrices; ++matrix_index)
        {
            
            double matrix_likelihood = 0.0;
            // iterate over the site rates
            for (size_t rate_index=0; rate_index<num_site_rates; ++rate_index)
            {
                matrix_likelihood += likelihoods[site_index][matrix_index][rate_index] * rates_probs[rate_index];
            }
            
            site_likelihood += matrix_likelihood * matrix_probs[matrix_index];
        }
        rv[site_index] = site_likelihood;
    }
    
}



template<class charType>
void RevBayesCore::CTMCProcess<charType>::computeSiteLikelihoodsPerSiteMatrix( std::vector< std::vector< double > > &rv ) const
{
    std::vector< std::vector< std::vector< double > > > likelihoods = std::vector< std::vector< std::vector< double > > >( num_patterns, std::vector< std::vector< double > >( num_site_matrices, std::vector<double>(num_site_rates, 0.0) ) );
    computeSiteLikelihoodsPerSiteRateAndMatrix( likelihoods );
    
    
    std::vector<double> rates_probs(num_site_rates, 1.0/num_site_rates);

    if ( site_rates_probs != NULL )
    {
        rates_probs = site_rates_probs->getValue();
    }
    
    // iterate over the number of sites
    for (size_t site_index=0; site_index<num_patterns; ++site_index)
    {
        
        // iterate over the site matrices
        for (size_t matrix_index=0; matrix_index<num_site_matrices; ++matrix_index)
        {
            
            double matrix_likelihood = 0.0;
            // iterate over the site rates
            for (size_t rate_index=0; rate_index<num_site_rates; ++rate_index)
            {
                matrix_likelihood += likelihoods[site_index][matrix_index][rate_index] * rates_probs[rate_index];
            }
            
            rv[site_index][matrix_index] = matrix_likelihood;
        }
    }
    
}

template<class charType>
void RevBayesCore::CTMCProcess<charType>::computeSiteLikelihoodsPerSiteRate( std::vector< std::vector< double > > &rv ) const
{
    std::vector< std::vector< std::vector< double > > > likelihoods = std::vector< std::vector< std::vector< double > > >( num_patterns, std::vector< std::vector< double > >( num_site_matrices, std::vector<double>(num_site_rates, 0.0) ) );
    computeSiteLikelihoodsPerSiteRateAndMatrix( likelihoods );
    
    
    std::vector<double> matrix_probs(num_site_matrices, 1.0/num_site_matrices);


    if ( site_matrix_probs != NULL )
    {
        matrix_probs = site_matrix_probs->getValue();
    }
    
    // iterate over the number of sites
    for (size_t site_index=0; site_index<num_patterns; ++site_index)
    {

        // iterate over the site rates
        for (size_t rate_index=0; rate_index<num_site_rates; ++rate_index)
        {
            
            double rate_likelihood = 0.0;
            // iterate over the site matrices
            for (size_t matrix_index=0; matrix_index<num_site_matrices; ++matrix_index)
            {
                rate_likelihood += likelihoods[site_index][matrix_index][rate_index] * matrix_probs[matrix_index];
            }
            rv[site_index][rate_index] = rate_likelihood;
        }
        
    }
    
}

template<class charType>
void RevBayesCore::CTMCProcess<charType>::computeSiteLikelihoodsPerSiteRateAndMatrix( std::vector< std::vector< std::vector< double > > > &rv ) const
{
    
    // compute the transition probabilities
    this->updateTransitionProbabilities();
    
    std::vector< std::vector<double> > rf;
    getRootFrequencies( rf );

    // iterate over the site matrices
    for (size_t matrix_index=0; matrix_index<num_site_matrices; ++matrix_index)
    {
        
        // iterate over the site rates
        for (size_t rate_index=0; rate_index<num_site_rates; ++rate_index)
        {
            
            const double* tp_begin = this->transition_prob_matrices[matrix_index][rate_index].theMatrix;

            // iterate over the number of sites
            for (size_t site_index=0; site_index<num_patterns; ++site_index)
            {
                
                // is this site a gap?
                if ( gap_vector[site_index] )
                {
                    // since this is a gap we need to assume that the actual state could have been any state
                    
                    // store the likelihood
                    rv[site_index][matrix_index][rate_index] = 1.0;
                }
                else // we have observed a character
                {

                    double this_site_prob = 0.0;
                    // iterate over all possible initial states
                    for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                    {

                        if ( this->using_ambiguous_characters == true && this->using_weighted_characters == false)
                        {
                            // compute the likelihood that we had a transition from state c1 to the observed state org_val
                            // note, the observed state could be ambiguous!
                            const RbBitSet &val = ambiguous_char_vector[site_index];

                            // get the pointer to the transition probabilities for the terminal states
                            const double* d  = tp_begin+(this->num_chars*c1);

                            double tmp = 0.0;

                            for ( size_t i=0; i<val.size(); ++i )
                            {
                                // check whether we observed this state
                                if ( val.test(i) == true )
                                {
                                    // add the probability
                                    tmp += *d;
                                }

                                // increment the pointer to the next transition probability
                                ++d;
                            } // end-for over all observed states for this character

                            // store the likelihood
                            this_site_prob += tmp * rf[matrix_index][c1];
                                        
                        }
//                        else if ( this->using_weighted_characters == true )
//                        {
//                            // compute the likelihood that we had a transition from state c1 to the observed state org_val
//                            // note, the observed state could be ambiguous!
////                            const RbBitSet &val = amb_char_node[site];
//                            size_t this_site_index = site_indices[site];
//                            const RbBitSet &val = this->value->getCharacter(0, this_site_index).getState();
//
//                            // get the pointer to the transition probabilities for the terminal states
//                            const double* d = tp_begin+(this->num_chars*c1);
//
//                            double tmp = 0.0;
//                            const std::vector< double >& weights = this->value->getCharacter(char_data_node_index, this_site_index).getWeights();
//                            for ( size_t i=0; i<val.size(); ++i )
//                            {
//                                // check whether we observed this state
//                                if ( val.test(i) == true )
//                                {
//                                    // add the probability
//                                    tmp += *d * weights[i] ;
//                                }
//
//                                // increment the pointer to the next transition probability
//                                ++d;
//                            } // end-for over all observed states for this character
//
//                            // store the likelihood
//                            this_site_prob += tmp * rf[matrix_index][c1];
//
//                        }
                        else // no ambiguous characters in use
                        {
                            unsigned long org_val = char_vector[site_index];

                            // store the likelihood
                            this_site_prob += tp_begin[c1*this->num_chars+org_val] * rf[matrix_index][c1];

                        }

                    } // end-for over all possible initial character for the branch
                    rv[site_index][matrix_index][rate_index] = this_site_prob;

                } // end-if a gap state vs observed state
                
            } // end for over all sites
            
        } // end for over all site rates
        
    } // end for over all Q matrices
    
}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::compress( void )
{

    // only if the value has been set
    if ( this->value == NULL )
    {
        return;
    }

    ambiguous_char_vector.clear();
    char_vector.clear();
    gap_vector.clear();
    pattern_counts.clear();
    num_patterns = 0;

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> site_indices = getIncludedSiteIndices();

    // check whether there are ambiguous characters (besides gaps)
    bool ambiguous_characters = false;

    AbstractDiscreteTaxonData& seq = *(this->value);
    // find the unique site patterns and compute their respective frequencies
    for (size_t site = 0; site < num_sites; ++site)
    {

        DiscreteCharacterState &c = seq.getCharacter(site_indices[site]);

        // if we treat unknown characters as gaps and this is an unknown character then we change it
        // because we might then have a pattern more
        if ( c.isGapState() == false && (c.isAmbiguous() || c.isMissingState()) )
        {
            ambiguous_characters = true;
            break;
        }
    }

    // set the global variable if we use ambiguous characters
    using_ambiguous_characters = ambiguous_characters;

    // check whether there are weighted characters
    bool weighted_characters = false;

    for (size_t site = 0; site < num_sites; ++site)
    {
        
        DiscreteCharacterState &c = seq.getCharacter(site_indices[site]);

        if ( c.isWeighted() )
        {
            weighted_characters = true;
            break;
        }
    }

    // set the global variable if we use ambiguous characters
    using_weighted_characters = weighted_characters;

    std::vector<bool> unique(num_sites, true);
    std::vector<size_t> index_of_site_pattern;

    // compress the character matrix if we're asked to
    if ( compressed == true )
    {
        // find the unique site patterns and compute their respective frequencies
        std::map<std::string,size_t> patterns;
        for (size_t site = 0; site < num_sites; ++site)
        {
            // create the site pattern
            std::string pattern = "";
            
            DiscreteCharacterState &c = seq.getCharacter(site_indices[site]);
            pattern += c.getStringValue();
            
            // check if we have already seen this site pattern
            std::map<std::string, size_t>::const_iterator index = patterns.find( pattern );
            if ( index != patterns.end() )
            {
                // we have already seen this pattern
                // increase the frequency counter
                pattern_counts[ index->second ]++;

                // obviously this site isn't unique nor the first encounter
                unique[site] = false;

                // remember which pattern this site uses
                site_pattern[site] = index->second;
            }
            else
            {
                // create a new pattern frequency counter for this pattern
                pattern_counts.push_back(1);

                // insert this pattern with the corresponding index in the map
                patterns.insert( std::pair<std::string,size_t>(pattern,num_patterns) );

                // remember which pattern this site uses
                site_pattern[site] = num_patterns;

                // increase the pattern counter
                num_patterns++;

                // add the index of the site to our pattern-index vector
                index_of_site_pattern.push_back( site );

                // flag that this site is unique (or the first occurence of this pattern)
                unique[site] = true;
            }
        }
    }
    else
    {
        // we do not compress
        num_patterns = num_sites;
        pattern_counts     = std::vector<size_t>(num_sites,1);
        index_of_site_pattern = std::vector<size_t>(num_sites,1);
        for (size_t i = 0; i < this->num_sites; i++)
        {
            index_of_site_pattern[i] = i;
        }
    }


    // compute which block of the data this process needs to compute
    pattern_block_start = size_t(floor( (double(this->pid   - this->active_PID) / this->num_processes ) * num_patterns) );
    pattern_block_end   = size_t(floor( (double(this->pid+1 - this->active_PID) / this->num_processes ) * num_patterns) );
    pattern_block_size  = pattern_block_end - pattern_block_start;


    std::vector<size_t> process_pattern_counts = std::vector<size_t>(pattern_block_size,0);
    // allocate and fill the cells of the vectors
    // resize the vectors
    ambiguous_char_vector.resize(pattern_block_size);
    char_vector.resize(pattern_block_size);
    gap_vector.resize(pattern_block_size);
    for (size_t pattern_index = 0; pattern_index < pattern_block_size; ++pattern_index)
    {
        // set the counts for this patter
        process_pattern_counts[pattern_index] = pattern_counts[pattern_index+pattern_block_start];

        charType &c = static_cast<charType &>( seq.getCharacter(site_indices[index_of_site_pattern[pattern_index+pattern_block_start]]) );
        gap_vector[pattern_index] = c.isGapState();

        if ( using_ambiguous_characters == true )
        {
            // we use the actual state
            ambiguous_char_vector[pattern_index] = c.getState();
        }
        else if ( c.isGapState() == false )
        {
            // we use the index of the state
            char_vector[pattern_index] = c.getStateIndex();
            if ( c.getStateIndex() >= this->num_chars )
                throw RbException("Problem with state index in CTMC-Process!");
        }
        else
        {
            // just to be safe
            char_vector[pattern_index] = -1;
        }

    }

    // now copy back the pattern count vector
    pattern_counts = process_pattern_counts;

}
        
        

template<class charType>
void RevBayesCore::CTMCProcess<charType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{

    if ( n == "siteLikelihoods" )
    {

        // make sure the likelihoods are updated
        const_cast<CTMCProcess<charType> *>( this )->computeLnProbability();

        // get the per site likelihood
        RbVector<double> tmp = RbVector<double>(num_patterns, 0.0);
        computeSiteLikelihoods( tmp );

        // now match it back to the actual sites
        if ( this->compressed == true && num_sites > num_patterns )
        {
            rv = RbVector<double>(num_sites, 0.0);

            for (size_t i=0; i<num_sites; ++i)
            {
                size_t pattern_index = site_pattern[i];
                rv[i] = tmp[pattern_index] / pattern_counts[pattern_index];
            }
        }
        else
        {
            rv = tmp;
        }

    }
    else if ( n == "siteRates" )
    {

        // make sure the likelihoods are updated
        const_cast<CTMCProcess<charType> *>( this )->computeLnProbability();

        // get the site rates
        std::vector<double> r;
        if ( this->rate_variation_across_sites == true )
        {
            r = this->site_rates->getValue();
        }
        else
        {
            r.push_back(1.0);
        }


        // get the per site rate likelihood
        MatrixReal tmp = MatrixReal(num_patterns, num_site_rates, 0.0);
//        computeSiteLikelihoodsPerSiteRate( tmp );

        // now match it back to the actual sites
        MatrixReal siteRateConditionalProb = MatrixReal(num_sites, num_site_rates, 0.0);
        for (size_t i=0; i<num_sites; ++i)
        {
            size_t pattern_index = site_pattern[i];
            double siteLikelihoods = 0.0;

            for (size_t j=0; j<num_site_rates; ++j)
            {
                siteRateConditionalProb[i][j] = exp(tmp[pattern_index][j] / pattern_counts[pattern_index]);
                siteLikelihoods += siteRateConditionalProb[i][j];
            }
            for (size_t j=0; j<num_site_rates; ++j)
            {
                siteRateConditionalProb[i][j] = siteRateConditionalProb[i][j] / siteLikelihoods;
            }
        }

        rv = RbVector<double>(num_sites, 0.0);

        // now either sample the site rates according to the conditional posterior probability of each rate category
        // or compute them as the posterior mean
        std::string method = static_cast<const TypedDagNode<std::string>* >( args[0] )->getValue();

        if (method == "sampling")
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            for (size_t i=0; i<num_sites; ++i)
            {
                size_t rateIndex = 0;
                double u = rng->uniform01();

                while ( true )
                {
                    u -= siteRateConditionalProb[i][rateIndex];

                    if ( u > 0.0 )
                    {
                        ++rateIndex;
                    }
                    else
                    {
                        break;
                    }
                }

                rv[i] = r[rateIndex];

            }

        }
        else if (method == "weightedAverage")
        {
            for (size_t i=0; i<num_sites; ++i)
            {
                for (size_t rateIndex=0; rateIndex < num_site_rates; ++rateIndex)
                {
                    rv[i] += r[rateIndex] * siteRateConditionalProb[i][rateIndex];
                }
            }
        }

    }
    else
    {
        throw RbException() << "The CTMC process does not have a member method called '" << n << "'.";
    }

}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, MatrixReal &rv) const
{

    if ( n == "siteRateLikelihoods" )
    {

        // make sure the likelihoods are updated
        const_cast<CTMCProcess<charType> *>( this )->computeLnProbability();

        // get the per site rate likelihood
        MatrixReal tmp = MatrixReal(num_patterns, num_site_rates, 0.0);
//        computeSiteLikelihoodsPerSiteRate( tmp );

        // now match it back to the actual sites
        rv = MatrixReal(num_sites, num_site_rates, 0.0);
        for (size_t i=0; i<num_sites; ++i)
        {
            size_t pattern_index = site_pattern[i];
            for (size_t j=0; j<num_site_rates; ++j)
            {
                rv[i][j] = tmp[pattern_index][j] / pattern_counts[pattern_index];
            }
        }

    }
    else if ( n == "siteMatrixLikelihoods" )
    {

        // make sure the likelihoods are updated
        const_cast<CTMCProcess<charType> *>( this )->computeLnProbability();

        // get the per site rate likelihood
        MatrixReal tmp = MatrixReal(num_patterns, num_site_matrices, 0.0);
//        computeSiteLikelihoodsPerSiteMixture( tmp );

        // now match it back to the actual sites
        rv = MatrixReal(num_sites, num_site_matrices, 0.0);
        for (size_t i=0; i<num_sites; ++i)
        {
            size_t pattern_index = site_pattern[i];
            for (size_t j=0; j<num_site_matrices; ++j)
            {
                rv[i][j] = tmp[pattern_index][j] / pattern_counts[pattern_index];
            }
        }

    }
    else
    {
        throw RbException() << "The CTMC process does not have a member method called '" << n << "'.";
    }

}




template<class charType>
std::vector<size_t> RevBayesCore::CTMCProcess<charType>::getIncludedSiteIndices( void )
{
//    return this->value->getIncludedSiteIndices();
    return std::vector<size_t>(true,this->value->getNumberOfCharacters());
}



template<class charType>
void RevBayesCore::CTMCProcess<charType>::getRootFrequencies( std::vector<std::vector<double> >& rf ) const
{
    if ( root_frequencies != NULL )
    {
        std::vector<double> f = root_frequencies->getValue();
        rf.push_back( f );
    }
    else if (site_rate_matrices !=  NULL)
    {
        for (size_t matrix = 0; matrix < this->num_site_matrices; matrix++)
        {
            const RateMatrix *rm = dynamic_cast<const RateMatrix *>(&site_rate_matrices->getValue()[matrix]);
            if ( rm != NULL )
            {
                rf.push_back( rm->getStationaryFrequencies() );
            }
            else
            {
                throw RbException("If you want to use RateGenerators that are not RateMatrices then you need to specify the root frequencies directly.");
            }
            
        }
    }
    else if (rate_matrix != NULL)
    {
        const RateMatrix *rm = dynamic_cast<const RateMatrix *>(&rate_matrix->getValue());
        if ( rm != NULL )
        {
            rf.push_back( rm->getStationaryFrequencies() );
        }
        else
        {
            throw RbException("If you want to use RateGenerators that are not RateMatrices then you need to specify the root frequencies directly.");
        }

    }
    else
    {
        rf.push_back( std::vector<double>(num_chars, 1.0/num_chars) );
    }
}

template<class charType>
std::vector<double> RevBayesCore::CTMCProcess<charType>::getRootFrequencies( size_t mixture ) const
{
    if (mixture > this->num_site_matrices)
    {
        throw RbException("Site mixture index out of bounds");
    }

    std::vector<std::vector<double> > rf;
    getRootFrequencies(rf);

    return rf[mixture];
}


template<class charType>
std::vector<double> RevBayesCore::CTMCProcess<charType>::getMixtureProbs( void ) const
{
    std::vector<double> rates_probs(num_site_rates, 1.0/num_site_rates);
    std::vector<double> matrix_probs(num_site_matrices, 1.0/num_site_matrices);
    size_t num_site_mixtures = num_site_rates * num_site_matrices;
    std::vector<double> probs(num_site_mixtures, 1.0/num_site_mixtures);

    if ( site_rates_probs != NULL )
    {
        rates_probs = site_rates_probs->getValue();
    }

    if ( site_matrix_probs != NULL )
    {
        matrix_probs = site_matrix_probs->getValue();
    }

    for (size_t matrix = 0; matrix < num_site_matrices; ++matrix)
    {
        for (size_t j = 0; j < this->num_site_rates; ++j)
        {
            probs[j * num_site_matrices + matrix] = matrix_probs[matrix] * rates_probs[j];
        }
    }

    return probs;
}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::redrawValue( void )
{
    
    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new DiscreteTaxonData<charType>( Taxon("") );
    
    // prepare the rate matrices
    updateTransitionProbabilities();
    
    // get the global random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    std::vector<std::vector<double> > freqs;
    getRootFrequencies(freqs);
    // simulate the root sequence
    DiscreteTaxonData< charType > root = DiscreteTaxonData<charType>( Taxon("") );
    for ( size_t i = 0; i < num_sites; ++i )
    {
        size_t matrix_index = 0;
        size_t rate_index   = 0;
        
        // draw the matrix index
        double u = 0.0;
        if ( num_site_matrices > 1 )
        {
            // draw the rate for this site
            u = rng->uniform01();

            std::vector<double> matrix_probs = getMixtureProbs();

            std::vector< double >::const_iterator freq = matrix_probs.begin();
            while ( true )
            {
                u -= *freq;

                if ( u > 0.0 )
                {
                    ++matrix_index;
                    ++freq;
                }
                else
                {
                    break;
                }
            }
        }
        
        // draw the rate index
        if ( num_site_rates > 1 )
        {
            // draw the rate for this site
            u = rng->uniform01();

            const std::vector<double>& rates_probs = site_rates_probs->getValue();

            std::vector< double >::const_iterator freq = rates_probs.begin();
            while ( true )
            {
                u -= *freq;

                if ( u > 0.0 )
                {
                    ++rate_index;
                    ++freq;
                }
                else
                {
                    break;
                }
            }

        }
        
        
        const std::vector< double > &stationary_freqs = freqs[matrix_index % freqs.size()];

        // create the character
        size_t root_state = 0;

        // draw the state
        u = rng->uniform01();
        std::vector< double >::const_iterator freq = stationary_freqs.begin();
        while ( true )
        {
            u -= *freq;
            if ( u > 0.0 )
            {
                ++root_state;
                ++freq;
            }
            else
            {
                break;
            }

        }
        
        TransitionProbabilityMatrix& p = transition_prob_matrices[matrix_index][rate_index];
        double *freqs = p[ root_state ];

        // create the character
        charType sim_c = charType( num_chars );
        sim_c.setToFirstState();
        size_t state_index = 0;
        
        // draw the state
        u = rng->uniform01();
        size_t stateIndex = 0;
        while ( true )
        {
            u -= *freqs;
            ++state_index;

            if ( u > 0.0 && state_index < this->num_chars)
            {
                ++sim_c;
                ++freqs;
            }
            else
            {
                break;
            }

        }

        // add the character to the sequence
        this->value->addCharacter( sim_c );
        
    }
}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::reInitialized( void )
{

    // we need to recompress because the tree may have changed
    compress();
}


template <class charType>
void RevBayesCore::CTMCProcess<charType>::setActivePIDSpecialized(size_t a, size_t n)
{

    // we need to recompress the data
    this->compress();
}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::setValue(charType *v, bool force)
{

    if ( v->getMaxObservedStateIndex() > this->num_chars - 1)
    {
        // We might use different sized matrices for different partitions depending on the observed number of states.
        std::stringstream ss;
        ss << "The number of observed states (" << v->getMaxObservedStateIndex() + 1 << ") is greater than the dimension of the Q matrix (" << this->num_chars << ")" << std::endl;
        throw RbException(ss.str());
    }

    if (v->getDataType() != this->template_state.getDataType() )
    {
      // There is a mismatch between the data type of the data matrix and the data type of the CTMC.
      std::stringstream ss;
      ss << "The data type of the data matrix ("<< v->getDataType() <<") differs from that of the CTMC object ("<< this->template_state.getDataType() <<")." << std::endl;
      throw RbException(ss.str());
    }

    // delegate to the parent class
    TypedDistribution< AbstractDiscreteTaxonData >::setValue(v, force);

    // reset the number of sites
    this->num_sites = v->getNumberOfIncludedCharacters();

    site_pattern.clear();
    site_pattern.resize(num_sites);

    // now compress the data and resize the likelihood vectors
    this->compress();

    // now we also set the template state
    template_state = charType( static_cast<const charType&>( this->value->getCharacter(0) ) );
    template_state.setToFirstState();
    template_state.setGapState( false );
    template_state.setMissingState( false );

}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::setRateMatrix(const TypedDagNode< RateGenerator > *rm)
{

    // remove the old parameter first
    if ( rate_matrix != NULL )
    {
        this->removeParameter( rate_matrix );
        rate_matrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( site_rate_matrices );
        site_rate_matrices = NULL;
    }

    // set the value
    rate_matrix = rm;
    num_site_matrices = 1;

    if (rm != NULL && rm->getValue().size() != this->num_chars)
    {
        std::stringstream ss;
        ss << "Rate generator dimensions (" << rm->getValue().size() << " do not match the number of character states (" << this->num_chars << ")";
        throw(RbException(ss.str()));
    }

    // add the new parameter
    this->addParameter( rate_matrix );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::setRateMatrix(const TypedDagNode< RbVector< RateGenerator > > *rm)
{

    // remove the old parameter first
    if ( rate_matrix != NULL )
    {
        this->removeParameter( rate_matrix );
        rate_matrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( site_rate_matrices );
        site_rate_matrices = NULL;
    }

    // set the value
    site_rate_matrices = rm;
    num_site_matrices = rm == NULL ? 1 : rm->getValue().size();

    if (rm != NULL && rm->getValue()[0].size() != this->num_chars)
    {
        std::stringstream ss;
        ss << "Rate generator dimensions (" << rm->getValue()[0].size() << " do not match the number of character states (" << this->num_chars << ")";
        throw(RbException(ss.str()));
    }

    // add the new parameter
    this->addParameter( site_rate_matrices );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::setRootFrequencies(const TypedDagNode< Simplex > *f)
{

    // remove the old parameter first
    if ( root_frequencies != NULL )
    {
        this->removeParameter( root_frequencies );
        root_frequencies = NULL;
    }

    if ( f != NULL )
    {
        // set the value
        root_frequencies = f;
    }

    // add the new parameter
    this->addParameter( root_frequencies );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::setSiteRates(const TypedDagNode< RbVector< double > > *r)
{

    // remove the old parameter first
    if ( site_rates != NULL )
    {
        this->removeParameter( site_rates );
        site_rates = NULL;
    }

    if ( r != NULL )
    {
        // set the value
        rate_variation_across_sites = true;
        site_rates = r;
        this->num_site_rates = r->getValue().size();
    }
    else
    {
        // set the value
        rate_variation_across_sites = false;
        site_rates = NULL;
        this->num_site_rates = 1;
    }

    // add the new parameter
    this->addParameter( site_rates );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}


template<class charType>
void RevBayesCore::CTMCProcess<charType>::setSiteRatesProbs(const TypedDagNode< Simplex > *rp)
{

    // remove the old parameter first
    if ( site_rates_probs != NULL )
    {
        this->removeParameter( site_rates_probs );
        site_rates_probs = NULL;
    }

    if (rp != NULL)
    {
        // set the value
        site_rates_probs = rp;
    }

    // add the new parameter
    this->addParameter( site_rates_probs );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}

template<class charType>
void RevBayesCore::CTMCProcess<charType>::setSiteMatricesProbs(const TypedDagNode< Simplex > *s)
{

    // remove the old parameter first
    if ( site_matrix_probs != NULL )
    {
        this->removeParameter( site_matrix_probs );
        site_matrix_probs = NULL;
    }

    // set the value
    site_matrix_probs = s;

    // add the new parameter
    this->addParameter( site_matrix_probs );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}


/** Swap a parameter of the distribution */
template<class charType>
void RevBayesCore::CTMCProcess<charType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == rate_matrix)
    {
        rate_matrix = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    }
    else if (oldP == site_rate_matrices)
    {
        site_rate_matrices = static_cast<const TypedDagNode< RbVector< RateGenerator > >* >( newP );
    }
    else if (oldP == root_frequencies)
    {
        root_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == site_matrix_probs)
    {
        site_matrix_probs = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == site_rates)
    {
        site_rates = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else if (oldP == site_rates_probs)
    {
        site_rates_probs = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

}


/*
 * Update the transition probability matrices for the branch attached to the given node index.
 */
template<class charType>
void RevBayesCore::CTMCProcess<charType>::updateTransitionProbabilities( void ) const
{

    double time = 1.0;
    
    if ( process_time != NULL )
    {
        time = process_time->getValue();
    }
    
    // first, get the rate matrix for this branch
    const RateGenerator *rm = NULL;

    for (size_t matrix_index = 0; matrix_index < this->num_site_matrices; ++matrix_index)
    {
        bool free_mem = false;
        if ( this->site_rate_matrices != NULL )
        {
            rm = &this->site_rate_matrices->getValue()[matrix_index];
        }
        else if ( this->rate_matrix != NULL )
        {
            rm = &this->rate_matrix->getValue();
        }
        else
        {
            RateMatrix_JC jc(this->num_chars);
            rm = new RateMatrix_JC(this->num_chars);
            
            free_mem = true;
        }

        for (size_t rate_index = 0; rate_index < this->num_site_rates; ++rate_index)
        {
            double r = 1.0;
            if ( this->rate_variation_across_sites == true )
            {
                r = this->site_rates->getValue()[rate_index];
            }
            rm->calculateTransitionProbabilities( time, 0.0,  r, this->transition_prob_matrices[matrix_index][rate_index] );
        }
        
        if (free_mem == true)
        {
            delete rm;
        }
    }
    
}

#endif
