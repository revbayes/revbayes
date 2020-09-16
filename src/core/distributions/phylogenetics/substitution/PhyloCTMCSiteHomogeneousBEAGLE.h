//---------------------------------------------------------------------------------------[ Header ]
//--{1

#ifndef PhyloCTMCSiteHomogeneousBEAGLE_H
#define PhyloCTMCSiteHomogeneousBEAGLE_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "TopologyNode.h"

#include "RbVector.h"

//-- TODO: are these needed?
#include "DnaState.h"
#include "RateMatrix.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"


namespace RevBayesCore
{

    template<class charType>
    class PhyloCTMCSiteHomogeneousBEAGLE : public AbstractPhyloCTMCSiteHomogeneous<charType>
    {

        public:

            //----====  Constructors  ====----

            //-- Default constructor
            PhyloCTMCSiteHomogeneousBEAGLE ( const TypedDagNode<Tree>* t
                                           , size_t nChars
                                           , bool c
                                           , size_t nSites
                                           , bool amb
                                           , bool internal
                                           , bool gapmatch
                                           );

            //-- Clone constructor
            PhyloCTMCSiteHomogeneousBEAGLE* clone ( void ) const;


            //-- Destructor
            virtual ~PhyloCTMCSiteHomogeneousBEAGLE ( void );



        protected:

            //----====  Protected Methods  ====----

            //-- Method used to call the appropriate likelihood function.
            virtual void   executeMethod                 ( const std::string &n
                                                         , const std::vector<const DagNode*> &args
                                                         , RbVector<double> &rv) const;

            //-- Method used to call the appropriate likelihood function.
            virtual void   executeMethod                 ( const std::string &n
                                                         , const std::vector<const DagNode*> &args
                                                         , MatrixReal &rv) const;

            //-- Calculate the lnLikelihood of a tree given the model.
            virtual double computeLnProbability          ( void );

            //-- BEAGLE compute lnLikelihood of a rooted tree.
            virtual void   computeRootLikelihood         ( size_t root
                                                         , size_t l
                                                         , size_t r
                                                         );

            //-- BEAGLE compute lnLikelihood of an unrooted tree.
            virtual void   computeRootLikelihood         ( size_t root
                                                         , size_t l
                                                         , size_t r
                                                         , size_t m
                                                         );

            //-- Collect a BEAGLE operation for an internal node into the computation queue.
            virtual void   computeInternalNodeLikelihood ( const TopologyNode &n
                                                         , size_t nIdx
                                                         , size_t l
                                                         , size_t r
                                                                      );

            //-- Collect a BEAGLE operation for a leaf node into the computation queue.
            virtual void   computeTipLikelihood          ( const TopologyNode &node
                                                         , size_t nIdx
                                                         );


        private:

            //----====  Private Variables  ====----

            //-- Keep track if a BEAGLE instance exists.
            bool isBeagleInitialized = false;

            //-- Accumulate BEAGLE lnLikelihood across all models.
            double ln_beagle_probability;


            //----====  Private Methods  ====----

            //-- Helper function to pretty print BEAGLE error codes.
        std::string parseReturnCode           ( int code );

            //-- Initialize a BEAGLE instance.
            void        initializeBeagleInstance  ( void );

            //-- Set the tree tip states for the BEAGLE instance.
            void        initializeBeagleTips      ( void );

            //-- Set BEAGLE buffers for site-model mappings.
            void        setSitePatterns           ( void );

            //-- Set BEAGLE buffers for the the stationary distribution of a model.
            void        setStationaryDistribution ( void );

            //-- Set BEAGLE buffers for the ASRV distribution for a model.
            void        setSiteRateDistribution   ( void );

            //-- Set BEAGLE buffers for the substitution models.
            void        setSubstitutionModels     ( void );

            //-- Calculate the tree branch length between adjacent tree nodes..
            double      calculateBranchLength     ( const TopologyNode &node
                                                  , size_t node_index
                                                  );

    };

} //-- End namespace

//--}

//-----------------------------------------------------------------------[ Imports / Constructors ]
//--{1

#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"

#include <cmath>
#include <cstring>


//-- Try to keep a clean(ish) file structure... Refer to the namespace of this class as 'This'.
//   Remember to '#undef This' at the end of header file and NEVER have '#include'
//   statements after this line.
#define This RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>



template<class charType>
This::PhyloCTMCSiteHomogeneousBEAGLE ( const TypedDagNode<Tree>* t
                                     , size_t nChars
                                     , bool c
                                     , size_t nSites
                                     , bool amb
                                     , bool internal
                                     , bool gapmatch
                                     )
    : AbstractPhyloCTMCSiteHomogeneous<charType> ( t
                                                 , nChars
                                                 , 1
                                                 , c
                                                 , nSites
                                                 , amb
                                                 , internal
                                                 , gapmatch
                                                 )
{

}



template<class charType>
This::~PhyloCTMCSiteHomogeneousBEAGLE ( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!

    // remove myself from the tree listeners
    if ( this->tau != NULL )
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener(this);
    }

    if ( RbSettings::userSettings().getUseBeagle() == true &&
         this->in_mcmc_mode == true && this->beagle_instance >= 0 )
    {
        #if defined ( RB_BEAGLE_DEBUG )
                RBOUT ( "Finalizing BEAGLE" );
        #endif /* RB_BEAGLE_DEBUG */

        beagleFinalizeInstance(this->beagle_instance);
        this->beagle_instance = -1;
    }
}



template<class charType>
This* This::clone ( void ) const
{
    return new PhyloCTMCSiteHomogeneousBEAGLE<charType>(*this);
}

//--}

//-----------------------------------------------------------------------------------[ Initialize ]
//--{1

template<class charType>
void This::initializeBeagleInstance ( void )
{
    std::stringstream ss;

    if ( RbSettings::userSettings().getUseBeagle() == true )
    {
        ss << std::endl;
        ss << "Using BEAGLE library v" << beagleGetVersion();
        ss << " for parallel likelihood evaluation (https://beagle-dev.github.io/)";
        ss << std::endl;

        int  b_resource            = RbSettings::userSettings().getBeagleResource();
        bool b_use_cpu_threading   = RbSettings::userSettings().getBeagleMaxCPUThreads() != 1
                                   ? true : false;
        bool b_use_scaling         = RbSettings::userSettings().getBeagleScalingMode() != "none"
                                   ? true : false;

        int  b_tipCount            = this->tau->getValue().getNumberOfTips();
        int  b_partialsBufferCount = 2 * this-> num_nodes
                                   + ( this->using_ambiguous_characters
                                     ? this->tau->getValue().getNumberOfTips()
                                     : 0
                                     );
        int  b_compactBufferCount  = this->tau->getValue().getNumberOfTips()
                                   - ( this->using_ambiguous_characters
                                     ? this->tau->getValue().getNumberOfTips()
                                     : 0
                                     );
        int  b_stateCount          = this->num_chars;
        int  b_patternCount        = this->pattern_block_size;
        int  b_eigenBufferCount    = this->num_site_mixtures * 2;
        int  b_matrixBufferCount   = this->num_nodes * 2;
        int  b_categoryCount       = this->num_site_rates;
        int  b_scaleBufferCount    = (b_use_scaling ? (this->num_nodes * 2) : 0);
        int* b_resourceList        = &b_resource;
        int  b_resourceCount       = 1;
        long b_preferenceFlags     = ( RbSettings::userSettings().getBeagleUseDoublePrecision()
                                     ? BEAGLE_FLAG_PRECISION_DOUBLE
                                     : BEAGLE_FLAG_PRECISION_SINGLE
                                     )
                                   | ( b_use_cpu_threading
                                     ? BEAGLE_FLAG_THREADING_CPP
                                     : 0
                                     );
        long b_requirementFlags    = 0;

        BeagleInstanceDetails b_return_info;

        if ( RbSettings::userSettings().getBeagleAuto() == true )
        {
            ss << "Running benchmarks to automatically select fastest BEAGLE resource... ";
            ss << std::endl;

            // select fastest resource
            BeagleBenchmarkedResourceList* rBList;
            rBList = beagleGetBenchmarkedResourceList( b_tipCount
                                                     , b_compactBufferCount
                                                     , b_stateCount
                                                     , b_patternCount
                                                     , b_categoryCount
                                                     , NULL  // resourceList
                                                     , 0     // resourceCount
                                                     , b_preferenceFlags
                                                     , b_requirementFlags
                                                     , b_eigenBufferCount
                                                     , 1     // partitionCount
                                                     , 0     // calculateDerivatives
                                                     , 0     // benchmarkFlags
                                                     );

            if (rBList != NULL)
            {
                b_resource = rBList->list[0].number;

                ss << "Using resource "
                   << rBList->list[0].number << ": " << rBList->list[0].name;
                if ( rBList->list[0].number != 0 )
                {
                    ss << " (" << rBList->list[0].performanceRatio << "x CPU)";
                }
                ss << std::endl;
            } else {
                ss << "Resource benchmarking failed, using resource "
		   << b_resource << ": " << rBList->list[0].name;
                ss << std::endl;
            }
        }

        this->beagle_instance = beagleCreateInstance( b_tipCount
                                                    , b_partialsBufferCount
                                                    , b_compactBufferCount
                                                    , b_stateCount
                                                    , b_patternCount
                                                    , b_eigenBufferCount
                                                    , b_matrixBufferCount
                                                    , b_categoryCount
                                                    , b_scaleBufferCount
                                                    , b_resourceList
                                                    , b_resourceCount
                                                    , b_preferenceFlags
                                                    , b_requirementFlags
                                                    , &b_return_info
                                                    );

        #if defined ( RB_BEAGLE_DEBUG )
            ss << "BEAGLE parameters"                                       << std::endl;
            ss << "\t" << "tipCount            : " << b_tipCount            << std::endl;
            ss << "\t" << "partialsBufferCount : " << b_partialsBufferCount << std::endl;
            ss << "\t" << "compactBufferCount  : " << b_compactBufferCount  << std::endl;
            ss << "\t" << "stateCount          : " << b_stateCount          << std::endl;
            ss << "\t" << "patternCount        : " << b_patternCount        << std::endl;
            ss << "\t" << "eigenBufferCount    : " << b_eigenBufferCount    << std::endl;
            ss << "\t" << "matrixBufferCount   : " << b_matrixBufferCount   << std::endl;
            ss << "\t" << "categoryCount       : " << b_categoryCount       << std::endl;
            ss << "\t" << "scaleBufferCount    : " << b_scaleBufferCount    << std::endl;
            ss << "\t" << "resource            : " << b_resource            << std::endl;
            ss << std::endl;
            ss << "BEAGLE instance: " << this->beagle_instance << std::endl;
        #endif /* RB_BEAGLE_DEBUG */

        if ( this->beagle_instance < 0 )
        {
            ss << "Failed to start BEAGLE instance. "
               << "Reverting to RevBayes likelihood calculator." << std::endl;
            RbSettings::userSettings().setUseBeagle(false);
        }
        else
        {
            ss << "Using BEAGLE resource " << b_return_info.resourceNumber << std::endl;
            ss << "\t" << "Rsrc Name : "   << b_return_info.resourceName   << std::endl;
            ss << "\t" << "Impl Name : "   << b_return_info.implName       << std::endl;
            ss << "\t" << "Flags :";
            ss << BeagleUtilities::printBeagleFlags(b_return_info.flags);
            ss << std::endl;

            if ( b_use_cpu_threading )
            {
                beagleSetCPUThreadCount( this->beagle_instance
                                       , RbSettings::userSettings().getBeagleMaxCPUThreads()
                                       );
            }

        }
    }
    RBOUT(ss.str());
}


template<class charType>
void This::initializeBeagleTips ( void )
{
    int*    b_inStates   = NULL;
    double* b_inPartials = NULL;
    int     b_tipIndex   = 0;

    if ( this->using_ambiguous_characters == true )
    {
        b_inPartials = new double[this->pattern_block_size * this->num_chars];
    }
    else
    {
        b_inStates = new int[this->pattern_block_size];
    }

    //-- Iterate over the all leaf nodes in the tree.
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for ( std::vector<TopologyNode*>::iterator it = nodes.begin() ; it != nodes.end() ; ++it )
    {
        if ( (*it)->isTip() )
        {
            // TODO: change this code here to allow for ambiguity codes
            b_tipIndex                                  = (*it)->getIndex();
            size_t data_tip_index                       = this->taxon_name_2_tip_index_map[(*it)->getName()];

            const std::vector<bool> &gap_node           = this->gap_matrix[data_tip_index];
            const std::vector<unsigned long> &char_node = this->char_matrix[data_tip_index];
            const std::vector<RbBitSet> &amb_char_node  = this->ambiguous_char_matrix[data_tip_index];

            // iterate over all sites
            for ( size_t b_pattern = 0; b_pattern < this->pattern_block_size; ++b_pattern )
            {
                // is this site a gap?
                if ( gap_node[b_pattern] )
                {
                    if ( this->using_ambiguous_characters == true )
                    {
                        for ( size_t c = 0; c < this->num_chars; ++c )
                        {
                            b_inPartials[b_pattern * this->num_chars + c] = 1.0;
                        }
                    }
                    else
                    {
                        b_inStates[b_pattern] = (int) this->num_chars;
                    }
                }
                else
                {
                    if ( this->using_ambiguous_characters == true )
                    {
                        for ( size_t c = 0; c < this->num_chars; ++c )
                        {
                            b_inPartials[b_pattern * this->num_chars + c] =
                                (amb_char_node[b_pattern].isSet(c) ? 1.0 : 0.0);
                        }
                    }
                    else
                    {
                        b_inStates[b_pattern] = (int) char_node[b_pattern];
                    }
                }
            }

            if ( this->using_ambiguous_characters == true )
            {
            int b_code_tip_partials =
                    beagleSetTipPartials( this->beagle_instance
                                        , b_tipIndex
                                        , b_inPartials
                                        );
                if ( b_code_tip_partials != 0 )
                {
                    throw RbException("Could not set tip partials for model"
                                     + this->parseReturnCode(b_code_tip_partials));
                }
            }
            else
            {
            int b_code_tip_states =
                    beagleSetTipStates( this->beagle_instance
                                      , b_tipIndex
                                      , b_inStates
                                      );
                if ( b_code_tip_states != 0 )
                {
                    throw RbException("Could not set tip states for model"
                                     + this->parseReturnCode(b_code_tip_states));
                }
            }

            #if defined ( RB_BEAGLE_DEBUG )
            std::stringstream ss;
            ss << "BEAGLE setTipStates, tip = " << b_tipIndex
               << ", block size = " << this->pattern_block_size
               << ", states = ";
            ss << std::endl;
            for ( size_t b_pattern = 0; b_pattern < this->pattern_block_size; ++b_pattern )
            {
                if ( this->using_ambiguous_characters == true )
                {
                    for ( size_t c = 0; c < this->num_chars; ++c )
		    {
                        ss << b_inPartials[b_pattern * this->num_chars + c];
                    }
                    ss << " ";
                }
                else
                {
                    ss << b_inStates[b_pattern];
                }
            }
            ss << std::endl;
            RBOUT(ss.str());
            #endif /* RB_BEAGLE_DEBUG_TIP_STATES */
        }
    }
    delete[] b_inStates;
}

//--}

//-----------------------------------------------------------------------[ Private Helper Methods ]
//--{1

template<class charType>
std::string This::parseReturnCode ( int code )
{
    std::string result;
    switch (code) {
    case 0:
      result = "Success!";
      break;
    case -1:
      result = "Unspecified error!";
      break;
    case -2:
      result = "Not enough memory could be allocated!";
      break;
    case -3:
      result = "Unspecified exception!";
      break;
    case -4:
      result = "The instance index is out of range, or the instance has not been created!";
      break;
    case -5:
      result = "One of the indices specified exceeded the range of the array!";
      break;
    case -6:
      result = "No resource matches requirements!";
      break;
    case -7:
      result = "No implementation matches requirements!";
      break;
    case -8:
      result = "Floating-point range exceeded!";
      break;
    default:
      result = "Error code not found!";
    }
    return result;
}


template<class charType>
double This::calculateBranchLength ( const TopologyNode &node, size_t node_index )
{
    double branch_time = node.getBranchLength();
    double rate        = 1.0;

    if ( this->branch_heterogeneous_clock_rates == true )
    {
        rate = this->heterogeneous_clock_rates->getValue()[node_index];
    }
    else if ( this->homogeneous_clock_rate != NULL )
    {
        rate = this->homogeneous_clock_rate->getValue();
    }

    //-- TODO: maybe needs extra checks here (invariable sites, etc.)...

    double branch_len = rate * branch_time;
    if ( branch_len < 0 )
    {
      throw RbException("Error : Negative branch length!");
    }
    return branch_len;
}

//--}

//------------------------------------------------------------------------------[ BEAGLE Updaters ]
//--{1

template<class charType>
void This::setSitePatterns ( void )
{
    this->b_inPatternWeights = std::vector<double>(this->pattern_block_size);
    for ( size_t b_pattern = 0; b_pattern < this->pattern_block_size; ++b_pattern )
    {
        this->b_inPatternWeights[b_pattern] = (double) this->pattern_counts[b_pattern];
    }

    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss; ss << "Setting Site Patterns...\n";
    ss << "\tPatterns : \n\t\t";
    for ( auto x : this->b_inPatternWeights )
    {
      ss << std::to_string(x) << " ";
    }
    ss << "\n";
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    int b_code_pattern_weight =
        beagleSetPatternWeights( this->beagle_instance
                               , &this->b_inPatternWeights[0]
                               );
    if ( b_code_pattern_weight != 0 )
    {
        throw RbException("Could not set pattern weights for model"
                         + this->parseReturnCode(b_code_pattern_weight));
    }
}


template<class charType>
void This::setStationaryDistribution ( void )
{
    //-- TODO : Make this work for multiple models.

    // get the root frequencies
    std::vector<std::vector<double>> ff;
    this->getRootFrequencies(ff);
    const std::vector<double>& b_f                     = ff[0];
    const double*              b_inStateFrequencies    = &b_f[0];
    int                        b_stateFrequenciesIndex = 0;        //-- model index

    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << "Setting Stationary Distribution...\n";
    ss << "\tRoot frequencies : \n\t\t";
    for ( auto x : b_f )
    {
      ss << std::to_string(x) << " ";
    }
    ss << "\n";
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    this->b_stateFrequenciesIndex = 0;
    int b_code_state_freq =
        beagleSetStateFrequencies( this->beagle_instance
                                 , b_stateFrequenciesIndex
                                 , b_inStateFrequencies
                                 );
    if ( b_code_state_freq != 0 )
    {
        throw RbException("Could not set state frequencies for model: "
                         + this->parseReturnCode(b_code_state_freq));
    }
}


template<class charType>
void This::setSiteRateDistribution ( void )
{
    if (this->rate_variation_across_sites)
    {
        this->b_inCategoryWeights = this->site_rates_probs->getValue();
        this->b_inCategoryRates   = this->site_rates->getValue();
    }
    else
    {
        this->b_inCategoryWeights = { 1.0 };
        this->b_inCategoryRates   = { 1.0 };
    }

    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << "Setting Site Rates...\n";
    ss << "\tCategory weights : \n\t\t";
    for ( auto x : this->b_inCategoryWeights )
    {
      ss << std::to_string(x) << " ";
    }
    ss << "\n";
    ss << "\tCategory rates : \n\t\t";
    for ( auto x : this->b_inCategoryRates )
    {
      ss << std::to_string(x) << " ";
    }
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    int b_code_weights;
    int b_code_rates;
    for ( auto model : this->b_model_indices )
    {
        b_code_weights =
            beagleSetCategoryWeights( this->beagle_instance
                                    , (int) model
                                    , &this->b_inCategoryWeights[0]
                                    );
        if ( b_code_weights != 0 )
        {
            throw RbException("Could not set category weights for model '" + std::to_string(model)
                              + "'. " + this->parseReturnCode(b_code_weights));
        }

        b_code_rates =
            beagleSetCategoryRatesWithIndex( this->beagle_instance
                                           , (int) model
                                           , &this->b_inCategoryRates[0]
                                           );
        if ( b_code_rates != 0 )
        {
            throw RbException("Could not set category rates for model '" + std::to_string(model)
                              + "'. " + this->parseReturnCode(b_code_rates));
        }
    }
}


template<class charType>
void This::setSubstitutionModels ( void )
{
    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << "Setting Models...\n";
    #endif /* RB_BEAGLE_DEBUG */

    //-- Clear the models.
    this->b_model_indices.clear();

    //-- Switch the active eigen system
    for ( size_t i = 0; i < this->active_eigen_system.size(); ++i )
    {
        if ( this->touched_eigen_system[i] == false )
        {
            this->touched_eigen_system[0] = true;
            this->active_eigen_system[i]  = (this->active_eigen_system[i] == 0 ? 1 : 0);
        }
    }

    RbVector<RateGenerator> rate_matrices;
    if ( this->homogeneous_rate_matrix )
    {
        rate_matrices.push_back(this->homogeneous_rate_matrix->getValue());
    }
    else
    {
        rate_matrices = this->heterogeneous_rate_matrices->getValue();
    }

    size_t              model_idx;
    EigenSystem*        eigen_system;
    std::vector<double> my_eigen_values;
    std::vector<double> flat_eigen_vectors;
    std::vector<double> flat_inv_eigen_vectors;

    int b_code_eigen_decomp;
    size_t num_models = rate_matrices.size();
    for ( size_t i = 0; i < num_models; ++i )
    {
        eigen_system = rate_matrices[i].getEigenSystem();
        model_idx    = i + this->active_eigen_system[i] * num_models;

        my_eigen_values        = eigen_system->getRealEigenvalues();
        //flat_eigen_vectors     = eigen_system->getEigenvectors().flattenMatrix();
        //flat_inv_eigen_vectors = eigen_system->getInverseEigenvectors().flattenMatrix();
        flat_inv_eigen_vectors = eigen_system->getEigenvectors().flattenMatrix();
        flat_eigen_vectors     = eigen_system->getInverseEigenvectors().flattenMatrix();

        #if defined ( RB_BEAGLE_DEBUG )
        ss << "\tEigenvalues : \n\t\t";
        for ( size_t j = 0; j < my_eigen_values.size(); ++j )
        {
            ss << my_eigen_values[j] << " ";
        }
        ss << "\n";
        ss << "\tEigenvectors : ";
        for ( size_t j = 0; j < flat_eigen_vectors.size(); ++j )
        {
            if ( (j % this->num_chars) == 0) { ss << "\n\t\t"; };
            ss << flat_eigen_vectors[j] << " ";

        }
        ss << "\n";
        ss << "\tInverse Eigenvectors : ";
        for ( size_t j = 0; j < flat_inv_eigen_vectors.size(); ++j )
        {
            if ( (j % this->num_chars) == 0) { ss << "\n\t\t"; };
            ss << flat_inv_eigen_vectors[j] << " ";

        }
        RBOUT(ss.str());
        #endif /* RB_BEAGLE_DEBUG */

        b_code_eigen_decomp =
            beagleSetEigenDecomposition( this->beagle_instance
                                       , model_idx
                                       , &flat_eigen_vectors[0]
                                       , &flat_inv_eigen_vectors[0]
                                       , &my_eigen_values[0]
                                       );
        if ( b_code_eigen_decomp != 0 )
        {
            throw RbException("Could not set eigen decomposition for model '" + std::to_string(model_idx)
                              + "'. " + this->parseReturnCode(b_code_eigen_decomp));
        }

        this->b_model_indices.push_back(model_idx);
    }
}

//--}

//----------------------------------------------------------------------[ Likelihood Calculations ]
//--{1

template<class charType>
double This::computeLnProbability( void )
{
    //-- Reset the lnProbability.
    this->lnProb = 0.0;

    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( this->tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
        this->dirty_nodes = std::vector<bool>(this->num_nodes, true);
    }

    //-- Make sure a BEAGLE instance has been created.
    if ( !this->isBeagleInitialized )
    {
      this->initializeBeagleInstance();
      this->initializeBeagleTips();

      //-- Note: Substitution model must be set before everything else.
      this->setSubstitutionModels();
      this->setSiteRateDistribution();

      //-- Note: The stationary distribution and site patterns are set only once.
      this->setStationaryDistribution();
      this->setSitePatterns();

      this->isBeagleInitialized = true;
    }

    //-- TODO: Only update if the model params change.
    if ( true )
    {
        this->setSubstitutionModels();
    }
    if ( true )
    {
        this->setSiteRateDistribution();
    }

    //-- Get the 'root' node.
    const TopologyNode& root = this->tau->getValue().getRoot();
    if ( root.getNumberOfChildren() == 2 )
    {
        this->computeRootLikelihood( root.getIndex()
                                   , root.getChild(0).getIndex()
                                   , root.getChild(1).getIndex()
                                   );

    }
    else if ( root.getNumberOfChildren() == 3 )
    {
        this->computeRootLikelihood( root.getIndex()
                                   , root.getChild(0).getIndex()
                                   , root.getChild(1).getIndex()
                                   , root.getChild(2).getIndex()
                                   );
    }
    else
    {
        throw RbException("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
    }

    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << "\tCalculated : " << this->ln_beagle_probability << "\n";
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    return this->ln_beagle_probability;
}



template<class charType>
void This::computeRootLikelihood ( size_t root, size_t left, size_t right )
{
    //-- TODO : Calculate the lnLikelihood for rooted trees. Should be able to copy the code for the
    //          unroooted case, but must change the child indexes/partials buffers.
}


template<class charType>
void This::computeRootLikelihood ( size_t root, size_t left, size_t right, size_t middle )
{
    size_t num_taxa  = (this->num_nodes + 2) / 2;

    size_t root_idx  = root   + this->num_nodes * this->activeLikelihood[root];
    //this->b_node_indices.push_back(root_idx); //-- TESTING!

    size_t mid_idx   = middle + this->num_nodes * this->activeLikelihood[middle];
    size_t left_idx  = left   + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right  + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t mid_partials   = (middle < num_taxa) ? middle : mid_idx;
    size_t left_partials  = (left   < num_taxa) ? left   : left_idx;
    size_t right_partials = (right  < num_taxa) ? right  : right_idx;

    //-- Push the last operation onto the queue
    BeagleOperation b_operation =
        { .destinationPartials    = (int) root_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };
    this->b_ops.push_back(b_operation);

    int     b_parentBufferIndices     = (int) root_idx;
    int     b_childBufferIndices      = (int) mid_partials;
    int     b_probabilityIndices      = (int) mid_idx;
    int*    b_firstDerivativeIndices  = NULL;
    int*    b_secondDerivativeIndices = NULL;
    int     b_cumulativeScaleIndices  = BEAGLE_OP_NONE;
    int     b_count                   = 1;
    double  b_outSumLogLikelihood     = 0;
    double* b_outSumFirstDerivative   = NULL;
    double* b_outSumSecondDerivative  = NULL;

    int b_code_update_transitions;
    int b_code_update_partials;
    int b_code_calc_edges;

    //-- these are just for testing... need to figure out what they actually are.
    int stateFrequencyIndex  = 0;  //-- should just be same as model index?
    int categoryWeightsIndex = 0;  //-- TODO : figure out what this is...

    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    #endif /* RB_BEAGLE_DEBUG */

    this->ln_beagle_probability = 0.0;
    for ( auto model : this->b_model_indices )
    {
        #if defined ( RB_BEAGLE_DEBUG )
        ss << "Calculating model : " << model << "\n";
        #endif /* RB_BEAGLE_DEBUG */

        //-- Update all transition matrices.
        b_code_update_transitions =
            beagleUpdateTransitionMatrices( this->beagle_instance
                                          , (int) model
                                          , &this->b_node_indices[0]
                                          , NULL
                                          , NULL
                                          , &this->b_branch_lengths[0]
                                          , this->b_branch_lengths.size()
                                          );
        if ( b_code_update_transitions != 0 )
        {
            throw RbException("Could not update transition matrix for model '" +
                    std::to_string(model) + "'. " +
                    this->parseReturnCode(b_code_update_transitions));
        }

        //-- Calculate and update all partial likelihood buffers
        b_code_update_partials =
            beagleUpdatePartials( this->beagle_instance
                                , &this->b_ops[0]
                                , this->b_ops.size()
                                , BEAGLE_OP_NONE
                                );
        if ( b_code_update_partials != 0 )
        {
          throw RbException("Could not update partials for model '" + std::to_string(model) + "'. "
                    + this->parseReturnCode(b_code_update_partials));
        }

        //-- Calclulate the lnLikelihood of the model
        b_code_calc_edges =
            beagleCalculateEdgeLogLikelihoods( this->beagle_instance
                                             , &b_parentBufferIndices
                                             , &b_childBufferIndices
                                             , &b_probabilityIndices
                                             , b_firstDerivativeIndices
                                             , b_secondDerivativeIndices

                                             //, &this->b_inCategoryRates[0]
                                             //, &this->b_inCategoryWeights[0]
                                             , &categoryWeightsIndex

                                             //, &this->b_stateFrequenciesIndex
                                             , &stateFrequencyIndex
                                             //, (int) model

                                             , &b_cumulativeScaleIndices
                                             , b_count
                                             , &b_outSumLogLikelihood
                                             , b_outSumFirstDerivative
                                             , b_outSumSecondDerivative
                                             );
        if ( b_code_calc_edges != 0 )
        {
            throw RbException("Could not calculate edge log likelihood for model '"
                    + std::to_string(model) + "'. " + this->parseReturnCode(b_code_calc_edges));
        }

        this->ln_beagle_probability += b_outSumLogLikelihood;
    }

    #if defined ( RB_BEAGLE_DEBUG )
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    //-- Reset the beagle operations queues
    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();
}



template<class charType>
void This::computeInternalNodeLikelihood ( const TopologyNode &node
                                         , size_t node_index
                                         , size_t left
                                         , size_t right
                                         )
{
    double branch_length = calculateBranchLength(node, node_index);

    size_t node_idx  = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t left_idx  = left       + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right      + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials  = (node.getChild(0).isTip()) ? left  : left_idx;
    size_t right_partials = (node.getChild(1).isTip()) ? right : right_idx;

    BeagleOperation b_operation =
        { .destinationPartials    = (int) node_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(node_idx);
    this->b_branch_lengths.push_back(branch_length);
}



template<class charType>
void This::computeTipLikelihood ( const TopologyNode &node, size_t node_index )
{
    double branch_length = calculateBranchLength(node, node_index);
    this->b_branch_lengths.push_back(branch_length);

    size_t node_idx = node_index + this->num_nodes * this->activeLikelihood[node_index];
    this->b_node_indices.push_back(node_idx);
}

//--}

//--------------------------------------------------------------------------------------[ Execute ]
//--{1

//-- TODO : See if these are actually needed here or we can use from the Abstract class...

template<class charType>
void This::executeMethod( const std::string &n
                        , const std::vector<const DagNode*> &args
                        , MatrixReal &rv) const
{
    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << "EXEC MAT";
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    if ( n == "siteLikelihoods" )
    {
        // TODO: maybe move model updaters here...

        // make sure the likelihoods are updated
        const_cast<PhyloCTMCSiteHomogeneousBEAGLE<charType>*>(this)->computeLnProbability();
    }
    else if ( n == "siteMixtureLikelihoods" )
    {
        // TODO: maybe move model updaters here...

        // make sure the likelihoods are updated
        const_cast<PhyloCTMCSiteHomogeneousBEAGLE<charType>*>(this)->computeLnProbability();
    }
    else
    {
        throw RbException("The PhyloCTMC process does not have a member method called '" + n + "'.");
    }
}



template<class charType>
void This::executeMethod( const std::string &n
                        , const std::vector<const DagNode*> &args
                        , RbVector<double> &rv) const
{
    #if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << "EXEC VECT";
    RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    if ( n == "siteLikelihoods" )
    {
        // TODO: maybe move model updaters here...

        // make sure the likelihoods are updated
        const_cast<PhyloCTMCSiteHomogeneousBEAGLE<charType>*>(this)->computeLnProbability();
    }
    else if ( n == "siteRates" )
    {
        // TODO: maybe move model updaters here...

        // make sure the likelihoods are updated
        const_cast<PhyloCTMCSiteHomogeneousBEAGLE<charType>*>(this)->computeLnProbability();
    }
    else
    {
        throw RbException("The PhyloCTMC process does not have a member method called '" + n + "'.");
    }
}


//--}

//--------------------------------------------------------------------------------------[ Cleanup ]
//--{1

//-- Undefine the local namespace shortcut
#undef This

#endif

//--}
