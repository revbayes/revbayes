#ifndef DPPAllocateAuxGibbsMove_H
#define DPPAllocateAuxGibbsMove_H

#include "DirichletProcessPriorDistribution.h"
#include "AbstractGibbsMove.h"
#include "RbMathLogic.h"

namespace RevBayesCore {

    /**
     * @brief The DDP Move for updating the allocation structure of DPP-distributed elements *templated class
     *
     * This move follows the algorithm 8 in Neal (2000, Markov Chain Sampling Methods for Dirichlet Process Mixture Models)
     * It is a Gibbs sampling proposal with auxiliary parameters.
     * The move operates on:
     *  @param v a vector of elements
     *
     * The tuning parameter of the move is essentially determined by
     *  @param na the number of auxiliary parameters (i.e., auxiliary tables)
     *
     */

    template <class valueType>
    class DPPAllocateAuxGibbsMove : public AbstractGibbsMove {
    
    public:

        DPPAllocateAuxGibbsMove(StochasticNode< RbVector<valueType> >* v, int na, double w);                                      //!< Internal constructor
    
        // Basic utility functions
        DPPAllocateAuxGibbsMove<valueType>*                     clone(void) const;                                                //!< Clone object
        const std::string&                                      getMoveName(void) const;                                          //!< Get the name of the move for summary printing
        
    protected:

        void                                                    performGibbsMove(void);                                           //!< Perform move
        void                                                    swapNodeInternal(DagNode *oldN, DagNode *newN);
        
    private:
        
        double                                                  getLnProbabilityForMove(void);
        int                                                     findTableIDForVal(std::vector<valueType> tvs, valueType val);
        void                                                    normalizeVector(std::vector<double> &v);
        int                                                     findElementNewTable(double u, const std::vector<double> &prob);
        
        StochasticNode< RbVector<valueType> >*                  variable;
        int                                                     num_aux_cat;
    };
}




template <class valueType>
RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::DPPAllocateAuxGibbsMove(StochasticNode< RbVector<valueType> > *v, int na, double w) : AbstractGibbsMove( w ),
    variable( v )
{
    addNode( variable );
    
    num_aux_cat = na;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class valueType>
RevBayesCore::DPPAllocateAuxGibbsMove<valueType>* RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::clone( void ) const
{
    return new DPPAllocateAuxGibbsMove<valueType>( *this );
}


template <class valueType>
const std::string& RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::getMoveName( void ) const
{
    static std::string name = "DPPAllocateAuxGibbsMove";
    
    return name;
}


/** Perform the move */
template <class valueType>
void RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::performGibbsMove( void )
{
    
    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    DirichletProcessPriorDistribution<valueType>& dist = static_cast<DirichletProcessPriorDistribution<valueType> &>( variable->getDistribution() );
    
    dist.createRestaurantVectors();
    int num_tables = dist.getNumberOfCategories();
    int num_elements = dist.getNumberOfElements();
    std::vector<valueType> table_values = dist.getTableParameters();
    std::vector<int> allocation_vector = dist.getElementAllocationVec();
    std::vector<int> num_allocated_to_table = dist.getNumElemPerTable();
    std::vector<valueType>& element_values = variable->getValue();
    double alpha = dist.getConcentrationParam();
    TypedDistribution<valueType>* base_distribution = dist.getBaseDistribution();
        
    double ln_alpha_over_num_aux = log( alpha/( (double)num_aux_cat ) );

    // loop over elements
    for (size_t i=0; i < num_elements; ++i)
    {
        
        // vector to hold lnL after moving this element to a different table
        std::vector<double> ln_probs;
        std::vector<valueType> temp_tables;
        
        // remove element i from current table
        int current_table = allocation_vector[i];
        num_allocated_to_table[current_table] -= 1;

        // try element i at other existing tables
        for (size_t j = 0; j < num_tables; ++j)
        {
            int num_seated = num_allocated_to_table[j];
            if (num_seated > 0)
            {
                valueType new_value = table_values[j];
                temp_tables.push_back(new_value);
                element_values[i] = new_value;
                variable->touch();
                double table_lnL = getLnProbabilityForMove(); // get lnL after changing tables
                ln_probs.push_back( log(num_seated) + table_lnL );
            }
            else
            {
                table_values[j] = valueType();
                num_allocated_to_table[j] = 0;
            }
        }
        
        // try element i at other auxillary (new) tables
        for (size_t j = 0; j < num_aux_cat; ++j)
        {
            base_distribution->redrawValue();
            valueType new_value = base_distribution->getValue();
            temp_tables.push_back(new_value);
            element_values[i] = new_value;
            variable->touch();
            double aux_lnL = getLnProbabilityForMove(); // get lnL after changing tables
            ln_probs.push_back( ln_alpha_over_num_aux + aux_lnL );
        }
        
        // normalize ln_probs vector
        normalizeVector(ln_probs);
                
        // choose new table for element i
        double u = rng->uniform01();
        int new_table = findElementNewTable(u, ln_probs);
        if (new_table != -1)
        {
            const valueType &final_value = temp_tables[new_table];
            element_values[i] = final_value;
            variable->touch();
            int tID = findTableIDForVal(table_values, final_value);
            if (tID == -1)
            {
                // add i to a new table
                num_tables += 1;
                num_allocated_to_table.push_back(1);
                table_values.push_back(final_value);
            }
            else
            {
                // add i to an existing table
                num_allocated_to_table[tID] += 1;
            }
        }
        variable->keep();
        ln_probs.clear();
        temp_tables.clear();
    }
    
    allocation_vector.clear();
    table_values.clear();
    num_allocated_to_table.clear();
    dist.createRestaurantVectors();
    variable->touch();
    variable->keep();

}


template <class valueType>
void RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    variable = static_cast<StochasticNode< RbVector<valueType> > * >( newN );

}


template <class valueType>
double RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::getLnProbabilityForMove(void)
{
    
    RbOrderedSet<DagNode*> affected;
    variable->initiateGetAffectedNodes( affected );
    double ln_probs = 0.0;
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        DagNode *the_node = *it;
        double lp = the_node->getLnProbability();
        ln_probs += lp;
    }
    return ln_probs;
}


template <class valueType>
int RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::findTableIDForVal(std::vector<valueType> tvs, valueType val) {
    
    for (size_t j=0; j<tvs.size(); j++)
    {
        if (tvs[j] == val)
        {
            // this is an existing table
            return int(j);
        }
    }
    
    return -1; // new table
}


template <class valueType>
void RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::normalizeVector(std::vector<double> &v)
{

    size_t n = v.size();
    double max = v[0];
    for (size_t i = 1; i < n; ++i)
    {
        if (v[i] > max)
        {
            max = v[i];
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        v[i] -= max;
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        if ( v[i] < -300.0 || RbMath::isNan(v[i]) )
        {
            v[i] = 0.0;
        }
        else
        {
            v[i] = exp( v[i] );
        }
        sum += v[i];
    }

    if (sum == 0.0)
    {
        for (size_t i=0; i<n; ++i)
        {
            v[i] /= sum;
        }
    }    
}


template <class valueType>
int RevBayesCore::DPPAllocateAuxGibbsMove<valueType>::findElementNewTable(double u, const std::vector<double> &prob)
{
    
    for (size_t j = 0; j < prob.size(); ++j)
    {
        u -= prob[j];
        if (u < 0.0)
        {
            return int(j);
        }
    }
   
    // could not find new table for element
    return -1;
    
}


#endif
