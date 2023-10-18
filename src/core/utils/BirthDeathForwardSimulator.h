#ifndef BirthDeathForwardSimulator_H
#define BirthDeathForwardSimulator_H

#include <cstddef>
#include <vector>

namespace RevBayesCore {

    class Tree;
    class TopologyNode;

    /**
     * General case birth-death forward simulator.
     *
     * All parameters are vectors of vectors. The outer vector is for the epochs. The inner vector is for the state specific rates.
     * For example, lambda[2][1] is the speciation rate in the 3rd epoch for a lineage in state 2.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class BirthDeathForwardSimulator {
        
    public:
        
        enum SIM_CONDITION { TIME, SURVIVAL, ROOT };

        BirthDeathForwardSimulator();
        
        // setters
        void                                    setBurstProbability( const std::vector<std::vector< double > > &l );
        void                                    setCompleteTree( bool c );
        void                                    setExtinctionRate( const std::vector<std::vector< double > > &m );
        void                                    setMaxNumLineages( size_t m );
        void                                    setMassExtinctionProbability( const std::vector<std::vector< double > > &m );
        void                                    setRootCategoryProbabilities( const std::vector<double> &p );
        void                                    setSamplingProbability( const std::vector<std::vector< double > > &p );
        void                                    setSamplingExtinctionProbability( const std::vector<std::vector< double > > &r );
        void                                    setSamplingRate( const std::vector<std::vector< double > > &p );
        void                                    setSamplingExtinctionRate( const std::vector<std::vector< double > > &r );
        void                                    setSpeciationRate( const std::vector<std::vector< double > > &l );
        void                                    setTimeline( const std::vector< double > &t );

//        Tree*                                   simulateTreeConditionTaxa( size_t n ) const;
        Tree*                                   simulateTreeConditionTime( double t, SIM_CONDITION cdt ) const;

    private:
        
        bool                                    checkParameters(void) const;
        size_t                                  getNumberOfCategories( void ) const;
        double                                  getLambdaProbability( size_t i, size_t n ) const;
        std::vector<double>                     getLambdaRate( size_t i, size_t n ) const;
        double                                  getMuProbability( size_t i, size_t n ) const;
        std::vector<double>                     getMuRate( size_t i, size_t n ) const;
        double                                  getPhiProbability( size_t i, size_t n ) const;
        std::vector<double>                     getPhiRate( size_t i, size_t n ) const;
        double                                  getRProbability( size_t i, size_t n ) const;
        std::vector<double>                     getRRate( size_t i, size_t n ) const;
        std::vector<double>                     getRootCategoryProbabilities( size_t n ) const;
        bool                                    hasExtantSurvivor(const TopologyNode &n) const;

        std::vector<std::vector< double > >     lambda;
        std::vector<std::vector< double > >     Lambda;
        std::vector<std::vector< double > >     mu;
        std::vector<std::vector< double > >     Mu;
        std::vector<std::vector< double > >     phi;
        std::vector<std::vector< double > >     Phi;
        std::vector<std::vector< double > >     r;
        std::vector<std::vector< double > >     R;
        std::vector<double>                     root_cat_probability;
        std::vector<double>                     timeline;
        size_t                                  MAX_NUM_LINEAGES;
        bool                                    complete_tree;        
    };
    
}


#endif
