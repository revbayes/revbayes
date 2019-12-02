#ifndef BirthDeathForwardSimulator_H
#define BirthDeathForwardSimulator_H

#include "Tree.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    
    /**
     * This class represents the writer object of character data objects into files in Fasta format.
     *
     * This class currently has only one functionality,
     * to write character data objects into a file in Fasta format.
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
        void                                    setExtinctionRate( const std::vector<std::vector< double > > &m );
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
        
    };
    
}


#endif
