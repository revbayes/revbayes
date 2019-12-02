/*
 * File:   PhyloOrnsteinUhlenbeckProcess.h
 * Author: Bastien Boussau
 *
 * Created on 04/10/2014
 */

#ifndef PhyloOrnsteinUhlenbeckProcess_H
#define	PhyloOrnsteinUhlenbeckProcess_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class PhyloOrnsteinUhlenbeckProcess : public TypedDistribution< RbVector<double> > {
        
    public:
        // constructor(s)
        PhyloOrnsteinUhlenbeckProcess(const TypedDagNode< Tree > *t);
        
        // public member functions
        PhyloOrnsteinUhlenbeckProcess*                          clone(void) const;                                                                      //!< Create an independent clone
        double                                                  computeLnProbability(void);
        void                                                    redrawValue(void);
        
        void                                                    setMean( const TypedDagNode< double >* x);
        void                                                    setMean( const TypedDagNode< RbVector< double > >* x);
        void                                                    setPhi( const TypedDagNode< double >* x);
        void                                                    setPhi( const TypedDagNode< RbVector< double > >* x);
        void                                                    setSigma( const TypedDagNode< double >* x);
        void                                                    setSigma( const TypedDagNode< RbVector< double > >* x);
        void                                                    setRoot( const TypedDagNode< double >* x);
        
    protected:
        // Parameter management functions
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter

    private:
        // helper methods
        void                                                    simulate();
        double                                                  recursiveLnProb(const TopologyNode& n);
        void                                                    recursiveSimulate(const TopologyNode& n);
        
        // private members
        const TypedDagNode< Tree >*                             tau;
        const TypedDagNode< RbVector< double > >*               heterogeneousMean;
        const TypedDagNode< RbVector< double > >*               heterogeneousPhi;
        const TypedDagNode< RbVector< double > >*               heterogeneousSigma;
        const TypedDagNode< double >*                           homogeneousMean;
        const TypedDagNode< double >*                           homogeneousPhi;
        const TypedDagNode< double >*                           homogeneousSigma;
        const TypedDagNode< double >*                           root;

    };
    
}


#endif	/* PhyloOrnsteinUhlenbeckProcess_H */

