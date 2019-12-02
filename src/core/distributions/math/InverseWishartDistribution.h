/* 
 * File:   InverseInverseWishartDistribution.h
 * Author: nl
 *
 * Created on 15 juillet 2014, 10:08
 */

#ifndef INVERSEWISHARTDISTRIBUTION_H
#define	INVERSEWISHARTDISTRIBUTION_H


#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
    
    class InverseWishartDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // inverse InverseWishart distribution of parameter sigma0 and df degrees of freedom
        InverseWishartDistribution(const TypedDagNode<MatrixReal> *insigma0, const TypedDagNode<long>* indf);

        // inverse InverseWishart distribution of parameter sigma0 = Diagonal(kappaVector) and df degrees of freedom
        InverseWishartDistribution(const TypedDagNode<RbVector<double> > *inkappaVector, const TypedDagNode<long>* indf);

        // inverse InverseWishart distribution of parameter sigma0 kappa * Identitymatrix and df degrees of freedom
        InverseWishartDistribution(const TypedDagNode<long>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<long>* indf);

        virtual                                            ~InverseWishartDistribution(void) {}
        
        // public member functions

        InverseWishartDistribution*                         clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        int                                                 getDF() const {return (int)df->getValue();}
        
        const TypedDagNode<MatrixReal>*                     getSigma0() const {return sigma0;}
        const TypedDagNode<RbVector<double> >*              getKappaVector() const {return kappaVector;}
        const TypedDagNode<double>*                         getKappa() const {return kappa;}

        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

    private:

        const TypedDagNode<MatrixReal>*                     sigma0;
        const TypedDagNode<RbVector<double> >*              kappaVector;
        const TypedDagNode<double>*                         kappa;
        const TypedDagNode<long>*                            df;
        const TypedDagNode<long>*                            dim;
                
    };
    
}


#endif	/* INVERSEWISHARTDISTRIBUTION_H */

