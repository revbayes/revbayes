#ifndef UnconstrainedSBN_H
#define UnconstrainedSBN_H

#include "MatrixReal.h"
#include "Subsplit.h"
#include "SBNParameters.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    class UnconstrainedSBN : public TypedDistribution<Tree> {

    public:
        UnconstrainedSBN(void);
        UnconstrainedSBN(const SBNParameters parameters);
		virtual                                            ~UnconstrainedSBN(void);                                                                    //!< Virtual destructor

        // public member functions
        UnconstrainedSBN*                        clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              computeLnProbabilityUnrootedTopologyMarginalize(void);
        double                                              computeLnProbabilityBranchLengths( void );
        double                                              logSumExp(std::vector<double> &x); // Calculate log(sum(exp(x))) accounting for the fact that exp(x) can produce 0s for x << 0
        double                                              logSumExp(double &x1, double &x2); // Calculate log(sum(exp(x))) accounting for the fact that exp(x) can produce 0s for x << 0
        double                                              logSumExpWeights( std::vector<double> &x, std::vector<double> &p ); // As logSumExp(), but with weights for the non-log values
        double                                              logSumExpWeights( double &x1, double &x2, double &p1, double &p2 ); // As logSumExp(), but with weights for the non-log values
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                    //!< Set the current value, e.g. attach an observation (clamp)

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

    private:

        // helper functions
        void                                                simulateTree(void);

        // members
        SBNParameters                                       parameters;  // The parameters for this SBN
        std::vector<Taxon>                                  taxa;
    };

}

#endif
