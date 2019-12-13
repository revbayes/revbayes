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
