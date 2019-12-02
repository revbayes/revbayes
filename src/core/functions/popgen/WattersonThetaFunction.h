#ifndef WattersonThetaFunction_H
#define WattersonThetaFunction_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Computing Watterson's theta function.
     *
     * This function computes Watterson's theta, which is an estimator of theta = 4*N*mu.
     * This function can be used to compute the per site or per locus estimator.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-04-30
     *
     */
    class WattersonThetaFunction : public TypedFunction<double> {
        
    public:
        WattersonThetaFunction(const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a, bool ps, bool excl);
        virtual                                                            ~WattersonThetaFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        WattersonThetaFunction*                                             clone(void) const;                                                              //!< Create an independent clone
        void                                                                update(void);
        
    protected:
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >*      alignment;
        bool                                                                per_site;
        bool                                                                exclude_ambiguous;
    };
    
}

#endif
