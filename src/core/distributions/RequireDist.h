#ifndef RequireDist_H
#define RequireDist_H

#include "PseudoObservation.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"
#include "RbBoolean.h"

namespace RevBayesCore
{
    class RequireDist : public TypedDistribution< PseudoObservation >
    {
        // the parameter the pseudodata applies to
        const TypedDagNode<Boolean>*                              predicate = nullptr;
        const TypedDagNode<double>*                            weight = nullptr;

    protected:

        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            if (oldP == predicate)
            {
                this->predicate = static_cast<const TypedDagNode<Boolean>*>(newP);
            }
            else if (oldP == weight)
            {
                this->weight = static_cast<const TypedDagNode<double>*>(newP);
            }
        }

    public:
        RequireDist(const TypedDagNode<Boolean>* p, const TypedDagNode<double>* w)
            :TypedDistribution<PseudoObservation>(nullptr),
             predicate(p),
             weight(w)
        {
            this->addParameter(predicate);
            this->addParameter(weight);
        }

        RequireDist(const RequireDist& d) = default;

        RequireDist*                                         clone(void) const override {return new RequireDist(*this);}

        double                                              computeLnProbability(void) override
        {
            // Set the probability to 1 if not clamped.
            if (not this->value) return 1;

            if (predicate->getValue())
                return 0;
            else
                return -log(weight->getValue());
        }

        void                                                redrawValue(void) override
        {
            delete this->value;
            this->value = nullptr;
        }
    };
}




#endif
