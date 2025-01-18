#ifndef PseudoDist_H
#define PseudoDist_H

#include "PseudoData.h"

/*
 * This distribution handle situations where:
 * - the data is unspecified
 * - the distribution is unspecified
 * - the likelihood function Pr(data|parameter) is specified as a function f(parameter).
 *
 * In RevBayes this would look like:
 *
 *   x ~ dnPseudo(parameter)
 *   x.clamp(f)
 *
 * The likelihood Pr(x|parameter) is then specified to be f(parameter).
 *
 * This only makes sense if x is clamped, so if the user just writes:
 *
 *   x ~ dnPseudo(parameter)
 *
 * then we set x=E where E is the empty pseudodata.  Then Pr(x) = E(x) = 1.
 */

namespace RevBayesCore
{

    template <typename T>
    class PseudoDist : public TypedDistribution< PseudoData<T>>
    {
        // the base distribution
        const TypedDagNode<T>*                              parameter = nullptr;

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            if (oldP == parameter)
            {
                this->parameter = static_cast<const TypedDagNode<T>*>(newP);
            }
        }

    public:
        PseudoDist(const TypedDagNode<T>* d)
            :TypedDistribution<PseudoData<T>>(nullptr),
             parameter(d)
        {
            this->addParameter(d);

            redrawValue();
        }

        PseudoDist(const PseudoDist<T>& d) = default;

        PseudoDist*                                         clone(void) const override {return new PseudoDist(*this);}

        double                                              computeLnProbability(void) override
        {
            // We return the likelihood Pr(unspecified data | x)
            auto& pseudo_data = this->getValue();
            auto& x = parameter->getValue();
            return pseudo_data(x);
        }

        void                                                redrawValue(void) override
        {
            delete this->value;

            // Supply a PseudoData object with no observations, leading to Pr(unspecified data|parameter) = 1.
            // This will be overridden only if the variable simulated from this distribution is clamped.
            this->value = new PseudoDataNothing<T>();
        }

    };
}




#endif
