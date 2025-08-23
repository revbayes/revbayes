#ifndef PseudoDist_H
#define PseudoDist_H

#include "PseudoData.h"
#include "PseudoObservation.h"

/*
 * This distribution handles situations where:
 * - we have data o that affects some parameter x
 * - the data o is unspecified
 * - the distribution is unspecified
 * - the likelihood function Pr(o|x, a, b, ...) is specified as a function l(x, a, b,...).
 *
 * In RevBayes this would look like:
 *
 *   o ~ dnPseudo(x, pdL(a, b, ...))
 *   o.clamp( pseudoObservation() )
 *
 * This only makes sense if o is clamped, so if the user just writes:
 *
 *   o ~ dnPseudo(x, pdL(a, b, ...))
 *
 * then we set the PseudoObservation pointer to nullptr.  If the pointer is nullptr, then
 * we return a probability of 1.
 *
 * Clamping this should add a value, which triggers returning the likelihood.
 */

namespace RevBayesCore
{
    template <typename T>
    class PseudoDist : public TypedDistribution< PseudoObservation >
    {
        // the parameter the pseudodata applies to
        const TypedDagNode<T>*                              parameter = nullptr;
        // the pseudodata
        const TypedDagNode<PseudoData<T>>*                  pseudo_data = nullptr;

    protected:

        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            if (oldP == parameter)
            {
                this->parameter = static_cast<const TypedDagNode<T>*>(newP);
            }
            else if (oldP == pseudo_data)
            {
                this->pseudo_data = static_cast<const TypedDagNode<PseudoData<T>>*>(newP);
            }
        }

    public:
        PseudoDist(const TypedDagNode<T>* p, const TypedDagNode<PseudoData<T>>* pd)
            :TypedDistribution<PseudoObservation>(nullptr),
             parameter(p),
             pseudo_data(pd)
        {
            this->addParameter(parameter);
            this->addParameter(pseudo_data);
        }

        PseudoDist(const PseudoDist<T>& d) = default;

        PseudoDist*                                         clone(void) const override {return new PseudoDist(*this);}

        double                                              computeLnProbability(void) override
        {
            if (this->value)
            {
                // We return the likelihood Pr(unspecified data | x)
                auto& x = parameter->getValue();
                auto& f = pseudo_data->getValue();
                return f(x);
            }
            else
                return 0;
        }

        void                                                redrawValue(void) override
        {
            delete this->value;
            this->value = nullptr;
        }
    };
}




#endif
