#ifndef PseudoDist_H
#define PseudoDist_H

#include "PseudoDataLikelihood.h"
#include "PseudoObservation.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"

#include <iostream>

/*
 * This distribution handles situations where:
 * - we have data o that supplies evidence about some parameter x
 * - the data o is unspecified
 * - the distribution is unspecified
 * - the likelihood function Pr(o|x, a, b, ...) is specified as a function l(x, a, b,...).
 *
 * In RevBayes this would look like:
 *
 *   o ~ dnPseudo(pdL(x, a, b, ...))
 *   o.clamp( pseudoObservation() )
 *
 * This only makes sense if o is clamped, so if the user just writes:
 *
 *   o ~ dnPseudo(pdL(x, a, b, ...))
 *
 * then we set the PseudoObservation pointer to nullptr.  If the pointer is nullptr, then
 * we return a probability of 1.
 *
 * Clamping this should add a value, which triggers returning the likelihood.
 */

namespace RevBayesCore
{
    class PseudoDist : public TypedDistribution< PseudoObservation >
    {
        // How many not-clamped warnings have we printed?
        int nWarnings = 0;
        // the pseudo-data likelihood
        const TypedDagNode<PseudoDataLikelihood>*                     pseudo_data_likelihood = nullptr;

    protected:

        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            if (oldP == pseudo_data_likelihood)
            {
                this->pseudo_data_likelihood = static_cast<const TypedDagNode<PseudoDataLikelihood>*>(newP);
            }
        }

    public:
        PseudoDist(const TypedDagNode<PseudoDataLikelihood>* pdl)
            :TypedDistribution<PseudoObservation>(nullptr),
             pseudo_data_likelihood(pdl)
        {
            this->addParameter(pseudo_data_likelihood);
        }

        PseudoDist(const PseudoDist& d) = default;

        PseudoDist*                                         clone(void) const override {return new PseudoDist(*this);}

        double                                              computeLnProbability(void) override
        {
            if (this->value)
            {
                // We return the likelihood Pr(unspecified data | x , parameters)
                PseudoDataLikelihood L = pseudo_data_likelihood->getValue();
                return L.log();
            }
            else
            {
                if (nWarnings < 10)
                {
                    nWarnings++;
                    std::cerr<<"Variable sampled from dnPseudo is not clamped!  It will have no effect.\n";
                    // Sadly, we can't find the name of the variable to clamp from here.
                }
                return 0;
            }
        }

        void                                                redrawValue(void) override
        {
            delete this->value;
            this->value = nullptr;
        }
    };
}




#endif
