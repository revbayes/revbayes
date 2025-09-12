#ifndef AbstractTrace_H
#define AbstractTrace_H

#include <cstddef>
#include <optional>

#include "Cloneable.h"

namespace RevBayesCore {

    class AbstractTrace : public Cloneable {

    public:

        virtual                         ~AbstractTrace(void) {}

        virtual void                    addValueFromString(const std::string &s) = 0;
        virtual int                     isCoveredInInterval(const std::string &v, double i, bool verbose, std::optional<bool> stochastic = std::nullopt) = 0;
    };

}

#endif
