#ifndef AbstractTrace_H
#define AbstractTrace_H

#include <stddef.h>

#include "Cloneable.h"

namespace RevBayesCore {

    class AbstractTrace : public Cloneable {

    public:

        virtual                         ~AbstractTrace(void) {}

        virtual void                    addValueFromString(const std::string &s) = 0;
        virtual int                     isCoveredInInterval(const std::string &v, double i, bool verbose) = 0;

    };

}

#endif
