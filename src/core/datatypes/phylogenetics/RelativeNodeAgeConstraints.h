#ifndef RelativeNodeAgeConstraints_H
#define RelativeNodeAgeConstraints_H

#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>
#include "RbFileManager.h"

#include "Cloneable.h"

namespace RevBayesCore {
class RelativeNodeAgeConstraintsReader;

    class RelativeNodeAgeConstraints : public Cloneable {
        
    public:
        RelativeNodeAgeConstraints();
        RelativeNodeAgeConstraints(RelativeNodeAgeConstraintsReader* tadr);
        RelativeNodeAgeConstraints(const RelativeNodeAgeConstraints& a);
        RelativeNodeAgeConstraints&                     operator=(const RelativeNodeAgeConstraints& a);
        virtual RelativeNodeAgeConstraints*             clone(void) const;
        
        const std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >& getConstraint(size_t i) const;
        const std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > >& getConstraints( ) const;

        size_t                                          getNumberOfConstraints(void) const;
        const path&                                     getFilename(void) const;
       // std::string                                     getDatatype(void) const;
    protected:
        std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > >								        olderYoungerConstraints;
        
    private:
        path                                            filename;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const RelativeNodeAgeConstraints& x);
}


#endif
