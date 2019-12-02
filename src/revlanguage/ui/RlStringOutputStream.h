#ifndef RlStringOutputStream_h
#define RlStringOutputStream_h

#include "RlUserInterfaceOutputStream.h"

#include <stdio.h>
#include <string>

class StringOutputStream : public UserInterfaceOutputStream {
    
public:
    
    StringOutputStream(void) {}
    virtual ~StringOutputStream(void) {}
    
    void                    emptyOutputString(void);
    const std::string&      getOutputString(void) const;
    void                    output(const std::string &o) const;
    void                    outputEndOfLine(void) const;
    
private:
    
    mutable std::string     current_output;
    
};


#endif /* RlStringOutputStream_hpp */
