#ifndef RlCommandLineOutputStream_h
#define RlCommandLineOutputStream_h

#include <iosfwd>

#include "RlUserInterfaceOutputStream.h"

class CommandLineOutputStream : public UserInterfaceOutputStream {
    
public:
    
    CommandLineOutputStream(void) {}
    virtual ~CommandLineOutputStream(void) {}
    
    void        output(const std::string &o) const;
    void        outputEndOfLine(void) const;
    
};


#endif /* RlCommandLineOutputStream_hpp */
