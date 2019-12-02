#ifndef RlUserInterfaceOutputStream_h
#define RlUserInterfaceOutputStream_h

#include <string>

class UserInterfaceOutputStream {
    
public:
    
    UserInterfaceOutputStream(void) {}
    virtual ~UserInterfaceOutputStream(void) {}
    
    virtual void        output(const std::string &o) const = 0;
    virtual void        outputEndOfLine(void) const = 0;

};


#endif /* RlUserInterfaceOutputStream_h */
