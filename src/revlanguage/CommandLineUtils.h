#ifndef COMMANDLINEUTILS_H
#define	COMMANDLINEUTILS_H


#include <string>
#include <cstdint>
std::string nonblankSeparators = "+->~()[]{}.,\"";
std::string separators = nonblankSeparators + " ";
std::string openingSeparators = "\"(,";

class CommandLineUtils {
public:

    std::string getFunctionName(const char* buf) {
        std::string s = "";
        for (size_t i = getLastSeparatorPosition(buf, true); i < strlen(buf); i++) {
            if (!contains(buf[i], separators.c_str())) {
                s += buf[i];
            } else {
                return s;
            }
        }
        return s;
    }
    
    bool endsWithSeparator(const char* buf){
        return contains(buf[strlen(buf)-1], nonblankSeparators.c_str());
    }

    // todo: implement regex instead

    size_t getLastSeparatorPosition(const char* buf, bool discardTrailingSeparators = false) {
        if (strlen(buf) == 0) {
            return 0;
        }
        size_t pos = 0;

        // search from end of string to previous separator
        bool foundNonSeparator = false;
        for (size_t i = strlen(buf) - 1; i > 0; i--) {
            if (discardTrailingSeparators) {
                if (!foundNonSeparator) {
                    foundNonSeparator = !contains(buf[i], separators.c_str());
                } else if (contains(buf[i], separators.c_str())) { // start checking if char is a separator first after a non separator is found
                    pos = i +1;
                    break;
                }
            } else if (contains(buf[i], separators.c_str())) {
                pos = i+1;
                break;
            }
        }
        return pos;
    }
    
    bool lastEffectiveSeparatorIsBlank(const char* buf) {
        std::int64_t i = std::int64_t(strlen(buf))-1;
        if (i < 0){
            return false;
        }
        
        while(buf[i] == ' ') i--;        
        return !contains(buf[i], openingSeparators.c_str());
    }

    bool contains(char c, const char* buf) {
        for (unsigned int i = 0; i < strlen(buf); i++) {
            if (buf[i] == c) {
                return true;
            }
        }
        return false;
    }

private:

};

#endif	/* COMMANDLINEUTILS_H */

