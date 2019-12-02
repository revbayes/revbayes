#ifndef UTILS_H
#define	UTILS_H

#include <string>
#include <cctype>
#include <string>

#include <boost/foreach.hpp>

typedef std::vector<std::string> StringVector;

class LineEditUtils {
public:

    struct FindLastInfo {
        int position;
        std::string needle;
        //std::string nonBlankNeedle;
    };
    
    /**
     * Find the last position of the longest needle possible in haystack
     * 
     * @param needles
     * @param haystack
     * @return 
     */

    static FindLastInfo findLastOf(StringVector needles, std::string haystack, int startFrom) {
        
        std::string lastMatch = "";
        //std::string lastNonBlankMatch = "";
        //int specialMatchType = 0; // 0 = no match, 1 = reference assignment operator, 2 = other
        std::string stringBuf = haystack.substr(startFrom);
        size_t rPos = std::string::npos, _rPos, tmpInt;

        BOOST_FOREACH(std::string s, needles) {

            _rPos = stringBuf.find(s, rPos); // first match
            tmpInt = _rPos;
            while (tmpInt != std::string::npos) { // find last match
                _rPos = tmpInt + s.length();
                tmpInt = stringBuf.find(s, tmpInt + s.length());
            }

            if (_rPos > rPos) {
//                std::cout << "_rPos = '" << _rPos << "' ";
//                if (s == "<-&") {
//                    specialMatchType = 1;
//                    std::cout << ", extra special found ";
//                } else if (specialMatchType == 1 && s != " ") {
//                    specialMatchType = 2;
//                    std::cout << ", trivial extra found ";
//                }
//                if (s == " "){
//                    if (lastMatch != " "){
//                        lastNonBlankMatch = lastMatch;
//                        lastMatch = " ";
//                    }
//                }else{
//                    lastMatch = s;
//                    lastNonBlankMatch = s;
//                }
                lastMatch = s;
                rPos = _rPos;
                //std::cout << "'\n\r";
            }
        }
        //startPos += rPos;
        FindLastInfo f;
        f.needle = lastMatch;
        f.position = (int)rPos;
        return f;
    }
    
    
    static std::string trim(const std::string &s) {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && isspace(*it))
            it++;

        std::string::const_reverse_iterator rit = s.rbegin();
        while (rit.base() != it && isspace(*rit))
            rit++;

        return std::string(it, rit.base());
    }

    static std::string reverseString(std::string s) {
        std::string result = "";
        for (size_t i = 0; i < s.size(); i++) {
            result = s.at(i) + result;
        }
        return result;
    }

    static bool contains(char c, std::string buf) {
        std::string s(buf);
        for (unsigned int i = 0; i < s.size(); i++) {
            if (s.at(i) == c) {
                return true;
            }
        }
        return false;
    }

    static bool lastCharContains(std::string buf, std::string test) {
        if (buf.size() <= 0) {
            return false;
        }
        return contains(buf.at(buf.size() - 1), test);
    }

    static std::string extractSubject(std::string cmd) {
        std::string s = "";
        bool foundStart = false;

        for (int i = (int)cmd.size() - 1; i >= 0; i--) {
            char c = cmd.at(i);
            if (!foundStart && isalpha(c)) {
                foundStart = true;
            }

            if (foundStart) {
                if (isalnum(c)) {
                    s += c;
                } else {
                    break;
                }
            }
        }
        return trim(reverseString(s));
    }


};
#endif	/* UTILS_H */

