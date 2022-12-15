#ifndef REVCLIENT_H
#define	REVCLIENT_H

#include <set>
#include <string>
#include <vector>


//struct tabCompletionInfo{
//        int startPos;
//        unsigned int specialMatchType;
//        std::vector<std::string> completions;
//        std::vector<std::string> matchingCompletions;
//        std::string compMatch;
//    };

namespace RevClient
{

    void  startInterpreter(void);
    int   interpret(const std::string& command);
}
    

//    std::set<std::string>       getDefaultCompletions();
//    std::vector<std::string>    getFileList(const std::string &path);
//        tabCompletionInfo           getTabCompleteInfo(const char *buf);
//        void                        setTabCompletionInfo(const char *buf);
    

#endif	/* REVCLIENT_H */

