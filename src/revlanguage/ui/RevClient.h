#ifndef REVCLIENT_H
#define	REVCLIENT_H

#include <set>
#include <string>
#include <vector>
#include <filesystem>


//struct tabCompletionInfo{
//        int startPos;
//        unsigned int specialMatchType;
//        std::vector<std::string> completions;
//        std::vector<std::string> matchingCompletions;
//        std::string compMatch;
//    };

namespace RevClient
{
    int   interpret(const std::string& command);

    void  execute_file(const std::filesystem::path& filename);
    void  shutdown();
    void  startInterpreter();
    void  startJupyterInterpreter();
}
    

//    std::set<std::string>       getDefaultCompletions();
//    std::vector<std::string>    getFileList(const std::string &path);
//        tabCompletionInfo           getTabCompleteInfo(const char *buf);
//        void                        setTabCompletionInfo(const char *buf);
    

#endif	/* REVCLIENT_H */

