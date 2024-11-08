#include <sys/stat.h>
#include <cstdio>
#include <sys/types.h> // IWYU pragma: keep
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "RbFileManager.h"
#include "RbSettings.h"
#include "RbException.h"
#include "StringUtilities.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

// TODO: remove all these includes
#ifdef _WIN32
#	include <dirent.h>
#   include <unistd.h>
#   include <windows.h>
#   include <shlwapi.h>
#else
#	include <dirent.h>
#   include <unistd.h>
#endif

// TODO: fix setValueFromFile

namespace fs = boost::filesystem;

namespace RevBayesCore
{

void createDirectoryForFile(const path& p)
{
    if (not p.parent_path().empty())
        create_directories(p.parent_path());
}

path appendToStem(const path& p, const std::string& s)
{
    return p.parent_path() / ( p.stem().string() + s + p.extension().string() );
}

/** Get line while safely handling cross-platform line endings.
 *  Modified from: https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
 *
 * @param[in] is stream from which to read
 * @param[out] t string to store the line read
 *
 * @return the input stream
 * */
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}


/**
 * Portable code to get a full path by expanding the user home directory.
 *
 * @param path relative path
 * @return full path including user home directory
 */
path expandUserDir(std::string dir)
{
    if ( !dir.empty() && dir[0] == '~')
    {
        char const* home = getenv("HOME");
        
        if (home or ((home = getenv("USERPROFILE"))))
        {
            dir.replace(0, 1, home);
        }
    }
    else if ( dir.empty() == false )
    {
        char const *hdrive = getenv("HOMEDRIVE"), *hpath = getenv("HOMEPATH");
        if ( hdrive != NULL )
        {
// TODO: remove this ifdef.
# ifdef _WIN32
            dir = std::string(hdrive) + hpath + "\\" + dir;
# else
            // FIXME: there is no leading `~`, so this doesn't make sense.
            dir.replace(0, 1, std::string(hdrive) + hpath);
# endif
        }
        
    }

    return dir;
}


/** Format an error exception string for problems specifying the file/path name
 * @param[out] errorStr string to store the formatted error
*/
void formatError(const path& p, std::string& errorStr)
{
    bool file_nameProvided    = (not p.filename().empty() and not p.filename_is_dot());
    bool isfile_nameGood      = is_regular_file(p);
    bool isDirectoryNameGood = is_directory( p.parent_path() );
    
    if ( file_nameProvided == false && isDirectoryNameGood == false )
    {
        errorStr += "Could not read contents of directory \"" + p.parent_path().string() + "\" because the directory does not exist";
    }
    else if (file_nameProvided == true && (isfile_nameGood == false || isDirectoryNameGood == false))
    {
        errorStr += "Could not read file named \"" + p.filename().string() + "\" in directory named \"" + p.parent_path().string() + "\" ";
        if (isfile_nameGood == false && isDirectoryNameGood == true)
        {
            errorStr += "because the file does not exist";
        }
        else if (isfile_nameGood == true && isDirectoryNameGood == false)
        {
            errorStr += "because the directory does not exist";
        }
        else
        {
            errorStr += "because neither the directory nor the file exist";
        }
    }
}

/** Fills in a vector with the names of the files in a directory
 *
 * @param[in] dirpath path to the directory to be listed
 * @param[out] sv vector to store the names
 * @param[in] recursive whether to list recursively, default true
 *
 * @return true
*/
bool setStringWithNamesOfFilesInDirectory(const path& dirpath, std::vector<path>& sv, bool recursive)
{
    // FIXME: It should be converted to use boost:filesystem.
    //        This is a holdover from the days of RbFileManager
    //        We should try and remove the #ifdef _WIN32, and 
    std::string dirstring = dirpath.string();

    DIR* dir = opendir( dirstring.c_str() );
    if (dir)
    {
        struct dirent* entry;
        while ( (entry = readdir( dir )) )
        {
            struct stat entryinfo;
            std::string entryname = entry->d_name;
            std::string entrypath = (dirpath  / entryname).string();
            
            bool skip_me = false;

#ifdef _WIN32
            if (stat( entrypath.c_str(), &entryinfo )) {
              // if this returned a non-zero value, something is wrong
              skip_me = true;
            }
#else
            if (!lstat( entrypath.c_str(), &entryinfo ))
            {
                
                // avoid recursing into symlinks that point to a directory above us
                if ( S_ISLNK( entryinfo.st_mode ) ) {
                    char *linkname = (char*) malloc(entryinfo.st_size + 1);
                    ssize_t r = readlink(entrypath.c_str(), linkname, entryinfo.st_size + 1);
                    if (r < 0 || r > entryinfo.st_size) {
                        // error
                        skip_me = true;
                    } else {
                        linkname[entryinfo.st_size] = '\0';
                        if (strspn(linkname, "./") == entryinfo.st_size) {
                            // this symlink consists entirely of dots and dashes and is likely a reference to a directory above us
                            skip_me = true;
                        } else {
                            // replace entryinfo with info from stat
                            if ( stat( entrypath.c_str(), &entryinfo ) ) {
                                // if this returned a non-zero value, something is wrong
                                skip_me = true;
                            }
                        }
                    }
                    free(linkname);
                }

            } else {
              // if this returned a non-zero value, something is wrong
              skip_me = true;
            }
#endif

            if (!skip_me) {

                if (entryname == "..")
                {
                    ;
                }
                else if (entryname == "." || entryname[0] == '.')
                {
                    ;
                }
                else if ( recursive == true && S_ISDIR( entryinfo.st_mode ) )
                {
                    setStringWithNamesOfFilesInDirectory( entrypath, sv );
                }
                else
                {
                    sv.push_back( entrypath );
                }
            }
            
        }
        
        closedir( dir );
    }
    
    // make sure that the file names are sorted
    std::sort(sv.begin(), sv.end());
    
    return true;
}

std::stringstream readFileAsStringStream(const path& fname)
{
    std::ifstream inFile(fname.string());
    if ( not inFile )
        throw RbException()<<"Could not open file "<<fname<<".";

    std::stringstream strStream;
    strStream << inFile.rdbuf(); //read the file
    return strStream;
}

std::string readFileAsString(const path& fname)
{
    return readFileAsStringStream(fname).str();
}


}
