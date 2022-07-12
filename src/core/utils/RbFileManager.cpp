#include <sys/stat.h>
#include <stdio.h>
#include <sys/types.h> // IWYU pragma: keep
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "RbFileManager.h"
#include "RbSettings.h"
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
    bool file_nameProvided    = p.filename().string() != "." and p.filename().string() != "" and p.filename().string() != "..";
    bool isfile_nameGood      = exists(p);
    bool isDirectoryNameGood = is_directory( p.parent_path() );

    if ( file_nameProvided == false && isDirectoryNameGood == false )
    {
        errorStr += "Could not read contents of directory \"" + p.string() + "\" because the directory does not exist";
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


/** Format an error exception string for problems specifying the file/path name
 * @param[out] errorStr string to store the formatted error
*/
void RbFileManager::formatError(std::string& errorStr)
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


const std::string& RbFileManager::getFullFileName( void ) const
{
    return full_file_name;
}


/** Get absolute file path from file_path
 * @return absolute path
 */
std::string RbFileManager::getFullFilePath( void ) const
{
    path p = file_path;

    if (not p.is_absolute())
        p = fs::current_path() / p;

    p.make_preferred();

    return p.string();
}


/** Get the last path component of full_file_name
 * @note any trailing path separator is removed, so x/y/z/ will return z
 * @return last path component
 */
std::string RbFileManager::getLastPathComponent( void )
{
    return fs::path(full_file_name).parent_path().filename().string();
}


/** Get the last path component of a path
 * @note any trailing path separator is NOT removed, so x/y/z/ will return an empty string
 * @param s input path
 * @return last path component
 */
std::string getLastPathComponent(const std::string& s)
{
    auto ss = fs::path(s).filename().string();
    if (ss == ".")
        ss == "";
    return ss;
}


std::string getPathSeparator( void )
{
#   ifdef _WIN32
    return "\\";
#   else
    return "/";
#   endif
}


/** Removes the last path component from a path
 * @note any trailing path separator is NOT removed, so x/y/z/ will return x/y/z
 * @return string without the last path component
 */
std::string getStringByDeletingLastPathComponent(const std::string& s)
{
    return fs::path(s).parent_path().make_preferred().string();
}


/** Checks whether full_file_name is a path to an existing directory */
bool RbFileManager::isDirectory( void ) const
{
    return fs::is_directory(full_file_name);
}


/** Tests whether a directory is present (and is a directory)
 * @param mp path to check
 * @return result of the test
 */
bool isDirectoryPresent(const std::string &mp)
{
    return fs::is_directory(mp);
}


/** Checks whether the path given by file_path + file_name is a path to an existing file */
bool RbFileManager::isFile( void ) const
{
    auto f = fs::path(file_path) / fs::path(file_name);
    return fs::is_regular_file(f) and not fs::is_directory(f);
}


/** Checks whether the file name is non-empty */
bool RbFileManager::isFileNamePresent(void) const
{
    
    if ( file_name == "" )
    {
        return false;
    }
    
    return true;
}


/** Checks whether a file passed in as its path and file name components is present (and a file)
 * @param mp file path - if empty, set to the working directory
 * @param mf file name
 * @return whether the file exists
*/
bool isFilePresent(const std::string &mp, const std::string &mf)
{ 
    auto f = fs::path(mp) / fs::path(mf);
    
    return fs::is_regular_file(f) and not fs::is_directory(f);
}

/** Checks whether a file passed in as its full path is present (and is a file)
 * @param fn full file path
 * @return whether the file exists
*/
bool isFilePresent(const std::string &fn)
{
    fs::path f = fn;

    return fs::is_regular_file(f) and not fs::is_directory(f);
}

void RbFileManager::setFileName(std::string const &s)
{
    file_name = s;
}


void RbFileManager::setFilePath(std::string const &s)
{
    file_path = s;
#	ifdef _WIN32
    StringUtilities::replaceSubstring(file_path,"/","\\");
#   endif
    
}

/** Fills in a vector with the names of the files in the directory file_path
 *
 * @param[out] sv vector to store the names
 * @param[in] recursive whether to list recursively, default true
 *
 * @return true
*/
bool RbFileManager::setStringWithNamesOfFilesInDirectory(std::vector<std::string>& sv, bool recursive)
{
    std::vector<path> tmp;
    bool ok = RevBayesCore::setStringWithNamesOfFilesInDirectory(file_path, tmp, recursive);

    for(auto& f: tmp)
        sv.push_back( f.string() );
    
    return ok;
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
    auto dir = canonical(dirpath);

    for(auto& dir_entry: directory_iterator(dir))
    {
        auto entry_path = dir_entry.path();
            
        // Is this a symlink that points to something non-existant?
        if (not exists(entry_path)) continue;

        try
        {
            entry_path = canonical(entry_path);
        }
        catch(...)
        {
            continue;
        }

        auto entry_name = relative( entry_path, dir );

        if (entry_name.empty()) continue;

        // skip symlinks that point outside of the directory
        if (entry_name.begin()->filename_is_dot_dot()) continue;

        // can this happen?
        if (entry_name.begin()->filename_is_dot()) continue;

        if (is_directory(entry_path))
        {
            if (recursive)
                setStringWithNamesOfFilesInDirectory( entry_path, sv, recursive );
        }
        else
            sv.push_back( entry_path );
    }
    
    // make sure that the file names are sorted
    std::sort(sv.begin(), sv.end());
    
    return true;
}
}
