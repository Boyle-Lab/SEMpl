#include "common.h"
using namespace std;

string revCompDNA(string dna){
  string rev = "";
  for(int i = static_cast<int>(dna.size()) - 1; i >= 0; i++){
    rev += dna[i];
  }
  rev = regex_replace(dna, regex("ACGTacgt"), "TGCAtgca");
  return rev;
}

bool fileExists(const string &filename){
  struct stat buffer;
  return (stat (filename.c_str(), &buffer) == 0);
}

// int parse_wc{
//   string s;
//   while(getline(s, ))
// }

// grabbed from internet
void GetFilesInDirectory(std::vector<string> &out, const string &directory)
{
#ifdef WINDOWS
    HANDLE dir;
    WIN32_FIND_DATA file_data;

    if ((dir = FindFirstFile((directory + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
        return; /* No files found */

    do {
        const string file_name = file_data.cFileName;
        const string full_file_name = directory + "/" + file_name;
        const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

        if (file_name[0] == '.')
            continue;

        if (is_directory)
            continue;

        out.push_back(full_file_name);
    } while (FindNextFile(dir, &file_data));

    FindClose(dir);
#else
    DIR *dir;
    class dirent *ent;
    class stat st;

    dir = opendir(directory);
    while ((ent = readdir(dir)) != NULL) {
        const string file_name = ent->d_name;
        const string full_file_name = directory + "/" + file_name;

        if (file_name[0] == '.')
            continue;

        if (stat(full_file_name.c_str(), &st) == -1)
            continue;

        const bool is_directory = (st.st_mode & S_IFDIR) != 0;

        if (is_directory)
            continue;

        out.push_back(full_file_name);
    }
    closedir(dir);
#endif
} // GetFilesInDirectory


// function taken from internet that splits  a string by a character
void split(std::string str, std::string splitBy, std::vector<std::string>& tokens)
{
    /* Store the original string in the array, so we can loop the rest
     * of the algorithm. */
    tokens.push_back(str);

    // Store the split index in a 'size_t' (unsigned integer) type.
    size_t splitAt;
    // Store the size of what we're splicing out.
    size_t splitLen = splitBy.size();
    // Create a string for temporarily storing the fragment we're processing.
    std::string frag;
    // Loop infinitely - break is internal.
    while(true)
    {
        /* Store the last string in the vector, which is the only logical
         * candidate for processing. */
        frag = tokens.back();
        /* The index where the split is. */
        splitAt = frag.find(splitBy);
        // If we didn't find a new split point...
        if(splitAt == string::npos)
        {
            // Break the loop and (implicitly) return.
            break;
        }
        /* Put everything from the left side of the split where the string
         * being processed used to be. */
        tokens.back() = frag.substr(0, splitAt);
        /* Push everything from the right side of the split to the next empty
         * index in the vector. */
        tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
    }
}

// REQUIRES: index is within str
string grab_string_at_index(const string &str, size_t index){
    auto ptr = str.c_str();
    size_t val = 0;
    while(val < index){
            while(*ptr != ' ' && *ptr){
                ++ptr;
            }
            ++ptr;
            ++val;
    }

    auto index_ptr = ptr;
    while(*index_ptr != ' ' && *index_ptr){
        ++index_ptr;
    }
    return string(ptr, static_cast<size_t>(index_ptr - ptr));

}
