#include <string>
#include <regex>
#include <sys/stat.h>

std::string revCompDNA(std::string DNA);

bool fileExists(const std::string &filename);

void split(std::string str, std::string splitBy, std::vector,std::string>& tokens);
