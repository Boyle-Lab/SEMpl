#include "iterativeSEM.hpp"
#include <iostream>
#include <regex>
#include <sys/stat.h>

std::string revCompDNA(std::string DNA);

bool fileExists(const std::string &filename);

void split(std::string str, std::string splitBy, std::vector<std::string>& tokens);

std::string grab_string_at_index(const std::string &str, size_t index);
