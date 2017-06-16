#ifndef COMMON_HPP
#define COMMON_HPP

#include "iterativeSEM.hpp"
#include <iostream>
#include <regex>
#include <dirent.h>
#include <sys/stat.h>

std::string revCompDNA(std::string DNA);

bool fileExists(const std::string &filename);

void split_string(const std::string &str, const std::string &splitBy, 
		   std::vector<std::string>& tokens);

std::string grab_string_at_index(const std::string &str, const size_t index, 
								 const std::string &split);

void GetFilesInDirectory(std::vector<std::string> &out, 
						 const std::string &directory);

int getLength(const Dataset &data);

#endif /* COMMON_HPP */
