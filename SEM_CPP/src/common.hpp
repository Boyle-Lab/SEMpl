#ifndef COMMON_HPP
#define COMMON_HPP

#include "iterativeSEM.hpp"
#include <iostream>
#include <sstream>
#include <regex>
#include <dirent.h>
#include <sys/stat.h>
#include <stdbool.h>

std::string revCompDNA(std::string DNA);

bool fileExists(const std::string &filename);

void split_string(const std::string &str, const std::string &splitBy,
		   std::vector<std::string>& tokens);
void split_string_white(const std::string &str,
		   std::vector<std::string>& tokens);

std::string grab_string_at_index(const std::string &str, const size_t index,
								 const std::string &split);

std::string grab_string_at_index_white(const std::string &str, const size_t index);

std::string grab_string_last_index(const std::string &s);

void GetFilesInDirectory(std::vector<std::string> &out,
						 const std::string &directory);

int getLength(const Dataset &data);

void grab_string_at_index(const std::string s, std::string &out, const size_t idx);

void grab_string_3_index(std::string s, std::string &out);

void grab_string_4_index(std::string s, std::string &out);

std::stringstream exec(const char* cmd);

uint64_t encode2bit(const char *original);

#endif /* COMMON_HPP */
