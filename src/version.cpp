#include "../include/version.hpp"

// Defined on the command line by the Makefile
#ifndef COSPT_GIT_SHA
#define COSPT_GIT_SHA "unknown"
#endif
#ifndef COSPT_BUILD_TIME
#define COSPT_BUILD_TIME "unknown"
#endif

std::string build_git_sha = COSPT_GIT_SHA;
std::string build_git_time = COSPT_BUILD_TIME;
