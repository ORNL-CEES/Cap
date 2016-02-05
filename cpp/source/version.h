#ifndef CAP_VERSION_H
#define CAP_VERSION_H

#include <string>

namespace cap
{

std::string version();

std::string git_branch();

std::string git_commit_hash();
}

#endif
