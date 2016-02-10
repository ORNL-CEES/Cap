/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#ifndef CAP_VERSION_H
#define CAP_VERSION_H

#include <string>

namespace cap
{

std::string version();

std::string git_branch();

std::string git_commit_hash();

std::string git_remote_url();

}

#endif
