#include <cap/version.h>
#include <string>

namespace cap {

std::string version() { return CAP_VERSION; }

std::string git_commit_hash() { return GIT_COMMIT_HASH; }

}
