#include <cache/post_processor.templates.h>

namespace cache {

template class PostprocessorParameters<2>;
template class Postprocessor<2>;
template class PostprocessorParameters<3>;
template class Postprocessor<3>;

template class SuperCapacitorPostprocessorParameters<2>;
template class SuperCapacitorPostprocessor<2>;
template class SuperCapacitorPostprocessorParameters<3>;
template class SuperCapacitorPostprocessor<3>;

} // end namespace cache
