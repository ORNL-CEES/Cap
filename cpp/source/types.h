#ifndef CAP_TYPES_H
#define CAP_TYPES_H

namespace dealii
{
namespace TrilinosWrappers
{
}

namespace parallel
{
namespace distributed
{
}
}

/**
 * Shorten dealii::TrilinosWrappers::MPI to dealii::Trilinos::MPI and
 * dealii::TrilinosWrappers to dealii::Trilinos.
 */
namespace Trilinos = TrilinosWrappers;

/**
 * Shorten dealii::parallel::distributed::Triangulation to
 * dealii::parallel::Triangulation
 */
namespace distributed = parallel::distributed;
}

#ifdef WITH_DEAL_II

#include <deal.II/base/types.h>
#include <limits>

namespace type
{

/**
 * Sometimes dealii::numbers::invalid_boundary_id triggers an overflow warning.
 * Use a typedef so we don't trigger the warning.
 */
static dealii::types::boundary_id const invalid_boundary_id =
    std::numeric_limits<::dealii::types::boundary_id>::max();
}

#endif

#endif
