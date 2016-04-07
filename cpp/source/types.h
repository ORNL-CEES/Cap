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

#endif
