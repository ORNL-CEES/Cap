#ifndef CAP_TYPES.h
#define CAP_TYPES.h

/**
 * Shorten dealii::TrilinosWrappers::MPI to dealii::Trilinos::MPI and
 * dealii::TrilinosWrappers to dealii::Trilinos.
 */
namespace dealii
{
 namespace Trilinos = TrilinosWrappers;
}

#endif
