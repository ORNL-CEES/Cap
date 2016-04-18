# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from .PyCap import *
from .data_helpers import *
from .time_evolution import *
from .end_criterion import *
from .stage import *
from .charge_discharge import *
from .voltammetry import *
from .ragone_plot import *
from .impedance_spectroscopy import *
from .observer_pattern import *
from mpi4py import MPI

__all__ = ['PyCap', 'data_helpers', 'time_evolution', 'end_criterion',
           'stage', 'charge_discharge', 'voltammetry', 'ragone_plot',
           'impedance_spectroscopy', 'observer_pattern']

__version__ = PyCap.__version__
__git_branch__ = PyCap.__git_branch__
__git_commit_hash__ = PyCap.__git_commit_hash__
__git_remote_url__ = PyCap.__git_remote_url__
__doc__ = PyCap.__doc__


# Override EnergyStorageDevice.__init__(...) to add a default value to the
# ``comm`` parameter.
def build(self, ptree, comm=MPI.COMM_SELF):
    """
    Parameters
    ----------
    ptree : pycap.PropertyTree
        Property Tree
    comm : mpi4py.MPI.Comm
        Communicator (the default value is mpi4py.MPI.COMM_SELF).
    """
    return EnergyStorageDevice.build(self, ptree, comm)

EnergyStorageDevice.__init__ = build;
