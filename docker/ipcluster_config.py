import os
c.IPClusterEngines.engine_launcher_class = 'MPI'
if 'NPROC' in os.environ:
    c.IPClusterEngines.n = int(os.environ['NPROC'])
