import os
c.IPClusterEngines.engine_launcher_class = 'MPI'
if 'JPY_ENGINES' in os.environ:
    c.IPClusterEngines.n = int(os.environ['JPY_ENGINES'])
