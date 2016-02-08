Appendix: Weak Formulation
==========================

Strong formulation
------------------

- In the collector:
  :math:`i_1 = -\sigma \nabla \Phi_1`
  :math:`\nabla \cdot i_1 = 0` 

- In the electrode:
  :math:`i_1 = -\sigma \nabla \Phi_1`                
  :math:`i_2 = -\kappa \nabla \Phi_2`                
  :math:`\nabla \cdot i_1 = \nabla \cdot i_2 = a i_n`


- In the seperator:
  :math:`i_2 = -\kappa \nabla \Phi_2`
  :math:`\nabla \cdot i_2 = 0`       

Weak formulation
----------------

- In the collector:

.. math::

  -\int_{\Omega_c}dr \phi_{1,i} \sigma \Delta \Phi_{1,j}\phi_{1,j} =
  -\int_{\partial {\Omega_c}}dr \phi_{1,i} \sigma \nabla \Phi_{1,j}\phi_{1,j}  +
  \int_{\Omega_c}dr \nabla \phi_{1,i} \sigma \nabla \Phi_{1,j}\phi_{1,j}

- In the electrode:

.. math::

  \begin{pmatrix}
  \int_{\Omega_e}dr \phi_{1,i} (-\sigma) \Delta \Phi_{1,j}\phi_{1,j}  + a C
  \frac{\partial \Phi_{1,j}\phi_{1,j}}{\partial t} & 
  - \int_{\Omega_e}dr \phi_{1,i} aC \frac{\partial \Phi_{2,j}\phi_{2,j}}{\partial
    t} \\
  -\int_{\Omega_e} dr \phi_{2,i} ac \frac{\partial \phi_{1,j}}{\partial t}  &
  \int_{\Omega_2} dr \phi_{2,i} (-\kappa) \Delta \Phi_{2,j} \phi_{2,j} -
  \phi_{2,i} aC \frac{\partial \Phi_{2,j}\phi_{2,j}}{\partial t} 
  \end{pmatrix} 
  =
  \begin{pmatrix}
  -\int_{\partial \Omega_e} dr \phi_{1,i} \ sigma \nabla \Phi_{1,j}\phi_{1,j} +
  \int_{\Omega_e}dr \nabla \phi_{1,i} \sigma) \nabla \Phi_{1,j}\phi_{1,j}  + a C
  \frac{\partial \Phi_{1,j}\phi_{1,j}}{\partial t} & 
  - \int_{\Omega_e}dr \phi_{1,i} aC \frac{\partial \Phi_{2,j}\phi_{2,j}}{\partial
    t} \\
  -\int_{\Omega_e} dr \phi_{2,i} aC \frac{\partial \phi_{1,j}}{\partial t}  &
  -\int_{\partial \omega_e} dr \phi_{2,i}\kappa \nabla \Phi_{2,j} \phi_{2,j} +
  \int_{\Omega_e} dr \nabla \phi_{2,i} \kappa \nabla \Phi_{2,j} \phi_{2,j} -
  \phi_{2,i} aC \frac{\partial \Phi_{2,j}\phi_{2,j}}{\partial t} 
  \end{pmatrix}

- In the seperator:

.. math::

  -\int_{\Omega_s}dr \phi_{2,i} \kappa \Delta \Phi_{2,j}\phi_{2,j} =
  -\int_{\partial {\Omega_s}}dr \phi_{2,i} \kappa \nabla \Phi_{2,j}\phi_{2,j} +
  \int_{\Omega_s}dr \nabla \phi_{2,i} \sigma \nabla \Phi_{2,j}\phi_{2,j}` 
