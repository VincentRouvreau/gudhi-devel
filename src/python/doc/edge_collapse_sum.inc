.. table::
   :widths: 30 40 30

   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | .. figure::                                                     | Edge collapse is able to reduce any flag filtration to a smaller     | :Author: Siddharth Pritam                                        |
   |      ../../doc/Collapse/dominated_edge.png                      | flag filtration with the same persistence, using only the            |                                                                  |
   |      :figclass: align-center                                    | 1-skeletons of a simplicial complex. The reduction is exact and the  | :Since: GUDHI 3.3.0                                              |
   |                                                                 | persistence homology of the reduced sequence is identical to the     |                                                                  |
   |                                                                 | persistence homology of the input sequence. The resulting method is  | :License: MIT                                                    |
   |                                                                 | simple and extremely efficient.                                      |                                                                  |
   |                                                                 |                                                                      |                                                                  |
   |                                                                 | Computation of edge collapse and persistent homology of a filtered   |                                                                  |
   |                                                                 | flag complex via edge collapse as described in                       |                                                                  |
   |                                                                 | :cite:`edgecollapsesocg2020`.                                        | :Requires: `Eigen <installation.html#eigen>`_ :math:`\geq` 3.1.0 |
   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | * :doc:`edge_collapse_user`                                     |                                                                                                                                         |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
