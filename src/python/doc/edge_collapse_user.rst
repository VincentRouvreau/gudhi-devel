:orphan:

.. To get rid of WARNING: document isn't included in any toctree

.. include:: edge_collapse_sum.inc

This module implements edge collapse of a filtered flag complex, in particular it reduces a filtration of Vietoris-Rips
complex from its graph to another smaller flag filtration with same persistence. Where a filtration is a sequence of
simplicial (here Rips) complexes connected with inclusions.

Edge collapse definition
========================

An edge :math:`e` in a simplicial complex :math:`K` is called a **dominated edge** if the link of :math:`e` in
:math:`K`, :math:`lk_K(e)` is a simplicial cone, that is, there exists a vertex :math:`v^{\prime} \notin e` and a
subcomplex :math:`L` in :math:`K`, such that :math:`lk_K(e) = v^{\prime}L`. We say that the vertex :math:`v^{\prime}`
is {dominating} :math:`e` and :math:`e` is {dominated} by :math:`v^{\prime}`. 
An **elementary egde collapse** is the removal of a dominated edge :math:`e` from :math:`K`, 
which we denote with :math:`K` :math:`{\searrow\searrow}^1` :math:`K\setminus e`. 
The symbol :math:`\mathbf{K\setminus e}` (deletion of :math:`e` from :math:`K`) refers to the subcomplex of :math:`K`
which has all simplices of :math:`K` except :math:`e` and the ones containing :math:`e`.
There is an **edge collapse** from a simplicial complex :math:`K` to its subcomplex :math:`L`, 
if there exists a series of elementary edge collapses from :math:`K` to :math:`L`, denoted as :math:`K`
:math:`{\searrow\searrow}` :math:`L`.

An edge collapse is a homotopy preserving operation, and it can be further expressed as sequence of the classical
elementary simple collapse. 
A complex without any dominated edge is called a :math:`1`- minimal complex and the core :math:`K^1` of simplicial
complex is a minimal complex such that :math:`K` :math:`{\searrow\searrow}` :math:`K^1`.
Computation of a core (not unique) involves computation of dominated edges and the dominated edges can be easily
characterized as follows:

-- For general simplicial complex: An edge :math:`e \in K` is dominated by another vertex :math:`v^{\prime} \in K`,
<i>if and only if</i> all the maximal simplices of :math:`K` that contain :math:`e` also contain :math:`v^{\prime}`

-- For a flag complex: An edge :math:`e \in K` is dominated by another vertex :math:`v^{\prime} \in K`, <i>if and only
if</i> all the vertices in :math:`K` that has an edge with both vertices of :math:`e`  also has an edge with
:math:`v^{\prime}`.

The algorithm to compute the smaller induced filtration is described in Section 5 :cite:`edgecollapsesocg2020`.
Edge collapse can be successfully employed to reduce any given filtration of flag complexes to a smaller induced
filtration which preserves the persistent homology of the original filtration and is a flag complex as well.

The general idea is that we consider edges in the filtered graph and sort them according to their filtration value
giving them a total order.
Each edge gets a unique index denoted as :math:`i` in this order.  To reduce the filtration, we move forward with
increasing filtration value 
in the graph and check if the current edge :math:`e_i` is dominated in the current graph :math:`G_i := \{e_1, .. e_i\}`
or not. 
If the edge :math:`e_i` is dominated we remove it from the filtration and move forward to the next edge
:math:`e_{i+1}`.
If :math:`e_i` is non-dominated then we keep it in the reduced filtration and then go backward in the current graph
:math:`G_i` to look for new non-dominated edges that was dominated before but might become non-dominated at this
point. 
If an edge :math:`e_j, j < i` during the backward search is found to be non-dominated, we include :math:`e_j` in to the
reduced filtration and we set its new filtration value to be :math:`i` that is the index of :math:`e_i`.
The precise mechanism for this reduction has been described in Section 5 :cite:`edgecollapsesocg2020`. 
Here we implement this mechanism for a filtration of Rips complex.
After perfoming the reduction the filtration reduces to a flag-filtration with the same persistence as the original
filtration. 


.. autofunction:: gudhi.flag_complex_collapse_edges


Basic edge collapse
-------------------

This example calls :func:`~gudhi.flag_complex_collapse_edges` from a proximity graph represented as a list of
`tuples(int, int, double)`. Then it collapses edges and displays a new list of `tuples(int, int, double)` (with less
edges) that will preserve the persistence homology computation.

.. testcode::

    from gudhi import flag_complex_collapse_edges as collapse_edges

    # 1   2
    # o---o
    # |\ /|
    # | x |
    # |/ \|
    # o---o
    # 0   3
    graph = [(0, 1, 1.),
             (1, 2, 1.),
             (2, 3, 1.),
             (3, 0, 1.),
             (0, 2, 2.),
             (1, 3, 2.) ]
    print(collapse_edges(graph))

.. testoutput::

    [(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0), (0, 2, 2.0)]
