:orphan:

.. To get rid of WARNING: document isn't included in any toctree

====================================
TensorFlow layer for representations
====================================

.. include:: representations_sum.inc

`PersLay <http://proceedings.mlr.press/v108/carriere20a.html>`_ is a layer for neural network architectures that allows
to automatically learn the best representation to use for persistence diagrams in supervised machine learning during
training time. Its parameters allow to reproduce most of the known finite-dimensional representations (such as, e.g.,
landscapes and images), and can be combined to create even new ones, that are best suited for a given supervised
machine learning task with persistence diagrams. PersLay is implemented in TensorFlow.

This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-perslay-visu.ipynb>`__ explains how to use
PersLay.

Example
-------

PersLay
^^^^^^^


.. testcode:: perslay

    import numpy             as np
    import torch
    from sklearn.preprocessing import MinMaxScaler
    import gudhi.representations as gdr
    import gudhi.tensorflow.perslay as prsl

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = torch.tensor(diagrams, dtype=torch.float32)

    rho = torch.nn.Identity()
    phi = prsl.GaussianPerslayPhi((5, 5), ((-.5, 1.5), (-.5, 1.5)), .1)
    weight = prsl.PowerPerslayWeight(1.,0.)
    perm_op = torch.sum

    perslay = prsl.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)
    print(vectors)

.. testoutput:: perslay

    tensor([[[[1.7266e-16],
              [4.1715e-09],
              [8.0383e-06],
              [8.0269e-06],
              [9.0331e-13]],

             [[4.1706e-09],
              [1.0074e-01],
              [1.5803e+00],
              [1.3066e+00],
              [1.4955e-07]],

             [[1.1337e-08],
              [2.7384e-01],
              [8.2997e-01],
              [9.0923e+00],
              [1.5146e-04]],

             [[8.5739e-12],
              [3.0724e-02],
              [1.2395e+01],
              [6.1665e-02],
              [1.0205e-06]],

             [[2.1244e-14],
              [7.6157e-05],
              [3.0724e-02],
              [1.3949e-06],
              [7.8094e-16]]]], grad_fn=<StackBackward0>)

Perslay reference manual
------------------------

.. autoclass:: gudhi.tensorflow.perslay.Perslay
   :members:
   :special-members:
   :show-inheritance:

Weight functions
^^^^^^^^^^^^^^^^
.. autoclass:: gudhi.tensorflow.perslay.GaussianMixturePerslayWeight
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.GridPerslayWeight
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.PowerPerslayWeight
   :members:
   :special-members:
   :show-inheritance:

Phi functions
^^^^^^^^^^^^^
.. autoclass:: gudhi.tensorflow.perslay.FlatPerslayPhi
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.GaussianPerslayPhi
   :members:
   :special-members:
   :show-inheritance:

.. autoclass:: gudhi.tensorflow.perslay.TentPerslayPhi
   :members:
   :special-members:
   :show-inheritance:
