import numpy as np
from sklearn import datasets
from gudhi.sklearn import CubicalPersistence
import gudhi

digits = datasets.load_digits().images[10]

  # Sklearn
cp = CubicalPersistence(homology_dimensions=[0, 1], n_jobs=-2)
cp.fit_transform([digits])
# [[array([[ 0.,  8.],
#        [ 0.,  9.],
#        [ 0., inf]]), array([[11., 16.],
#        [11., 14.],
#        [10., 16.],
#        [10., 11.]])]]

  # tensorflow
from gudhi.tensorflow import CubicalLayer
import tensorflow as tf
X = tf.Variable(digits)
cl = CubicalLayer(homology_dimensions=[0, 1])
cl.call(X)
# [(<tf.Tensor: shape=(2, 2), dtype=float64, numpy=
# array([[0., 8.],
#        [0., 9.]])>, <tf.Tensor: shape=(1, 1), dtype=float64, numpy=array([[0.]])>), (<tf.Tensor: shape=(4, 2), dtype=float64, numpy=
# array([[10., 11.],
#        [11., 14.],
#        [11., 16.],
#        [10., 16.]])>, <tf.Tensor: shape=(0, 1), dtype=float32, numpy=array([], shape=(0, 1), dtype=float32)>)]

  # PyTorch Array API
from gudhi.array_api import CubicalLayer
import torch
X = torch.tensor(digits)
cl = CubicalLayer(homology_dimensions=[0, 1])
cl(X)
# [(tensor([[0., 8.],
#         [0., 9.]], dtype=torch.float64), tensor([[0.]], dtype=torch.float64)), (tensor([[10., 11.],
#         [11., 14.],
#         [11., 16.],
#         [10., 16.]], dtype=torch.float64), tensor([], size=(0, 1), dtype=torch.float64))]

X = torch.tensor([[0.,2.,2.],[2.,2.,2.],[2.,2.,1.]], requires_grad=True)
cl = CubicalLayer(homology_dimensions=[0, 1])
cl(X)
dgm = cl(X)[0][0]
loss = torch.sum(torch.square(0.5 * (dgm[:, 1] - dgm[:, 0])))
grads = torch.autograd.grad(loss, X)
print(grads)
# (tensor([[ 0.0000,  0.0000,  0.0000],
#         [ 0.0000,  0.5000,  0.0000],
#         [ 0.0000,  0.0000, -0.5000]]),)

from gudhi.array_api import cubical_persistence
dgm = cubical_persistence(X, homology_dimensions=[0, 1])
loss = torch.sum(torch.square(0.5 * (dgm[0][0][:, 1] - dgm[0][0][:, 0])))
grads = torch.autograd.grad(loss, X)
print(grads)
# (tensor([[ 0.0000,  0.0000,  0.0000],
#         [ 0.0000,  0.5000,  0.0000],
#         [ 0.0000,  0.0000, -0.5000]]),)

  # NumPy Array API
cl(digits)
# [(array([[0., 8.],
#        [0., 9.]]), array([[0.]])), (array([[10., 11.],
#        [11., 14.],
#        [11., 16.],
#        [10., 16.]]), array([], shape=(0, 1), dtype=float64))]

  # JAX Array API
import jax.numpy as jnp
X = jnp.array(digits)
cl(X)
# [(Array([[0., 8.],
#        [0., 9.]], dtype=float32), Array([[0.]], dtype=float32)), (Array([[10., 11.],
#        [11., 14.],
#        [11., 16.],
#        [10., 16.]], dtype=float32), Array([], shape=(0, 1), dtype=float32))]


### input_type="vertices"

cp = CubicalPersistence(homology_dimensions=[0, 1], input_type="vertices", n_jobs=-2)
cp.fit_transform([digits])
# [[array([[ 0.,  8.],
#        [ 0., 12.],
#        [ 0., inf]]), array([[10., 11.],
#        [15., 16.]])]]

from gudhi.array_api import CubicalLayer
import torch
X = torch.tensor(digits)
cl = CubicalLayer(homology_dimensions=[0, 1], input_type="vertices")
cl(X)
# [(tensor([[ 0.,  8.],
#         [ 0., 12.]], dtype=torch.float64), tensor([[0.]], dtype=torch.float64)), (tensor([[10., 11.],
#         [15., 16.]], dtype=torch.float64), tensor([], size=(0, 1), dtype=torch.float64))]
