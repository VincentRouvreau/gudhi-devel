import numpy as np
from gudhi.array_api import cubical_persistence, CubicalLayer
import torch
from array_api_compat import get_namespace, is_torch_array, is_numpy_array

# Taken from (and to avoid to import sklearn and download dataset):
# from sklearn import datasets
# digit = datasets.load_digits().images[10]
digit = np.array([[ 0.,  0.,  1.,  9., 15., 11.,  0.,  0.],
       [ 0.,  0., 11., 16.,  8., 14.,  6.,  0.],
       [ 0.,  2., 16., 10.,  0.,  9.,  9.,  0.],
       [ 0.,  1., 16.,  4.,  0.,  8.,  8.,  0.],
       [ 0.,  4., 16.,  4.,  0.,  8.,  8.,  0.],
       [ 0.,  1., 16.,  5.,  1., 11.,  3.,  0.],
       [ 0.,  0., 12., 12., 10., 10.,  0.,  0.],
       [ 0.,  0.,  1., 10., 13.,  3.,  0.,  0.]])

def np_sort_dgm_by_lifetime(dgm):
    np_dgm = np.asarray(dgm)
    sorted_indices = np.argsort(np_dgm[:, 0] - np_dgm[:, 1])
    return np_dgm[sorted_indices]
    
def test_cubical_persistence_function_array_api():
    for input_type in ["top_dimensional_cells", "vertices"]:
        dimensions = [0, 1]
        np_dgms = cubical_persistence(digit, homology_dimensions=dimensions, input_type = input_type)
        assert len(np_dgms) == len(dimensions)
    
        X = torch.tensor(digit)
        torch_dgms = cubical_persistence(X, homology_dimensions=dimensions, input_type = input_type)
        assert len(torch_dgms) == len(dimensions)
    
        for idx in range(len(dimensions)):
            np_dgm = np_dgms[idx]
            assert is_numpy_array(np_dgm)
            sorted_np_dgm = np_sort_dgm_by_lifetime(np_dgm)
            torch_dgm = torch_dgms[idx]
            assert is_torch_array(torch_dgm)
            sorted_np_form_torch_dgm = np_sort_dgm_by_lifetime(torch_dgm)
            np.testing.assert_array_almost_equal(sorted_np_dgm, sorted_np_form_torch_dgm, decimal=6)


def test_cubical_persistence_function_array_api_with_gradients():
    X = torch.tensor([[0.,2.,2.],[2.,2.,2.],[2.,2.,1.]], requires_grad=True)
    dgm_0 = cubical_persistence(X, homology_dimensions=[0])[0]
    # Remove inf values
    mask = ~torch.isinf(dgm_0).any(dim=1)
    finite_dgm_0 = dgm_0[mask]
    loss = torch.sum(torch.square(0.5 * (finite_dgm_0[:, 1] - finite_dgm_0[:, 0])))
    grads = torch.autograd.grad(loss, X)
    np_grads = np.asarray(grads[0])
    np_expected_grads = np.asarray([[ 0.0000,  0.0000,  0.0000],
                                 [ 0.0000,  0.5000,  0.0000],
                                 [ 0.0000,  0.0000, -0.5000]])
    np.testing.assert_array_almost_equal(np_grads, np_expected_grads, decimal=6)

# TODO:
#  # JAX Array API
# import jax.numpy as jnp
# X = jnp.array(digits)
# cl(X)
# [(Array([[0., 8.],
#        [0., 9.]], dtype=float32), Array([[0.]], dtype=float32)), (Array([[10., 11.],
#        [11., 14.],
#        [11., 16.],
#        [10., 16.]], dtype=float32), Array([], shape=(0, 1), dtype=float32))]
