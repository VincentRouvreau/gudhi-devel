# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2026 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from array_api_compat import get_namespace, is_torch_namespace, is_jax_namespace
import numpy as np
import warnings
from typing import Literal, Iterable, Optional

from ..cubical_complex import CubicalComplex
from .._cubical_complex_ext import (
    _Bitmap_cubical_complex_interface,
    _Cubical_complex_persistence_interface,
    _Bitmap_cubical_complex_interface_float,
    _Cubical_complex_persistence_interface_float,
)

from .._pers_cub_low_dim_ext import (
    _persistence_on_a_line,
    _persistence_on_rectangle_from_top_cells,
)

# jax is an optional dependency (torch/numpy users shouldn't need it installed), so this
# import is guarded rather than unconditional at module scope.
try:
    import jax
    from functools import partial

    _HAS_JAX = True
except ImportError:
    _HAS_JAX = False


def _persistence_from_cells_without_autodiff(
    array_api_namespace,
    cells,
    homology_dimensions: Iterable[int],
    input_is_from_top_cells: bool,
    min_persistence: float = 0.0,
    homology_coeff_field: int = 11,
):
    xp = array_api_namespace
    cells = xp.asarray(cells)
    if len(cells.shape) == 1 and min_persistence >= 0.0:
        res = xp.asarray(_persistence_on_a_line(cells))
        if min_persistence > 0.0:
            # It would be more efficient inside _persistence_on_a_line, but not worth it?
            res = res[res[:, 1] - res[:, 0] > min_persistence]
        # Wasteful if dim_list_ does not contain 0, but that seems unlikely.
        return [res if i == 0 else xp.empty((0, 2)) for i in homology_dimensions]

    if len(cells.shape) == 2 and input_is_from_top_cells and min_persistence >= 0.0:
        if cells.size == 0:
            diags = [xp.empty((0, 2)), xp.empty((0, 2))]
        elif cells.shape[0] == 1 or cells.shape[1] == 1:
            diags = [xp.asarray(_persistence_on_a_line(cells.reshape(-1))), xp.empty((0, 2))]
        elif cells.shape[0] == 2:
            diags = [xp.asarray(_persistence_on_a_line(cells.min(0))), xp.empty((0, 2))]
        elif cells.shape[1] == 2:
            diags = [xp.asarray(_persistence_on_a_line(cells.min(1))), xp.empty((0, 2))]
        else:
            diags = _persistence_on_rectangle_from_top_cells(cells, min_persistence)
        return [xp.asarray(diags[i]) if i in (0, 1) else xp.empty((0, 2)) for i in homology_dimensions]

    xp = array_api_namespace
    # Compute pixels associated to positive and negative simplices
    # Don't compute gradient for this operation
    Xflat = xp.reshape(cells, [-1])
    # We reverse the dimensions because CubicalComplex uses Fortran ordering
    Xdim = xp.asarray(cells.shape[::-1])

    if Xflat.dtype == xp.float64:
      cc = _Bitmap_cubical_complex_interface(Xdim, Xflat, input_is_from_top_cells)
      pers = _Cubical_complex_persistence_interface(cc, True)
    elif Xflat.dtype == xp.float32:
        cc = _Bitmap_cubical_complex_interface_float(Xdim, Xflat, input_is_from_top_cells)
        pers = _Cubical_complex_persistence_interface_float(cc, True)
    else:
        raise TypeError(f"Unknown cells type {Xflat.dtype}")
    pers._compute_persistence(homology_coeff_field, 0.0)

    return [xp.asarray(pers._intervals_in_dimension(dim)) for dim in homology_dimensions]


def _compute_cubical_persistence_pairs(array_api_namespace, Xflat, Xdim, input_is_from_top_cells, homology_coeff_field):
    """
    Runs the C++ backend directly on Xflat/Xdim and returns the raw, variable-length array
    of persistence pairs, as produced by `_cofaces_of_cubical_persistence_pairs` /
    `_vertices_of_cubical_persistence_pairs`.

    Xflat/Xdim must be *concrete*, buffer-backed arrays -- a numpy array, or a torch tensor
    (CPU, contiguous) -- which nanobind can read with no copy. They must NOT be a
    torch/jax autodiff Tracer, since those have no buffer nanobind can extract.

    `array_api_namespace` is only used to compare `Xflat.dtype` against `xp.float32` /
    `xp.float64` (dtype objects differ between namespaces, e.g. `torch.float64` vs
    `np.float64`, so this can't be hardcoded to numpy).
    """
    xp = array_api_namespace
    if Xflat.dtype == xp.float64:
        cc = _Bitmap_cubical_complex_interface(Xdim, Xflat, input_is_from_top_cells)
        pers = _Cubical_complex_persistence_interface(cc, True)
    elif Xflat.dtype == xp.float32:
        cc = _Bitmap_cubical_complex_interface_float(Xdim, Xflat, input_is_from_top_cells)
        pers = _Cubical_complex_persistence_interface_float(cc, True)
    else:
        raise TypeError(f"Unknown cells type {Xflat.dtype}")

    pers._compute_persistence(homology_coeff_field, 0.0)

    # TODO: verify the return type of cofaces_of_cubical_persistence_pairs() by nanobind
    # a copy is perhaps avoidable?
    if input_is_from_top_cells:
        return np.array(pers._cofaces_of_cubical_persistence_pairs())
    return np.array(pers._vertices_of_cubical_persistence_pairs())


def _compute_padded_pers_indices_numpy(
    cells_np, homology_dimensions, input_is_from_top_cells, homology_coeff_field, max_pairs
):
    """
    Same underlying computation as `_compute_cubical_persistence_pairs`, but repackaged
    into a *fixed-shape* array so it can be returned from a `jax.pure_callback`, whose
    output shape/dtype must be declared ahead of time -- before the number of persistence
    pairs is actually known. Only used on the JAX path: `jax.pure_callback` always hands
    its Python callback concrete numpy arrays, regardless of the caller's namespace, so
    this can be numpy-specific.

    Returns:
        indices: int32 array of shape (len(homology_dimensions), max_pairs, 2).
            indices[i, j] = [birth_flat_index, death_flat_index] for the j-th finite pair
            in homology_dimensions[i], or [-1, -1] for unused padding slots. int32 rather
            than int64 because jax.pure_callback's declared result dtype must be 32-bit
            unless the caller has opted into `jax_enable_x64` -- fine here since these are
            just flat-array indices, well within 32-bit range for any realistic input.
        index_essential: flat index of the global minimum cell (birth of the essential
            homology-dimension-0 feature).
    """
    Xflat_np = cells_np.reshape(-1)
    Xdim_np = np.asarray(cells_np.shape[::-1], dtype=np.int32)

    pers_pairs = _compute_cubical_persistence_pairs(
        np, Xflat_np, Xdim_np, input_is_from_top_cells, homology_coeff_field
    )

    # Get only finite persistence pairs - when top-dimensional coface of negative simplex is not -1
    finite_pers_pairs_ind = np.argwhere(pers_pairs[:, 2] != -1)[:, 0]
    finite_pers_pairs = pers_pairs[finite_pers_pairs_ind]

    indices = np.full((len(homology_dimensions), max_pairs, 2), -1, dtype=np.int32)
    for idx_dim, dimension in enumerate(homology_dimensions):
        hidxs = np.argwhere(finite_pers_pairs[:, 0] == dimension)[:, 0]
        pairs = finite_pers_pairs[hidxs][:, 1:]
        indices[idx_dim, : len(pairs)] = pairs

    index_essential = np.argmin(Xflat_np)
    return indices, np.int32(index_essential)


if _HAS_JAX:

    @partial(jax.custom_jvp, nondiff_argnums=(1, 2, 3, 4))
    def _pers_indices_via_callback(cells, homology_dimensions, input_is_from_top_cells, homology_coeff_field, n_cells):
        """
        Runs the C++ backend on 'cells' via jax.pure_callback and returns the padded index
        arrays from `_compute_padded_pers_indices_numpy`.

        jax.pure_callback has NO default differentiation rule: JAX raises "Pure callbacks
        do not support JVP" the instant it's asked to differentiate through one, even
        trivially -- there's no implicit "treat it as a constant" fallback, unlike some
        other frameworks. It must be given an explicit rule via jax.custom_jvp (see
        `_pers_indices_via_callback_jvp` below).

        `homology_dimensions` is a nondiff arg here and so must be hashable -- callers pass
        a tuple, not a list.
        """
        result_shape_dtypes = (
            jax.ShapeDtypeStruct((len(homology_dimensions), n_cells, 2), np.int32),
            jax.ShapeDtypeStruct((), np.int32),
        )

        def _callback(cells_concrete):
            return _compute_padded_pers_indices_numpy(
                np.asarray(cells_concrete), homology_dimensions, input_is_from_top_cells, homology_coeff_field, n_cells
            )

        return jax.pure_callback(_callback, result_shape_dtypes, cells)

    @_pers_indices_via_callback.defjvp
    def _pers_indices_via_callback_jvp(
        homology_dimensions, input_is_from_top_cells, homology_coeff_field, n_cells, primals, tangents
    ):
        """
        The JVP rule required by jax.custom_jvp above. `primals` here are already plain,
        concrete values (not Tracers) -- that's how the custom_jvp/pure_callback contract
        works -- so calling `_pers_indices_via_callback` again to get the primal output is
        safe and does not recurse into tracing.

        The outputs of `_pers_indices_via_callback` are integer flat-array *indices*: which
        positions in `cells` are birth/death cells. They have no meaningful derivative of
        their own -- the actual differentiable computation is the xp.take(Xflat, indices)
        gather that happens afterwards, outside this function, on the real (traced) `cells`.

        Since the outputs are integer-dtyped, their tangents must use JAX's special
        zero-size `float0` dtype (JAX's marker for "this value has no meaningful tangent"),
        not `zeros_like` in the same integer dtype -- `jax.custom_jvp` raises a dtype
        mismatch error otherwise ("expecting tangent float0[...]").
        """
        (cells,) = primals
        primal_out = _pers_indices_via_callback(cells, homology_dimensions, input_is_from_top_cells, homology_coeff_field, n_cells)
        tangent_out = jax.tree_util.tree_map(
            lambda x: jax.numpy.zeros(x.shape, dtype=jax.dtypes.float0), primal_out
        )
        return primal_out, tangent_out


def _dgms_from_flat_and_indices(xp, Xflat, homology_dimensions, indices_by_dim, index_essential, min_persistence):
    """
    Shared final step for both the torch/numpy and jax autodiff paths: turns per-dimension
    flat index arrays into persistence diagrams by gathering from `Xflat`. This is the part
    that must stay a *traced, differentiable* array-API operation, so that autodiff
    (torch.autograd / jax.grad) can backprop through the `xp.take` gathers.
    """
    dgms = []
    for idx_dim, dimension in enumerate(homology_dimensions):
        indices_flat = indices_by_dim[idx_dim]
        finite_dgm = xp.reshape(xp.take(Xflat, indices_flat), [-1, 2])

        if min_persistence >= 0:
            pers = xp.abs(finite_dgm[:, 1] - finite_dgm[:, 0])
            mask = pers > min_persistence
            sel = xp.nonzero(mask)[0]
            finite_dgm = xp.take(finite_dgm, sel, axis=0)

        if dimension == 0:
            essential_dgm = xp.reshape(xp.take(Xflat, index_essential), [-1, 1])
            # Extend with a +inf value at the end for essential diagram to return a bar code [birth, +inf]
            essential_dgm = xp.concat((essential_dgm, xp.asarray([[xp.inf]])), axis=1)
            dgm = xp.concat((finite_dgm, essential_dgm))
        else:
            dgm = finite_dgm
        dgms.append(dgm)
    return dgms


def _persistence_from_cells_with_autodiff(
    array_api_namespace,
    cells,
    homology_dimensions: Iterable[int],
    input_is_from_top_cells: bool,
    min_persistence: float = 0.0,
    homology_coeff_field: int = 11,
):
    if (is_torch_namespace(array_api_namespace) or is_jax_namespace(array_api_namespace)) == False:
        warnings.warn(
            """
            The input 'cells' is not a torch tensor, nor a jax numpy array and it is asked to preserve gradient.
            'preserve_gradient' is slower and not available for other namespaces than torch and jax
            """,
            UserWarning,
        )

    xp = array_api_namespace
    homology_dimensions = list(homology_dimensions)

    # Compute pixels associated to positive and negative simplices
    # Don't compute gradient for this operation
    Xflat = xp.reshape(cells, [-1])
    # We reverse the dimensions because CubicalComplex uses Fortran ordering
    Xdim = xp.asarray(cells.shape[::-1])

    if is_jax_namespace(xp):
        # Under jax.grad (or jax.jit), 'cells'/'Xflat' are abstract Tracers while this
        # function is traced -- a Tracer has no buffer nanobind can read, so the C++
        # backend can't be called on it directly (unlike torch below).
        # _pers_indices_via_callback forces 'cells' to materialize to real host memory
        # (even inside grad/jit) via jax.pure_callback, runs the C++ computation on that
        # concrete numpy array, and hands the (fixed-shape, padded) *indices* back into the
        # trace -- with an explicit jax.custom_jvp rule attached, since pure_callback has no
        # differentiation rule of its own (JAX errors outright otherwise, even for a
        # trivially-zero gradient). Gradients flow through the xp.take gathers in
        # _dgms_from_flat_and_indices instead, which are ordinary, differentiable JAX ops on
        # the still-traced Xflat. The callback necessarily makes a host-memory copy of
        # 'cells'; there's no way around that for a traced JAX array.
        if not _HAS_JAX:
            raise ImportError("preserve_gradient=True with a JAX array requires jax to be installed.")

        n_cells = int(np.prod(cells.shape))
        padded_indices, index_essential = _pers_indices_via_callback(
            cells, tuple(homology_dimensions), input_is_from_top_cells, homology_coeff_field, n_cells
        )

        indices_by_dim = []
        for idx_dim in range(len(homology_dimensions)):
            pair_indices = padded_indices[idx_dim]
            valid = pair_indices[:, 0] != -1
            indices_by_dim.append(xp.reshape(pair_indices[valid], [-1]))

        return _dgms_from_flat_and_indices(xp, Xflat, homology_dimensions, indices_by_dim, index_essential, min_persistence)

    # Torch (and, with the warning above, plain numpy) path: 'cells' is already a concrete,
    # buffer-backed tensor/array -- unlike the jax Tracer case, nanobind can read Xflat/Xdim
    # directly with no copy (for a contiguous CPU tensor), exactly as in the original code.
    # Gradients flow through the xp.take gathers in _dgms_from_flat_and_indices, which
    # torch.autograd differentiates natively; the C++ call below is simply not part of the
    # recorded autograd graph, which is fine since it's only used to determine *which*
    # indices are birth/death cells, not their values.
    index_essential = xp.argmin(Xflat)

    pers_pairs = _compute_cubical_persistence_pairs(xp, Xflat, Xdim, input_is_from_top_cells, homology_coeff_field)

    # Get only finite persistence pairs - when top-dimensional coface of negative simplex is not -1
    finite_pers_pairs_ind = np.argwhere(pers_pairs[:, 2] != -1)[:, 0]
    finite_pers_pairs = pers_pairs[finite_pers_pairs_ind]

    indices_by_dim = []
    for dimension in homology_dimensions:
        # xp.take requires a vector of indices with pytorch (but can be an array for numpy)
        hidxs = np.argwhere(finite_pers_pairs[:, 0] == dimension)[:, 0]
        indices_flat = np.reshape(finite_pers_pairs[hidxs][:, 1:], [-1])
        # Force dtype - maybe better somewhere else - maybe comes when empty
        indices_flat = xp.asarray(indices_flat, dtype=index_essential.dtype)
        indices_by_dim.append(indices_flat)

    return _dgms_from_flat_and_indices(xp, Xflat, homology_dimensions, indices_by_dim, index_essential, min_persistence)


def cubical_persistence(
    cells,
    homology_dimensions: Iterable[int],
    input_type: Literal["top_dimensional_cells", "vertices"] = "top_dimensional_cells",
    min_persistence: float = 0.0,
    homology_coeff_field: int = 11,
    preserve_gradient: bool = False,
    array_api_namespace=None,
):
    """
    Returns the persistent homology bar code from the cubical complex.
    Parameters:
        homology_dimensions: The returned persistence diagrams dimension(s).
        input_type: 'top_dimensional_cells' if the filtration values passed are those of the top-dimensional cells,
            'vertices' if they correspond to the vertices. Default is 'top_dimensional_cells'.
        homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11.
        min_persistence: The minimum persistence value to take into account (strictly greater than
            `min_persistence`). Default value is `0.0`. Set `min_persistence` to `-1.0` to see all values.
        preserve_gradient: Shall the function preserve the input gradient or not. Default value is `False` (faster).
            If the input is a PyTorch or JAX array with gradients, `preserve_gradient` can be set to `True` to allow
            persistence optimization.
    """
    if input_type not in ["top_dimensional_cells", "vertices"]:
        raise ValueError("input_type can only be 'top_dimensional_cells' or 'vertices'")

    if array_api_namespace is None:
        xp = get_namespace(cells)
    else:
        xp = array_api_namespace

    input_is_from_top_cells = input_type == "top_dimensional_cells"
    if preserve_gradient == True:
        return _persistence_from_cells_with_autodiff(
            xp, cells, homology_dimensions, input_is_from_top_cells, min_persistence, homology_coeff_field
        )
    else:
        return _persistence_from_cells_without_autodiff(
            xp, cells, homology_dimensions, input_is_from_top_cells, min_persistence, homology_coeff_field
        )


class CubicalLayer:
    """
    Layer for computing the persistent homology of a cubical complex
    """

    def __init__(
        self,
        homology_dimensions: Iterable[int],
        input_type: Literal["top_dimensional_cells", "vertices"] = "top_dimensional_cells",
        min_persistence: float = 0.0,
        homology_coeff_field: int = 11,
    ):
        """Constructor for the CubicalLayer class

        Parameters:
            homology_dimensions: list of homology dimensions.
            input_type: 'top_dimensional_cells' if the filtration values passed to `__call__()` are those of the
                top-dimensional cells, 'vertices' if they correspond to the vertices.
            min_persistence: minimum distance-to-diagonal of the points in the output persistence diagrams
                (default None, in which case `0.` is used for all dimensions)
            homology_coeff_field: homology field coefficient. Must be a prime number. Default value is 11.
        """
        self.homology_dimensions = homology_dimensions
        if input_type not in ["top_dimensional_cells", "vertices"]:
            raise ValueError("input_type can only be 'top_dimensional_cells' or 'vertices'")
        self.input_type = input_type
        self.input_top_cells = input_type == "top_dimensional_cells"
        self.min_persistence = min_persistence
        self.homology_coeff_field = homology_coeff_field

    def __call__(self, X):
        """Compute persistence diagram associated to a cubical complex filtered by some pixel values.

        Parameters:
            X (anything compatible to an Array API standard): pixel values of the cubical complex

        Returns:
            List[Tuple[tf.Tensor,tf.Tensor]]: List of cubical persistence diagrams. The length of this list is the same
                than that of dimensions, i.e., there is one persistence diagram per homology dimension provided in the
                input list dimensions. Moreover, the finite and essential parts of the persistence diagrams are
                provided separately: each element of this list is a tuple of size two that contains the finite and
                essential parts of the corresponding persistence diagram, of shapes [num_finite_points, 2] and
                [num_essential_points, 1] respectively. Note that the essential part is always empty in cubical
                persistence diagrams, except in homology dimension zero, where the essential part always contains a
                single point, with abscissa equal to the smallest value in the complex, and infinite ordinate.
        """
        self.dgms = cubical_persistence(
            X,
            homology_dimensions=self.homology_dimensions,
            input_type=self.input_type,
            min_persistence=self.min_persistence,
            homology_coeff_field=self.homology_coeff_field,
        )
        return self.dgms
