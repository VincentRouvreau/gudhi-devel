# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2026 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


from array_api_compat import get_namespace
import numpy as np
from ..cubical_complex import CubicalComplex
from typing import Iterable, Optional
from .._cubical_complex_ext import _Bitmap_cubical_complex_interface, _Cubical_complex_persistence_interface


class CubicalLayer():
    """
    Layer for computing the persistent homology of a cubical complex
    """
    def __init__(self, homology_dimensions: Iterable[int], input_type: str = "top_dimensional_cells",
                 min_persistence: Optional[Iterable[float]] = None, homology_coeff_field: int = 11):
        """Constructor for the CubicalLayer class

        Parameters:
            homology_dimensions: list of homology dimensions.
            input_type: 'top_dimensional_cells' if the filtration values passed to `__call__()` are those of the
                top-dimensional cells, 'vertices' if they correspond to the vertices.
            min_persistence: minimum distance-to-diagonal of the points in the output persistence diagrams
                (default None, in which case `0.` is used for all dimensions)
            homology_coeff_field: homology field coefficient. Must be a prime number. Default value is 11.
        """
        self.dimensions = homology_dimensions
        if input_type not in ["top_dimensional_cells", "vertices"]:
            raise ValueError("input_type can only be 'top_dimensional_cells' or 'vertices'")
        self.input_top_cells = (input_type == "top_dimensional_cells")
        self.min_persistence = min_persistence
        if min_persistence is None:
            self.min_persistence = np.zeros(len(self.dimensions), dtype=np.float64)
        self.hcf = homology_coeff_field
        if len(self.min_persistence) != len(self.dimensions):
            raise ValueError("'homology_dimensions' and 'min_persistence' do not have coherent sizes")

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
        xp = get_namespace(X)
        # Compute pixels associated to positive and negative simplices
        # Don't compute gradient for this operation
        Xflat = xp.reshape(X, [-1])
        # We reverse the dimensions because CubicalComplex uses Fortran ordering
        Xdim = xp.asarray(X.shape[::-1])
        
        # index of minimum pixel value for essential persistence diagram
        index_essential = xp.argmin(Xflat)

        cc = _Bitmap_cubical_complex_interface(Xdim, Xflat, self.input_top_cells)
        pers = _Cubical_complex_persistence_interface(cc, True)
        pers._compute_persistence(self.hcf, 0.)
        
        # TODO: verify the return type of cofaces_of_cubical_persistence_pairs() by nanobind
        # a copy is perhaps avoidable?
        if self.input_top_cells:
            pers_pairs = np.array(pers._cofaces_of_cubical_persistence_pairs())
        else:
            pers_pairs = np.array(pers._vertices_of_cubical_persistence_pairs())
        
        # Get only finite persistence pairs - when top-dimensional coface of negative simplex is not -1
        finite_pers_pairs_ind = np.argwhere(pers_pairs[:, 2] != -1)[:, 0]
        finite_pers_pairs = pers_pairs[finite_pers_pairs_ind]
        
        # Get persistence diagram by simply picking the corresponding entries in the image
        self.dgms = []
        for idx_dim, dimension in enumerate(self.dimensions):
            # xp.take requires a vector of indices with pytorch (but can be an array for numpy)
            hidxs = np.argwhere(finite_pers_pairs[:, 0] == dimension)[:, 0]
            indices_flat = np.reshape(finite_pers_pairs[hidxs][:, 1:], [-1])
            # Force dtype - maybe better somewhere else - maybe comes when empty
            indices_flat = xp.asarray(indices_flat, dtype=index_essential.dtype)
            finite_dgm = xp.reshape(xp.take(Xflat, indices_flat), [-1, 2])
            if dimension == 0:
                essential_dgm = xp.reshape(xp.take(Xflat, index_essential), [-1, 1])
            else:
                essential_dgm = xp.zeros([0, 1], dtype=X.dtype)

            min_pers = self.min_persistence[idx_dim]
            if min_pers >= 0:
                pers = xp.abs(finite_dgm[:, 1] - finite_dgm[:, 0])
                mask = pers > min_pers
                idx = xp.nonzero(mask)[0]
                finite_dgm = xp.take(finite_dgm, idx, axis=0)
                self.dgms.append((finite_dgm, essential_dgm))
            else:
                self.dgms.append((finite_dgm, essential_dgm))
        return self.dgms
