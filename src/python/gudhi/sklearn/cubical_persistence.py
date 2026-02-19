# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

__license__ = "MIT"


import numpy as np
from numpy.typing import ArrayLike
from typing import Union, Literal, Optional
from sklearn.base import BaseEstimator, TransformerMixin
from joblib import Parallel, delayed
from array_api_compat import numpy as numpy_namespace
from ..array_api.cubical_layer import cubical_persistence

# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#   participant USER
#   participant CP as CubicalPersistence
#   USER->>CP: fit_transform(X)
#   CP->>thread1: _tranform(X[0])
#   CP->>thread2: _tranform(X[1])
#   Note right of CP: ...
#   thread1->>CP: [array( H0(X[0]) ), array( H1(X[0]) )]
#   thread2->>CP: [array( H0(X[1]) ), array( H1(X[1]) )]
#   Note right of CP: ...
#   CP->>USER: [[array( H0(X[0]) ), array( H1(X[0]) )],<br/> [array( H0(X[1]) ), array( H1(X[1]) )],<br/> ...]


class CubicalPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for computing the persistence diagrams from a cubical complex.
    """

    def __init__(
        self,
        homology_dimensions: Union[int, ArrayLike],
        input_type: Literal["top_dimensional_cells", "vertices"] = "top_dimensional_cells",
        homology_coeff_field: int = 11,
        min_persistence: float = 0.0,
        n_jobs: Optional[int] = None,
    ):
        """
        Constructor for the CubicalPersistence class.

        Parameters:
            homology_dimensions: The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
                dimension matters (in other words, when `homology_dimensions` is an int).
            input_type: 'top_dimensional_cells' if the filtration values passed to `transform()` are those of the
                top-dimensional cells, 'vertices' if they correspond to the vertices.
            homology_coeff_field: The homology coefficient field. Must be a prime number. Default value is 11.
            min_persistence: The minimum persistence value to take into account (strictly greater than
                `min_persistence`). Default value is `0.0`. Set `min_persistence` to `-1.0` to see all values.
            n_jobs: cf. https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html
        """
        self.homology_dimensions = homology_dimensions
        self.input_type = input_type
        self.homology_coeff_field = homology_coeff_field
        self.min_persistence = min_persistence
        self.n_jobs = n_jobs

        # Done twice (in __init__ and fit), but exception is better the sooner
        dim_list = np.asarray(self.homology_dimensions, dtype=int)
        if dim_list.ndim not in [0, 1]:
            raise ValueError(f"Invalid dimension. Got {self.homology_dimensions=}, expected type=int|ArrayLike[int].")

    def fit(self, X, Y=None):
        """
        Fit the `CubicalPersistence` class in function of `homology_dimensions` type.
        """
        # Must be in the `fit` part, as `transform` should be const and as `__init__` is not called on a parallel grid
        # search for instance
        self._dim_list = np.asarray(self.homology_dimensions, dtype=int)
        self._unwrap = False
        if self._dim_list.ndim == 0:
            self._unwrap = True
            self._dim_list = self._dim_list.reshape(1)
        return self

    def transform(self, X, Y=None):
        """Compute all the cubical complexes and their associated persistence diagrams.

        :param X: Filtration values of the top-dimensional cells or vertices for each complex.
        :type X: list of array-like

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]`
              - If `homology_dimensions` was set to `[i, j]`:
                `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of (,2) array_like or list of list of (,2) array_like
        """
        # threads is preferred as cubical construction and persistence computation releases the GIL
        res = Parallel(n_jobs=self.n_jobs, prefer="threads")(
            delayed(cubical_persistence)(
                np.asarray(cells),
                homology_dimensions=self._dim_list,
                input_type=self.input_type,
                min_persistence=self.min_persistence,
                homology_coeff_field=self.homology_coeff_field,
                array_api_namespace=numpy_namespace,
            )
            for cells in X
        )
        # cf. `fit`
        if self._unwrap:
            res = [d[0] for d in res]
        return res
