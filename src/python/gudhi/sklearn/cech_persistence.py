# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2024 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from .. import DelaunayComplex
from sklearn.base import BaseEstimator, TransformerMixin

# joblib is required by scikit-learn
from joblib import Parallel, delayed

# Mermaid sequence diagram - https://mermaid-js.github.io/mermaid-live-editor/
# sequenceDiagram
#   participant USER
#   participant R as CechPersistence
#   USER->>R: fit_transform(X)
#   Note right of R: homology_dimensions=[i,j]
#   R->>thread1: _tranform(X[0])
#   R->>thread2: _tranform(X[1])
#   Note right of R: ...
#   thread1->>R: [array( Hi(X[0]) ), array( Hj(X[0]) )]
#   thread2->>R: [array( Hi(X[1]) ), array( Hj(X[1]) )]
#   Note right of R: ...
#   R->>USER: [[array( Hi(X[0]) ), array( Hj(X[0]) )],<br/> [array( Hi(X[1]) ), array( Hj(X[1]) )],<br/>...]


class CechPersistence(BaseEstimator, TransformerMixin):
    """
    This is a class for constructing Čech complexes and computing the persistence diagrams from them.
    """

    def __init__(
        self,
        homology_dimensions,
        input_type="point cloud",
        filtration=True,
        precision="safe",
        max_alpha_square=float("inf"),
        homology_coeff_field=11,
        n_jobs=None,
    ):
        """
        Constructor for the CechPersistence class.

        Parameters:
            homology_dimensions (int or list of int): The returned persistence diagrams dimension(s).
                Short circuit the use of :class:`~gudhi.representations.preprocessing.DimensionSelector` when only one
                dimension matters (in other words, when `homology_dimensions` is an int).
            input_type (str): Can be 'point cloud' when inputs are point clouds, or 'weighted point cloud', when
                inputs are point clouds plus the weight (as the last column value). Default is 'point cloud'.
            filtration (bool): Can be True, the filtration value of each simplex is computed, or False, the filtration
                value of each simplex is not computed (set to NaN). Default is True.
            precision (str): Delaunay complex precision can be 'fast', 'safe' or 'exact'. Default is 'safe'.
            max_alpha_square (float): The maximum alpha square threshold the simplices shall not exceed. Default is set
                to infinity, and there is very little point using anything else since it does not save time.
            homology_coeff_field (int): The homology coefficient field. Must be a prime number. Default value is 11.
            n_jobs (int): Number of jobs to run in parallel. `None` (default value) means `n_jobs = 1` unless in a
                joblib.parallel_backend context. `-1` means using all processors. cf.
                https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for more details.
        """
        self.homology_dimensions = homology_dimensions
        self.input_type = input_type
        self.filtration = filtration
        self.precision = precision
        self.max_alpha_square = max_alpha_square
        self.homology_coeff_field = homology_coeff_field
        self.n_jobs = n_jobs
        if self.input_type not in ["point cloud", "weighted point cloud"]:
            raise ValueError("Unknown input type")
        if self.precision not in ["safe", "fast", "exact"]:
            raise ValueError("Unknown precision")

    def fit(self, X, Y=None):
        """
        Nothing to be done, but useful when included in a scikit-learn Pipeline.
        """
        return self

    def __transform(self, inputs):
        max_dimension = max(self.dim_list_) + 1

        # Default filtration value
        fltr = None
        if self.input_type == "point cloud":
            pts = inputs
            cech = DelaunayComplex(points=pts, precision=self.precision)
            if self.filtration:
                fltr = "cech"

        elif self.input_type == "weighted point cloud":
            wgts = inputs[:, -1]
            pts = inputs[:, :-1]
            cech = DelaunayComplex(points=pts, weights=wgts, precision=self.precision)
            if self.filtration:
                fltr = "alpha"

        stree = cech.create_simplex_tree(max_alpha_square=self.max_alpha_square, filtration=fltr)

        persistence_dim_max = False
        # Specific case where, despite expansion(max_dimension), stree has a lower dimension
        if max_dimension > stree.dimension():
            persistence_dim_max = True

        stree.compute_persistence(
            homology_coeff_field=self.homology_coeff_field, persistence_dim_max=persistence_dim_max
        )

        return [stree.persistence_intervals_in_dimension(dim) for dim in self.dim_list_]

    def transform(self, X, Y=None):
        """Compute all the Čech complexes and their associated persistence diagrams.

        :param X: list of point clouds as Euclidean coordinates, plus the weight (as the last column value) if
            :paramref:`~gudhi.sklearn.cech_persistence.CechPersistence.input_type` was set to 'weighted point cloud'.
        :type X: list of list of float OR list of numpy.ndarray

        :return: Persistence diagrams in the format:

              - If `homology_dimensions` was set to `n`: `[array( Hn(X[0]) ), array( Hn(X[1]) ), ...]`
              - If `homology_dimensions` was set to `[i, j]`:
                `[[array( Hi(X[0]) ), array( Hj(X[0]) )], [array( Hi(X[1]) ), array( Hj(X[1]) )], ...]`
        :rtype: list of numpy ndarray of shape (,2) or list of list of numpy ndarray of shape (,2)
        """
        # Depends on homology_dimensions is an integer or a list of integer (else case)
        if isinstance(self.homology_dimensions, int):
            unwrap = True
            self.dim_list_ = [self.homology_dimensions]
        else:
            unwrap = False
            self.dim_list_ = self.homology_dimensions

        # threads is preferred as Rips construction and persistence computation releases the GIL
        res = Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(inputs) for inputs in X)

        if unwrap:
            res = [d[0] for d in res]
        return res

    def get_feature_names_out(self):
        """Provide column names for implementing sklearn's set_output API."""
        return [f"H{i}" for i in self.dim_list_]