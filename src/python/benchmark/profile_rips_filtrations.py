from gudhi.datasets.generators import points
from gudhi import read_lower_triangular_matrix_from_csv_file
from gudhi.filtrations import rips_complex
from scipy.spatial.distance import cdist
import numpy as np
from timeit import default_timer as timer
from datetime import timedelta

pts = points.sphere(n_samples = 5000, ambient_dim = 3, radius = 1)
start = timer(); d = cdist(pts, pts); print(timedelta(seconds=timer()-start))
# 0:00:00.079720

f = '5000pts_on_2sphere.csv'
delimiter = ','
np.savetxt(f, d, delimiter=delimiter)
# To get the sequence format that is also accepted
sd = read_lower_triangular_matrix_from_csv_file(f, separator=delimiter)

start = timer(); st1 = rips_complex(points=pts, max_edge_length=float('inf')); print(timedelta(seconds=timer()-start))
# 0:00:01.056092
start = timer(); st2 = rips_complex(distance_matrix=d, max_edge_length=float('inf')); print(timedelta(seconds=timer()-start))
# 0:00:00.116341
start = timer(); st3 = rips_complex(distance_matrix=sd, max_edge_length=float('inf')); print(timedelta(seconds=timer()-start))
# 0:00:01.086323

max_edge_length = 1.5 # almost 0.75 * np.max(d)
start = timer(); st1 = rips_complex(points=pts, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.687482
start = timer(); st2 = rips_complex(distance_matrix=d, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.193933
start = timer(); st3 = rips_complex(distance_matrix=sd, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.729296

max_edge_length = 1. # almost 0.5 * np.max(d)
start = timer(); st1 = rips_complex(points=pts, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.343175
start = timer(); st2 = rips_complex(distance_matrix=d, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.201040
start = timer(); st3 = rips_complex(distance_matrix=sd, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.393524

max_edge_length = 0.5 # almost 0.25 * np.max(d)
start = timer(); st1 = rips_complex(points=pts, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.111799
start = timer(); st2 = rips_complex(distance_matrix=d, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.151324
start = timer(); st3 = rips_complex(distance_matrix=sd, max_edge_length=max_edge_length); print(timedelta(seconds=timer()-start))
# 0:00:00.161206
