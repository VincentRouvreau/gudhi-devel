#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__license__ = "MIT"


from gudhi.filtrations import rips_complex


print("#####################################################################")
print("Vietoris-Rips complex creation from points")
cplx = rips_complex(points=[[0, 0], [1, 0], [0, 1], [1, 1]], max_edge_length=42, max_dimension=1)

print("filtrations=")
for simplex_with_filtration in cplx.get_filtration():
    print("(%s, %.2f)" % tuple(simplex_with_filtration))

print("star([0])=", cplx.get_star([0]))
print("coface([0], 1)=", cplx.get_cofaces([0], 1))
