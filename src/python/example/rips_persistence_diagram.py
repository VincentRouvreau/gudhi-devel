#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__license__ = "MIT"


import matplotlib.pyplot as plot
import gudhi as gd


print("#####################################################################")
print("Vietoris-Rips complex creation from points")
cplx = gd.filtrations.rips_complex(points=[[0, 0], [1, 0], [0, 1], [1, 1]], max_edge_length=42, max_dimension=1)

diag = cplx.persistence(homology_coeff_field=2, min_persistence=0)
print(f"diag={diag}")

gd.plot_persistence_diagram(diag)
plot.show()
