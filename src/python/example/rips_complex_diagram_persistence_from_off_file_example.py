#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__license__ = "MIT"


import argparse
import errno
import os
import gudhi as gd


parser = argparse.ArgumentParser(
    description="Vietoris-Rips complex creation from " "points read in a OFF file.",
    epilog="Example: "
    "example/rips_complex_diagram_persistence_from_off_file_example.py "
    "-f ../data/points/tore3D_300.off -a 0.6"
    "- Constructs a Rips complex with the "
    "points from the given OFF file.",
)
parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-e", "--max_edge_length", type=float, default=0.5)
parser.add_argument("-d", "--max_dimension", type=int, default=1)
parser.add_argument("-b", "--band", type=float, default=0.0)
parser.add_argument(
    "--no-diagram",
    default=False,
    action="store_true",
    help="Flag for not to display the diagrams",
)

args = parser.parse_args()

with open(args.file) as f:
    first_line = f.readline()
    if (first_line == "OFF\n") or (first_line == "nOFF\n"):
        print("##############################################################")
        print("Vietoris-Rips complex creation from points read in a OFF file")

        message = "Vietoris-Rips complex with max_edge_length=" + repr(args.max_edge_length)
        print(message)

        point_cloud = gd.read_points_from_off_file(off_file=args.file)
        cplx = gd.filtrations.rips_complex(
            points=point_cloud, max_edge_length=args.max_edge_length, max_dimension=args.max_dimension
        )

        print(f"Number of simplices={cplx.num_simplices()}")

        diag = cplx.persistence()

        print(f"betti_numbers()={cplx.betti_numbers()}")

        if args.no_diagram == False:
            import matplotlib.pyplot as plot

            gd.plot_persistence_diagram(diag, band=args.band)
            plot.show()
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.file)

    f.close()
