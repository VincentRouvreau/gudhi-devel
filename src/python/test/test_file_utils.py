""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2022 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

import gudhi as gd
import numpy as np
import pytest
from tempfile import NamedTemporaryFile


def test_off_rw():
    for dim in range(2, 6):
        X = np.random.rand(123, dim)
        gd.write_points_to_off_file("rand.off", X)
        Y = gd.read_points_from_off_file("rand.off")
        assert Y == pytest.approx(X)


def test_human_off():
    pts = gd.read_points_from_off_file("human.off")
    # Should not try to read faces
    assert pts.shape == (4706, 3)


def test_invalid_off_file():
    name = NamedTemporaryFile().name
    with open(name, "w") as f:
        f.write("nOFL\n1000 200 2\n1.2 1.3 1.4")
    with pytest.raises(ValueError):
        gd.read_points_from_off_file(name)
    # Try to open a non-existing file - a new temp file name should not exist
    with pytest.raises(FileNotFoundError):
        gd.read_points_from_off_file(NamedTemporaryFile().name)


def test_off_file_header():
    name = NamedTemporaryFile().name
    print(name)
    # OFF
    ## classic case
    random_nb_points = np.random.randint(0, 100)
    with open(name, "w") as f:
        f.write(f"OFF\n{random_nb_points} 200 2\n1.2 1.3 1.4")
    with open(name) as f:
        assert gd._file_utils._read_off_file_header(f) == (3, random_nb_points)
    ## with comments
    random_nb_points = np.random.randint(0, 100)
    with open(name, "w") as f:
        f.write(
            f"# comment on the first line\nOFF\n# comment on the third line\n{random_nb_points} 200 2\n"
            "# comment before points\n1.2 1.3 1.4"
        )
    with open(name) as f:
        assert gd._file_utils._read_off_file_header(f) == (3, random_nb_points)
    # nOFF
    ## when 'dim nb_vertices nb_faces nb_edges' on the same line
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(f"nOFF\n{random_dim} {random_nb_points} 200 2\n1.2 1.3 1.4")
    with open(name) as f:
        assert gd._file_utils._read_off_file_header(f) == (random_dim, random_nb_points)
    ## when 'dim nb_vertices nb_faces nb_edges' on the same line + comments
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(
            f"# comment on the first line\nnOFF\n# comment on the third line\n{random_dim} {random_nb_points} 200 2\n"
            "# comment before points\n1.2 1.3 1.4"
        )
    with open(name) as f:
        assert gd._file_utils._read_off_file_header(f) == (random_dim, random_nb_points)
    ## when 'dim' and 'nb_vertices nb_faces nb_edges' on the separated lines
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(f"nOFF\n{random_dim}\n{random_nb_points} 200 2\n1.2 1.3 1.4")
    with open(name) as f:
        assert gd._file_utils._read_off_file_header(f) == (random_dim, random_nb_points)
    ## when 'dim' and 'nb_vertices nb_faces nb_edges' on the separated lines + comments
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(
            f"# first comment\nnOFF\n# second comment\n{random_dim}\n# third comment\n{random_nb_points} 200 2\n"
            "# another comment\n1.2 1.3 1.4"
        )
    with open(name) as f:
        assert gd._file_utils._read_off_file_header(f) == (random_dim, random_nb_points)

def test_get_next_line():
    for end_of_line in ['\r\n', '\n']:
        for comment_char in ['#', '@']:
            name = NamedTemporaryFile().name
            print(name)
            with open(name, "w") as f:
                f.write(f"{comment_char}{end_of_line}{end_of_line}{comment_char}{end_of_line}First{end_of_line}")
                f.write(f"{end_of_line}{comment_char}{end_of_line}{end_of_line}Second{end_of_line}{end_of_line}")
            with open(name) as f:
                assert gd._file_utils._get_next_line(f, comment=comment_char).split() == ["First"]
                assert gd._file_utils._get_next_line(f, comment=comment_char).split() == ["Second"]

def test_non_existing_csv_file():
    # Try to open a non existing file
    matrix = gd.read_lower_triangular_matrix_from_csv_file(
        csv_file=NamedTemporaryFile().name
    )
    assert matrix == []


def test_full_square_distance_matrix_csv_file():
    # Create test file
    test_file = open("full_square_distance_matrix.csv", "w")
    test_file.write("0;1;2;3;\n1;0;4;5;\n2;4;0;6;\n3;5;6;0;")
    test_file.close()
    matrix = gd.read_lower_triangular_matrix_from_csv_file(
        csv_file="full_square_distance_matrix.csv", separator=";"
    )
    assert matrix == [[], [1.0], [2.0, 4.0], [3.0, 5.0, 6.0]]


def test_lower_triangular_distance_matrix_csv_file():
    # Create test file
    test_file = open("lower_triangular_distance_matrix.csv", "w")
    test_file.write("\n1,\n2,3,\n4,5,6,\n7,8,9,10,")
    test_file.close()
    matrix = gd.read_lower_triangular_matrix_from_csv_file(
        csv_file="lower_triangular_distance_matrix.csv", separator=","
    )
    assert matrix == [[], [1.0], [2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0, 10.0]]


def test_non_existing_persistence_file():
    # Try to open a non existing file
    persistence = gd.read_persistence_intervals_grouped_by_dimension(
        persistence_file=NamedTemporaryFile().name
    )
    assert persistence == []
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file=NamedTemporaryFile().name, only_this_dim=1
    )
    np.testing.assert_array_equal(persistence, [])


def test_read_persistence_intervals_without_dimension():
    # Create test file
    test_file = open("persistence_intervals_without_dimension.pers", "w")
    test_file.write(
        "# Simple persistence diagram without dimension\n2.7 3.7\n9.6 14.\n34.2 34.974\n3. inf"
    )
    test_file.close()
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_without_dimension.pers"
    )
    np.testing.assert_array_equal(
        persistence, [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float("Inf"))]
    )
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_without_dimension.pers", only_this_dim=0
    )
    np.testing.assert_array_equal(persistence, [])
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_without_dimension.pers", only_this_dim=1
    )
    np.testing.assert_array_equal(persistence, [])
    persistence = gd.read_persistence_intervals_grouped_by_dimension(
        persistence_file="persistence_intervals_without_dimension.pers"
    )
    assert persistence == {
        -1: [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float("Inf"))]
    }


def test_read_persistence_intervals_with_dimension():
    # Create test file
    test_file = open("persistence_intervals_with_dimension.pers", "w")
    test_file.write(
        "# Simple persistence diagram with dimension\n0 2.7 3.7\n1 9.6 14.\n3 34.2 34.974\n1 3. inf"
    )
    test_file.close()
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers"
    )
    np.testing.assert_array_equal(
        persistence, [(2.7, 3.7), (9.6, 14.0), (34.2, 34.974), (3.0, float("Inf"))]
    )
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=0
    )
    np.testing.assert_array_equal(persistence, [(2.7, 3.7)])
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=1
    )
    np.testing.assert_array_equal(persistence, [(9.6, 14.0), (3.0, float("Inf"))])
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=2
    )
    np.testing.assert_array_equal(persistence, [])
    persistence = gd.read_persistence_intervals_in_dimension(
        persistence_file="persistence_intervals_with_dimension.pers", only_this_dim=3
    )
    np.testing.assert_array_equal(persistence, [(34.2, 34.974)])
    persistence = gd.read_persistence_intervals_grouped_by_dimension(
        persistence_file="persistence_intervals_with_dimension.pers"
    )
    assert persistence == {
        0: [(2.7, 3.7)],
        1: [(9.6, 14.0), (3.0, float("Inf"))],
        3: [(34.2, 34.974)],
    }
