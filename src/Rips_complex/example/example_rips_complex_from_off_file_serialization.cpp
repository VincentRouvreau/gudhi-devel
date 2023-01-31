#include <gudhi/Rips_complex.h>
// to construct Rips_complex from a OFF file of points
#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Clock.h>

// include headers that implement a archive in simple text format
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off threshold archive_file.txt\n";
  std::cerr << "       i.e.: " << progName << " ../../data/points/alphacomplexdoc.off 60.0 archive.txt\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 4) usage(argc, (argv[0] - 1));

  std::string off_file_name(argv[1]);
  double threshold = atof(argv[2]);
  std::string archive_file_name(argv[3]);

  // Type definitions
  using Point = std::vector<double>;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<double>;

  // ----------------------------------------------------------------------------
  // Init of a Rips complex from an OFF file
  // ----------------------------------------------------------------------------
  Gudhi::Points_off_reader<Point> off_reader(off_file_name);
  Gudhi::Clock rips_clock("Rips ctor");
  rips_clock.begin();
  Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), threshold, Gudhi::Euclidean_distance());
  rips_clock.end();
  std::clog << rips_clock;

  Gudhi::Clock serial_clock("serialize Rips");
  serial_clock.begin();
  {
    std::ofstream ofs(archive_file_name);
    boost::archive::binary_oarchive oarchive(ofs);
    oarchive << rips_complex_from_file;
  }
  serial_clock.end();
  std::clog << serial_clock;

  std::vector<std::vector<double>> empty{};
  Rips_complex rips_complex_from_archive(empty, 0.);

  Gudhi::Clock deserial_clock("deserialize Rips");
  deserial_clock.begin();
  {
    std::ifstream ifs(archive_file_name);
    boost::archive::binary_iarchive iarchive(ifs);
    iarchive >> rips_complex_from_archive;
  }
  deserial_clock.end();
  std::clog << deserial_clock;

  return 0;
}

// 4000 vertices on 3-torus - Rips ctor: 0.262s - serialize Rips: 0.801s - deserialize Rips: 0.673s - 123 Mo
// 8000 vertices on 3-torus - Rips ctor: 1.032s - serialize Rips: 3.649s - deserialize Rips: 2.442s - 489 Mo
