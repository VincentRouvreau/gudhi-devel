/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
  
/** \brief Concept of the template parameter for the classes
 * `Gudhi::Hasse_diagram::Hasse_diagram_persistence<HasseDiagramOptions>` and
 * `Gudhi::Hasse_diagram::Hasse_diagram<HasseDiagramOptions>`.
 *
 * One model for this is `Gudhi::Hasse_diagram::Hasse_diagram_options_full_featured`.
 * If you want to provide your own, it is recommended that you derive from it and override some parts instead of
 * writing a class from scratch.
 */
struct HasseDiagramOptions {
  /// Must be an integer type.
  typedef IncidenceType Incidence_type;
  /// Must be comparable with operator<.
  typedef FiltrationValue Filtration_value;
  /// Must be an integer type.
  typedef SimplexKey Simplex_key;
  /// Can be of any type to store additional information
  typedef AdditionalInformation Additional_information;
  /// If true, each cell has extra storage for one `Simplex_key`. Necessary for `Persistent_cohomology`.
  static const bool store_key;
  /// If true, each cell has extra storage for one `Filtration_value`. Without it, `Persistent_cohomology` degenerates to computing usual (non-persistent) cohomology.
  static const bool store_filtration;
  /// If true, each cell has extra storage for one `Additional_information`, and this value is propagated by operations like `Gudhi::Simplex_tree::expansion`. Without it, `Persistent_cohomology` degenerates to computing usual (non-persistent) cohomology.
  static const bool store_additional_information = true;
};
