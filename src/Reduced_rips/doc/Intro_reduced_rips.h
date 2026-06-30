/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *    - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_REDUCED_RIPS_INTRO_REDUCED_RIPS_H_
#define DOC_REDUCED_RIPS_INTRO_REDUCED_RIPS_H_

namespace Gudhi {

namespace reduced_rips {

/**  @defgroup reduced_rips Reduced Vietoris-Rips degree-1 persistent homology
 *
 * @author    Thomas Burnett, Musashi Koyama
 *
 * @{
 *
 * @section reducedripsdefinition Definition
 *
 * This module computes the degree-1 (i.e. @f$H_1@f$) Vietoris-Rips persistent homology of a point cloud or a
 * distance matrix, without ever building the full Vietoris-Rips complex. It implements the <em>Reduced
 * Vietoris-Rips filtration</em> of Koyama, M&eacute;moli, Robins and Turner @cite
 * koyama2026computationdegree1persistenthomology.
 *
 * The reduction works on a small part of the complex. The 1-skeleton that matters for @f$H_1@f$ is the
 * relative neighborhood graph, and the only 2-simplices that enter the boundary matrix are the ones a lune
 * test certifies. This lets it handle much larger inputs than building the full complex would.
 *
 * @section reducedripsinput Input
 *
 * Two inputs are accepted:
 * - a <b>Euclidean point cloud</b>, for which the relative neighborhood graph is built from a Delaunay
 *   triangulation in dimension 2 and 3 and a direct @f$O(n^2)@f$ construction otherwise.
 * - an <b>arbitrary symmetric distance matrix</b> (Gudhi::reduced_rips::Reduced_rips::from_distance_matrix),
 *   for which the same reduction is driven purely by the supplied distances.
 *
 * The reduced filtration is exact for any symmetric matrix of non-negative dissimilarities.
 *
 * @section reducedripsscope Scope and limitations
 *
 * - Only homological dimension 1 is computed.
 *
 * @section reducedripsexample Example
 *
 * @include Reduced_rips/example_reduced_rips_from_off.cpp
 *
 * @}
 */

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // DOC_REDUCED_RIPS_INTRO_REDUCED_RIPS_H_
