
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/*! \file triangulation/dim3/tetrahedron3.h
 *  \brief Internal header for tetrahedra in a 3-manifold triangulation.
 *
 *  This file is automatically included from triangulation/dim3.h; there is
 *  no need for end users to include this header explicitly.
 */

// We include triangulation3.h before the header guard, to ensure that the
// various face/component headers are processed in exactly the right order.
// This header will be re-included from within triangulation3.h.
#include "triangulation/dim3/triangulation3.h"

#ifndef __TETRAHEDRON3_H
#ifndef __DOXYGEN
#define __TETRAHEDRON3_H
#endif

namespace regina {

/**
 * \weakgroup triangulation
 * @{
 */

/**
 * Represents a tetrahedron within a 3-manifold triangulation.
 *
 * This is a specialisation of the generic Simplex class template; see the
 * generic Simplex documentation for an overview of how this class works.
 *
 * This 3-dimensional specialisation contains some extra functionality.
 * In particular, each tetrahedron stores additional details on how this
 * tetrahedron and its sub-faces integrate into the overall skeletal
 * structure of the triangulation.
 *
 * An implementation note: the Triangulation<3> class is responsible for
 * creating, maintaining and destroying this extra skeletal information.
 */
template <>
class REGINA_API Simplex<3> : public detail::SimplexBase<3> {
    private:
        int tmpOrientation_[4];
            /**< Temporary array used to represent orientations of triangles
                 when calculating orientability of boundary components.
                 Each orientation will be +/-1.
                 The array should only be used within these orientability
                 routines, and its contents afterwards are unpredictable. */

    public:
        /**
         * A dimension-specific alias for adjacentSimplex().
         *
         * See adjacentSimplex() for further information.
         */
        Simplex* adjacentTetrahedron(int face) const;
        /**
         * A dimension-specific alias for adjacentFacet().
         *
         * See adjacentFacet() for further information.
         */
        int adjacentFace(int face) const;

    private:
        /**
         * Creates a new tetrahedron with empty description and no
         * faces joined to anything.
         *
         * @param tri the triangulation to which the new tetrahedron belongs.
         */
        Simplex(Triangulation<3>* tri);
        /**
         * Creates a new tetrahedron with the given description and
         * no edges joined to anything.
         *
         * @param desc the description to give the new tetrahedron.
         * @param tri the triangulation to which the new tetrahedron belongs.
         */
        Simplex(const std::string& desc, Triangulation<3>* tri);

    friend class Triangulation<3>;
    friend class detail::TriangulationBase<3>;
        /**< Allow access to private members. */
};

/**
 * A convenience typedef for Simplex<3>.
 */
typedef Simplex<3> NTetrahedron;

/*@}*/

// Inline functions for Simplex<3>

inline Simplex<3>* Simplex<3>::adjacentTetrahedron(int face) const {
    return adjacentSimplex(face);
}

inline int Simplex<3>::adjacentFace(int face) const {
    return adjacentFacet(face);
}

inline Simplex<3>::Simplex(Triangulation<3>* tri) : detail::SimplexBase<3>(tri) {
}

inline Simplex<3>::Simplex(const std::string& desc, Triangulation<3>* tri) :
        detail::SimplexBase<3>(desc, tri) {
}

} // namespace regina

#endif

