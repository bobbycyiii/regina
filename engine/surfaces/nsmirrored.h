
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

/*! \file surfaces/nsmirrored.h
 *  \brief Provides a normal surface vector that is mirrored in another
 *  coordinate system to avoid frequent lengthy calculations.
 */

#ifndef __NSMIRRORED_H
#ifndef __DOXYGEN
#define __NSMIRRORED_H
#endif

#include "regina-core.h"
#include "surfaces/nnormalsurface.h"

namespace regina {

/**
 * \weakgroup surfaces
 * @{
 */

/**
 * A normal surface vector that is mirrored in another coordinate system
 * to avoid frequent lengthy calculations.  When it is difficult to
 * convert from the native coordinate system to standard tri-quad-oct
 * coordinates, use this as a base class.  The conversion of the entire
 * vector will be done once only, and future coordinate lookups will be
 * performed through the pre-converted mirror vector.
 *
 * Subclasses need not implement any of the coordinate lookup routines.
 * The default implementation will be to pass the lookup to the mirror.
 * If any particular lookup can be done more efficiently in the native
 * coordinate system, the corresponding routine should be overridden.
 *
 * Subclasses must however implement two variants of makeMirror():
 *
 * - a static function
 *   NNormalSurfaceVector* makeMirror(const Ray&, const NTriangulation*),
 *   which builds the mirror vector from the native vector;
 * - a virtual function
 *   NNormalSurfaceVector* makeMirror(const NTriangulation*) const,
 *   which simply calls the static function described above.
 *
 * \ifacespython Not present.
 */
class REGINA_API NNormalSurfaceVectorMirrored : public NNormalSurfaceVector {
    private:
        NNormalSurfaceVector* mirror;
            /**< The mirror vector. */

    public:
        /**
         * Creates a new vector all of whose entries are initialised to
         * zero.
         *
         * @param length the number of elements in the new vector.
         */
        NNormalSurfaceVectorMirrored(size_t length);
        /**
         * Creates a new vector that is a clone of the given vector.
         *
         * @param cloneMe the vector to clone.
         */
        NNormalSurfaceVectorMirrored(const Vector<LargeInteger>& cloneMe);
        /**
         * Creates a new vector that is a clone of the given vector.
         *
         * @param cloneMe the vector to clone.
         */
        NNormalSurfaceVectorMirrored(const NNormalSurfaceVectorMirrored&
            cloneMe);
        /**
         * Destroys this vector and its mirror if appropriate.
         */
        virtual ~NNormalSurfaceVectorMirrored();

        /**
         * Creates a new mirror vector corresponding to this vector.
         *
         * This function must return a newly created vector that describes
         * this normal surface in the corresponding "standard" coordinate
         * system (i.e., NNormalSurfaceVectorStandard if this subclass is
         * NNormalSurfaceVectorQuad, or NNormalSurfaceVectorANStandard if
         * this subclass is NNormalSurfaceVectorQuadOct, and so on.)
         *
         * Typically subclasses should implement this function by calling the
         * corresponding static makeMirror() function, which does the real work.
         *
         * @param triang the triangulation in which this normal surface lives.
         * @return a newly created vector that describes the same normal
         * surface in a "standard" coordinate system, as described above.
         */
        virtual NNormalSurfaceVector* makeMirror(const NTriangulation* triang)
            const = 0;

#ifdef __DOXYGEN
        /**
         * Creates a new mirror vector corresponding to the given vector.
         * This function should be implemented by each subclass of
         * NNormalSurfaceVectorMirrored.
         *
         * This function must return a newly created vector that describes
         * the given normal surface in a different coordinate system:
         *
         * - The given vector must be in the coordinate system described
         *   by the subclass in which this function is implemented.
         *
         * - The return vector must be in the corresponding "standard"
         *   coordinate system (i.e., NNormalSurfaceVectorStandard if
         *   this subclass is NNormalSurfaceVectorQuad, or
         *   NNormalSurfaceVectorANStandard if this subclass is
         *   NNormalSurfaceVectorQuadOct, and so on.)
         *
         * @param original a vector in the coordinate system corresponding to
         * the subclass in which this function is implemented.
         * @param triang the triangulation in which the corresponding
         * normal surface lives.
         * @return a newly created vector that describes the same normal
         * surface in a "standard" coordinate system, as described above.
         */
        static NNormalSurfaceVector* makeMirror(const Ray& original,
            const NTriangulation* triang);
#endif

        virtual void setElement(size_t index, const LargeInteger& value)
            override;
        virtual void operator += (const NNormalSurfaceVector& other) override;
        virtual void scaleDown() override;

        virtual LargeInteger triangles(size_t tetIndex,
            int vertex, const NTriangulation* triang) const override;
        virtual LargeInteger orientedTriangles(size_t tetIndex,
            int vertex, const NTriangulation* triang, bool orientation) const
            override;
        virtual LargeInteger quads(size_t tetIndex,
            int quadType, const NTriangulation* triang) const override;
        virtual LargeInteger orientedQuads(size_t tetIndex,
            int quadType, const NTriangulation* triang, bool orientation) const
            override;
        virtual LargeInteger octs(size_t tetIndex,
            int octType, const NTriangulation* triang) const override;
        virtual LargeInteger edgeWeight(size_t edgeIndex,
            const NTriangulation* triang) const override;
        virtual LargeInteger arcs(size_t triIndex,
            int triVertex, const NTriangulation* triang) const override;
};

/*@}*/

// Inline functions for NNormalSurfaceVectorMirrored

inline NNormalSurfaceVectorMirrored::NNormalSurfaceVectorMirrored(
        size_t length) : NNormalSurfaceVector(length), mirror(0) {
}
inline NNormalSurfaceVectorMirrored::NNormalSurfaceVectorMirrored(
        const Vector<LargeInteger>& cloneMe) :
        NNormalSurfaceVector(cloneMe), mirror(0) {
}
inline NNormalSurfaceVectorMirrored::NNormalSurfaceVectorMirrored(
        const NNormalSurfaceVectorMirrored& cloneMe) :
        NNormalSurfaceVector(cloneMe), mirror(0) {
}
inline NNormalSurfaceVectorMirrored::~NNormalSurfaceVectorMirrored() {
    delete mirror;
}

inline void NNormalSurfaceVectorMirrored::setElement(size_t index,
        const LargeInteger& value) {
    coords_.setElement(index, value);
    delete mirror;
}

inline void NNormalSurfaceVectorMirrored::operator += (
        const NNormalSurfaceVector& other) {
    coords_ += other.coords();
    delete mirror;
}

inline void NNormalSurfaceVectorMirrored::scaleDown() {
    coords_.scaleDown();
    delete mirror;
}

inline LargeInteger NNormalSurfaceVectorMirrored::triangles(
        size_t tetIndex, int vertex, const NTriangulation* triang) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->triangles(tetIndex, vertex, triang);
}
inline LargeInteger NNormalSurfaceVectorMirrored::orientedTriangles(
        size_t tetIndex, int vertex, const NTriangulation* triang,
        bool orientation) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->orientedTriangles(tetIndex, vertex, triang,
        orientation);
}
inline LargeInteger NNormalSurfaceVectorMirrored::quads(
        size_t tetIndex, int quadType, const NTriangulation* triang) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->quads(tetIndex, quadType, triang);
}
inline LargeInteger NNormalSurfaceVectorMirrored::orientedQuads(
        size_t tetIndex, int quadType, const NTriangulation* triang,
        bool orientation) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->orientedQuads(tetIndex, quadType, triang,
        orientation);
}
inline LargeInteger NNormalSurfaceVectorMirrored::octs(
        size_t tetIndex, int octType, const NTriangulation* triang) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->octs(tetIndex, octType, triang);
}
inline LargeInteger NNormalSurfaceVectorMirrored::edgeWeight(
        size_t edgeIndex, const NTriangulation* triang) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->edgeWeight(edgeIndex, triang);
}
inline LargeInteger NNormalSurfaceVectorMirrored::arcs(
        size_t triIndex, int triVertex,
        const NTriangulation* triang) const {
    if (! mirror)
        const_cast<NNormalSurfaceVectorMirrored*>(this)->mirror =
            makeMirror(triang);
    return mirror->arcs(triIndex, triVertex, triang);
}

} // namespace regina

#endif

