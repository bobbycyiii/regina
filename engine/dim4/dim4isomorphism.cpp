
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2011, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
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

/* end stub */

#include "dim4/dim4isomorphism.h"
#include "dim4/dim4triangulation.h"
#include "triangulation/ngeneralisomorphism.tcc"

namespace regina {

// Instatiate all templates from the .tcc file.
template void NGeneralIsomorphism<4>::writeTextShort(std::ostream&) const;
template void NGeneralIsomorphism<4>::writeTextLong(std::ostream&) const;
template bool NGeneralIsomorphism<4>::isIdentity() const;
template NGeneralIsomorphism<4>::NGeneralIsomorphism(
    const NGeneralIsomorphism<4>&);
template Dim4Isomorphism* NGeneralIsomorphism<4>::random(unsigned);
template Dim4Triangulation* NGeneralIsomorphism<4>::apply(
        const Dim4Triangulation*) const;
template void NGeneralIsomorphism<4>::applyInPlace(Dim4Triangulation*) const;

} // namespace regina

