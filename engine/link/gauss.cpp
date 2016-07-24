
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

#include "link/link.h"
#include "utilities/stringutils.h"
#include <iterator>

namespace regina {

Link* Link::fromOrientedGauss(const std::string& s) {
    std::vector<std::string> terms;
    basicTokenise(std::back_inserter(terms), s);
    return fromOrientedGauss(terms.begin(), terms.end());
}

bool Link::parseOrientedGaussTerm(const std::string& s,
        size_t nCross, size_t& crossing, int& strand, int& sign) {
    if (s.length() < 3)
        return false;

    if (s[0] == '+')
        strand = 1;
    else if (s[0] == '-')
        strand = 0;
    else
        return false;

    if (s[1] == '<')
        sign = (strand == 1 ? 1 : -1);
    else if (s[1] == '>')
        sign = (strand == 1 ? -1 : 1);
    else
        return false;

    if (! valueOf(s.substr(2), crossing))
        return false;

    return (crossing > 0 && crossing <= nCross);
}

bool Link::parseOrientedGaussTerm(const char* s,
        size_t nCross, size_t& crossing, int& strand, int& sign) {
    if (! (s[0] && s[1] && s[2]))
        return false;

    if (s[0] == '+')
        strand = 1;
    else if (s[0] == '-')
        strand = 0;
    else
        return false;

    if (s[1] == '<')
        sign = (strand == 1 ? 1 : -1);
    else if (s[1] == '>')
        sign = (strand == 1 ? -1 : 1);
    else
        return false;

    char* endPtr;
    crossing = static_cast<size_t>(strtol(s + 2, &endPtr, 10));
    return (*endPtr == 0 && crossing > 0 && crossing <= nCross);
}

} // namespace regina

