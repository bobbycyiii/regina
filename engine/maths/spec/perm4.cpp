
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2021, Ben Burton                                   *
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

#include <sstream>
#include "maths/perm.h"

namespace regina {

Perm<4>::Perm(int a0, int a1, int b0, int b1, int c0, int c1, int d0, int d1) {
    int image[4];
    image[a0] = a1;
    image[b0] = b1;
    image[c0] = c1;
    image[d0] = d1;
    code_ = static_cast<Code>(S4Index(image[0], image[1], image[2], image[3]));
}

Perm<4>::Perm(const int* a, const int* b) {
    int image[4];
    image[a[0]] = b[0];
    image[a[1]] = b[1];
    image[a[2]] = b[2];
    image[a[3]] = b[3];
    code_ = static_cast<Code>(S4Index(image[0], image[1], image[2], image[3]));
}

std::string Perm<4>::str() const {
    char ans[5];
    for (int i = 0; i < 4; i++)
        ans[i] = static_cast<char>('0' + imageTable[code_][i]);
    ans[4] = 0;

    return ans;
}

std::string Perm<4>::trunc(unsigned len) const {
    char ans[5];
    for (unsigned i = 0; i < len; ++i)
        ans[i] = static_cast<char>('0' + imageTable[code_][i]);
    ans[len] = 0;
    return ans;
}

std::string Perm<4>::trunc2() const {
    char ans[3];
    ans[0] = static_cast<char>('0' + imageTable[code_][0]);
    ans[1] = static_cast<char>('0' + imageTable[code_][1]);
    ans[2] = 0;
    return ans;
}

std::string Perm<4>::trunc3() const {
    char ans[4];
    ans[0] = static_cast<char>('0' + imageTable[code_][0]);
    ans[1] = static_cast<char>('0' + imageTable[code_][1]);
    ans[2] = static_cast<char>('0' + imageTable[code_][2]);
    ans[3] = 0;
    return ans;
}

} // namespace regina

