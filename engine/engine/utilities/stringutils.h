
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2002, Ben Burton                                   *
 *  For further details contact Ben Burton (benb@acm.org).                *
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
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,        *
 *  MA 02111-1307, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

/*! \file stringutils.h
 *  \brief Provides various routines for use with C++ strings.
 */

#ifndef __STRINGUTILS_H
#ifndef __DOXYGEN
#define __STRINGUTILS_H
#endif

#include <string>

namespace regina {

/**
 * Creates a new C string that is a duplicate of the given C++ string.
 *
 * The deallocation of the new C string is the responsibility of
 * the caller of this routine.
 *
 * \ifaces Not present.
 *
 * @param str the C++ string to duplicate.
 * @return the new duplicate C string.
 */
char* duplicate(const std::string& str);

/**
 * Decomposes the given string into tokens.
 * This is an extremely simple tokeniser; tokens are defined to be
 * separated by whitespace.
 *
 * @param results the output iterator to which the resulting tokens will
 * be written; this must accept objects of type <tt>const std::string&</tt>.
 * @param str the string to decompose.
 * @return the number of tokens found.
 */
template <class OutputIterator>
unsigned basicTokenise(OutputIterator results, const std::string& str);

} // namespace regina

// Template definitions

#include "utilities/stringutils.tcc"

#endif

