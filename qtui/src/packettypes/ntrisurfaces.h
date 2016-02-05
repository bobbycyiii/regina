
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2014, Ben Burton                                   *
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

/*! \file ntrisurfaces.h
 *  \brief Provides a normal surface properties viewer for triangulations.
 */

#ifndef __NTRISURFACES_H
#define __NTRISURFACES_H

#include "../packettabui.h"

class QAbstractButton;
class QLabel;

namespace regina {
    class NPacket;
    template <int> class Triangulation;
    typedef Triangulation<3> NTriangulation;
};

/**
 * A triangulation page for viewing normal surface properties.
 */
class NTriSurfacesUI : public QObject, public PacketViewerTab {
    Q_OBJECT

    private:
        /**
         * Packet details
         */
        regina::NTriangulation* tri;

        /**
         * Internal components
         */
        QWidget* ui;
        QLabel* zeroEff;
        QLabel* titleZeroEff;
        QLabel* splitting;
        QLabel* titleSplitting;
        QLabel* threeSphere;
        QLabel* titleThreeSphere;
        QLabel* threeBall;
        QLabel* titleThreeBall;
        QLabel* solidTorus;
        QLabel* titleSolidTorus;
        QLabel* irreducible;
        QLabel* titleIrreducible;
        QLabel* haken;
        QLabel* titleHaken;
        QLabel* strict;
        QLabel* titleStrict;
        QLabel* hyperbolic;
        QLabel* titleHyperbolic;
        QAbstractButton* btnZeroEff;
        QAbstractButton* btnSplitting;
        QAbstractButton* btnThreeSphere;
        QAbstractButton* btnThreeBall;
        QAbstractButton* btnSolidTorus;
        QAbstractButton* btnIrreducible;
        QAbstractButton* btnHaken;
        QAbstractButton* btnStrict;
        QLabel* manifold;
        QLabel* census;

    public:
        /**
         * Constructor.
         */
        NTriSurfacesUI(regina::NTriangulation* packet,
            PacketTabbedUI* useParentUI);

        /**
         * PacketViewerTab overrides.
         */
        regina::NPacket* getPacket();
        QWidget* getInterface();
        void refresh();

    public slots:
        /**
         * Calculate properties.
         */
        void calculateZeroEff();
        void calculateSplitting();
        void calculateThreeSphere();
        void calculateThreeBall();
        void calculateSolidTorus();
        void calculateIrreducible();
        void calculateHaken();
        void calculateStrict();

        /**
         * Notify that preferences have changed.
         */
        void updatePreferences();
};

inline void NTriSurfacesUI::updatePreferences() {
    refresh();
}

#endif
