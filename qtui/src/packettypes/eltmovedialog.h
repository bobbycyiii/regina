
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
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

/* end stub */

/*! \file eltmovedialog.h
 *  \brief Provides a dialog through which the user can select and
 *  perform an elementary move.
 */

#ifndef __ELTMOVEDIALOG_H
#define __ELTMOVEDIALOG_H

#include <QDialog>
#include <vector>

class QButtonGroup;
class QRadioButton;
class TetrahedronChooser;
class FaceChooser;
class EdgeChooser;
class EdgeIntChooser;
class VertexChooser;

namespace regina {
    class NTriangulation;
};

/**
 * A dialog used to select and perform an elementary move on a
 * triangulation.
 */
class EltMoveDialog : public QDialog {
    Q_OBJECT

    private:
        /**
         * Internal components:
         */
        EdgeChooser* box32;
        FaceChooser* box23;
        EdgeIntChooser* box44;
        EdgeChooser* box20e;
        VertexChooser* box20v;
        EdgeIntChooser* box21;
        FaceChooser* boxOpenBook;
        EdgeChooser* boxCloseBook;
        TetrahedronChooser* boxShellBdry;
        EdgeChooser* boxCollapseEdge;
        QRadioButton* use32;
        QRadioButton* use23;
        QRadioButton* use44;
        QRadioButton* use20e;
        QRadioButton* use20v;
        QRadioButton* use21;
        QRadioButton* useOpenBook;
        QRadioButton* useCloseBook;
        QRadioButton* useShellBdry;
        QRadioButton* useCollapseEdge;
        QButtonGroup* moveTypes;

        /**
         * Packet tree structure:
         */
        regina::NTriangulation* tri;

    public:
        /**
         * Constructor and destructor.
         */
        EltMoveDialog(QWidget* parent, regina::NTriangulation* useTri);
        ~EltMoveDialog();

    protected slots:
        /**
         * Ok has been clicked.
         */
        virtual void slotOk();
};

#endif
