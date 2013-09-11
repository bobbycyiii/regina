
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

#include "triangulation/ntriangulation.h"

#include "eltmovedialog.h"
#include "reginasupport.h"
#include "choosers/edgechooser.h"
#include "choosers/edgeintchooser.h"
#include "choosers/facechooser.h"
#include "choosers/tetrahedronchooser.h"
#include "choosers/vertexchooser.h"

#include <QButtonGroup>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QLayout>
#include <QRadioButton>
#include <QWhatsThis>

namespace {
    const int ID_32 = 0;
    const int ID_23 = 1;
    const int ID_44 = 2;
    const int ID_20E = 3;
    const int ID_20V = 4;
    const int ID_21 = 5;
    const int ID_OPENBOOK = 6;
    const int ID_CLOSEBOOK = 7;
    const int ID_SHELLBDRY = 8;
    const int ID_COLLAPSEEDGE = 9;

    bool has20v(regina::NVertex* v) {
        return v->getTriangulation()->twoZeroMove(v, true, false);
    }

    bool has32(regina::NEdge* e) {
        return e->getTriangulation()->threeTwoMove(e, true, false);
    }

    bool has20e(regina::NEdge* e) {
        return e->getTriangulation()->twoZeroMove(e, true, false);
    }

    bool hasCloseBook(regina::NEdge* e) {
        return e->getTriangulation()->closeBook(e, true, false);
    }

    bool hasCollapseEdge(regina::NEdge* e) {
        return e->getTriangulation()->collapseEdge(e, true, false);
    }

    bool has44(regina::NEdge* e, int axis) {
        return e->getTriangulation()->fourFourMove(e, axis, true, false);
    }

    bool has21(regina::NEdge* e, int end) {
        return e->getTriangulation()->twoOneMove(e, end, true, false);
    }

    bool has23(regina::NFace* f) {
        return f->getTriangulation()->twoThreeMove(f, true, false);
    }

    bool hasOpenBook(regina::NFace* f) {
        return f->getTriangulation()->openBook(f, true, false);
    }

    bool hasShellBoundary(regina::NTetrahedron* t) {
        return t->getTriangulation()->shellBoundary(t, true, false);
    }
}

EltMoveDialog::EltMoveDialog(QWidget* parent, regina::NTriangulation* useTri) :
        QDialog(parent), // tr("Elementary Move"), Ok|Cancel, Ok, parent),
        tri(useTri) {
    setWindowTitle(tr("Elementary Move"));
    QVBoxLayout *dialogLayout = new QVBoxLayout(this);

    QGridLayout* layout = new QGridLayout();
      //, 10, 2, 0 /* margin */, spacingHint());
    dialogLayout->addLayout(layout);

    use32 = new QRadioButton(tr("&3-2"), this);
    use32->setWhatsThis( tr("<qt>Perform a 3-2 move on this "
        "triangulation.<p>"
        "A <i>3-2 move</i> involves replacing three tetrahedra joined along "
        "an edge of degree three with two tetrahedra joined along a "
        "single face.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(use32, 0, 0);
    use23 = new QRadioButton(tr("&2-3"), this);
    use23->setWhatsThis( tr("<qt>Perform a 2-3 move on this "
        "triangulation.<p>"
        "A <i>2-3 move</i> involves replacing two tetrahedra joined along "
        "a single face with three tetrahedra joined along an edge of "
        "degree three.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(use23, 1, 0);
    use44 = new QRadioButton(tr("&4-4"), this);
    use44->setWhatsThis( tr("<qt>Perform a 4-4 move on this "
        "triangulation.<p>"
        "A <i>4-4 move</i> involves replacing four tetrahedra joined along "
        "an edge of degree four with four new tetrahedra joined along a "
        "different edge in a different position.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(use44, 2, 0);
    use20e = new QRadioButton(tr("2-0 (&edge)"), this);
    use20e->setWhatsThis( tr("<qt>Perform a 2-0 edge move on this "
        "triangulation.<p>"
        "A <i>2-0 edge move</i> involves taking two tetrahedra joined along "
        "an edge of degree two and squashing them flat.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(use20e, 3, 0);
    use20v = new QRadioButton(tr("2-0 (&vertex)"), this);
    use20v->setWhatsThis(tr("<qt>Perform a 2-0 vertex move on this "
        "triangulation.<p>"
        "A <i>2-0 vertex move</i> involves taking two tetrahedra meeting at "
        "a vertex of degree two and squashing them together.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(use20v, 4, 0);
    use21 = new QRadioButton(tr("2-&1"), this);
    use21->setWhatsThis( tr("<qt>Perform a 2-1 move on this "
        "triangulation.<p>"
        "A <i>2-1 move</i> involves taking a tetrahedron joined to itself "
        "about an edge of degree one and merging it with an adjacent "
        "tetrahedron.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(use21, 5, 0);
    useOpenBook = new QRadioButton(tr("&Open book"), this);
    useOpenBook->setWhatsThis( tr("<qt>Perform a book opening "
        "move on this triangulation.<p>"
        "A <i>book opening move</i> involves taking an internal face that "
        "meets the boundary of the triangulation along at least one edge "
        "and ungluing the tetrahedra along that face, thereby &quot;opening "
        "out&quot; that face and exposing two more tetrahedron faces to "
        "the boundary.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(useOpenBook, 6, 0);
    useCloseBook = new QRadioButton(tr("C&lose book"), this);
    useCloseBook->setWhatsThis( tr("<qt>Perform a book closing "
        "move on this triangulation.<p>"
        "A <i>book closing move</i> involves taking an edge on the boundary "
        "of the triangulation and folding together the two boundary faces "
        "on either side.  The aim of this move is to simplify the boundary "
        "of the triangulation.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(useCloseBook, 7, 0);
    useShellBdry = new QRadioButton(tr("&Shell boundary"), this);
    useShellBdry->setWhatsThis( tr("<qt>Perform a boundary shelling "
        "move on this triangulation.<p>"
        "A <i>boundary shelling move</i> simply involves removing a "
        "tetrahedron that meets the triangulation boundary along one or "
        "more faces.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(useShellBdry, 8, 0);
    useCollapseEdge = new QRadioButton(tr("&Collapse edge"), this);
    useCollapseEdge->setWhatsThis( tr("<qt>Collapse an edge in this "
        "triangulation.<p>"
        "<i>Collapsing an edge</i> involves taking an edge between two "
        "distinct vertices and collapsing that edge to a point.  Any "
        "tetrahedra containing that edge will be flattened into faces.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered in the adjacent drop-down list.</qt>"));
    layout->addWidget(useCollapseEdge, 9, 0);

    box32 = new EdgeChooser(tri, &has32, this);
    box32->setWhatsThis( tr("<qt>Select the degree three edge about which "
        "the 3-2 move will be performed.  The edge numbers in this list "
        "correspond to the edge numbers seen when viewing the "
        "triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(box32, 0, 1);
    box23 = new FaceChooser(tri, &has23, this);
    box23->setWhatsThis( tr("<qt>Select the face about which "
        "the 2-3 move will be performed.  The face numbers in this list "
        "correspond to the face numbers seen when viewing the "
        "triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(box23, 1, 1);
    box44 = new EdgeIntChooser(tri, 0, 1, tr("axis"), &has44, this);
    box44->setWhatsThis( tr("<qt>Select the degree four edge about which "
        "the 4-4 move will be performed.  You must also select the axis "
        "along which the four new tetrahedra will be inserted (there are "
        "two different ways in which this can be done).<p>"
        "The edge numbers in this list correspond to the edge numbers seen "
        "when viewing the triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(box44, 2, 1);
    box20e = new EdgeChooser(tri, &has20e, this);
    box20e->setWhatsThis( tr("<qt>Select the degree two edge about which "
        "the 2-0 edge move will be performed.  The edge numbers in this list "
        "correspond to the edge numbers seen when viewing the "
        "triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(box20e, 3, 1);
    box20v = new VertexChooser(tri, &has20v, this);
    box20v->setWhatsThis( tr("<qt>Select the degree two vertex about "
        "which the 2-0 vertex move will be performed.  The vertex numbers "
        "in this list correspond to the vertex numbers seen when viewing the "
        "triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(box20v, 4, 1);
    box21 = new EdgeIntChooser(tri, 0, 1, tr("end"), &has21, this);
    box21->setWhatsThis( tr("<qt>Select the degree one edge about which "
        "the 2-1 move will be performed.  You must also select at which "
        "end of the edge the surrounding tetrahedron will be merged with "
        "its neighbour.<p>"
        "The edge numbers in this list correspond to the edge numbers seen "
        "when viewing the triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(box21, 5, 1);
    boxOpenBook = new FaceChooser(tri, &hasOpenBook, this);
    boxOpenBook->setWhatsThis( tr("<qt>Select the internal face "
        "that should be opened out.  The face numbers in this list "
        "correspond to the face numbers seen when viewing the "
        "triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(boxOpenBook, 6, 1);
    boxCloseBook = new EdgeChooser(tri, &hasCloseBook, this);
    boxCloseBook->setWhatsThis( tr("<qt>Select the boundary edge "
        "around which the book will be closed.  The edge numbers in this list "
        "correspond to the edge numbers seen when viewing the "
        "triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(boxCloseBook, 7, 1);
    boxShellBdry = new TetrahedronChooser(tri, &hasShellBoundary, this);
    boxShellBdry->setWhatsThis( tr("<qt>Select the boundary tetrahedron "
        "that should be removed.  The tetrahedron numbers in this list "
        "are the usual tetrahedron numbers seen in the gluings editor.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(boxShellBdry, 8, 1);
    boxCollapseEdge = new EdgeChooser(tri, &hasCollapseEdge, this);
    boxCollapseEdge->setWhatsThis( tr("<qt>Select the edge joining "
        "two distinct vertices that should be collapsed.  "
        "The edge numbers in this list correspond to the edge numbers seen "
        "when viewing the triangulation skeleton.<p>"
        "Only moves that do not change the underlying 3-manifold are "
        "offered.</qt>"));
    layout->addWidget(boxCollapseEdge, 9, 1);

    use32->setEnabled(box32->count() > 0);
    use23->setEnabled(box23->count() > 0);
    use44->setEnabled(box44->count() > 0);
    use20e->setEnabled(box20e->count() > 0);
    use20v->setEnabled(box20v->count() > 0);
    use21->setEnabled(box21->count() > 0);
    useOpenBook->setEnabled(boxOpenBook->count() > 0);
    useCloseBook->setEnabled(boxCloseBook->count() > 0);
    useShellBdry->setEnabled(boxShellBdry->count() > 0);
    useCollapseEdge->setEnabled(boxCollapseEdge->count() > 0);

    box32->setEnabled(box32->count() > 0);
    box23->setEnabled(box23->count() > 0);
    box44->setEnabled(box44->count() > 0);
    box20e->setEnabled(box20e->count() > 0);
    box20v->setEnabled(box20v->count() > 0);
    box21->setEnabled(box21->count() > 0);
    boxOpenBook->setEnabled(boxOpenBook->count() > 0);
    boxCloseBook->setEnabled(boxCloseBook->count() > 0);
    boxShellBdry->setEnabled(boxShellBdry->count() > 0);
    boxCollapseEdge->setEnabled(boxCollapseEdge->count() > 0);

    moveTypes = new QButtonGroup();
    moveTypes->addButton(use32, ID_32);
    moveTypes->addButton(use23, ID_23);
    moveTypes->addButton(use44, ID_44);
    moveTypes->addButton(use20e, ID_20E);
    moveTypes->addButton(use20v, ID_20V);
    moveTypes->addButton(use21, ID_21);
    moveTypes->addButton(useOpenBook, ID_OPENBOOK);
    moveTypes->addButton(useCloseBook, ID_CLOSEBOOK);
    moveTypes->addButton(useShellBdry, ID_SHELLBDRY);
    moveTypes->addButton(useCollapseEdge, ID_COLLAPSEEDGE);

    QDialogButtonBox *buttonBox = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    dialogLayout->addWidget(buttonBox);

    connect(buttonBox, SIGNAL(accepted()), this, SLOT(slotOk()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
}

EltMoveDialog::~EltMoveDialog() {
    delete moveTypes;
}

void EltMoveDialog::slotOk() {
    if (use32->isChecked())
        tri->threeTwoMove(box32->selected());
    else if (use23->isChecked())
        tri->twoThreeMove(box23->selected());
    else if (use44->isChecked()) {
        std::pair<regina::NEdge*, int> s = box44->selected();
        tri->fourFourMove(s.first, s.second);
    } else if (use20e->isChecked())
        tri->twoZeroMove(box20e->selected());
    else if (use20v->isChecked())
        tri->twoZeroMove(box20v->selected());
    else if (use21->isChecked()) {
        std::pair<regina::NEdge*, int> s = box21->selected();
        tri->twoOneMove(s.first, s.second);
    } else if (useOpenBook->isChecked())
        tri->openBook(boxOpenBook->selected());
    else if (useCloseBook->isChecked())
        tri->closeBook(boxCloseBook->selected());
    else if (useShellBdry->isChecked())
        tri->shellBoundary(boxShellBdry->selected());
    else if (useCollapseEdge->isChecked())
        tri->collapseEdge(boxCollapseEdge->selected());
    else {
        ReginaSupport::info(this, tr("Please select a move."));
        return;
    }

    accept();
}

