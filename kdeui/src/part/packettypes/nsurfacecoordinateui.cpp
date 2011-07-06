
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2009, Ben Burton                                   *
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

// Regina core includes:
#include "surfaces/nnormalsurfacelist.h"
#include "surfaces/nsurfacefilter.h"
#include "triangulation/ntriangulation.h"

// UI includes:
#include "coordinatechooser.h"
#include "coordinates.h"
#include "nsurfacecoordinateitem.h"
#include "nsurfacecoordinateui.h"
#include "../packetchooser.h"
#include "../packetfilter.h"
#include "../reginapart.h"

#include <KActionCollection>
#include <KComponentData>
#include <klocale.h>
#include <kmessagebox.h>
#include <QHeaderView>
#include <QHelpEvent>
#include <qlabel.h>
#include <qlayout.h>
#include <qstyle.h>
#include <QTableWidget>
#include <QTableWidgetItem>

#define DEFAULT_COORDINATE_COLUMN_WIDTH 40

using regina::NNormalSurfaceList;
using regina::NPacket;

NSurfaceCoordinateUI::NSurfaceCoordinateUI(regina::NNormalSurfaceList* packet,
        PacketTabbedUI* useParentUI, bool readWrite) :
        PacketEditorTab(useParentUI), surfaces(packet), appliedFilter(0),
        newName(0), isReadWrite(readWrite), currentlyResizing(false) {
    // Prepare the array of modified surface names.

    if (surfaces->getNumberOfSurfaces() > 0)
        newName = new QString[surfaces->getNumberOfSurfaces()];

    // Set up the UI.

    ui = new QWidget();
    uiLayout = new QVBoxLayout(ui);
    uiLayout->addSpacing(5);

    QBoxLayout* hdrLayout = new QHBoxLayout();
    uiLayout->addLayout(hdrLayout);
    hdrLayout->setSpacing(5);
    hdrLayout->addSpacing(5);

    // Set up the coordinate selector.
    QLabel* label = new QLabel(i18n("Display coordinates:"), ui);
    hdrLayout->addWidget(label);
    coords = new CoordinateChooser(ui);
    coords->insertAllViewers(surfaces);
    coords->setCurrentSystem(surfaces->getFlavour());
    connect(coords, SIGNAL(activated(int)), this, SLOT(refreshLocal()));
    hdrLayout->addWidget(coords);
    QString msg = i18n("Allows you to view these normal surfaces in a "
        "different coordinate system.");
    label->setWhatsThis(msg);
    coords->setWhatsThis(msg);

    hdrLayout->addStretch(1);

    // Set up the filter selector.
    label = new QLabel(i18n("Apply filter:"), ui);
    hdrLayout->addWidget(label);
    filter = new PacketChooser(surfaces->getTreeMatriarch(),
        new SingleTypeFilter<regina::NSurfaceFilter>(), true, 0, ui);
    filter->setAutoUpdate(true);
    connect(filter, SIGNAL(activated(int)), this, SLOT(refreshLocal()));
    hdrLayout->addWidget(filter);
    msg = i18n("<qt>Allows you to filter this list so that only normal "
        "surfaces satisfying particular properties are displayed.<p>"
        "To use this feature you need a separate surface filter.  You "
        "can create new surface filters through the <i>Packet Tree</i> "
        "menu.</qt>");
    label->setWhatsThis(msg);
    filter->setWhatsThis(msg);

    hdrLayout->addSpacing(5);
    uiLayout->addSpacing(5);

    // And leave space for the table.
    // We won't actually set up the table until we refresh.
    tableWhatsThis = i18n("Displays details of the individual normal "
        "surfaces in this list.<p>"
        "Each row represents a single normal (or almost normal) surface.  "
        "As well as various properties of the surface, each row contains "
        "a detailed representation the surface in the currently selected "
        "coordinate system.<p>"
        "For details on what each property means or what each coordinate "
        "represents, hover the mouse over the column header (or refer "
        "to the users' handbook).</qt>");

    // Set up the surface list actions.
    surfaceActions = new KActionCollection(0, ReginaPart::factoryInstance());

    actCutAlong = surfaceActions->addAction("surface_cutalong");
    actCutAlong->setText(i18n("Cu&t Along Surface"));
    actCutAlong->setToolTip(i18n("Cut the triangulation along the "
        "selected surface"));
    actCutAlong->setEnabled(false);
    actCutAlong->setWhatsThis(i18n("<qt>Cuts open the surround triangulation "
        "along the selected surface.  This triangulation will not "
        "be changed; instead a new cut-open triangulation will be created.<p>"
        "This operation will never change the topology of the underlying "
        "3-manifold beyond just cutting along the surface (as opposed to "
        "the related <i>crushing</i> operation, which might).  However, "
        "because the new surface boundaries are created from real "
        "boundary faces, the resulting number of tetrahedra might be very "
        "large.</qt>"));
    connect(actCutAlong, SIGNAL(triggered()), this, SLOT(cutAlong()));
    surfaceActionList.append(actCutAlong);

    actCrush = surfaceActions->addAction("surface_crush");
    actCrush->setText("Crus&h Surface");
    actCrush->setToolTip(i18n("Crush the selected surface to a point"));
    actCrush->setEnabled(false);
    actCrush->setWhatsThis(i18n("<qt>Crushes the selected surface to a point "
        "within the surrounding triangulation.  This triangulation will not "
        "be changed; instead a new crushed triangulation will be created.<p>"
        "<b>Warning:</b> This routine simply removes all tetrahedra "
        "containing quadrilateral discs and rejoins the others "
        "appropriately.  In some circumstances this might change the "
        "topology of the underlying 3-manifold beyond just slicing along "
        "the surface and shrinking the resulting boundary/boundaries "
        "to points.</qt>"));
    surfaceActionList.append(actCrush);
    connect(actCrush, SIGNAL(triggered()), this, SLOT(crush()));

    // Tidy up.
    refresh();
}

NSurfaceCoordinateUI::~NSurfaceCoordinateUI() {
    if (newName)
        delete[] newName;

    // Make sure the actions, including separators, are all deleted.
    surfaceActionList.clear();
    delete surfaceActions;
}

const QLinkedList<KAction*>& NSurfaceCoordinateUI::getPacketTypeActions() {
    return surfaceActionList;
}

regina::NPacket* NSurfaceCoordinateUI::getPacket() {
    return surfaces;
}

QWidget* NSurfaceCoordinateUI::getInterface() {
    return ui;
}

void NSurfaceCoordinateUI::commit() {
    for (unsigned long i = 0; i < surfaces->getNumberOfSurfaces(); i++)
        const_cast<regina::NNormalSurface*>(surfaces->getSurface(i))->
            setName(newName[i].toAscii().constData());

    setDirty(false);
}

void NSurfaceCoordinateUI::refreshLocal() {
    // Update the current filter.
    filter->refreshContents();

    if (filter->selectedPacket() != appliedFilter) {
        if (appliedFilter)
            appliedFilter->unlisten(this);
        appliedFilter = dynamic_cast<regina::NSurfaceFilter*>(
            filter->selectedPacket());
        if (appliedFilter)
            appliedFilter->listen(this);
    }

    // Remove the old table.
    table.reset(0);

    // Set up the new table.
    table.reset(new QTableWidget(ui));
    //table->setAllColumnsShowFocus(true);
    table->setSortingEnabled(false);
    table->setSelectionMode(QAbstractItemView::SingleSelection);
    //table->setDefaultRenameAction(QListView::Accept); TODO
    table.get()->setWhatsThis(tableWhatsThis);
    uiLayout->addWidget(table.get(), 1);

    // Add table columns.
    int coordSystem = coords->getCurrentSystem();
    regina::NTriangulation* tri = surfaces->getTriangulation();

    bool embeddedOnly = surfaces->isEmbeddedOnly();
    bool almostNormal = surfaces->allowsAlmostNormal();
    int propCols = NSurfaceCoordinateItem::propertyColCount(embeddedOnly,
        almostNormal);
    long coordCols = Coordinates::numColumns(coordSystem, tri);

    long i;
    QStringList columnHeaders;
    for (i = 0; i < propCols; i++)
        columnHeaders << NSurfaceCoordinateItem::propertyColName(i,
            embeddedOnly, almostNormal);
    for (i = 0; i < coordCols; i++)
        columnHeaders << Coordinates::columnName(coordSystem, i, tri);
   
    table->setHorizontalHeaderLabels(columnHeaders);

    headerTips.reset(new SurfaceHeaderToolTip(surfaces, coordSystem,
        table->horizontalHeader()));
    table->viewport()->installEventFilter(headerTips.get());
    //connect(table->header(), SIGNAL(sizeChange(int, int, int)),
    //    this, SLOT(columnResized(int, int, int)));

    // Insert surfaces into the table.
    const regina::NNormalSurface* s;
    for (i = surfaces->getNumberOfSurfaces() - 1; i >= 0; i--) {
        s = surfaces->getSurface(i);
        if (appliedFilter && ! appliedFilter->accept(*s))
            continue;
        (new NSurfaceCoordinateItem(table.get(), surfaces, i, newName[i],
            coordSystem)); //->setRenameEnabled(1, isReadWrite);
    }

    //for (i = 0; i < table->columns(); i++)
    //    table->adjustColumn(i);

    // Hook up the cut and crush actions to the new table.
    actCutAlong->setEnabled(false);
    actCrush->setEnabled(false);
    connect(table.get(), SIGNAL(selectionChanged()),
        this, SLOT(updateActionStates()));

    // Final tidying up.
    connect(table.get(), SIGNAL(itemRenamed(QListViewItem*, int,
        const QString&)), this, SLOT(notifySurfaceRenamed()));
    table->show();
}

void NSurfaceCoordinateUI::refresh() {
    // Refresh the surface names from the underlying packet.
    for (unsigned long i = 0; i < surfaces->getNumberOfSurfaces(); i++)
        newName[i] = surfaces->getSurface(i)->getName().c_str();

    // Refresh the table of surfaces.
    refreshLocal();

    setDirty(false);
}

void NSurfaceCoordinateUI::setReadWrite(bool readWrite) {
    isReadWrite = readWrite;

    if (table.get()) {
        QList<QTableWidgetItem*> children = table->findItems("",
            Qt::MatchWrap|Qt::MatchRecursive);
        for ( int i=0; i < children.count() ; i++) {
            QTableWidgetItem* item = children[i];
            Qt::ItemFlags flags = item->flags();
            if (readWrite)
                flags = flags | Qt::ItemIsEditable;
            else 
                flags = flags & (~Qt::ItemIsEditable);
            item->setFlags(flags);
        }
    }
    updateActionStates();
}

void NSurfaceCoordinateUI::packetToBeDestroyed(NPacket*) {
    // Our currently applied filter is about to be destroyed.
    filter->setCurrentItem(0); // (i.e., None)
    refreshLocal();
}

void NSurfaceCoordinateUI::cutAlong() {
    QList<QTableWidgetItem*> items = table->selectedItems();
    if ( items.count() == 0) {
        KMessageBox::error(ui,
            i18n("No normal surface is currently selected to cut along."));
        return;
    }
    QTableWidgetItem* item = items[0];

    const regina::NNormalSurface* toCutAlong =
        dynamic_cast<NSurfaceCoordinateItem*>(item)->getSurface();
    if (! toCutAlong->isCompact()) {
        KMessageBox::error(ui, i18n("The selected surface is non-compact "
            "and so cannot be cut along."));
        return;
    }

    // Go ahead and cut along the surface.
    // Be nice and simplify the triangulation, which could be very large.
    regina::NTriangulation* ans = toCutAlong->cutAlong();
    ans->intelligentSimplify();
    ans->setPacketLabel(surfaces->makeUniqueLabel(i18n("Cut-open %1").arg(
        surfaces->getTriangulation()->getPacketLabel().c_str())
          .toAscii().constData()));
    surfaces->insertChildLast(ans);

    enclosingPane->getPart()->packetView(ans, true);
}

void NSurfaceCoordinateUI::crush() {
    QList<QTableWidgetItem*> items = table->selectedItems();
    if ( items.count() == 0) {
        KMessageBox::error(ui,
            i18n("No normal surface is currently selected to crush."));
        return;
    }
    QTableWidgetItem* item = items[0];

    const regina::NNormalSurface* toCrush =
        dynamic_cast<NSurfaceCoordinateItem*>(item)->getSurface();
    if (! toCrush->isCompact()) {
        KMessageBox::error(ui, i18n("The selected surface is non-compact "
            "and so cannot be crushed."));
        return;
    }

    // Go ahead and crush it.
    regina::NTriangulation* ans = toCrush->crush();
    ans->setPacketLabel(surfaces->makeUniqueLabel(i18n("Crushed %1").arg(
        surfaces->getTriangulation()->getPacketLabel().c_str())
          .toAscii().constData()));
    surfaces->insertChildLast(ans);

    enclosingPane->getPart()->packetView(ans, true);
}

void NSurfaceCoordinateUI::updateActionStates() {
    bool canCrushOrCut = isReadWrite && table.get() &&
        table->selectedItems().count() != 0 && (! surfaces->allowsAlmostNormal())
        && surfaces->isEmbeddedOnly();

    actCutAlong->setEnabled(canCrushOrCut);
    actCrush->setEnabled(canCrushOrCut);
}

void NSurfaceCoordinateUI::columnResized(int section, int, int newSize) {
    int nNonCoordSections = NSurfaceCoordinateItem::propertyColCount(
        surfaces->isEmbeddedOnly(), surfaces->allowsAlmostNormal());
    if (currentlyResizing || section < nNonCoordSections)
        return;

    // A coordinate column has been resized.
    // Resize all coordinate columns.
    currentlyResizing = true;
    //for (long i = nNonCoordSections; i < table->columnCount(); i++)
    //    table->setColumnWidth(i, newSize);
    currentlyResizing = false;
}

void NSurfaceCoordinateUI::notifySurfaceRenamed() {
    setDirty(true);
}

SurfaceHeaderToolTip::SurfaceHeaderToolTip(
        regina::NNormalSurfaceList* useSurfaces, int useCoordSystem,
        QHeaderView* header) : header(header), surfaces(useSurfaces),
    coordSystem(useCoordSystem) {
}

bool SurfaceHeaderToolTip::eventFilter(QObject *obj, QEvent *event) {
    if ( event->type() == QEvent::ToolTip) {
        QHelpEvent *helpEvent = static_cast<QHelpEvent *>(event);
        QPoint p = helpEvent->pos();
        int section = header->logicalIndexAt(p.x());
        if (section < 0)
            return false;

        int propertyCols = NSurfaceCoordinateItem::propertyColCount(
            surfaces->isEmbeddedOnly(), surfaces->allowsAlmostNormal());

        QString tipString;
        if (section < propertyCols)
            tipString = NSurfaceCoordinateItem::propertyColDesc(section,
                surfaces->isEmbeddedOnly(), surfaces->allowsAlmostNormal());
        else
            tipString = Coordinates::columnDesc(coordSystem,
                section - propertyCols, surfaces->getTriangulation());
        QToolTip::showText(helpEvent->globalPos(), tipString);
        return true;
    }
    return false;
}

