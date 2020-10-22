
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Test Suite                                                            *
 *                                                                        *
 *  Copyright (c) 2020, Robert C. Haraway, III.                           *
 *  For further details contact Robert Haraway (bobbycyiii@fastmail.com). *
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

#include <cppunit/extensions/HelperMacros.h>
#include "surfaces/normalcoords.h"
#include "surfaces/normalflags.h"
#include "surfaces/normalsurface.h"
#include "surfaces/normalsurfaces.h"
#include "triangulation/dim3.h"
#include "triangulation/example3.h"
#include "testsuite/dim3/testtriangulation.h"
#include "triangulation/generic/triangulation.h"


using regina::Example;
using regina::NormalSurfaces;
using regina::Triangulation;
using regina::NormalListFlags;
using regina::NS_VERTEX;
using regina::NS_EMBEDDED_ONLY;

class FaultTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(FaultTest);

    CPPUNIT_TEST(separates);
    CPPUNIT_TEST(isEssentialSphere);
    CPPUNIT_TEST(isEssentialTorus);
    CPPUNIT_TEST(isSolidTorusAnnulus);
    
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp() {
    }

    void tearDown() {
    }

    Triangulation<3>* verifyAllSeparating(Triangulation<3>* tri, 
            const std::string& triName) {
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, regina::NS_EMBEDDED_ONLY);
        bool nonseparating_known = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (!(s->surface(i)->separates())){
                nonseparating_known = true;
                break;
            }
        delete s;
        if (nonseparating_known)
            CPPUNIT_FAIL(("A surface in " + triName +
                          " is computed to be nonseparating.").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasNonSeparating(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, regina::NS_EMBEDDED_ONLY);
        bool all_known_surfaces_separating = true;
        for (unsigned i = 0; i < s->size(); ++i)
            if (!(s->surface(i)->separates())){
                all_known_surfaces_separating = false;
                break;
            }
        delete s;
        if (all_known_surfaces_separating)
            CPPUNIT_FAIL(("No surfaces in " + triName +
                          " were computed to be nonseparating.").c_str());
        return tri;
    }

    void separates() {
        Triangulation<3>* tri;

        // Manifolds without nonseparating surfaces

        tri = Example<3>::threeSphere();
        delete verifyAllSeparating(tri, "Minimal 3-sphere");
        tri = Example<3>::simplicialSphere();
        delete verifyAllSeparating(tri, "Pentachoron boundary 3-sphere");
        tri = Example<3>::ball();
        delete verifyAllSeparating(tri, "One-tetrahedron ball");

        int p = 3;
        int q = 2;
        while (p < 1000){
            if (p % 2 != 0){
                tri = Example<3>::lens(p,q);
                delete verifyAllSeparating(tri, "Fibonacci lens space with odd p");
            }
            p = p + q;
            q = p - q;
        }
        
        tri = Example<3>::poincareHomologySphere();
        delete verifyAllSeparating(tri, "Poincare homology sphere");
        tri = Example<3>::weeks();
        delete verifyAllSeparating(tri, "Weeks-Matveev-Fomenko manifold");


        // Manifolds with nonseparating surfaces
        
        tri = Example<3>::s2xs1();
        delete verifyHasNonSeparating(tri, "S2xS1");
        tri = Example<3>::rp2xs1();
        delete verifyHasNonSeparating(tri, "RP2xS1");
        tri = Example<3>::rp3rp3();
        delete verifyHasNonSeparating(tri, "RP3#RP3");
        tri = Example<3>::smallClosedNonOrblHyperbolic();
        delete verifyHasNonSeparating(tri, "Smallest known closed nonorientable hyperbolic");

        p = 3;
        q = 2;
        while (p < 1000){
            tri = Example<3>::lst(p,q);
            delete verifyHasNonSeparating(tri, "Solid torus");
            if (p % 2 == 0){
                tri = Example<3>::lens(p,q);
                delete verifyHasNonSeparating(tri, "Fibonacci lens space with even p");
            }
            p = p + q;
            q = p - q;
        }

        tri = Example<3>::solidKleinBottle();
        delete verifyHasNonSeparating(tri, "Solid Klein bottle");

        tri = Example<3>::figureEight();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasNonSeparating(tri, "Figure eight knot exterior");
        tri = Example<3>::trefoil();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasNonSeparating(tri, "Trefoil knot exterior");
        tri = Example<3>::whiteheadLink();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasNonSeparating(tri, "Whitehead link");
        tri = Example<3>::gieseking();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasNonSeparating(tri, "Gieseking manifold");
        tri = Example<3>::cuspedGenusTwoTorus();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasNonSeparating(tri, "Genus two surface x I");
    }

    Triangulation<3>* verifyAllFundamentalSpheresTrivial(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, regina::NS_EMBEDDED_ONLY);
        bool found_essential_sphere = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isEssentialSphere()){
                found_essential_sphere = true;
                break;
            }
        delete s;
        if (found_essential_sphere)
            CPPUNIT_FAIL(("The irreducible manifold " + triName +
                          " was computed to have an essential sphere.").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasEssentialSphere(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, NS_EMBEDDED_ONLY);
        bool found_essential_sphere = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isEssentialSphere()){
                found_essential_sphere = true;
                break;
            }
        delete s;
        if (!found_essential_sphere)
            CPPUNIT_FAIL(("No essential spheres were found in the composite manifold " +
                          triName + ".").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasQVertexEssentialSphere(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_QUAD,
            static_cast<NormalListFlags>(NS_VERTEX & NS_EMBEDDED_ONLY));
        bool found_essential_sphere = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isEssentialSphere()){
                found_essential_sphere = true;
                break;
            }
        delete s;
        if (!found_essential_sphere)
            CPPUNIT_FAIL(("No essential spheres were found in the composite manifold " +
                          triName + ".").c_str());
        return tri;
    }

    void isEssentialSphere() {
        Triangulation<3>* tri;

        // Closed irreducible manifolds
        tri = Example<3>::threeSphere();
        delete verifyAllFundamentalSpheresTrivial(tri, "Minimal 3-sphere");
        tri = Example<3>::simplicialSphere();
        delete verifyAllFundamentalSpheresTrivial(tri, "Pentachoron boundary 3-sphere");
        tri = Example<3>::poincareHomologySphere();
        delete verifyAllFundamentalSpheresTrivial(tri, "Poincare homology sphere");
        tri = Example<3>::weeks();
        delete verifyAllFundamentalSpheresTrivial(tri, "Weeks-Matveev-Fomenko manifold");
        int p = 3;
        int q = 2;
        while (p < 1000){
            tri = Example<3>::lens(p,q);
            delete verifyAllFundamentalSpheresTrivial(tri, "Fibonacci lens space");
            p = p + q;
            q = p - q;
        }

        // Bounded irreducible manifolds
        tri = Example<3>::ball();
        delete verifyAllFundamentalSpheresTrivial(tri, "One tetrahedron ball");
        tri = Example<3>::cuspedGenusTwoTorus();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(
            tri, "Trivial I-bundle over genus two surface");
        tri = Example<3>::figureEight();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Figure eight knot exterior");
        tri = Example<3>::trefoil();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Trefoil knot exterior");
        tri = Example<3>::whiteheadLink();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Whitehead link exterior");
        tri = Example<3>::gieseking();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Gieseking manifold");
        tri = Example<3>::solidKleinBottle();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Solid Klein bottle");
        p = 3;
        q = 2;
        while (p < 1000){
            tri = Example<3>::lst(p,q);
            delete verifyAllFundamentalSpheresTrivial(tri, "Solid torus");
            p = p + q;
            q = p - q;
        }

        // The prime reducible manifold S^2 x S^1
        tri = Example<3>::s2xs1();
        delete verifyHasEssentialSphere(tri, "S2xS1");
        
        // RP3#RP3 has a reducing sphere,
        // but in the given triangulation, there is no fundamental reducing sphere.
        // Instead, there are two fundamental embedded one-sided RP2s
        // that double to reducing spheres.

        tri = Example<3>::rp3rp3();
        delete verifyAllFundamentalSpheresTrivial(tri, "RP3#RP3");
        
        // Reducible manifolds require larger triangulations.
        // To reduce testing time, we use quad-vertex surfaces.
        
        p = 3;
        q = 2;
        while (p < 100){
            tri = Example<3>::lens(p,q);
            tri->connectedSumWith(*tri);
            tri->intelligentSimplify();
            delete verifyHasQVertexEssentialSphere(tri, "Lens space # itself");
            p = p + q;
            q = p - q;
        }
    }


    Triangulation<3>* verifyAllFundamentalToriInessential(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, regina::NS_EMBEDDED_ONLY);
        bool found_essential_torus = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isEssentialTorus()){
                found_essential_torus = true;
                break;
            }
        delete s;
        if (found_essential_torus)
            CPPUNIT_FAIL(("The atoroidal manifold " + triName +
                          " was computed to have an essential torus.").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasEssentialTorus(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, NS_EMBEDDED_ONLY);
        bool found_essential_torus = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isEssentialTorus()){
                found_essential_torus = true;
                break;
            }
        delete s;
        if (!found_essential_torus)
            CPPUNIT_FAIL(("No essential tori were found in the toroidal manifold " +
                          triName + ".").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasQVertexEssentialTorus(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_QUAD,
            static_cast<NormalListFlags>(NS_VERTEX & NS_EMBEDDED_ONLY));
        bool found_essential_torus = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isEssentialTorus()){
                found_essential_torus = true;
                break;
            }
        delete s;
        if (!found_essential_torus)
            CPPUNIT_FAIL(("No essential tori were found in the toroidal manifold " +
                          triName + ".").c_str());
        return tri;
    }
    
    void isEssentialTorus() {
        Triangulation<3>* tri;

        // Closed atoroidal manifolds
        // Closed irreducible manifolds
        tri = Example<3>::threeSphere();
        delete verifyAllFundamentalSpheresTrivial(tri, "Minimal 3-sphere");
        tri = Example<3>::simplicialSphere();
        delete verifyAllFundamentalSpheresTrivial(tri, "Pentachoron boundary 3-sphere");
        tri = Example<3>::poincareHomologySphere();
        delete verifyAllFundamentalSpheresTrivial(tri, "Poincare homology sphere");
        tri = Example<3>::weeks();
        delete verifyAllFundamentalSpheresTrivial(tri, "Weeks-Matveev-Fomenko manifold");
        int p = 3;
        int q = 2;
        while (p < 1000){
            tri = Example<3>::lens(p,q);
            delete verifyAllFundamentalSpheresTrivial(tri, "Fibonacci lens space");
            p = p + q;
            q = p - q;
        }

        // Bounded atoroidal manifolds
        tri = Example<3>::ball();
        delete verifyAllFundamentalSpheresTrivial(tri, "One tetrahedron ball");
        tri = Example<3>::cuspedGenusTwoTorus();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(
            tri, "Trivial I-bundle over genus two surface");
        tri = Example<3>::figureEight();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Figure eight knot exterior");
        tri = Example<3>::trefoil();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Trefoil knot exterior");
        tri = Example<3>::whiteheadLink();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Whitehead link exterior");
        tri = Example<3>::gieseking();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Gieseking manifold");
        tri = Example<3>::solidKleinBottle();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalSpheresTrivial(tri, "Solid Klein bottle");
        p = 3;
        q = 2;
        while (p < 1000){
            tri = Example<3>::lst(p,q);
            delete verifyAllFundamentalSpheresTrivial(tri, "Solid torus");
            p = p + q;
            q = p - q;
        }

        // Toroidal manifolds require larger triangulations.
        // To reduce testing time, we use quad-vertex surfaces.
        tri = new Triangulation<3>("uLLvPAPAzzvQPQccdeghiihjjlmqspstrstthshgbvrndhakecbcqvndm");
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasQVertexEssentialTorus(tri, "Doubled figure eight knot exterior");

        tri = new Triangulation<3>("iLALzQcbccefhgghlpkkucjjs");
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasQVertexEssentialTorus(tri, "Doubled trefoil knot exterior");
    }

    void isSolidTorusAnnulus() {
        CPPUNIT_FAIL("Not implemented yet");
    }

};

void addFaultFinding(CppUnit::TextUi::TestRunner& runner) {
    runner.addTest(FaultTest::suite());
}
