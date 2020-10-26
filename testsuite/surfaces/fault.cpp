
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
#include "manifold/sfs.h"
#include "surfaces/normalcoords.h"
#include "surfaces/normalflags.h"
#include "surfaces/normalsurface.h"
#include "surfaces/normalsurfaces.h"
#include "triangulation/dim3.h"
#include "triangulation/example3.h"
#include "testsuite/dim3/testtriangulation.h"
#include "triangulation/generic/triangulation.h"
#include "manifold/sfs.h"

#include <list>

using regina::Example;
using regina::NormalListFlags;
using regina::NormalSurfaces;
using regina::NS_EMBEDDED_ONLY;
using regina::NS_QUAD;
using regina::NS_STANDARD;
using regina::NS_VERTEX;
using regina::SFSFibre;
using regina::SFSpace;
using regina::Triangulation;


class FaultTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(FaultTest);

    CPPUNIT_TEST(separates);
    CPPUNIT_TEST(isEssentialSphere);
    CPPUNIT_TEST(isEssentialTorus);
    CPPUNIT_TEST(isSolidTorusAnnulus);
    CPPUNIT_TEST(findFault);
    
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
        tri = Example<3>::threeSphere();
        delete verifyAllFundamentalToriInessential(tri, "Minimal 3-sphere");
        tri = Example<3>::simplicialSphere();
        delete verifyAllFundamentalToriInessential(tri, "Pentachoron boundary 3-sphere");
        tri = Example<3>::poincareHomologySphere();
        delete verifyAllFundamentalToriInessential(tri, "Poincare homology sphere");
        tri = Example<3>::weeks();
        delete verifyAllFundamentalToriInessential(tri, "Weeks-Matveev-Fomenko manifold");
        int p = 3;
        int q = 2;
        while (p < 1000){
            tri = Example<3>::lens(p,q);
            delete verifyAllFundamentalToriInessential(tri, "Fibonacci lens space");
            p = p + q;
            q = p - q;
        }

        // Bounded atoroidal manifolds
        tri = Example<3>::ball();
        delete verifyAllFundamentalToriInessential(tri, "One tetrahedron ball");
        tri = Example<3>::cuspedGenusTwoTorus();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(
            tri, "Trivial I-bundle over genus two surface");
        tri = Example<3>::figureEight();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(tri, "Figure eight knot exterior");
        tri = Example<3>::trefoil();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(tri, "Trefoil knot exterior");
        tri = Example<3>::whiteheadLink();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(tri, "Whitehead link exterior");
        tri = Example<3>::gieseking();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(tri, "Gieseking manifold");
        tri = Example<3>::solidKleinBottle();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(tri, "Solid Klein bottle");
        p = 3;
        q = 2;
        while (p < 1000){
            tri = Example<3>::lst(p,q);
            delete verifyAllFundamentalToriInessential(tri, "Solid torus");
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

        Triangulation<3>* verifyAllFundamentalAnnuli(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, regina::NS_EMBEDDED_ONLY);
        bool found_solid_torus_annulus = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isSolidTorusAnnulus()){
                found_solid_torus_annulus = true;
                break;
            }
        delete s;
        if (found_solid_torus_annulus)
            CPPUNIT_FAIL(("The manifold " + triName +
                          " was computed to split into two solid tori along an annulus.").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasSolidTorusAnnulus(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_STANDARD, NS_EMBEDDED_ONLY);
        bool found_solid_torus_annulus = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isSolidTorusAnnulus()){
                found_solid_torus_annulus = true;
                break;
            }
        delete s;
        if (!found_solid_torus_annulus)
            CPPUNIT_FAIL(("The fundamental search did not find an annulus that splits " +
                          triName + " into two solid tori.").c_str());
        return tri;
    }

    Triangulation<3>* verifyHasQVertexSolidTorusAnnulus(Triangulation<3>* tri, 
            const std::string& triName){
        NormalSurfaces* s = NormalSurfaces::enumerate(
            tri, regina::NS_QUAD,
            static_cast<NormalListFlags>(NS_VERTEX & NS_EMBEDDED_ONLY));
        bool found_solid_torus_annulus = false;
        for (unsigned i = 0; i < s->size(); ++i)
            if (s->surface(i)->isSolidTorusAnnulus()){
                found_solid_torus_annulus = true;
                break;
            }
        delete s;
        if (!found_solid_torus_annulus)
            CPPUNIT_FAIL(("The quad-vertex search did not find an annulus that splits " +
                          triName + " into two solid tori.").c_str());
        return tri;
    }

    void isSolidTorusAnnulus() {
        Triangulation<3>* tri;
        
        // Closed anannular manifolds
        tri = Example<3>::threeSphere();
        delete verifyAllFundamentalToriInessential(tri, "Minimal 3-sphere");
        tri = Example<3>::simplicialSphere();
        delete verifyAllFundamentalToriInessential(tri, "Pentachoron boundary 3-sphere");
        tri = Example<3>::poincareHomologySphere();
        delete verifyAllFundamentalToriInessential(tri, "Poincare homology sphere");
        tri = Example<3>::weeks();
        delete verifyAllFundamentalToriInessential(tri, "Weeks-Matveev-Fomenko manifold");
        int p = 3;
        int q = 2;
        while (p < 1000){
            tri = Example<3>::lens(p,q);
            delete verifyAllFundamentalToriInessential(tri, "Fibonacci lens space");
            p = p + q;
            q = p - q;
        }
        
        
        // Bounded anannular manifolds
        tri = Example<3>::ball();
        delete verifyAllFundamentalToriInessential(tri, "One tetrahedron ball");
        tri = Example<3>::cuspedGenusTwoTorus();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(
            tri, "Trivial I-bundle over genus two surface");
        tri = Example<3>::figureEight();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyAllFundamentalToriInessential(tri, "Figure eight knot exterior");

        // Annular manifolds

        // Quad vertex coordinates aren't necessary for these, seemingly.
        tri = Example<3>::trefoil();
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasSolidTorusAnnulus(tri, "Trefoil knot complement");

        tri = new Triangulation<3>("gLKjMcacdeffjkxnqs");
        delete verifyHasSolidTorusAnnulus(tri, "SFS D: (3,1) (5,-2)");
        
        tri = new Triangulation<3>("mHLyMzMzMcbcdefghijklldwuxhqxhqxhw");
        delete verifyHasSolidTorusAnnulus(tri, "SFS D: (2,1) (144,-89)");
    }

    Triangulation<3>* verifyHasFault(Triangulation<3>* tri, 
            const std::string& triName){
        if (!tri->findFault(NS_QUAD, NS_VERTEX)){
            CPPUNIT_FAIL(("No fault found in " + triName + ".").c_str());
        }
        return tri;
    }

    Triangulation<3>* verifyFaultless(Triangulation<3>* tri, 
            const std::string& triName){
        if (tri->findFault(NS_QUAD, NS_VERTEX)){
            CPPUNIT_FAIL(("Supposed \"fault\" found in " + triName + ".").c_str());
        }
        return tri;
    }

    Triangulation<3>* small(long alpha0, long beta0,
                   long alpha1, long beta1,
                   long alpha2, long beta2){
        const SFSFibre* fiber0 = new SFSFibre(alpha0, beta0);
        const SFSFibre* fiber1 = new SFSFibre(alpha1, beta1);
        const SFSFibre* fiber2 = new SFSFibre(alpha2, beta2);
        SFSpace* sfs = new SFSpace(
            SFSpace::o1, 0, 0, 0, 0, 0);
        sfs->insertFibre(*fiber0);
        sfs->insertFibre(*fiber1);
        sfs->insertFibre(*fiber2);
        Triangulation<3>* sfs_tri = sfs->construct();
        delete sfs;
        delete fiber0;
        delete fiber1;
        delete fiber2;
        return sfs_tri;
    }
    
    void findFault(){
        Triangulation<3>* tri;
        
        // Faultless
        // Faultless: closed manifolds

        // Faultless: Small SFSes with positive orbifold Euler char.
        tri = Example<3>::threeSphere();
        delete verifyFaultless(tri, "Minimal 3-sphere");
        tri = Example<3>::simplicialSphere();
        delete verifyFaultless(tri, "Pentachoron boundary 3-sphere");

        int p = 3;
        int q = 2;
        while (p < 1000){
            tri = Example<3>::lens(p,q);
            delete verifyFaultless(tri, "Fibonacci lens space");
            p = p + q;
            q = p - q;
        }
        tri = Example<3>::poincareHomologySphere();
        delete verifyFaultless(tri, "Poincare homology sphere");

#define VERIFY_SMALL_FAULTLESS(a0,b0,a1,b1,a2,b2) \
        tri = small(a0,b0,a1,b1,a2,b2);           \
        tri->intelligentSimplify();               \
        delete verifyFaultless(tri, "SFS [S2: (a0,b0) (a1,b1) (a2,b2)]")

        VERIFY_SMALL_FAULTLESS(2,1,3,1,3,-2);
        
        VERIFY_SMALL_FAULTLESS(2,1,3,1,3,-1);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,4,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,5,-4);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,3,-1);
        
        VERIFY_SMALL_FAULTLESS(2,1,3,1,3,1);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,3,2);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,4,-1);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,5,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,5,-2);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,3,2);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,4,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,4,-1);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,5,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,5,-2);

        // Faultless: Small SFSes with negative orbifold Euler char.
        VERIFY_SMALL_FAULTLESS(2,1,3,1,7,-6);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,7,-5);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,7,-5);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,7,-4);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,7,-2);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,8,-5);
        VERIFY_SMALL_FAULTLESS(2,1,3,1,8,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,7,-5);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,7,-4);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,7,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,7,-2);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,8,-5);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,8,-3);
        VERIFY_SMALL_FAULTLESS(2,1,3,2,8,-3);
        VERIFY_SMALL_FAULTLESS(2,1,4,1,5,-4);
        VERIFY_SMALL_FAULTLESS(2,1,4,1,5,-3);
        VERIFY_SMALL_FAULTLESS(2,1,4,1,5,-2);
        VERIFY_SMALL_FAULTLESS(2,1,4,3,5,-3);
        VERIFY_SMALL_FAULTLESS(2,1,4,3,5,-2);
        VERIFY_SMALL_FAULTLESS(2,1,5,2,5,-3);
        VERIFY_SMALL_FAULTLESS(2,1,5,2,5,-2);
        VERIFY_SMALL_FAULTLESS(2,1,5,3,5,-2);

        VERIFY_SMALL_FAULTLESS(3,1,3,1,4,-3);
        VERIFY_SMALL_FAULTLESS(3,1,3,1,4,-1);
        VERIFY_SMALL_FAULTLESS(3,1,3,1,5,-3);
        VERIFY_SMALL_FAULTLESS(3,1,3,1,5,-2);
        VERIFY_SMALL_FAULTLESS(3,1,3,2,4,-3);
        VERIFY_SMALL_FAULTLESS(3,1,3,2,4,-1);
        VERIFY_SMALL_FAULTLESS(3,1,3,2,5,-3);
        VERIFY_SMALL_FAULTLESS(3,1,3,2,5,-2);
        VERIFY_SMALL_FAULTLESS(3,2,3,2,4,-3);
        VERIFY_SMALL_FAULTLESS(3,2,3,2,4,-1);
        VERIFY_SMALL_FAULTLESS(3,2,3,2,5,-3);
        VERIFY_SMALL_FAULTLESS(3,2,3,2,5,-2);

        // Faultless: Closed hyperbolic
        tri = Example<3>::weeks();
        delete verifyFaultless(tri, "Weeks-Matveev-Fomenko manifold");
        tri = new Triangulation<3>("jLvALQQadfgfihhiijaqjvuxxec");
        delete verifyFaultless(tri, "Vol2");
        tri = new Triangulation<3>("jLLLzQQcdfeghihiiqceqgvlkuq");
        delete verifyFaultless(tri, "Vol3");
        tri = new Triangulation<3>("jLLvQMQcdfeghiihihsqaggqtni");
        delete verifyFaultless(tri, "Vol4");

        // Faultless: bounded, < 10 tets 
        tri = new Triangulation<3>("ifvLQQiegghfhfhlxhixgon");
        delete verifyFaultless(tri, "truncated m003");
        tri = new Triangulation<3>("jLLvQMIacdfeghiiihsqaggqti");
        delete verifyFaultless(tri, "truncated m006");
        tri = new Triangulation<3>("jLLLzQkacdfeghiihhsmdfvtaw");
        delete verifyFaultless(tri, "truncated m010");
        
        // Faulty

        // Faulty: Example<3> triangulations
        tri = Example<3>::rp3rp3();
        delete verifyHasFault(tri, "RP3#RP3");
        tri = Example<3>::rp2xs1();
        delete verifyHasFault(tri, "RP2xS1");
        tri = Example<3>::s2xs1();
        delete verifyHasFault(tri, "S2xS1");
        tri = Example<3>::solidKleinBottle();
        delete verifyHasFault(tri, "Solid Klein bottle");
        tri = Example<3>::ballBundle();
        delete verifyHasFault(tri, "Solid torus");
        tri = Example<3>::trefoil();
        delete verifyHasFault(tri, "Trefoil knot exterior");
        tri = Example<3>::cuspedGenusTwoTorus();
        delete verifyHasFault(tri, "Genus two handlebody");

        // Faulty: Others
        tri = Example<3>::trefoil();
        tri->connectedSumWith(*tri);
        tri->idealToFinite();
        tri->intelligentSimplify();
        delete verifyHasFault(tri, "Trefoil # trefoil");

        tri = new Triangulation<3>("jLLLAAQbdeffghiiixniqxruawo");
        delete verifyHasFault(tri, "Trefoil cable");

        tri = new Triangulation<3>(
            "uLLLzAQwzwLAMQccdefhhjihjlnrqqpssttthsdmsvhdkbbmaqdawwiei");
        delete verifyHasFault(tri, "Trefoil knot exterior doubled along boundary");
    }
};

void addFaultFinding(CppUnit::TextUi::TestRunner& runner) {
    runner.addTest(FaultTest::suite());
}
