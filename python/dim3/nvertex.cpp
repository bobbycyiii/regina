
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
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

#include <boost/python.hpp>
#include "triangulation/dim2.h"
#include "triangulation/dim3.h"
#include "../globalarray.h"
#include "../helpers.h"
#include "../safeheldtype.h"

using namespace boost::python;
using namespace regina::python;
using regina::Face;
using regina::FaceEmbedding;
using regina::NVertex;
using regina::VertexEmbedding;

namespace {
    boost::python::list vertex_embeddings_list(const NVertex* v) {
        boost::python::list ans;
        for (auto& emb: *v)
            ans.append(emb);
        return ans;
    }

    regina::Triangulation<2>* vertex_buildLink(const NVertex* v) {
        return new regina::Triangulation<2>(*(v->buildLink()));
    }

    boost::python::tuple vertex_buildLinkDetail_bool(const NVertex* v,
            bool labels = true) {
        regina::NIsomorphism* iso;
        regina::Triangulation<2>* link = new regina::Triangulation<2>(
            *(v->buildLinkDetail(labels, &iso)));
        return boost::python::make_tuple(
            boost::python::object(boost::python::handle<>(
                regina::python::to_held_type<>::
                apply<regina::Triangulation<2>*>::type()(link))),
            boost::python::object(boost::python::handle<>(
                boost::python::manage_new_object::
                apply<regina::NIsomorphism*>::type()(iso))));
    }

    boost::python::tuple vertex_buildLinkDetail_void(const NVertex* v) {
        return vertex_buildLinkDetail_bool(v);
    }
}

void addNVertex() {
    class_<FaceEmbedding<3, 0>>("FaceEmbedding3_0",
            init<regina::Tetrahedron<3>*, int>())
        .def(init<const VertexEmbedding<3>&>())
        .def("simplex", &VertexEmbedding<3>::simplex,
            return_value_policy<reference_existing_object>())
        .def("tetrahedron", &VertexEmbedding<3>::tetrahedron,
            return_value_policy<reference_existing_object>())
        .def("face", &VertexEmbedding<3>::face)
        .def("vertex", &VertexEmbedding<3>::vertex)
        .def("vertices", &VertexEmbedding<3>::vertices)
        .def(regina::python::add_output())
        .def(regina::python::add_eq_operators())
    ;

    {
        scope s = class_<Face<3, 0>, std::auto_ptr<Face<3, 0>>,
                boost::noncopyable>("Face3_0", no_init)
            .def("index", &NVertex::index)
            .def("embeddings", vertex_embeddings_list)
            .def("embedding", &NVertex::embedding,
                return_internal_reference<>())
            .def("front", &NVertex::front,
                return_internal_reference<>())
            .def("back", &NVertex::back,
                return_internal_reference<>())
            .def("triangulation", &NVertex::triangulation,
                return_value_policy<to_held_type<> >())
            .def("component", &NVertex::component,
                return_value_policy<reference_existing_object>())
            .def("boundaryComponent", &NVertex::boundaryComponent,
                return_value_policy<reference_existing_object>())
            .def("degree", &NVertex::degree)
            .def("link", &NVertex::link)
            .def("buildLink", &vertex_buildLink,
                return_value_policy<to_held_type<>>())
            .def("buildLinkDetail", vertex_buildLinkDetail_void)
            .def("buildLinkDetail", vertex_buildLinkDetail_bool)
            .def("isLinkClosed", &NVertex::isLinkClosed)
            .def("isIdeal", &NVertex::isIdeal)
            .def("isBoundary", &NVertex::isBoundary)
            .def("isStandard", &NVertex::isStandard)
            .def("isValid", &NVertex::isValid)
            .def("hasBadIdentification", &NVertex::hasBadIdentification)
            .def("hasBadLink", &NVertex::hasBadLink)
            .def("isLinkOrientable", &NVertex::isLinkOrientable)
            .def("linkEulerChar", &NVertex::linkEulerChar)
            .def("ordering", &NVertex::ordering)
            .def("faceNumber", &NVertex::faceNumber)
            .def("containsVertex", &NVertex::containsVertex)
            .def(regina::python::add_output())
            .def(regina::python::add_eq_operators())
            .staticmethod("ordering")
            .staticmethod("faceNumber")
            .staticmethod("containsVertex")
        ;

        enum_<regina::NVertex::LinkType>("LinkType")
            .value("SPHERE", regina::NVertex::SPHERE)
            .value("DISC", regina::NVertex::DISC)
            .value("TORUS", regina::NVertex::TORUS)
            .value("KLEIN_BOTTLE", regina::NVertex::KLEIN_BOTTLE)
            .value("NON_STANDARD_CUSP", regina::NVertex::NON_STANDARD_CUSP)
            .value("INVALID", regina::NVertex::INVALID)
            ;

        s.attr("SPHERE") = NVertex::SPHERE;
        s.attr("DISC") = NVertex::DISC;
        s.attr("TORUS") = NVertex::TORUS;
        s.attr("KLEIN_BOTTLE") = NVertex::KLEIN_BOTTLE;
        s.attr("NON_STANDARD_CUSP") = NVertex::NON_STANDARD_CUSP;
        s.attr("INVALID") = NVertex::INVALID;
    }

    scope().attr("NVertexEmbedding") = scope().attr("FaceEmbedding3_0");
    scope().attr("VertexEmbedding3") = scope().attr("FaceEmbedding3_0");
    scope().attr("NVertex") = scope().attr("Face3_0");
    scope().attr("Vertex3") = scope().attr("Face3_0");
}
