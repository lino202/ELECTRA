/*
 * IMP. Image to Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "IMP/engine/meshing/adjacency.hpp"

namespace IMP {


Adjacency::Adjacency() : faces_to_cells_(), edges_to_faces_(), verts_to_edges_()
{}


Adjacency::~Adjacency()
{}


void Adjacency::MapFacesToCells(const std::vector<Triangle> &unique_faces, const std::vector<Tetrahedron> &cells)
{   
    // Clear face to cell adjacency.
    this->faces_to_cells_.clear();
    this->faces_to_cells_.resize(unique_faces.size());

    // Initialize containers for faces with sorted connectivity.
    Triangle sort_face;
    std::vector<Triangle> sorted_faces;  sorted_faces.reserve(unique_faces.size());
    std::vector<int> sorted_face_ids;    sorted_face_ids.reserve(unique_faces.size());
    std::vector<std::pair<Triangle,int>> sorted_faces_pairs; sorted_faces_pairs.reserve(unique_faces.size());

    // Sort in container the faces with their ids.
    int fid = 0;
    for (const auto &face : unique_faces) {
        sort_face = face;
        sort_face.SortVertIds();
        sorted_faces_pairs.emplace_back(std::make_pair(sort_face,fid++));
    }
    std::sort(std::begin(sorted_faces_pairs), std::end(sorted_faces_pairs));

    // Retrieve the sorted faces.
    std::transform(std::begin(sorted_faces_pairs), std::end(sorted_faces_pairs),
                   std::back_inserter(sorted_faces),[](auto const& p){ return p.first; });

    // Retrieve the sorted face ids.
    std::transform(std::begin(sorted_faces_pairs), std::end(sorted_faces_pairs),
                   std::back_inserter(sorted_face_ids),[](auto const& p){ return p.second; });

    // Construct the face to cell adjacency.
    int cid = 0; int pos = 0;
    std::array<Triangle,4> cell_faces;
    for (const auto &cell : cells) {
        cell_faces[0] = Triangle(cell.V0(), cell.V2(), cell.V1(), cid);
        cell_faces[1] = Triangle(cell.V0(), cell.V3(), cell.V2(), cid);
        cell_faces[2] = Triangle(cell.V1(), cell.V2(), cell.V3(), cid);
        cell_faces[3] = Triangle(cell.V0(), cell.V1(), cell.V3(), cid);
        cid++;

        for (const auto &cell_face : cell_faces) {
            // Find position of the similar sorted face.
            sort_face = cell_face;
            sort_face.SortVertIds();
            auto range = std::equal_range(sorted_faces.begin(), sorted_faces.end(), sort_face);

            if (range.first == range.second) {
                std::string error_str = "Could not construct face to cell adjacency. A cell face did not match with any of the unique faces.";
                throw std::runtime_error(Logger::Error(error_str));
            }

            // Index of the face in the sorted faces container.
            pos = range.first-sorted_faces.begin();
            fid = sorted_face_ids[pos];

            // Add parent cell of the current face to the adjacency.
            this->faces_to_cells_[fid].emplace_back(cell_face.ParentCellId());
        }
    }
} 


void Adjacency::MapEdgesToFaces(const std::vector<Edge> &unique_edges, const std::vector<Triangle> &faces)
{   
    // Clear edges to faces adjacency.
    this->edges_to_faces_.clear();
    this->edges_to_faces_.resize(unique_edges.size());

    // Initialize containers for edges with sorted connectivity.
    Edge sort_edge;
    std::vector<Edge> sorted_edges;  sorted_edges.reserve(unique_edges.size());
    std::vector<int> sorted_edge_ids;    sorted_edge_ids.reserve(unique_edges.size());
    std::vector<std::pair<Edge,int>> sorted_edges_pairs; sorted_edges_pairs.reserve(unique_edges.size());

    // Sort in container the edges with their ids.
    int eid = 0;
    for (const auto &edge : unique_edges) {
        sort_edge = edge;
        sort_edge.SortVertIds();
        sorted_edges_pairs.emplace_back(std::make_pair(sort_edge,eid++));
    }
    std::sort(std::begin(sorted_edges_pairs), std::end(sorted_edges_pairs));

    // Retrieve the sorted edges.
    std::transform(std::begin(sorted_edges_pairs), std::end(sorted_edges_pairs),
                   std::back_inserter(sorted_edges),[](auto const& p){ return p.first; });

    // Retrieve the sorted edge ids.
    std::transform(std::begin(sorted_edges_pairs), std::end(sorted_edges_pairs),
                   std::back_inserter(sorted_edge_ids),[](auto const& p){ return p.second; });

    // Construct the edges to faces adjacency.
    int fid = 0; int pos = 0;
    std::array<Edge,3> face_edges;
    for (const auto &face : faces) {
        face_edges[0] = Edge(face.V0(), face.V1(), fid);
        face_edges[1] = Edge(face.V1(), face.V2(), fid);
        face_edges[2] = Edge(face.V2(), face.V0(), fid++);

        for (const auto &face_edge : face_edges) {
            // Find position of the similar sorted edge.
            sort_edge = face_edge;
            sort_edge.SortVertIds();
            auto range = std::equal_range(sorted_edges.begin(), sorted_edges.end(), sort_edge);

            if (range.first == range.second) {
                std::string error_str = "Could not construct edges to faces adjacency. A face edge did not match with any of the unique edges.";
                throw std::runtime_error(Logger::Error(error_str));
            }

            // Index of the edge in the sorted edges container.
            pos = range.first-sorted_edges.begin();
            eid = sorted_edge_ids[pos];

            // Add parent face of the current edge to the adjacency.
            this->edges_to_faces_[eid].emplace_back(face_edge.ParentFaceId());
        }
    }

} 


void Adjacency::MapVertsToEdges(const std::vector<Vertex> &unique_verts, const std::vector<Edge> &edges)
{
    // Clear vertices to edges adjacency.
    this->verts_to_edges_.clear();
    this->verts_to_edges_.resize(unique_verts.size());

    // Create vertices to edges adjacency map.
    int eid=0;
    for (const auto &edge : edges) {
        this->verts_to_edges_[edge.V0()].emplace_back(eid);
        this->verts_to_edges_[edge.V1()].emplace_back(eid++);
    }
}
    

} // End of namespace IMP.
