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


#include "IMP/engine/meshing/tet_to_poly.hpp"


namespace IMP {

TetToPoly::TetToPoly() :
primal_map_(), dual_map_(),  primal_verts_(), dual_verts_(), primal_cells_(), primal_faces_(),
dual_cells_(), dual_faces_(), primal_edges_(), dual_edges_(), dual_node_sets(), vertex_to_dual_face_id_(), edge_to_dual_face_id_(),
cell_to_dual_point_id_(), face_to_dual_point_id_(), edge_to_dual_point_id_(), vertex_to_dual_point_id_(), feat_angle_(90.)
{}


TetToPoly::~TetToPoly()
{}


void TetToPoly::ClearPrimalMesh()
{
    this->primal_verts_.clear(); this->primal_verts_.shrink_to_fit();
    this->primal_edges_.clear(); this->primal_edges_.shrink_to_fit();
    this->primal_faces_.clear(); this->primal_faces_.shrink_to_fit();
    this->primal_cells_.clear(); this->primal_cells_.shrink_to_fit();
}


void TetToPoly::LoadPrimalMesh(const Mesh<3,4> &primal_mesh)
{   
    // Clear primal mesh.
    this->ClearPrimalMesh();

    // Load primal mesh verts.
    this->primal_verts_.resize(primal_mesh.NodesNum());
    int nid = 0;
    for (const auto &node : primal_mesh.Nodes()) {
        this->primal_verts_[nid++].SetCoords(node[0], node[1], node[2]);
    }

    // Load primal mesh cells.
    this->primal_cells_.resize(primal_mesh.CellsNum());
    int cid = 0;
    for (const auto &cell : primal_mesh.Cells()) {
        this->primal_cells_[cid++].SetVertIds(cell.N(0), cell.N(1), cell.N(2), cell.N(3));
    }

    if (primal_mesh.NodeSetsNum() != 0) this->dual_node_sets = primal_mesh.NodeSets();
}


void TetToPoly::OrientPrimalCells()
{
    for (auto &cell : this->primal_cells_) {
        // Connectivity of the cell.
        auto v0 = cell.V0();
        auto v1 = cell.V1();
        auto v2 = cell.V2();
        auto v3 = cell.V3();

        // Vertices of the cell.
        auto p1 = this->primal_verts_[v0];
        auto p2 = this->primal_verts_[v1];
        auto p3 = this->primal_verts_[v2];
        auto p4 = this->primal_verts_[v3];

        // Edges starting from first vertex to the other vertices of the cell.
        auto e1 = Vertex({p2.X()-p1.X(), p2.Y()-p1.Y(), p2.Z()-p1.Z()});
        auto e2 = Vertex({p3.X()-p1.X(), p3.Y()-p1.Y(), p3.Z()-p1.Z()});
        auto e3 = Vertex({p4.X()-p1.X(), p4.Y()-p1.Y(), p4.Z()-p1.Z()});

        // Compute the determinant of the three vectors.
        auto det = Det(e1, e2, e3);

        // Swap cell connectivity if not in counter-clockwise order.
        if (det < 2*std::numeric_limits<double>::epsilon())  cell.SetVertIds(v2,v1,v0,v3);
    }
}


void TetToPoly::ExtractPrimalFaces()
{
    if (this->primal_cells_.size() == 0) {
        std::string error_str = "Could not extract primal faces. Primal cells are not set.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Reset primal faces container.
    this->primal_faces_.clear();
    this->primal_faces_.reserve(4*this->primal_cells_.size());

    // Get unique faces of the primal cells.
    auto cell_faces = std::array<Triangle,4>{};
    auto unique_faces = std::unordered_set<Triangle, Triangle::HashFunc>{};
    auto sort_face = Triangle{};
    auto cid = int{0};
    for (const auto &cell : this->primal_cells_) {

        cell_faces[0] = Triangle(cell.V0(), cell.V2(), cell.V1(), cid);
        cell_faces[1] = Triangle(cell.V0(), cell.V3(), cell.V2(), cid);
        cell_faces[2] = Triangle(cell.V1(), cell.V2(), cell.V3(), cid);
        cell_faces[3] = Triangle(cell.V0(), cell.V1(), cell.V3(), cid++);

        // Store only the unique faces of the cell.
        for (const auto &face : cell_faces) {
            sort_face = face;
            sort_face.SortVertIds();
            if (unique_faces.insert(sort_face).second == true) {
                this->primal_faces_.emplace_back(face);
            }
        }
    }
    this->primal_faces_.shrink_to_fit();

    // Map the primal faces to all their parent cells.
    this->primal_map_.MapFacesToCells(this->primal_faces_, this->primal_cells_);

    // Set the boundary flag of the faces.
    int fid = 0;
    for (const auto &parent_cell_ids : this->primal_map_.FacesToCells()) {
        if (parent_cell_ids.size() == 1) {
            this->primal_faces_[fid].SetOnBoundary(true);
        }
        fid++;
    }

}


void TetToPoly::ExtractPrimalEdges()
{
    if (this->primal_faces_.size() == 0) {
        std::string error_str = "Could not extract primal edges. Primal faces are not set.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Reset primal edges container.
    this->primal_edges_.clear();
    this->primal_edges_.reserve(3*this->primal_faces_.size());

    // Get unique edges of the primal faces.
    std::array<Edge,3> face_edges;
    std::unordered_set<Edge, Edge::HashFunc> unique_edges;
    Edge sort_edge;
    int fid = 0;
    for (const auto &face : this->primal_faces_) {

        face_edges[0] = Edge(face.V0(), face.V1(), fid);
        face_edges[1] = Edge(face.V1(), face.V2(), fid);
        face_edges[2] = Edge(face.V2(), face.V0(), fid);
        fid++;

        // Store only the unique edges of the cell.
        for (const auto &edge : face_edges) {
            sort_edge = edge;
            sort_edge.SortVertIds();
            if (unique_edges.insert(sort_edge).second == true) {
                this->primal_edges_.emplace_back(edge);
            }
        }
    }
    this->primal_edges_.shrink_to_fit();

    // Map the primal edges to all their parent faces.
    this->primal_map_.MapEdgesToFaces(this->primal_edges_, this->primal_faces_);

    // Set the boundary flag of the edges.
    int eid = 0;
    for (const auto &parent_face_ids : this->primal_map_.EdgesToFaces()) {
        for (const auto &fid : parent_face_ids) {
            if (this->primal_faces_[fid].OnBoundary()) {
                this->primal_edges_[eid].SetOnBoundary(true);
                break;
            }
        }
        eid++;
    }

}


void TetToPoly::FindEdgesOnFeature()
{
    if (this->primal_map_.EdgesToFaces().size() == 0) {
        std::string error_str = "Could not find edges on feature. Extract the primar edges first.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    auto eid = 0;
    for (auto &edge : this->primal_edges_) {
        if (edge.OnBoundary()) {
            // Get the indices of the attached faces to the edge.
            auto face_ids = this->primal_map_.EdgesToFaces()[eid];

            // Get the indices of the two boundary faces.
            auto bnd_fids = std::vector<int>{};
            for (const auto &fid : face_ids) {
                if (this->primal_faces_[fid].OnBoundary()) {
                    bnd_fids.emplace_back(fid);
                }
            }

            if (bnd_fids.size() == 2) {

                // Compute the normals to the faces.
                auto a = this->primal_faces_[bnd_fids[0]].Normal(this->primal_verts_);
                a /= a.Norm();
                auto b = this->primal_faces_[bnd_fids[1]].Normal(this->primal_verts_);
                b /= b.Norm();
                // The edge is on a feature if the normals' angle between is smaller than the feature angle.
                auto dot1 = Dot(a, b);
                if (std::abs(dot1) < std::cos(ALGORITHMS::DegToRad(this->feat_angle_))) edge.SetOnFeature(true);
            }
        }
        eid++;
    }

}


void TetToPoly::FindVertsOnBoundaryAndFeature()
{
    this->primal_map_.MapVertsToEdges(this->primal_verts_, this->primal_edges_);

    int vid = 0; int feat_edge_num = 0;
    for (const auto &edges_on_vert : this->primal_map_.VertsToEdges()) {
        for (const auto &eid : edges_on_vert) {
            if (this->primal_edges_[eid].OnBoundary()) this->primal_verts_[vid].SetOnBoundary(true);
            if (this->primal_edges_[eid].OnFeature()) feat_edge_num++;
        }

        if (feat_edge_num == 2) this->primal_verts_[vid].SetOnFeatureEdge(true);
        if (feat_edge_num > 2) this->primal_verts_[vid].SetOnFeatureCorner(true);
        vid++;
        feat_edge_num = 0;
    }
}


void TetToPoly::GenerateDualVerts()
{
    // Check that all adjacency information on the primal mesh is available.
    if (!this->primal_map_.IsReady()) {
        std::string error_str = "Could not generate dual vertices. Adjacency information of primal map is missing.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    /* RANGES ARE NOT WELL SET */
    this->dual_verts_.clear();

    this->cell_to_dual_point_id_.clear();
    this->face_to_dual_point_id_.clear();
    this->edge_to_dual_point_id_.clear();
    this->vertex_to_dual_point_id_.clear();

    // Create dual vertices in the cells of the primal mesh.
    int dvid = 0; int cid = 0;
    for (const auto &cell : this->primal_cells_) {
        this->dual_verts_.emplace_back(cell.Centroid(this->primal_verts_));
        this->cell_to_dual_point_id_[cid] = dvid;
        dvid++;
        cid++;
    }

    // Create dual vertices on the boundary faces of the primal mesh.
    int fid = 0;
    for (const auto &face : this->primal_faces_) {
        if (face.OnBoundary()) {
            this->dual_verts_.emplace_back(face.Centroid(this->primal_verts_));
            this->face_to_dual_point_id_[fid] = dvid;
            dvid++;
        }
        fid++;
    }

    // Create dual vertices on the feature edges of the primal mesh.
    int eid = 0;
    for (const auto &edge : this->primal_edges_) {
        if (edge.OnFeature()) {
            this->dual_verts_.emplace_back(edge.Centroid(this->primal_verts_));
            this->edge_to_dual_point_id_[eid] = dvid;
            dvid++;
        }
        eid++;
    }

    // Create dual vertices at the feature corners of the primal mesh.
    int vid = 0;
    for (const auto &vert : this->primal_verts_) {
        if (vert.OnFeatureCorner()) {
            this->dual_verts_.emplace_back(vert);
            this->vertex_to_dual_point_id_[vid] = dvid;
            dvid++;
        }
        vid++;
    }
}


void TetToPoly::GenerateDualFaces()
{
    this->dual_faces_.clear();
    this->dual_faces_.reserve(this->primal_edges_.size()+this->primal_verts_.size());
    this->vertex_to_dual_face_id_.clear();
    this->vertex_to_dual_face_id_.reserve(this->primal_verts_.size());
    this->edge_to_dual_face_id_.clear();
    this->edge_to_dual_face_id_.reserve(this->primal_edges_.size());

    // Generate dual faces corresponding to primal edges.
    std::vector<int> edge_faces_ids;                // The indices of the adjacent faces to the edge.
    std::vector<int> bound_face_ids;                // The indices of the adjacent faces to the edge that are on the boundary.
    std::unordered_set<int> edge_cells_ids;         // The indices of the adjacent cells to the edge.
    std::vector<int> dual_face_conn;
    std::vector<bool> visited;
    int eid = 0; int dual_face_id = 0;
    for (const auto &edge : this->primal_edges_) {

        // Get the indices of the faces and cells that are attached to the edge.
        edge_faces_ids = this->primal_map_.EdgesToFaces()[eid];
        for (const auto &fid : edge_faces_ids) {
            // Collect the unique indices of the connected cells to the edge. 
            std::copy(std::begin(this->primal_map_.FacesToCells()[fid]), 
                      std::end(this->primal_map_.FacesToCells()[fid]),
                      std::inserter(edge_cells_ids, std::end(edge_cells_ids)));

            // Collect the indices of the boundary faces connected to the edge, if any.
            if (this->primal_faces_[fid].OnBoundary())  bound_face_ids.emplace_back(fid);
        }

        // Add the index of the first cell in the dual face connectivity.
        dual_face_conn.resize(edge_cells_ids.size());
        visited.resize(edge_cells_ids.size(), false);
        if (bound_face_ids.size() == 2) {
            // Get the index of the attached cell to the first boundary face.    
            int id = this->primal_map_.FacesToCells()[bound_face_ids[0]][0];

            int iter = 0;
            for (const auto &cid : edge_cells_ids) {
                if (cid == id) { dual_face_conn[0] = this->cell_to_dual_point_id_.at(cid);  visited[iter] = true;  break; }
                iter++;
            }
        } else {
            dual_face_conn[0] = this->cell_to_dual_point_id_.at(*edge_cells_ids.begin());
            visited[0] = true;
        }

        // Add the indices of the rest cells in the dual face connectivity.
        std::size_t pos = 0;
        while (pos != dual_face_conn.size()-1) {
            int i = 0;
            for (const auto &cid : edge_cells_ids) {
                if (this->primal_cells_[dual_face_conn[pos]].IsSharingFace(this->primal_cells_[cid]) && !visited[i]) {
                    visited[i] = true;
                    dual_face_conn[++pos] = this->cell_to_dual_point_id_.at(cid);
                    break;
                }
                i++;
            }
        }

        // Add the index of the dual points corresponding to the boundary faces.
        if (bound_face_ids.size() == 2) {
            int fid1 = bound_face_ids[0];
            int fid2 = bound_face_ids[1];
            dual_face_conn.insert(std::begin(dual_face_conn), this->face_to_dual_point_id_.at(fid1));
            dual_face_conn.emplace_back(this->face_to_dual_point_id_.at(fid2));
        }

        // Add the index of the dual point corresponding to the edge.
        if (edge.OnFeature()) { 
            dual_face_conn.emplace_back(this->edge_to_dual_point_id_.at(eid));
        }

        // Add the dual face in the container.
        this->dual_faces_.emplace_back(Polygon(dual_face_conn));
        this->edge_to_dual_face_id_[eid] = dual_face_id++;

        // Clear containers.
        edge_faces_ids.clear();
        bound_face_ids.clear();
        edge_cells_ids.clear();
        visited.clear();
        dual_face_conn.clear();
        eid++;
    } // End of Generate dual faces corresponding to primal edges.


    // Generate dual faces corresponding to primal vertices.
    std::vector<int> vert_edges_ids;               // The indices of the adjacent edges to the vertex.
    std::vector<int> feat_edge_ids;                // The indices of the adjacent edges to the vertex that are on a feature.
    std::unordered_set<int> vert_faces_ids;        // The indices of the adjacent faces to the vertex.
    std::vector<int> vert_dual_faces_ids;

    int vid = 0;
    for (const auto &vert : this->primal_verts_) {

        if (vert.OnBoundary()) {
            // Clear containers.
            vert_edges_ids.clear();  vert_faces_ids.clear();  
            feat_edge_ids.clear();   vert_dual_faces_ids.clear();

            // Get the indices of the edges and faces that are attached to the vertex.
            vert_edges_ids = this->primal_map_.VertsToEdges()[vid];
            for (const auto &edg_id : vert_edges_ids) {
                for (const auto &fid : this->primal_map_.EdgesToFaces()[edg_id]) {
                    if (this->primal_faces_[fid].OnBoundary()) vert_faces_ids.insert(fid);
                }
                // Collect the indices of the feature edges connected to the vertex, if any.
                if (this->primal_edges_[edg_id].OnFeature())  feat_edge_ids.emplace_back(edg_id);
            }

            // Generate dual faces for boundary vertices on features. 
            if (feat_edge_ids.size() > 1) {
                std::unordered_map<int,int> faces_pos;
                std::vector<std::vector<int>> feat_faces(vert_faces_ids.size());
                std::vector<Vertex> face_normals(vert_faces_ids.size());
                int i = 0; 
                for (const auto &fid : vert_faces_ids) {
                    // Compute normal of face.
                    face_normals[i] = this->primal_faces_[fid].Normal(this->primal_verts_);

                    // Search the feature edges of the face.
                    for (const auto &feat_eid : feat_edge_ids) {
                        if (this->primal_faces_[fid].HasEdge(this->primal_edges_[feat_eid])) 
                            feat_faces[i].emplace_back(feat_eid);
                    }
                    // Map the index of the face to its position in the container.
                    faces_pos[fid] = i++;
                }

                // Set all the faces not visited.
                visited.clear();  visited.resize(vert_faces_ids.size(), false);

                // Generate dual faces for different feature faces.
                int ff_pos = 0; int ff_id = 0;
                for (const auto &ff : feat_faces) {
                    ff_pos = &ff - &feat_faces[0];
                    ff_id = *std::next(vert_faces_ids.begin(),ff_pos);
                    if (!visited[faces_pos.at(ff_id)]) {
                        visited[faces_pos.at(ff_id)] = true;

                        if (ff.size() == 2) {
                            // Create dual face connectivity.
                            dual_face_conn.clear();
                            dual_face_conn.emplace_back(this->edge_to_dual_point_id_.at(ff[0]));
                            dual_face_conn.emplace_back(this->face_to_dual_point_id_.at(ff_id));
                            dual_face_conn.emplace_back(this->edge_to_dual_point_id_.at(ff[1]));
                            if (vert.OnFeatureCorner()) dual_face_conn.emplace_back(this->vertex_to_dual_point_id_.at(vid));

                            // Add the dual face in the container.
                            this->dual_faces_.emplace_back(Polygon(dual_face_conn));
                            vert_dual_faces_ids.emplace_back(dual_face_id++);

                        } else if (ff.size() == 1) {
                            // Create dual face connectivity.
                            dual_face_conn.clear();
                            dual_face_conn.emplace_back(this->edge_to_dual_point_id_.at(ff[0]));
                            dual_face_conn.emplace_back(this->face_to_dual_point_id_.at(ff_id));

                            // Get all the faces that are in the same feature.
                            std::vector<int> same_feat_face_ids;
                            auto a = face_normals[ff_pos];
                            a /= a.Norm();
                            for (const auto &fid : vert_faces_ids) {
                                auto b = face_normals[faces_pos.at(fid)];
                                b /= b.Norm();
                                auto dot1 = Dot(a,b);
                                if (fid != ff_id && std::abs(dot1) > std::cos(ALGORITHMS::DegToRad(this->feat_angle_))) {
                                    same_feat_face_ids.emplace_back(fid);
                                }
                            }

                            int prev_fid = ff_id;
                            if (same_feat_face_ids.size() == 1) {
                                dual_face_conn.emplace_back(this->face_to_dual_point_id_.at(same_feat_face_ids[0]));
                                prev_fid = same_feat_face_ids[0];
                            } else {
                                std::size_t pass = 0; int safe_guard = 0;
                                std::vector<bool> processed(same_feat_face_ids.size(), false);
                                while (pass != same_feat_face_ids.size()) {
                                    if (safe_guard == 100) { std::cout << "During dual face generation: reach safeguard\n"; break; }
                                    int it = 0;
                                    for (const auto &sfid : same_feat_face_ids) {
                                        if (!processed[it] && this->primal_faces_[prev_fid].IsSharingEdge(this->primal_faces_[sfid])) {
                                            dual_face_conn.emplace_back(this->face_to_dual_point_id_.at(sfid));
                                            prev_fid = sfid;
                                            processed[it] = true;
                                            visited[faces_pos.at(sfid)] = true;
                                            pass++;
                                        }
                                        it++;
                                        safe_guard++;
                                    }
                                }
                            }

                            bool found_end_edge = false;
                            for (const auto &feat_eid : feat_edge_ids) {
                                if (this->primal_faces_[prev_fid].HasEdge(this->primal_edges_[feat_eid])) {
                                    dual_face_conn.emplace_back(this->edge_to_dual_point_id_.at(feat_eid));
                                    found_end_edge = true;
                                }
                            }

                            if (!found_end_edge) std::cout << "During dual face generation: not found final edge\n";

                            if (vert.OnFeatureCorner()) dual_face_conn.emplace_back(this->vertex_to_dual_point_id_.at(vid));

                            // Add the dual face in the container.
                            this->dual_faces_.emplace_back(Polygon(dual_face_conn));
                            vert_dual_faces_ids.emplace_back(dual_face_id++);


                        }
                    }
                } // End Generate dual faces for different feature faces.

              
            } else { // Generate dual faces for the rest boundary vertices.
                dual_face_conn.clear(); dual_face_conn.resize(vert_faces_ids.size());
                visited.clear(); visited.resize(vert_faces_ids.size(), false);
                int prev_fid = *vert_faces_ids.begin();
                dual_face_conn[0] = this->face_to_dual_point_id_.at(prev_fid);
                visited[0] = true;

                // Add the indices of the rest faces in the dual face connectivity.
                std::size_t pos = 0;
                while (pos != dual_face_conn.size()-1) {
                    int i = 0;
                    for (const auto &fid : vert_faces_ids) {
                        if (this->primal_faces_[prev_fid].IsSharingEdge(this->primal_faces_[fid]) && !visited[i]) {
                            visited[i] = true;
                            dual_face_conn[++pos] = this->face_to_dual_point_id_.at(fid);
                            prev_fid = fid;
                            break;
                        }
                        i++;
                    }
                }
                // Add the dual face in the container.
                this->dual_faces_.emplace_back(Polygon(dual_face_conn));
                vert_dual_faces_ids.emplace_back(dual_face_id++);
            }
            this->vertex_to_dual_face_id_[vid] = vert_dual_faces_ids;
        }
        vid++;
    }

}


void TetToPoly::GenerateDualCells()
{   
    this->dual_cells_.clear();
    this->dual_cells_.resize(this->primal_verts_.size());

    // Constuct a polyhedral cell for each vertex of the primal mesh.
    int vid = 0;
    std::vector<int> polygon_ids;
    for (const auto &primal_vert : this->primal_verts_) {
        
        // Clear poly_faces container.
        polygon_ids.clear();

        // Collect dual faces of connected edges.
        auto vert_edges_ids = this->primal_map_.VertsToEdges()[vid];
        polygon_ids.clear();  polygon_ids.reserve(vert_edges_ids.size());
        int fid = 0;
        for (const auto &eid : vert_edges_ids) {
            fid = &this->dual_faces_[this->edge_to_dual_face_id_.at(eid)] - &this->dual_faces_[0];
            polygon_ids.emplace_back(fid);
        }

        // Add capping dual faces.
        if (primal_vert.OnBoundary()) {
            for (const auto &cap_df : this->vertex_to_dual_face_id_.at(vid)) {
                fid = &this->dual_faces_[cap_df] - &this->dual_faces_[0];
                polygon_ids.emplace_back(fid);
            }
        }

        // Construct dual cell from its poly faces.
        this->dual_cells_[vid++] = Polyhedron(polygon_ids);
    }
}


void TetToPoly::SavePrimalFaces(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << this->primal_verts_.size() << " float\n";
    for (const auto &vert : this->primal_verts_) {
        out << vert.X() << " " << vert.Y() << " " << vert.Z() << std::endl;
    }
    
    out << "CELLS " << this->primal_faces_.size() << " " << 4*this->primal_faces_.size() << "\n";
    for (const auto &face : this->primal_faces_) {
        out << "3 " << face.V0() << " " << face.V1() << " " << face.V2() << "\n";
    }
    out << "CELL_TYPES " << this->primal_faces_.size() << "\n";
    std::size_t i = 0;
    while (i != this->primal_faces_.size()) {
        out << "5\n";
        i++;
    }

    // Distinguish faces on the boundary.
    out << "CELL_DATA " << this->primal_faces_.size() << "\n";
    out << "SCALARS Boundary int 1\n";
    out << "LOOKUP_TABLE default\n";
    int pad = 0;
    for (const auto &face : this->primal_faces_) {
        face.OnBoundary() ? out << "1" : out << "0";
        
        if (pad == 20) { out << "\n"; pad = 0; }
        else { out << " "; pad++; }
    }
    out << "\n";

    out << "NORMALS cell_normals float\n";
    for (const auto &face : this->primal_faces_) {
        auto normal = face.Normal(this->primal_verts_);
        out << normal.X() << " " << normal.Y() << " " << normal.Z() << std::endl;
    }
    
    out.close();
}


void TetToPoly::SavePrimalEdges(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << this->primal_verts_.size() << " float\n";
    for (const auto &vert : this->primal_verts_) {
        out << vert.X() << " " << vert.Y() << " " << vert.Z() << std::endl;
    }
    
    out << "CELLS " << this->primal_edges_.size() << " " << 3*this->primal_edges_.size() << "\n";
    for (const auto &edge : this->primal_edges_) {
        out << "2 " << edge.V0() << " " << edge.V1() << "\n";
    }

    out << "CELL_TYPES " << this->primal_edges_.size() << "\n";
    std::size_t i = 0;
    while (i != this->primal_edges_.size()) {
        out << "3\n";
        i++;
    }

    // Distinguish edges on the boundary and on features.
    out << "CELL_DATA " << this->primal_edges_.size() << "\n";
    out << "SCALARS Boundary int 1\n";
    out << "LOOKUP_TABLE default\n";
    int pad = 0;
    for (const auto &edge : this->primal_edges_) {
        edge.OnBoundary() ? out << "1" : out << "0";
        
        if (pad == 20) { out << "\n"; pad = 0; }
        else { out << " "; pad++; }
    }
    out << "\n";
    out << "SCALARS Features int 1\n";
    out << "LOOKUP_TABLE default\n";
    pad = 0;
    for (const auto &edge : this->primal_edges_) {
        edge.OnFeature() ? out << "1" : out << "0";
        
        if (pad == 20) { out << "\n"; pad = 0; }
        else { out << " "; pad++; }
    }
    
    out.close();
}


void TetToPoly::SavePrimalVerts(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << this->primal_verts_.size() << " float\n";
    for (const auto &vert : this->primal_verts_) {
        out << vert.X() << " " << vert.Y() << " " << vert.Z() << std::endl;
    }
    
    out << "CELLS " << this->primal_verts_.size() << " " << 2*this->primal_verts_.size() << "\n";
    std::size_t vid = 0;
    while (vid != this->primal_verts_.size()) {
        out << "1 " << vid++ << "\n";
    }

    out << "CELL_TYPES " << this->primal_verts_.size() << "\n";
    vid = 0;
    while (vid != this->primal_verts_.size()) {
        out << "1\n";
        vid++;
    }

    // Distinguish vertices on the boundary and on feature corners.
    out << "CELL_DATA " << this->primal_verts_.size() << "\n";
    out << "SCALARS Boundary int 1\n";
    out << "LOOKUP_TABLE default\n";
    int pad = 0;
    for (const auto &vert : this->primal_verts_) {
        vert.OnBoundary() ? out << "1" : out << "0";
        
        if (pad == 20) { out << "\n"; pad = 0; }
        else { out << " "; pad++; }
    }
    out << "\n";
    out << "SCALARS Corners int 1\n";
    out << "LOOKUP_TABLE default\n";
    pad = 0;
    for (const auto &vert : this->primal_verts_) {
        vert.OnFeatureCorner() ?  out << "1" : out << "0";
        
        if (pad == 20) { out << "\n"; pad = 0; }
        else { out << " "; pad++; }
    }
    
    out.close();
}


void TetToPoly::SaveDualVerts(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << this->dual_verts_.size() << " float\n";
    for (const auto &vert : this->dual_verts_) {
        out << vert.X() << " " << vert.Y() << " " << vert.Z() << std::endl;
    }
    
    out << "CELLS " << this->dual_verts_.size() << " " << 2*this->dual_verts_.size() << "\n";
    std::size_t vid = 0;
    while (vid != this->dual_verts_.size()) {
        out << "1 " << vid++ << "\n";
    }

    out << "CELL_TYPES " << this->dual_verts_.size() << "\n";
    vid = 0;
    while (vid != this->dual_verts_.size()) {
        out << "1\n";
        vid++;
    }
    
    out.close();
}


void TetToPoly::SaveDualFaces(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    out << "POINTS " << this->dual_verts_.size() << " float\n";
    for (const auto &vert : this->dual_verts_) {
        out << vert.X() << " " << vert.Y() << " " << vert.Z() << std::endl;
    }
    
    int dual_face_conn_entries = 0;
    for (const auto &face : this->dual_faces_) dual_face_conn_entries += (face.VertIds().size()+1);

    out << "POLYGONS " << this->dual_faces_.size() << " " << dual_face_conn_entries << "\n";
    for (const auto &face : this->dual_faces_) {
        out << face.VertIds().size();
        for (const auto &id : face.VertIds()) out << " " << id;
        out << "\n";
    }
    
    out.close();
}


void TetToPoly::SaveDualCells(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    out << "POINTS " << this->dual_verts_.size() << " float\n";
    for (const auto &vert : this->dual_verts_) {
        out << vert.X() << " " << vert.Y() << " " << vert.Z() << std::endl;
    }

    int dual_faces_num = 0; int dual_face_conn_entries = 0;
    for (const auto &cell : this->dual_cells_) {
        for (const auto &fid : cell.FaceIds()) {
            dual_faces_num++;
            dual_face_conn_entries += (this->dual_faces_[fid].VertIds().size()+1);
        }
    }
    out << "POLYGONS " << dual_faces_num << " " << dual_face_conn_entries << "\n";
    for (const auto &cell : this->dual_cells_) {
        for (const auto &fid : cell.FaceIds()) {
            out << this->dual_faces_[fid].VertIds().size();
            for (const auto &id : this->dual_faces_[fid].VertIds()) out << " " << id;
            out << "\n";
        }
    }
    
    out.close();
}


void TetToPoly::SaveDualMesh(const std::string &outfile)
{
    std::ofstream out(outfile);

    out << "*\n** ELECTRA Voronoi Tesselation file\n";
    out << "** Generated by: IMP library\n**\n";
    
    // Save primal mesh vertices (nodes).
    out << "*Nodes " << this->primal_verts_.size() << "\n";
    int id=0;
    for (const auto &n : this->primal_verts_) {
        out << ++id << ", " << std::setprecision(15) << n.X() << ", " << n.Y() << ", " << n.Z() << "\n";
    }

    // Save dual mesh vertices (points).
    out << "*Points " << this->dual_verts_.size() << "\n";
    id = 0;
    for (const auto &p : this->dual_verts_) {
        out << ++id << ", " << std::setprecision(15) << p.X() << ", " << p.Y() << ", " << p.Z() << "\n";
    }

    // Save dual mesh faces (polygons).
    out << "*Facets " << this->dual_faces_.size() << "\n";
    id = 0;
    for (const auto &polygon : this->dual_faces_) {
        out << ++id << ", " << polygon.VertIds().size();
        for (const auto &pid : polygon.VertIds())  out << ", " << pid+1;
        out << "\n";
    }

    // Save dual mesh cells (polyhedra)
    out << "*Cells " << this->dual_cells_.size() << "\n";
    id = 0;
    for (const auto &polyhedron : this->dual_cells_) {
        out << ++id << ", " << polyhedron.FaceIds().size();
        for (const auto &fid : polyhedron.FaceIds())  out << ", " << fid+1;
        out << "\n";
    }

    // Save dual mesh node sets if any.
    if (this->dual_node_sets.size() != 0) {
        int pad , cnt;
        for (const auto &nset : this->dual_node_sets) {
            out << "*Nset, nset=" << nset.first << "\n";
            pad = 0;
            cnt = 0;
            for (const auto &nid : nset.second.NodeIds()) {
                out << nid+1; pad++; cnt++;
                if (pad == 20 || cnt == static_cast<int>(nset.second.NodeIds().size())) { pad = 0; out << "\n"; }
                else { out << ", "; }    
            }
        }
    }

    out << "*END\n**\n";    
    out.close();
}


} // End of namespace IMP.
