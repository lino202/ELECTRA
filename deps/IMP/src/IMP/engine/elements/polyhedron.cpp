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


#include "IMP/engine/elements/polyhedron.hpp"

namespace IMP {


Vertex Polyhedron::Centroid(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const
{
    std::unordered_set<int> unique_vert_ids;
    for (const auto &fid : this->face_ids_) {
        for (const auto &vid : faces[fid].VertIds()) {
            unique_vert_ids.insert(vid);
        }
    }

    Vertex centroid_coords;
    for (const auto &uvid : unique_vert_ids) centroid_coords += verts[uvid];
    return centroid_coords / unique_vert_ids.size();
}


TetraMesh Polyhedron::SymmetricDecomposition(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const
{   
    // Compute the polyhedron's centroid.
    Vertex cell_center = this->Centroid(verts, faces);

    // Map unique ids of vertices and faces of the polyhedron.
    int new_vid = 0, new_fid = 0;
    std::unordered_map<int,int> unique_vert_ids; unique_vert_ids.reserve(this->face_ids_.size());
    std::unordered_map<int,int> unique_face_ids; unique_face_ids.reserve(this->face_ids_.size());
    for (const auto &fid : this->face_ids_) {   
        auto added_face_id = unique_face_ids.insert({fid, new_fid});
        if (added_face_id.second) new_fid++;

        for (const auto &vid : faces[fid].VertIds()) {
            auto added_vert_id = unique_vert_ids.insert({vid, new_vid});
            if (added_vert_id.second) new_vid++;
        }
    }
    auto unique_verts_num = unique_vert_ids.size();
    auto unique_faces_num = unique_face_ids.size();

    // Get the vertices of the decomposition.
    std::vector<Vertex> decomp_verts(unique_vert_ids.size()+unique_face_ids.size()+1);
    
    // Decomposition vertices at edge corners.
    for (const auto &vid : unique_vert_ids) decomp_verts[vid.second] = verts[vid.first];
    
    // Decomposition vertices at faces centroids.
    for (const auto &fid : unique_face_ids) decomp_verts[unique_verts_num+fid.second] = faces[fid.first].Centroid(verts);
    
    // Decomposition vertex at cell centroid.
    decomp_verts[decomp_verts.size()-1] = cell_center;

    // Establish the connectivities of the symmetric decomposition tetrahedra.
    std::vector<Tetrahedron> decomp_tetras; decomp_tetras.reserve(5*this->face_ids_.size());
    std::size_t vnum = 0;
    int v0 = 0, v1 = 0, v2 = 0, v3 = 0;
    Tetrahedron dctet;
    for (const auto &fid : this->face_ids_) {       
        // Get the number of the vertices in the polyhedron's face.
        vnum = faces[fid].VertIds().size();

        // Iterate edges of the polyhedron's face.
        v1 = unique_vert_ids.at(faces[fid].VertIds()[vnum-1]);
        for (std::size_t i=0; i != vnum; ++i) {
            v0 = v1;
            v1 = unique_vert_ids.at(faces[fid].VertIds()[i]);
            v2 = unique_verts_num + unique_face_ids.at(fid);
            v3 = unique_verts_num + unique_faces_num;
            dctet.SetVertIds(v0,v1,v3,v2);
            if (dctet.SignedVolume(decomp_verts) < 0.) dctet.SetVertIds(v1,v0,v3,v2);
            decomp_tetras.emplace_back(dctet);
        }
    }
    decomp_tetras.shrink_to_fit();

    // std::string filename = "/home/mood/Desktop/decomp/decomp_cell.vtk";
    // std::ofstream out(filename);
    // out << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    // out << "POINTS " << decomp_verts.size() << " float\n";
    // for (const auto &vert : decomp_verts) { out << vert.X() << " " << vert.Y() << " " << vert.Z() << "\n"; }
    // out << "CELLS " << decomp_tetras.size() << " " << 5*decomp_tetras.size() << "\n";
    // for (const auto &tet : decomp_tetras) { out << "4 " << tet.V0() << " " << tet.V1() << " " << tet.V2() << " " << tet.V3() << "\n"; }
    // out << "CELL_TYPES " << decomp_tetras.size() << "\n";
    // std::size_t i = 0;
    // while (i != decomp_tetras.size()) { out << "10\n"; i++; }    
    // out.close();

    return TetraMesh(decomp_verts, decomp_tetras);
}


double Polyhedron::SignedVolume(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const
{
    double volume = 0.;
    TetraMesh decomp = this->SymmetricDecomposition(verts, faces);
    for (const auto &tet : decomp.Tetras()) volume += tet.SignedVolume(decomp.Verts());
    return volume;
}


double Polyhedron::AbsVolume(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const
{
    return std::abs(this->SignedVolume(verts, faces));
}


double Polyhedron::ConditionNumber(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const
{
    double cond = 0.;
    Vertex e1, e2, e3;
    Vertex cross12, cross23, cross31;

    // Get symmetric decomposition of the polyhedron.
    auto decomp = this->SymmetricDecomposition(verts,faces);

    // Compute the condition number of each trivalent corner.
    double A = 0., V = 0., L = 0.;
    int corner_id = 0;
    for (const auto &tet : decomp.Tetras()) {
        // Compute condition number for both trivalent corners on polyhedron edges.
        for (int i = 0; i != 2; ++i) {
            // Get edges emanating from the trivalent corner.
            if (i==0) {
                corner_id = tet.V0();
                e1 = (decomp.Verts()[tet.V1()] - decomp.Verts()[corner_id]);
            } else {
                corner_id = tet.V1();
                e1 = (decomp.Verts()[tet.V0()] - decomp.Verts()[corner_id]);
            }
            e2 = decomp.Verts()[tet.V2()] - decomp.Verts()[corner_id];
            e3 = decomp.Verts()[tet.V3()] - decomp.Verts()[corner_id];

            cross12 = Cross(e1,e2);
            cross23 = Cross(e2,e3);
            cross31 = Cross(e3,e1);

            L = std::sqrt(e1.Norm()*e1.Norm() + e2.Norm()*e2.Norm() + e3.Norm()*e3.Norm());
            A = std::sqrt(cross12.Norm()*cross12.Norm() + cross23.Norm()*cross23.Norm() + cross31.Norm()*cross31.Norm());
            V = std::abs(Dot(Cross(e1,e2),e3));

            // Add in the condition number of the polyhedron.
            cond += (L*A)/V;
        }
    }

    // Return normalized condition number.
    return cond / (6*decomp.Tetras().size());

}


bool Polyhedron::IsValid(const std::vector<Vertex> &verts, const std::vector<Polygon> &faces) const
{
    TetraMesh decomp = this->SymmetricDecomposition(verts, faces);
    for (const auto &tet : decomp.Tetras()) {
        if (!(tet.SignedVolume(decomp.Verts()) > 0.)) return false;
    }
    return true;
}


bool operator == (const Polyhedron &poly1, const Polyhedron &poly2)
{
    return (poly1.FaceIds() == poly2.FaceIds());    
}

    
bool operator != (const Polyhedron &poly1, const Polyhedron &poly2)
{ 
    return !(poly1 == poly2); 
}


} // End of namespace IMP.