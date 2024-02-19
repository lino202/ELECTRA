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


#include "IMP/engine/tesselations/vertex_optimizer.hpp"


namespace IMP {


double VertexOptimizer::ComputeObjectiveFunction(const Vertex &v)
{
    double fx = 0.;
    Vertex e1, e2, e3;
    Vertex cross12, cross23, cross31;
    double A = 0., V = 0., L = 0., d = 0.;
    for (const auto &corner : this->trival_corners_) {
        e1 = corner[0] - v;
        e2 = corner[1] - v;
        e3 = corner[2] - v;

        cross12 = Cross(e1,e2);
        cross23 = Cross(e2,e3);
        cross31 = Cross(e3,e1);

        L = std::sqrt(e1.Norm()*e1.Norm() + e2.Norm()*e2.Norm() + e3.Norm()*e3.Norm());
        A = std::sqrt(cross12.Norm()*cross12.Norm() + cross23.Norm()*cross23.Norm() + cross31.Norm()*cross31.Norm());
        V = std::abs(Dot(Cross(e1,e2),e3));
        d = 1. / (1. + 50.*std::exp(V/this->neigh_volume_));

        // Add the quality metric of the corner in the objective function's value.
        fx += (2*A*L) / (V + std::sqrt(V*V + d*d));
    }
    return fx/static_cast<double>(this->trival_corners_.size());
}


void VertexOptimizer::ComputeGradient(const Vertex &v, Eigen::Vector3d &grad)
{
    double h = std::sqrt(std::numeric_limits<double>::epsilon());
    Vertex vxp(v.X()+h, v.Y(),   v.Z());
    Vertex vyp(v.X(),   v.Y()+h, v.Z());
    Vertex vzp(v.X(),   v.Y(),   v.Z()+h);

    Vertex vxm(v.X()-h, v.Y(),   v.Z());
    Vertex vym(v.X(),   v.Y()-h, v.Z());
    Vertex vzm(v.X(),   v.Y(),   v.Z()-h);

    double gfxp = this->ComputeObjectiveFunction(vxp);
    double gfyp = this->ComputeObjectiveFunction(vyp);
    double gfzp = this->ComputeObjectiveFunction(vzp);

    double gfxm = this->ComputeObjectiveFunction(vxm);
    double gfym = this->ComputeObjectiveFunction(vym);
    double gfzm = this->ComputeObjectiveFunction(vzm); 

    grad.coeffRef(0) = (gfxp-gfxm) / (2.*h);
    grad.coeffRef(1) = (gfyp-gfym) / (2.*h);
    grad.coeffRef(2) = (gfzp-gfzm) / (2.*h);

}


void VertexOptimizer::ComputeHessian(const Eigen::Vector3d &grad, Eigen::Matrix3d &hess)
{
    Eigen::Vector3d grad1, grad2, grad3;
    this->ComputeGradient(grad.coeff(0), grad1);
    this->ComputeGradient(grad.coeff(1), grad2);
    this->ComputeGradient(grad.coeff(2), grad3);

    hess.col(0) = grad1;
    hess.col(1) = grad1;
    hess.col(2) = grad1;

}


bool VertexOptimizer::IsPositiveDefinite(const Eigen::Matrix3d &hess)
{   
    if (hess.coeffRef(0,0) < 0.) return false;

    Eigen::Matrix2d subhess;
    subhess.coeffRef(0,0) = hess.coeff(0,0);
    subhess.coeffRef(0,1) = hess.coeff(0,1);
    subhess.coeffRef(1,0) = hess.coeff(1,0);
    subhess.coeffRef(1,1) = hess.coeff(1,1);

    if (hess.coeffRef(0,0) > 0. &&
        subhess.determinant() > 0. &&
        hess.determinant() > 0.) return true;

    return false;
}


void VertexOptimizer::MakePositiveDefinite(Eigen::Matrix3d &hess)
{   
    double beta = 0.001;
    double min_diag = 0., tau = 0.;
    for (int i = 0; i != hess.rows(); ++i) {
        if (hess.coeffRef(i,i) < min_diag) min_diag = hess.coeffRef(i,i);
    }
    if (min_diag < 0.) tau = -min_diag + beta;

    while (!this->IsPositiveDefinite(hess)) {
        hess = hess + beta*Eigen::Matrix3d::Identity();
        tau = std::max(2*tau, beta);
    }
}


bool VertexOptimizer::SatisfiesArmijo(const Vertex &v, const Eigen::Vector3d &grad, const Eigen::Vector3d &p, double a)
{
    Vertex vnew(v.X()+a*p.coeff(0), v.Y()+a*p.coeff(1), v.Z()+a*p.coeff(2));

    double lhs = this->ComputeObjectiveFunction(vnew);
    double rhs = this->ComputeObjectiveFunction(v) + 0.0001*a*grad.transpose()*p;    
    if (lhs <= rhs) return true;

    return false;

}


void VertexOptimizer::Set(const Vertex &vert, const std::vector<std::vector<Vertex>> &trival_corners, double neigh_volume) 
{
    this->vert_ = vert;
    this->trival_corners_ = trival_corners;
    this->neigh_volume_ = neigh_volume;
}


void VertexOptimizer::Optimize()
{
    // Compute objective function.
    Vertex v = this->vert_;
    double a = 1.;
    int max_iter = 100;
    int it = 0;
    Eigen::Vector3d grad;
    Eigen::Matrix3d hess;
    while (it != max_iter) {

        this->ComputeGradient(v, grad);
        if (std::abs(grad.sum()) < 0.00001) break;

        this->ComputeHessian(grad, hess);
        this->MakePositiveDefinite(hess);

        Eigen::Vector3d p = -hess.inverse()*grad;
        while (!this->SatisfiesArmijo(v, grad, p, a))  a *= 0.8;
        
        p *= a;
        v.SetCoords(v.X()+p.coeff(0), v.Y()+p.coeff(1), v.Z()+p.coeff(2)); 
        it++;
    }

    this->vert_ = v;
}



} // End of namespace IMP.
