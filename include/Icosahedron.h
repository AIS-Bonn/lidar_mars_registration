/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration
Copyright (c) 2021, Computer Science Institute VI, University of Bonn.

This file contains adapted code from
https://github.com/caosdoar/spheres/blob/master/src/spheres.cpp
Copyright (c) 2015, Oscar Sebio Cajaraville

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#pragma once
#include <Eigen/Dense>
#include <vector>
// https://github.com/caosdoar/spheres/blob/master/src/spheres.cpp
struct Icosahedron
{
    int subdivided = 0;
    Eigen::Matrix3Xf vertices_mat;
    Eigen::Matrix3Xi triangles_mat;
    std::vector<Eigen::Vector3f> vertices;
    std::vector<Eigen::Vector3i> triangles;
    uint32_t triangleCount() const { return triangles.size(); }

    void addTriangle(uint32_t a, uint32_t b, uint32_t c)
    {
        triangles.emplace_back(Eigen::Vector3i(a,b,c));
    }
    inline Eigen::Vector3f normalize(const Eigen::Vector3f &a) const
    {
        return a.normalized();
    }
    void toMatrices()
    {
        vertices_mat = Eigen::Matrix3Xf(3,vertices.size());
        const int num_vertices = vertices.size();
        for ( int col = 0; col < num_vertices; ++col)
        {
            vertices_mat.col(col) = vertices[col];
        }
        triangles_mat = Eigen::Matrix3Xi(3,triangles.size());
        const int num_triangles = triangles.size();
        for ( int col = 0; col < num_triangles; ++col)
        {
            triangles_mat.col(col) = triangles[col];
        }
    }

    Icosahedron( const int num_subdivision = 0 )
    {
        subdivided = 0;
        const double t = (1.0 + sqrt(5.0)) / 2.0;

        // Vertices
        vertices.emplace_back(normalize(Eigen::Vector3f(-1.0,  t, 0.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f( 1.0,  t, 0.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f(-1.0, -t, 0.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f( 1.0, -t, 0.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f(0.0, -1.0,  t)));
        vertices.emplace_back(normalize(Eigen::Vector3f(0.0,  1.0,  t)));
        vertices.emplace_back(normalize(Eigen::Vector3f(0.0, -1.0, -t)));
        vertices.emplace_back(normalize(Eigen::Vector3f(0.0,  1.0, -t)));
        vertices.emplace_back(normalize(Eigen::Vector3f( t, 0.0, -1.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f( t, 0.0,  1.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f(-t, 0.0, -1.0)));
        vertices.emplace_back(normalize(Eigen::Vector3f(-t, 0.0,  1.0)));

        // Faces
        addTriangle(0, 11, 5);
        addTriangle(0, 5, 1);
        addTriangle(0, 1, 7);
        addTriangle(0, 7, 10);
        addTriangle(0, 10, 11);
        addTriangle(1, 5, 9);
        addTriangle(5, 11, 4);
        addTriangle(11, 10, 2);
        addTriangle(10, 7, 6);
        addTriangle(7, 1, 8);
        addTriangle(3, 9, 4);
        addTriangle(3, 4, 2);
        addTriangle(3, 2, 6);
        addTriangle(3, 6, 8);
        addTriangle(3, 8, 9);
        addTriangle(4, 9, 5);
        addTriangle(2, 4, 11);
        addTriangle(6, 2, 10);
        addTriangle(8, 6, 7);
        addTriangle(9, 8, 1);
        toMatrices();

        while ( subdivided < num_subdivision ){ subdivide(); };
    }

    void subdivide()
    {
        ++subdivided;
        const std::vector<Eigen::Vector3i> faces = triangles;
        triangles.clear();
        for ( const Eigen::Vector3i & tri : faces )
        {
            const Eigen::Vector3f & v0 = vertices[tri(0)];
            const Eigen::Vector3f & v1 = vertices[tri(1)];
            const Eigen::Vector3f & v2 = vertices[tri(2)];
            const Eigen::Vector3f v3 = normalize(0.5 * (v0 + v1));
            const Eigen::Vector3f v4 = normalize(0.5 * (v1 + v2));
            const Eigen::Vector3f v5 = normalize(0.5 * (v2 + v0));
            const int v3i = vertices.size();vertices.emplace_back(v3);
            const int v4i = vertices.size();vertices.emplace_back(v4);
            const int v5i = vertices.size();vertices.emplace_back(v5);
            addTriangle(tri(0),v3i,v5i);
            addTriangle(v3i,tri(1),v4i);
            addTriangle(v5i,v4i,tri(2));
            addTriangle(v3i,v4i,v5i);
        }
        toMatrices();
    }

    void getTriangles ( Eigen::Matrix3Xi & indices ) const
    {
        indices = triangles_mat;
    }
    void getEllipsoid ( const Eigen::Matrix3f & evec, const Eigen::Vector3f & evals, const Eigen::Vector3f & mean, Eigen::Matrix3Xf & vert) const
    {
        const Eigen::Matrix3f rm = evec*evals.asDiagonal();
        const Eigen::Matrix3Xf tv = (rm*vertices_mat).colwise() + mean;
        vert = tv;
    }
};
