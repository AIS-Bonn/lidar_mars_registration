/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration

Copyright (c) 2021, Computer Science Institute VI, University of Bonn.
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
#ifdef USE_EASY_PBR

#include <memory>
#include <Eigen/Geometry>
#include "MarsFwdDec.h"

namespace easy_pbr{
    class Mesh;
    typedef std::shared_ptr<Mesh> MeshSharedPtr;
}
typedef easy_pbr::Mesh Cloud;
typedef easy_pbr::MeshSharedPtr CloudPtr;

enum class SurfelType : int
{
    TriangleEllipsoid = 0, Quad, Ellipse
};

class VisMesh
{
private:
    std::vector<Eigen::Vector3f> vertices;
    std::vector<Eigen::Vector3f> colors;
    std::vector<float> intensities;
    std::vector<Eigen::Vector3f> reflectivities;
    std::vector<Eigen::Vector3f> vertex_normals;
    std::vector<Eigen::Vector3f> tangent_direction;
    std::vector<float> tangent_other_length;
    std::vector<float> tangent_normal_length;
    std::vector<Eigen::Vector2i> edges;
    std::vector<Eigen::Vector3i> faces;
    std::vector<Eigen::Vector3f> face_normals;

    CloudPtr prev_mesh;

public:
    typedef std::shared_ptr<VisMesh> Ptr;
    static Ptr create();
    static Eigen::Vector3f randColor();
    static Eigen::Vector3f colorFromNormal( const Eigen::Vector3f & normal );
    static Eigen::Vector3f getOneOfTenColors( const int & idx );
    static Eigen::Vector3f getPlasma ( const float & s);

    static Eigen::Vector3f black;// = Eigen::Vector3f (0,0,0);
    static Eigen::Vector3f red;// = Eigen::Vector3f  (1,0,0);
    static Eigen::Vector3f green;//= Eigen::Vector3f (0,1,0);
    static Eigen::Vector3f blue;// = Eigen::Vector3f (0,0,1);
    static Eigen::Vector3f yellow;// = Eigen::Vector3f (1,1,0);
    static Eigen::Vector3f white;// = Eigen::Vector3f (1,1,1);

    size_t size() const;

    void resetForAggregation();
    void showAggregatedEdgeMesh ( const Sophus::SE3d & oldLastScenePose, const std::string & meshName, const bool & show = true );
    void showAggregatedPointMesh ( const Sophus::SE3d & oldLastScenePose, const std::string & meshName, const float & randomSampleFactor = 0., const bool & show = true );
    void showAggregatedSurfelTriangleEllipsoidMesh ( const Sophus::SE3d & oldLastScenePose, const std::string & meshName );

    void addPoint( const Eigen::Vector3f & vertex_one,
                   const Eigen::Vector3f & color,
                   const float * intensity = nullptr,
                   const uint16_t * reflectivity = nullptr );
    void addEdge( const Eigen::Vector3f & vertex_one,
                  const Eigen::Vector3f & vertex_two,
                  const Eigen::Vector3f & color );
    void addEdge( const Eigen::Vector3f & vertex_one,
                  const Eigen::Vector3f & vertex_two,
                  const Eigen::Vector3f & color_one,
                  const Eigen::Vector3f & color_two );


    void addSurfel( const Eigen::MatrixXf & evecs,
                    const Eigen::VectorXf & evals,
                    const Eigen::Vector3f & mean,
                    const Eigen::Vector3f & color,
                    const float & sigma = 1,
                    const SurfelType & type = SurfelType::TriangleEllipsoid );

    void addSurfelQuad( const Eigen::MatrixXf & evecs,
                        const Eigen::VectorXf & evals,
                        const Eigen::Vector3f & mean,
                        const Eigen::Vector3f & color,
                        const float & sigma = 3);
    void addSurfelEllipsoid( const Eigen::MatrixXf & evecs,
                             const Eigen::VectorXf & evals,
                             const Eigen::Vector3f & mean,
                             const Eigen::Vector3f & color,
                             const float & sigma = 3);
    void addSurfelTriangleEllipsoid( const Eigen::MatrixXf & evecs,
                                     const Eigen::VectorXf & evals,
                                     const Eigen::Vector3f & mean,
                                     const Eigen::Vector3f & color,
                                     const float & sigma = 1 );

    void addNormal( const Eigen::Vector3f & mean,
                    const Eigen::Vector3f & normal,
                    const Eigen::Vector3f & color );

    CloudPtr getPointMesh( ) const;
    CloudPtr getCubeMesh( const float & cube_size ) const;
    CloudPtr getEdgeMesh( ) const;
    CloudPtr getSurfelMesh( ) const;
    CloudPtr getSurfelEllipsoidMesh( ) const;
    CloudPtr getSurfelTriangleEllipsoidMesh( ) const;
    CloudPtr getEdgeMeshFromPoints( ) const;


    void reservePoints( const int num_points );
    void reserveSurfelTriangleEllipsoid ( const int num_surfels );

    CloudPtr showPointMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d(), const int & point_size = 2 ) const;
    CloudPtr showCubeMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d(), const float & cube_size = 2 ) const;
    CloudPtr showEdgeMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d() ) const;
    CloudPtr showEdgeMeshFromPoints( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d() ) const;
    CloudPtr showSurfelMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d(), const SurfelType & type = SurfelType::TriangleEllipsoid ) const;
    CloudPtr showSurfelQuadMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d() ) const;
    CloudPtr showSurfelEllipsoidMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d() ) const;
    CloudPtr showSurfelTriangleEllipsoidMesh( const std::string & meshName, const bool & rotate90deg = false, const Sophus::SE3d & tf = Sophus::SE3d() ) const;
};
#endif
