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

#ifdef USE_EASY_PBR
#include "VisUtils.h"
#include "Icosahedron.h"

#include "easy_pbr/Mesh.h"
#include "easy_pbr/Scene.h"
using namespace easy_pbr;
using namespace radu::utils;

#define LOGURU_REPLACE_GLOG 1
#include <loguru.hpp>

Icosahedron iso ( 1 );

Eigen::Vector3f VisMesh::black = Eigen::Vector3f (0,0,0);
Eigen::Vector3f VisMesh::red = Eigen::Vector3f  (1,0,0);
Eigen::Vector3f VisMesh::green= Eigen::Vector3f (0,1,0);
Eigen::Vector3f VisMesh::blue = Eigen::Vector3f (0,0,1);
Eigen::Vector3f VisMesh::yellow = Eigen::Vector3f (1,1,0);
Eigen::Vector3f VisMesh::white = Eigen::Vector3f (1,1,1);

Eigen::VectorXf vec2eigen ( const std::vector<float> & vec )
{
    if ( vec.empty() ) return Eigen::VectorXf();
    Eigen::VectorXf mat = Eigen::VectorXf::Zero(vec.size());
    for ( int row = 0; row < int(vec.size()); ++row )
        mat(row) = vec[row];
    return mat;
}
Eigen::MatrixXf vec2eigen ( const std::vector<Eigen::VectorXf> & vec )
{
    if ( vec.empty() ) return Eigen::MatrixXf();
    Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(vec.size(),vec.front().rows());
    for ( int row = 0; row < int(vec.size()); ++row )
        mat.row(row) = vec[row].transpose();
    return mat;
}
Eigen::MatrixXf vec2eigen ( const std::vector<Eigen::Vector3f> & vec )
{
    if ( vec.empty() ) return Eigen::MatrixXf(0,3);
    Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(vec.size(),3);
    for ( int row = 0; row < int(vec.size()); ++row )
        mat.row(row) = vec[row].transpose();
    return mat;
}

Eigen::MatrixXi vec2eigen ( const std::vector<Eigen::VectorXi> & vec )
{
    if ( vec.empty() ) return Eigen::MatrixXi();
    Eigen::MatrixXi mat = Eigen::MatrixXi::Zero(vec.size(),vec.front().rows());
    for ( int row = 0; row < int(vec.size()); ++row )
        mat.row(row) = vec[row].transpose();
    return mat;
}

Eigen::MatrixXi vec2eigen ( const std::vector<Eigen::Vector2i> & vec )
{
    if ( vec.empty() ) return Eigen::MatrixXi();
    Eigen::MatrixXi mat = Eigen::MatrixXi::Zero(vec.size(),2);
    for ( int row = 0; row < int(vec.size()); ++row )
        mat.row(row) = vec[row].transpose();
    return mat;
}

Eigen::MatrixXi vec2eigen ( const std::vector<Eigen::Vector3i> & vec )
{
    if ( vec.empty() ) return Eigen::MatrixXi();
    Eigen::MatrixXi mat = Eigen::MatrixXi::Zero(vec.size(),3);
    for ( int row = 0; row < int(vec.size()); ++row )
        mat.row(row) = vec[row].transpose();
    return mat;
}


VisMesh::Ptr VisMesh::create()
{
    return std::make_shared<VisMesh>();
}

Eigen::Vector3f VisMesh::randColor()
{
    return Eigen::Vector3f::Random().cwiseAbs();
}

Eigen::Vector3f VisMesh::getPlasma ( const float & s)
{
    // plasma map from https://github.com/yuki-koyama/tinycolormap/blob/master/include/tinycolormap.hpp
    static Eigen::Matrix3Xf plasma = (Eigen::Matrix<float,256,3>() <<
         0.050383, 0.029803, 0.527975 ,
         0.063536, 0.028426, 0.533124 ,
         0.075353, 0.027206, 0.538007 ,
         0.086222, 0.026125, 0.542658 ,
         0.096379, 0.025165, 0.547103 ,
         0.105980, 0.024309, 0.551368 ,
         0.115124, 0.023556, 0.555468 ,
         0.123903, 0.022878, 0.559423 ,
         0.132381, 0.022258, 0.563250 ,
         0.140603, 0.021687, 0.566959 ,
         0.148607, 0.021154, 0.570562 ,
         0.156421, 0.020651, 0.574065 ,
         0.164070, 0.020171, 0.577478 ,
         0.171574, 0.019706, 0.580806 ,
         0.178950, 0.019252, 0.584054 ,
         0.186213, 0.018803, 0.587228 ,
         0.193374, 0.018354, 0.590330 ,
         0.200445, 0.017902, 0.593364 ,
         0.207435, 0.017442, 0.596333 ,
         0.214350, 0.016973, 0.599239 ,
         0.221197, 0.016497, 0.602083 ,
         0.227983, 0.016007, 0.604867 ,
         0.234715, 0.015502, 0.607592 ,
         0.241396, 0.014979, 0.610259 ,
         0.248032, 0.014439, 0.612868 ,
         0.254627, 0.013882, 0.615419 ,
         0.261183, 0.013308, 0.617911 ,
         0.267703, 0.012716, 0.620346 ,
         0.274191, 0.012109, 0.622722 ,
         0.280648, 0.011488, 0.625038 ,
         0.287076, 0.010855, 0.627295 ,
         0.293478, 0.010213, 0.629490 ,
         0.299855, 0.009561, 0.631624 ,
         0.306210, 0.008902, 0.633694 ,
         0.312543, 0.008239, 0.635700 ,
         0.318856, 0.007576, 0.637640 ,
         0.325150, 0.006915, 0.639512 ,
         0.331426, 0.006261, 0.641316 ,
         0.337683, 0.005618, 0.643049 ,
         0.343925, 0.004991, 0.644710 ,
         0.350150, 0.004382, 0.646298 ,
         0.356359, 0.003798, 0.647810 ,
         0.362553, 0.003243, 0.649245 ,
         0.368733, 0.002724, 0.650601 ,
         0.374897, 0.002245, 0.651876 ,
         0.381047, 0.001814, 0.653068 ,
         0.387183, 0.001434, 0.654177 ,
         0.393304, 0.001114, 0.655199 ,
         0.399411, 0.000859, 0.656133 ,
         0.405503, 0.000678, 0.656977 ,
         0.411580, 0.000577, 0.657730 ,
         0.417642, 0.000564, 0.658390 ,
         0.423689, 0.000646, 0.658956 ,
         0.429719, 0.000831, 0.659425 ,
         0.435734, 0.001127, 0.659797 ,
         0.441732, 0.001540, 0.660069 ,
         0.447714, 0.002080, 0.660240 ,
         0.453677, 0.002755, 0.660310 ,
         0.459623, 0.003574, 0.660277 ,
         0.465550, 0.004545, 0.660139 ,
         0.471457, 0.005678, 0.659897 ,
         0.477344, 0.006980, 0.659549 ,
         0.483210, 0.008460, 0.659095 ,
         0.489055, 0.010127, 0.658534 ,
         0.494877, 0.011990, 0.657865 ,
         0.500678, 0.014055, 0.657088 ,
         0.506454, 0.016333, 0.656202 ,
         0.512206, 0.018833, 0.655209 ,
         0.517933, 0.021563, 0.654109 ,
         0.523633, 0.024532, 0.652901 ,
         0.529306, 0.027747, 0.651586 ,
         0.534952, 0.031217, 0.650165 ,
         0.540570, 0.034950, 0.648640 ,
         0.546157, 0.038954, 0.647010 ,
         0.551715, 0.043136, 0.645277 ,
         0.557243, 0.047331, 0.643443 ,
         0.562738, 0.051545, 0.641509 ,
         0.568201, 0.055778, 0.639477 ,
         0.573632, 0.060028, 0.637349 ,
         0.579029, 0.064296, 0.635126 ,
         0.584391, 0.068579, 0.632812 ,
         0.589719, 0.072878, 0.630408 ,
         0.595011, 0.077190, 0.627917 ,
         0.600266, 0.081516, 0.625342 ,
         0.605485, 0.085854, 0.622686 ,
         0.610667, 0.090204, 0.619951 ,
         0.615812, 0.094564, 0.617140 ,
         0.620919, 0.098934, 0.614257 ,
         0.625987, 0.103312, 0.611305 ,
         0.631017, 0.107699, 0.608287 ,
         0.636008, 0.112092, 0.605205 ,
         0.640959, 0.116492, 0.602065 ,
         0.645872, 0.120898, 0.598867 ,
         0.650746, 0.125309, 0.595617 ,
         0.655580, 0.129725, 0.592317 ,
         0.660374, 0.134144, 0.588971 ,
         0.665129, 0.138566, 0.585582 ,
         0.669845, 0.142992, 0.582154 ,
         0.674522, 0.147419, 0.578688 ,
         0.679160, 0.151848, 0.575189 ,
         0.683758, 0.156278, 0.571660 ,
         0.688318, 0.160709, 0.568103 ,
         0.692840, 0.165141, 0.564522 ,
         0.697324, 0.169573, 0.560919 ,
         0.701769, 0.174005, 0.557296 ,
         0.706178, 0.178437, 0.553657 ,
         0.710549, 0.182868, 0.550004 ,
         0.714883, 0.187299, 0.546338 ,
         0.719181, 0.191729, 0.542663 ,
         0.723444, 0.196158, 0.538981 ,
         0.727670, 0.200586, 0.535293 ,
         0.731862, 0.205013, 0.531601 ,
         0.736019, 0.209439, 0.527908 ,
         0.740143, 0.213864, 0.524216 ,
         0.744232, 0.218288, 0.520524 ,
         0.748289, 0.222711, 0.516834 ,
         0.752312, 0.227133, 0.513149 ,
         0.756304, 0.231555, 0.509468 ,
         0.760264, 0.235976, 0.505794 ,
         0.764193, 0.240396, 0.502126 ,
         0.768090, 0.244817, 0.498465 ,
         0.771958, 0.249237, 0.494813 ,
         0.775796, 0.253658, 0.491171 ,
         0.779604, 0.258078, 0.487539 ,
         0.783383, 0.262500, 0.483918 ,
         0.787133, 0.266922, 0.480307 ,
         0.790855, 0.271345, 0.476706 ,
         0.794549, 0.275770, 0.473117 ,
         0.798216, 0.280197, 0.469538 ,
         0.801855, 0.284626, 0.465971 ,
         0.805467, 0.289057, 0.462415 ,
         0.809052, 0.293491, 0.458870 ,
         0.812612, 0.297928, 0.455338 ,
         0.816144, 0.302368, 0.451816 ,
         0.819651, 0.306812, 0.448306 ,
         0.823132, 0.311261, 0.444806 ,
         0.826588, 0.315714, 0.441316 ,
         0.830018, 0.320172, 0.437836 ,
         0.833422, 0.324635, 0.434366 ,
         0.836801, 0.329105, 0.430905 ,
         0.840155, 0.333580, 0.427455 ,
         0.843484, 0.338062, 0.424013 ,
         0.846788, 0.342551, 0.420579 ,
         0.850066, 0.347048, 0.417153 ,
         0.853319, 0.351553, 0.413734 ,
         0.856547, 0.356066, 0.410322 ,
         0.859750, 0.360588, 0.406917 ,
         0.862927, 0.365119, 0.403519 ,
         0.866078, 0.369660, 0.400126 ,
         0.869203, 0.374212, 0.396738 ,
         0.872303, 0.378774, 0.393355 ,
         0.875376, 0.383347, 0.389976 ,
         0.878423, 0.387932, 0.386600 ,
         0.881443, 0.392529, 0.383229 ,
         0.884436, 0.397139, 0.379860 ,
         0.887402, 0.401762, 0.376494 ,
         0.890340, 0.406398, 0.373130 ,
         0.893250, 0.411048, 0.369768 ,
         0.896131, 0.415712, 0.366407 ,
         0.898984, 0.420392, 0.363047 ,
         0.901807, 0.425087, 0.359688 ,
         0.904601, 0.429797, 0.356329 ,
         0.907365, 0.434524, 0.352970 ,
         0.910098, 0.439268, 0.349610 ,
         0.912800, 0.444029, 0.346251 ,
         0.915471, 0.448807, 0.342890 ,
         0.918109, 0.453603, 0.339529 ,
         0.920714, 0.458417, 0.336166 ,
         0.923287, 0.463251, 0.332801 ,
         0.925825, 0.468103, 0.329435 ,
         0.928329, 0.472975, 0.326067 ,
         0.930798, 0.477867, 0.322697 ,
         0.933232, 0.482780, 0.319325 ,
         0.935630, 0.487712, 0.315952 ,
         0.937990, 0.492667, 0.312575 ,
         0.940313, 0.497642, 0.309197 ,
         0.942598, 0.502639, 0.305816 ,
         0.944844, 0.507658, 0.302433 ,
         0.947051, 0.512699, 0.299049 ,
         0.949217, 0.517763, 0.295662 ,
         0.951344, 0.522850, 0.292275 ,
         0.953428, 0.527960, 0.288883 ,
         0.955470, 0.533093, 0.285490 ,
         0.957469, 0.538250, 0.282096 ,
         0.959424, 0.543431, 0.278701 ,
         0.961336, 0.548636, 0.275305 ,
         0.963203, 0.553865, 0.271909 ,
         0.965024, 0.559118, 0.268513 ,
         0.966798, 0.564396, 0.265118 ,
         0.968526, 0.569700, 0.261721 ,
         0.970205, 0.575028, 0.258325 ,
         0.971835, 0.580382, 0.254931 ,
         0.973416, 0.585761, 0.251540 ,
         0.974947, 0.591165, 0.248151 ,
         0.976428, 0.596595, 0.244767 ,
         0.977856, 0.602051, 0.241387 ,
         0.979233, 0.607532, 0.238013 ,
         0.980556, 0.613039, 0.234646 ,
         0.981826, 0.618572, 0.231287 ,
         0.983041, 0.624131, 0.227937 ,
         0.984199, 0.629718, 0.224595 ,
         0.985301, 0.635330, 0.221265 ,
         0.986345, 0.640969, 0.217948 ,
         0.987332, 0.646633, 0.214648 ,
         0.988260, 0.652325, 0.211364 ,
         0.989128, 0.658043, 0.208100 ,
         0.989935, 0.663787, 0.204859 ,
         0.990681, 0.669558, 0.201642 ,
         0.991365, 0.675355, 0.198453 ,
         0.991985, 0.681179, 0.195295 ,
         0.992541, 0.687030, 0.192170 ,
         0.993032, 0.692907, 0.189084 ,
         0.993456, 0.698810, 0.186041 ,
         0.993814, 0.704741, 0.183043 ,
         0.994103, 0.710698, 0.180097 ,
         0.994324, 0.716681, 0.177208 ,
         0.994474, 0.722691, 0.174381 ,
         0.994553, 0.728728, 0.171622 ,
         0.994561, 0.734791, 0.168938 ,
         0.994495, 0.740880, 0.166335 ,
         0.994355, 0.746995, 0.163821 ,
         0.994141, 0.753137, 0.161404 ,
         0.993851, 0.759304, 0.159092 ,
         0.993482, 0.765499, 0.156891 ,
         0.993033, 0.771720, 0.154808 ,
         0.992505, 0.777967, 0.152855 ,
         0.991897, 0.784239, 0.151042 ,
         0.991209, 0.790537, 0.149377 ,
         0.990439, 0.796859, 0.147870 ,
         0.989587, 0.803205, 0.146529 ,
         0.988648, 0.809579, 0.145357 ,
         0.987621, 0.815978, 0.144363 ,
         0.986509, 0.822401, 0.143557 ,
         0.985314, 0.828846, 0.142945 ,
         0.984031, 0.835315, 0.142528 ,
         0.982653, 0.841812, 0.142303 ,
         0.981190, 0.848329, 0.142279 ,
         0.979644, 0.854866, 0.142453 ,
         0.977995, 0.861432, 0.142808 ,
         0.976265, 0.868016, 0.143351 ,
         0.974443, 0.874622, 0.144061 ,
         0.972530, 0.881250, 0.144923 ,
         0.970533, 0.887896, 0.145919 ,
         0.968443, 0.894564, 0.147014 ,
         0.966271, 0.901249, 0.148180 ,
         0.964021, 0.907950, 0.149370 ,
         0.961681, 0.914672, 0.150520 ,
         0.959276, 0.921407, 0.151566 ,
         0.956808, 0.928152, 0.152409 ,
         0.954287, 0.934908, 0.152921 ,
         0.951726, 0.941671, 0.152925 ,
         0.949151, 0.948435, 0.152178 ,
         0.946602, 0.955190, 0.150328 ,
         0.944152, 0.961916, 0.146861 ,
         0.941896, 0.968590, 0.140956 ,
         0.940015, 0.975158, 0.131326).finished().transpose();
    return plasma.col(static_cast<int>(std::round( std::max<float>(0.f,std::min<float>(s,1.f)) * 255.0)));
}

Eigen::Vector3f VisMesh::colorFromNormal( const Eigen::Vector3f & normal )
{
    const float ngfx = (normal.x()+1)*.5;
    const float ngfy = (normal.y()+1)*.5;
    const float ngfz = (normal.z()+1)*.5;
    return Eigen::Vector3f(ngfx,ngfy,ngfz);
}

size_t VisMesh::size() const
{
    return vertices.size();
}

void VisMesh::addPoint( const Eigen::Vector3f & vertex_one,
                        const Eigen::Vector3f & color,
                        const float * intensity,
                        const uint16_t * reflectivity
                      )
{
    if ( !vertex_one.allFinite() ) return;
    vertices.emplace_back(vertex_one); colors.emplace_back(color);
    if ( intensity ) intensities.emplace_back(*intensity);
    if ( reflectivity ) reflectivities.emplace_back(Eigen::Vector3f::Constant(*reflectivity/65535.));
}

CloudPtr VisMesh::getPointMesh( ) const
{
    CloudPtr mesh = Cloud::create();
    if ( vertices.empty() ) { mesh->V.resize(vertices.size(),3); return mesh; }
    mesh->V = vec2eigen(vertices).cast<double>();
    mesh->C = vec2eigen(colors).cast<double>();
    if ( intensities.size() == vertices.size() ) mesh->I = vec2eigen(intensities).cast<double>() * 3;
    //if ( reflectivities.size() == vertices.size() ) mesh->C = vec2eigen(reflectivities).cast<double>() * 1.5;
    mesh->m_vis.m_show_points=true;
    mesh->m_vis.m_color_type=MeshColorType::PerVertColor;
    return mesh;
}

CloudPtr VisMesh::showPointMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf, const int & point_size ) const
{
    CloudPtr mesh = getPointMesh();
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    mesh->m_vis.m_point_size = std::max(point_size,1);
    Scene::show(mesh, meshName);
    return mesh;
}

CloudPtr VisMesh::getCubeMesh( const float & cube_size ) const
{
    CloudPtr mesh = Cloud::create();
    if ( vertices.empty() ) { mesh->V.resize(vertices.size(),3); return mesh; }
    if (vertices.empty()) return mesh;

    mesh->C.resize(6*4 * vertices.size(),3);
    mesh->NV.resize(6*4 * vertices.size(),3);
    mesh->F.resize(6*2 * vertices.size(),3);
    mesh->NF.resize(6*2 * vertices.size(),3);

    static const Eigen::Vector3f origins[6] =
    {
        Eigen::Vector3f(-1.0, -1.0, -1.0),
        Eigen::Vector3f(1.0, -1.0, -1.0),
        Eigen::Vector3f(1.0, -1.0, 1.0),
        Eigen::Vector3f(-1.0, -1.0, 1.0),
        Eigen::Vector3f(-1.0, 1.0, -1.0),
        Eigen::Vector3f(-1.0, -1.0, 1.0)
    };
    static const Eigen::Vector3f rights[6] =
    {
        Eigen::Vector3f(2.0, 0.0, 0.0),
        Eigen::Vector3f(0.0, 0.0, 2.0),
        Eigen::Vector3f(-2.0, 0.0, 0.0),
        Eigen::Vector3f(0.0, 0.0, -2.0),
        Eigen::Vector3f(2.0, 0.0, 0.0),
        Eigen::Vector3f(2.0, 0.0, 0.0)
    };
    static const Eigen::Vector3f ups[6] =
    {
        Eigen::Vector3f(0.0, 2.0, 0.0),
        Eigen::Vector3f(0.0, 2.0, 0.0),
        Eigen::Vector3f(0.0, 2.0, 0.0),
        Eigen::Vector3f(0.0, 2.0, 0.0),
        Eigen::Vector3f(0.0, 0.0, 2.0),
        Eigen::Vector3f(0.0, 0.0, -2.0)
    };
    Eigen::Matrix3Xf normals = Eigen::Matrix3Xf(3,6);
    Eigen::Matrix3Xf colors = Eigen::Matrix3Xf(3,6);
    std::vector<Eigen::Matrix3Xf> faces;
    for ( int faceIdx = 0; faceIdx < 6; ++faceIdx )
    {
        faces.emplace_back(Eigen::Matrix3Xf(3,4));
        Eigen::Matrix3Xf & pts = faces[faceIdx];
        pts.col(0) = origins[faceIdx];
        pts.col(1) = origins[faceIdx] + rights[faceIdx];
        pts.col(2) = origins[faceIdx] + rights[faceIdx] + ups[faceIdx];
        pts.col(3) = origins[faceIdx] + ups[faceIdx];
        normals.col(faceIdx) = rights[faceIdx].cross(ups[faceIdx]);
        colors.col(faceIdx) = colorFromNormal(normals.col(faceIdx));
    }

    for ( size_t i = 0; i < (vertices.size()); ++i )
    {
        for ( int faceIdx = 0; faceIdx < 6; ++faceIdx )
        {
            const Eigen::Matrix3Xf & pts = faces[faceIdx];
            const int pointOffset = i*6*4+faceIdx*4;
            const int faceOffset = i*6*2+faceIdx*2;

            mesh->F.row(faceOffset) = Eigen::Vector3i(pointOffset, pointOffset+1, pointOffset+2).transpose();
            mesh->F.row(faceOffset+1) = Eigen::Vector3i(pointOffset+2, pointOffset+3, pointOffset).transpose();
            mesh->NF.row(faceOffset) = normals.col(faceIdx).transpose().cast<double>();
            mesh->NF.row(faceOffset+1) = normals.col(faceIdx).transpose().cast<double>();
            for ( int vertIdx = 0; vertIdx < 4; ++vertIdx)
            {
                mesh->V.row(pointOffset + vertIdx) = ((pts.col(vertIdx))*(cube_size/2) + vertices[i]).transpose().cast<double>();
                mesh->NV.row(pointOffset + vertIdx) = normals.col(faceIdx).transpose().cast<double>();
                mesh->C.row(pointOffset + vertIdx) = colors.col(faceIdx).transpose().cast<double>();
            }
        }
    }
    //LOG(INFO) << "F:\n" << mesh->F;
    LOG(INFO) << "occAll: F: " << mesh->F.rows() << " " << mesh->F.maxCoeff() << " V: " << mesh->V.rows() << " v: " << vertices.size();

//    mesh->V = vec2eigen(vertices);
//    mesh->C = vec2eigen(colors);
//    if ( intensities.size() == vertices.size() ) mesh->I = vec2eigen(intensities) * 3;
//    //if ( reflectivities.size() == vertices.size() ) mesh->C = vec2eigen(reflectivities) * 1.5;
    mesh->m_vis.m_show_points=true;
    mesh->m_vis.m_show_mesh=true;
    mesh->m_vis.m_color_type=MeshColorType::PerVertColor;
    return mesh;
}

CloudPtr VisMesh::showCubeMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf, const float & cube_size ) const
{
    CloudPtr mesh = getCubeMesh( cube_size );
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    Scene::show(mesh, meshName);
    return mesh;
}


CloudPtr VisMesh::getEdgeMeshFromPoints( ) const
{
    CloudPtr mesh = Cloud::create();
    if ( vertices.empty() ) { mesh->V.resize(vertices.size(),3); return mesh; }
    mesh->V = vec2eigen(vertices).cast<double>();
    mesh->C = vec2eigen(colors).cast<double>();
    mesh->E = Eigen::MatrixXi(std::max<int>(0,mesh->V.rows()-1), 2);
    for ( int row = 0; row < mesh->V.rows()-1; ++row )
    {
        mesh->E.row(row) << row, row+1;
    }
    mesh->m_vis.m_show_points=true;
    mesh->m_vis.m_show_lines=true;
    mesh->m_vis.m_color_type=MeshColorType::PerVertColor;
    return mesh;
}

CloudPtr VisMesh::showEdgeMeshFromPoints( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf ) const
{
    CloudPtr mesh = getEdgeMeshFromPoints();
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    Scene::show(mesh, meshName);
    return mesh;
}


void VisMesh::addEdge( const Eigen::Vector3f & vertex_one,
                       const Eigen::Vector3f & vertex_two,
                       const Eigen::Vector3f & color
                     )
{
    addEdge(vertex_one,vertex_two,color,color);
}
void VisMesh::addEdge( const Eigen::Vector3f & vertex_one,
                       const Eigen::Vector3f & vertex_two,
                       const Eigen::Vector3f & color_one,
                       const Eigen::Vector3f & color_two
                     )
{
    if ( !vertex_one.allFinite() || ! vertex_two.allFinite() ) return;
    const int firstIdx = vertices.size();
    const int secondIdx = firstIdx+1;
    vertices.emplace_back(vertex_one); colors.emplace_back(color_one);
    vertices.emplace_back(vertex_two); colors.emplace_back(color_two);
    edges.emplace_back((Eigen::Vector2i()<<firstIdx,secondIdx).finished());
}

CloudPtr VisMesh::getEdgeMesh( ) const
{
    CloudPtr mesh = Cloud::create();
    if ( vertices.empty() ) { mesh->V.resize(vertices.size(),3); return mesh; }
    mesh->V = vec2eigen(vertices).cast<double>();
    mesh->C = vec2eigen(colors).cast<double>();
    mesh->E = vec2eigen(edges);
    mesh->m_vis.m_show_points=true;
    mesh->m_vis.m_show_lines=true;
    mesh->m_vis.m_color_type=MeshColorType::PerVertColor;
    return mesh;
}
CloudPtr VisMesh::showEdgeMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf) const
{
    CloudPtr mesh = getEdgeMesh();
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    //mesh->m_vis.m_point_size = 24;
    //mesh->m_vis.m_show_lines = false;
    Scene::show(mesh, meshName);
    return mesh;
}

void VisMesh::addNormal( const Eigen::Vector3f & mean,
                         const Eigen::Vector3f & normal,
                         const Eigen::Vector3f & color
                       )
{
    addEdge( mean, mean+0.1*normal, color);
}

void VisMesh::addSurfel( const Eigen::MatrixXf & evecs,
                         const Eigen::VectorXf & evals,
                         const Eigen::Vector3f & mean,
                         const Eigen::Vector3f & color,
                         const float & sigma,
                         const SurfelType & type )
{
    switch(type)
    {
        default:
        case SurfelType::TriangleEllipsoid: addSurfelTriangleEllipsoid( evecs, evals, mean, color, sigma ); break;
        case SurfelType::Quad: addSurfelQuad( evecs, evals, mean, color ); break;
        case SurfelType::Ellipse: addSurfelEllipsoid( evecs, evals, mean, color); break;
    }
}

void VisMesh::addSurfelQuad( const Eigen::MatrixXf & evecs,
                             const Eigen::VectorXf & evals,
                             const Eigen::Vector3f & mean,
                             const Eigen::Vector3f & color,
                             const float & sigma
                           )
{
    const Eigen::Vector3f normal = evecs.col(0).cast<float>().normalized();
    const float r_v = sigma*sqrt(evals.y());
    const float r_u = sigma*sqrt(evals.z());
    const Eigen::Vector3f u = evecs.col(2).cast<float>().normalized(); // eigenvector to biggest eigenvalue
    const Eigen::Vector3f v = normal.normalized().cross(u.normalized()).normalized();
    const int firstIdx = vertices.size();
    const int secondIdx = firstIdx+1;
    const int thirdIdx  = firstIdx+2;
    const int fourthIdx = firstIdx+3;
    vertices.emplace_back((mean + u*r_u + v*r_v)); vertex_normals.emplace_back(normal); colors.emplace_back(color);
    vertices.emplace_back((mean + v*r_v - u*r_u)); vertex_normals.emplace_back(normal); colors.emplace_back(color);
    vertices.emplace_back((mean - v*r_v + u*r_u)); vertex_normals.emplace_back(normal); colors.emplace_back(color);
    vertices.emplace_back((mean - u*r_u - v*r_v)); vertex_normals.emplace_back(normal); colors.emplace_back(color);
    faces.emplace_back((Eigen::VectorXi(3,1)<< firstIdx, secondIdx, thirdIdx).finished()); face_normals.emplace_back(normal);
    faces.emplace_back((Eigen::VectorXi(3,1)<< secondIdx, fourthIdx, thirdIdx).finished()); face_normals.emplace_back(normal);
}

CloudPtr VisMesh::getSurfelMesh( ) const
{
    CloudPtr mesh = Cloud::create();
    if ( vertices.empty() ) { mesh->V.resize(vertices.size(),3); return mesh; }
    mesh->V = vec2eigen(vertices).cast<double>();
    mesh->C = vec2eigen(colors).cast<double>();
    mesh->F = vec2eigen(faces);
    mesh->NV = vec2eigen(vertex_normals).cast<double>();
    mesh->NF = vec2eigen(face_normals).cast<double>();
    mesh->m_vis.m_show_points=false;
    mesh->m_vis.m_show_mesh=true;
    mesh->m_vis.m_color_type=MeshColorType::PerVertColor;
    return mesh;
}

CloudPtr VisMesh::showSurfelQuadMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf ) const
{
    CloudPtr mesh = getSurfelMesh();
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    Scene::show(mesh, meshName);
    return mesh;
}

void VisMesh::addSurfelEllipsoid(  const Eigen::MatrixXf & evecs,
                                   const Eigen::VectorXf & evals,
                                   const Eigen::Vector3f & mean,
                                   const Eigen::Vector3f & color,
                                   const float & sigma
                                 )
{
    const Eigen::Vector3f normal = evecs.col(0).cast<float>().normalized();
    const Eigen::Vector3f u = evecs.col(2).cast<float>().normalized(); // eigenvector to biggest eigenvalue
    const float r_w = sigma*sqrt(evals.x());
    const float r_v = sigma*sqrt(evals.y());
    const float r_u = sigma*sqrt(evals.z());
    tangent_direction.emplace_back(r_u * u);
    tangent_other_length.emplace_back(r_v);
    tangent_normal_length.emplace_back(r_w);
    vertices.emplace_back(mean); vertex_normals.emplace_back(normal); colors.emplace_back(color);
}

CloudPtr VisMesh::getSurfelEllipsoidMesh( ) const
{
    CloudPtr mesh = Cloud::create();
    if ( vertices.empty() ) { mesh->V.resize(vertices.size(),3); return mesh; }
    mesh->V = vec2eigen(vertices).cast<double>();
    mesh->C = vec2eigen(colors).cast<double>();
    mesh->NV = vec2eigen(vertex_normals).cast<double>();
    //for surfel rendering each vertex has a 2 vectors that are tangent defining the span of the elipsoid. For memory usage we don't store the 2 vectors directly because we alreayd have a normal vector, rather we store one tangent vector in full (vec3) and the other one we store only the norm of it because it's direction can be inferred as the cross product between the normal and the first tangent vector
    mesh->V_tangent_u  = vec2eigen(tangent_direction).cast<double>();
    mesh->V_length_v = vec2eigen(tangent_other_length).cast<double>();
    mesh->m_vis.m_show_points=false;
    mesh->m_vis.m_show_mesh=false;
    mesh->m_vis.m_show_surfels = true;
    mesh->m_vis.m_color_type=MeshColorType::PerVertColor;
    return mesh;
}

CloudPtr VisMesh::showSurfelEllipsoidMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf ) const
{
    CloudPtr mesh = getSurfelEllipsoidMesh();
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    Scene::show(mesh, meshName);
    return mesh;
}

void VisMesh::reserveSurfelTriangleEllipsoid ( const int num_surfels )
{
    const int num_faces = iso.triangles.size();
    const int num_verts = iso.vertices.size();
    faces.reserve( num_surfels * num_faces);
    face_normals.reserve(num_surfels * num_faces);
    vertices.reserve(num_surfels * num_verts);
    vertex_normals.reserve(num_surfels * num_verts);
    colors.reserve(num_surfels * num_verts);
}

void VisMesh::addSurfelTriangleEllipsoid(  const Eigen::MatrixXf & evecs,
                                           const Eigen::VectorXf & evals,
                                           const Eigen::Vector3f & mean,
                                           const Eigen::Vector3f & color,
                                           const float & sigma
                                        )
{

    //while ( iso.subdivided < 1 ){ iso.subdivide(); };
    static Eigen::Matrix3Xi indices;
    if ( indices.cols() <= 0 )
        iso.getTriangles ( indices );

    const Eigen::Vector3f normal = evecs.col(0).cast<float>().normalized();
    const Eigen::Vector3f u = evecs.col(2).cast<float>().normalized(); // eigenvector to biggest eigenvalue
    const Eigen::Vector3f v = normal.normalized().cross(u.normalized()).normalized();
    const float r_w = sigma*sqrt(evals.x());
    const float r_v = sigma*sqrt(evals.y());
    const float r_u = sigma*sqrt(evals.z());
    Eigen::Matrix3f evec = Eigen::Matrix3f::Identity();
    evec.col(0) = normal; evec.col(1) = v; evec.col(2) = u;
    const Eigen::Vector3f scaled_evals ( r_w, r_v, r_u );
    Eigen::Matrix3Xf verts;
    iso.getEllipsoid ( evec, scaled_evals, mean, verts );

    const int firstIdx = vertices.size();
    const Eigen::Vector3i firstIdxVec = Eigen::Vector3i::Constant(firstIdx);
    const int num_indices = indices.cols();
    for ( int idx = 0; idx < num_indices; ++idx )
    {
        const Eigen::Vector3i & face = indices.col(idx);
        faces.emplace_back(face+firstIdxVec); face_normals.emplace_back(normal);
    }
    const int num_verts = verts.cols();
    for ( int idx = 0; idx < num_verts; ++idx )
    {
        const Eigen::Vector3f & vert = verts.col(idx);
        vertices.emplace_back(vert); vertex_normals.emplace_back(normal); colors.emplace_back(color);
    }
}

CloudPtr VisMesh::getSurfelTriangleEllipsoidMesh( ) const
{
    return getSurfelMesh();
}

CloudPtr VisMesh::showSurfelTriangleEllipsoidMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf ) const
{
    CloudPtr mesh = getSurfelTriangleEllipsoidMesh();
    if ( ! tf.params().isApprox(Sophus::SE3d().params()) )
        mesh->transform_vertices_cpu(Eigen::Affine3d(tf.matrix()));
    if ( rotate90deg )
        mesh->worldROS2worldGL();
    Scene::show(mesh, meshName);
    return mesh;
}
void VisMesh::showAggregatedSurfelTriangleEllipsoidMesh ( const Sophus::SE3d & oldLastScenePose, const std::string & meshName )
{
    CloudPtr new_mesh = getSurfelTriangleEllipsoidMesh();
    if ( ! oldLastScenePose.params().isApprox(Sophus::SE3d().params()) )
        new_mesh->transform_vertices_cpu(Eigen::Affine3d(oldLastScenePose.matrix()));
    new_mesh->worldROS2worldGL();
    if ( ! prev_mesh ) prev_mesh = new_mesh;
    else prev_mesh->add(*new_mesh);
    Scene::show(prev_mesh, meshName);
}

CloudPtr VisMesh::showSurfelMesh( const std::string & meshName, const bool & rotate90deg, const Sophus::SE3d & tf, const SurfelType & type ) const
{
    switch( type )
    {
        default: case SurfelType::TriangleEllipsoid: return showSurfelTriangleEllipsoidMesh( meshName, rotate90deg, tf ); break;
        case SurfelType::Quad: return showSurfelQuadMesh( meshName, rotate90deg, tf ); break;
        case SurfelType::Ellipse: return showSurfelEllipsoidMesh( meshName, rotate90deg, tf ); break;
    }
}


void VisMesh::resetForAggregation()
{
    vertices.clear();
    colors.clear();
    vertex_normals.clear();
    tangent_direction.clear();
    tangent_other_length.clear();
    tangent_normal_length.clear();
    edges.clear();
    faces.clear();
    face_normals.clear();
    intensities.clear();
    reflectivities.clear();
    // dont clear prev_mesh
}

void VisMesh::showAggregatedPointMesh ( const Sophus::SE3d & oldLastScenePose, const std::string & meshName, const float & randomSampleFactor, const bool & show )
{
    CloudPtr new_mesh = getPointMesh();
    //LOG(INFO) << "NewMesh: " << new_mesh->V.rows() << " rsf: " << randomSampleFactor;
    if ( randomSampleFactor > 0. && randomSampleFactor < 1. )
        new_mesh->random_subsample(randomSampleFactor);
    //LOG(INFO) << "NewMesh after RSS: " << new_mesh->V.rows() << " rsf: " << randomSampleFactor;
    if ( ! oldLastScenePose.params().isApprox(Sophus::SE3d().params()) )
        new_mesh->transform_vertices_cpu(Eigen::Affine3d(oldLastScenePose.matrix()));
    new_mesh->worldROS2worldGL();
    if ( ! prev_mesh ) prev_mesh = new_mesh;
    else prev_mesh->add(*new_mesh);
    //prev_mesh->m_vis.m_color_type = MeshColorType::Height;
    //prev_mesh->m_min_max_y_for_plotting = Eigen::Vector2f(-2,2);
    if ( show && prev_mesh )
        Scene::show(prev_mesh, meshName);
}

void VisMesh::showAggregatedEdgeMesh ( const Sophus::SE3d & oldLastScenePose, const std::string & meshName, const bool & show )
{
    CloudPtr new_mesh = getEdgeMesh();
    if ( ! oldLastScenePose.params().isApprox(Sophus::SE3d().params()) )
        new_mesh->transform_vertices_cpu(Eigen::Affine3d(oldLastScenePose.matrix()));
    new_mesh->worldROS2worldGL();
    if ( ! prev_mesh ) prev_mesh = new_mesh;
    else prev_mesh->add(*new_mesh);
    if ( show )
        Scene::show(prev_mesh, meshName);
}

void VisMesh::reservePoints( const int num_points )
{
    vertices.reserve(num_points);
    colors.reserve(num_points);
    intensities.reserve(num_points);
    reflectivities.reserve(num_points);
    vertex_normals.reserve(num_points);
}
#endif
