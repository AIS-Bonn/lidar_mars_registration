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
#include "EasyPBRwrapper.h"

#include <fstream>
#include <chrono>
#include <thread>

//loguru
//#define LOGURU_NO_DATE_TIME 1
//#define LOGURU_NO_UPTIME 1
//#define LOGURU_REPLACE_GLOG 1
#include <loguru.hpp>

//My stuff 
#include "UtilsGL.h"
#include "easy_pbr/Viewer.h"
#include "easy_pbr/Gui.h"
#include "easy_pbr/MeshGL.h"
#include "easy_pbr/Mesh.h"
#include "easy_pbr/Scene.h"
#include "RandGenerator.h"

//Add this header after we add all opengl stuff because we need the profiler to have glFinished defined
//#define ENABLE_GL_PROFILING 1
//#include "Profiler.h"

//#include "imgui_plot.h"

//configuru
#define CONFIGURU_WITH_EIGEN 1
#define CONFIGURU_IMPLICIT_CONVERSIONS 1
#include <configuru.hpp>
using namespace configuru;

//boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "BagAdaptor.h"
#include "MarsSplineRegistrationAdaptor.h"

using namespace easy_pbr;

// SyntheticGenerator::SyntheticGenerator(const std::string& config_file):
EasyPBRwrapper::EasyPBRwrapper(const std::string& config_file, const std::shared_ptr<Viewer>& view):
    #ifdef WITH_DIR_WATCHER
    dir_watcher(std::string(PROJECT_SOURCE_DIR)+"/shaders/",5),
    #endif
    m_view(view),
    m_fullscreen_quad(MeshGL::create()),
    m_iter(0),
    m_start(false)
{

    init_params(config_file);
    compile_shaders();
    init_opengl();
    install_callbacks(view);

    //creates a mesh and adds it to the scene
    //create_mesh();
}

EasyPBRwrapper::~EasyPBRwrapper(){

}

void EasyPBRwrapper::init_params(const std::string config_file){
    //read all the parameters
    std::string config_file_abs;
    if (fs::path(config_file).is_relative()){
        config_file_abs=(fs::path(PROJECT_SOURCE_DIR) / config_file).string();
    }else{
        config_file_abs=config_file;
    }
    Config cfg = configuru::parse_file(config_file_abs, CFG);
    Config wrapper_config=cfg["wrapper"];

    m_init_auto_start = wrapper_config["auto_start"];
    m_init_play_num = wrapper_config["play_num"];
    m_init_skip_num = wrapper_config["skip_num"];

    unsigned int seed=0;
    m_rand_gen=std::make_shared<radu::utils::RandGenerator> (seed);

}


void EasyPBRwrapper::compile_shaders(){

    // m_detect_balloon_shader.compile( std::string(PROJECT_SOURCE_DIR)+"/shaders/detect_balloon_vert.glsl", std::string(PROJECT_SOURCE_DIR)+"/shaders/detect_balloon_frag.glsl" ) ;
    // m_detect_copter_shader.compile( std::string(PROJECT_SOURCE_DIR)+"/shaders/detect_copter_vert.glsl", std::string(PROJECT_SOURCE_DIR)+"/shaders/detect_copter_frag.glsl" ) ;
}

void EasyPBRwrapper::init_opengl(){
    //create a fullscreen quad which we will use for rendering full screen images
    m_fullscreen_quad->m_core->create_full_screen_quad();
    m_fullscreen_quad->upload_to_gpu();
}

void EasyPBRwrapper::hotload_shaders(){
#ifdef WITH_DIR_WATCHER
    std::vector<std::string> changed_files=dir_watcher.poll_files();
    if(changed_files.size()>0){
        compile_shaders();
    }
#endif
}

void EasyPBRwrapper::install_callbacks(const std::shared_ptr<Viewer>& view){
    //pre draw functions (can install multiple functions and they will be called in order)
    //view->add_callback_pre_draw( [this]( Viewer& v ) -> void{ this->pre_draw_animate_mesh(v); }  );
    //view->add_callback_pre_draw( [this]( Viewer& v ) -> void{ this->pre_draw_colorize_mesh(v); }  );

    //post draw functions
    view->add_callback_post_draw( [this]( Viewer& v ) -> void{ this->post_draw(v); }  );
}

void EasyPBRwrapper::create_mesh(){
    MeshSharedPtr mesh= Mesh::create();
    mesh->V.resize(4,3);
    mesh->V.row(0)<< -1,-1,0;
    mesh->V.row(1)<< 1,-1,0;
    mesh->V.row(2)<< -1,1,0;
    mesh->V.row(3)<< 1,1,0;

    mesh->F.resize(2,3);
    mesh->F.row(0)<< 0,1,2;
    mesh->F.row(1)<< 2,1,3;

    mesh->m_vis.m_show_wireframe=true;
    mesh->m_vis.m_line_width=5.0;
    Scene::show(mesh, "mesh");
}


void EasyPBRwrapper::pre_draw_animate_mesh(Viewer& view){
    //get a handle to the mesh from the scene and move it in the x direction
    MeshSharedPtr mesh= Scene::get_mesh_with_name("mesh");
    mesh->model_matrix().translation().x() = std::sin( m_iter/100.0 );

}

void EasyPBRwrapper::pre_draw_colorize_mesh(Viewer& view){
    //get a handle to the mesh from the scene and modify its color
    MeshSharedPtr mesh= Scene::get_mesh_with_name("mesh");
    mesh->m_vis.m_solid_color.x() = std::fabs(std::sin( m_iter/50.0 ));

}

void EasyPBRwrapper::post_draw(Viewer& view){
    //get the final render as a opencv mat
    //cv::Mat mat = view.rendered_tex_no_gui().download_to_cv_mat();

    //the opencv mat can now be written to disk or even rendered in the GUI as a texture
    //cv::flip(mat, mat, 0);
    //cv::cvtColor(mat, mat, cv::COLOR_BGR2RGB);
    //Gui::show(mat, "mat");

    static bool auto_start = m_init_auto_start;
    static int play_num = m_init_play_num;
    static int skip_num = m_init_skip_num;
    static int last_num = -1;

    if ( m_init_auto_start )
    {
        if ( skip_num > 0 )
        {
            m_skip_num = skip_num;
            m_skip = true;
        }
        if ( play_num > 0 )
        {
            m_play_num = play_num;
            m_play = true;
        }
        m_init_auto_start = false;
    }

    m_should_skip = false;
    m_should_play = false;

    ImGuiWindowFlags window_flags = 0;
    ImGui::Begin("Menu", nullptr, window_flags);
    ImGui::Separator();
    if ( ImGui::Button("Start") )
    {
        m_start = true;
    }
    if (ImGui::CollapsingHeader("BagAdaptor", ImGuiTreeNodeFlags_::ImGuiTreeNodeFlags_DefaultOpen)) {
        const int cur_num = m_bag_adaptor->get_current_index();
        const int total_num = m_bag_adaptor->get_total_num();
        const int num_diff = cur_num-last_num;
        const bool was_keyframe = m_registration_adaptor->was_keyframe();
        if ( m_skip )
        {
            if ( m_skip_num > 0 )
            {
                if ( num_diff > 0 ) m_skip_num -= num_diff;
                if ( m_skip_num > 0 ) m_should_skip = true;
                else m_skip = false;
            }
        }
        else
        {
            if ( m_play_till_keyframe && was_keyframe && play_num != m_play_num )
            {
                m_play_till_keyframe = false;
                m_play = false;
                m_skip = false;
                m_play_num = 0;
                m_skip_num = 0;
            }
            if ( m_play )
            {
                if ( num_diff > 0 ) m_play_num -= num_diff;
                if ( m_play_num > 0 ) m_should_play = true;
                else m_play = false;
            }
        }
        last_num = cur_num;

        ImGui::Checkbox("AutoStart", &auto_start);
        if ( ImGui::Button("skip:") )
        {
            m_skip_num = skip_num;
            m_skip = true;
        }
        ImGui::SameLine();
        ImGui::SliderInt("skipN", &skip_num, 0,1000);
        if ( ImGui::Button("play:") )
        {
            m_play_num = play_num;
            m_play = true;
        }
        ImGui::SameLine();
        ImGui::SliderInt("playN", &play_num, 1,1000);
        if ( ImGui::Button("play next") )
        {
            m_play_num = 1;
            m_play = true;
        }
        ImGui::SameLine();
        if (ImGui::Button ("play2kf") )
        {
            m_play_num = play_num;
            m_play = true;
            m_play_till_keyframe = true;
        }
        if ( ImGui::Button("stop") )
        {
            m_play_num = 0;
            m_play = false;
        }

        ImGui::Text("skipping: %d",m_skip_num);
        //ImGui::Text(("playing: " + std::to_string(m_play_num)).c_str());
        ImGui::Text("playing: %d",m_play_num);
        //ImGui::Text(( std::to_string(cur_num) + " / " + std::to_string(total_num) ).c_str());
        ImGui::Text("%d / %d", cur_num, total_num);
    }
    ImGui::End();
    //m_iter++;
}

void EasyPBRwrapper::register_module(const std::shared_ptr<BagAdaptor> bagAdaptor){
    m_bag_adaptor = bagAdaptor;
}
void EasyPBRwrapper::register_module(const std::shared_ptr<MarsSplineRegistrationAdaptor> regAdaptor){
    m_registration_adaptor = regAdaptor;
}
