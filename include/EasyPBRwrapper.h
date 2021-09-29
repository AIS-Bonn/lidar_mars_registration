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

//c++
#include <iostream>
#include <memory>
#include <set>

// #include <glad/glad.h> // Initialize with gladLoadGL()
// Include glfw3.h after our OpenGL definitions
// #include <GLFW/glfw3.h>

//opencv
// #include "opencv2/opencv.hpp"

//dir watcher
#ifdef WITH_DIR_WATCHER
    #include "dir_watcher/dir_watcher.hpp"
#endif

//gl
#include "Texture2D.h"
#include "Shader.h"


//forward declarations
namespace radu { namespace utils {
    class RandGenerator;
    }
}
namespace easy_pbr{
    class Mesh;
    class MeshGL;
    class Viewer;
    typedef std::shared_ptr<Mesh> MeshSharedPtr;
}
class BagAdaptor;
class MarsSplineRegistrationAdaptor;

//in order to dissalow building on the stack and having only ptrs https://stackoverflow.com/a/17135547
class EasyPBRwrapper;

class EasyPBRwrapper: public std::enable_shared_from_this<EasyPBRwrapper>
{
public:
    //https://stackoverflow.com/questions/29881107/creating-objects-only-as-shared-pointers-through-a-base-class-create-method
    template <class ...Args>
    static std::shared_ptr<EasyPBRwrapper> create( Args&& ...args ){
        return std::shared_ptr<EasyPBRwrapper>( new EasyPBRwrapper(std::forward<Args>(args)...) );
    }
    ~EasyPBRwrapper();

    bool hasStarted() const{return m_start;}
    bool isRunning() const{return m_skip || m_play;}
    bool shouldSkipOne() const{return m_should_skip; }
    bool shouldPlayOne() const{return m_should_play; }
    void register_module(const std::shared_ptr<BagAdaptor> bagAdaptor);
    void register_module(const std::shared_ptr<MarsSplineRegistrationAdaptor> regAdaptor);

private:
    EasyPBRwrapper(const std::string& config_file, const std::shared_ptr<easy_pbr::Viewer>& view); // we put the constructor as private so as to dissalow creating the object on the stack because we want to only used shared ptr for it

    std::shared_ptr<easy_pbr::Viewer> m_view;
    std::shared_ptr<radu::utils::RandGenerator> m_rand_gen;

    std::shared_ptr<easy_pbr::MeshGL> m_fullscreen_quad; //we store it here because we precompute it and then we use for composing the final image after the deffered geom pass
    int m_iter = 0;

    bool m_start = false;
    bool m_play = false;
    bool m_skip = false;
    bool m_play_till_keyframe = false;
    bool m_should_skip = false;
    bool m_should_play = false;
    int m_skip_num = 0;
    int m_play_num = 0;

    bool m_init_auto_start = false;
    int m_init_play_num = 0;
    int m_init_skip_num = 0;

    std::shared_ptr<BagAdaptor> m_bag_adaptor = nullptr;
    std::shared_ptr<MarsSplineRegistrationAdaptor> m_registration_adaptor = nullptr;
    //params

    void init_params(const std::string config_file);
    void compile_shaders();
    void init_opengl();
    void hotload_shaders();
    void install_callbacks(const std::shared_ptr<easy_pbr::Viewer>& view); //installs some callbacks that will be called by the viewer after it finishes an update

    void create_mesh();

    //pre draw callbacks
    void pre_draw_animate_mesh(easy_pbr::Viewer& view);
    void pre_draw_colorize_mesh(easy_pbr::Viewer& view);

    //post draw callbacks
    void post_draw(easy_pbr::Viewer& view);

};
