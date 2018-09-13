// Copyright (c) 2009, V. Lepetit, EPFL
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies, 
//   either expressed or implied, of the FreeBSD Project.

#include <iostream>
#include <fstream>
using namespace std;

#include "epnp.h"

const float uc = 320;
const float vc = 240;
const float fu = 800;
const float fv = 800;


// MtM takes more time than 12x12 opencv SVD with about 180 points and more:

const int n = 10;
const float noise = 10;

float rand(float min, float max)
{
    return min + (max - min) * float(rand()) / RAND_MAX;
}

void random_pose(float R[3][3], float t[3])
{
    const float range = 1;

    float phi   = rand(0, range * 3.14159 * 2);
    float theta = rand(0, range * 3.14159);
    float psi   = rand(0, range * 3.14159 * 2);

    R[0][0] = cos(psi) * cos(phi) - cos(theta) * sin(phi) * sin(psi);
    R[0][1] = cos(psi) * sin(phi) + cos(theta) * cos(phi) * sin(psi);
    R[0][2] = sin(psi) * sin(theta);

    R[1][0] = -sin(psi) * cos(phi) - cos(theta) * sin(phi) * cos(psi);
    R[1][1] = -sin(psi) * sin(phi) + cos(theta) * cos(phi) * cos(psi);
    R[1][2] = cos(psi) * sin(theta);

    R[2][0] = sin(theta) * sin(phi);
    R[2][1] = -sin(theta) * cos(phi);
    R[2][2] = cos(theta);

    t[0] = 0.0f;
    t[1] = 0.0f;
    t[2] = 6.0f;
}

void random_point(float & Xw, float & Yw, float & Zw)
{
    float theta = rand(0, 3.14159), phi = rand(0, 2 * 3.14159), R = rand(0, +2);

    Xw =  sin(theta) * sin(phi) * R;
    Yw = -sin(theta) * cos(phi) * R;
    Zw =  cos(theta) * R;
}

void project_with_noise(float R[3][3], float t[3],
float Xw, float Yw, float Zw,
float & u, float & v)
{
    float Xc = R[0][0] * Xw + R[0][1] * Yw + R[0][2] * Zw + t[0];
    float Yc = R[1][0] * Xw + R[1][1] * Yw + R[1][2] * Zw + t[1];
    float Zc = R[2][0] * Xw + R[2][1] * Yw + R[2][2] * Zw + t[2];

    //  float nu = rand(-noise, +noise);//with noise
    //  float nv = rand(-noise, +noise);
    //  u = uc + fu * Xc / Zc + nu;
    //  v = vc + fv * Yc / Zc + nv;
    u = uc + fu * Xc / Zc;//without noise
    v = vc + fv * Yc / Zc;
}

//#include "opencv2/opencv.hpp"

int main(int /*argc*/, char ** /*argv*/)
{
    epnp PnP;

    //srand(time(0));

    PnP.set_internal_parameters(uc, vc, fu, fv);
    PnP.set_maximum_number_of_correspondences(n);

    float R_true[3][3]={  0.0774722,-0.925266,-0.371322,
                           -0.985687,-0.015148,-0.167906,
                           0.149733,0.379015,-0.913196};
    float t_true[3]={0,0,6};

    //random_pose(R_true, t_true);
    //  for(int i=0;i<3;i++){
    //      t_true[i]=(float)i/10;
    //      for(int j=0;j<3;j++)
    //          R_true[i][j]=(float)i/10+(float)j/5;
    //  }

    PnP.reset_correspondences();

//    for(int i = 0; i < n; i++) {
//        float Xw, Yw, Zw, u, v;

//        random_point(Xw, Yw, Zw);
//        cout<<i<<":"<<Xw<<","<<Yw<<","<<Zw<<endl;

//        project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);

//        PnP.add_correspondence(Xw, Yw, Zw, u, v);
//    }

    float Xw, Yw, Zw, u, v;

    float te[3]={0.30323,0.728222,0.922768};
    Xw=te[0];Yw=te[1];Zw=te[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te1[3]={1.78068,0.295685,-0.283669};
    Xw=te1[0];Yw=te1[1];Zw=te1[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te2[3]={-0.099836,0.0399842,1.54034};
    Xw=te2[0];Yw=te2[1];Zw=te2[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te3[3]={-0.275392,0.088704,-0.506886};
    Xw=te3[0];Yw=te3[1];Zw=te3[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te4[3]={-0.485978,1.88505,0.251724};
    Xw=te4[0];Yw=te4[1];Zw=te4[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te5[3]={0.188951,0.169818,-1.04954};
    Xw=te5[0];Yw=te5[1];Zw=te5[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te6[3]={-0.772051,0.717037,0.412913};
    Xw=te6[0];Yw=te6[1];Zw=te6[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te7[3]={-0.337377,-1.6338,0.933979};
    Xw=te7[0];Yw=te7[1];Zw=te7[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te8[3]={0.0838386,-0.0543909,-1.50634};
    Xw=te8[0];Yw=te8[1];Zw=te8[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);

    float te9[3]={-0.0189287,-0.109389,-0.163641};
    Xw=te9[0];Yw=te9[1];Zw=te9[2];
    project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);






    float R_est[3][3], t_est[3];
    float err2 = PnP.compute_pose(R_est, t_est);
    float rot_err, transl_err;

    PnP.relative_error(rot_err, transl_err, R_true, t_true, R_est, t_est);
    cout << ">>> Reprojection error: " << err2 << endl;
    cout << ">>> rot_err: " << rot_err << ", transl_err: " << transl_err << endl;
    cout << endl;
    cout << "'True reprojection error':"
         << PnP.reprojection_error(R_true, t_true) << endl;
    cout << endl;
    cout << "True pose:" << endl;
    PnP.print_pose(R_true, t_true);
    cout << endl;
    cout << "Found pose:" << endl;
    PnP.print_pose(R_est, t_est);

	system("Pause");

    //  cv::Mat img;
    //  img=cv::imread("1.bmp");
    //  cv::imshow("Pi",img);
    //  cv::waitKey();


    return 0;
}
