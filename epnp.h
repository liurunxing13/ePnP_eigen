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

#ifndef epnp_h
#define epnp_h

//#include <opencv/cv.h>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Dense>

using namespace std;
//using namespace cv;

class epnp {
 public:
  epnp(void);
  ~epnp();

  void set_parameter(int n, const float * point_3d,const float * point_2d, const float * inter_matrix);

  void set_internal_parameters(const float uc, const float vc,
			       const float fu, const float fv);

  void set_maximum_number_of_correspondences(const int n);
  void reset_correspondences(void);
  void add_correspondence(const float X, const float Y, const float Z,
			  const float u, const float v);

  float compute_pose(float R[3][3], float T[3]);

  void relative_error(float & rot_err, float & transl_err,
		      const float Rtrue[3][3], const float ttrue[3],
		      const float Rest[3][3],  const float test[3]);

  void print_pose(const float R[3][3], const float t[3]);
  float reprojection_error(const float R[3][3], const float t[3]);

 private:
  void choose_control_points(void);
  void compute_barycentric_coordinates(void);
  void fill_M(Eigen::MatrixXf * M, const int row, const float * alphas, const float u, const float v);
  void compute_ccs(const float * betas, const float * ut);
  void compute_pcs(void);

  void solve_for_sign(void);

  void find_betas_approx_1(const Eigen::MatrixXf * L_6x10, const Eigen::Map<Eigen::MatrixXf> * Rho, float * betas);
  void find_betas_approx_2(const Eigen::MatrixXf * L_6x10, const Eigen::Map<Eigen::MatrixXf> * Rho, float * betas);
  void find_betas_approx_3(const Eigen::MatrixXf * L_6x10, const Eigen::Map<Eigen::MatrixXf> * Rho, float * betas);
  //void qr_solve(CvMat * A, CvMat * b, CvMat * X);

  float dot(const float * v1, const float * v2);
  float dist2(const float * p1, const float * p2);

  void compute_rho(float * rho);
  void compute_L_6x10(const float * ut, float * l_6x10);

  void gauss_newton(float * l_6x10, float * rho, float current_betas[4]);
  void compute_A_and_b_gauss_newton(const float * l_6x10, const float * rho,
	  float betas[4], Eigen::Map<Eigen::MatrixXf> * A, Eigen::Map<Eigen::MatrixXf> * b);

  float compute_R_and_t(const float * ut, const float * betas,
			 float R[3][3], float t[3]);

  void estimate_R_and_t(float R[3][3], float t[3]);

  void copy_R_and_t(const float R_dst[3][3], const float t_dst[3],
		    float R_src[3][3], float t_src[3]);

  void mat_to_quat(const float R[3][3], float q[4]);

  //相机内参
  float uc, vc, fu, fv;
  //pws 表示世界坐标系下的3D参考点坐标，大小为3*n
  //us 表示图像坐标系下的2D参考点坐标，大小为2*n
  //alphas 表示齐次重心坐标，大小为4*n
  //pcs 表示摄像机坐标系下的3D点坐标，摄像机坐标系下的控制点坐标的线性组合得到。大小为3^n
  float * pws, * us, * alphas, * pcs;
  int maximum_number_of_correspondences;
  int number_of_correspondences;

  //cws为世界坐标系下的4个控制点(基)的3D坐标
  //ccs为摄像机坐标系下的4个控制点(基)的3D坐标
  float cws[4][3], ccs[4][3];
  float cws_determinant;
};

#endif
