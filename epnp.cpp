// implement of class epnp.h

#include "epnp.h"

epnp::epnp(void)
{
	maximum_number_of_correspondences = 0;
	number_of_correspondences = 0;

	pws = 0;
	us = 0;
	alphas = 0;
	pcs = 0;
}

epnp::~epnp()
{
	delete[] pws;
	delete[] us;
	delete[] alphas;
	delete[] pcs;
}

void epnp::set_internal_parameters(float uc, float vc, float fu, float fv)
{
	this->uc = uc;
	this->vc = vc;
	this->fu = fu;
	this->fv = fv;
}

void epnp::set_maximum_number_of_correspondences(int n)
{
	if (maximum_number_of_correspondences < n) {
		if (pws != 0) delete[] pws;
		if (us != 0) delete[] us;
		if (alphas != 0) delete[] alphas;
		if (pcs != 0) delete[] pcs;

		maximum_number_of_correspondences = n;
		pws = new float[3 * maximum_number_of_correspondences];
		us = new float[2 * maximum_number_of_correspondences];
		alphas = new float[4 * maximum_number_of_correspondences];
		pcs = new float[3 * maximum_number_of_correspondences];
	}
}

void epnp::reset_correspondences(void)
{
	number_of_correspondences = 0;
}

void epnp::add_correspondence(float X, float Y, float Z, float u, float v)
{
	pws[3 * number_of_correspondences] = X;
	pws[3 * number_of_correspondences + 1] = Y;
	pws[3 * number_of_correspondences + 2] = Z;

	us[2 * number_of_correspondences] = u;
	us[2 * number_of_correspondences + 1] = v;

	number_of_correspondences++;
}

void epnp::choose_control_points(void)
{
	// Take C0 as the reference points centroid:
	cws[0][0] = cws[0][1] = cws[0][2] = 0;
	for (int i = 0; i < number_of_correspondences; i++)
		for (int j = 0; j < 3; j++)
			cws[0][j] += pws[3 * i + j];

	for (int j = 0; j < 3; j++)
		cws[0][j] /= number_of_correspondences;

	Eigen::MatrixX3f PW0(number_of_correspondences, 3);

	float pw0tpw0[3 * 3], dc[3], uct[3 * 3];
	Eigen::Map<Eigen::MatrixXf> PW0tPW0(pw0tpw0, 3, 3);
	Eigen::Map<Eigen::MatrixXf> DC(dc, 3, 1);
	Eigen::Map<Eigen::MatrixXf> UCt(uct, 3, 3);

	for (int i = 0; i < number_of_correspondences; i++) {
		for (int j = 0; j < 3; j++) {
			PW0(i, j) = pws[3 * i + j] - cws[0][j];
		}
	}
	PW0tPW0 = PW0.transpose()*PW0;
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(PW0tPW0, Eigen::ComputeThinU | Eigen::ComputeThinV);
	DC = svd.singularValues();//dc������Ϊ�Ǹ�����
	UCt = svd.matrixU();
	Eigen::Matrix3f VCt_t = svd.matrixV();

	for (int i = 1; i < 4; i++) {
		float k = sqrt(DC(i - 1) / number_of_correspondences);
		for (int j = 0; j < 3; j++) {
			cws[i][j] = cws[0][j] + k * UCt((i - 1), j);
			//cout << "cws:\n" << cws[i][j] << endl;
		}
	}

}

void epnp::compute_barycentric_coordinates(void)
{
	Eigen::Matrix3f CC, CC_inv;


	for (int i = 0; i < 3; i++)
		for (int j = 1; j < 4; j++)
			CC(i, j - 1) = cws[j][i] - cws[0][i];

	CC_inv = CC.inverse();

	for (int i = 0; i < number_of_correspondences; i++) {
		float * pi = pws + 3 * i;
		float * a = alphas + 4 * i;

		for (int j = 0; j < 3; j++) {
			a[1 + j] =
				CC_inv(j, 0) * (pi[0] - cws[0][0]) +
				CC_inv(j, 1) * (pi[1] - cws[0][1]) +
				CC_inv(j, 2) * (pi[2] - cws[0][2]);
			//cout << "a:\n" << a[1 + j] << endl;

		}
		a[0] = 1.0f - a[1] - a[2] - a[3];
	}
}

void epnp::fill_M(Eigen::MatrixXf * M,
	const int row, const float * as, const float u, const float v)
{
	for (int i = 0; i < 4; i++)
	{
		(*M)(row, 3 * i) = as[i] * fu;
		(*M)(row, 3 * i + 1) = 0.0;
		(*M)(row, 3 * i + 2) = as[i] * (uc - u);

		(*M)(row + 1, 3 * i) = 0.0;
		(*M)(row + 1, 3 * i + 1) = as[i] * fv;
		(*M)(row + 1, 3 * i + 2) = as[i] * (vc - v);
	}
}

void epnp::compute_ccs(const float * betas, const float * ut)
{
	for (int i = 0; i < 4; i++)
		ccs[i][0] = ccs[i][1] = ccs[i][2] = 0.0f;

	for (int i = 0; i < 4; i++) {
		const float * v = ut + 12 * (11 - i);
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 3; k++)
				ccs[j][k] += betas[i] * v[3 * j + k];
	}
}

void epnp::compute_pcs(void)
{
	for (int i = 0; i < number_of_correspondences; i++) {
		float * a = alphas + 4 * i;
		float * pc = pcs + 3 * i;

		for (int j = 0; j < 3; j++)
			pc[j] = a[0] * ccs[0][j] + a[1] * ccs[1][j] + a[2] * ccs[2][j] + a[3] * ccs[3][j];
	}
}

float epnp::compute_pose(float R[3][3], float t[3])
{
	choose_control_points();
	compute_barycentric_coordinates();
	Eigen::MatrixXf M(2 * number_of_correspondences, 12);//M size n*12;

														 //cout<<number_of_correspondences<<endl;

	for (int i = 0; i < number_of_correspondences; i++)
	{
		//cout << "alphas:" << *(alphas + 4 * i) << "us:" << us[2 * i] << endl;
		fill_M(&M, 2 * i, alphas + 4 * i, us[2 * i], us[2 * i + 1]);
	}
	//cout << "M" << M << endl;

	float ut[12 * 12];
	Eigen::MatrixXf MtM(12, 12);
	Eigen::MatrixXf D(12, 1);
	Eigen::Map<Eigen::MatrixXf> Ut(ut, 12, 12);
	MtM = M.transpose()*M;
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(MtM, Eigen::ComputeThinU | Eigen::ComputeThinV);
	D = svd.singularValues();
	Ut = svd.matrixU();


	float l_6x10[6 * 10], rho[6];
	compute_L_6x10(ut, l_6x10);
	compute_rho(rho);

	Eigen::Map<Eigen::MatrixXf> L_6x10_t(l_6x10, 10, 6);
	Eigen::Map<Eigen::MatrixXf> Rho(rho, 6, 1);

	Eigen::MatrixXf L_6x10 = L_6x10_t.transpose();
	//cout << "L_6x10:" << L_6x10 << endl;
	//cout << "Rho:" << Rho << endl;


	float Betas[4][4], rep_errors[4];
	float Rs[4][3][3], ts[4][3];

	find_betas_approx_1(&L_6x10, &Rho, Betas[1]);
	gauss_newton(l_6x10, rho, Betas[1]);
	rep_errors[1] = compute_R_and_t(ut, Betas[1], Rs[1], ts[1]);

	find_betas_approx_2(&L_6x10, &Rho, Betas[2]);
	gauss_newton(l_6x10, rho, Betas[2]);
	rep_errors[2] = compute_R_and_t(ut, Betas[2], Rs[2], ts[2]);

	find_betas_approx_3(&L_6x10, &Rho, Betas[3]);
	gauss_newton(l_6x10, rho, Betas[3]);
	rep_errors[3] = compute_R_and_t(ut, Betas[3], Rs[3], ts[3]);

	int N = 1;
	if (rep_errors[2] < rep_errors[1]) N = 2;
	if (rep_errors[3] < rep_errors[N]) N = 3;

	copy_R_and_t(Rs[N], ts[N], R, t);

	return rep_errors[N];
}

void epnp::copy_R_and_t(const float R_src[3][3], const float t_src[3],
	float R_dst[3][3], float t_dst[3])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			R_dst[i][j] = R_src[i][j];
		t_dst[i] = t_src[i];
	}
}

float epnp::dist2(const float * p1, const float * p2)
{
	return
		(p1[0] - p2[0]) * (p1[0] - p2[0]) +
		(p1[1] - p2[1]) * (p1[1] - p2[1]) +
		(p1[2] - p2[2]) * (p1[2] - p2[2]);
}

float epnp::dot(const float * v1, const float * v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

float epnp::reprojection_error(const float R[3][3], const float t[3])
{
	float sum2 = 0.0;

	for (int i = 0; i < number_of_correspondences; i++) {
		float * pw = pws + 3 * i;
		float Xc = dot(R[0], pw) + t[0];
		float Yc = dot(R[1], pw) + t[1];
		float inv_Zc = 1.0 / (dot(R[2], pw) + t[2]);
		float ue = uc + fu * Xc * inv_Zc;
		float ve = vc + fv * Yc * inv_Zc;
		float u = us[2 * i], v = us[2 * i + 1];

		sum2 += sqrt((u - ue) * (u - ue) + (v - ve) * (v - ve));
	}

	return sum2 / number_of_correspondences;
}

void epnp::estimate_R_and_t(float R[3][3], float t[3])
{
	float pc0[3], pw0[3];

	pc0[0] = pc0[1] = pc0[2] = 0.0;
	pw0[0] = pw0[1] = pw0[2] = 0.0;

	for (int i = 0; i < number_of_correspondences; i++) {
		const float * pc = pcs + 3 * i;
		const float * pw = pws + 3 * i;

		for (int j = 0; j < 3; j++) {
			pc0[j] += pc[j];
			pw0[j] += pw[j];
		}
	}
	for (int j = 0; j < 3; j++) {
		pc0[j] /= number_of_correspondences;
		pw0[j] /= number_of_correspondences;
	}

	float abt[3 * 3], abt_d[3], abt_u[3 * 3], abt_v[3 * 3];
	Eigen::Map<Eigen::MatrixXf> ABt_t(abt, 3, 3);
	Eigen::Map<Eigen::MatrixXf> ABt_D(abt_d, 3, 1);
	Eigen::Map<Eigen::MatrixXf> ABt_U(abt_u, 3, 3);
	Eigen::Map<Eigen::MatrixXf> ABt_V(abt_v, 3, 3);
	ABt_t.setZero();

	for (int i = 0; i < number_of_correspondences; i++) {
		float * pc = pcs + 3 * i;
		float * pw = pws + 3 * i;
		for (int j = 0; j < 3; j++) {
			abt[3 * j] += (pc[j] - pc0[j]) * (pw[0] - pw0[0]);
			abt[3 * j + 1] += (pc[j] - pc0[j]) * (pw[1] - pw0[1]);
			abt[3 * j + 2] += (pc[j] - pc0[j]) * (pw[2] - pw0[2]);
		}
	}

	Eigen::JacobiSVD<Eigen::MatrixXf> svd(ABt_t.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
	ABt_D = svd.singularValues();
	ABt_U.transpose() = svd.matrixU();
	ABt_V.transpose() = svd.matrixV();

	//cout << "ABt:" << ABt_t.transpose()<<endl;
	//cout << "ABt_U:" << ABt_U << endl;
	//cout << "ABt_V:" << ABt_V << endl;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = dot(abt_u + 3 * i, abt_v + 3 * j);

	const float det =
		R[0][0] * R[1][1] * R[2][2] + R[0][1] * R[1][2] * R[2][0] + R[0][2] * R[1][0] * R[2][1] -
		R[0][2] * R[1][1] * R[2][0] - R[0][1] * R[1][0] * R[2][2] - R[0][0] * R[1][2] * R[2][1];

	if (det < 0) {
		R[2][0] = -R[2][0];
		R[2][1] = -R[2][1];
		R[2][2] = -R[2][2];
	}

	t[0] = pc0[0] - dot(R[0], pw0);
	t[1] = pc0[1] - dot(R[1], pw0);
	t[2] = pc0[2] - dot(R[2], pw0);
}

void epnp::print_pose(const float R[3][3], const float t[3])
{
	cout << R[0][0] << " " << R[0][1] << " " << R[0][2] << " " << t[0] << endl;
	cout << R[1][0] << " " << R[1][1] << " " << R[1][2] << " " << t[1] << endl;
	cout << R[2][0] << " " << R[2][1] << " " << R[2][2] << " " << t[2] << endl;
}

void epnp::solve_for_sign(void)
{
	if (pcs[2] < 0.0) {
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 3; j++)
				ccs[i][j] = -ccs[i][j];

		for (int i = 0; i < number_of_correspondences; i++) {
			pcs[3 * i] = -pcs[3 * i];
			pcs[3 * i + 1] = -pcs[3 * i + 1];
			pcs[3 * i + 2] = -pcs[3 * i + 2];
		}
	}
}

float epnp::compute_R_and_t(const float * ut, const float * betas,
	float R[3][3], float t[3])
{
	compute_ccs(betas, ut);
	compute_pcs();

	solve_for_sign();

	estimate_R_and_t(R, t);

	return reprojection_error(R, t);
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_1 = [B11 B12     B13         B14]

void epnp::find_betas_approx_1(const Eigen::MatrixXf * L_6x10, const Eigen::Map<Eigen::MatrixXf> * Rho,
	float * betas)
{
	float l_6x4[6 * 4], b4[4];
	Eigen::Map<Eigen::MatrixXf> L_6x4(l_6x4, 6, 4);
	Eigen::Map<Eigen::MatrixXf> B4(b4, 4, 1);

	for (int i = 0; i < 6; i++) {
		L_6x4(i, 0) = (*L_6x10)(i, 0);
		L_6x4(i, 1) = (*L_6x10)(i, 1);
		L_6x4(i, 2) = (*L_6x10)(i, 3);
		L_6x4(i, 3) = (*L_6x10)(i, 6);
	}

	B4 = L_6x4.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(*Rho);

	//cout << "L_6x4:" << L_6x4 << "B4:" << B4<<endl;
	//cout << "L_610:" << (*L_6x10) << endl;

	if (b4[0] < 0) {
		betas[0] = sqrt(-b4[0]);
		betas[1] = -b4[1] / betas[0];
		betas[2] = -b4[2] / betas[0];
		betas[3] = -b4[3] / betas[0];
	}
	else {
		betas[0] = sqrt(b4[0]);
		betas[1] = b4[1] / betas[0];
		betas[2] = b4[2] / betas[0];
		betas[3] = b4[3] / betas[0];
	}
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_2 = [B11 B12 B22                            ]

void epnp::find_betas_approx_2(const Eigen::MatrixXf * L_6x10, const Eigen::Map<Eigen::MatrixXf> * Rho,
	float * betas)
{
	float l_6x3[6 * 3], b3[3];
	Eigen::Map<Eigen::MatrixXf> L_6x3(l_6x3, 6, 3);
	Eigen::Map<Eigen::MatrixXf> B3(b3, 3, 1);
	L_6x3 = (*L_6x10).block(0, 0, 6, 3);


	B3 = L_6x3.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(*Rho);

	if (b3[0] < 0) {
		betas[0] = sqrt(-b3[0]);
		betas[1] = (b3[2] < 0) ? sqrt(-b3[2]) : 0.0;
	}
	else {
		betas[0] = sqrt(b3[0]);
		betas[1] = (b3[2] > 0) ? sqrt(b3[2]) : 0.0;
	}

	if (b3[1] < 0) betas[0] = -betas[0];

	betas[2] = 0.0;
	betas[3] = 0.0;
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_3 = [B11 B12 B22 B13 B23                    ]

void epnp::find_betas_approx_3(const Eigen::MatrixXf * L_6x10, const Eigen::Map<Eigen::MatrixXf> * Rho,
	float * betas)
{
	float l_6x5[6 * 5], b5[5];
	Eigen::Map<Eigen::MatrixXf> L_6x5(l_6x5, 6, 5);
	Eigen::Map<Eigen::MatrixXf> B5(b5, 5, 1);
	L_6x5 = (*L_6x10).block(0, 0, 6, 5);


	B5 = L_6x5.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(*Rho);

	if (b5[0] < 0) {
		betas[0] = sqrt(-b5[0]);
		betas[1] = (b5[2] < 0) ? sqrt(-b5[2]) : 0.0;
	}
	else {
		betas[0] = sqrt(b5[0]);
		betas[1] = (b5[2] > 0) ? sqrt(b5[2]) : 0.0;
	}
	if (b5[1] < 0) betas[0] = -betas[0];
	betas[2] = b5[3] / betas[0];
	betas[3] = 0.0;
}

void epnp::compute_L_6x10(const float * ut, float * l_6x10)
{
	const float * v[4];

	v[0] = ut + 12 * 11;
	v[1] = ut + 12 * 10;
	v[2] = ut + 12 * 9;
	v[3] = ut + 12 * 8;

	float dv[4][6][3];

	for (int i = 0; i < 4; i++) {
		int a = 0, b = 1;
		for (int j = 0; j < 6; j++) {
			dv[i][j][0] = v[i][3 * a] - v[i][3 * b];
			dv[i][j][1] = v[i][3 * a + 1] - v[i][3 * b + 1];
			dv[i][j][2] = v[i][3 * a + 2] - v[i][3 * b + 2];

			b++;
			if (b > 3) {
				a++;
				b = a + 1;
			}
		}
	}

	for (int i = 0; i < 6; i++) {
		float * row = l_6x10 + 10 * i;

		row[0] = dot(dv[0][i], dv[0][i]);
		row[1] = 2.0f * dot(dv[0][i], dv[1][i]);
		row[2] = dot(dv[1][i], dv[1][i]);
		row[3] = 2.0f * dot(dv[0][i], dv[2][i]);
		row[4] = 2.0f * dot(dv[1][i], dv[2][i]);
		row[5] = dot(dv[2][i], dv[2][i]);
		row[6] = 2.0f * dot(dv[0][i], dv[3][i]);
		row[7] = 2.0f * dot(dv[1][i], dv[3][i]);
		row[8] = 2.0f * dot(dv[2][i], dv[3][i]);
		row[9] = dot(dv[3][i], dv[3][i]);
	}
}

void epnp::compute_rho(float * rho)
{
	rho[0] = dist2(cws[0], cws[1]);
	rho[1] = dist2(cws[0], cws[2]);
	rho[2] = dist2(cws[0], cws[3]);
	rho[3] = dist2(cws[1], cws[2]);
	rho[4] = dist2(cws[1], cws[3]);
	rho[5] = dist2(cws[2], cws[3]);
}

void epnp::compute_A_and_b_gauss_newton(const float * l_6x10, const float * rho,
	float betas[4], Eigen::Map<Eigen::MatrixXf> * A, Eigen::Map<Eigen::MatrixXf> * b)
{
	//cout << endl << "l_6x10:" << endl;
	//for (int i = 0; i<60; i++)
	//{
	//	cout << l_6x10[i] << "  ";
	//}
	//cout << endl << "rho:" << endl;
	//for (int i = 0; i<6; i++)
	//{
	//	cout << rho[i] << "  ";
	//}

	for (int i = 0; i < 6; i++) {
		const float * rowL = l_6x10 + i * 10;
		float * rowA = A->data() + i * 4;

		rowA[0] = 2 * rowL[0] * betas[0] + rowL[1] * betas[1] + rowL[3] * betas[2] + rowL[6] * betas[3];
		rowA[1] = rowL[1] * betas[0] + 2 * rowL[2] * betas[1] + rowL[4] * betas[2] + rowL[7] * betas[3];
		rowA[2] = rowL[3] * betas[0] + rowL[4] * betas[1] + 2 * rowL[5] * betas[2] + rowL[8] * betas[3];
		rowA[3] = rowL[6] * betas[0] + rowL[7] * betas[1] + rowL[8] * betas[2] + 2 * rowL[9] * betas[3];

		(*b)(i) = rho[i] -
			(
				rowL[0] * betas[0] * betas[0] +
				rowL[1] * betas[0] * betas[1] +
				rowL[2] * betas[1] * betas[1] +
				rowL[3] * betas[0] * betas[2] +
				rowL[4] * betas[1] * betas[2] +
				rowL[5] * betas[2] * betas[2] +
				rowL[6] * betas[0] * betas[3] +
				rowL[7] * betas[1] * betas[3] +
				rowL[8] * betas[2] * betas[3] +
				rowL[9] * betas[3] * betas[3]
				);
	}
}

void epnp::gauss_newton(float * l_6x10, float * rho, float current_betas[4])
{
	//for (int i = 0; i < 4; i++)
	//	cout << endl << "betes:" << current_betas[i];

	const int iterations_number = 5;

	float a[6 * 4], b[6], x[4];
	Eigen::Map<Eigen::MatrixXf> A_t(a, 4, 6);//sotage as colcum;
	Eigen::Map<Eigen::MatrixXf> B(b, 6, 1);
	Eigen::Map<Eigen::MatrixXf> X(x, 4, 1);

	for (int k = 0; k < iterations_number; k++) {
		compute_A_and_b_gauss_newton(l_6x10, rho,
			current_betas, &A_t, &B);		
		X = A_t.transpose().colPivHouseholderQr().solve(B);

		for (int i = 0; i < 4; i++)
			current_betas[i] += x[i];
	}

	//for (int i = 0; i < 4; i++)
	//	cout<<endl<<"After iterate betes:"<<current_betas[i];	
}


void epnp::relative_error(float & rot_err, float & transl_err,
	const float Rtrue[3][3], const float ttrue[3],
	const float Rest[3][3], const float test[3])
{
	float qtrue[4], qest[4];

	mat_to_quat(Rtrue, qtrue);
	mat_to_quat(Rest, qest);

	float rot_err1 = sqrt((qtrue[0] - qest[0]) * (qtrue[0] - qest[0]) +
		(qtrue[1] - qest[1]) * (qtrue[1] - qest[1]) +
		(qtrue[2] - qest[2]) * (qtrue[2] - qest[2]) +
		(qtrue[3] - qest[3]) * (qtrue[3] - qest[3])) /
		sqrt(qtrue[0] * qtrue[0] + qtrue[1] * qtrue[1] + qtrue[2] * qtrue[2] + qtrue[3] * qtrue[3]);

	float rot_err2 = sqrt((qtrue[0] + qest[0]) * (qtrue[0] + qest[0]) +
		(qtrue[1] + qest[1]) * (qtrue[1] + qest[1]) +
		(qtrue[2] + qest[2]) * (qtrue[2] + qest[2]) +
		(qtrue[3] + qest[3]) * (qtrue[3] + qest[3])) /
		sqrt(qtrue[0] * qtrue[0] + qtrue[1] * qtrue[1] + qtrue[2] * qtrue[2] + qtrue[3] * qtrue[3]);

	rot_err = min(rot_err1, rot_err2);

	transl_err =
		sqrt((ttrue[0] - test[0]) * (ttrue[0] - test[0]) +
		(ttrue[1] - test[1]) * (ttrue[1] - test[1]) +
			(ttrue[2] - test[2]) * (ttrue[2] - test[2])) /
		sqrt(ttrue[0] * ttrue[0] + ttrue[1] * ttrue[1] + ttrue[2] * ttrue[2]);
}

void epnp::mat_to_quat(const float R[3][3], float q[4])
{
	float tr = R[0][0] + R[1][1] + R[2][2];
	float n4;

	if (tr > 0.0f) {
		q[0] = R[1][2] - R[2][1];
		q[1] = R[2][0] - R[0][2];
		q[2] = R[0][1] - R[1][0];
		q[3] = tr + 1.0f;
		n4 = q[3];
	}
	else if ((R[0][0] > R[1][1]) && (R[0][0] > R[2][2])) {
		q[0] = 1.0f + R[0][0] - R[1][1] - R[2][2];
		q[1] = R[1][0] + R[0][1];
		q[2] = R[2][0] + R[0][2];
		q[3] = R[1][2] - R[2][1];
		n4 = q[0];
	}
	else if (R[1][1] > R[2][2]) {
		q[0] = R[1][0] + R[0][1];
		q[1] = 1.0f + R[1][1] - R[0][0] - R[2][2];
		q[2] = R[2][1] + R[1][2];
		q[3] = R[2][0] - R[0][2];
		n4 = q[1];
	}
	else {
		q[0] = R[2][0] + R[0][2];
		q[1] = R[2][1] + R[1][2];
		q[2] = 1.0f + R[2][2] - R[0][0] - R[1][1];
		q[3] = R[0][1] - R[1][0];
		n4 = q[2];
	}
	float scale = 0.5f / float(sqrt(n4));

	q[0] *= scale;
	q[1] *= scale;
	q[2] *= scale;
	q[3] *= scale;
}

void epnp::set_parameter(int n, const float * point_3d, const float * point_2d, const float * inter_matrix)
{
	/*
	n is point number	
	point_3d storage as XW1,YW1,ZW1,XW2,YW2,ZW2.....Size 3*n
	point_2d storage as u1, v1, u2, v2.....Size 2*n
	inter_matrix storage as fx,fy,u0,v0 Size 4*1
	*/
	set_internal_parameters(inter_matrix[0], inter_matrix[1], inter_matrix[2], inter_matrix[3]);
	set_maximum_number_of_correspondences(n);
	reset_correspondences();
	float Xw, Yw, Zw, u, v;
	for (int i = 0; i < n; i++)
	{
		Xw = point_3d[3*i + 0] ;
		Yw = point_3d[3*i + 1] ;
		Zw = point_3d[3*i + 2] ;
		u = point_2d[2*i + 0] ;
		v = point_2d[2*i + 1] ;
		add_correspondence(Xw, Yw, Zw, u, v);
	}
}
