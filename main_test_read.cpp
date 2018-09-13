//epnp using example

#include "epnp.h"
#include <time.h>
#include "fstream"

//#pragma comment(lib, "Epnp_Eigen.lib") //

int main(int /*argc*/, char ** /*argv*/)
{

	const int n = 10;
	float P3[30], P2[20], iM[4];
	float R_est[3][3], t_est[3];
    epnp PnP;

	ifstream f_3d, f_2d, i_M;
	f_3d.open("point_3d.txt");
	f_2d.open("point_2d.txt");
	i_M.open("inter_matrix.txt");
	for (int i = 0; i < 4; i++)
		i_M >> iM[i];
	for (int i = 0; i < 30; i++)
		f_3d >> P3[i];
	for (int i = 0; i < 20; i++)
		f_2d >> P2[i];

	PnP.set_parameter(n, P3, P2, iM);

	clock_t start, stop;
	start = clock();
	int num = 1000;
	while(num--)
		float err2 = PnP.compute_pose(R_est, t_est);
	stop = clock();
	cout << "Programm run 1000 times:" << (float)(stop - start) << "ms" << endl;

	

    //float rot_err, transl_err;
    //PnP.relative_error(rot_err, transl_err, R_true, t_true, R_est, t_est);
    //cout << ">>> Reprojection error: " << err2 << endl;
    //cout << ">>> rot_err: " << rot_err << ", transl_err: " << transl_err << endl;
    //cout << endl;
    //cout << "'True reprojection error':"
    //     << PnP.reprojection_error(R_true, t_true) << endl;
    //cout << endl;
    //cout << "True pose:" << endl;
    //PnP.print_pose(R_true, t_true);
    //cout << endl;
    cout << "Found pose:" << endl;
    PnP.print_pose(R_est, t_est);

	system("Pause");

    return 0;
}
