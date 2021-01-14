#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {

	Eigen::Matrix3d identI, crossA, crossAT, RT;
	Eigen::Matrix36d copB;
	Eigen::Matrix66d copC;

	identI.setIdentity();
	crossA.setZero();
	copC.setZero();
	crossA(0, 1) = -x[2];
	crossA(0, 2) = x[1];
	crossA(1, 0) = x[2];
	crossA(1, 2) = -x[0];
	crossA(2, 0) = -x[1];
	crossA(2, 1) = x[0];
	crossAT = crossA.transpose();
	RT = R.transpose();

	for (int m = 0; m < 3; m++) {
		for (int n = 0; n < 3; n++) {
			copB(m, n) = crossAT(m, n);
		}
	}

	for (int m = 0; m < 3; m++) {
		for (int n = 3; n < 6; n++) {
			copB(m, n) = identI(m, n - 3);
		}
	}

	for (int m = 0; m < 3; m++) {
		for (int n = 0; n < 3; n++) {
			copC(m, n) = RT(m, n);
		}
	}

	for (int m = 3; m < 6; m++) {
		for (int n = 3; n < 6; n++) {
			copC(m, n) = RT(m - 3, n - 3);
		}
	}

	J = R * copB * copC;
    
}