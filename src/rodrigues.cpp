#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
	Eigen::Vector3d a;
	Eigen::Matrix3d identI, crossA;
	double absOmega;

	identI.setIdentity();

	absOmega = sqrt(pow(omega[0], 2) + pow(omega[1], 2) + pow(omega[2], 2));

	a = omega * (1.0 / absOmega);

	crossA.setZero();
	crossA(0, 1) = -a[2];
	crossA(0, 2) = a[1];
	crossA(1, 0) = a[2];
	crossA(1, 2) = -a[0];
	crossA(2, 0) = -a[1];
	crossA(2, 1) = a[0];

	R = identI + sin(absOmega) * crossA + (1 - cos(absOmega)) * crossA * crossA;
    

}