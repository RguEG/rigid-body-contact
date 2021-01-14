#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>
#include <rigid_body_jacobian.h>
#include <inverse_rigid_body.h>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//  n - list of collision normals
//  x - list of world space collision points
//  obj - list of collision object ids 
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
inline void exponential_euler_lcp_contact(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces,
                            std::vector<Eigen::Vector3d> &n, std::vector<Eigen::Vector3d> &x, std::vector<std::pair<int,int> > &obj) {
	for (int i = 0; i < q.rows() / 12; i++) {
		Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12.0 * i).data());

		Eigen::Matrix3d inertiaM, expR;

		inertiaM.setZero();
		expR.setZero();

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				inertiaM(m, n) = masses.at(i)(m, n);
			}
		}

		Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12 * i + 9.0).data());

		Eigen::Vector3d w = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6.0 * i).data());
		Eigen::Vector3d pdot = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6.0 * i + 3.0).data());
		Eigen::Vector3d text = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6.0 * i).data());
		Eigen::Vector3d fext = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6.0 * i + 3.0).data());
		Eigen::Vector6d qdoti = Eigen::Map<const Eigen::Vector6d>(qdot.segment<6>(6.0 * i).data());

		rodrigues(expR, dt * w);
		Eigen::SparseMatrixd A;
		Eigen::Vector3d b, u, k;
		k = (R * inertiaM * R.transpose()) * w;
		k = (dt * w).cross(k);
		A = (R * inertiaM * R.transpose()).sparseView();
		b = (R * inertiaM * R.transpose()) * w + k + dt * text;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		u = solver.solve(b);
		w = u;
		pdot = pdot + (1.0 / masses.at(i)(5, 5)) * dt * fext;

		qdot[6.0 * i] = w[0];
		qdot[6.0 * i + 1.0] = w[1];
		qdot[6.0 * i + 2.0] = w[2];

		qdot[6.0 * i + 3.0] = pdot[0];
		qdot[6.0 * i + 4.0] = pdot[1];
		qdot[6.0 * i + 5.0] = pdot[2];

		Eigen::MatrixXd RM;
		RM.resize(masses.at(i).rows(), masses.at(i).cols());
		RM = masses.at(i);
		Eigen::Matrix3d inertiaRM;
		inertiaRM = R * inertiaM * R.transpose();

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				RM(m, n) = inertiaRM(m, n);
			}
		}

		Eigen::Vector3d lastfTop;
		lastfTop = w.cross(inertiaRM * w) + text;
		Eigen::Vector6d lastf;

		for (int m = 0; m < 3; m++) {
			lastf[m] = lastfTop[m];
		}

		for (int n = 3; n < 6; n++) {
			lastf[n] = fext[n - 3];
		}

		Eigen::VectorXd alpha;
		alpha.resize(x.size());
		alpha.setZero();
		Eigen::Vector6d fTotal;
		fTotal = lastf;

		for (int k = 0; k < 10; k++) {
			if (x.size() == 0) {
				break;
			}
			for (int z = 0; z < x.size(); z++) {

				Eigen::VectorXd qdotStar;

				qdotStar = qdot + dt * RM.inverse() * fTotal;

				double alphai = 0;
				Eigen::VectorXd delta, gamma;
				Eigen::MatrixXd g;
				Eigen::Matrix36d j;
				Eigen::Vector3d x_undeformed;
				inverse_rigid_body(x_undeformed, x.at(z), R, p);
				rigid_body_jacobian(j, R, p, x.at(z));
				g = j.transpose() * n.at(z);

				delta = dt * (g.transpose() * RM.inverse() * g);
				gamma = g.transpose() * (qdotStar + dt * RM.inverse() * fTotal);

				alphai = -gamma[0] / delta[0];

				if (alphai < 0) {
					alphai = 0;
				}
				alpha[z] = alphai;
			}
			Eigen::VectorXd alphaTog;
			alphaTog.resize(qdot.size());
			alphaTog.setZero();
			for (int z = 0; z < x.size(); z++) {

				Eigen::Matrix36d j;
				rigid_body_jacobian(j, R, p, x.at(z));
				alphaTog += alpha[z] * j.transpose() * n.at(z);
			}

			qdot = qdot + dt * RM.inverse() * fTotal + dt * RM.inverse() * alphaTog;
			w = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6.0 * i).data());
			pdot = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6.0 * i + 3.0).data());

			lastfTop = w.cross(inertiaRM * w) + text;

			for (int m = 0; m < 3; m++) {
				lastf[m] = lastfTop[m];
			}
			fTotal = lastf;
		}

		R = expR * R;

		p = p + dt * pdot;

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				q[12 * i + 3.0 * m + n] = R(n, m);
			}
		}
		q[12.0 * i + 9.0] = p[0];
		q[12.0 * i + 10.0] = p[1];
		q[12.0 * i + 11.0] = p[2];
	}
}