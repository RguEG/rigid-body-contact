#include <inverse_rigid_body.h>
//#include <rigid_body_jacobian.h>

void inverse_rigid_body(Eigen::Vector3d &x, Eigen::Ref<const Eigen::Vector3d> x_world, 
                        Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p) {
    //Input:
//  x - world space position
//  R - rotation from undeformed to world space
//  p - center-of-mass translation 
//Output:
//  X - undeformed position of point x
//  qdot - updated generalized velocities 

	//Eigen::Matrix36d J;
	Eigen::SparseMatrixd A;
	Eigen::Vector3d b;
	A = R.sparseView();
	b = x_world - p;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(b);

	//rigid_body_jacobian(J,R,p,x_world);
    
	//qdot = 
}