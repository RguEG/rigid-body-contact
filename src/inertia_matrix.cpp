#include <inertia_matrix.h>
#include <cassert>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
    Eigen::VectorXd N;
    Eigen::Vector3d a, b, c, temC, crossA, crossB;
    double absN, integralX, area;

    double side[3];

    temC.setZero();
    I.setZero();

    for (int i = 0; i < F.rows(); i++) {

        a = V.row(F(i, 0));
        b = V.row(F(i, 1));
        c = V.row(F(i, 2));

        side[0] = sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
        side[1] = sqrt(pow(a[0] - c[0], 2) + pow(a[1] - c[1], 2) + pow(a[2] - c[2], 2));
        side[2] = sqrt(pow(c[0] - b[0], 2) + pow(c[1] - b[1], 2) + pow(c[2] - b[2], 2));

        double p = (side[0] + side[1] + side[2]) / 2;
        area = sqrt(p * (p - side[0]) * (p - side[1]) * (p - side[2]));

        crossA = b - a;
        crossB = c - a;
        N = ((crossA).cross(crossB));
        absN = pow(N[0], 2.0) + pow(N[1], 2.0) + pow(N[2], 2.0);
        N = (crossA).cross(crossB) * (1.0 / sqrt(absN));

        integralX = 2 * area * density * (1.0 / 6.0) * N[0] * (a[0] + b[0] + c[0]);

        mass += integralX;

        I(0, 0) += 2 * area * density * N[1] * (1.0 / 3.0) * ((a[1] * (b[1] * b[1])) / 2.0E+1 + ((a[1] * a[1]) * b[1]) / 2.0E+1 + (a[1] * (c[1] * c[1])) / 2.0E+1 + ((a[1] * a[1]) * c[1]) / 2.0E+1 + (b[1] * (c[1] * c[1])) / 2.0E+1 + ((b[1] * b[1]) * c[1]) / 2.0E+1 + (a[1] * a[1] * a[1]) / 2.0E+1 + (b[1] * b[1] * b[1]) / 2.0E+1 + (c[1] * c[1] * c[1]) / 2.0E+1 + (a[1] * b[1] * c[1]) / 2.0E+1) +
            2 * area * density * N[2] * (1.0 / 3.0) * ((a[2] * (b[2] * b[2])) / 2.0E+1 + ((a[2] * a[2]) * b[2]) / 2.0E+1 + (a[2] * (c[2] * c[2])) / 2.0E+1 + ((a[2] * a[2]) * c[2]) / 2.0E+1 + (b[2] * (c[2] * c[2])) / 2.0E+1 + ((b[2] * b[2]) * c[2]) / 2.0E+1 + (a[2] * a[2] * a[2]) / 2.0E+1 + (b[2] * b[2] * b[2]) / 2.0E+1 + (c[2] * c[2] * c[2]) / 2.0E+1 + (a[2] * b[2] * c[2]) / 2.0E+1);

        I(1, 1) += 2 * area * density * N[0] * (1.0 / 3.0) * ((a[0] * (b[0] * b[0])) / 2.0E+1 + ((a[0] * a[0]) * b[0]) / 2.0E+1 + (a[0] * (c[0] * c[0])) / 2.0E+1 + ((a[0] * a[0]) * c[0]) / 2.0E+1 + (b[0] * (c[0] * c[0])) / 2.0E+1 + ((b[0] * b[0]) * c[0]) / 2.0E+1 + (a[0] * a[0] * a[0]) / 2.0E+1 + (b[0] * b[0] * b[0]) / 2.0E+1 + (c[0] * c[0] * c[0]) / 2.0E+1 + (a[0] * b[0] * c[0]) / 2.0E+1) +
            2 * area * density * N[2] * (1.0 / 3.0) * ((a[2] * (b[2] * b[2])) / 2.0E+1 + ((a[2] * a[2]) * b[2]) / 2.0E+1 + (a[2] * (c[2] * c[2])) / 2.0E+1 + ((a[2] * a[2]) * c[2]) / 2.0E+1 + (b[2] * (c[2] * c[2])) / 2.0E+1 + ((b[2] * b[2]) * c[2]) / 2.0E+1 + (a[2] * a[2] * a[2]) / 2.0E+1 + (b[2] * b[2] * b[2]) / 2.0E+1 + (c[2] * c[2] * c[2]) / 2.0E+1 + (a[2] * b[2] * c[2]) / 2.0E+1);

        I(2, 2) += 2 * area * density * N[1] * (1.0 / 3.0) * ((a[1] * (b[1] * b[1])) / 2.0E+1 + ((a[1] * a[1]) * b[1]) / 2.0E+1 + (a[1] * (c[1] * c[1])) / 2.0E+1 + ((a[1] * a[1]) * c[1]) / 2.0E+1 + (b[1] * (c[1] * c[1])) / 2.0E+1 + ((b[1] * b[1]) * c[1]) / 2.0E+1 + (a[1] * a[1] * a[1]) / 2.0E+1 + (b[1] * b[1] * b[1]) / 2.0E+1 + (c[1] * c[1] * c[1]) / 2.0E+1 + (a[1] * b[1] * c[1]) / 2.0E+1) +
            2 * area * density * N[0] * (1.0 / 3.0) * ((a[0] * (b[0] * b[0])) / 2.0E+1 + ((a[0] * a[0]) * b[0]) / 2.0E+1 + (a[0] * (c[0] * c[0])) / 2.0E+1 + ((a[0] * a[0]) * c[0]) / 2.0E+1 + (b[0] * (c[0] * c[0])) / 2.0E+1 + ((b[0] * b[0]) * c[0]) / 2.0E+1 + (a[0] * a[0] * a[0]) / 2.0E+1 + (b[0] * b[0] * b[0]) / 2.0E+1 + (c[0] * c[0] * c[0]) / 2.0E+1 + (a[0] * b[0] * c[0]) / 2.0E+1);

        I(1, 0) += -2 * area * density * N[0] * (1.0 / 2.0) * ((a[1] * (a[0] * a[0])) / 2.0E+1 + (a[1] * (b[0] * b[0])) / 6.0E+1 + (b[1] * (a[0] * a[0])) / 6.0E+1 + (a[1] * (c[0] * c[0])) / 6.0E+1 + (b[1] * (b[0] * b[0])) / 2.0E+1 + (c[1] * (a[0] * a[0])) / 6.0E+1 + (b[1] * (c[0] * c[0])) / 6.0E+1 + (c[1] * (b[0] * b[0])) / 6.0E+1 + (c[1] * (c[0] * c[0])) / 2.0E+1 + (a[1] * a[0] * b[0]) / 3.0E+1 + (a[1] * a[0] * c[0]) / 3.0E+1 + (b[1] * a[0] * b[0]) / 3.0E+1 + (a[1] * b[0] * c[0]) / 6.0E+1 + (b[1] * a[0] * c[0]) / 6.0E+1 + (c[1] * a[0] * b[0]) / 6.0E+1 + (b[1] * b[0] * c[0]) / 3.0E+1 + (c[1] * a[0] * c[0]) / 3.0E+1 + (c[1] * b[0] * c[0]) / 3.0E+1);
        I(0, 1) += -2 * area * density * N[0] * (1.0 / 2.0) * ((a[1] * (a[0] * a[0])) / 2.0E+1 + (a[1] * (b[0] * b[0])) / 6.0E+1 + (b[1] * (a[0] * a[0])) / 6.0E+1 + (a[1] * (c[0] * c[0])) / 6.0E+1 + (b[1] * (b[0] * b[0])) / 2.0E+1 + (c[1] * (a[0] * a[0])) / 6.0E+1 + (b[1] * (c[0] * c[0])) / 6.0E+1 + (c[1] * (b[0] * b[0])) / 6.0E+1 + (c[1] * (c[0] * c[0])) / 2.0E+1 + (a[1] * a[0] * b[0]) / 3.0E+1 + (a[1] * a[0] * c[0]) / 3.0E+1 + (b[1] * a[0] * b[0]) / 3.0E+1 + (a[1] * b[0] * c[0]) / 6.0E+1 + (b[1] * a[0] * c[0]) / 6.0E+1 + (c[1] * a[0] * b[0]) / 6.0E+1 + (b[1] * b[0] * c[0]) / 3.0E+1 + (c[1] * a[0] * c[0]) / 3.0E+1 + (c[1] * b[0] * c[0]) / 3.0E+1);

        I(2, 0) += -2 * area * density * N[0] * (1.0 / 2.0) * ((a[2] * (a[0] * a[0])) / 2.0E+1 + (a[2] * (b[0] * b[0])) / 6.0E+1 + (b[2] * (a[0] * a[0])) / 6.0E+1 + (a[2] * (c[0] * c[0])) / 6.0E+1 + (b[2] * (b[0] * b[0])) / 2.0E+1 + (c[2] * (a[0] * a[0])) / 6.0E+1 + (b[2] * (c[0] * c[0])) / 6.0E+1 + (c[2] * (b[0] * b[0])) / 6.0E+1 + (c[2] * (c[0] * c[0])) / 2.0E+1 + (a[2] * a[0] * b[0]) / 3.0E+1 + (a[2] * a[0] * c[0]) / 3.0E+1 + (b[2] * a[0] * b[0]) / 3.0E+1 + (a[2] * b[0] * c[0]) / 6.0E+1 + (b[2] * a[0] * c[0]) / 6.0E+1 + (c[2] * a[0] * b[0]) / 6.0E+1 + (b[2] * b[0] * c[0]) / 3.0E+1 + (c[2] * a[0] * c[0]) / 3.0E+1 + (c[2] * b[0] * c[0]) / 3.0E+1);
        I(0, 2) += -2 * area * density * N[0] * (1.0 / 2.0) * ((a[2] * (a[0] * a[0])) / 2.0E+1 + (a[2] * (b[0] * b[0])) / 6.0E+1 + (b[2] * (a[0] * a[0])) / 6.0E+1 + (a[2] * (c[0] * c[0])) / 6.0E+1 + (b[2] * (b[0] * b[0])) / 2.0E+1 + (c[2] * (a[0] * a[0])) / 6.0E+1 + (b[2] * (c[0] * c[0])) / 6.0E+1 + (c[2] * (b[0] * b[0])) / 6.0E+1 + (c[2] * (c[0] * c[0])) / 2.0E+1 + (a[2] * a[0] * b[0]) / 3.0E+1 + (a[2] * a[0] * c[0]) / 3.0E+1 + (b[2] * a[0] * b[0]) / 3.0E+1 + (a[2] * b[0] * c[0]) / 6.0E+1 + (b[2] * a[0] * c[0]) / 6.0E+1 + (c[2] * a[0] * b[0]) / 6.0E+1 + (b[2] * b[0] * c[0]) / 3.0E+1 + (c[2] * a[0] * c[0]) / 3.0E+1 + (c[2] * b[0] * c[0]) / 3.0E+1);

        I(2, 1) += -2 * area * density * N[1] * (1.0 / 2.0) * ((a[2] * (a[1] * a[1])) / 2.0E+1 + (a[2] * (b[1] * b[1])) / 6.0E+1 + (b[2] * (a[1] * a[1])) / 6.0E+1 + (a[2] * (c[1] * c[1])) / 6.0E+1 + (b[2] * (b[1] * b[1])) / 2.0E+1 + (c[2] * (a[1] * a[1])) / 6.0E+1 + (b[2] * (c[1] * c[1])) / 6.0E+1 + (c[2] * (b[1] * b[1])) / 6.0E+1 + (c[2] * (c[1] * c[1])) / 2.0E+1 + (a[2] * a[1] * b[1]) / 3.0E+1 + (a[2] * a[1] * c[1]) / 3.0E+1 + (b[2] * a[1] * b[1]) / 3.0E+1 + (a[2] * b[1] * c[1]) / 6.0E+1 + (b[2] * a[1] * c[1]) / 6.0E+1 + (c[2] * a[1] * b[1]) / 6.0E+1 + (b[2] * b[1] * c[1]) / 3.0E+1 + (c[2] * a[1] * c[1]) / 3.0E+1 + (c[2] * b[1] * c[1]) / 3.0E+1);
        I(1, 2) += -2 * area * density * N[1] * (1.0 / 2.0) * ((a[2] * (a[1] * a[1])) / 2.0E+1 + (a[2] * (b[1] * b[1])) / 6.0E+1 + (b[2] * (a[1] * a[1])) / 6.0E+1 + (a[2] * (c[1] * c[1])) / 6.0E+1 + (b[2] * (b[1] * b[1])) / 2.0E+1 + (c[2] * (a[1] * a[1])) / 6.0E+1 + (b[2] * (c[1] * c[1])) / 6.0E+1 + (c[2] * (b[1] * b[1])) / 6.0E+1 + (c[2] * (c[1] * c[1])) / 2.0E+1 + (a[2] * a[1] * b[1]) / 3.0E+1 + (a[2] * a[1] * c[1]) / 3.0E+1 + (b[2] * a[1] * b[1]) / 3.0E+1 + (a[2] * b[1] * c[1]) / 6.0E+1 + (b[2] * a[1] * c[1]) / 6.0E+1 + (c[2] * a[1] * b[1]) / 6.0E+1 + (b[2] * b[1] * c[1]) / 3.0E+1 + (c[2] * a[1] * c[1]) / 3.0E+1 + (c[2] * b[1] * c[1]) / 3.0E+1);

        temC[0] += 2 * area * density * N[0] * (1.0 / 2.0) * ((a[0] * b[0]) / 1.2E+1 + (a[0] * c[0]) / 1.2E+1 + (b[0] * c[0]) / 1.2E+1 + (a[0] * a[0]) / 1.2E+1 + (b[0] * b[0]) / 1.2E+1 + (c[0] * c[0]) / 1.2E+1);
        temC[1] += 2 * area * density * N[1] * (1.0 / 2.0) * ((a[1] * b[1]) / 1.2E+1 + (a[1] * c[1]) / 1.2E+1 + (b[1] * c[1]) / 1.2E+1 + (a[1] * a[1]) / 1.2E+1 + (b[1] * b[1]) / 1.2E+1 + (c[1] * c[1]) / 1.2E+1);
        temC[2] += 2 * area * density * N[2] * (1.0 / 2.0) * ((a[2] * b[2]) / 1.2E+1 + (a[2] * c[2]) / 1.2E+1 + (b[2] * c[2]) / 1.2E+1 + (a[2] * a[2]) / 1.2E+1 + (b[2] * b[2]) / 1.2E+1 + (c[2] * c[2]) / 1.2E+1);
    }

    center = (1.0 / mass) * temC;
    
}