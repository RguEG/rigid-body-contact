#include <collision_box_floor.h>
#include <iostream>

void collision_box_floor(std::vector<Eigen::Vector3d> &x, std::vector<Eigen::Vector3d> &n, std::vector<std::pair<int,int>> &objs,
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, unsigned int obj_id,
                         Eigen::Ref<Eigen::MatrixXd> V, 
                         Eigen::Ref<const Eigen::Vector3d> dir, Eigen::Ref<const Eigen::Vector3d> pos) {
    //Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  obj_id - id of object being checked for collision with the floor
//  V - the nx3 matrix of undeformed vertex positions
//  dir - the outward facing normal for the floor
//  pos - the world space position of a point on the floor plane
//Output:
//  x - world space, per-vertex collision points. Computed as any vertex that is on the "wrong side" of the floor plane
//  n - collision normals, one for each collision point. These point away from the floor. 
//  objs - Pairs of ids for objects involved in collisions. The first id is for the object, away from which the normal points. The second id 
//  is for the object towards which the normal points. The floor has an id of -1. 

    double dis;

    for (int i = 0; i < V.rows();i++) {
        Eigen::Vector3d someV;
        someV = R * V.row(i).transpose() +p;
        if (dir.transpose() * (someV - pos) <= 0) {
            x.push_back(V.row(i));
            n.push_back(dir);
            objs.push_back(std::pair<int,int> (obj_id,-1));
        }
    }
}