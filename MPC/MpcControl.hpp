#ifndef _MPC_CONTROL_H_
#define _MPC_CONTROL_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <qpOASES.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Timer.hpp"

extern "C" {
#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_aux_ext_dep.h>

#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp_utils.h>
#include <hpipm_d_part_cond.h>

#include <hpipm_d_dense_qp.h>
#include <hpipm_d_dense_qp_sol.h>
#include <hpipm_d_dense_qp_ipm.h>
#include <hpipm_d_dense_qp_dim.h>
#include <hpipm_d_dense_qp_res.h>
#include <hpipm_timing.h>

#include <acados/utils/print.h>
#include <acados_c/ocp_qp_interface.h>

#include <acados/dense_qp/dense_qp_common.h>
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/utils/print.h>
#include <acados/utils/mem.h>
#include <acados_c/dense_qp_interface.h>
}

#include <casadi/casadi.hpp>
#include <OsqpEigen/OsqpEigen.h>

using S1 = casadi::Slice;
using casadi::MX;
using casadi::DM;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using std::cout;
using std::endl;

typedef Matrix<double,3,1> Vec3;
typedef Matrix<double,6,1> Vec6;
typedef Eigen::Matrix3d Mat3;

enum CoordinateAxis {X, Y, Z};
template <typename T>
using templateMat3 = typename Eigen::Matrix<T, 3, 3>;
template <typename T>
templateMat3<T> coordinateRotation(CoordinateAxis axis, const T& theta){
  T s = std::sin(theta);
  T c = std::cos(theta);
  templateMat3<T> R;
  if(axis == CoordinateAxis::X){
    R << 1, 0, 0, 0, c, -s, 0, s, c;
  }
  else if(axis == CoordinateAxis::Y){
    R << c, 0, s, 0, 1, 0, -s, 0, c;
  }
  else if(axis == CoordinateAxis::Z){
    R << c, -s, 0, s, c, 0, 0, 0, 1;
  }
  
  return R;
}

#define MPC_INFINITY 5e10
#define MPC_CT 0.03
#define HORIZON 15
#define QP_MAX_TIME 0.02
class MpcAction
{
    public:
        MpcAction();
        ~MpcAction();
        Mat3 vecCross(Vec3 v);
        void matrix_to_real(qpOASES::real_t* dst, Matrix<double,Dynamic,Dynamic> src, int rows, int cols);
        void matrix_to_doubleArray(double* dst, Matrix<double,Dynamic,Dynamic> src, int rows, int cols);
        void matrix_to_intArray(int* dst, Matrix<int,Dynamic,Dynamic> src, int rows, int cols);
        void real_to_matrix(Matrix<double,Dynamic,1>& dst, qpOASES::real_t* src, int n_items);

        void getInitState();
        void getConStateMatrix();//continuous state space model
        void getDisStateMatrix();//discrete state space model
        void getQpMatrix();
        void setRefTrajectory();
        void getMpcMatrix();

        Vec6 solveMpcQpoasesDense();
        Vec6 solveMpcHpipmOcp();
        Vec6 solveMpcHpipmOcpParDense();
        Vec6 solveMpcHpipmDense();

        Vec6 solveMpcCasadiOcp();
        Vec6 solveMpcCasadiDense();

        //for symbolic computation of casadi
        MX crossProductMx(MX vec);
        MX EulerZYX2RotMatMx(MX rpy);
        MX getAngVelTransMatrixMx(MX rpy);

        Vec6 solveMpcOsqpDense();

        Vec6 solveMpcAcadosOcp();
        Vec6 solveMpcAacadosDense();

    private:
        Matrix<double,3,1> rpy, xyz, wI, vI, r;
        Matrix<double, 6, 1> desiredVel;
        Matrix<double,3,3> BodyWorld, angVelTrans, invAngVelTrans, Ib, Iw, invIw;
        Matrix<double,12,1> X0;
        Matrix<double,12,1> C, Cdt;
        Matrix<double,12,6> B, Bdt;
        Matrix<double,12,12> A, Adt;
        Matrix<double,12,12> L;
        Matrix<double,12,1> stateWeight;
        Matrix<double,6,1> controlWeight;
        Matrix<double,6,6> K;
        Matrix<double,8,6> f_block;
        Matrix<double,8,1> inequa_lb, inequa_ub;

        Eigen::DiagonalMatrix<double, Dynamic> Lqp;
        Eigen::DiagonalMatrix<double, Dynamic> Kqp;
        Matrix<double,Dynamic, 12> Aqp;
        Matrix<double,Dynamic, Dynamic> Bqp;
        Matrix<double,Dynamic, 1> Cqp;
        Matrix<double,Dynamic,1> Xdqp;
        Matrix<double,Dynamic, 1> Lbqp, Ubqp;
        Matrix<double,Dynamic, Dynamic> fmat;
        Matrix<double,Dynamic, Dynamic> Hqp;
        Matrix<double,Dynamic, 1> Gqp;

        Matrix<double,Dynamic,1> U_sol;
        double solve_time;
        
        qpOASES::real_t* H_qpoases;
        qpOASES::real_t* g_qpoases;
        qpOASES::real_t* A_qpoases;
        qpOASES::real_t* lb_qpoases;
        qpOASES::real_t* ub_qpoases;

        qpOASES::real_t* q_soln;
    public:
        const double m = 10; 
        const double  Ixx = 1;
        const double  Ixy =  0;
        const double  Ixz =  0;
        const double  Iyy = 1;
        const double  Iyz =  0;
        const double  Izz = 1;

        const int horizon = HORIZON;
        const double miu = 0.1;
        const double mu = 1.0 / miu;
        const double fzmin = 0;
        const double fzmax = 300;
        const double torXmin = -50;
        const double torXmax = 50;
        const double torYmin = -50;
        const double torYmax = 50;
        const double torZmin = -50;
        const double torZmax = 50;
        const double CT = MPC_CT;
        qpOASES::real_t cputime = QP_MAX_TIME;
        Matrix<double,3,1> g;
        Matrix<double,3,1> RPY, pPos, pVel, dRPY, comToFoot;
        Matrix<double,3,1> desired_rpy,desired_drpy,desired_pos,desired_vel;
        
        double t_qpoaDen, t_hpDen, t_hpOcp, t_hpOcpPar, t_casDen, t_osqpDen, t_casOcp,
                t_acados0cp, t_acadosDense;
};



#endif