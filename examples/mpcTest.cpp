#include "MPC/MpcControl.hpp"
#include <iostream>
#include <casadi/casadi.hpp>

using std::cout;
using std::endl;
using casadi::DM;
using Eigen::Matrix;

int main(void){
    MpcAction mpc;
    mpc.RPY << 0.1, 0.1, 0;
    mpc.pPos << 0.02, 0, 0.79;
    mpc.dRPY << 0.2, 0.2, 0.;
    mpc.pVel << 0.4, 0.01, 0.02;
    mpc.desired_rpy << 0, 0, 0;
    mpc.desired_drpy << 0, 0, 0;
    mpc.desired_pos << 0, 0, 0.8;
    mpc.desired_vel << 0.5, 0, 0;
    mpc.comToFoot << 0.01, 0.01, -0.77;

    Matrix<double, 6, 1> sol_qpoases, sol_hpipmOcp, sol_hpipmDense, 
        sol_hpipmOcpParCond, sol_casadiOcp, sol_casadiDense,
        sol_osqpDense, sol_acadosOcp, sol_acadosOcpDense;
    sol_qpoases = mpc.solveMpcQpoasesDense();
    sol_hpipmOcp = mpc.solveMpcHpipmOcp();
    sol_hpipmOcpParCond = mpc.solveMpcHpipmOcpParDense();
    sol_hpipmDense = mpc.solveMpcHpipmDense();
    sol_casadiOcp = mpc.solveMpcCasadiOcp();
    sol_casadiDense = mpc.solveMpcCasadiDense();
    sol_osqpDense = mpc.solveMpcOsqpDense();
    sol_acadosOcp = mpc.solveMpcAcadosOcp();
    
    sol_acadosOcpDense = mpc.solveMpcAacadosDense();
    cout << "qpoases solve time  " << mpc.t_qpoaDen << " ms\n" << sol_qpoases.transpose() << 
        "\n\ncasadi-dense(solver:qrqp, similar to qpoases, active set) solve time  " << mpc.t_casDen << " ms\n"<< sol_casadiDense.transpose() <<
        "\n\nosqp-dense(ADMM) solve time  " <<mpc.t_osqpDen <<  " ms\n" << sol_osqpDense.transpose() <<
        "\n\ncasadi-ocp(solver:ipopt) solve  " << mpc.t_casOcp << " ms\n"<<sol_casadiOcp.transpose() << 
        "\n\nhpipm-dense solve time  " << mpc.t_hpDen << " ms\n" << sol_hpipmDense.transpose() <<
        "\n\nacados-dense(solver:hpipm, interior point) solve  " << mpc.t_acadosDense << " ms\n"<<sol_acadosOcpDense.transpose() <<
        "\n\nhpipm-ocp solve time  " << mpc.t_hpOcp << " ms\n" << sol_hpipmOcp.transpose() << 
        "\n\nhpipm-ocp-part-dense solve  "<<mpc.t_hpOcpPar<<  " ms\n"<< sol_hpipmOcpParCond.transpose() << 
        "\n\nacados-ocp solve  " << mpc.t_acados0cp << " ms\n"<<sol_acadosOcp.transpose() << endl;
}