#include "MpcControl.hpp"


MpcAction::MpcAction() 
{
    stateWeight.setZero();
    controlWeight.setZero();
    
    stateWeight << 1, 1, 1,   1, 1, 1,  1, 1, 1,   1, 1, 1;
    controlWeight << 1, 1, 1,   1, 1, 1;


    RPY.setZero();
    dRPY.setZero();
    pPos.setZero();
    pVel.setZero();
    comToFoot.setZero();
    Lqp.resize(12*HORIZON);
    Kqp.resize(6*HORIZON);
    Lqp.setZero();
    Kqp.setZero();
    Aqp = MatrixXd::Zero(12*HORIZON, 12);
    Bqp = MatrixXd::Zero(12*HORIZON, 6*HORIZON);
    Cqp = MatrixXd::Zero(12*HORIZON, 1);
    Xdqp = MatrixXd::Zero(12*HORIZON, 1);
    Lbqp = MatrixXd::Zero(8*HORIZON, 1);
    Ubqp = MatrixXd::Zero(8*HORIZON, 1);
    fmat = MatrixXd::Zero(8*HORIZON, 6*HORIZON);
    Hqp = MatrixXd::Zero(6*HORIZON, 6*HORIZON);
    Gqp = MatrixXd::Zero(6*HORIZON, 1);
    g << 0, 0, -9.81;
    Ib << Ixx,Ixy,Ixz,
         Ixy,Iyy,Iyz,
         Ixz,Iyz,Izz;
    BodyWorld << coordinateRotation<double>(CoordinateAxis::Z, rpy[2]) * 
                 coordinateRotation<double>(CoordinateAxis::Y, rpy[1]) *
                 coordinateRotation<double>(CoordinateAxis::X, rpy[0]);

    angVelTrans << cos(rpy[1])*cos(rpy[2]), -sin(rpy[2]), 0,
                    cos(rpy[1])*sin(rpy[2]), cos(rpy[2]), 0,
                    -sin(rpy[1]), 0, 1;
    invAngVelTrans = angVelTrans.inverse();
    Iw = BodyWorld * Ib * BodyWorld.transpose(); 
    invIw = Iw.inverse();

    Lqp.diagonal() = stateWeight.replicate(horizon, 1);
    Kqp.diagonal() = controlWeight.replicate(horizon, 1);
    
    f_block <<   mu, 0, 1, 0, 0, 0,
                -mu, 0, 1, 0, 0, 0,
                0, mu, 1, 0, 0, 0,
                0, -mu, 1, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1;
    inequa_lb << 0, 0, 0, 0, fzmin, torXmin, torYmin, torZmin;
    inequa_ub << MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, fzmax, torXmax, torYmax, torZmax;

    fmat.setZero();
    for(int i = 0; i < horizon; i++){
        fmat.block(i*8, i*6, 8, 6) = f_block;
    }
    Lbqp.setZero();
    Ubqp.setZero();
    for(int i = 0; i < horizon; i++){
        Lbqp.block(i*8, 0, 8, 1) = inequa_lb;
        Ubqp.block(i*8, 0, 8, 1) = inequa_ub;
    }

    H_qpoases = (qpOASES::real_t*)malloc(6*horizon*6*horizon*sizeof(qpOASES::real_t));
    g_qpoases = (qpOASES::real_t*)malloc(6*horizon*sizeof(qpOASES::real_t));
    A_qpoases = (qpOASES::real_t*)malloc(8*horizon*6*horizon*sizeof(qpOASES::real_t));
    lb_qpoases = (qpOASES::real_t*)malloc(8*horizon*sizeof(qpOASES::real_t));
    ub_qpoases = (qpOASES::real_t*)malloc(8*horizon*sizeof(qpOASES::real_t));
    q_soln = (qpOASES::real_t*)malloc(6*horizon*sizeof(qpOASES::real_t));

    matrix_to_real(A_qpoases, fmat, 8*horizon, 6*horizon);
    matrix_to_real(lb_qpoases, Lbqp, 8*horizon, 1);
    matrix_to_real(ub_qpoases, Ubqp, 8*horizon, 1);
}

MpcAction::~MpcAction()
{
}

Matrix<double,3,3> MpcAction::vecCross(Matrix<double,3,1> v)
{
    Matrix<double,3,3> R;
    R << 0,-v[2],v[1],
         v[2],0,-v[0],
         -v[1],v[0],0;
    return R;
}
void MpcAction::matrix_to_real(qpOASES::real_t* dst, Matrix<double,Dynamic,Dynamic> src, int rows, int cols)
{
  int a = 0;
  for(int r = 0; r < rows; r++)
  {
    for(int c = 0; c < cols; c++)
    {
      dst[a] = (double)src(r,c);
      a++;
    }
  }
}
void MpcAction::matrix_to_doubleArray(double* dst, Matrix<double,Dynamic,Dynamic> src, int rows, int cols)
{
  int a = 0;
  for(int r = 0; r < rows; r++)
  {
    for(int c = 0; c < cols; c++)
    {
      dst[a] = (double)src(r,c);
      a++;
    }
  }
}
void MpcAction::matrix_to_intArray(int* dst, Matrix<int,Dynamic,Dynamic> src, int rows, int cols)
{
  int a = 0;
  for(int r = 0; r < rows; r++)
  {
    for(int c = 0; c < cols; c++)
    {
      dst[a] = (int)src(r,c);
      a++;
    }
  }
}

void MpcAction::real_to_matrix(Matrix<double,Dynamic,1>& dst, qpOASES::real_t* src, int n_items)
{
    dst.resize(n_items,1);
    dst.setZero();
    for(int i=0;i<n_items;i++)
        dst(i) = (double)src[i];
}

void MpcAction::getInitState()
{
    rpy = RPY;
    xyz = pPos;
    vI = pVel;
    angVelTrans << cos(rpy[1]) * cos(rpy[2]), -sin(rpy[2]), 0,
        cos(rpy[1])* sin(rpy[2]), cos(rpy[2]), 0,
        -sin(rpy[1]), 0, 1;
    wI = angVelTrans * dRPY;
    r = comToFoot;
    BodyWorld << coordinateRotation<double>(CoordinateAxis::Z, rpy[2]) * 
                 coordinateRotation<double>(CoordinateAxis::Y, rpy[1]) *
                 coordinateRotation<double>(CoordinateAxis::X, rpy[0]);
    invAngVelTrans = angVelTrans.inverse();
    Iw = BodyWorld * Ib * BodyWorld.transpose();
    invIw = Iw.inverse();

    X0 << rpy, xyz, Iw * wI, m * vI;
}


void MpcAction::getConStateMatrix()
{

    A  <<   MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), angVelTrans.inverse()*invIw, MatrixXd::Zero(3,3),
            MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), MatrixXd::Identity(3,3)/m,
            MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), MatrixXd::Zero(3,3),
            MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), MatrixXd::Zero(3,3), MatrixXd::Zero(3,3);
    B << MatrixXd::Zero(3,3), MatrixXd::Zero(3,3),
            MatrixXd::Zero(3,3), MatrixXd::Zero(3,3),
            vecCross(r), MatrixXd::Identity(3,3),
            MatrixXd::Identity(3,3), MatrixXd::Zero(3,3);
    C << MatrixXd::Zero(3,1), MatrixXd::Zero(3,1), MatrixXd::Zero(3,1), m*g;

//    cout<<"A:"<<A<<endl;
//    cout<<"B:"<<B.transpose()<<endl;
//    cout<<"C:"<<C.transpose()<<endl;
}

void MpcAction::getDisStateMatrix()
{
    Adt = MatrixXd::Identity(12,12) + CT * A;
    Bdt = CT * B;
    Cdt = CT * C;
    Cdt(5,0) = 0.5*(-9.81)*CT*CT;
}

void MpcAction::getQpMatrix()
{
    Aqp.setZero();
    Bqp.setZero();
    Cqp.setZero();

    Matrix<double, 12, 12> powerMatsA[HORIZON+1];
    Matrix<double, 12, 1> powerMatsC[HORIZON+1];
    powerMatsA[0].setIdentity();
    powerMatsC[0].setZero();
    for(int i=1; i<horizon + 1; i++){
        powerMatsA[i] = Adt * powerMatsA[i-1];
        powerMatsC[i] = powerMatsA[i-1] * Cdt + powerMatsC[i-1];
    }

    for(int r=0;r<horizon;r++){
        Aqp.block(12*r, 0, 12, 12) = powerMatsA[r+1];
        Cqp.block(12*r, 0, 12, 1) = powerMatsC[r+1];
        for(int c=0;c<horizon;c++){
            if(r >= c){
                int a = r-c;
                Bqp.block(12*r, 6*c, 12, 6) = powerMatsA[a] * Bdt;
            }
        }
    }
}

void MpcAction::setRefTrajectory(){
    Xdqp.setZero();

    Matrix<double, 12, 1> powerMatsTraj[HORIZON];
    Matrix<double, 12, 1> velCmd;
    Matrix<double, 3, 1> desiredLcom, desiredPcom;

    desiredVel << desired_drpy, desired_vel;//rpy + xyz
    velCmd << desiredVel, MatrixXd::Zero(6, 1);
    desiredLcom = Iw * angVelTrans * desiredVel.head(3);
    desiredPcom = m * desiredVel.tail(3);
    powerMatsTraj[0] << 0, 0, rpy[2], xyz[0], xyz[1], desired_pos[2], desiredLcom, desiredPcom;
    for (int i = 1; i < horizon; ++i) {
        powerMatsTraj[i] = powerMatsTraj[i - 1] + velCmd * CT;
    }
    for (int r = 0; r < horizon; ++r) {
        Xdqp.block(12 * r, 0, 12, 1) = powerMatsTraj[r];
    }
}

void MpcAction::getMpcMatrix()
{
   //Timer t1;
    getInitState();
    getConStateMatrix();
    getDisStateMatrix();
    getQpMatrix();
    setRefTrajectory();
    Hqp.setZero();
    Gqp.setZero();
    
    Hqp = Bqp.transpose()*Lqp*Bqp;
    Hqp.diagonal() += Kqp.diagonal();
    Gqp = Bqp.transpose()*Lqp*(Aqp*X0 + Cqp - Xdqp);
    
    //cout<<"H && G time: "<< t2.getMs() << " ms" <<endl;
}

Vec6 MpcAction::solveMpcQpoasesDense()
{
    Timer cal_time;
    getMpcMatrix();
    //Timer trans;
     matrix_to_real(H_qpoases, Hqp, 6*horizon, 6*horizon);
     matrix_to_real(g_qpoases, Gqp, 6*horizon, 1);
    //cout<<"matrix_to_real time: "<< trans.getMs() << " ms" <<endl;
    int num_constraints = fmat.rows();
    int num_variables = fmat.cols();
    qpOASES::int_t nWSR = 200;//1000
    qpOASES::QProblem problem_red(num_variables, num_constraints);
    qpOASES::Options op;
    op.setToMPC();
    //op.setToReliable();
    op.printLevel = qpOASES::PL_NONE;
    problem_red.setOptions(op);
//    Timer sol;
    int rval2 = problem_red.init(H_qpoases, g_qpoases, A_qpoases, NULL, NULL, lb_qpoases, ub_qpoases, nWSR, &cputime);
    //int rval2 = problem_red.init(H_qpoases, g_qpoases, A_qpoases, NULL, NULL, NULL, NULL, nWSR, &cputime);
    (void) rval2;
     if(rval2 != qpOASES::SUCCESSFUL_RETURN){
         cout<<"failed to init QP!"<<endl;
          if(rval2 == qpOASES::RET_MAX_NWSR_REACHED) cout<<"RET_MAX_NWSR_REACHED"<<endl;
          else if (rval2 == qpOASES::RET_INIT_FAILED) cout<<"RET_INIT_FAILED"<<endl;
          else if(rval2 == qpOASES::RET_INIT_FAILED_HOTSTART) cout<<"QP Cannot Be Solved!"<<endl;
          cout<<"rval2: "<<rval2<<endl;
     }
    int rval = problem_red.getPrimalSolution(q_soln);
    if(rval != qpOASES::SUCCESSFUL_RETURN)
        cout<<"failed to solve!"<<endl;
    //cout<<"solve time: "<< sol.getMs() << " ms" <<endl;
    // for(int i= 0;i<HORIZON;++i)
    //     std::cout << q_soln[i*6 + 1] << std::endl;
    real_to_matrix(U_sol, q_soln, 6);
    solve_time = cal_time.getMs();
    t_qpoaDen = cal_time.getMs();
    //std::cout << "qpoases solve_time   " <<solve_time<<std::endl;
    return U_sol;//只取第一步的解
}

Vec6 MpcAction::solveMpcHpipmOcp(){
    Timer cal_time;

    getInitState();
    getConStateMatrix();
    getDisStateMatrix();
    setRefTrajectory();

    static int N = HORIZON;
    std::vector<int> nx_(N + 1, 12);
    std::vector<int> nu_(N + 1, 6);
    std::fill_n(nu_.begin() + N, 1, 0);
    // for(std::vector<int>::iterator it = nu_.begin(); it != nu_.end(); ++it)
    //     std::cout << *it << std::endl;
    std::vector<int> nbx_(N + 1, 0);
    nbx_.at(0) = 12;
    //std::fill_n(nbx_.begin() + N, 1, 12);

    std::vector<int> nbu_(N + 1, 4);
    std::fill_n(nbu_.begin() + N, 1, 0);

    std::vector<int> ng_(N + 1, 4);
    std::fill_n(ng_.begin() + N, 1, 0);

    std::vector<int> nsbx_(N + 1, 0);
    std::vector<int> nsbu_(N + 1, 0);
    std::vector<int> nsg_(N + 1, 0);
    
    //ocp qp dim
    int dim_size = d_ocp_qp_dim_memsize(N);
    void *dim_mem = malloc(dim_size);
    if(!dim_mem) std::cout << "malloc dim memory failed" << std::endl;

    struct d_ocp_qp_dim dim;
	d_ocp_qp_dim_create(N, &dim, dim_mem);
    d_ocp_qp_dim_set_all(nx_.data(), nu_.data(), nbx_.data(), nbu_.data(), ng_.data(), nsbx_.data(), nsbu_.data(), nsg_.data(), &dim);
    //ocp qp
    int qp_size = d_ocp_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);
    if(!qp_mem) std::cout << "malloc qp memory failed" << std::endl;

    struct d_ocp_qp qp;
	d_ocp_qp_create(&dim, &qp, qp_mem);

    static double hA_[12 * 12];
    static double hB_[12 * 6];
    static double hb_[12 * 1];

    static double Q_[12 * 12];
    static double S_[6 * 12] = {0};
    static double R_[6 * 6];
    static double q_[12 * 1];
    static double r_[6 * 1] = {0};
    
    static int idxbu_[4 * 1];
    static double lbu_[4 * 1];
    static double ubu_[4 * 1];

    static int idxbx0_[12 * 1];
    static double x0_[12 * 1];

    static double C_[4 * 12] = {0};
    static double D_[4 * 6];
    static double lg_[4 * 1];
    static double ug_[4 * 1];

    Eigen::Matrix<int, 4, 1> Jbu;
    Eigen::Matrix<int, 12, 1> Jbx;
    Eigen::Matrix<double, 4, 6> Dn;
    Eigen::Matrix<double, 4, 1> Ulow, Uup;
    Eigen::Matrix<double, 12, 1> X_0;
    Eigen::Matrix<double, 4, 1> dlow, dup;
    Eigen::Matrix<double, 12, 1> qn;

    Eigen::DiagonalMatrix<double, 12> Q_temp;Q_temp.setZero();
    Eigen::DiagonalMatrix<double, 6> R_temp;R_temp.setZero();
    Eigen::Matrix<double, 12, 12> Q;
    Eigen::Matrix<double, 6, 6> R;
    Q_temp.diagonal() = stateWeight.replicate(1, 1);
    R_temp.diagonal() = controlWeight.replicate(1, 1);
    
    Q = Q_temp.diagonal().asDiagonal();
    R = R_temp.diagonal().asDiagonal();
    //std::cout << Q <<"\n\n\n" << R << std::endl;
    qn = -2 * Q * Xdqp.head(12);
    Jbu << 2, 3, 4, 5; 
    Jbx << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;

    Dn << mu, 0, 1, 0, 0, 0,
        -mu, 0, 1, 0, 0, 0,
        0, mu, 1, 0, 0, 0,
        0, -mu, 1, 0, 0, 0;
    Ulow << fzmin, torXmin, torYmin, torZmin;
    Uup << fzmax, torXmax, torYmax, torZmax;
    X_0 = X0;
    dlow << 0, 0, 0, 0;
    dup << MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, MPC_INFINITY;
    //std::cout << Adt << std::endl;
    matrix_to_doubleArray(hA_, Adt, 12, 12);
    matrix_to_doubleArray(hB_, Bdt, 12, 6);
    matrix_to_doubleArray(hb_, Cdt, 12, 1);
    matrix_to_doubleArray(Q_, Q, 12, 12);
    matrix_to_doubleArray(R_, R, 6, 6);
    matrix_to_doubleArray(q_, qn, 12, 1);
    matrix_to_intArray(idxbu_, Jbu, 4, 1);
    matrix_to_doubleArray(lbu_, Ulow, 4, 1);
    matrix_to_doubleArray(ubu_, Uup, 4, 1);
    matrix_to_doubleArray(D_, Dn, 4, 6);
    matrix_to_doubleArray(lg_, dlow, 4, 1);
    matrix_to_doubleArray(ug_, dup, 4, 1);
    matrix_to_intArray(idxbx0_, Jbx, 12, 1);
    matrix_to_doubleArray(x0_, X_0, 12, 1);

    std::vector<double*> AA(N, hA_);
    std::vector<double*> BB(N, hB_);
    std::vector<double*> bb(N, hb_);
    std::vector<double*> hQ(N + 1, Q_);
    std::vector<double*> hR(N + 1, R_);
    std::vector<double*> hq(N + 1, q_);
    std::vector<double*> hS(N + 1, S_);
    std::vector<double*> hr(N + 1, r_);

    std::vector<int*> hidxbx(N + 1, idxbx0_);
    std::vector<double*> hlbx(N + 1, x0_);
    std::vector<double*> hubx(N + 1, x0_);
    std::fill_n(hidxbx.begin() + 1, N, nullptr);
    std::fill_n(hlbx.begin() + 1, N, nullptr);
    std::fill_n(hubx.begin() + 1, N, nullptr);

    std::vector<int*> hidxbu(N + 1, idxbu_);
    std::vector<double*> hlbu(N + 1, lbu_);
    std::vector<double*> hubu(N + 1, ubu_);

    std::vector<double*> hC(N + 1, C_);
    std::vector<double*> hD(N + 1, D_);
    std::vector<double*> hlg(N + 1, lg_);
    std::vector<double*> hug(N + 1, ug_);

    std::vector<double*> hZl(N + 1, nullptr);
    std::vector<double*> hZu(N + 1, nullptr);
    std::vector<double*> hzl(N + 1, nullptr);
    std::vector<double*> hzu(N + 1, nullptr);
    std::vector<int*> hidxs(N + 1, nullptr);
    std::vector<double*> hlls(N + 1, nullptr);
    std::vector<double*> hlus(N + 1, nullptr);

    d_ocp_qp_set_all(AA.data(), BB.data(), bb.data(),
        hQ.data(), hS.data(), hR.data(), hq.data(), hr.data(),
        hidxbx.data(), hlbx.data(), hubx.data(),
        hidxbu.data(), hlbu.data(), hubu.data(),
        hC.data(), hD.data(), hlg.data(), hug.data(),
        hZl.data(), hZu.data(), hzl.data(), hzu.data(),
        hidxs.data(), hlls.data(), hlus.data(),
        &qp
        );
    //ocp sol
    int qp_sol_size = d_ocp_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    if(!qp_sol_mem) std::cout << "malloc qp solve memory failed" << std::endl;

	struct d_ocp_qp_sol qp_sol;
	d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);
    //ipm arg
    int ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);

    hpipm_mode mode = hpipm_mode::SPEED;//BALANCE;//ROBUST;//SPEED;
    int iter_max = 200;
    double alpha_min = 1e-4;
    double mu0 = 10;
    double tol_stat = 1e-5;
    double tol_eq = 1e-5;
    double tol_ineq = 1e-5;
    double tol_comp = 1e-5;
    double reg_prim = 1e-8;
    int warm_start = 0;
    int pred_corr = 1;
    int ric_alg = 0;
    int split_step = 1;
	d_ocp_qp_ipm_arg_set_default(mode, &arg);
	d_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
	d_ocp_qp_ipm_arg_set_alpha_min(&alpha_min, &arg);
	d_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_ocp_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
	d_ocp_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
	d_ocp_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
	d_ocp_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
	d_ocp_qp_ipm_arg_set_reg_prim(&reg_prim, &arg);
	d_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);
	d_ocp_qp_ipm_arg_set_pred_corr(&pred_corr, &arg);
	d_ocp_qp_ipm_arg_set_ric_alg(&ric_alg, &arg);
	d_ocp_qp_ipm_arg_set_split_step(&split_step, &arg);
    //ipm workspace
    int ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_ws workspace;
	d_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);
    //ipm solver
    int nrep = 1;
    int hpipm_status;
    for(int rep = 0; rep < nrep; rep++){
		// call solver
		d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
		d_ocp_qp_ipm_get_status(&workspace, &hpipm_status);
	}
    if(hpipm_status == 0)
		{
        std::cout << "hpipm ocp qp solved\n";
		}
	else if(hpipm_status==1)
		{
        std::cout << "hpipm ocp qp failed! Maximum number of iterations reached\n";
		}
	else if(hpipm_status==2)
		{
        std::cout << "hpipm ocp qp failed! Minimum step lenght reached\n";
		}
	else if(hpipm_status==3)
		{
        std::cout << "hpipm ocp qp failed! NaN in computations\n";
		}
    
    std::vector<Eigen::Matrix<double, 6, 1>> u_sol(N);
    std::vector<Eigen::Matrix<double, 12, 1>> x_sol(N + 1);
    for(int i = 0; i < N; ++i){
        d_ocp_qp_sol_get_u(i, &qp_sol, u_sol[i].data());
    }
    for(int i = 0; i < N + 1; ++i){
        d_ocp_qp_sol_get_x(i, &qp_sol, x_sol[i].data());
    }
    Matrix<double,6,1> Usol_hpipm;
    Usol_hpipm = u_sol[0];
    // for(int i= 0;i<N;++i)
    //     std::cout << u_sol[i][1] << std::endl;
    t_hpOcp = cal_time.getMs();
    //std::cout << "hpipm ocp solve time:  " << cal_time.getMs() << std::endl;
    free(dim_mem);
    free(qp_mem);
	free(qp_sol_mem);
	free(ipm_arg_mem);
	free(ipm_mem);

    return Usol_hpipm;
}
Vec6 MpcAction::solveMpcHpipmDense(){
    const int N = HORIZON;

    int nv = 6 * N;
    int ne = 0;
    int nb = 4 * N;
    int ng = 4 * N;
    int ns = 0;
    int nsb = 0;
    int nsg = 0;

    double hH[6 * N * 6 * N];
    double hg[6 * N];

    double hA[6 * N * 6 * N] = {0};
    double hb[6 * N] = {0};

    int hidv[4 * N];
    double hlv[4 * N];
    double huv[4 * N];
    
    double hC[4 * N * 6 * N];
    double hlg[4 * N];
    double hug[4 * N];
    
    Eigen::Matrix<int, 4, 1> idv;
    Eigen::Matrix<int, 4 * N, 1> idv_temp;
    idv << 2, 3, 4, 5;
    idv_temp = idv.replicate(N, 1);
    //std::cout << idv_temp << std::endl;
    Eigen::Matrix<double, 4, 1> lv, uv, lg, ug;
    Eigen::Matrix<double, 4, 6> hc_;
    Eigen::Matrix<double, 4 * N, 1> lv_temp, uv_temp, lg_temp, ug_temp;
    Eigen::Matrix<double, 4 * N, 6 * N> hc_temp;
    hc_temp.setZero();

    lv << fzmin, torXmin, torYmin, torZmin;
    uv << fzmax, torXmax, torYmax, torZmax;
    lg << 0, 0, 0, 0;
    ug << MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, MPC_INFINITY;
    hc_ << mu, 0, 1, 0, 0, 0,
        -mu, 0, 1, 0, 0, 0,
        0, mu, 1, 0, 0, 0,
        0, -mu, 1, 0, 0, 0;

    lv_temp = lv.replicate(N, 1); 
    uv_temp = uv.replicate(N, 1);
    lg_temp = lg.replicate(N, 1);
    ug_temp = ug.replicate(N, 1);

    for(int i=0;i<N;++i)
        hc_temp.block(4 * i, 6 * i, 4, 6) = hc_;
    //std::cout << hc_temp << std::endl;
    matrix_to_intArray(hidv, idv_temp, 4 * N, 1);
    matrix_to_doubleArray(hlv, lv_temp, 4 * N, 1);
    matrix_to_doubleArray(huv, uv_temp, 4 * N, 1);
    matrix_to_doubleArray(hlg, lg_temp, 4 * N, 1);
    matrix_to_doubleArray(hug, ug_temp, 4 * N, 1);
    matrix_to_doubleArray(hC, hc_temp, 4 * N, 6 * N);

    int dense_qp_dim_size = d_dense_qp_dim_memsize();
    void *qp_dim_mem = malloc(dense_qp_dim_size);
    if(!qp_dim_mem) std::cout << "dense qp dim memory malloc failed" << std::endl;

	struct d_dense_qp_dim qp_dim;
	d_dense_qp_dim_create(&qp_dim, qp_dim_mem);
	d_dense_qp_dim_set_all(nv, ne, nb, ng, nsb, nsg, &qp_dim);

    int qp_size = d_dense_qp_memsize(&qp_dim);
    void *qp_mem = malloc(qp_size);
    if(!qp_mem) std::cout << "dense qp memory malloc failed" << std::endl;

	struct d_dense_qp qp;
	d_dense_qp_create(&qp_dim, &qp, qp_mem);

    Timer solve_time;
    getMpcMatrix();
    matrix_to_doubleArray(hH, Hqp, 6 * N, 6 * N);
    matrix_to_doubleArray(hg, Gqp, 6 * N, 1);
    
    d_dense_qp_set_all(hH, hg, 
        nullptr, nullptr, 
        hidv, hlv, huv, 
        hC, hlg, hug, 
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 
        &qp
    );
    int qp_sol_size = d_dense_qp_sol_memsize(&qp_dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    if(!qp_sol_mem) std::cout << "dense qp solve dim memory malloc failed" << std::endl;

	struct d_dense_qp_sol qp_sol;
	d_dense_qp_sol_create(&qp_dim, &qp_sol, qp_sol_mem);

    int ipm_arg_size = d_dense_qp_ipm_arg_memsize(&qp_dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    if(!ipm_arg_mem) std::cout << "dense qp argument dim memory malloc failed" << std::endl;

	struct d_dense_qp_ipm_arg arg;
	d_dense_qp_ipm_arg_create(&qp_dim, &arg, ipm_arg_mem);
    hpipm_mode mode = hpipm_mode::SPEED;//BALANCE;//ROBUST;//SPEED;
    //https://github1s.com/giaf/hpipm/blob/HEAD/include/hpipm_d_dense_qp_ipm.h#L62
    int iter_max = 20; //25;
	double mu0 = 0.1;
    double alpha_min = 1e-12;
	int comp_res_exit = 0;
	double tol_stat = 1e-4;
	double tol_eq = 1e-4;
	double tol_ineq = 1e-4;
	double tol_comp = 1e-4;
	int kkt_fact_alg = 0;
	int remove_lin_dep_eq = 0;
	double tau_min = 1e-3;
    int split_step = 0;
    
    d_dense_qp_ipm_arg_set_default(mode, &arg);
    d_dense_qp_ipm_arg_set_alpha_min(&alpha_min, &arg);
    d_dense_qp_ipm_arg_set_iter_max(&iter_max, &arg);
	d_dense_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_dense_qp_ipm_arg_set_comp_res_exit(&comp_res_exit, &arg);
	d_dense_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
	d_dense_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
	d_dense_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
	d_dense_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
	d_dense_qp_ipm_arg_set_kkt_fact_alg(&kkt_fact_alg, &arg);
	d_dense_qp_ipm_arg_set_remove_lin_dep_eq(&remove_lin_dep_eq, &arg); 
    d_dense_qp_ipm_arg_set_split_step(&split_step, &arg);

    int ipm_size = d_dense_qp_ipm_ws_memsize(&qp_dim, &arg);
    void *ipm_mem = malloc(ipm_size);
	struct d_dense_qp_ipm_ws workspace;
	d_dense_qp_ipm_ws_create(&qp_dim, &arg, &workspace, ipm_mem);
    
    int hpipm_dense_status, nrep = 1;
    for(int rep=0; rep<nrep; rep++){
		d_dense_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
		d_dense_qp_ipm_get_status(&workspace, &hpipm_dense_status);
	}
    if(hpipm_dense_status == 0)
		{
        std::cout << "hpipm dense QP solved\n";
		}
	else if(hpipm_dense_status==1)
		{
        std::cout << "hpipm dense QP failed! Maximum number of iterations reached\n";
		}
	else if(hpipm_dense_status==2)
		{
        std::cout << "hpipm dense QP failed! Minimum step lenght reached\n";
		}
	else if(hpipm_dense_status==3)
		{
        std::cout << "hpipm dense QP failed! NaN in computations\n";
		}
    
    std::vector<double> u_sol(6 * N);
    d_dense_qp_sol_get_v(&qp_sol, u_sol.data());
    // for(int i=0;i<N;++i)
    //     std::cout << u_sol[i*6 + 2] << std::endl;
    Eigen::Matrix<double,6,1> Usol_hpipmDense;
    for(int i=0;i<6;++i)
        Usol_hpipmDense[i] = u_sol[i];

    t_hpDen = solve_time.getMs();
    //std::cout << "hpipm dense qp solve time:  " << solve_time.getMs() << std::endl;

    return Usol_hpipmDense;
}

Vec6 MpcAction::solveMpcHpipmOcpParDense(){
    static int N = HORIZON;
    std::vector<int> nx_(N + 1, 12);
    std::vector<int> nu_(N + 1, 6);
    std::fill_n(nu_.begin() + N, 1, 0);
    // for(std::vector<int>::iterator it = nu_.begin(); it != nu_.end(); ++it)
    //     std::cout << *it << std::endl;
    std::vector<int> nbx_(N + 1, 0);
    nbx_.at(0) = 12;
    //std::fill_n(nbx_.begin() + N, 1, 12);

    std::vector<int> nbu_(N + 1, 4);
    std::fill_n(nbu_.begin() + N, 1, 0);

    std::vector<int> ng_(N + 1, 4);
    std::fill_n(ng_.begin() + N, 1, 0);

    std::vector<int> nsbx_(N + 1, 0);
    std::vector<int> nsbu_(N + 1, 0);
    std::vector<int> nsg_(N + 1, 0);
    
    
    //ocp qp dim
    int dim_size = d_ocp_qp_dim_memsize(N);
    void *dim_mem = malloc(dim_size);
    if(!dim_mem) std::cout << "malloc dim memory failed" << std::endl;

    struct d_ocp_qp_dim dim;
	d_ocp_qp_dim_create(N, &dim, dim_mem);
    d_ocp_qp_dim_set_all(nx_.data(), nu_.data(), nbx_.data(), nbu_.data(), ng_.data(), nsbx_.data(), nsbu_.data(), nsg_.data(), &dim);
    
    //ocp qp part condense
    int N2 = N;// horizon length of partially condensed OCP QP
    int dim_size2 = d_ocp_qp_dim_memsize(N2);
    void *dim_mem2 = malloc(dim_size2);
    if(!dim_mem2) std::cout << "malloc dim2 memory failed" << std::endl;

	struct d_ocp_qp_dim dim2;
	d_ocp_qp_dim_create(N2, &dim2, dim_mem2);

	int *block_size = (int*)malloc((N+1)*sizeof(int));
    if(!block_size) std::cout << "malloc block memory failed" << std::endl;
	d_part_cond_qp_compute_block_size(N, N2, block_size);
    d_part_cond_qp_compute_dim(&dim, block_size, &dim2);

    //ocp qp
    int qp_size = d_ocp_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);
    if(!qp_mem) std::cout << "malloc qp memory failed" << std::endl;

    struct d_ocp_qp qp;
	d_ocp_qp_create(&dim, &qp, qp_mem);
    
    static double hA_[12 * 12];
    static double hB_[12 * 6];
    static double hb_[12 * 1];

    static double Q_[12 * 12];
    static double S_[6 * 12] = {0};
    static double R_[6 * 6];
    static double q_[12 * 1];
    static double r_[6 * 1] = {0};
    
    static int idxbu_[4 * 1];
    static double lbu_[4 * 1];
    static double ubu_[4 * 1];

    static int idxbx0_[12 * 1];
    static double x0_[12 * 1];

    static double C_[4 * 12] = {0};
    static double D_[4 * 6];
    static double lg_[4 * 1];
    static double ug_[4 * 1];

    Eigen::Matrix<int, 4, 1> Jbu;
    Eigen::Matrix<int, 12, 1> Jbx;
    Eigen::Matrix<double, 4, 6> Dn;
    Eigen::Matrix<double, 4, 1> Ulow, Uup;
    Eigen::Matrix<double, 12, 1> X_0;
    Eigen::Matrix<double, 4, 1> dlow, dup;
    Eigen::Matrix<double, 12, 1> qn;

    Eigen::DiagonalMatrix<double, 12> Q_temp;Q_temp.setZero();
    Eigen::DiagonalMatrix<double, 6> R_temp;R_temp.setZero();
    Eigen::Matrix<double, 12, 12> Q;
    Eigen::Matrix<double, 6, 6> R;
    Q_temp.diagonal() = stateWeight.replicate(1, 1);
    R_temp.diagonal() = controlWeight.replicate(1, 1);
    
    Q = Q_temp.diagonal().asDiagonal();
    R = R_temp.diagonal().asDiagonal();
    Jbu << 2, 3, 4, 5; 
    Jbx << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;

    Dn << mu, 0, 1, 0, 0, 0,
        -mu, 0, 1, 0, 0, 0,
        0, mu, 1, 0, 0, 0,
        0, -mu, 1, 0, 0, 0;
    Ulow << fzmin, torXmin, torYmin, torZmin;
    Uup << fzmax, torXmax, torYmax, torZmax;
    X_0 = X0;
    dlow << 0, 0, 0, 0;
    dup << MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, MPC_INFINITY;
    //std::cout << Q <<"\n\n\n" << R << std::endl;
    std::vector<Eigen::Matrix<double, 12, 1>> qn_con(N + 1);
    std::vector<double*> hqn_con;
    for(int i=0;i< N +1;++i)
        hqn_con.push_back(new double[12]);
    
    //std::cout << Adt << std::endl;
    
    matrix_to_doubleArray(Q_, Q, 12, 12);
    matrix_to_doubleArray(R_, R, 6, 6);
    //matrix_to_doubleArray(q_, qn, 12, 1);
    matrix_to_intArray(idxbu_, Jbu, 4, 1);
    matrix_to_doubleArray(lbu_, Ulow, 4, 1);
    matrix_to_doubleArray(ubu_, Uup, 4, 1);
    matrix_to_doubleArray(D_, Dn, 4, 6);
    matrix_to_doubleArray(lg_, dlow, 4, 1);
    matrix_to_doubleArray(ug_, dup, 4, 1);
    matrix_to_intArray(idxbx0_, Jbx, 12, 1);
    matrix_to_doubleArray(x0_, X_0, 12, 1);
    std::vector<double*> hS(N + 1, S_);
    std::vector<double*> hr(N + 1, r_);
    std::vector<int*> hidxbx(N + 1, idxbx0_);
    std::vector<double*> hlbx(N + 1, x0_);
    std::vector<double*> hubx(N + 1, x0_);
    std::fill_n(hidxbx.begin() + 1, N, nullptr);
    std::fill_n(hlbx.begin() + 1, N, nullptr);
    std::fill_n(hubx.begin() + 1, N, nullptr);

    std::vector<double*> hQ(N + 1, Q_);
    std::vector<double*> hR(N + 1, R_);
    std::vector<int*> hidxbu(N + 1, idxbu_);
    std::vector<double*> hlbu(N + 1, lbu_);
    std::vector<double*> hubu(N + 1, ubu_);

    std::vector<double*> hC(N + 1, C_);
    std::vector<double*> hD(N + 1, D_);
    std::vector<double*> hlg(N + 1, lg_);
    std::vector<double*> hug(N + 1, ug_);

    std::vector<double*> hZl(N + 1, nullptr);
    std::vector<double*> hZu(N + 1, nullptr);
    std::vector<double*> hzl(N + 1, nullptr);
    std::vector<double*> hzu(N + 1, nullptr);
    std::vector<int*> hidxs(N + 1, nullptr);
    std::vector<double*> hlls(N + 1, nullptr);
    std::vector<double*> hlus(N + 1, nullptr);

    Timer cal_time;
    getInitState();
    getConStateMatrix();
    getDisStateMatrix();
    setRefTrajectory();
    qn_con[0] = -2 * Q * Xdqp.head(12);
    for(int i = 1;i < N + 1;++i){
        qn_con[i] = -2 * Q * Xdqp.segment(12 * (i - 1), 12);
    }
    for(int i =0;i<N + 1;++i){
        matrix_to_doubleArray(hqn_con[i], qn_con[i], 12, 1);
    }
    matrix_to_doubleArray(hA_, Adt, 12, 12);
    matrix_to_doubleArray(hB_, Bdt, 12, 6);
    matrix_to_doubleArray(hb_, Cdt, 12, 1);
    std::vector<double*> AA(N, hA_);
    //cout << AA[1][0] << endl;
    std::vector<double*> BB(N, hB_);
    std::vector<double*> bb(N, hb_);
    // std::vector<double*> hq(N + 1, q_);
    std::vector<double*> hq(hqn_con);
    

    d_ocp_qp_set_all(AA.data(), BB.data(), bb.data(),
        hQ.data(), hS.data(), hR.data(), hq.data(), hr.data(),
        hidxbx.data(), hlbx.data(), hubx.data(),
        hidxbu.data(), hlbu.data(), hubu.data(),
        hC.data(), hD.data(), hlg.data(), hug.data(),
        hZl.data(), hZu.data(), hzl.data(), hzu.data(),
        hidxs.data(), hlls.data(), hlus.data(),
        &qp
        );

    //ocp qp part condense
    int qp_size2 = d_ocp_qp_memsize(&dim2);
	void *qp_mem2 = malloc(qp_size2);

	struct d_ocp_qp qp2;
	d_ocp_qp_create(&dim2, &qp2, qp_mem2);

    //ocp sol
    int qp_sol_size = d_ocp_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    if(!qp_sol_mem) std::cout << "malloc qp solve memory failed" << std::endl;

	struct d_ocp_qp_sol qp_sol;
	d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    //ocp qp sol part cond
    int qp_sol_size2 = d_ocp_qp_sol_memsize(&dim2);
	void *qp_sol_mem2 = malloc(qp_sol_size2);

	struct d_ocp_qp_sol qp_sol2;
	d_ocp_qp_sol_create(&dim2, &qp_sol2, qp_sol_mem2);

    //part condense arg
    int part_cond_arg_size = d_part_cond_qp_arg_memsize(dim2.N);
	void *part_cond_arg_mem = malloc(part_cond_arg_size);

	struct d_part_cond_qp_arg part_cond_arg;
	d_part_cond_qp_arg_create(dim2.N, &part_cond_arg, part_cond_arg_mem);

	d_part_cond_qp_arg_set_default(&part_cond_arg);

	d_part_cond_qp_arg_set_ric_alg(0, &part_cond_arg);

    //ipm arg
    int ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);

    hpipm_mode mode = hpipm_mode::SPEED;//BALANCE;//ROBUST;//SPEED;
    int iter_max = 200;
    double alpha_min = 1e-4;
    double mu0 = 10;
    double tol_stat = 1e-5;
    double tol_eq = 1e-5;
    double tol_ineq = 1e-5;
    double tol_comp = 1e-5;
    double reg_prim = 1e-8;
    int warm_start = 0;
    int pred_corr = 1;
    int ric_alg = 0;
    int split_step = 1;
	d_ocp_qp_ipm_arg_set_default(mode, &arg);
	d_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
	d_ocp_qp_ipm_arg_set_alpha_min(&alpha_min, &arg);
	d_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_ocp_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
	d_ocp_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
	d_ocp_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
	d_ocp_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
	d_ocp_qp_ipm_arg_set_reg_prim(&reg_prim, &arg);
	d_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);
	d_ocp_qp_ipm_arg_set_pred_corr(&pred_corr, &arg);
	d_ocp_qp_ipm_arg_set_ric_alg(&ric_alg, &arg);
	d_ocp_qp_ipm_arg_set_split_step(&split_step, &arg);

    //part cond workspace
    int part_cond_size = d_part_cond_qp_ws_memsize(&dim, block_size, &dim2, &part_cond_arg);
	void *part_cond_mem = malloc(part_cond_size);

	struct d_part_cond_qp_ws part_cond_ws;
	d_part_cond_qp_ws_create(&dim, block_size, &dim2, &part_cond_arg, &part_cond_ws, part_cond_mem);

    //ipm workspace
    int ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_ws workspace;
	d_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    int nrep = 1;
    int hpipm_status;
    //part cond
    for(int rep = 0; rep < nrep; rep++)
		{
		d_part_cond_qp_cond(&qp, &qp2, &part_cond_arg, &part_cond_ws);
		}
    //ipm solver
    for(int rep = 0; rep < nrep; rep++){
		// call solver
		d_ocp_qp_ipm_solve(&qp2, &qp_sol2, &arg, &workspace);
		d_ocp_qp_ipm_get_status(&workspace, &hpipm_status);
	}
    if(hpipm_status == 0)
		{
        std::cout << "hpipm partial ocp qp solved\n";
		}
	else if(hpipm_status==1)
		{
        std::cout << "hpipm partial ocp qp failed! Maximum number of iterations reached\n";
		}
	else if(hpipm_status==2)
		{
        std::cout << "hpipm partial ocp qp failed! Minimum step lenght reached\n";
		}
	else if(hpipm_status==3)
		{
        std::cout << "hpipm partial ocp qp failed! NaN in computations\n";
		}
    //part expand
    for(int rep = 0; rep < nrep; rep++)
		d_part_cond_qp_expand_sol(&qp, &qp2, &qp_sol2, &qp_sol, &part_cond_arg, &part_cond_ws);


    std::vector<Eigen::Matrix<double, 6, 1>> u_sol(N);
    std::vector<Eigen::Matrix<double, 12, 1>> x_sol(N + 1);
    for(int i = 0; i < N; ++i){
        d_ocp_qp_sol_get_u(i, &qp_sol, u_sol[i].data());
    }
    for(int i = 0; i < N + 1; ++i){
        d_ocp_qp_sol_get_x(i, &qp_sol, x_sol[i].data());
    }
    // for(int i=0;i<N;++i)
    //     std::cout << u_sol[i][2] << std::endl;
    Matrix<double,6,1> Usol_hpipm;
    Usol_hpipm = u_sol[0];

    t_hpOcpPar = cal_time.getMs();
    //std::cout << "hpipm ocp part condense solve time:  " << cal_time.getMs() << std::endl;
    free(dim_mem);
    free(qp_mem);
	free(qp_sol_mem);
	free(ipm_arg_mem);
	free(ipm_mem);
    //cout << "\n\nAdt[15]\n" << Aqp.bottomRightCorner(12, 12) << "\n\n";
    return Usol_hpipm;

}

Vec6 MpcAction::solveMpcCasadiOcp(){
    //casadi::Opti opti = casadi::Opti("conic");
    casadi::Opti opti = casadi::Opti();

    MX theta = MX::sym("theta", 3); 
    MX p = MX::sym("p", 3);  
    MX AngMomen = MX::sym("angularmomentum", 3); 
    MX LinMomen = MX::sym("linearmomentum", 3); 
    
    MX states = MX::vertcat({theta, p, AngMomen, LinMomen});
    int n_states = states.size1();

    MX force = MX::sym("force",3); 
    MX torque = MX::sym("torque", 3); 
    MX controls = MX::vertcat({force,torque});
    int n_controls = controls.size1();

    MX U = opti.variable(n_controls, HORIZON);
    MX X = opti.variable(n_states, HORIZON+1);
    MX dState = opti.parameter(n_states * HORIZON);
    MX currState = opti.parameter(n_states);
   
    MX com_to_foot = opti.parameter(3);
    MX r1_mat = crossProductMx(com_to_foot); 
    DM I = DM({{Ixx, 0, 0},{0, Iyy, 0},{0, 0, Izz}}); 

    MX RM = EulerZYX2RotMatMx(currState(S1(0,3)));
    MX gI = MX::mtimes(MX::mtimes(RM,I),RM.T());
    MX invIw = MX::inv(gI);
    MX angVelWorld2Body = getAngVelTransMatrixMx(currState(S1(0,3))); 
    
    DM g = DM({0,0,-9.8});
    
    MX theta_next = theta + MPC_CT * MX::mtimes(angVelWorld2Body, MX::mtimes(invIw,AngMomen));
    
    MX p_next = p + MPC_CT * LinMomen/m;
    MX angmomen_next = AngMomen + MPC_CT * MX::mtimes(r1_mat,force) + MPC_CT * torque;
    MX linmomen_next = LinMomen + MPC_CT * force + MPC_CT * m * g;
   
    MX states_next = MX::vertcat({theta_next, p_next, angmomen_next, linmomen_next});
    casadi::Function function_states_next = casadi::Function("function_states_next", {states, controls}, {states_next});

    DM Q_vec = DM::zeros(n_states); 
    DM R_vec = DM::zeros(n_controls); 
    std::copy(stateWeight.data(), stateWeight.data() + stateWeight.size(), Q_vec.ptr());
    std::copy(controlWeight.data(), controlWeight.data() + controlWeight.size(), R_vec.ptr());
    DM Q_mat = DM::diag(Q_vec);
    DM R_mat = DM::diag(R_vec);


    MX delta_x0 = X(S1(), 0)-dState(S1(0, n_states));
    MX obj = MX::mtimes(MX::mtimes(delta_x0.T(),Q_mat),delta_x0);
    //MX obj = 0;

    for(int i =0;i<HORIZON;i++) 
    {
      MX delta_x = X(S1(),i+1)-dState(S1(i*n_states, (i+1)*n_states));
      MX delta_u = U(S1(),i);
      obj = obj + MX::mtimes(MX::mtimes(delta_x.T(),Q_mat),delta_x) + MX::mtimes(MX::mtimes(delta_u.T(),R_mat),delta_u);
    }
    opti.minimize(obj);

    for(int i=0;i<HORIZON;i++){                                   
      opti.subject_to(X(S1(),i+1)==function_states_next(std::vector<MX>{X(S1(),i),U(S1(),i)}).at(0));
    }
    for(int i=0;i<HORIZON;i++){
      opti.subject_to( fzmin <= U(2,i) <= fzmax );
      opti.subject_to( torXmin <= U(3,i) <= torXmax );
      opti.subject_to( torYmin <= U(4,i) <= torYmax );
      opti.subject_to( torZmin <= U(5,i) <= torZmax );

      opti.subject_to(-miu*U(2,i) <= U(0,i) <= miu*U(2,i));
      opti.subject_to(-miu*U(2,i) <= U(1,i) <= miu*U(2,i));
    }
    opti.subject_to(X(S1(),0) == currState);  
    casadi::Dict opts_setting = {{"ipopt.max_iter", 500}, {"ipopt.print_level", 0}, {"print_time", 0},{"ipopt.acceptable_tol",1e-5},{"ipopt.acceptable_obj_change_tol",1e-6}};
    opti.solver("ipopt", opts_setting);
    
    // casadi::Dict opts_setting = {{"printLevel", "none"},{"sparse",false},{"schur",false},{"max_schur", 20},{"hessian_type", "posdef"}};
    // opti.solver("qpoases", opts_setting);

    // casadi::Dict opts_setting = {{"print_iter", false},{"print_header", false},{"max_iter", 20}};
    // opti.solver("qrqp", opts_setting);

    Timer sol_time;
    getInitState();
    setRefTrajectory();

    DM dm_dState = DM::zeros(12 * HORIZON, 1);
    DM dm_com2foot = DM::zeros(3, 1);
    DM dm_currState = DM::zeros(12, 1);
    std::copy(Xdqp.data(), Xdqp.data() + Xdqp.size(), dm_dState.ptr());
    std::copy(comToFoot.data(), comToFoot.data() + comToFoot.size(), dm_com2foot.ptr());
    std::copy(X0.data(), X0.data() + X0.size(), dm_currState.ptr());

    opti.set_value(currState, dm_currState);
    opti.set_value(dState, dm_dState);
    opti.set_value(com_to_foot,dm_com2foot);
   
    std::unique_ptr<casadi::OptiSol> sol;
    sol = std::make_unique<casadi::OptiSol>(opti.solve());
    
    DM u_sol = sol->value(U);
    DM x_sol = sol->value(X);

    std::vector<double> vec_u = static_cast<std::vector<double>>(u_sol(S1(),0));
    // for(int i=0;i<HORIZON;++i)
    //  cout << vec_u[i*6 + 2] << endl;
    Matrix<double, 6, 1> sol_casadi(vec_u.data());

    t_casOcp = sol_time.getMs();
    //cout << "casadi ocp solve time: " << sol_time.getMs() << endl;

    return sol_casadi;
}
Vec6 MpcAction::solveMpcCasadiDense(){
    casadi::Opti opti = casadi::Opti("conic");

    MX U = opti.variable(6* HORIZON, 1);
    DM H_mat = DM::zeros(6* HORIZON, 6* HORIZON);
    DM g_mat = DM::zeros(6* HORIZON, 1);

    Timer start_time;
    getMpcMatrix();
    Matrix<double, 6* HORIZON, 6* HORIZON, Eigen::RowMajor> H_temp;
    Matrix<double, 6* HORIZON, 1> g_temp;
    H_temp = Hqp;
    g_temp = 2 * Gqp;
    std::copy(H_temp.data(), H_temp.data() + H_temp.size(), H_mat.ptr());
    std::copy(g_temp.data(), g_temp.data() + g_temp.size(), g_mat.ptr());

    MX obj = MX::mtimes(MX::mtimes(U.T(), H_mat), U) + MX::mtimes(U.T(), g_mat);
    opti.minimize(obj);

    DM A_mat = DM::zeros(8* HORIZON, 6* HORIZON);
    DM lbA_mat = DM::zeros(8* HORIZON, 1);
    DM ubA_mat = DM::zeros(8* HORIZON, 1);
    Matrix<double, 8* HORIZON, 6* HORIZON, Eigen::RowMajor> A_temp;

    A_temp = fmat;   
    std::copy(A_temp.data(), A_temp.data() + A_temp.size(), A_mat.ptr());
    std::copy(Lbqp.data(), Lbqp.data() + Lbqp.size(), lbA_mat.ptr());
    std::copy(Ubqp.data(), Ubqp.data() + Ubqp.size(), ubA_mat.ptr());
    //opti.subject_to(lbA_mat <= MX::mtimes(A_mat, U) <= ubA_mat);
    
    // for(int i=0;i<HORIZON;i++){ 
    //   opti.set_initial(U(S1(i*6, (i+1)*6)), DM({0, 0, 5, 0, 0, 0}) );  
    // }

    //casadi::Dict opts_setting = {{"ipopt.max_iter", 500}, {"ipopt.print_level", 0}, {"print_time", 0},{"ipopt.acceptable_tol",1e-5},{"ipopt.acceptable_obj_change_tol",1e-6}};
    //opti.solver("ipopt", opts_setting);
    
    // casadi::Dict opts_setting = {{"printLevel", "none"},{"sparse",false},{"schur",false},{"max_schur", 50},{"hessian_type", "posdef"}};
    // opti.solver("qpoases", opts_setting);

    casadi::Dict opts_setting = {{"print_iter", false},{"print_header", false},{"max_iter", 20}};
    opti.solver("qrqp", opts_setting);
    std::unique_ptr<casadi::OptiSol> sol;
    sol = std::make_unique<casadi::OptiSol>(opti.solve());
    
    DM u_sol = sol->value(U);

    std::vector<double> vec_u = static_cast<std::vector<double>>(u_sol(S1(0, 6)));
    Matrix<double, 6, 1> sol_casadi_dense(vec_u.data());
    // for(int i=0;i<HORIZON;++i)
    //     cout << u_sol(i*6 + 2) << endl;
    t_casDen = start_time.getMs();
    //cout << "casadi dense solve time:  " << start_time.getMs() << endl;
    
    return sol_casadi_dense;
}

MX MpcAction::crossProductMx(MX vec)
{
    MX mat = MX::zeros(3,3);
    mat(0,1) = -vec(2);mat(0,2) = vec(1);
    mat(1,0) = vec(2); mat(1,2) = -vec(0);
    mat(2,0) = -vec(1);mat(2,1) = vec(0);
    return mat;
}
MX MpcAction::getAngVelTransMatrixMx(MX rpy)// theta_dot = R * w;
{
    MX RM = MX::zeros(3,3);
    MX su = sin(rpy(0));MX cu = cos(rpy(0));
    MX sv = sin(rpy(1));MX cv = cos(rpy(1));
    MX sw = sin(rpy(2));MX cw = cos(rpy(2));
    RM(0,0) = cw / cv;RM(0,1) =   sw/cv;RM(0,2) = 0;
    RM(1,0) =     -sw;RM(1,1) =      cw;RM(1,2) = 0;
    RM(2,0) =cw*sv/cv;RM(2,1) =sw*sv/cv;RM(2,2) = 1;
    return RM;
}
MX MpcAction::EulerZYX2RotMatMx(MX rpy) 
{
    MX RM = MX::zeros(3,3);
    MX su = sin(rpy(0));MX cu = cos(rpy(0));
    MX sv = sin(rpy(1));MX cv = cos(rpy(1));
    MX sw = sin(rpy(2));MX cw = cos(rpy(2));
    RM(0,0)=cv * cw;RM(0,1)=-cu * sw + su * cw * sv;RM(0,2)=su * sw + cu * cw * sv;
    RM(1,0)=cv * sw;RM(1,1)= cu * cw + su * sv * sw;RM(1,2)=-su * cw + cu * sv * sw;
    RM(2,0)=-sv;RM(2,1)=su * cv;RM(2,2)= cu * cv;
    return RM;
}

Vec6 MpcAction::solveMpcOsqpDense(){
    OsqpEigen::Solver solver;

    solver.settings()->setWarmStart(false);
    solver.settings()->setMaxIteration(50);
    solver.settings()->setAlpha(1.5);

    solver.data()->setNumberOfVariables(6 * HORIZON);
    solver.data()->setNumberOfConstraints(8 * HORIZON);

    Eigen::SparseMatrix<double> hessian(6*HORIZON, 6*HORIZON);
    Eigen::SparseMatrix<double> linearMatrix(8*HORIZON, 6*HORIZON);

    Timer osqp_time;
    getMpcMatrix();
    for(int row = 0;row<6 * HORIZON;++row){
        for(int col=0;col < 6*HORIZON;++col){
            hessian.insert(row, col) = Hqp(row, col);
        }
    }
    //cout << hessian <<endl;
    //cout << Hqp << endl;
    for(int row = 0;row<8 * HORIZON;++row){
        for(int col=0;col < 6*HORIZON;++col){
            linearMatrix.insert(row, col) = fmat(row, col);
        }
    }

    Eigen::VectorXd gradient;
    gradient = Gqp;

    
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(Lbqp);
    solver.data()->setUpperBound(Ubqp);
    // instantiate the solver
    solver.initSolver();

    if(solver.solve())
     cout << "osqp solved" << endl;
    Eigen::VectorXd QPSolution;
    QPSolution = solver.getSolution();
    // for(int i=0;i<HORIZON;++i)
    //  cout << QPSolution[i*6 + 2] << endl;
    Matrix<double, 6, 1> sol_osqp;
    sol_osqp = QPSolution.head(6);
    t_osqpDen = osqp_time.getMs();
    //cout << "osqp dense qp solve time  " << osqp_time.getMs() << endl;

    return sol_osqp;
}

Vec6 MpcAction::solveMpcAcadosOcp(){
    Timer cal_time;

    getInitState();
    getConStateMatrix();
    getDisStateMatrix();
    setRefTrajectory();

    static const int N = HORIZON;
    
    static double hA_[12 * 12];
    static double hB_[12 * 6];
    static double hb_[12 * 1];

    static double Q_[12 * 12];
    static double S_[6 * 12] = {0};
    static double R_[6 * 6];
    static double q_[12 * 1];
    static double r_[6 * 1] = {0};
    
    static int idxbu_[4 * 1];
    static double lbu_[4 * 1];
    static double ubu_[4 * 1];

    static int idxbx_[12 * 1];
    static double x0[12 * 1];
    static double xn[12 * 1];

    static double C_[4 * 12] = {0};
    static double D_[4 * 6];
    static double lg_[4 * 1];
    static double ug_[4 * 1];

    Eigen::Matrix<int, 4, 1> Jbu;
    Eigen::Matrix<int, 12, 1> Jbx;
    Eigen::Matrix<double, 12, 1> X_0, X_N;
    Eigen::Matrix<double, 4, 6> Dn;
    Eigen::Matrix<double, 4, 1> Ulow, Uup;
    Eigen::Matrix<double, 4, 1> dlow, dup;
    Eigen::Matrix<double, 12, 1> qn;

    Eigen::DiagonalMatrix<double, 12> Q_temp;Q_temp.setZero();
    Eigen::DiagonalMatrix<double, 6> R_temp;R_temp.setZero();
    Eigen::Matrix<double, 12, 12> Q;
    Eigen::Matrix<double, 6, 6> R;
    Q_temp.diagonal() = stateWeight.replicate(1, 1);
    R_temp.diagonal() = controlWeight.replicate(1, 1);
    
    Q = Q_temp.diagonal().asDiagonal();
    R = R_temp.diagonal().asDiagonal();

    //std::cout << Q <<"\n\n\n" << R << std::endl;
    qn = -2 * Q * Xdqp.head(12);
    Jbu << 2, 3, 4, 5; 
    Jbx << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;
    Dn << mu, 0, 1, 0, 0, 0,
        -mu, 0, 1, 0, 0, 0,
        0, mu, 1, 0, 0, 0,
        0, -mu, 1, 0, 0, 0;
    Ulow << fzmin, torXmin, torYmin, torZmin;
    Uup << fzmax, torXmax, torYmax, torZmax;
    dlow << 0, 0, 0, 0;
    dup << MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, MPC_INFINITY;
    X_0 = X0;

    //std::cout << Adt << std::endl;
    matrix_to_doubleArray(hA_, Adt, 12, 12);
    matrix_to_doubleArray(hB_, Bdt, 12, 6);
    matrix_to_doubleArray(hb_, Cdt, 12, 1);
    matrix_to_doubleArray(Q_, Q, 12, 12);
    matrix_to_doubleArray(R_, R, 6, 6);
    matrix_to_doubleArray(q_, qn, 12, 1);
    matrix_to_intArray(idxbu_, Jbu, 4, 1);
    matrix_to_doubleArray(lbu_, Ulow, 4, 1);
    matrix_to_doubleArray(ubu_, Uup, 4, 1);
    matrix_to_intArray(idxbx_, Jbx, 12, 1);
    matrix_to_doubleArray(x0, X_0, 12, 1);
    matrix_to_doubleArray(D_, Dn, 4, 6);
    matrix_to_doubleArray(lg_, dlow, 4, 1);
    matrix_to_doubleArray(ug_, dup, 4, 1);

    ocp_qp_solver_plan plan;
    plan.qp_solver = FULL_CONDENSING_HPIPM; //FULL_CONDENSING_QPOASES; // FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM
    // PARTIAL_CONDENSING_HPIPM, PARTIAL_CONDENSING_HPMPC, PARTIAL_CONDENSING_OOQP,
    //  PARTIAL_CONDENSING_OSQP, PARTIAL_CONDENSING_QPDUNES, FULL_CONDENSING_QORE, FULL_CONDENSING_OOQP,

    ocp_qp_xcond_solver_config *config = ocp_qp_xcond_solver_config_create(plan);
    ocp_qp_dims *dims = ocp_qp_dims_create(N);
    int nx = 12;
    int nu = 6;
    int nbu = 4;
    int ng = 4;
    int nbx = 12;
    int nzero = 0;
    for (int i = 0; i < N+1; i++){
        ocp_qp_dims_set(config, dims, i, "nx", &nx);
    }
    for (int i = 0; i < N; i++){
        ocp_qp_dims_set(config, dims, i, "nu", &nu);
        ocp_qp_dims_set(config, dims, i, "nbu", &nbu);
        ocp_qp_dims_set(config, dims, i, "ng", &ng);
    }
    ocp_qp_dims_set(config, dims, 0, "nbx", &nbx);
    ocp_qp_dims_set(config, dims, N, "nu", &nzero);

    ocp_qp_in *qp_in = ocp_qp_in_create(dims);

    for (int i = 0; i < N; i++)
    {
        ocp_qp_in_set(config, qp_in, i, "A", hA_);
        ocp_qp_in_set(config, qp_in, i, "B", hB_);
        ocp_qp_in_set(config, qp_in, i, "b", hb_);
        ocp_qp_in_set(config, qp_in, i, "Q", Q_);
        ocp_qp_in_set(config, qp_in, i, "S", S_);
        ocp_qp_in_set(config, qp_in, i, "R", R_);
        ocp_qp_in_set(config, qp_in, i, "q", q_);
        ocp_qp_in_set(config, qp_in, i, "r", r_);

        ocp_qp_in_set(config, qp_in, i, "idxbu", idxbu_);
        ocp_qp_in_set(config, qp_in, i, "lbu", lbu_);
        ocp_qp_in_set(config, qp_in, i, "ubu", ubu_);

        ocp_qp_in_set(config, qp_in, i, "C", C_);
        ocp_qp_in_set(config, qp_in, i, "D", D_);
        ocp_qp_in_set(config, qp_in, i, "lg", lg_);
        ocp_qp_in_set(config, qp_in, i, "ug", ug_);
    }
    ocp_qp_in_set(config, qp_in, 0, "idxbx", idxbx_);
    ocp_qp_in_set(config, qp_in, 0, "lbx", x0);
    ocp_qp_in_set(config, qp_in, 0, "ubx", x0);

    //std::cout << X_N.transpose() << std::endl;

    ocp_qp_xcond_solver_dims *solver_dims =
                ocp_qp_xcond_solver_dims_create_from_ocp_qp_dims(config, dims);


    void *opts = ocp_qp_xcond_solver_opts_create(config, solver_dims);


    // set partial condensing option
    if (plan.qp_solver == PARTIAL_CONDENSING_HPIPM ||
        // plan.qp_solver == PARTIAL_CONDENSING_HPMPC ||
        // plan.qp_solver == PARTIAL_CONDENSING_OOQP ||
        plan.qp_solver == PARTIAL_CONDENSING_OSQP ||
        plan.qp_solver ==  PARTIAL_CONDENSING_QPDUNES)
    {
        int N2 = N;
        ocp_qp_xcond_solver_opts_set(config, (ocp_qp_xcond_solver_opts*)opts, "cond_N", &N2);
    }

    ocp_qp_out *qp_out = ocp_qp_out_create(dims);

    ocp_qp_solver *qp_solver = ocp_qp_create(config, solver_dims, opts);

    int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

    if (acados_return == ACADOS_SUCCESS)
        std::cout << "acados ocp qp success" << std::endl;
    else if(acados_return == ACADOS_FAILURE)
        std::cout << "acados ocp qp failure" << std::endl;
    else if(acados_return == ACADOS_MAXITER)
        std::cout << "acados ocp qp maximum number of iterations reached" << std::endl;
    else if(acados_return == ACADOS_MINSTEP)
        std::cout << "acados ocp qp  minimum step size in QP solver reached" << std::endl;
    else if(acados_return == ACADOS_QP_FAILURE)
        std::cout << "acados ocp qp solver failed" << std::endl;

    // printf("\nqp output:\n");
    //print_ocp_qp_out(qp_out);
    Matrix<double, 6*N, 1> u;
    for(int r=0;r<N;++r){
        for(int c=0;c<6;++c){
            /*
            ocp_qp_out = d_ocp_qp_sol https://github1s.com/acados/acados/blob/HEAD/acados/ocp_qp/ocp_qp_common.h#L54
            d_ocp_qp_sol https://github1s.com/giaf/hpipm/blob/7f8e8de828be86b6ae3d373e8e56f9f3b1fd45aa/include/hpipm_d_ocp_qp_sol.h
            blasfeo_dvec https://usermanual.wiki/Pdf/guide.1256587877/view
            */
            u[6*r + c] = qp_out->ux->pa[(12+6)*r + 12 + c];
        }
    }
    t_acados0cp = cal_time.getMs();
    //print_ocp_qp_out(qp_out);
    return u.head(6);  
}

Vec6 MpcAction::solveMpcAacadosDense(){
    Timer cal_time;
    getMpcMatrix();

    const int N = HORIZON;
    double hH[6 * N * 6 * N];
    double hg[6 * N];

    // double hA[6 * N * 6 * N] = {0};
    // double hb[6 * N] = {0};

    int hidv[4 * N] = {0};
    double hlv[4 * N] = {0};
    double huv[4 * N] = {0};
    
    double hC[4 * N * 6 * N] = {0};
    double hlg[4 * N] = {0};
    double hug[4 * N] = {0};
    
    Eigen::Matrix<int, 4, 1> idv;
    Eigen::Matrix<int, 4 * N, 1> idv_temp;
    idv << 2, 3, 4, 5;
    idv_temp = idv.replicate(N, 1);
    //std::cout << idv_temp << std::endl;
    Eigen::Matrix<double, 4, 1> lv, uv, lg, ug;
    Eigen::Matrix<double, 4, 6> hc_;
    Eigen::Matrix<double, 4 * N, 1> lv_temp, uv_temp, lg_temp, ug_temp;
    Eigen::Matrix<double, 4 * N, 6 * N> hc_temp;
    hc_temp.setZero();

    lv << fzmin, torXmin, torYmin, torZmin;
    uv << fzmax, torXmax, torYmax, torZmax;
    lg << 0, 0, 0, 0;
    ug << MPC_INFINITY, MPC_INFINITY, MPC_INFINITY, MPC_INFINITY;
    hc_ << mu, 0, 1, 0, 0, 0,
        -mu, 0, 1, 0, 0, 0,
        0, mu, 1, 0, 0, 0,
        0, -mu, 1, 0, 0, 0;

    lv_temp = lv.replicate(N, 1); 
    uv_temp = uv.replicate(N, 1);
    lg_temp = lg.replicate(N, 1);
    ug_temp = ug.replicate(N, 1);

    matrix_to_doubleArray(hH, Hqp, 6 * N, 6 * N);
    matrix_to_doubleArray(hg, Gqp, 6 * N, 1);
    for(int i=0;i<N;++i)
        hc_temp.block(4 * i, 6 * i, 4, 6) = hc_;
    //std::cout << hc_temp << std::endl;
    matrix_to_intArray(hidv, idv_temp, 4 * N, 1);
    matrix_to_doubleArray(hlv, lv_temp, 4 * N, 1);
    matrix_to_doubleArray(huv, uv_temp, 4 * N, 1);
    matrix_to_doubleArray(hlg, lg_temp, 4 * N, 1);
    matrix_to_doubleArray(hug, ug_temp, 4 * N, 1);
    matrix_to_doubleArray(hC, hc_temp, 4 * N, 6 * N);

    dense_qp_solver_plan plan;
    plan.qp_solver = DENSE_QP_HPIPM;//DENSE_QP_HPIPM, DENSE_QP_QORE, DENSE_QP_QPOASES, DENSE_QP_OOQP
    qp_solver_config *config = dense_qp_config_create(&plan);
    dense_qp_dims dims;
    dims.nv = 6 * N;
    dims.ne = 0; 
    dims.nb = 4 * N;
    dims.ng = 4 * N;
    // dims.nb = 0;
    // dims.ng = 0;
    dims.ns = 0;
    dims.nsb = 0;
    dims.nsg = 0;
    dense_qp_in *qp_in = dense_qp_in_create(config, &dims);

    d_dense_qp_set_all(hH, hg, 
        nullptr, nullptr, 
        hidv, hlv, huv, 
        hC, hlg, hug, 
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 
        qp_in
    );

    // for(int i=0;i<4*N;++i)
    //   std::cout << hlv[i] << "  ";
    void *opts = dense_qp_opts_create(config, &dims);
    dense_qp_out *qp_out = dense_qp_out_create(config, &dims);
    dense_qp_solver *qp_solver = dense_qp_create(config, &dims, opts);
    int acados_return = dense_qp_solve(qp_solver, qp_in, qp_out);
    if (acados_return == ACADOS_SUCCESS)
        std::cout << "acados dense qp success" << std::endl;
    else if(acados_return == ACADOS_FAILURE)
        std::cout << "acados dense qp failure" << std::endl;
    else if(acados_return == ACADOS_MAXITER)
        std::cout << "acados dense qp maximum number of iterations reached" << std::endl;
    else if(acados_return == ACADOS_MINSTEP)
        std::cout << "acados dense qp  minimum step size in QP solver reached" << std::endl;
    else if(acados_return == ACADOS_QP_FAILURE)
        std::cout << "acados dense qp solver failed" << std::endl;

    Matrix<double, 6*N, 1> u;
    for(int r=0;r<N;++r){
        for(int c=0;c<6;++c){
            /*
            dense_qp_out = d_dense_qp_sol https://github1s.com/acados/acados/blob/HEAD/acados/dense_qp/dense_qp_common.h#L51
            d_dense_qp_sol https://github1s.com/giaf/hpipm/blob/7f8e8de828be86b6ae3d373e8e56f9f3b1fd45aa/include/hpipm_d_dense_qp_sol.h#L56
            blasfeo_dvec https://usermanual.wiki/Pdf/guide.1256587877/view
            */
            u[6*r + c] = qp_out->v->pa[6*r + c];
        }
    }
    t_acadosDense = cal_time.getMs();
    //print_ocp_qp_out(qp_out);
    return u.head(6);  
}
