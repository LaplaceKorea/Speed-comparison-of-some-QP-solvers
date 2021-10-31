#include <Eigen/Dense>
#include <iostream>
#include <vector>
extern "C" {
#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_timing.h>
}

void matrix_to_doubleArray(double dst[], Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> src, int rows, int cols);
void matrix_to_intArray(int* dst, Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> src, int rows, int cols);

void matrix_to_intArray(int* dst, Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> src, int rows, int cols)
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
void matrix_to_doubleArray(double dst[], Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> src, int rows, int cols)
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
int main(void){
    // QP size

// horizon lenght
int N = 5;
// number of input
static int nnu[5] = {1, 1, 1, 1, 1};
// number of states
static int nnx[6] = {2, 2, 2, 2, 2, 2};
// number of input box constraints
static int nnbu[5] = {0, 0, 0, 0, 0};
// number of states box constraints
static int nnbx[6] = {1, 0, 0, 0, 0, 0};
// number of general constraints
static int nng[6] = {2, 2, 2, 2, 2, 2};
// number of softed constraints on state box constraints
static int nnsbx[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on input box constraints
static int nnsbu[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on general constraints
static int nnsg[6] = {0, 0, 0, 0, 0, 0};
// number of input box constraints considered as equality
static int nnbue[6] = {0, 0, 0, 0, 0, 0};
// number of states box constraints considered as equality
static int nnbxe[6] = {0, 0, 0, 0, 0, 0};
// number of general constraints considered as equality
static int nnge[6] = {0, 0, 0, 0, 0, 0};


// QP data

//
static double A[] = {1, 0, 1, 1};
//
static double B[] = {0, 1};
//
static double b[] = {0, 0};

//
static double Q[] = {1, 0, 0, 1};
//
static double R[] = {1};
//
static double S[] = {0, 0};
//
static double q[] = {1, 1};
//
static double r[] = {0};

//
static double lbx0[] = {7, 4};
//
static double ubx0[] = {7, 4};
//
static int idxbx0[] = {0, 1};

//
static double u_guess[] = {0};
//
static double x_guess[] = {0, 0};
//
static double sl_guess[] = {};
//
static double su_guess[] = {};

// array of pointers

//
static double *AA[5] = {A, A, A, A, A};
//
static double *BB[5] = {B, B, B, B, B};
//
static double *bb[5] = {b, b, b, b, b};
//
static double *QQ[6] = {Q, Q, Q, Q, Q, Q};
//
static double *RR[6] = {R, R, R, R, R, R};
//
static double *SS[6] = {S, S, S, S, S, S};
//
static double *qq[6] = {q, q, q, q, q, q};
//
static double *rr[6] = {r, r, r, r, r, r};
//
static int *iidxbx[6] = {idxbx0, NULL, NULL, NULL, NULL, NULL};
//
static double *llbx[6] = {lbx0, NULL, NULL, NULL, NULL, NULL};
//
static double *uubx[6] = {ubx0, NULL, NULL, NULL, NULL, NULL};
//
static int *iidxbu[6] = {};
//
static double *llbu[6] = {};
//
static double *uubu[6] = {};
//
static double *CC[6] = {};
//
static double *DD[6] = {};
//
static double *llg[6] = {};
//
static double *uug[6] = {};
//
static double *ZZl[6] = {};
//
static double *ZZu[6] = {};
//
static double *zzl[6] = {};
//
static double *zzu[6] = {};
//
static int *iidxs[6] = {};
//
static double *llls[6] = {};
//
static double *llus[6] = {};
//
static double *iidxe[6] = {};

//
static double *uu_guess[6] = {u_guess, u_guess, u_guess, u_guess, u_guess, u_guess};
//
static double *xx_guess[6] = {x_guess, x_guess, x_guess, x_guess, x_guess, x_guess};
//
static double *ssl_guess[6] = {sl_guess, sl_guess, sl_guess, sl_guess, sl_guess, sl_guess};
//
static double *ssu_guess[6] = {su_guess, su_guess, su_guess, su_guess, su_guess, su_guess};



// export as global data

int *nu = nnu;
int *nx = nnx;
int *nbu = nnbu;
int *nbx = nnbx;
int *ng = nng;
int *nsbx = nnsbx;
int *nsbu = nnsbu;
int *nsg = nnsg;
int *nbue = nnbue;
int *nbxe = nnbxe;
int *nge = nnge;

double **hA = AA;
double **hB = BB;
double **hb = bb;
double **hQ = QQ;
double **hR = RR;
double **hS = SS;
double **hq = qq;
double **hr = rr;
int **hidxbx = iidxbx;
double **hlbx = llbx;
double **hubx = uubx;
int **hidxbu = iidxbu;
double **hlbu = llbu;
double **hubu = uubu;
double **hC = CC;
double **hD = DD;
double **hlg = llg;
double **hug = uug;
double **hZl = ZZl;
double **hZu = ZZu;
double **hzl = zzl;
double **hzu = zzu;
int **hidxs = iidxs;
double **hlls = llls;
double **hlus = llus;
double **hidxe = iidxe;

double **hu_guess = uu_guess;
double **hx_guess = xx_guess;
double **hsl_guess = ssl_guess;
double **hsu_guess = ssu_guess;

// arg
hpipm_mode mode = hpipm_mode::SPEED;
int iter_max = 30;
double alpha_min = 1e-8;
double mu0 = 1e4;
double tol_stat = 1e-4;
double tol_eq = 1e-5;
double tol_ineq = 1e-5;
double tol_comp = 1e-5;
double reg_prim = 1e-12;
int warm_start = 0;
int pred_corr = 1;
int ric_alg = 0;
int split_step = 1;

int ii, jj;

	int hpipm_status;

	int rep, nrep=10;

	hpipm_timer timer;

/************************************************
* ocp qp dim
************************************************/

	hpipm_size_t dim_size = d_ocp_qp_dim_memsize(N);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_ocp_qp_dim_create(N, &dim, dim_mem);

	d_ocp_qp_dim_set_all(nx, nu, nbx, nbu, ng, nsbx, nsbu, nsg, &dim);

//	d_ocp_qp_dim_codegen("examples/c/data/test_d_ocp_data.c", "w", &dim);

/************************************************
* ocp qp
************************************************/

	hpipm_size_t qp_size = d_ocp_qp_memsize(&dim);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_ocp_qp_create(&dim, &qp, qp_mem);

    Eigen::Matrix<double, 2, 2> A_temp, C_temp;
    Eigen::Matrix<double, 2, 1> B_temp, b_temp, lg_temp, ug_temp;
     Eigen::Matrix<double, 1, 1> ubx_temp, lbx_temp;
    Eigen::Matrix<int, 1, 1> idxbx_temp;
    A_temp << 1, 0, 1, 1;
    B_temp << 0, 1;
    b_temp << 0, 0;
    idxbx_temp << 0;
    ubx_temp << 7;
    lbx_temp << 7;

    C_temp << 1, 3, -1, 3;
    lg_temp << 1, 0;
    ug_temp << 7, 9;
    static double hA_temp[2 * 2];
    // for(int i = 0; i< 4;++i)
    //     std::cout << hA_temp[i] << std::endl;
    static double hB_temp[2 * 1];
    static double hb_temp[2 * 1];
    static int hidxbx_temp[1 * 1];
    static double hlbx_temp[1 * 1];
    static double hubx_temp[1 * 1];
    static double hC_temp[2 * 2];
    static double hD_temp[2 * 1] = {0};
    static double hlg_temp[2 * 1];
    static double hug_temp[2 * 1];
    matrix_to_doubleArray(hA_temp, A_temp, 2, 2);
    matrix_to_doubleArray(hB_temp, B_temp, 2, 1);
    matrix_to_doubleArray(hb_temp, b_temp, 2, 1);
    matrix_to_intArray(hidxbx_temp, idxbx_temp, 1, 1);
    matrix_to_doubleArray(hlbx_temp, lbx_temp, 1, 1);
    matrix_to_doubleArray(hubx_temp, ubx_temp, 1, 1);
    matrix_to_doubleArray(hC_temp, C_temp, 2, 2);
    matrix_to_doubleArray(hlg_temp, lg_temp, 2, 1);
    matrix_to_doubleArray(hug_temp, ug_temp, 2, 1);
    std::vector<double*> hA_(N, hA_temp);
    std::vector<double*> hB_(N, hB_temp);
    std::vector<double*> hb_(N, hb_temp);
     std::vector<int*>hidxbx_(N + 1, hidxbx_temp);
    std::vector<double*>hlbx_(N + 1, hlbx_temp);
    std::vector<double*>hubx_(N + 1, hubx_temp);
    std::vector<int*>hidxbu_(N + 1, nullptr);
    std::vector<double*>hlbu_(N + 1, nullptr);
    std::vector<double*>hubu_(N + 1, nullptr);
    std::vector<double*>hC_(N + 1, hC_temp);
    std::vector<double*>hD_(N + 1, hD_temp);
    std::vector<double*>hlg_(N + 1, hlg_temp);
    std::vector<double*>hug_(N + 1, hug_temp);


    std::vector<double*> hZl_(N + 1, nullptr);
    std::vector<double*> hZu_(N + 1, nullptr);
    std::vector<double*> hzl_(N + 1, nullptr);
    std::vector<double*> hzu_(N + 1, nullptr);
    std::vector<int*> hidxs_(N + 1, nullptr);
    std::vector<double*> hlls_(N + 1, nullptr);
    std::vector<double*> hlus_(N + 1, nullptr);
    //std::vector<double*> hr_(N + 1, nullptr);
	 d_ocp_qp_set_all(hA_.data(), hB_.data(), hb_.data(), 
     hQ, hS, hR, hq, hr, 
     hidxbx_.data(), hlbx_.data(), hubx_.data(), 
     hidxbu_.data(), hlbu_.data(), hubu_.data(), 
     hC_.data(), hD_.data(), hlg_.data(), hug_.data(), 
     hZl_.data(), hZu_.data(), hzl_.data(), hzu_.data(),
     hidxs_.data(), hlls_.data(), hlus_.data(), 
     &qp
     );
//d_ocp_qp_set_all(hA, hB, hb, hQ, hS, hR, hq, hr, hidxbx, hlbx, hubx, hidxbu, hlbu, hubu, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hlls, hlus, &qp);
//	d_ocp_qp_codegen("examples/c/data/test_d_ocp_data.c", "a", &dim, &qp);

/************************************************
* ocp qp sol
************************************************/

	hpipm_size_t qp_sol_size = d_ocp_qp_sol_memsize(&dim);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

/************************************************
* ipm arg
************************************************/

	hpipm_size_t ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);

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

//	d_ocp_qp_ipm_arg_codegen("examples/c/data/test_d_ocp_data.c", "a", &dim, &arg);

/************************************************
* ipm workspace
************************************************/

	hpipm_size_t ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_ws workspace;
	d_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

/************************************************
* ipm solver
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		// call solver
		d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
		d_ocp_qp_ipm_get_status(&workspace, &hpipm_status);
		}

	double time_ipm = hpipm_toc(&timer) / nrep;

/************************************************
* print solution info
************************************************/

    printf("\nHPIPM returned with flag %i.\n", hpipm_status);
    if(hpipm_status == 0)
		{
        printf("\n -> QP solved!\n");
		}
	else if(hpipm_status==1)
		{
        printf("\n -> Solver failed! Maximum number of iterations reached\n");
		}
	else if(hpipm_status==2)
		{
        printf("\n -> Solver failed! Minimum step lenght reached\n");
		}
	else if(hpipm_status==3)
		{
        printf("\n -> Solver failed! NaN in computations\n");
		}
	else
		{
        printf("\n -> Solver failed! Unknown return flag\n");
		}
    printf("\nAverage solution time over %i runs: %e [s]\n", nrep, time_ipm);
	printf("\n\n");

/************************************************
* extract and print solution
************************************************/
    std::vector<Eigen::Matrix<double, 1, 1>> u_sol(N);
    std::vector<Eigen::Matrix<double, 2, 1>> x_sol(N+1);
    for(int i = 0; i < N; ++i){
        d_ocp_qp_sol_get_u(i, &qp_sol, u_sol[i].data());
    }
    for(int i = 0; i < N + 1; ++i){
        d_ocp_qp_sol_get_x(i, &qp_sol, x_sol[i].data());
    }

    std::cout << "first step x:  " << x_sol[0][0] << "  "<<  x_sol[0][1]<< std::endl;

	return 0;


}