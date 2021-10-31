#include <stdio.h>
#include <stdlib.h>

#include "acados/utils/print.h"
#include "acados_c/ocp_qp_interface.h"

#define QP_HORIZON 5

int main() {

    double A[] = {1, 0, 1, 1};
    double B[] = {0, 1};
    double b[] = {0, 0};

    double Q[] = {1, 0, 0, 1};
    double S[] = {0, 0};
    double R[] = {1};
    double q[] = {1, 1};
    double r[] = {0};

    double x0[] = {1, 1};
    int idxb0[] = {0, 1};

    double xn[] = {-0.8, 0.8};
    int idxbn[] = {0, 1};

    ocp_qp_solver_plan plan;
    plan.qp_solver = FULL_CONDENSING_HPIPM; //FULL_CONDENSING_QPOASES; // FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM
    // PARTIAL_CONDENSING_HPIPM, PARTIAL_CONDENSING_HPMPC, PARTIAL_CONDENSING_OOQP,
    //  PARTIAL_CONDENSING_OSQP, PARTIAL_CONDENSING_QPDUNES, FULL_CONDENSING_QORE, FULL_CONDENSING_OOQP,

    ocp_qp_xcond_solver_config *config = ocp_qp_xcond_solver_config_create(plan);

    int N = QP_HORIZON;
    ocp_qp_dims *dims = ocp_qp_dims_create(N);

    int nx = 2;
    int nu = 1;
    // here: no general linear constraints (ng), soft constraints (ns, nsbx, nsbu, nsg)

    for (int i = 0; i < N+1; i++)
    {
        ocp_qp_dims_set(config, dims, i, "nx", &nx);
    }
    for(int i = 0; i < N; ++i){
         ocp_qp_dims_set(config, dims, i, "nu", &nu);
    }
    // initial value for x
    ocp_qp_dims_set(config, dims, 0, "nbx", &nx);
    // final value for x
    ocp_qp_dims_set(config, dims, N, "nbx", &nx);

    // printf("\nqp dimensions:\n");
    // print_ocp_qp_dims(dims);

    ocp_qp_in *qp_in = ocp_qp_in_create(dims);

    for (int i = 0; i < N; i++)
    {
        ocp_qp_in_set(config, qp_in, i, "A", A);
        ocp_qp_in_set(config, qp_in, i, "B", B);
        ocp_qp_in_set(config, qp_in, i, "b", b);
        ocp_qp_in_set(config, qp_in, i, "Q", Q);
        ocp_qp_in_set(config, qp_in, i, "S", S);
        ocp_qp_in_set(config, qp_in, i, "R", R);
        ocp_qp_in_set(config, qp_in, i, "q", q);
        ocp_qp_in_set(config, qp_in, i, "r", r);
    }
    ocp_qp_in_set(config, qp_in, 0, "idxbx", idxb0);
    ocp_qp_in_set(config, qp_in, 0, "lbx", x0);
    ocp_qp_in_set(config, qp_in, 0, "ubx", x0);

    ocp_qp_in_set(config, qp_in, N, "idxbx", idxbn);
    ocp_qp_in_set(config, qp_in, N, "lbx", xn);
    ocp_qp_in_set(config, qp_in, N, "ubx", xn);

    // printf("\nqp input:\n");
    // print_ocp_qp_in(qp_in);

    ocp_qp_xcond_solver_dims *solver_dims =
                ocp_qp_xcond_solver_dims_create_from_ocp_qp_dims(config, dims);


    void *opts = ocp_qp_xcond_solver_opts_create(config, solver_dims);


    // set partial condensing option
    if (plan.qp_solver == PARTIAL_CONDENSING_HPIPM)
    // are not defined by default
        // plan.qp_solver == PARTIAL_CONDENSING_HPMPC ||
        // plan.qp_solver == PARTIAL_CONDENSING_OOQP ||
        // plan.qp_solver == PARTIAL_CONDENSING_OSQP ||
        // plan.qp_solver ==  PARTIAL_CONDENSING_QPDUNES)
    {
        int N2 = 2;
        ocp_qp_xcond_solver_opts_set(config, (ocp_qp_xcond_solver_opts*)opts, "cond_N", &N2);
    }

    ocp_qp_out *qp_out = ocp_qp_out_create(dims);

    ocp_qp_solver *qp_solver = ocp_qp_create(config, solver_dims, opts);

    int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);



    // printf("\nqp output:\n");
    print_ocp_qp_out(qp_out);

    return 0;
}