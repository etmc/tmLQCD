
extern int geometric(int k,int l, double q2, double eps_sq);
extern int solve_cg(int k,int l, double q2, double eps_sq);
extern int bicg(int k,int l, double q2, double eps_sq);
extern int eva(double *lambda, int k, double q_off, double eps_sq);
extern int evamax(double *lambda, int k, double q_off, double eps_sq);
extern int evamax0(double *lambda, int k, double q_off, double eps_sq);
