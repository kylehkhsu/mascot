/*
 * RungeKutta4.hh
 *
 *  created on: 22.04.2015
 *      author: rungger
 */

#ifndef RUNGEKUTTA4_HH_
#define RUNGEKUTTA4_HH_

/* class: RungeKutta4 
 * Fixed step size ode solver implementing a Runge Kutta scheme of order 4 */
class OdeSolver {
private:
  /* dimension */
  const double dim_;
  /* number of intermediate steps */
  const int nint_;
  /* intermidiate step size = tau_/nint_ */
  const double h_;

public:
  const double tau_;
  /* function: OdeSolver
   *
   * Input:
   * dim - state space dimension 
   * nint - number of intermediate steps
   * tau - sampling time
   */
  OdeSolver(const int dim, const int nint, const double tau) : dim_(dim), nint_(nint), h_(tau/nint), tau_(tau) { }

  /* function: ()
   *
   * Input:
   * rhs - rhs of ode  dx/dt = rhs(x,u)
   * x - current state x
   * u - current input u
   */
  template<class RHS, class X, class U>
  inline void operator()(RHS rhs, X &x, U &u) {
		X k[4];
		X tmp;

		for(int t=0; t<nint_; t++) {
			rhs(k[0],x,u);
			for(int i=0;i<dim_;i++)
				tmp[i]=x[i]+h_/2*k[0][i];

			rhs(k[1],tmp, u);
			for(int i=0;i<dim_;i++)
				tmp[i]=x[i]+h_/2*k[1][i];

			rhs(k[2],tmp, u);
			for(int i=0;i<dim_;i++)
				tmp[i]=x[i]+h_*k[2][i];

			rhs(k[3],tmp, u);
			for(int i=0; i<dim_; i++)
				x[i] = x[i] + (h_/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
		}
	}
    
    /* function: ()
     *
     * Input:
     * rhs - rhs of ode  dx/dt = rhs(x,u)
     * x - current state x
     * u - current input u
     * obstacles - set of hyper-rectangular obstacles as vector. Each element of the vector correspond to one obstacle, whose elements are arranged as: {-lb_x1, ub_x1, -lb_x2, ub_x2, ...} where x1, x2, ... are the state variables, and lb, ub represent the lower and upper bound respectively
     * unsafeAt - point of collision if collision is detected
     */
    template<class RHS, class X, class U, class T>
    inline void operator()(RHS rhs, X &x, U &u, const std::vector<std::vector<T>> obstacles, std::vector<T>& unsafeAt) {
        X k[4];
        X tmp;
        /* arrange the obstacles in suitable form */
        std::vector<std::vector<double>> lb,ub;
        for (int j=0; j<obstacles.size(); j++) {
            std::vector<double> l,u;
            for (size_t k=0; k<dim_; k++) {
                l.push_back(-obstacles[j][2*k]);
                u.push_back(obstacles[j][2*k+1]);
            }
            lb.push_back(l);
            ub.push_back(u);
        }
        
        bool col_flag = false;
        for(int t=0; t<nint_; t++) {
            rhs(k[0],x,u);
            for(int i=0;i<dim_;i++)
                tmp[i]=x[i]+h_/2*k[0][i];
            
            rhs(k[1],tmp, u);
            for(int i=0;i<dim_;i++)
                tmp[i]=x[i]+h_/2*k[1][i];
            
            rhs(k[2],tmp, u);
            for(int i=0;i<dim_;i++)
                tmp[i]=x[i]+h_*k[2][i];
            
            rhs(k[3],tmp, u);
            for(int i=0; i<dim_; i++)
                x[i] = x[i] + (h_/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
            
            if (!col_flag) { /* not detected a collision yet */
                /* collision check */
                unsafeAt.clear();
                /* iterate over all obstacles (boxes) */
                for (int j=0; j<obstacles.size(); j++) {
                    col_flag = true;
                    for (int k=0; k<dim_; k++) {
                        if (x[k]<lb[j][k] || x[k]>ub[j][k]) {
                            col_flag = false;
                            break;
                        }
                    }
                    if (col_flag) {
                        /* end collision check, but finish the computation up to tau */
                        for (int i=0; i<dim_; i++) {
                            unsafeAt.push_back(x[i]);
                        }
                        break;
                    }
                }
            }
        }
    }
    // kaushik
    inline double getStepSize() {
        return h_;
    }
};

#endif /* RUNGEKUTTA4_HH_ */
