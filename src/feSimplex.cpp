#include "feSimplex.h"

// #include "soplex.h"

// using namespace soplex;
using namespace std;

int n, m;
double A[MAX_M][MAX_N], b[MAX_M], c[MAX_N], v;
int N[MAX_N], B[MAX_M]; // nonbasic & basic

// pivot yth variable around xth constraint
inline void pivot(int x, int y)
{
  // printf("Pivoting variable %d around constraint %d.\n", y, x);

  // first rearrange the x-th row
  for(int j = 0; j < n; j++) {
    if(j != y) {
      A[x][j] /= -A[x][y];
    }
  }
  b[x] /= -A[x][y];
  A[x][y] = 1.0 / A[x][y];

  // now rearrange the other rows
  for(int i = 0; i < m; i++) {
    if(i != x) {
      for(int j = 0; j < n; j++) {
        if(j != y) {
          A[i][j] += A[i][y] * A[x][j];
        }
      }
      b[i] += A[i][y] * b[x];
      A[i][y] *= A[x][y];
    }
  }

  // now rearrange the objective function
  for(int j = 0; j < n; j++) {
    if(j != y) {
      c[j] += c[y] * A[x][j];
    }
  }
  v += c[y] * b[x];
  c[y] *= A[x][y];

  // finally, swap the basic & nonbasic variable
  swap(B[x], N[y]);
}

// Run a single iteration of the simplex algorithm.
// Returns: 0 if OK, 1 if STOP, -1 if UNBOUNDED
inline int iterate_simplex()
{
  // printf("--------------------\n");
  // printf("State:\n");
  // printf("Maximise: ");
  // for (int j=0;j<n;j++) printf("%lfx_%d + ", c[j], N[j]);
  // printf("%lf\n", v);
  // printf("Subject to:\n");
  // for (int i=0;i<m;i++)
  // {
  //     for (int j=0;j<n;j++) printf("%lfx_%d + ", A[i][j], N[j]);
  //     printf("%lf = x_%d\n", b[i], B[i]);
  // }

  // getchar(); // uncomment this for debugging purposes!

  int ind = -1, best_var = -1;
  for(int j = 0; j < n; j++) {
    if(c[j] > 0) {
      if(best_var == -1 || N[j] < ind) {
        ind = N[j];
        best_var = j;
      }
    }
  }
  if(ind == -1) return 1;

  double max_constr = INFINITY;
  int best_constr = -1;
  for(int i = 0; i < m; i++) {
    if(A[i][best_var] < 0) {
      double curr_constr = -b[i] / A[i][best_var];
      if(curr_constr < max_constr) {
        max_constr = curr_constr;
        best_constr = i;
      }
    }
  }
  if(isinf(max_constr))
    return -1;
  else
    pivot(best_constr, best_var);

  return 0;
}

// (Possibly) converts the LP into a slack form with a feasible basic solution.
// Returns 0 if OK, -1 if INFEASIBLE
inline int initialise_simplex()
{
  int k = -1;
  double min_b = -1;
  for(int i = 0; i < m; i++) {
    if(k == -1 || b[i] < min_b) {
      k = i;
      min_b = b[i];
    }
  }

  if(b[k] >= 0) // basic solution feasible!
  {
    for(int j = 0; j < n; j++) N[j] = j;
    for(int i = 0; i < m; i++) B[i] = n + i;
    return 0;
  }

  // generate auxiliary LP
  n++;
  for(int j = 0; j < n; j++) N[j] = j;
  for(int i = 0; i < m; i++) B[i] = n + i;

  // store the objective function
  double c_old[MAX_N];
  for(int j = 0; j < n - 1; j++) c_old[j] = c[j];
  double v_old = v;

  // aux. objective function
  c[n - 1] = -1;
  for(int j = 0; j < n - 1; j++) c[j] = 0;
  v = 0;
  // aux. coefficients
  for(int i = 0; i < m; i++) A[i][n - 1] = 1;

  // perform initial pivot
  pivot(k, n - 1);

  // now solve aux. LP
  int code;
  while(!(code = iterate_simplex()))
    ;

  // assert(code == 1); // aux. LP cannot be unbounded!!!

  if(v != 0) return -1; // infeasible!

  int z_basic = -1;
  for(int i = 0; i < m; i++) {
    if(B[i] == n - 1) {
      z_basic = i;
      break;
    }
  }

  // if x_n basic, perform one degenerate pivot to make it nonbasic
  if(z_basic != -1) pivot(z_basic, n - 1);

  int z_nonbasic = -1;
  for(int j = 0; j < n; j++) {
    if(N[j] == n - 1) {
      z_nonbasic = j;
      break;
    }
  }
  assert(z_nonbasic != -1);

  for(int i = 0; i < m; i++) {
    A[i][z_nonbasic] = A[i][n - 1];
  }
  swap(N[z_nonbasic], N[n - 1]);

  n--;
  for(int j = 0; j < n; j++)
    if(N[j] > n) N[j]--;
  for(int i = 0; i < m; i++)
    if(B[i] > n) B[i]--;

  for(int j = 0; j < n; j++) c[j] = 0;
  v = v_old;

  for(int j = 0; j < n; j++) {
    bool ok = false;
    for(int jj = 0; jj < n; jj++) {
      if(j == N[jj]) {
        c[jj] += c_old[j];
        ok = true;
        break;
      }
    }
    if(ok) continue;
    for(int i = 0; i < m; i++) {
      if(j == B[i]) {
        for(int jj = 0; jj < n; jj++) {
          c[jj] += c_old[j] * A[i][jj];
        }
        v += c_old[j] * b[i];
        break;
      }
    }
  }

  return 0;
}

// Runs the simplex algorithm to optimise the LP.
// Returns a vector of -1s if unbounded, -2s if infeasible.
// pair<vector<double>, double> simplex()
// {
//     if (initialise_simplex() == -1)
//     {
//         return make_pair(vector<double>(n + m, -2), INFINITY);
//     }

//     int code;
//     while (!(code = iterate_simplex()));

//     if (code == -1) return make_pair(vector<double>(n + m, -1), INFINITY);

//     vector<double> ret;
//     ret.resize(n + m);
//     for (int j=0;j<n;j++)
//     {
//         ret[N[j]] = 0;
//     }
//     for (int i=0;i<m;i++)
//     {
//         ret[B[i]] = b[i];
//     }

//     return make_pair(ret, v);
// }

pair<vector<double>, double> simplex(int nInput, int mInput, vector<double> AInput,
                                     vector<double> bInput, vector<double> cInput, double vInput)
{
  n = nInput;
  m = mInput;

  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      A[i][j] = AInput[i * n + j];
    }
  }

  for(int i = 0; i < n; ++i) {
    c[i] = cInput[i];
  }

  for(int i = 0; i < m; ++i) {
    b[i] = bInput[i];
  }

  v = vInput;

  if(initialise_simplex() == -1) {
    return make_pair(vector<double>(n + m, -2), INFINITY);
  }

  int code;
  while(!(code = iterate_simplex()))
    ;

  if(code == -1) return make_pair(vector<double>(n + m, -1), INFINITY);

  vector<double> ret;
  ret.resize(n + m);
  for(int j = 0; j < n; j++) {
    ret[N[j]] = 0;
  }
  for(int i = 0; i < m; i++) {
    ret[B[i]] = b[i];
  }

  return make_pair(ret, v);
}