#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// Calculates the lower bound to the number of false hypotheses
// [[Rcpp::export]]
int findDiscSum(int                  s,   // size of I
                int                  m,   // number of all p-values
                Rcpp::NumericVector& u,   // [min(u,v), u, max(u,v)*(m-s+1)]
                Rcpp::NumericVector& v,   // [min(u,v), v, max(u,v)*(s+1)]
                Rcpp::NumericVector& cs)  // critical values for transformed p-value (1:m)
{
    int    k = 1;
    int    r = -1;
    double T = 0;

    for (int i = 1; i <= m; i++)
    {
        if (i == 1 || u[k+r+1] <= v[i-k-r])
        {
            T = T + u[k+r+1];
            r = r + 1;
        }
        else
        {
            T = T + v[i-k-r];
        }

        while (k <= std::min(s,i) && T < cs[i-1])
        {
            if (r > 0)
            {
                r = r - 1;
            }
            else
            {
                T = T + u[k+1] - v[i-k];
            }
            k = k + 1;
        }
    }

    return (s-k+1);
}
