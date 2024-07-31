#include <RcppArmadillo.h>
#include <omp.h>
#include <progress.hpp>


using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
// Hamming Dist between string of the same length and where the Ns in the reference do not count
arma::mat hammingDistAll(const arma::mat& m, const arma::mat& m2, int ncores=1,bool verbose=true,int th=6)
{
  arma::mat hdist(m.n_rows, m2.n_rows);
  hdist.fill(0);
  Progress p(m.n_rows*m2.n_rows, verbose);
#pragma omp parallel for num_threads(ncores) shared(hdist)
  for(int i=0;i<m.n_rows;i++) {
    for(int j=0;j<m2.n_rows;j++)
    {
      int mismatchs = 0;
      for(int k=0;k<std::min(m.n_cols,m2.n_cols) && mismatchs<=th;k++) {
        if( m(i,k) != m2(j,k) ) {
          mismatchs++;
        }
      }
      hdist(i,j) = mismatchs;
    }
  }
  return hdist;
}
