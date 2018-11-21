
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

// [[Rcpp::export]]
SEXP calc_M2(int nSize, int nSpecies, arma::mat N, arma::vec scLinf, arma::mat ration, arma::mat wgt, arma::cube suitability ,double phiMin, double otherFood ) {

  arma::mat M2 = arma::zeros(nSize,nSpecies);

  for (int msp = 0; msp<nSpecies; msp++) { //prey species
    for (int nsz = 0; nsz < scLinf(msp); nsz++){ // prey size
      for (int isp=0; isp<nSpecies; isp++) { //predator species
        for (int jsz = 0; jsz < scLinf(isp); jsz++) { //predator size
          int index = (isp * nSize) + jsz;
          double numerator = ration(jsz,isp)*N(jsz,isp)*suitability(index,nsz,msp);
          // additional loops to sum over suitabilities
          arma::mat tempSuit = arma::zeros(nSize,nSpecies);
          for (int lsp=0 ; lsp < nSpecies ; lsp++) {
            for (int ksz=0 ; ksz < scLinf(lsp) ; ksz++) {
              tempSuit(ksz,lsp) = suitability(index,ksz,lsp);
            }
          }
          double denominator = accu(tempSuit % wgt % N) + otherFood;

          M2(nsz,msp) += numerator/denominator;

        }

      }
    }

  }

  M2 = 10 * nSize;
  Rcpp::Rcout << "This is nSize = " << nSize << std::endl;

  return Rcpp::wrap(M2);
}


// this checks to see if user has hit ctrl + c to abort script
// Rcpp::checkUserInterrupt()

// calc_M2 <- function(nSize,nSpecies,N,ration,suitability,phiMin,otherFood){
//
//   M2 <- matrix(data=0,nrow = nSize,ncol = nSpecies)
//   for (msp in 1:nSpecies) {  # prey species
//     for (nsz in 1:ration$scLinf[msp]) { # prey size class
//       for (isp in 1:nSpecies) { # predator
//         for (jsz in 1:ration$scLinf[isp]) { # predator size
//           numerator <- ration$ration[jsz,isp]*N[jsz,isp]*suitability[jsz,isp,nsz,msp]
//           denominator <- sum(suitability[jsz,isp,,] * ration$wgt * N) + otherFood
//           M2[nsz,msp] <- M2[nsz,msp] + numerator/denominator
//         }
//       }
//     }
//   }
//   M2 <- M2*phiMin
//
//     return(M2)
// }
