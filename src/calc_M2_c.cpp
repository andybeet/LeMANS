// [[Rccp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//' Calculates predation mortality, M2 (C++)
//'
//' Calculates the predation mortality for each species in each size class
//'
//'@param nSize Number of size class intervals species can grow through
//'@param nSpecies Number of species in the model
//'@param N A matrix (nSize x nSpecies) of abundance (number of individuals)
//'@param scLinf A vector (nSpecies) of the size class at which each species reaches L_inf (maximum length)
//'@param ration A matrix (nSize x nSpecies) of (growth in time interval)/growth efficiency values. See \code{\link{calc_ration}}
//'@param wgt A matrix (nSize x nSpecies) of species weight at the mid point of each size class (Units: grams). See \code{\link{calc_ration}}
//'@param suitability 3D array of predator size preference for prey size. See \code{\link{calc_sizepref_suitability}}
//'@param phiMin Scalar. Model timestep (years). See \code{\link{calc_phi}}
//'@param otherFood Scalar. Amount of other food available not accounted for in the model (g)
//'
//'@return A matrix is returned
//'
//'   \item{M2}{A matrix (nSize x nSpecies) of M2 (predation mortality) values.
//'    Note: M2 = 0 for size classes in which a species is not preyed upon.}
//'
//'@seealso \code{\link{calc_sizepref_suitability}}, \code{\link{calc_ration}} \code{\link{calc_phi}}
//'@section References:
//'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
//'
//'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
//'alc_
//' @export
// [[Rcpp::export]]
arma::mat calc_M2_c(int nSize, int nSpecies, arma::mat N, arma::vec scLinf, arma::mat ration, arma::mat wgt,
                  arma::cube suitability, double phiMin, double otherFood) {

  arma::mat M2 = arma::zeros(nSize,nSpecies);


  for (int msp = 0; msp<nSpecies; msp++) { //prey species
    // this checks to see if user has hit ctrl + c to abort script. Add it to the body of code
    Rcpp::checkUserInterrupt();
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

          M2(nsz,msp) += (numerator/denominator);

        }
      }
    }
  }

  M2 = M2*phiMin;
  return M2;
}


//std::ofstream ofs;
//ofs.open ("test.txt", std::ofstream::out | std::ofstream::app);

//          ofs << msp<< nsz<< isp<< jsz <<" "<<numerator << "\n";
//          ofs << msp<< nsz<< isp<< jsz <<" "<<denominator << "\n";
//ofs.close();
