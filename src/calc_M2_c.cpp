// calculates M2 using c++
//
// @param loads of parameters. document later
//@export

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <fstream>


// this checks to see if user has hit ctrl + c to abort script. Add it to the body of code
// Rcpp::checkUserInterrupt()

using namespace std;

// [[Rcpp::export]]
arma::mat calc_M2_c(int nSize, int nSpecies, arma::mat N, arma::vec scLinf, arma::mat ration, arma::mat wgt,
                  arma::cube suitability, double phiMin, double otherFood) {

  arma::mat M2 = arma::zeros(nSize,nSpecies);


  for (int msp = 0; msp<nSpecies; msp++) { //prey species
    for (int nsz = 0; nsz < scLinf(msp); nsz++){ // prey size
      for (int isp=0; isp<nSpecies; isp++) { //predator species
        for (int jsz = 0; jsz < scLinf(isp); jsz++) { //predator size
          int index = (isp * nSize) + jsz;
          double numerator = ration(jsz,isp)*N(jsz,isp)*suitability(index,nsz,msp);

          //Rcpp::Rcout << "numerator = " << numerator << std::endl;
          // additional loops to sum over suitabilities
          arma::mat tempSuit = arma::zeros(nSize,nSpecies);
          for (int lsp=0 ; lsp < nSpecies ; lsp++) {
            for (int ksz=0 ; ksz < scLinf(lsp) ; ksz++) {
              tempSuit(ksz,lsp) = suitability(index,ksz,lsp);
            }
          }
          double denominator = accu(tempSuit % wgt % N) + otherFood;
          //Rcpp::Rcout << "denominator = " << denominator << std::endl;

          M2(nsz,msp) += (numerator/denominator);

        }
      }
    }
  }


  M2 = M2/phiMin;
  //M2 = M2/4.58;

  //Rcpp::Rcout << "This is dim of M2 = " << M2.size() << std::endl;

  return M2;
}


//std::ofstream ofs;
//ofs.open ("test.txt", std::ofstream::out | std::ofstream::app);

//          ofs << msp<< nsz<< isp<< jsz <<" "<<numerator << "\n";
//          ofs << msp<< nsz<< isp<< jsz <<" "<<denominator << "\n";
//ofs.close();
