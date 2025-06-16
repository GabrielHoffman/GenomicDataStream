
/***********************************************************************
 * @file        PCAOne.h
 * @author      Gabriel Hoffman, Zilong Li
 * @email       gabriel.hoffman@mssm.edu
 * @brief       PCA One algorithm
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef PCA_ONE_H
#define PCA_ONE_H

#include <iostream>

#include "GenomicDataStream.h"
#include "Rand.hpp"

using namespace std;
using namespace arma;
using namespace gds;

struct PCA {

	PCA( Eigen::MatrixXd u):
		u(u) {}

	PCA( const Eigen::MatrixXd &u, const Eigen::MatrixXd &v, const Eigen::VectorXd &d, const vector<string> featureIds):
		u(u), v(v), d(d), featureIds(featureIds) {}

	Eigen::MatrixXd u, v;
	Eigen::VectorXd d;
	const vector<string> featureIds;
};


PCA pcaone(	const shared_ptr<GenomicDataStream> gds,
						const string &region,
						int m, 
						int k, 
						int nchunks,
						int s = 20, 
						int p = 7, 
						int B = 64, 
						int threads = 4,
						const bool verbose = true,
						const bool scaleAndCenter = true,
            ostream& outstrm = Rcpp::Rcout) {
  
  auto regions = splitRegionString( region ); 

  int n = gds->n_samples();

  const int l = k + s;
  auto randomEngine = std::default_random_engine{};
  Eigen::MatrixXd Omg = StandardNormalRandom<Eigen::MatrixXd, std::default_random_engine>(n, l, randomEngine);
  Eigen::MatrixXd Omg2 = Omg;
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero(n, l);
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero(n, l);
  Eigen::MatrixXd H(n, l), G(m, l), R(l, l), Rt(l, l);
  
  // Test: Read in Serial
  shared_ptr<GenomicDataStream> gdsStream = createFileView( gds->getParam() );
  gdsStream->setRegions( regions );
  DataChunk<Eigen::MatrixXd> ch;
  
  // Parallel part using Thread Building Blocks
  GenomicDataStreamParallel<Eigen::MatrixXd> gsp(gds->getParam(), regions, nchunks, threads);

  size_t band = ceil((double) nchunks / B);
  vector<string> featureIds;
  std::mutex pcaMutex;  

  for (int pi = 0; pi <= p; pi++) {

    if( verbose ){
      outstrm << "\rEpoch " << pi << " / " << p << "  ";
    }

    if (std::pow(2, pi) >= B) {
      // reset H1, H2 to zero
      H1.setZero();
      H2.setZero();
    }
    band = std::fmin(band * 2, nchunks);

    size_t i{1},  start{0};
    gsp.processChunks([&](const gds::DataChunk<Eigen::MatrixXd> &chunk, size_t b) {
      auto Ab = chunk.getData();

      if( scaleAndCenter ){
        standardize(Ab);
        Ab /= sqrt(n-1);    
      }
      {
        std::lock_guard<std::mutex> lock(pcaMutex);

        if( pi == 0){
          auto tmp = chunk.getInfo<VariantInfo>()->getFeatureNames();
          featureIds.insert(featureIds.end(), 
              tmp.begin(), tmp.end());
        }

        G.middleRows(start, Ab.cols()).noalias() = Ab.transpose() * Omg;
        if (i <= band / 2)
          H1.noalias() += Ab * G.middleRows(start, Ab.cols());
        else
          H2.noalias() += Ab * G.middleRows(start, Ab.cols());
        bool adjacent =
          (pi > 0 && (b + 1) == std::pow(2, pi - 1) && std::pow(2, pi) < B);
        if((!((b + 1) < band && !adjacent)) && ((i == band) || (i == band / 2) || adjacent)){
          H = H1 + H2;
          Eigen::HouseholderQR<Eigen::MatrixXd> qr(H);
          Omg.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(n, l);
          flipOmg(Omg2, Omg);
          if (i == band) {
            H1.setZero();
            i = 0;
          } else {
            H2.setZero();
          }
        }
        start += Ab.cols();
        i++;
      }
    }); // end processChunks
  } // end for

  // get USV
  if( verbose ){
		outstrm << "\rFinal decompositions";
	}
  {
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(G);
    R.noalias() = Eigen::MatrixXd::Identity(l, m) * qr.matrixQR().triangularView<Eigen::Upper>();
    G.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(m, l);
  }
  {
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(G);
    Rt.noalias() = Eigen::MatrixXd::Identity(l, m) * qr.matrixQR().triangularView<Eigen::Upper>();
    G.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(m, l);
  }

  R = Rt * R;

  Eigen::MatrixXd out = R.transpose().fullPivHouseholderQr().solve(H.transpose());
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(out, Eigen::ComputeThinU | Eigen::ComputeThinV);

  if( verbose ){
    outstrm << "\rCompleted            " << std::endl;
  }

  return PCA(svd.matrixV().leftCols(k),
			  	G * svd.matrixU().leftCols(k),
			  	svd.singularValues().head(k),
			  	featureIds);
}

#endif
