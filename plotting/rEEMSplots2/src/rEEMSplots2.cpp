
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

void rcppcompute_Nm(const Eigen::VectorXd &mtiles, const Eigen::VectorXd &mrates,
		    const Eigen::MatrixXd &mseeds, const Eigen::VectorXd &qtiles,
		    const Eigen::VectorXd &qrates, const Eigen::MatrixXd &qseeds,
		    const double grandmean_m, const double grandmean_N, const double grandmean_Nm,
		    const Eigen::MatrixXd &marks, const std::string &distm, const double weight,
		    Eigen::VectorXd &means_m , Eigen::VectorXd &prGTx_m , Eigen::VectorXd &prLTx_m,
		    Eigen::VectorXd &means_N , Eigen::VectorXd &prGTx_N , Eigen::VectorXd &prLTx_N,
		    Eigen::VectorXd &means_Nm, Eigen::VectorXd &prGTx_Nm, Eigen::VectorXd &prLTx_Nm
		    );
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y);
Eigen::MatrixXd greatcirc_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y);


RcppExport SEXP rEEMSplots2__rcppcompute_Nm(SEXP mtiles_, SEXP mrates_, SEXP mseeds_,
					    SEXP qtiles_, SEXP qrates_, SEXP qseeds_,
					    SEXP grandmean_m_,
					    SEXP grandmean_N_,
					    SEXP grandmean_Nm_,
					    SEXP marks_, SEXP distm_,
					    SEXP weight_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mtiles(mtiles_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mrates(mrates_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type mseeds(mseeds_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type qtiles(qtiles_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type qrates(qrates_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type qseeds(qseeds_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type marks(marks_);
    Rcpp::traits::input_parameter< const std::string& >::type distm(distm_);
    Rcpp::traits::input_parameter< const double& >::type grandmean_m(grandmean_m_);
    Rcpp::traits::input_parameter< const double& >::type grandmean_N(grandmean_N_);
    Rcpp::traits::input_parameter< const double& >::type grandmean_Nm(grandmean_Nm_);
    Rcpp::traits::input_parameter< const double& >::type weight(weight_);
    Eigen::VectorXd means_m , prGTx_m , prLTx_m;
    Eigen::VectorXd means_N , prGTx_N , prLTx_N;
    Eigen::VectorXd means_Nm, prGTx_Nm, prLTx_Nm;
    rcppcompute_Nm(mtiles, mrates, mseeds, qtiles, qrates, qseeds,
		   grandmean_m, grandmean_N, grandmean_Nm,
		   marks, distm, weight,
		   means_m , prGTx_m , prLTx_m,
		   means_N , prGTx_N , prLTx_N,
		   means_Nm, prGTx_Nm, prLTx_Nm);
    return Rcpp::List::create(Rcpp::Named("means.m")  = means_m,
			      Rcpp::Named("prGTx.m")  = prGTx_m,
			      Rcpp::Named("prLTx.m")  = prLTx_m,
			      Rcpp::Named("means.N")  = means_N,
			      Rcpp::Named("prGTx.N")  = prGTx_N,
			      Rcpp::Named("prLTx.N")  = prLTx_N,
			      Rcpp::Named("means.Nm") = means_Nm,
			      Rcpp::Named("prGTx.Nm") = prGTx_Nm,
			      Rcpp::Named("prLTx.Nm") = prLTx_Nm);
    END_RCPP
}
RcppExport SEXP rEEMSplots2__euclidean_dist( SEXP X_, SEXP Y_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(X_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(Y_);
    __result = Rcpp::wrap(euclidean_dist(X, Y));
    return __result;
    END_RCPP    
}
RcppExport SEXP rEEMSplots2__greatcirc_dist( SEXP X_, SEXP Y_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(X_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(Y_);
    __result = Rcpp::wrap(greatcirc_dist(X, Y));
    return __result;
    END_RCPP    
}
