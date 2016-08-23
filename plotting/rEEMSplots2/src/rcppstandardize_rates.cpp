
// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

const double pi_180 = M_PI / 180.0;
const double Earth_radiusX2 = 2.0 * 6378137.0;

// Compute the effective resistance matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd resistance_distance(const Eigen::MatrixXd &M) {
    Eigen::MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum(); Hinv.array() += 1.0;
    Eigen::MatrixXd H = Hinv.inverse();
    Eigen::MatrixXd R = - 2.0 * H;
    Eigen::VectorXd u = Eigen::VectorXd::Ones(M.rows());
    R.noalias() += H.diagonal() * u.transpose();
    R.noalias() += u * H.diagonal().transpose();
    return R;
}
// Squared Euclidean distance: Taking the square root is not necessary
// because EEMS uses the distances to find closest points. For example,
//   pairwise_distance(X,Y).col(0).minCoeff( &closest )
// finds the row/point in X that is closest to the first row/point in Y
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y) {
    return(  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
                 + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
                 - 2.0 * X * Y.transpose() );
}
// Great circle distance, using the haversine formula
Eigen::MatrixXd greatcirc_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y) {
    int nr = X.rows();
    int nc = Y.rows();
    // Convert from degrees to radians
    Eigen::ArrayXXd lon1 = X.col(0).replicate(1, nc).array() * pi_180;
    Eigen::ArrayXXd lat1 = X.col(1).replicate(1, nc).array() * pi_180;
    Eigen::ArrayXXd lon2 = Y.col(0).transpose().replicate(nr, 1).array() * pi_180;
    Eigen::ArrayXXd lat2 = Y.col(1).transpose().replicate(nr, 1).array() * pi_180;
    // The haversine function is hav(theta) = (1 - cos(theta)) / 2
    Eigen::ArrayXXd hav_lon = 0.5 * (1.0 - (lon2 - lon1).cos());
    Eigen::ArrayXXd hav_lat = 0.5 * (1.0 - (lat2 - lat1).cos());
    Eigen::ArrayXXd h = hav_lat + hav_lon * lat1.cos() * lat2.cos();
    // The great circle distance d is given by d = 2 * radius * arcsine( sqrt(h) )
    Eigen::MatrixXd d = (h < 1.0).select(h.sqrt(), 1.0).asin(); // Avoid numerical issues by ensuring h <= 1.0
    return (Earth_radiusX2 * d);
}
// Compute pairwise distances between the rows of X and the rows of Y.
// Choose either Euclidean or great circle distance
Eigen::MatrixXd pairwise_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y, const std::string &distm) {
    if (!distm.compare("greatcirc")) {
        return (greatcirc_dist(X, Y));
    } else {
        return (euclidean_dist(X, Y));
    }
}
// Compute one contour, by filling in each of the pixels/marks
Eigen::VectorXd compute_contour_values(const Eigen::MatrixXd &marks,
				       const Eigen::VectorXd &now_rates, const Eigen::MatrixXd &now_seeds,
				       const std::string &distm) {
    int nmrks = marks.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Zero(nmrks);
    Eigen::MatrixXd distances = pairwise_dist(marks, now_seeds, distm);
    for ( int row = 0, closest = 0 ; row < nmrks ; row++ ) {
        distances.row(row).minCoeff( &closest );
        values(row) = now_rates( closest );
    }
    return values;
}
// Compute the average contours, by calling compute_contour_vals repeatedly.
// Do this for the m rates, N rates and N*m rates, at the same time.
//
// [[Rcpp::export]]
void rcppcompute_Nm(const Eigen::VectorXd &mtiles, const Eigen::VectorXd &mrates,
		    const Eigen::MatrixXd &mseeds, const Eigen::VectorXd &qtiles,
		    const Eigen::VectorXd &qrates, const Eigen::MatrixXd &qseeds,
		    const double grandmean_m, const double grandmean_N, const double grandmean_Nm,
		    const Eigen::MatrixXd &marks, const std::string &distm, const double weight,
		    Eigen::VectorXd &means_m , Eigen::VectorXd &prGTx_m , Eigen::VectorXd &prLTx_m,
		    Eigen::VectorXd &means_N , Eigen::VectorXd &prGTx_N , Eigen::VectorXd &prLTx_N,
		    Eigen::VectorXd &means_Nm, Eigen::VectorXd &prGTx_Nm, Eigen::VectorXd &prLTx_Nm) {
    int nmrks = marks.rows();
    Eigen::VectorXd values_m  = Eigen::VectorXd::Zero(nmrks);
    Eigen::VectorXd values_N  = Eigen::VectorXd::Zero(nmrks);
    Eigen::VectorXd values_Nm = Eigen::VectorXd::Zero(nmrks);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(nmrks);
    // Vectors to store the migration rates m, population sizes N and their product N * m.
    means_m  = Eigen::VectorXd::Zero(nmrks);
    means_N  = Eigen::VectorXd::Zero(nmrks);
    means_Nm = Eigen::VectorXd::Zero(nmrks);
    // Vectors to store the number of times (m > X) and (m < X), and so on.
    prGTx_m  = Eigen::VectorXd::Zero(nmrks);
    prLTx_m  = Eigen::VectorXd::Zero(nmrks);
    prGTx_N  = Eigen::VectorXd::Zero(nmrks);
    prLTx_N  = Eigen::VectorXd::Zero(nmrks);
    prGTx_Nm = Eigen::VectorXd::Zero(nmrks);
    prLTx_Nm = Eigen::VectorXd::Zero(nmrks);
    int mstart = 0, qstart = 0;
    for ( int iter = 0 ; iter < mtiles.size() ; iter++ ) {
        int now_mtiles = (int)mtiles(iter);
	int now_qtiles = (int)qtiles(iter);
        Eigen::VectorXd now_mrates = mrates.segment(mstart, now_mtiles);
        Eigen::VectorXd now_qrates = qrates.segment(qstart, now_qtiles);
        Eigen::MatrixXd now_mseeds = mseeds.block(mstart, 0, now_mtiles, 2);
        Eigen::MatrixXd now_qseeds = qseeds.block(qstart, 0, now_qtiles, 2);
        values_m  = compute_contour_values(marks, now_mrates, now_mseeds, distm);
        values_N  = compute_contour_values(marks, now_qrates, now_qseeds, distm);
	values_Nm = values_N.cwiseProduct(values_m);
        means_m  += weight * values_m;
	means_N  += weight * values_N;
	means_Nm += weight * values_Nm;
        prGTx_m  += weight * (values_m.array()  > grandmean_m ).select(ones, 0.0);
        prLTx_m  += weight * (values_m.array()  < grandmean_m ).select(ones, 0.0);
        prGTx_N  += weight * (values_N.array()  > grandmean_N ).select(ones, 0.0);
        prLTx_N  += weight * (values_N.array()  < grandmean_N ).select(ones, 0.0);
        prGTx_Nm += weight * (values_Nm.array() > grandmean_Nm).select(ones, 0.0);
        prLTx_Nm += weight * (values_Nm.array() < grandmean_Nm).select(ones, 0.0);
        mstart += now_mtiles;
	qstart += now_qtiles;
    }
}
