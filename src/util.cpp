
#include "util.hpp"

// NEED TO DELETE THE SIGMA SCALE STUFF???

extern string dist_metric;

Params::Params( ) { }
Params::~Params( ) { }
Params::Params(const string &params_file, const long seed_from_command_line) {
    try {
        po::options_description eems_options("EEMS options from parameter file");
        eems_options.add_options()
        ("seed", po::value<long>(&seed)->default_value(seed_from_command_line), "Random seed")
        ("datapath", po::value<string>(&datapath)->required(), "Path to coord/sims/outer files")
        ("mcmcpath", po::value<string>(&mcmcpath)->required(), "Path to output directory")
        ("prevpath", po::value<string>(&prevpath)->default_value(""), "Path to previous output directory")
        ("gridpath", po::value<string>(&gridpath)->default_value(""), "Path to demes/edges/ipmap files")
        ("nIndiv", po::value<int>(&nIndiv)->required(), "nIndiv")
        //("genomeSize", po::value<double>(&genomeSize)->default_value(3000000000), "genomeSize")
        //("cutOff", po::value<double>(&cutOff)->default_value(3000000), "cutOff")
        ("nDemes", po::value<int>(&nDemes)->default_value(1), "nDemes")
        ("diploid", po::value<bool>(&diploid)->default_value(true), "diploid")
        ("distance", po::value<string>(&distance)->default_value("euclidean"), "distance")
        ("numMCMCIter", po::value<int>(&numMCMCIter)->default_value(1), "numMCMCIter")
        ("numBurnIter", po::value<int>(&numBurnIter)->default_value(0), "numBurnIter")
        ("numThinIter", po::value<int>(&numThinIter)->default_value(0), "numThinIter")
        ("mSeedsProposalS2", po::value<double>(&mSeedsProposalS2)->default_value(0.01), "mSeedsProposalS2")
        ("qSeedsProposalS2", po::value<double>(&qSeedsProposalS2)->default_value(0.1), "qSeedsProposalS2")
        ("mEffctProposalS2", po::value<double>(&mEffctProposalS2)->default_value(0.1), "mEffctProposalS2")
        ("qEffctProposalS2", po::value<double>(&qEffctProposalS2)->default_value(0.001), "qEffctProposalS2")
        ("mrateMuProposalS2", po::value<double>(&mrateMuProposalS2)->default_value(0.01), "mrateMuProposalS2")
        ("qrateMuProposalS2", po::value<double>(&qrateMuProposalS2)->default_value(0.01), "qrateMuProposalS2")
        ("qVoronoiPr", po::value<double>(&qVoronoiPr)->default_value(0.5), "qVoronoiPr")
        ("mrateShape", po::value<double>(&mrateShape_2)->default_value(0.001), "mrateShape")
        ("qrateShape", po::value<double>(&qrateShape_2)->default_value(0.001), "qrateShape")
        ("sigmaShape", po::value<double>(&sigmaShape_2)->default_value(0.001), "sigmaShape")
        ("qrateScale", po::value<double>(&qrateScale_2)->default_value(1.0), "qrateScale")
        ("mrateScale", po::value<double>(&mrateScale_2)->default_value(1.0), "mrateScale")
        ("sigmaScale", po::value<double>(&sigmaScale_2)->default_value(1.0), "sigmaScale")
        ("negBiProb", po::value<double>(&negBiProb)->default_value(0.67), "negBiProb")
        ("negBiSize", po::value<int>(&negBiSize)->default_value(10), "negBiSize") ;
        ifstream instrm(params_file.c_str());
        po::variables_map vm;
        po::store(po::parse_config_file(instrm,eems_options,true),vm);
        po::notify(vm);
        instrm.close();
    } catch(exception& e) {
        cerr << "[EEMS::Params] Error parsing input parameters in " << params_file << ": " << endl;
        cerr << e.what() << endl; exit(1);
    }
    mrateShape_2 /= 2.0;
    qrateShape_2 /= 2.0;
    sigmaShape_2 /= 2.0;
    mrateScale_2 /= 2.0;
    qrateScale_2 /= 2.0;
    sigmaScale_2 /= 2.0;
    
    dfmin = nIndiv;
    dfmax = 1e6;
    testing = false;

    
    
    // let's assume a maximum population size of 2N = 500
    // and maximum migration rate of m = 0.1. Remember, rates are paramterized on the log scale
    mrateMuUpperBound = -0.301; // log10(0.5)
    qrateMuUpperBound = -2.3; // log10(0.005)
    
    mrateMuLowerBound = -10.0;
    qrateMuLowerBound = -10.0;

    // Ensure that mrateMuUpperBound + mEffectHalfInterval <= 0 so rates are between 0 and 1.
    mEffctHalfInterval = 0.301;
    qEffctHalfInterval = 2.0;
}
ostream& operator<<(ostream& out, const Params& params) {
    out << "               datapath = " << params.datapath << endl
    << "               mcmcpath = " << params.mcmcpath << endl
    << "               prevpath = " << params.prevpath << endl
    << "               gridpath = " << params.gridpath << endl
    << "               distance = " << params.distance << endl
    << "                 nIndiv = " << params.nIndiv << endl
    //<< "             genomeSize = " << params.genomeSize << endl
    //<< "                 cutOff = " << params.cutOff << endl
    << "                 nDemes = " << params.nDemes << endl
    << "                   seed = " << params.seed << endl
    << "            numMCMCIter = " << params.numMCMCIter << endl
    << "            numBurnIter = " << params.numBurnIter << endl
    << "            numThinIter = " << params.numThinIter << endl
    << "              negBiSize = " << params.negBiSize << endl
    << fixed << setprecision(10)
    << "              negBiProb = " << params.negBiProb << endl
    << "             qVoronoiPr = " << params.qVoronoiPr << endl
    << "             mrateShape = " << 2.0*params.mrateShape_2 << endl
    << "             qrateShape = " << 2.0*params.qrateShape_2 << endl
    << "             sigmaShape = " << 2.0*params.sigmaShape_2 << endl
    << "             qrateScale = " << 2.0*params.qrateScale_2 << endl
    << "             mrateScale = " << 2.0*params.mrateScale_2 << endl
    << "             sigmaScale = " << 2.0*params.sigmaScale_2 << endl
    << "       mSeedsProposalS2 = " << params.mSeedsProposalS2 << endl
    << "       qSeedsProposalS2 = " << params.qSeedsProposalS2 << endl
    << "       mEffctProposalS2 = " << params.mEffctProposalS2 << endl
    << "       qEffctProposalS2 = " << params.qEffctProposalS2 << endl
    << "      mrateMuProposalS2 = " << params.mrateMuProposalS2 << endl
    << "      qrateMuProposalS2 = " << params.qrateMuProposalS2 << endl;
    return out;
}
bool Params::check_input_params( ) const {
    bool error = false;
    boost::filesystem::path mcmcdir(mcmcpath.c_str());
    boost::filesystem::path prevdir(prevpath.c_str());
    cerr << "Using Boost " << BOOST_LIB_VERSION
    << " and Eigen " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << endl
    << "  EEMS was tested with Boost 1_57 and Eigen 3.2.4" << endl << endl;
    if (!numeric_limits<double>::has_infinity || !numeric_limits<double>::infinity()) {
        cerr << "  Infinity not supported on this platform" << endl;
        error = true;
    }
    if (!boost::filesystem::create_directory(mcmcdir) && !exists(mcmcdir)) {
        cerr << "  Failed to create output directory " << mcmcpath << endl;
        error = true;
    }
    if (!(!distance.compare("euclidean") || !distance.compare("greatcirc"))) {
        cerr << "  Choose either 'euclidean' or 'greatcirc' distance metric" << endl;
        error = true;
    }
    if (!boost::filesystem::exists(datapath + ".coord") ||
        !boost::filesystem::exists(datapath + ".sims") ||
        !boost::filesystem::exists(datapath + ".outer")) {
        cerr << "  Failed to find input files " << datapath << ".coord/sims/outer" << endl;
        error = true;
    }
    if (!gridpath.empty() &&
        (!boost::filesystem::exists(gridpath + ".demes") ||
         !boost::filesystem::exists(gridpath + ".edges") ||
         !boost::filesystem::exists(gridpath + ".ipmap"))) {
            // Path to input population grid is specified
            // but one of the required files is missing
            cerr << "  Failed to find graph files " << gridpath << ".demes/edges/ipmap" << endl;
            error = true;
        }
    if (!prevpath.empty() && !exists(prevdir)) {
        cerr << "  Failed to find directory " << prevpath << " to previous EEMS output" << endl;
        error = true;
    }
    if (!(mSeedsProposalS2>0) || !(mEffctProposalS2>0) || !(mrateMuProposalS2>0) ||
        !(qSeedsProposalS2>0) || !(qEffctProposalS2>0) || !(qrateMuProposalS2>0)) {
        cerr << "  Choose positive variance parameters for the proposal distributions:" << endl
        << "  mrateMuProposalS2 = " << mrateMuProposalS2 << ", qrateMuProposalS2 = " << qrateMuProposalS2 << endl
        << "   mSeedsProposalS2 = " << mSeedsProposalS2 << ", mEffctProposalS2 = " << mEffctProposalS2 << endl
        << "   qSeedsProposalS2 = " << qSeedsProposalS2 << ", qEffctProposalS2 = " << qEffctProposalS2 << endl;
        error = true;
    }
    
    /*if (genomeSize > 3.3e9){
        cerr << "  Error with genome size: " << endl
        << " genomeSize = " << genomeSize << endl;
        error = true;
    }
     
    
    if (cutOff > genomeSize){
        cerr << "  Error with IBD cut off: " << endl
        << " cutOff = " << cutOff << endl;
        error = true;
    }
     */
    
    if (!(numMCMCIter>0) || !(numBurnIter>=0) || !(numThinIter>=0) || !(numMCMCIter>numBurnIter) ||
        !(numMCMCIter>numBurnIter+numThinIter)) {
        cerr << "  Error with the MCMC parameters:" << endl
        << "  numMCMCIter = " << numMCMCIter << ", numBurnIter = " << numBurnIter << ", numThinIter " << numThinIter << endl;
        error = true;
    }
    if (!(qrateShape_2>0) || !(mrateShape_2>0) || !(sigmaShape_2) ||
        !(qrateScale_2>0) || !(mrateScale_2>0) || !(sigmaScale_2)) {
        cerr << "  Error with the Inverse Gamma hyperparameters:" << endl
        << "  qrateShape = " << 2.0*qrateShape_2 << ", qrateScale = " << 2.0*qrateScale_2 << endl
        << "  mrateShape = " << 2.0*mrateShape_2 << ", mrateScale = " << 2.0*mrateScale_2 << endl
        << "  sigmaShape = " << 2.0*sigmaShape_2 << ", sigmaScale = " << 2.0*sigmaScale_2 << endl;
        error = true;
    }
    if (!(negBiSize>0) || !( (negBiProb>0) && (negBiProb<1) )) {
        cerr << "  Error with the Negative Binomial hyperparameters:" << endl
        << "  negBiSize = " << negBiSize << ", negBiProb = " << negBiProb << endl;
        error = true;
    }
    return(error);
}
VectorXd split(const string &line) {
    istringstream in(line);
    vector<double> numbers;
    double number;
    while (!in.eof()) { in >> number;
        if (!in.fail()) { numbers.push_back(number); } else { break; }
    }
    if (in.fail()||in.bad()) { return (VectorXd::Zero(0)); }
    // Since we don't know how many numbers there are,
    // first we store the data in a std::vector of doubles,
    // and then typecast it to Eigen::VectorXd
    return (VectorXd::Map(&numbers[0],numbers.size()));
}
bool isposdef(const MatrixXd &A) {
    SelfAdjointEigenSolver<MatrixXd> eig(A,EigenvaluesOnly);
    double minval = eig.eigenvalues().minCoeff();
    return (minval>0);
}

double logdet(const MatrixXd &A) {
    return (A.selfadjointView<Lower>().ldlt().vectorD().array().log().sum());
}
double pseudologdet(const MatrixXd &A, const int rank) {
    SelfAdjointEigenSolver<MatrixXd> x(A);
    return (x.eigenvalues().reverse().array().head(rank).log().sum());
}

double poisln(const MatrixXd &expectedIBD, const MatrixXd &observedIBD, const MatrixXd &cMatrix){
    double ll = 0;
    double epsilon = 1e-8;
    int n = expectedIBD.rows();
    double lamda;
    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++){
            if (expectedIBD(i,j) < epsilon){
                lamda = epsilon;
            } else{
                lamda = expectedIBD(i,j);
            }
            ll += observedIBD(i,j)*log(lamda)-cMatrix(i,j)*lamda;
        }
    }
    return(ll);
}

double mvgammaln(const double a, const int p) {
    double val = 0.25*log_pi*p*(p-1);
    for (int i = 0 ; i < p ; i++ ) {
        val += lgamma(a-0.5*i);
    }
    return (val);
}
/*
 Squared Euclidean distance: Taking the square root is not necessary
 because EEMS uses the distances to find closest points. For example,
 pairwise_distance(X,Y).col(0).minCoeff( &closest )
 finds the row/point in X that is closest to the first row/point in Y
 */
MatrixXd euclidean_dist(const MatrixXd &X, const MatrixXd &Y) {
    return (  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
            + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
            - 2.0*X*Y.transpose() );
}
// Great circle distance, up to a constant of proportionality equal to 2*R
// where R is the earth's radius
MatrixXd greatcirc_dist(const MatrixXd &X, const MatrixXd &Y) {
    int nr = X.rows();
    int nc = Y.rows();
    MatrixXd lon1 = X.col(0).replicate(1,nc) * pi_180;
    MatrixXd lat1 = X.col(1).replicate(1,nc) * pi_180;
    MatrixXd lon2 = Y.col(0).transpose().replicate(nr,1) * pi_180;
    MatrixXd lat2 = Y.col(1).transpose().replicate(nr,1) * pi_180;
    MatrixXd dlon = 0.5*(lon2 - lon1);
    MatrixXd dlat = 0.5*(lat2 - lat1);
    MatrixXd a = dlat.array().sin().square().matrix() +
    (dlon.array().sin().square() * lat1.array().cos() * lat2.array().cos()).matrix();
    MatrixXd c = (a.array()<1.0).select(a.array().sqrt(),MatrixXd::Ones(nr,nc)).array().asin();
    return (c); // Instead of (2*R*c) where R = 6378137 is the Earth's radius.
}
// Compute pairwise distances between the rows of X and the rows of Y.
// Choose either Euclidean or great circle distance
MatrixXd pairwise_distance(const MatrixXd &X, const MatrixXd &Y) {
    if (!dist_metric.compare("greatcirc")) {
        return (greatcirc_dist(X,Y));
    } else {
        return (euclidean_dist(X,Y));
    }
}
// Read a matrix, with unknown dimensions, from a text file
// Return an empty matrix (0 rows, 0 columns) if there is an error
MatrixXd readMatrixXd(const string &filename) {
    ifstream instrm(filename.c_str(), ios::in);
    if (!instrm.is_open( )) { return (MatrixXd::Zero(0,0)); }
    string line; getline(instrm,line);
    boost::algorithm::trim(line);
    // Split the first line into numbers
    VectorXd row = split(line);
    // This tells us the number of columns
    int cols = row.size();
    if (!cols) { return (MatrixXd::Zero(0,0)); }
    MatrixXd mat(1,cols); mat.row(0) = row;
    while (getline(instrm,line)) {
        boost::algorithm::trim(line);
        row = split(line);
        if (row.size()==cols) {
            // Resize the matrix to add each row (done once at the start of EEMS, so okay)
            insertRow(mat,row);
        } else {
            return (MatrixXd::Zero(0,0));
        }
    }
    return(mat);
}
double trace_AxB(const MatrixXd &A, const MatrixXd &B) {
    return (A.cwiseProduct(B).sum());
}
// If there are two draws and the first has two tiles and the second -- three tiles,
// then sizes = c(2,3) and array = c(m_{1t_1},m_{1t_2},m_{2t_1},m_{2t_2},m_{2t_3})
bool dlmcell(const string &filename, const VectorXd &sizes, const vector<double> &array) {
    bool error = false;
    if (array.size()!=sizes.sum()) { error = true; return(error); }
    ofstream out(filename.c_str(),ofstream::out);
    if (!out.is_open()) { error = true; return(error); }
    vector<double>::const_iterator it = array.begin();
    for ( int i = 0 ; i < sizes.size() ; i++ ) {
        for ( int j = 0 ; j < sizes(i) ; j++ ) {
            out << fixed << setprecision(14) << *it << " "; it++;
        }
        out << endl;
    }
    out.close( );
    return(error);
}
void removeElem(VectorXd &vec, const int elemToRemove)
{
    int elems = vec.size() - 1;
    if (elemToRemove <= elems) {
        // Copy the last element into the element to remove
        if (elemToRemove < elems) {
            vec(elemToRemove) = vec(elems);
        }
        vec.conservativeResize(elems);
    }
}
void removeRow(MatrixXd &mat, const int rowToRemove)
{
    int rows = mat.rows() - 1;
    int cols = mat.cols();
    if( rowToRemove <= rows ) {
        // Copy the last row into the row to remove
        if ( rowToRemove < rows ) {
            mat.row(rowToRemove) = mat.row(rows);
        }
        mat.conservativeResize(rows,cols);
    }
}
void insertElem(VectorXd &vec, const double &elem)
{
    int elems = vec.size();
    vec.noalias() = (VectorXd(elems+1) << vec, elem).finished();
}
void insertRow(MatrixXd &mat, const VectorXd &row)
{
    int rows = mat.rows();
    int cols = mat.cols();
    mat.noalias() = (MatrixXd(rows+1,cols) << mat, row.transpose()).finished();
}
/*
 Log probability mass function of the negative binomial distribution, with
 parameters size (number of failures) and prob (probability of success)
 * up to a constant of proportionality that depends on size and prob,
 but not on the random variable k (number of successes)
 Also -- log pmf of the zero-truncated negative binomial, with the same parameters,
 as long as k>0, since truncating only changes the constant of proportionality
 Uses gamma(n) = (n-1)!
 */
double dnegbinln(const int k, const int size, const double prob) {
    double pln = -Inf;
    if ( (k>=0) && (size>0) && (prob>0) && (prob<1) ) {
        pln = lgamma(size+k) - lgamma(k+1) + k*log(prob);
    }
    return (pln);
}
/*
 Log probability density function of the inverse gamma distribution,
 with parameters shape and scale
 * Up to a constant of proportionality that depends on the shape and the scale,
 but not on the random variable x
 */
double dinvgamln(const double x, const double shape, const double scale) {
    double pln = -Inf;
    if ( (x>0.0) && (shape>0.0) && (scale>0.0) ) {
        pln = - (shape+1.0)*log(x) - scale/x;
    }
    return (pln);
}
/*
 Multivariate normal probability density function, on the log scale,
 up to a constant of proportionality
 */
double dmvnormln(const VectorXd &x, const VectorXd &mu, const MatrixXd &sigma) {
    double pln = -Inf;
    if ( isposdef(sigma) ) {
        LLT<MatrixXd> Chol(sigma);
        MatrixXd L = Chol.matrixL();
        pln = - L.diagonal().array().log().sum();
        pln -= 0.5 * (x-mu).transpose() * Chol.solve(x-mu);
    }
    return (pln);
}
/*
 Truncated normal probability density function, on the log scale,
 with support [0,bnd], including the normalizing constant
 */
double dtrnormln(const double x, const double mu, const double sigma2, const double bnd) {
    double pln = -Inf;
    if ( (sigma2>0) && (x>=-bnd) && (x<=bnd) ) {
        boost::math::normal pnorm(mu,sqrt(sigma2));
        pln = - 0.5 * log(sigma2) - 0.5 * (x-mu) * (x-mu) / sigma2
        - log(cdf(pnorm,bnd) - cdf(pnorm,-bnd));
    }
    return (pln);
}
VectorXd slice(const VectorXd &A, const VectorXi &I) {
    int elems = I.size();
    VectorXd B(elems);
    assert(I.minCoeff() >= 0);
    assert(I.maxCoeff() < A.size());
    for ( int i = 0 ; i < elems ; i++ ) {
        B(i) = A(I(i));
    }
    return (B);
}
MatrixXd slice(const MatrixXd &A, const VectorXi &R, const VectorXi &C) {
    int rows = R.size();
    int cols = C.size();
    MatrixXd B(rows,cols);
    assert(R.minCoeff() >= 0);
    assert(C.minCoeff() >= 0);
    assert(R.maxCoeff() < A.rows());
    assert(C.maxCoeff() < A.cols());
    for ( int i = 0 ; i < rows ; i++ ) {
        for ( int j = 0 ; j< cols ; j++ ) {
            B(i,j) = A(R(i),C(j));
        } }
    return (B);
}

double max(double a, double b){
    if (a > b){
        return(a);
    }
    return(b);
}

void getWeights(VectorXd &w, VectorXd &x){
    // REQUIRES: w and x vectors of length 30, r is recombination rate, and L (in base pairs) of cutoff.
    // MODIFIES: w and x
    // EFFECTS: x will contain the x-values telling you where to evaluate P(T_mrca = x); w will contains the weights
    // This function allows user to compute an integral by computing \sum_i P(T_mrca = x_i) * w_i
    
    x << 0.118440697736960550688, 0.3973475034735802657556, 0.8365549141880933313119, 1.437175158191620443607,
    2.200789508440616292336, 3.129448303166859096349, 4.225699164493802071261, 5.492626704368934083587,
    6.933903364122364597039, 8.553853192793023779194, 10.35753137020864105106, 12.35082332811269876439,
    14.54056869943518703492, 16.93471724415800802837, 19.54252664684054185266, 22.37481610233449499411,
    25.44429563058376261798, 28.76600031447167014762, 32.35787326932856805551, 36.24156497875364752439,
    40.44355691460364227197, 44.99678841355200250088, 49.94309754094208987181, 55.33704611950810443499,
    61.25224904369593075136, 67.79260716731075303985, 75.11420274687672563149, 83.47405073153149030595,
    93.36359463048878316735, 106.0462505962874034422;
    
    w << 0.02093564741472521761, 0.09585049298017654367, 0.18833296435057945936, 0.23281944819987904471,
    0.2060782293528492151, 0.138528960450616358, 0.07293919110208096649, 0.030605607903988887905,
    0.010333948458420042431, 0.002821608083735993584, 6.2402663742264620427E-4, 1.1168849922460852198E-4,
    1.6129719270580565631E-5, 1.87044426274856472768E-6, 1.72995513372709914535E-7, 1.26506996496773906645E-8,
    7.2352574135703022224E-10, 3.19320138447436406004E-11, 1.069761647687436460972E-12, 2.66597906070505518515E-14,
    4.82019019925788439097E-16, 6.12740480626441608041E-18, 5.26125812567892365789E-20, 2.89562589607893296815E-22,
    9.51695437836864011982E-25, 1.69046847745875738033E-27, 1.39738002075239812243E-30, 4.20697826929603166432E-34,
    2.89826026866498969507E-38, 1.411587124593531584E-43;
    
}
