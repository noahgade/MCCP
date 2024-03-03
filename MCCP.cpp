#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mean_nrmse(arma::mat data, arma::mat prediction) {
  arma::mat sse = arma::pow((prediction - data), 2);
  arma::mat mse = arma::mean(sse, 0);
  arma::mat norm = arma::var(data, 0, 0);
  arma::mat nrmse = arma::mean(arma::pow((mse / norm), 0.5), 1);
  return nrmse;
}

// [[Rcpp::export]]
arma::mat scaleCpp(arma::mat data) {
  float dim = data.n_cols;
  for(int i = 0; i < dim; i++) {
    data.col(i) = (data.col(i) - arma::mean(data.col(i))) / arma::stddev(data.col(i)) / 2;
  }
  return data;
}

// [[Rcpp::export]]
arma::cube genW1(float RNN_dim, float data_dim, float n_init) {
  arma::cube W1 = arma::cube(RNN_dim, 1 + data_dim + RNN_dim, n_init);
  for(int init = 0; init < n_init; init++) {
    arma::mat W1_bias = arma::randn(RNN_dim, 1);
    arma::mat W1_input = arma::randn(RNN_dim, data_dim);
    arma::sp_mat W1_hidden = arma::sprandn(RNN_dim, RNN_dim, 10 / RNN_dim);
    arma::mat W1_hid(W1_hidden);
    float eigval = abs(eig_gen(W1_hid)).max() / 0.8;
    W1.slice(init) = arma::join_rows(W1_bias, W1_input, W1_hid / eigval);
  }
  return W1;
}

// [[Rcpp::export]]
arma::cube RNN1(arma::mat data, arma::cube W1, float RSInitial, float bscale, float iscale) {
  float RNN_dim = W1.n_rows;
  float n_init = W1.n_slices;
  float data_dim = data.n_cols;
  float L = data.n_rows;
  arma::mat dat = scaleCpp(data);
  arma::cube RS = RSInitial * arma::cube(RNN_dim, n_init, 1, arma::fill::ones);
  arma::cube ResStates = arma::join_slices(RS, arma::cube(RNN_dim, n_init, L));
  for(int t = 0; t < L; t++) {
    for(int init = 0; init < n_init; init++) {
      ResStates.slice(t + 1).col(init) = tanh(bscale * W1.slice(init).col(0) + iscale * W1.slice(init).cols(1, data_dim) * trans(dat.row(t)) + W1.slice(init).cols(data_dim + 1, data_dim + RNN_dim) * ResStates.slice(t).col(init));
    }
  }
  return ResStates.slices(1, L);
}

// [[Rcpp::export]]
arma::cube WoutCalc(arma::mat data, arma::cube ResStates, float regular) {
  float RNN_dim = ResStates.n_rows;
  float n_init = ResStates.n_cols;
  float data_dim = data.n_cols;
  arma::cube Wout = arma::cube(RNN_dim, data_dim, n_init);
  for(int j = 0; j < n_init; j++) {
    arma::mat RS = ResStates.col(j);
    arma::mat WS = inv(RS * trans(RS) + regular * arma::mat(RNN_dim, RNN_dim, arma::fill::eye)) * RS * data;
    Wout.slice(j) = WS;
  }
  return Wout;
}

// [[Rcpp::export]]
arma::cube OutCalc(arma::cube Wout, arma::cube ResStates) {
  float n_init = Wout.n_slices;
  float L = ResStates.n_slices;
  float data_dim = Wout.n_cols;
  arma::cube output = arma::cube(L, data_dim, n_init);
  for(int j = 0; j < n_init; j++) {
    arma::mat WS = Wout.slice(j);
    arma::mat RS = ResStates.col(j);
    arma::mat os = trans(WS) * RS;
    output.slice(j) = trans(os);
  }
  return output;
}

// [[Rcpp::export]]
arma::vec RNNParamFit(arma::mat data, int washL, int trainL, float n_init) {
  float L = data.n_rows;
  float data_dim = data.n_cols;
  float RNN_dim = floor(0.9 * trainL);
  arma::vec bscales = arma::linspace(0.1, 0.5, 3);
  arma::vec iscales = arma::linspace(0.6, 1.4, 3);
  arma::mat errors = arma::mat(3, 3);
  float regular = 1e-6;
  for(int b = 0; b < 3; b++) {
    for(int i = 0; i < 3; i++) {
      arma::cube W1 = genW1(RNN_dim, data_dim, n_init);
      arma::cube RS = RNN1(data, W1, 0, bscales(b), iscales(i));
      arma::cube Wout = WoutCalc(data.rows(washL, washL + trainL - 1), RS.slices(washL, washL + trainL - 1), regular);
      float endpt = 2 * washL + trainL - 1;
      if(endpt > (L - 1)) {
        endpt = L - 1;
      }
      arma::cube data_hat = OutCalc(Wout, RS.slices(washL + trainL, endpt));
      arma::vec err = arma::vec(n_init);
      for(int j = 0; j < n_init; j++) {
        arma::mat e = mean_nrmse(data_hat.slice(j), data.rows(washL + trainL, endpt));
        err(j) = e(0, 0);
      }
      errors(b, i) = arma::mean(err);
    }
  }
  float bscale;
  float iscale;
  float min_err = errors.min();
  for(int b = 0; b < 3; b++) {
    for(int i = 0; i < 3; i++) {
      if(errors(b, i) == errors.min()) {
        bscale = bscales(b);
        iscale = iscales(i);
      }
    }
  }
  arma::vec out {RNN_dim, bscale, iscale, min_err};
  return out;
}

// [[Rcpp::export]]
arma::cube CCalc(arma::cube ResStates, float aperture) {
  float RNN_dim = ResStates.n_rows;
  float n_init = ResStates.n_cols;
  float L = ResStates.n_slices;
  arma::cube C = arma::cube(RNN_dim, RNN_dim, n_init);
  for(int i = 0; i < n_init; i++) {
    arma::mat RStemp = ResStates.col(i);
    arma::mat R = RStemp * trans(RStemp) / L;
    arma::mat U, V, O;
    arma::vec D;
    arma::svd(U, D, V, R);
    arma::mat S = arma::diagmat(D) * inv(arma::diagmat(D) + pow(aperture, -2) * O.eye(RNN_dim, RNN_dim));
    C.slice(i) = U * S * trans(U);
  }
  return C;
}

// [[Rcpp::export]]
arma::cube CRNN1(arma::mat data, arma::cube W1, arma::cube C, float RSInitial, float washL, float bscale, float iscale) {
  float RNN_dim = C.n_rows;
  float n_init = C.n_slices;
  float data_dim = data.n_cols;
  float L = data.n_rows;
  arma::mat dat = scaleCpp(data);
  arma::cube RS = RSInitial * arma::cube(RNN_dim, n_init, 1, arma::fill::ones);
  arma::cube ResStates = arma::join_slices(RS, arma::cube(RNN_dim, n_init, L));
  arma::cube CStates = arma::join_slices(RS, arma::cube(RNN_dim, n_init, L));
  for(int t = 0; t < washL; t++) {
    for(int init = 0; init < n_init; init++) {
      ResStates.slice(t + 1).col(init) = tanh(bscale * W1.slice(init).col(0) + iscale * W1.slice(init).cols(1, data_dim) * trans(dat.row(t)) + W1.slice(init).cols(data_dim + 1, data_dim + RNN_dim) * ResStates.slice(t).col(init));
      CStates.slice(t + 1).col(init) = ResStates.slice(t + 1).col(init);
    }
  }
  for(int t = washL; t < L; t++) {
    for(int init = 0; init < n_init; init++) {
      ResStates.slice(t + 1).col(init) = tanh(bscale * W1.slice(init).col(0) + iscale * W1.slice(init).cols(1, data_dim) * trans(dat.row(t)) + W1.slice(init).cols(data_dim + 1, data_dim + RNN_dim) * ResStates.slice(t).col(init));
      CStates.slice(t + 1).col(init) = C.slice(init) * ResStates.slice(t + 1).col(init);
    }
  }
  arma::cube CRNN = arma::join_slices(ResStates.slices(1, L), CStates.slices(1, L));
  return CRNN;
}

// [[Rcpp::export]]
arma::mat angleCalc(arma::cube RS, arma::cube CRS) {
  float L = RS.n_slices;
  float n_init = RS.n_cols;
  arma::mat angles = arma::mat(L, n_init);
  for(int i = 0; i < n_init; i++) {
    arma::mat CRSt = CRS.col(i);
    arma::mat RSt = RS.col(i);
    angles.col(i) = diagvec(trans(CRSt) * RSt)/ sqrt(diagvec(trans(CRSt) * CRSt)) / sqrt(diagvec(trans(RSt) * RSt));
  }
  return angles;
}

// [[Rcpp::export]]
arma::vec KSCalc(arma::vec seq, int washL, int trainL) {
  float L = seq.n_elem;
  arma::vec KSseq = arma::zeros(L);
  arma::vec left_len = arma::regspace(1, L - washL - trainL - 1);
  arma::vec right_endpt = left_len - 1 + 2 * washL + 2 * trainL;
  for(int ii = 0; ii < left_len.n_elem; ii++) {
    if(right_endpt(ii) > (L - 1)) {
      right_endpt(ii) = L - 1;
    }
  }
  arma::vec right_len = right_endpt - washL - trainL - left_len + 1;
  arma::vec len = left_len + right_len;
  for(int tpt = 0; tpt < (L - washL - trainL - 1); tpt++) {
    arma::vec sortseq = arma::sort(seq.rows(washL + trainL, right_endpt(tpt)));
    arma::vec leftseq = seq.rows(washL + trainL, washL + trainL + left_len(tpt) - 1);
    arma::vec rightseq = seq.rows(washL + trainL + left_len(tpt), right_endpt(tpt));
    arma::vec left = arma::vec(len(tpt));
    arma::vec right = arma::vec(len(tpt));
    for(int t = 0; t < len(tpt); t++) {
      float left_count = 0;
      for(int lt = 0; lt < left_len(tpt); lt++) {
        if(leftseq(lt) <= sortseq(t)) {
          left_count = left_count + 1;
        }
      }
      left(t) = left_count / left_len(tpt);
      float right_count = 0;
      for(int rt = 0; rt < right_len(tpt); rt++) {
        if(rightseq(rt) <= sortseq(t)) {
          right_count = right_count + 1;
        }
      }
      right(t) = right_count / right_len(tpt);
    }
    float coef1 = left_len(tpt) * right_len(tpt) / pow(len(tpt), 2);
    arma::vec qscale = {pow((left_len(tpt) / len(tpt)), 0.5) * pow(1 - (left_len(tpt) / len(tpt)), 0.5), 0.01};
    float coef2 = arma::max(qscale);
    KSseq(washL + trainL + tpt) = coef1 * arma::max(abs(left - right)) / coef2;
  }
  return KSseq;
}

// [[Rcpp::export]]
arma::mat select(arma::vec KSseq, int washL, int trainL) {
  float L = KSseq.n_rows;
  arma::vec points = arma::zeros(L);
  arma::vec stats = arma::zeros(L);
  float start = washL + trainL;
  float end = L - washL - trainL;
  for(int t = start; t < end; t++) {
    float startpt = t - washL - trainL;
    float endpt = t + washL + trainL;
    if(startpt < 0) {
      startpt = 0;
    }
    if(endpt > (L - 1)) {
      endpt = L - 1;
    }
    if(KSseq(t) >= KSseq.rows(startpt, endpt).max()) {
      points(t) = t + 1;
      stats(t) = KSseq(t);
    }
  }
  arma::mat output = arma::join_horiz(points, stats);
  arma::mat ret = output.rows(find(output.col(0) > 0));
  return ret;
}

// [[Rcpp::export]]
arma::cube genBD(arma::mat data, int washL, int trainL, int BL, int n_boot) {
  float L = data.n_rows;
  float n_block = ceil((L - washL - trainL) / BL);
  arma::mat wrap_data = arma::join_vert(data.rows(washL + trainL, L - 1), data.rows(washL + trainL, L - 1));
  arma::mat ind = floor(arma::randu(n_block, n_boot) * (L - washL - trainL));
  arma::cube bd_temp = arma::cube(washL + trainL + n_block * BL, data.n_cols, n_boot);
  arma::cube bd = arma::cube(L, data.n_cols, n_boot);
  for(int boot = 0; boot < n_boot; boot++) {
    bd_temp.slice(boot).rows(0, washL + trainL - 1) = data.rows(0, washL + trainL - 1);
    for(int block = 0; block < n_block; block++) {
      bd_temp.slice(boot).rows(block * BL + washL + trainL, block * BL + washL + trainL + BL - 1) = wrap_data.rows(ind(block, boot), ind(block, boot) + BL - 1);
    }
    bd.slice(boot) = bd_temp.slice(boot).rows(0, L - 1);
  }
  return bd;
}

// [[Rcpp::export]]
arma::vec MBB(arma::cube bd, arma::cube W1, arma::cube C, int washL, int trainL, arma::vec params) {
  float L = bd.n_rows;
  float n_boot = bd.n_slices;
  arma::vec mbbKS = arma::vec(n_boot);
  for(int boot = 0; boot < n_boot; boot++) {
    arma::cube CRS = CRNN1(bd.slice(boot), W1, C, 0, washL, params(1), params(2));
    arma::mat angles = angleCalc(CRS.slices(0, L - 1), CRS.slices(L, 2 * L - 1));
    arma::vec seq = arma::mean(angles, 1);
    arma::vec seq_out = arma::vec(seq);
    arma::vec KSseq = KSCalc(seq_out, washL, trainL);
    mbbKS(boot) = max(KSseq);
  }
  return mbbKS;
}

// [[Rcpp::export]]
List ESM(arma::mat data, int washL, int trainL, float threshold) {
  float L = data.n_rows;
  float data_dim = data.n_cols;
  arma::vec params = RNNParamFit(data, washL, trainL, 5);
  arma::cube W1 = genW1(params(0), data_dim, 100);
  arma::cube RS = RNN1(data, W1, 0, params(1), params(2));
  arma::cube C = CCalc(RS.slices(washL, washL + trainL - 1), 100);
  arma::cube CRS = CRNN1(data, W1, C, 0, washL, params(1), params(2));
  arma::mat angles = angleCalc(CRS.slices(0, L - 1), CRS.slices(L, 2 * L - 1));
  arma::mat seq = arma::mean(angles, 1);
  arma::vec seq_out = arma::vec(seq);
  arma::vec KSseq = KSCalc(seq_out, washL, trainL);
  arma::mat sel = select(KSseq, washL, trainL);
  float n_boot = 239;
  float flags = sel.n_rows;
  arma::mat mbbKS = arma::mat(n_boot, flags);
  float quant = 1;
  float BL = ceil(pow(L, 1/3));
  arma::vec zero = arma::zeros(flags);
  arma::mat POI = arma::join_horiz(sel, zero);
  for(int flag = 0; flag < flags; flag++) {
    float bdL = sel(flag, 0) + washL + trainL;
    if(bdL > (L - 1)) {
      bdL = L - 1;
    }
    arma::cube boot_data = genBD(data.rows(0, bdL - 1), washL, trainL, BL, n_boot);
    arma::mat mbbKStemp = MBB(boot_data, W1, C, washL, trainL, params);
    float count = 1;
    for(int counter = 0; counter < n_boot; counter++) {
      if(sel(flag, 1) <= mbbKStemp(counter)) {
        count = count + 1;
      }
    }
    quant = count / (n_boot + 1);
    POI(flag, 2) = quant;
    mbbKS.col(flag) = mbbKStemp;
    if(quant <= threshold) {
      break;
    }
  }
  arma::mat retPOI = POI.rows(find(POI.col(2) > 0));
  arma::mat retmbbKS = mbbKS.cols(find(mbbKS.row(0) > 0));
  return List::create(_["POI"] = retPOI, _["MBBnull"] = retmbbKS);
}

// [[Rcpp::export]]
List SCCP(arma::mat data, int washL, int trainL, float threshold, bool FWER) {
  int L = data.n_rows;
  int Pmax = floor((L - washL - trainL) / (washL + trainL + 1));
  arma::vec ret_estimates = arma::zeros(Pmax);
  arma::vec ret_stats = arma::zeros(Pmax);
  arma::vec ret_quants = arma::zeros(Pmax);
  arma::vec cutoff;
  if(FWER) {
    cutoff = threshold / arma::regspace(Pmax, -1, 1);
  } else {
    cutoff = threshold * arma::vec(Pmax, arma::fill::ones);
  }
  int n_rejects = 0;
  List fwd0 = ESM(data, washL, trainL, cutoff(0));
  arma::mat POI = fwd0[0];
  int points = POI.n_rows;
  int estimate = POI(points - 1, 0);
  float quant = POI(points - 1, 2);
  arma::mat MBBnull = fwd0[1];
  while(estimate < (L - washL - trainL)) {
    if(quant > cutoff(n_rejects)) {
      estimate = L;
    } else {
      ret_estimates(n_rejects) = estimate;
      ret_stats(n_rejects) = POI(points - 1, 1);
      ret_quants(n_rejects) = quant;
      n_rejects = n_rejects + 1;
      if(estimate >= (L - 2 * washL - 2* trainL)) {
        estimate = L;
      } else {
        arma::mat new_data = data.rows(estimate, L - 1);
        List fwd = ESM(new_data, washL, trainL, cutoff(n_rejects));
        arma::mat tempPOI = fwd[0];
        if(tempPOI.n_rows != 0) {
          tempPOI.col(0) = tempPOI.col(0) + estimate;
          arma::mat tempMBB = fwd[1];
          POI = arma::join_vert(POI, tempPOI);
          MBBnull = arma::join_horiz(MBBnull, tempMBB);
          points = POI.n_rows;
          estimate = POI(points - 1, 0);
          quant = POI(points - 1, 2);
        } else {
          estimate = L;
        }
      }
    }
  }
  arma::mat return_mat = arma::join_horiz(arma::join_horiz(ret_estimates, ret_stats), ret_quants);
  return_mat = return_mat.rows(find(return_mat.col(0) > 0));
  return List::create(_["estCPs"] = return_mat, _["POI"] = POI, _["MBB"] = MBBnull);
}

// [[Rcpp::export]]
arma::vec Reconcile(arma::mat fwd_estimates, arma::mat bwd_estimates, int gammastar, float L) {
  arma::vec taufb = arma::sort(arma::join_vert(fwd_estimates.col(0), bwd_estimates.col(0)));
  int nfb = taufb.n_rows;
  arma::vec tau;
  if(nfb == 0) {
    tau = arma::vec{0, L};;
  } else {
    arma::mat mathcalB = arma::zeros(nfb, 2);
    mathcalB(0, 0) = taufb(0);
    mathcalB(0, 1) = taufb(0);
    int nB = 1;
    int ntau = 1;
    while(ntau < nfb) {
      if((taufb(ntau) - mathcalB(nB - 1, 1)) < gammastar) {
        mathcalB(nB - 1, 1) = taufb(ntau);
        ntau = ntau + 1;
      } else {
        mathcalB(nB, 0) = taufb(ntau);
        mathcalB(nB, 1) = taufb(ntau);
        ntau = ntau + 1;
        nB = nB + 1;
      }
    }
    tau = arma::zeros(nfb);
    int nt = 0;
    for(int nb = 0; nb < nB; nb ++) {
      if((mathcalB(nb, 1) - mathcalB(nb, 0)) < gammastar) {
        tau(nt) = floor((mathcalB(nb, 0) + mathcalB(nb, 1)) / 2);
        nt = nt + 1;
      } else {
        tau(nt) = mathcalB(nb, 0);
        tau(nt + 1) = mathcalB(nb, 1);
        nt = nt + 2;
      }
    }
  }
  return tau;
}

// [[Rcpp::export]]
List MCCP(arma::mat data, int washL, int trainL, float threshold, bool FWER) {
  int L = data.n_rows;
  List fwd = SCCP(data, washL, trainL, threshold / 2, FWER);
  arma::mat fwd_estimates = fwd[0];
  arma::mat fwd_POI = fwd[1];
  arma::mat fwd_MBB = fwd[2];
  arma::mat bwddata = arma::flipud(data);
  List bwd = SCCP(bwddata, washL, trainL, threshold / 2, FWER);
  arma::mat bwd_estimates = bwd[0];
  arma::mat bwd_POI = bwd[1];
  arma::mat bwd_MBB = bwd[2];
  bwd_estimates.col(0) = L - bwd_estimates.col(0);
  bwd_POI.col(0) = L - bwd_POI.col(0);
  arma::vec tau = Reconcile(fwd_estimates.col(0), bwd_estimates.col(0), washL + trainL + 1, L);
  tau = tau(find(tau > 0));
  tau = tau(find(tau < L));
  return List::create(_["estCPs"] = tau, _["fwd_estCPs"] = fwd_estimates, _["bwd_estCPs"] = bwd_estimates, _["fwd_POI"] = fwd_POI, _["bwd_POI"] = bwd_POI, _["fwd_MBB"] = fwd_MBB, _["bwd_MBB"] = bwd_MBB);
}