//	Logconcave Cox PH model code
//	Vector of derivatives will be very costly
	//use Active set algorithm with inner/outer loops to reduce number of times we calculate them
//	Need to have a step where we update nu = exp( x^t b)

//Allow for augmented estimator
	// If data is uncensored, expected missing is 1/n on each side
	// If right most data is uncensored, this is motivated for the same reason
	// If data is right censored, is there any motivation for this?

//h(t) = h_0(t) * nu
//f(t) = h(t) * S(t)
//f(t) = h_0(t) * nu * S_0(t)^nu
//f(t) = nu * f_0(t)/S_0(t) * S_0(t)^nu
//f(t) = nu * f_0(t) * S_0(t)^(nu-1)
//log(f(t)) = log(nu) + log(f_(t)) + (nu-1) * log(S_0(t))

//#include "Cpp2Source/baseClass.cpp"
//#include "logconcave.cpp"
#include "Cpp2Source/logconcave.cpp"

extern "C" SEXP LC_CoxPH(SEXP R_L, SEXP R_R, SEXP R_x,
	SEXP R_b, SEXP R_actInds, SEXP move_X,
	SEXP AugL, SEXP AugR, SEXP covariates){
	bool c_move_x = LOGICAL(move_X)[0] == TRUE;
	int n = LENGTH(R_L);
	int k = LENGTH(R_x);
	vector<double> x(k);
	vector<double> b(k);
	for(int i = 0; i < k; i++){
		x[i] = REAL(R_x)[i];
		b[i] = REAL(R_b)[i];
	}
	vector<int> cL(n);
	vector<int> cR(n);
	vector<double> cRepVec(n);
	bool leftDone,rightDone;
	double leftval, rightval;
	for(int i = 0; i < n; i++){
		leftval = REAL(R_L)[i];
		rightval = REAL(R_R)[i];
		leftDone = false;
		rightDone = false;
		for(int j = 0; j < k; j++){
			if(leftDone == false){
				if(leftval == x[j]){
					cL[i] = j;
					leftDone = true;
				}	
			}
			if(rightDone == false){
				if(rightval == x[j]){
					cR[i] = j;
					rightDone = true;
				}
			}
		}
		if(rightDone == false | leftDone == false){
			Rprintf("problem: data does not match x! Quiting\n");
			return(R_NilValue);
		}
	}			
	SEXP dims = getAttrib(covariates, R_DimSymbol);
	PROTECT(dims);
	int num_cov = INTEGER(dims)[1];
	int cov_n = INTEGER(dims)[0];
	UNPROTECT(1);
	
	if(cov_n != n){
		Rprintf("Number of covariates does not match number of observations\n");
		return(R_NilValue);
	}
	

	QuadProgPP::Matrix<double> Covars(cov_n,num_cov);
	for(int i = 0; i < cov_n; i++){
		for(int j = 0; j < num_cov; j++)
			Covars[i][j] = REAL(covariates)[i + j*cov_n];
	}

	
	int ak = LENGTH(R_actInds);
	vector<int> actInds(ak);	
	for(int i = 0; i < ak; i++)	
		actInds[i] = INTEGER(R_actInds)[i] - 1;
	
	
	double augl = REAL(AugL)[0];
	double augr = REAL(AugR)[0];

	
	LogConCenPH optObj(0, x, b, cL, cR, actInds, c_move_x, augl, augr, Covars, num_cov);



	optObj.updateNu();

	optObj.updateRegress();

	optObj.VEMstep();	
	optObj.ICMstep();
	
//	double old_llk = -INFINITY;
		
//	double new_llk = optObj.llk();
	
	int inner_it = 0;
	int outer_it = 0;
	int max_inner_its = 2*sqrt(optObj.x.size()) + 5;
	int max_outer_its = 100;
	
		
	double outer_tol = pow(10.0, -10);
	double inner_tol = pow(10.0, -12);
	double outer_error = outer_tol + 1;
	double reg_error = outer_tol + 1;
	double inner_error;
	double inner_llk, outer_llk;
	int loopcount = 0;
	bool start_move_x = false;
	while( (outer_it < max_outer_its) & ( (outer_error > outer_tol) | (loopcount < 2))){
		outer_it++;
		outer_llk = optObj.llk();
		optObj.VEMstep();
		inner_error = inner_tol + 1;
		reg_error = outer_tol + 1;
		inner_it = 0;
		while(inner_it < max_inner_its & inner_error > inner_tol){
			inner_it++;
			if(start_move_x & c_move_x)
				optObj.updateXs();
			if(reg_error > outer_tol){
				reg_error = optObj.updateRegress();
				}
			inner_llk = optObj.llk();
			optObj.ICMstep();
			inner_error = optObj.llk() - inner_llk;
			if(inner_error != inner_error)
				break;
			optObj.checkEnds();			
		}
		optObj.recenterBeta();
		outer_error = optObj.llk() - outer_llk;	
		if(outer_error != outer_error){
			Rprintf("Warning: undefined error. Algorithm terminated at invalid estimate\n");
			break;
		}
		if(outer_error < 0.0001)
			start_move_x = true;
		if(outer_error < outer_tol){
			loopcount++;
			}
		else 
			loopcount = 0;
	}
	
	optObj.makePropDist();

	ak = optObj.getAK();
	SEXP output = PROTECT(allocVector(VECSXP, 5) );
	SEXP regressP = PROTECT(allocVector(REALSXP, num_cov) );
	for(int i = 0; i < num_cov; i++)
		REAL(regressP)[i] = optObj.cov_b[i];
	SEXP x_out = PROTECT(allocVector(REALSXP, ak) );
	SEXP b_out = PROTECT(allocVector(REALSXP, ak) );
	SEXP llk_out = PROTECT(allocVector(REALSXP, 1) );
	for(int i = 0; i < ak; i++){
		REAL(x_out)[i] = optObj.x[optObj.actIndex[i]];
		REAL(b_out)[i] = optObj.b[optObj.actIndex[i]];
	}
	SEXP Hess_out = PROTECT(allocMatrix(REALSXP, num_cov, num_cov));
	for(int i = 0; i < num_cov; i++){
		for(int j = 0; j < num_cov; j++)
			REAL(Hess_out)[i + num_cov*j] = optObj.Hess_b[i][j];
	}
	
	REAL(llk_out)[0] = optObj.llk();
	SET_VECTOR_ELT(output, 0, x_out);
	SET_VECTOR_ELT(output, 1, b_out);
	SET_VECTOR_ELT(output, 2, llk_out);
	SET_VECTOR_ELT(output, 3, regressP);
	SET_VECTOR_ELT(output, 4, Hess_out);
	UNPROTECT(6);
	return(output);
}

LogConCenPH::LogConCenPH(int MoveX, 
				   	vector<double> X, 
				   	vector<double> B, 
				   	vector<int> L_ind, 
				   	vector<int> R_ind,
				   	vector<int> actInds,
				   	bool allow_x_move,
				   	double augl,
					double augr,
					QuadProgPP::Matrix<double> Covars,
				   	int num_cov
				   	){
				   	x = X;
				   	b = B;
				   	Lindex = L_ind;
				   	Rindex = R_ind;
				   	actIndex = actInds;
				   	allowMoveX = MoveX;
				   	h = 0.0001;
				   	slpTol = pow(10.0, -14.5);
				   	cur_Err = INFINITY;
				   	endTol = -10;
				   	int k = X.size();
				   	dx.resize(k);
				   	for(int i = 0; i < (k-1); i++)
				   		dx[i] = x[i+1] - x[i];
					allActDervs = vector<double>(k);
					s.resize(k);
					p.resize(k);
					int n = Lindex.size();
					augLeft = augl;
					augRight = augr;
					covars = Covars;
					nu.resize(n);
					for(int i = 0; i < n ; i ++)
						nu[i] = 1.0;
					cov_b.resize(num_cov);
					for(int i = 0; i < num_cov; i++)
						cov_b[i] = 0;
					cov_b_d.resize(num_cov);
					Hess_b.resize(num_cov, num_cov);
					test_Hess_b.resize(num_cov, num_cov);
					b_h.resize(num_cov);
					b_l.resize(num_cov);
				   }
				   
	
double LogConCenPH::llk(){
	int k = x.size();
	double db;
	s[0] = 0;
//	double der;
	if(x[0] > -INFINITY){
		if(b[1] == -INFINITY | b[0] == -INFINITY){
				s[1] = s[0];
			} else {
			db = b[1] - b[0];
			if(db < 0.00001 && db > -0.00001)	
				s[1] = exp(b[0]) * dx[0] *(1 + db/2);
			else			
				s[1] = dx[0]/db * (exp(b[1]) - exp(b[0]) ); 
		}
	}
	if(x[0] == -INFINITY){
		if(b[1] == -INFINITY){
			s[1] = 0;
		} else {
			db = (b[2] - b[1]) / dx[1];
			if(db <= 0){
				return(-INFINITY);
			}
			s[1] = exp(b[1])/db;
			}
		}
	for(int i = 1; i < k-2; i++){
		if(b[i+1] == -INFINITY | b[i] == -INFINITY){
			s[i+1] = s[i];
			continue;
		} 
		db = b[i+1] - b[i];
		if(db < 0.00001 && db > -0.00001)	
			s[i+1] = s[i] + exp(b[i]) * dx[i] *(1 + db/2);
		else			
			s[i+1] = s[i] + dx[i]/db * (exp(b[i+1]) - exp(b[i]) ); 				
		}
	if(x[k-1] < INFINITY){
		if(b[k-1] == -INFINITY | b[k-2] == -INFINITY){
			s[k-1] = s[k-2];
			} 
		else{
			db = b[k-1] - b[k-2];
			if(db < 0.00001 && db > -0.00001)	
				s[k-1] = s[k-2] + exp(b[k-1]) * dx[k-2] *(1 + db/2);
			else			
				s[k-1] = s[k-2] + dx[k-2]/db * (exp(b[k-1]) - exp(b[k-2]) ); 
			}
		}
	if(x[k-1] == INFINITY){
		if(b[k-2] == -INFINITY){
			s[k-1] = s[k-2];
			}	 
		else {
			db = (b[k-2] - b[k-3]) / dx[k-3];//(x[k-2] - x[k-3]);
			if(db >= 0){
				return(-INFINITY);
			}
			s[k-1] = s[k-2] - exp(b[k-2])/db;
		}					
	}
	scaleValue = (1 - augLeft - augRight) / s[k-1];
	double logScale = log(scaleValue);
			
	//Should be scaleValue * b[i] for uncensored, I think??

	for(int i = 0; i < k; i++)
		s[i] = s[i] * scaleValue + augLeft;		
	int n_row = Lindex.size();
//	double tot_rvec = 0;
	double log_sum = 0;
	double p_ob = 0;
	int hi_ind;
	int lo_ind;
	for(int i = 0; i < n_row; i++){
		hi_ind = Rindex[i];
		lo_ind = Lindex[i];
		if(hi_ind == lo_ind){
			if(abs(nu[i]-1) >  0.00000001)
				log_sum += log(nu[i]) + (b[lo_ind] + logScale) + (nu[i]-1) * log(1-s[hi_ind]);	
			else
				log_sum += (b[lo_ind] + logScale);
			continue;
		}
		p_ob  =  pow(1-s[lo_ind], nu[i]) - pow(1-s[hi_ind], nu[i]);
					
		if(p_ob == 0) {
			 return(-INFINITY);
		}
		log_sum = log_sum + log(p_ob);
	}
	if(log_sum == INFINITY | log_sum == -INFINITY){
		return(-INFINITY);
	}
	if(log_sum != log_sum){
		return(-INFINITY);
	}
	return (log_sum);
}

double LogConCenPH::nullk(){
//	int k = x.size();
	updateNu();
//	double scaleValue = (1 - augLeft - augRight) / s[k-1];
	double logScale = log(scaleValue);
//	for(int i = 0; i < k; i++)
//		s[i] = s[i] * scaleValue + augLeft;		
	int n_row = Lindex.size();
//	double tot_rvec = 0;
	double log_sum = 0;
	double p_ob = 0;
	int hi_ind;
	int lo_ind;
	for(int i = 0; i < n_row; i++){
		hi_ind = Rindex[i];
		lo_ind = Lindex[i];
		if(hi_ind == lo_ind){
			if( abs(nu[i] - 1) >  0.000000001 )
				log_sum += log(nu[i]) + (b[lo_ind] + logScale) + (nu[i]-1) * log(1-s[hi_ind]);	
			else
				log_sum += (b[lo_ind] + logScale);
			continue;
		}
		p_ob  =  pow(1-s[lo_ind], nu[i]) - pow(1-s[hi_ind], nu[i]);
		if(!(p_ob > 0) ) {
			 return(-INFINITY);
		}
		log_sum += log(p_ob);
	}
	if(log_sum == INFINITY | log_sum == -INFINITY){
		return(-INFINITY);
	}
	if(log_sum != log_sum){
		return(-INFINITY);
	}
	return (log_sum);
}


void LogConCenPH::calcDervVec(){
/*	int k = x.size();
	for(int i = 0; i < k; i++){
		allActDervs[i] = numericDervs(i)[0];		
	}	*/
	fastNumActDers();
}

void LogConCenPH::updateNu(){
	int k = cov_b.size();
	int n = nu.size();
	for(int i = 0; i < n ; i++){
		nu[i] = 0;
		for(int j = 0; j < k; j++)
			nu[i] = nu[i] + covars[i][j] * cov_b[j];
		nu[i] = exp(nu[i] );
	}	
}

void LogConCenPH::calcNuDerv(){
	double llk_0 = llk();
	llk_0 = nullk();
	double llk_ll, llk_hh;
	double pd;
	int k = cov_b.size();
	for(int i = 0; i < k; i++){
		cov_b[i] = cov_b[i] + h;
		b_h[i] = nullk();
		cov_b[i] = cov_b[i] - 2*h;
		b_l[i] = nullk();
		cov_b[i] = cov_b[i] + h;
		cov_b_d[i] = (b_h[i] - b_l[i])/(2*h);
		Hess_b[i][i] = (b_h[i] + b_l[i] - 2 * llk_0)/pow(h,2.0);
	}
	
	
	for(int i = 0; i < k; i++){
		for(int j = i+1; j < k; j++){
			cov_b[i] = cov_b[i] + h;
			cov_b[j] = cov_b[j] + h;
			llk_hh = nullk();
			cov_b[i] = cov_b[i] - 2*h;
			cov_b[j] = cov_b[j] - 2*h;
			llk_ll = nullk();
			cov_b[i] = cov_b[i] + h;
			cov_b[j] = cov_b[j] + h;
			pd = (llk_hh + llk_ll + 2 * llk_0 - b_l[i] - b_h[i] - b_l[j] - b_h[j] ) / (4 * pow(h , 2.0) );
			Hess_b[i][j] = pd;
			Hess_b[j][i] = pd;
		}
	}	
}

double LogConCenPH::updateRegress(){
	double llk_begin = llk();
	calcNuDerv();
//	double deter;
//	bool can_inv;

	QuadProgPP::Vector<double> propStep(cov_b.size());
	cov_b_d = -cov_b_d;
	Hess_b = -Hess_b;	

	for(int i = 0; i < cov_b.size(); i++){
		for(int j = 0; j < cov_b.size(); j++)
			test_Hess_b[i][j] = Hess_b[i][j];
	}
//	QuadProgPP::Matrix<double> test_Hess_b = Hess_b; //I don't think this is actually copying!!
	bool decompOK = QuadProgPP::cholesky_decomposition_safe(test_Hess_b);
	
	int tries = 0;
	double diagAdd = 1;
	
	while(decompOK == false & tries < 10){
		test_Hess_b = Hess_b;
		for(int i = 0; i< cov_b_d.size(); i++)
			test_Hess_b[i][i] = test_Hess_b[i][i] + diagAdd;
		diagAdd = diagAdd*2;
		tries++;
		decompOK = QuadProgPP::cholesky_decomposition_safe(test_Hess_b);
	}
	if(decompOK == false){
		Rprintf("warning: Hessian for covariates not positive definite\n");
		return(1.0);
	}
	cholesky_solve(test_Hess_b, propStep, cov_b_d);
	if(propStep[0] != propStep[0]){
		Rprintf("warning: propstep undefined\n");
		return(1.0);
	}
	for (int i = 0; i < cov_b_d.size(); i++)
	    propStep[i] = -propStep[i];

	for(int i = 0; i < cov_b_d.size(); i++)
		cov_b[i] = cov_b[i] + propStep[i];
	double llk_new = nullk();
	if(llk_new < llk_begin){
		propStep = -propStep;
		int it = 0;
		while(it < 5 & llk_new < llk_begin){
			it++;
			for(int i = 0; i < propStep.size(); i++)
				propStep[i] = propStep[i]/2;
			for(int i = 0; i < cov_b_d.size(); i++)
				cov_b[i] = cov_b[i] + propStep[i];
			llk_new = nullk();
		}
		if(llk_new < llk_begin){
			for(int i = 0; i < cov_b_d.size(); i++)
				cov_b[i] = cov_b[i] + propStep[i];
		llk_new = nullk();	
		}
	}
	return(llk_new - llk_begin);
}


void LogConCenPH::update_p(int index1, int index2){
	int k = x.size();
	p[0] = 0;
	double slp = (b[index2] - b[index1])/(x[index2] - x[index1]);
	if(index1 == 0){
		if(x[0] > -INFINITY){
		if(b[1] == -INFINITY | b[0] == -INFINITY){
				p[1] = 0;
			} else {
			double db = b[1] - b[0];
			if(db < 0.00001 && db > -0.00001)	
				p[1] = exp(b[0]) * dx[0] *(1 + db/2);
			else			
				p[1] = 1/slp * (exp(b[1]) - exp(b[0]) ); 
		}
	}
	if(x[0] == -INFINITY){
		if(b[1] == -INFINITY){
			p[1] = 0;
		} else {
			double db = (b[2] - b[1]) / dx[1];
			p[1] = exp(b[1])/db;
			}
		}
	index1 = 1;
	}
	if(slp == 0){
		for(int i = index1; i < index2; i++)
			p[i+1] = exp(b[i]) * dx[i]; 
	}
	else{
		for(int i = index1; i < index2; i++)
			p[i+1] = 1/slp * (exp(b[i+1]) - exp(b[i]) ); 				
	}
	
	if(index2 == k-1){
		if(x[k-1] < INFINITY){
			if(b[k-1] == -INFINITY | b[k-2] == -INFINITY){
				p[k-1] = 0;
				} 
			else{
				double db = b[k-1] - b[k-2];
				if(db < 0.00001 && db > -0.00001)	
					p[k-1] = exp(b[k-1]) * dx[k-2] *(1 + db/2);
				else			
					p[k-1] = dx[k-2]/db * (exp(b[k-1]) - exp(b[k-2]) ); 
				}
			}
		if(x[k-1] == INFINITY){
			if(b[k-2] == -INFINITY)
				p[k-1] = 0;
			else {
				double db = (b[k-2] - b[k-3]) / dx[k-3];//(x[k-2] - x[k-3]);
				p[k-1] = -exp(b[k-2])/db;
			}					
		}
	}
}

void LogConCenPH::p2s(){
	int k = x.size();
	s[0] = 0;
	for(int i = 1; i < k; i++)
		s[i] = p[i]+s[i-1];
}

double LogConCenPH::fastBasellk(){
	int k = x.size();
	p2s();
	scaleValue = (1 - augLeft - augRight) / s[k-1];
	double logScale = log(scaleValue);
			
	for(int i = 0; i < k; i++)
		s[i] = s[i] * scaleValue + augLeft;		
		
	int n_row = Lindex.size();
//	double tot_rvec = 0;
	double log_sum = 0;
	double p_ob = 0;
	int hi_ind;
	int lo_ind;
	for(int i = 0; i < n_row; i++){
		hi_ind = Rindex[i];
		lo_ind = Lindex[i];
		if(hi_ind == lo_ind){
			if(abs(nu[i] - 1) > 0.00000001)
				log_sum += log(nu[i]) + (b[lo_ind] + logScale) + (nu[i]-1) * log(1-s[hi_ind]);	
			else
				log_sum += (b[lo_ind] + logScale);
			continue;
		}
		p_ob  =  pow(1-s[lo_ind], nu[i]) - pow(1-s[hi_ind], nu[i]);
					
		if(p_ob == 0) {
			 return(-INFINITY);
		}
		log_sum += log(p_ob);
	}
	if(log_sum == INFINITY | log_sum == -INFINITY){
		return(-INFINITY);
	}
	if(log_sum != log_sum){
		return(-INFINITY);
	}
	return (log_sum);
}

void LogConCenPH::fastNumActDers(){
	double llk_l, llk_h;// = llk();
	double h_fast = h/100;	
	int k = x.size();
	for(int i = 0; i < k; i++){
		p[i] = 0;
		allActDervs[i] = 0;
	}
	for(int act_i = 0; act_i < getAK() - 1; act_i++)
		update_p(actIndex[act_i], actIndex[act_i+1]);
	for(int act_i = 0; act_i < getAK() - 1; act_i++){
		for(int i = actIndex[act_i]+1; i < actIndex[act_i+1]-1; i++){
			move_act_b(i, h_fast);
			update_p(actIndex[act_i], i);
			update_p(i, actIndex[act_i+1]);
			llk_h = fastBasellk();
			move_act_b(i, -2*h_fast);
			update_p(actIndex[act_i], i);
			update_p(i, actIndex[act_i+1]);
			llk_l = fastBasellk();
			move_act_b(i, h_fast);
			update_p(actIndex[act_i], i);
			update_p(i, actIndex[act_i+1]);
			allActDervs[i] = (llk_h-llk_l)/(h_fast);
		}
	}
}

