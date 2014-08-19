#include "baseClass.cpp"

#include <R.h>
#include <Rinternals.h>

extern "C" SEXP uniVarLCCens(SEXP R_L, SEXP R_R, SEXP R_x,
	SEXP R_b, SEXP R_repVec, SEXP R_actInds, SEXP move_X){
	bool c_move_x = LOGICAL(move_X)[0] == TRUE;
	int n = LENGTH(R_L);
	vector<int> cL(n);
	vector<int> cR(n);
	vector<double> cRepVec(n);
	for(int i = 0; i < n; i++){
		cL[i] = INTEGER(R_L)[i] - 1;
		cR[i] = INTEGER(R_R)[i] - 1;
		cRepVec[i] = REAL(R_repVec)[i];
	}		
	
	int k = LENGTH(R_x);
	vector<double> x(k);
	vector<double> b(k);
	for(int i = 0; i < k; i++){
		x[i] = REAL(R_x)[i];
		b[i] = REAL(R_b)[i];
	}
	

	int ak = LENGTH(R_actInds);
	vector<int> actInds(ak);	
	for(int i = 0; i < ak; i++)	
		actInds[i] = INTEGER(R_actInds)[i] - 1;
	
	LogConCen optObj(0, x, b, cL, cR, actInds, cRepVec, c_move_x);


	int max_it = 500;
	int it = 0;
	double tol = pow(10.0, -10.0);
	double oldLike = -INFINITY;
	double Err = tol + 1;
	double newLike = 0;
	bool startMoveX = false;
	
	while(it < max_it & Err > tol){
		it ++;		
		optObj.VEMstep();
		
		optObj.ICMstep();		
		
		optObj.recenterBeta();
			
		optObj.checkEnds();

		if(startMoveX & c_move_x){
			for(int i = 1; i < optObj.getAK()- 1; i++){
				if(i < optObj.getAK() - 1 & optObj.actIndex[i] % 2 == 1)
					optObj.updateX(i);
			}
		}
			
		newLike = optObj.llk();
		Err = newLike - oldLike;
		oldLike = newLike;
		if((Err < tol) & (!startMoveX)){
			startMoveX = true;
			Err = tol + 1;
			}
	}
	optObj.makePropDist();

	ak = optObj.getAK();
	SEXP output = PROTECT(allocVector(VECSXP, 3) );
	SEXP x_out = PROTECT(allocVector(REALSXP, ak) );
	SEXP b_out = PROTECT(allocVector(REALSXP, ak) );
	SEXP llk_out = PROTECT(allocVector(REALSXP, 1) );
	for(int i = 0; i < ak; i++){
		REAL(x_out)[i] = optObj.x[optObj.actIndex[i]];
		REAL(b_out)[i] = optObj.b[optObj.actIndex[i]];
	}
	REAL(llk_out)[0] = optObj.llk();
	SET_VECTOR_ELT(output, 0, x_out);
	SET_VECTOR_ELT(output, 1, b_out);
	SET_VECTOR_ELT(output, 2, llk_out);
	UNPROTECT(4);
	return(output);
}

LogConCen::LogConCen(int MoveX, 
				   	vector<double> X, 
				   	vector<double> B, 
				   	vector<int> L_ind, 
				   	vector<int> R_ind,
				   	vector<int> actInds,
				   	vector<double> repVec,
				   	bool allow_x_move
				   	){
				   	x = X;
				   	b = B;
				   	Lindex = L_ind;
				   	Rindex = R_ind;
				   	actIndex = actInds;
				   	allowMoveX = MoveX;
				   	rep_vec = repVec;
				   	h = 0.0001;
				   	slpTol = pow(10.0, -14.5);
				   	cur_Err = INFINITY;
				   	endTol = -15;
				   	int k = X.size();
				   	dx.resize(k);
				   	for(int i = 0; i < (k-1); i++)
				   		dx[i] = x[i+1] - x[i];
				  	d_bl = vector<double>(k);
				  	d_bu = vector<double>(k);
				   	allBaseDervs = vector<double>(k);
					allActDervs = vector<double>(k);
					s.resize(k);
				   }

double LogConCen::llk() {
	int k = x.size();
//	double dx;
	double db;
//	vector<double> s;
//	s.resize(k);
	s[0] = 0;
//	double der;
	if(x[0] > -INFINITY){
		if(b[1] == -INFINITY | b[0] == -INFINITY){
				s[1] = s[0];
			} else {
//			dx = x[1] - x[0];
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
		//dx = x[i+1] - x[i];
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
		//	dx = x[k-1] - x[k-2];
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
		
	int n_row = Lindex.size();
	double tot_rvec = 0;
	double log_sum = 0;
	double p_ob = 0;
	int hi_ind;
	int lo_ind;
	for(int i = 0; i < n_row; i++){
		hi_ind = Rindex[i];
		lo_ind = Lindex[i];
		if(hi_ind == lo_ind){
			log_sum = log_sum + b[lo_ind] * rep_vec[i];
			tot_rvec = tot_rvec + rep_vec[i];
			continue;
		}
		p_ob  = s[hi_ind] - s[lo_ind];
			
		if(p_ob == 0) {
			if(rep_vec[i] == 0)
			 	p_ob = 1;
			if(rep_vec[i] != 0){
			 	return(-INFINITY);
			 }
		}
		log_sum = log_sum + log(p_ob) * rep_vec[i];
		tot_rvec = tot_rvec + rep_vec[i];
	}
	if(log_sum == INFINITY | log_sum == -INFINITY){
		return(-INFINITY);
	}
	double output = log_sum - tot_rvec * log(s[k-1]) ;
	if(output != output){
		Rprintf("Warning: likelihood function undefined!\n");
		return(-INFINITY);
	}
	return (output);
}


void LCBase::checkEnds(){
	double cur_llk = llk();
	double new_llk;
	double old_b;
	bool exit_loop = false;
	double limit = std::pow(10.0, -10.0);
	int it = 0;
	while(exit_loop == false & it < 5){
		it++;
		if(b[actIndex[0] ] < endTol){
			old_b = b[actIndex[0] ];
			b[actIndex[0]] = -INFINITY;
			new_llk = llk();
			if(new_llk > (cur_llk - limit) ){
				int moveInd = actIndex[0] ;
				addActive(moveInd+1);
				removeActive(moveInd);
	//			moveInd++;
			} else {
			b[actIndex[0]] = old_b;		
			exit_loop = true;
			}
		}
	}
	
	it = 0;
	exit_loop = false;
	int endIndex;
	while(exit_loop == false & it < 5){
		it++;
		endIndex = actIndex[getAK() - 1];
		if(b[endIndex] < endTol){
			old_b = b[endIndex];
			b[endIndex] = -INFINITY;
			new_llk = llk();
			if(new_llk > (cur_llk - limit) ){
				addActive(endIndex-1);
				removeActive(endIndex);
			} else {
			b[endIndex] = old_b;
			exit_loop = true;
			}
		}
	}
}

vector<double> LCBase::getLimits(int index){
	double minD = -INFINITY;
	double max_l = INFINITY;
	double max_r = INFINITY;
	double maxD = INFINITY;
	double l_slope = INFINITY;
	double r_slope = -INFINITY;
	double m_slope;
	int ak = getAK();
	int act_l = -1;
	int act_r = actIndex[ak-1] + 1;
	for(int i = 0; i < ak; i++){
		if(actIndex[i] > index){
			act_r = i;
			break;
		}
	}
	for(int i = ak-1; i >= 0; i--){
		if(actIndex[i] < index){
			act_l = i;
			break;
		}
	}
	
	 
	if(ak > 2){
		if(act_l > 0){
			l_slope = (b[actIndex[act_l]] - b[actIndex[act_l-1] ] )/(x[actIndex[act_l]] - x[actIndex[act_l-1]]);
			if(l_slope < INFINITY)
				max_l = b[actIndex[act_l]] + l_slope * (x[index] - x[actIndex[act_l]]) - b[index];
		}
		if(act_r < ak - 1){
			r_slope = (b[actIndex[act_r+1]] - b[actIndex[act_r] ] )/(x[actIndex[act_r+1]] - x[actIndex[act_r]]);
			if(r_slope > -INFINITY)
				max_r = b[actIndex[act_r]] + r_slope * (x[index] - x[actIndex[act_r]]) - b[index];
		}
	maxD = min(max_l, max_r);
	if(act_l >= 0 & act_r <= actIndex[ak-1]){
		m_slope = (b[actIndex[act_r]] - b[actIndex[act_l]]) / (x[actIndex[act_r]] - x[actIndex[act_l]]);
		minD = b[actIndex[act_l]] + m_slope * (x[index] - x[actIndex[act_l]]) - b[index];
		}
	}
	vector<double> output(2);
	output[0] = minD;
	output[1] = maxD;
	return(output);
}

void LogConCen::calcBaseDervs(){
	int k = x.size();
	double der = 0, db = 0;//, dx = 0;
	d_bl[k-1] = 0;
	for(int i = 0; i < k-1; i++)
		{
		d_bl[i] = 0;
		if(b[i] == -INFINITY | b[i+1] == -INFINITY)
			{
			continue;
			}
	//	dx = x[i+1] - x[i];
		db = b[i+1] - b[i];
		if(db < 0.00001 && db > -0.00001)
			{
			der = exp(b[i])/3 * dx[i];
			d_bl[i] = exp(b[i])/2 * dx[i] - der * db;
			}
		else
			{
			d_bl[i] = -(exp(b[i]) * db - (exp(b[i+1]) - exp(b[i]))) * dx[i] / pow(db, 2);
			}
		}
	if(x[0] == -INFINITY)
		d_bl[0] = 0;	
/*	if(x[0] == -INFINITY & x[1] > -INFINITY)
		{
		dx = x[1] - x[2];
		db = betas[1] - x[2];
		d_bl[1] = (exp(betas[1]) * dx * (1 - db)  )  / (db * db);
		} */
		
	if(x[k-1] == INFINITY & b[k-2] > -INFINITY)
		{
		//dx = x[k-2] - x[k-3];
		db = b[k-2] - b[k-3];
		d_bl[k-2] = (exp(b[k-2]) * dx[k-3] * (1 - db)  )  / (db * db);
		d_bl[k-3] = d_bl[k-3] - exp(b[k-2]) * dx[k-3] / (db * db);
		}
//	double d_bu[k];	
	d_bu[0] = 0;
	for(int i = 1; i < k; i++)
		{
		d_bu[i] = 0;
		if(b[i-1] == -INFINITY | b[i] == -INFINITY)
			{
			d_bu[i] = 0;
			continue;
			}
	//	dx = x[i] - x[i-1];
		db = b[i] - b[i-1];
		if(db < 0.00001 && db > -0.00001)
			{
			der = exp(b[i-1])/3 * dx[i-1];
			d_bu[i] = exp(b[i-1])/2 * dx[i-1] + der * db;
			}
		else
			{
			d_bu[i] = (exp(b[i]) * db - (exp(b[i]) - exp(b[i-1]))) * dx[i-1] / pow(db, 2);
			}
		}
		if(x[0] == -INFINITY & b[1] > -INFINITY)
			{
//			dx = x[2] - x[1];
			db = b[2] - b[1];
			d_bu[1] = (exp(b[1]) * dx[1] * (db + 1) ) / (db * db);
			d_bl[2] = d_bl[2] - exp(b[1]) * dx[1] / (db * db);
			}
	 
//		double s[k]; 
		s[0] = 0;
		if(x[0] > -INFINITY)
			{
			if(b[1] == -INFINITY | b[0] == -INFINITY)
				{
				s[1] = s[0];
				} else
				{ 
			//	dx = x[1] - x[0];
				db = b[1] - b[0];
				if(db < 0.00001 && db > -0.00001)	
					s[1] = exp(b[0]) * dx[0] *(1 + db/2);
				else			
					s[1] = dx[0]/db * (exp(b[1]) - exp(b[0]) ); 
				}
			} 
		if(x[0] == -INFINITY)
			{
			if(b[1] == -INFINITY)
				{
				s[1] = 0;
				} else
				{
				db = (b[2] - b[1]) / (x[2] - x[1]);
				if(db <= 0)
					s[1] = INFINITY;
				s[1] = exp(b[1])/db;
				}
			}
		for(int i = 1; i < k-2; i++)
			{
			if(b[i+1] == -INFINITY | b[i] == -INFINITY)
				{
				s[i+1] = s[i];
				continue;
				} 
		//	dx = x[i+1] - x[i];
			db = b[i+1] - b[i];
			if(db < 0.00001 && db > -0.00001)	
				s[i+1] = s[i] + exp(b[i]) * dx[i] *(1 + db/2);
			else			
				s[i+1] = s[i] + dx[i]/db * (exp(b[i+1]) - exp(b[i]) ); 
			}
		if(x[k-1] < INFINITY)
			{
			if(b[k-1] == -INFINITY | b[k-2] == -INFINITY)
				{
				s[k-1] = s[k-2];
				} else
				{
			//	dx = x[k-1] - x[k-2];
				db = b[k-1] - b[k-2];
				if(db < 0.00001 && db > -0.00001)	
					s[k-1] = s[k-2] + exp(b[k-1]) * dx[k-2] *(1 + db/2);
				else			
					s[k-1] = s[k-2] + dx[k-2]/db * (exp(b[k-1]) - exp(b[k-2]) ); 
				}
			}
		if(x[k-1] == INFINITY)
			{
			if(b[k-2] == -INFINITY)
				{
				s[k-1] = s[k-2];
				} else
				{
				db = (b[k-2] - b[k-3]) / dx[k-3];//(x[k-2] - x[k-3]);
				if(db >= 0)
					s[k-1] = INFINITY;
				s[k-1] = s[k-2] - exp(b[k-2])/db;
				}
			}		
		int n_row = Lindex.size();
		double tot_rvec = 0;
		for(int i = 0; i < k; i++)
			allBaseDervs[i] = 0;
		double p_ob = 0;
		int hi_ind;
		int lo_ind;
//		double lOnly;
//		double both;
//		double rOnly;
		for(int i = 0; i < n_row; i++)
			{
			if(rep_vec[i] < 0.0000000001 )
				continue;
			hi_ind = Rindex[i];
			lo_ind = Lindex[i];
			if(hi_ind == lo_ind)
				{
				allBaseDervs[lo_ind] = allBaseDervs[lo_ind] + rep_vec[i];
				tot_rvec = tot_rvec + rep_vec[i];
				continue;
				}
			p_ob  = (s[hi_ind] - s[lo_ind])/rep_vec[i];
			for(int j = lo_ind; j < hi_ind; j++)
				allBaseDervs[j] = allBaseDervs[j] + d_bl[j]/p_ob;
			for(int j = lo_ind + 1; j < hi_ind + 1; j++)
				allBaseDervs[j] = allBaseDervs[j] + d_bu[j]/p_ob;
			tot_rvec = tot_rvec + rep_vec[i];
			}
		for(int i = 0; i < k; i ++ )
			allBaseDervs[i] = allBaseDervs[i] - tot_rvec * (d_bl[i] + d_bu[i])/s[k-1];	
		if(x[0] == -INFINITY)
			allBaseDervs[0] = 0;
		if(x[k-1] == -INFINITY)
			allBaseDervs[k-1] = 0;
}	

void LCBase::checkAllActive(){
	int ak = getAK();
	vector<double> lims; 
	for(int i = 0; i < ak; i++){
		lims = getLimits(actIndex[i]);
		if(lims[0] > -slpTol){
			removeActive(actIndex[i]);
			}
	}
}

int LCBase::findMaxError(){
//	int k = x.size();
	double max_Err = 0;
	int max_Index = 0;
	int curr_act = 0;
	double mark_Err = 0;
	int begin = actIndex[0];
	int end = actIndex[getAK()-1];
	for(int i = begin; i < end; i++){
		if(i > actIndex[curr_act])
			curr_act++;
		if(i == actIndex[curr_act])
			mark_Err = abs(allActDervs[i]);
		else
			mark_Err = max(allActDervs[i], 0.0);
		if(mark_Err > max_Err){
			max_Err = mark_Err;
			max_Index = i;
		}
	}	

	cur_Err = max_Err;
	return(max_Index);
}

vector<int> LCBase::findMaxIntError(){
//	int k = x.size();
	int numPoints = getAK() - 1;
	vector<int> points(numPoints);
	double max_Err;
	int max_Index;
	int curr_act;
	double mark_Err;
	int begin;
	int end;
	for(int cur_ai = 0; cur_ai < numPoints; cur_ai++){
		max_Err = 0;
		max_Index = 0;
		curr_act = 0;
		mark_Err = 0;
		begin = actIndex[cur_ai];
		end = actIndex[cur_ai+1];

		for(int i = begin; i < end; i++){
			if(i > actIndex[curr_act])
				curr_act++;
			if(i == actIndex[curr_act])
				mark_Err = abs(allActDervs[i]);
			else
				mark_Err = max(allActDervs[i], 0.0);
			if(mark_Err > max_Err){
				max_Err = mark_Err;
				max_Index = i;
			}
		}	
	cur_Err = max_Err;
	points[cur_ai] = max_Index;
	}
	return(points);
}


void LCBase::qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec){
	int ak = getAK();
	vector<double> dx(ak - 1);
	vector<double> db(ak - 1);
	double setMatValue;
	for(int i = 0; i < ak-1; i++){
		dx[i] = x[actIndex[i+1] ] - x[actIndex[i] ];
		db[i] = b[actIndex[i+1] ] - b[actIndex[i] ];
		}
	for(int i = 0; i < ak - 2; i++){
		conVec[i] = - db[i+1]/dx[i+1] + db[i]/dx[i];
		setMatValue = 1/dx[i];
		Amat[i][i] = setMatValue;
		setMatValue = -1/dx[i] - 1/dx[i+1];
		Amat[i+1][i] = setMatValue;
		setMatValue  = 1 / dx[i+1];
		Amat[i+2][i] = setMatValue;
	}
}

void LCBase::recenterBeta(){
	double maxVal = -INFINITY;
	int k = b.size();
	for(int i = 0; i < k; i++)
		maxVal = max(maxVal, b[i]);
	for(int i = 0; i < k; i++)
		b[i] = b[i] - maxVal;
}

void LCBase::makePropDist(){
	int k = x.size();
//	double dx;
	double db;
//	vector<double> s;
//	s.resize(k);
	s[0] = 0;
//	double der;
	if(x[0] > -INFINITY){
		if(b[1] == -INFINITY | b[0] == -INFINITY){
				s[1] = s[0];
			} else {
//			dx = x[1] - x[0];
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
			s[1] = exp(b[1])/db;
			}
		}
	for(int i = 1; i < k-2; i++){
		if(b[i+1] == -INFINITY | b[i] == -INFINITY){
			s[i+1] = s[i];
			continue;
		} 
		//dx = x[i+1] - x[i];
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
		//	dx = x[k-1] - x[k-2];
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
			s[k-1] = s[k-2] - exp(b[k-2])/db;
		}					
	}		
	
	double logIntegral = log(s[k-1]);
	for(int i = 0; i < k; i++)
		b[i] = b[i] - logIntegral;	
}



void LogConCen::baseDervs2ActDervs() {
	int k = x.size();
	int a_k = getAK();
	int cur_l = actIndex[0];
	int cur_r = actIndex[1];
	int act_cnt = 1; 
	for(int i = 0; i < k; i++)
		allActDervs[i] = 0;
	double l_sum = 0;
	double r_sum = 0;
	bool New_Act = true;
	

	for(int i = cur_l+1; i < cur_r; i++)
		r_sum = r_sum + allBaseDervs[i] * (x[cur_r] - x[i])/(x[cur_r] - x[cur_l]);
	allActDervs[cur_l] = r_sum + allBaseDervs[cur_l];
	for(int i = actIndex[0] + 1; i < actIndex[a_k-1]; i++)
		{

		l_sum = (l_sum + allBaseDervs[i-1]) * (x[i-1] - x[cur_l])/(x[i]-x[cur_l]) ;
		r_sum = r_sum  * (x[cur_r] - x[i-1])/(x[cur_r] - x[i]) - allBaseDervs[i];
		if(New_Act == true)
			{
			New_Act = false;
			cur_l = actIndex[act_cnt-1];
			l_sum = 0;
			if(cur_r == i)
				r_sum = 0;
			} 
		if(i == cur_r)
			{
			act_cnt++;
			if(act_cnt == a_k)
				Rprintf("Warning: act_cnt == a_k in allActDervs!\n");
			cur_r = actIndex[act_cnt];
			r_sum = 0;
			for(int j = i + 1; j < cur_r; j++)
				{
					r_sum = r_sum + allBaseDervs[j] * (x[cur_r] - x[j])/(x[cur_r] - x[i]);
				}	
			New_Act = true;
			}
			allActDervs[i] = l_sum + r_sum + allBaseDervs[i];	
		}
	int k_last = actIndex[a_k-1];
	l_sum = (l_sum + allBaseDervs[k_last-1]) * (x[k_last-1] - x[cur_l])/(x[k_last]-x[cur_l]);
	allActDervs[k_last] = l_sum + allBaseDervs[k_last];
	
	
	if(x[0] == -INFINITY)
		allActDervs[0] = 0;
	if(x[k-1] == INFINITY)
		allActDervs[k-1] = 0;
}

void LogConCen::calcDervVec(){
	calcBaseDervs();
	baseDervs2ActDervs();
}
