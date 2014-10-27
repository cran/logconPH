#include "baseClass.h"

/*
double abs(double x){
	if(x < 0)
		return(-x);
	return(x);
}
*/

double min(double x, double y){
	if(x < y)
		return(x);
	return(y);
}

bool checkQPconstraints(QuadProgPP::Matrix<double> Amat, QuadProgPP::Vector<double> conVec, 
						QuadProgPP::Vector<double> propVec){
	int I = conVec.size();
	int J = propVec.size();
	double tot;
	for(int i = 0; i < I; i++){
		tot = 0;
		for(int j = 0; j < J; j++)
			tot = tot + Amat[j][i] * propVec[j];
		if(tot + conVec[i] < 0){
			return(false);
		}
	}
	return(true);
}

template <class T>
void printVec (vector<T> &vec){
//	Rprintf("C code called printVec, but not currently supported\n");
}

template <class T>
void printVec (vector<T> &vec, int begin, int end){
//	Rprintf("C code called printVec, but not currently supported\n");
}



/*
actSetBase::actSetBase(int MoveX, 
				   vector<double> X, 
				   vector<double> B, 
				   vector<int> L_ind, 
				   vector<int> R_ind,
				   vector<int> actInds){
				   x = X;
				   b = B;
				   Lindex = L_ind;
				   Rindex = R_ind;
				   actIndex = actInds;
				   allowMoveX = MoveX;
				   slpTol = pow(10, -13);
				   }
*/

void actSetBase::move_act_b(int index, double h) {
	int a_k = actIndex.size();
//	int k = b.size();
	b[index] = b[index] + h;
	int ind_l = -1;
	int ind_r = -1;
	double leftSlp = R_PosInf;
	double rightSlp = R_NegInf;
	if(a_k == 2 && index == actIndex[0])
		{
		leftSlp = (b[actIndex[1]] - b[actIndex[0]]) / (x[actIndex[1]] - x[actIndex[0]]);
		for(int i = actIndex[1]; i >= actIndex[0]; i--)
			b[i] = b[actIndex[1] ] + leftSlp * (x[i] - x[actIndex[1]]);
		return;
		}
	
	if(a_k == 2 && index == actIndex[1])
		{
		leftSlp = (b[actIndex[1] ] - b[actIndex[0]]) / (x[actIndex[1]] - x[actIndex[0]]);
		for(int i = actIndex[0]; i <= actIndex[1]; i ++)		
			b[i] = b[actIndex[0]] + leftSlp * (x[i] - x[actIndex[0] ]);
		return;
		}
	
	
	if(index > actIndex[0])
		{
		for(int i = 1; i < a_k; i++)
			{
			if(actIndex[i]	>= index )
				{
				ind_l = actIndex[i-1];
				break;
				}
			}
		if(ind_l < index && ind_l > -1)
			{
			leftSlp = (b[index] - b[ind_l])/(x[index] - x[ind_l]);
			
			for(int i = ind_l; i < index; i++) 
				{
				b[i] = b[ind_l] + leftSlp * (x[i] - x[ind_l]);
				}
			}
		}
	if(index < actIndex[a_k-1])
		{
		for(int i = a_k - 2; i >= 0; i--)
			{
			if(actIndex[i] <= index)
				{
				ind_r = actIndex[i+1];
				break;
				}
			}
	 	if(ind_r > index)
	 		{
	 		rightSlp = (b[ind_r] - b[index])/(x[ind_r] - x[index]);
	 		for(int i = index; i < ind_r; i++)
	 			{
	 			b[i] = b[index] + rightSlp * (x[i] - x[index]);
	 			}
	 		}
		}
}	

void actSetBase::addActive(int index){
	int ak = getAK();
	vector<int>::iterator iter;
	for(int i = 0 ; i < ak ; i++){
		if(index == actIndex[i]){
			return;
			}
		if(index < actIndex[i]){
			iter = actIndex.begin() + i;
			actIndex.insert(iter, index);
			return;
		}
	}
}

void actSetBase::removeActive(int index){
	int ak = getAK();
	vector<int>::iterator iter;
	for(int i = 0 ; i < ak  ; i++){
		if(actIndex[i] > index)
			return;
		if(actIndex[i] == index){
			iter = actIndex.begin() + i;
			actIndex.erase(iter);
			return;
		}	
	}
}

void actSetBase::flatten(int beg, int end){
	double slp = (b[end] - b[beg]) / (x[end] - x[beg]);
	double inter = b[beg];
	for(int i = beg + 1; i < end; i++)
		b[i] = inter + slp * (x[i] - x[beg]);
}

bool actSetBase::local3OK(int a_index){
	vector<double> limits(2);
	if(a_index > 0){
		limits = getLimits(a_index-1);
		if(limits[0] > limits[1]){
			return false;
		}
	}
	limits = getLimits(a_index);
		if(limits[0] > limits[1]){
			return false;
		}
	if(a_index < getAK() - 1){
		if(limits[0] > limits[1]){
			return false;
		}
	}
	return true;
}
				   
vector<double> actSetBase::numericDervs(int index){
	vector<double> ders;
	ders.resize(2);
	double lk_l, lk_h, lk_0;
	move_act_b(index, h);
	lk_h = llk();
	move_act_b(index, -2*h);
	lk_l = llk();
	move_act_b(index, h);
	lk_0 = llk();
	ders[0] = (lk_h - lk_l)/(2.0*h);
	ders[1] = (lk_h + lk_l - 2.0*lk_0)/(h*h);
	return(ders);
}

double actSetBase::numericParitialDerv(int index1, int index2){
	double lk_hh, lk_ll, lk_lh, lk_hl, pder;
	move_act_b(index1, h);
	move_act_b(index2, h);
	lk_hh = llk();
	move_act_b(index1, -2.0*h);
	lk_lh = llk();
	move_act_b(index2, -2.0*h);
	lk_ll = llk();
	move_act_b(index1, 2.0*h);
	lk_hl = llk();
	move_act_b(index1, -h);
	move_act_b(index2, h);
	pder = (lk_hh + lk_ll - lk_lh - lk_hl)/(4*h);
	return(pder);
}


// Note: this function only allows movement within an interval (although it will jump to end point if the exact value is given
// Will not move actIndex[0] or actIndex[ak-1]
// Will not move inactive points 
// Will not move points with odd index
void actSetBase::moveX(int index, double h){
	if(index % 2 == 0){
	//	Rprintf("Attempting to move fixed point \n");
		return;
	}
	int ak = getAK();
	if(index == actIndex[0] || index == actIndex[ak-1]){
	//	Rprintf("Attempting to move endpoints \n");
		return;
	}
	bool isActive = false;
	for(int i = 0; i < ak; i++){
		if(actIndex[i] == index){
			isActive = true;
			break;
		}
	}
	if(isActive == false){
		return;
	}
	if(h > x[index+1] - x[index] || h < x[index-1] - x[index]){
		return;
	}
	if(h < x[index+1] - x[index] && h > x[index-1] - x[index]) {
		x[index] = x[index] + h;
		move_act_b(index, 0);
		return;
	}
	if(h == x[index+1] - x[index]){
		b[index+1] = b[index];
		x[index] = (x[index-1] + x[index+1])/2; 
		removeActive(index);
		addActive(index+1);
		move_act_b(index + 1, 0);
		return;
	}
	if(h == x[index-1] - x[index]){
		b[index-1] = b[index];
		x[index] = (x[index-1] + x[index+1])/2; 
		removeActive(index);
		addActive(index-1);
		move_act_b(index-1, 0);
	}
}

double actSetBase::lk_xb(int index, double x_h, double b_h){
	move_act_b(index, b_h);
	moveX(index, x_h);
	double output = llk();
	move_act_b(index, -b_h);
	moveX(index, -x_h);
	return(output);	
}

vector<double> actSetBase::dervMoveX(int index){
	vector<double> output;
	output.resize(5);
	double h_xl = min(h, (x[index] - x[index-1])/2);
	double h_xh = min(h, (x[index+1] - x[index])/2);


	double lk_xh_bh = lk_xb(index, h_xh, h);
	double lk_xh_bl = lk_xb(index, h_xh, -h);
	double lk_xl_bh = lk_xb(index, -h_xl, h);
	double lk_xl_bl = lk_xb(index, -h_xl, -h);
	double lk_xh_b0 = lk_xb(index, h_xh, 0);
	double lk_xl_b0 = lk_xb(index, -h_xl, 0);
	double lk_x0_bh = lk_xb(index, 0, h);
	double lk_x0_bl = lk_xb(index, 0, -h);
	double lk_00 = llk();


	double b_d = (lk_x0_bh - lk_x0_bl) / (2 * h);
	double b_d2 = (lk_x0_bh + lk_x0_bl - 2 * lk_00) / pow(h,2);
	
	double x_d = (lk_xh_b0 - lk_xl_b0) / (h_xh + h_xl);
	double x_d2 = (lk_xh_b0 + lk_xl_b0 - 2 * lk_00) / (h_xl * h_xh);
	
	double dxdb = (lk_xh_bh + lk_xl_bl - lk_xh_bl - lk_xl_bh)/ (2 * h *(h_xl+h_xh) );
	//	Need dxdb...
	
	
	output[0] = b_d;
	output[1] = b_d2;
	output[2] = x_d;
	output[3] = x_d2;
	output[4] = dxdb;
	return(output);
}


void actSetBase::updateXs(){
	for(int i = 1; i < (getAK()- 1); i++){
		if(i < getAK() - 1 && actIndex[i] % 2 == 1)
			updateX(i);
	}
	int ai;
	for(int i = 1; i < (getAK() - 1); i++){
		if(i < getAK() - 1 && actIndex[i] % 2 == 1){
			ai = actIndex[i];
			if(x[ai] - x[ai -1] < pow(10.0, -8.0) ){
				x[ai] = (x[ai+1] + x[ai-1])/2;
				addActive(ai-1);
				removeActive(ai);
				move_act_b(actIndex[i], 0);
			}
			if(x[ai+1] - x[ai] < pow(10.0, -8.0) ){
				x[ai] = (x[ai+1] + x[ai-1])/2;
				addActive(ai+1);
				removeActive(ai);
				move_act_b(actIndex[i], 0);			
			}
		}
	}
}


void actSetBase::updateX(int a_index){
	int index = actIndex[a_index];
	if(index % 2 == 0){
	//	Rprintf("Warning: attempted to update odd x position\n");
	}	
	vector<double> dervs = dervMoveX(index);
	QuadProgPP::Matrix<double> Hess(2,2);
	Hess[0][0] = dervs[1];
	Hess[1][1] = dervs[3];
	Hess[1][0] = dervs[4];
	Hess[0][1] = dervs[4];
	double llk_old = llk();
	double llk_new;
	dervs[1] = dervs[2];
	dervs.resize(2);
	if(Hess[0][0] < 0 && Hess[1][1] < 0 && Hess[0][0] * Hess[1][1] > Hess[1][0]*Hess[0][1]){
		cur_Err = max(cur_Err, abs(dervs[0]) );
		double det = Hess[0][0] * Hess[1][1] - Hess[1][0] * Hess[0][1];
		QuadProgPP::Matrix<double> InvHess(2,2);
		InvHess[0][0] = 1/det * Hess[1][1];
		InvHess[1][1] = 1/det * Hess[0][0];
		InvHess[0][1] = -1/det * Hess[1][0];
		InvHess[1][0] = -1/det * Hess[0][1];
		vector<double> propVec(2);
		propVec[0] = -(dervs[0] * InvHess[0][0] + dervs[1] * InvHess[1][0]);	// prop step for b
		propVec[1] = -(dervs[0] * InvHess[0][1] + dervs[1] * InvHess[1][1]);	// prop step for x
		double scale = 0;
		if(propVec[1] > 0)
			scale = (x[index+1] - x[index])/propVec[1] * 4/5;
	 	if(propVec[1] < 0)
			scale = abs((x[index-1] - x[index])/propVec[1] * 4/5);
		
		if(scale < 1){
			propVec[0] = propVec[0]*scale;
			propVec[1] = propVec[1]*scale;
		}
		
		if(propVec[0] != propVec[0]){
//			Rprintf("propVec[0] undefined. quiting updateX\n");
			return;
		}
		if(propVec[1] != propVec[1]){
//			Rprintf("propVec[1] undefined. quitting updateX\n");
			return;
		}
		move_act_b(index, propVec[0]);
		moveX(index, propVec[1]);
		
		llk_new = llk();
		if( (llk_new < llk_old) || (!local3OK(a_index)) ){
			propVec[0] = -propVec[0]/2;
			propVec[1] = -propVec[1]/2;
			int it = 0;
			while( (it < 10) && (llk_new < llk_old || !local3OK(a_index) ) ){
				it++;
				move_act_b(index, propVec[0]);
				moveX(index, propVec[1]);
				propVec[0] = propVec[0]/2;
				propVec[1] = propVec[1]/2;
				llk_new = llk();
			}
			
			if( (llk_new < llk_old) || (!local3OK(a_index)) ){
		//		Rprintf("note: something bad happened in updateX\n");
		//		Rprintf("propVec[0] = %f, propVec[1] = %f\n", propVec[0], propVec[1]);
				move_act_b(index, propVec[0] * 2);
				moveX(index, propVec[1] * 2);
			}
		}
	}	
	else{	
		cur_Err = max(cur_Err, abs(dervs[0]) );
		double delta = -dervs[0]/Hess[0][0];
		delta = min(delta, 9/10 * (x[index + 1] - x[index]));
		delta = max(delta, 9/10 * (x[index - 1] - x[index]));
		moveX(index, delta);
		llk_new = llk();
		if( (llk_new < llk_old) || !local3OK(a_index)){
			delta = -delta/2;
			int it = 0;
			while(it < 5 && (llk_new < llk_old || !local3OK(a_index))){
				it++;
				moveX(index, delta);
				delta = delta/2;
				llk_new = llk();
				}
			if(llk_new > llk_old){
				moveX(index, delta * 2);
			}
		}
		if(x[index]-x[index - 1] < 0.0001){
			if(b[index-1] == R_NegInf)
				b[index-1] = b[index];
			addActive(index-1);
			removeActive(index);
			x[index] = (x[index-1] + x[index+1])/2;
			b[index] = (b[index-1] + b[index+1])/2;
		}
		if(x[index+1]-x[index] < 0.0001){
			if(b[index+1] == R_NegInf)
				b[index+1] = b[index];
			addActive(index+1);
			removeActive(index);
			x[index] = (x[index+1] + x[index+1])/2;
			b[index] = (b[index+1] + b[index+1])/2;
		}		
	}
}


void actSetBase::update1Var(int index) {
	int minInd = actIndex[0];
	int maxInd = actIndex[getAK()-1];
	if(index > maxInd){
		Rprintf("Warning: invalid (large) index selected in update1Var\n");
		return;
	}
	if(index < minInd){
		Rprintf("Warning: invalid (small) index selected in update1Var\n");
		return;
	}
	addActive(index);
	vector<double> ders = numericDervs(index);
	vector<double> lims = getLimits(index);
	int its = 0;
	double loglike = llk();
	if(ders[1] < 0){
		double delta = min( max( lims[0], -ders[0]/ders[1] ), lims[1] );
		move_act_b(index, delta);
		double new_loglike = llk();
		
		delta = delta * -1;
		while( (new_loglike < loglike) && (its < 10)){
			its++;
			delta = delta/2;
			move_act_b(index, delta);
			new_loglike = llk();
		}
		if(new_loglike < loglike){
			move_act_b(index, delta);
			}
	} else {
		double delta;
		if(ders[0] < 0 )
			delta = max(lims[0], -1.0);
		else 
			delta = min(lims[1],1.0);
		
		move_act_b(index, delta);
		double new_loglike = llk();
		delta = delta * - 1;
		while((new_loglike < loglike) && (its < 25)) {
			its ++;
			delta = delta/2;
			move_act_b(index, delta);
			new_loglike = llk();
		}
	}
	checkAllActive();
}




void actSetBase::VEMstep(){
/*	calcDervVec();
	int updateIndex = findMaxError();
	if(updateIndex == 0)
		return;
	update1Var(updateIndex);
*/
	calcDervVec();
	vector<int> points = findMaxIntError();
	int pointsSize = points.size();
	int thisIndex;
	for(int i = 0; i < pointsSize;i++){
		thisIndex = points[i];
		if(thisIndex > 0)
			update1Var(thisIndex);
	}	
}

void actSetBase::ICMstep(){
	double parD;
	int ak = getAK();
	
	checkAllActive();
	
	QuadProgPP::Vector<double> d1(ak);
	QuadProgPP::Vector<double> propStep(ak);
	QuadProgPP::Vector<double> blankEqs(0);
	QuadProgPP::Matrix<double> blankMat(ak, 0);
	QuadProgPP::Matrix<double> Amat(ak, ak-2);
	QuadProgPP::Vector<double> conVec(ak-2);
	QuadProgPP::Matrix<double> ParHess(ak, ak);
	
	for(int i = 0; i < ak; i++){
		if(i < ak-2)
			conVec[i] = 0.0;
		propStep[i] = 0.0;
		for(int j = 0; j < ak; j++){
			if(j < ak - 2)
				Amat[i][j] = 0.0;
			ParHess[i][j] = 0.0;
		}
	}
	vector<double> currentDervs(2);

		if(ak <= 2){
			vector<int> actIndsCopy = actIndex;
			for(int i = 0; i < ak; i++)
				update1Var(actIndsCopy[i]);
//			checkAllActive();
			return;
			}
		for(int i = 0; i < ak; i++){
			currentDervs = numericDervs(actIndex[i]);
			d1[i] = -currentDervs[0];
			if(currentDervs[1] >= 0){
				vector<int> curInds = actIndex;
				int indSize = curInds.size();
				for(int j = 0; j < indSize; j++){
					update1Var(curInds[j]);
//					checkAllActive();
					}
				return;
				}
			ParHess[i][i] = -currentDervs[1];
			}
		for(int i = 0; i<(ak-1); i++){
			parD = numericParitialDerv(actIndex[i], actIndex[i+1]);
			if(parD * parD < ParHess[i][i] * ParHess[i+1][i+1]){
				ParHess[i][i+1] = -parD;
				ParHess[i+1][i] = -parD;
			}
		}	
		qpLimMatrix(Amat, conVec);	
		double inf = 1.0E300;

		double QP_result = QuadProgPP::solve_quadprog(ParHess, d1, blankMat, blankEqs, Amat, conVec, propStep);	
		if(QP_result == inf){
		//	Rprintf("Warning: no feasible solution to QP problem! Skipping ICM step\n");
			return;
		}	
		double str_llk = llk();
		vector<double> limits(2);
		bool propOK = true;
		for(int i = 0; i < ak; i++){
			propOK = false;
			if(propStep[i] < 10000.0 && propStep[i] > -10000.0)
				propOK = true;
		}
		if(propOK == false){
		//	Rprintf("note: propStep is really crazy. skipping\n");
			return;
		}
		for(int i = 0; i < ak; i++)
			move_act_b(actIndex[i], propStep[i]);
		for(int j = 0; j < ak; j++){
			limits = getLimits(actIndex[j]);
			if(limits[0] > limits[1]){
				propOK = false;
				break;
				}
			}	
		if(propOK){
			double new_llk = llk();
			if(new_llk < str_llk){
				for(int i = 0; i < ak; i++)
					propStep[i] = propStep[i] * -0.5;
				int it = 0;
				while(str_llk > new_llk && it < 15){
					it++;
					for(int i = 0; i < ak; i++){
						move_act_b(actIndex[i], propStep[i]);
						propStep[i] = propStep[i]/2;
					}
					new_llk = llk();	
					}
				}
			if(new_llk < str_llk){
				for(int i = 0; i < ak; i++){
					move_act_b(actIndex[i], propStep[i] * 2);
				}
			}
		checkAllActive();
		}
		if(!propOK){
			for(int i = 0; i<ak ;i++)
				move_act_b(actIndex[i], -propStep[i]);
		}
		vector<int> curInds = actIndex;
		int curIndSize = curInds.size();
		for(int i = 0; i < curIndSize; i++){
			update1Var(curInds[i]);
			checkAllActive();
		}
}

	
