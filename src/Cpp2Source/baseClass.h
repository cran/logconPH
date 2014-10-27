#include <iostream>
#include <vector>
#include "QuadProg++.cc"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>



using namespace std;


class actSetBase {
	public:
	vector<double> x;
	vector<double> dx;
	vector<double> b;
	vector<int> Lindex;
	vector<int> Rindex;
	vector<int> actIndex;
	vector<double> rep_vec;
	vector<double> allActDervs;
	vector<double> s;
	
	virtual void calcDervVec() = 0;
	
	int getAK() {return(actIndex.size() ); };
	double cur_llk;
	double cur_Err;
	double inner_Err;
	double lk_xb(int index, double x_h, double b_h);
	void move_act_b(int index, double h);
	vector<double> numericDervs(int index);
	double numericParitialDerv(int index1, int index2);
	void VEMstep();
	virtual int findMaxError() = 0;
	virtual vector<int> findMaxIntError() = 0;
	int allowMoveX;
	void moveX(int index, double h);
	vector<double> dervMoveX(int index);
	void addActive(int index);
	void removeActive(int index);
	double slpTol;
	double h;
	double h2;
	void update1Var(int index);
	void updateX(int index);
	void updateXs();
//	actSetBase::actSetBase(int, vector<double>, vector<double>, vector<int>, vector<int>, vector<int>);
	virtual double llk() = 0;
	virtual vector<double> getLimits(int index) = 0;
	void flatten(int beg, int end);
	virtual void checkAllActive() = 0;
	bool local3OK(int a_index);
	virtual void recenterBeta() = 0;
	virtual void makePropDist() = 0;
	virtual void qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec) = 0;
//	void updateActDervs(){calcBaseDervs(); baseDervs2ActDervs();};	
	void ICMstep();
	};

class LCBase : public actSetBase{
	public:
	double endTol;
	vector<double> getLimits(int index);
	void checkEnds();
	vector<double> d_bl;
	vector<double> d_bu;
	void checkAllActive();
	int findMaxError();
	vector<int> findMaxIntError();
	void qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec);
	void recenterBeta();
	void makePropDist();
	vector<double> allBaseDervs;	
	virtual void calcDervVec() = 0;
};

class LogConCen : public LCBase {
	public:
	double llk();
	LogConCen(int, vector<double>, vector<double>, vector<int>, vector<int>, vector<int>, vector<double>, bool);
	vector<double> allBaseDervs;
	void calcBaseDervs();
	void baseDervs2ActDervs();	
	void calcDervVec();
};		


class LogConCenPH : public LCBase {
	public:
	double llk();
	double nullk();
	LogConCenPH(int MoveX, vector<double> X, vector<double> B, 
				vector<int> L_ind, vector<int> R_ind, vector<int> actInds,
				bool allow_x_move, double augl, double augr, QuadProgPP::Matrix<double> Covars,
				int num_cov);
//	vector <double> nu;
    void updateNu();
	double	augLeft, augRight;
	void calcDervVec();
    void calcNuDerv();
    double updateRegress();
    QuadProgPP::Matrix<double> covars;
    QuadProgPP::Vector<double> cov_b;
    QuadProgPP::Vector<double> nu;
    
    QuadProgPP::Vector<double> cov_b_d;
    QuadProgPP::Vector<double> b_h;
    QuadProgPP::Vector<double> b_l;
    QuadProgPP::Matrix<double> Hess_b;
    QuadProgPP::Matrix<double> test_Hess_b;

	vector<double> p;
	void update_p(int index1, int index2);
	void p2s();
	double fastBasellk();
	void fastNumActDers();
	double scaleValue;
	void moveCovOrBase(int i, double delta);
	double partialDerCovOrBase(int i, int j);
	
	double updateOnFlyLoop();
};
