#ifndef MODELSTATS_HPP
#define MODELSTATS_HPP
#include <cstdlib>
#include <cmath>


#if (! NOGSL)
#include "gsl/gsl_sf.h"
#else
double lgamma(double x) {
	double ans;
	ans = x*log(x)-x;
	return ans;
};
double gsl_sf_lnbeta(double a, double b) {
	double ans;
	ans = lgamma(a)+lgamma(b) -lgamma(a+b);
	return ans;
};
double gsl_sf_gamma(double a) {
	return exp(lgamma(a));
};

#endif

float lnBetaFunction(float a, float b) {
	return gsl_sf_lnbeta(a,b);
};



float mySafeLog(float x) {
	if (x<1e-8) return 1e10;
	else if (std::isnan(x)) return 0;
	else return log(x);
};


class ModelPairStatsBase;

class ModelSelfStatsBase {
	public:
		virtual ModelSelfStatsBase* Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair)=0;
		virtual void AddOutGoing(float x) {}
};

class BinomialSelfStats: public ModelSelfStatsBase {
	public:
		float nV;
		BinomialSelfStats() {
			this->nV = 1;
		}
		virtual BinomialSelfStats* Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair);
};

class WSelfStats: public ModelSelfStatsBase {
	public:
		float degree;
		float selfMissing;
		virtual WSelfStats* Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair);
		virtual void AddOutGoing(float x);
		WSelfStats() {
			this->degree = 0;
			this->selfMissing = 0;
		}
};

class PoissonSelfStats: public ModelSelfStatsBase {
	public:
		float nV;
		PoissonSelfStats() {
			nV =1;
		}
		virtual PoissonSelfStats* Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair);
};

class ModelPairStatsBase {
	public:
		virtual float simple()=0;
		virtual ModelPairStatsBase* Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3)=0;
		virtual ModelPairStatsBase* Add2(ModelPairStatsBase* b2)=0;
		virtual void AddEdge(float x)=0;
		virtual float MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb)=0;
		virtual float MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx)=0;
		virtual float FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb)=0;
		virtual float FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx)=0;
};

class BinomialPairStats: public ModelPairStatsBase {
	public:
		float nE;
		BinomialPairStats() {
			this->nE = 0;
		}
		virtual float simple() { return nE; }
		virtual BinomialPairStats* Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3);
		virtual BinomialPairStats* Add2(ModelPairStatsBase* b2);
		virtual void AddEdge(float x);
		virtual float MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb);
		virtual float MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx);
		virtual float FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb);
		virtual float FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx);
};

class PoissonPairStats: public ModelPairStatsBase {
	public:
		float sumE;
		float sumlgam;
		PoissonPairStats() {
			this->sumE = 0; this->sumlgam = 0;
		}
		virtual float simple() { return sumE; }
		virtual PoissonPairStats* Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3);
		virtual PoissonPairStats* Add2(ModelPairStatsBase* b2);
		virtual void AddEdge(float x);
		virtual float MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb);
		virtual float MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx);
		virtual float FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb);
		virtual float FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx);
};

class WPairStats: public ModelPairStatsBase {
	public:
		float nE;
		WPairStats() {
			this->nE = 0;
		}
		virtual float simple() { return nE; }
		virtual WPairStats* Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3);
		virtual WPairStats* Add2(ModelPairStatsBase* b2);
		virtual void AddEdge(float x);
		virtual float MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb);
		virtual float MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx);
		virtual float FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb);
		virtual float FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx);
};

BinomialSelfStats* BinomialSelfStats::Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair) {
	BinomialSelfStats* p = new BinomialSelfStats;
	p->nV = this->nV + ((BinomialSelfStats*) b2)->nV;
	return p;
}


PoissonSelfStats* PoissonSelfStats::Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair) {
	PoissonSelfStats* p = new PoissonSelfStats;
	p->nV = this->nV + ((PoissonSelfStats*) b2)->nV;
	return p;
}


WSelfStats* WSelfStats::Add2(ModelSelfStatsBase* b2,ModelPairStatsBase* pair) {
	WSelfStats* p = new WSelfStats;
	p->degree = this->degree + ((WSelfStats*) b2)->degree;
	p->selfMissing = this->selfMissing + ((WSelfStats*) b2)->selfMissing + (this->degree * ((WSelfStats*) b2)->degree);
	if (pair != NULL) {
		p->selfMissing  -= ((WPairStats* )pair)->nE;
	}
	return p;
}

void WSelfStats::AddOutGoing(float x) {
	this->degree = this->degree + x;
}

float BinomialPairStats::FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb) {
	float Hab, Eab, Eaa, Ebb, Haa, Hbb, ans;
	Eab = this->nE;
	if (aa==NULL) Eaa = 0;
	else Eaa = ((BinomialPairStats* )aa)->nE;
	if (bb==NULL) Ebb = 0;
	else Ebb = ((BinomialPairStats* )bb)->nE;
	Hab = ((BinomialSelfStats*) sa)->nV * ((BinomialSelfStats*) sb)->nV - Eab;
	Haa = (((BinomialSelfStats*) sa)->nV * ((BinomialSelfStats*) sa)->nV - 1)/2.0f;
	Hbb = (((BinomialSelfStats*) sb)->nV * ((BinomialSelfStats*) sb)->nV - 1)/2.0f;
	ans = lnBetaFunction(Eaa+Ebb+Eab+1,Haa+Hbb+Hab+1)
	    - lnBetaFunction(Eaa+1,Haa+1)
	    - lnBetaFunction(Ebb+1,Hbb+1)
	    - lnBetaFunction(Eab+1,Hab+1);
	return ans;
}

float BinomialPairStats::MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb) {
	float Tab, Eab, Eaa, Ebb, Taa, Tbb, ans;
	Tab = ((BinomialSelfStats*) sa)->nV * ((BinomialSelfStats*) sb)->nV;
	Taa = (((BinomialSelfStats*) sa)->nV * ((BinomialSelfStats*) sa)->nV - 1)/2.0f;
	Tbb = (((BinomialSelfStats*) sb)->nV * ((BinomialSelfStats*) sb)->nV - 1)/2.0f;
	Eab = this->nE;
	if (aa==NULL) Eaa = 0;
	else Eaa = ((BinomialPairStats* )aa)->nE;
	if (bb==NULL) Ebb = 0;
	else Ebb = ((BinomialPairStats* )bb)->nE;
	float Etot = Eaa+Ebb+Eab;
	float Ttot = Taa+Tbb+Tab;
	float thetaTot = Etot/Ttot;
	float thetaA = Eaa/Taa;
	float thetaB = Ebb/Tbb;
	float thetaAB = Eab/Tab;
	ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
		- Eaa * mySafeLog(thetaA) - (Taa - Eaa) * mySafeLog(1-thetaA)
		- Ebb * mySafeLog(thetaB) - (Tbb - Ebb) * mySafeLog(1-thetaB)
		- Eab * mySafeLog(thetaAB) - (Tab - Eab) * mySafeLog(1-thetaAB);
	return ans;
}


float BinomialPairStats::MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx) {
	float Eax, Ebx, Tax, Tbx, ans, da, db, dx;
	if (pax==NULL) Eax = 0;
	else Eax = ((BinomialPairStats*) pax)->nE;
	if (pbx==NULL) Ebx = 0;
	else Ebx = ((BinomialPairStats*) pbx)->nE;
	da = ((BinomialSelfStats*) sa)->nV;
	db = ((BinomialSelfStats*) sb)->nV;
	dx = ((BinomialSelfStats*) sx)->nV;
	Tax = da * dx;
	Tbx = db * dx;
	float Etot = Eax+Ebx;
	float Ttot = Tax+Tbx;
	float thetaTot = Etot/Ttot;
	float thetaA = Eax/Tax;
	float thetaB = Ebx/Tbx;
	ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
		- Eax * mySafeLog(thetaA) - (Tax - Eax) * mySafeLog(1-thetaA)
		- Ebx * mySafeLog(thetaB) - (Tbx - Ebx) * mySafeLog(1-thetaB);
	return ans;
}


float BinomialPairStats::FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx) {
	float Eax, Ebx, Tax, Tbx, ans, da, db, dx;
	if (pax==NULL) Eax = 0;
	else Eax = ((BinomialPairStats*) pax)->nE;
	if (pbx==NULL) Ebx = 0;
	else Ebx = ((BinomialPairStats*) pbx)->nE;
	da = ((BinomialSelfStats*) sa)->nV;
	db = ((BinomialSelfStats*) sb)->nV;
	dx = ((BinomialSelfStats*) sx)->nV;
	Tax = da * dx;
	Tbx = db * dx;
	float Hax = Tax - Eax;
	float Hbx = Tbx - Ebx;
	ans = lnBetaFunction(Eax+Ebx+1.0f,Hax+Hbx+1.0f)
		- lnBetaFunction(Eax+1.0f,Hax+1.0f)
		- lnBetaFunction(Ebx+1.0f,Hbx+1.0f);
	return ans;
}

void BinomialPairStats::AddEdge(float x) {
	this->nE += x;
}


BinomialPairStats* BinomialPairStats::Add2(ModelPairStatsBase* b2) {
	BinomialPairStats* p = new BinomialPairStats;
	p->nE = this->nE;
	if (b2 != NULL) p->nE += ((BinomialPairStats*) b2)->nE;
	return p;
}


BinomialPairStats* BinomialPairStats::Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3) {
	BinomialPairStats* p = new BinomialPairStats;
	if (b2 != NULL) p->nE += ((BinomialPairStats*) b2)->nE;
	if (b3 != NULL) p->nE += ((BinomialPairStats*) b3)->nE;
	return p;
}

float PoissonPairStats::FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb) {
	float Tab, Eab, Eaa, Ebb, lEab, lEbb, lEaa, Taa, Tbb, ans;
	Tab = ((PoissonSelfStats*) sa)->nV * ((PoissonSelfStats*) sb)->nV;
	Taa = (((PoissonSelfStats*) sa)->nV * ((PoissonSelfStats*) sa)->nV - 1);
	Tbb = (((PoissonSelfStats*) sb)->nV * ((PoissonSelfStats*) sb)->nV - 1);
	Eab = this->sumE;
	lEab = this->sumlgam;
	if (aa==NULL) { Eaa = 0; lEaa = 0; }
	else { Eaa = ((PoissonPairStats* )aa)->sumE; lEaa = ((PoissonPairStats* )aa)->sumlgam; }
	if (bb==NULL) { Ebb = 0; lEbb = 0; }
	else { Ebb = ((PoissonPairStats* )bb)->sumE; lEbb = ((PoissonPairStats* )bb)->sumlgam; }
	float Etot = Eaa+Ebb+Eab;
	float lEtot = lEaa + lEbb + lEab;
	float Ttot = Taa+Tbb+Tab;
	ans = Etot*(mySafeLog(Etot/Ttot)-1) - lEtot
		- Eaa*(mySafeLog(Eaa/Taa)-1) + lEaa
		- Eab*(mySafeLog(Eab/Tab)-1) + lEab
		- Ebb*(mySafeLog(Ebb/Tbb)-1) + lEbb;
	return ans;
}



float PoissonPairStats::MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb) {
	float Tab, Eab, Eaa, Ebb, lEab, lEbb, lEaa, Taa, Tbb, ans;
	Tab = ((PoissonSelfStats*) sa)->nV * ((PoissonSelfStats*) sb)->nV;
	Taa = (((PoissonSelfStats*) sa)->nV * ((PoissonSelfStats*) sa)->nV - 1);
	Tbb = (((PoissonSelfStats*) sb)->nV * ((PoissonSelfStats*) sb)->nV - 1);
	Eab = this->sumE;
	lEab = this->sumlgam;
	if (aa==NULL) { Eaa = 0; lEaa = 0; }
	else { Eaa = ((PoissonPairStats* )aa)->sumE; lEaa = ((PoissonPairStats* )aa)->sumlgam; }
	if (bb==NULL) { Ebb = 0; lEbb = 0; }
	else { Ebb = ((PoissonPairStats* )bb)->sumE; lEbb = ((PoissonPairStats* )bb)->sumlgam; }
	float Etot = Eaa+Ebb+Eab;
	float lEtot = lEaa + lEbb + lEab;
	float Ttot = Taa+Tbb+Tab;
	ans = Etot*(mySafeLog(Etot/Ttot) -1) - lEtot
		- Eaa*(mySafeLog(Eaa/Taa) -1) + lEaa
		- Ebb*(mySafeLog(Ebb/Tbb) -1) + lEbb
		- Eab*(mySafeLog(Eab/Tab) -1) + lEab;
	return ans;
}

float PoissonPairStats::MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx) {
	float Eax, Ebx, lEax, lEbx, Tax, Tbx, ans, da, db, dx;
	if (pax==NULL) { Eax = 0; lEax = 0;}
	else {
		Eax = ((PoissonPairStats*) pax)->sumE;
		lEax = ((PoissonPairStats*) pax)->sumlgam;
	}
	if (pbx==NULL) { Ebx = 0; lEbx = 0;}
	else {
		Ebx = ((PoissonPairStats*) pbx)->sumE;
		lEbx = ((PoissonPairStats*) pbx)->sumlgam;
	}
	da = ((PoissonSelfStats*) sa)->nV;
	db = ((PoissonSelfStats*) sb)->nV;
	dx = ((PoissonSelfStats*) sx)->nV;
	Tax = da * dx;
	Tbx = db * dx;
	float Etot = Eax+Ebx;
	float lEtot = lEax + lEbx;
	float Ttot = Tax+Tbx;
	ans = Etot*(mySafeLog(Etot/Ttot)-1) - lEtot
		- Eax*(mySafeLog(Eax/Tax)-1) + lEax
		- Ebx*(mySafeLog(Ebx/Tbx)-1) + lEbx;
	return ans;
}


float PoissonPairStats::FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx) {
	float Eax, Ebx, lEax, lEbx, Tax, Tbx, ans, da, db, dx;
	if (pax==NULL) { Eax = 0; lEax = 0;}
	else {
		Eax = ((PoissonPairStats*) pax)->sumE;
		lEax = ((PoissonPairStats*) pax)->sumlgam;
	}
	if (pbx==NULL) { Ebx = 0; lEbx = 0;}
	else {
		Ebx = ((PoissonPairStats*) pbx)->sumE;
		lEbx = ((PoissonPairStats*) pbx)->sumlgam;
	}
	da = ((PoissonSelfStats*) sa)->nV;
	db = ((PoissonSelfStats*) sb)->nV;
	dx = ((PoissonSelfStats*) sx)->nV;
	Tax = da * dx;
	Tbx = db * dx;
	float Etot = Eax+Ebx;
	float lEtot = lEax + lEbx;
	float Ttot = Tax+Tbx;
	ans = - lEtot + lgamma(1+Etot) - (1+Etot)*mySafeLog(Ttot)
	      + lEax  - lgamma(1+Eax)  + (1+Eax )*mySafeLog(Tax)
	      + lEbx  - lgamma(1+Ebx)  + (1+Ebx )*mySafeLog(Tbx);
	return ans;
}




void PoissonPairStats::AddEdge(float x) {
	this->sumE += x;
	this->sumlgam += lgamma(x+1);
}


PoissonPairStats* PoissonPairStats::Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3) {
	PoissonPairStats* p = new PoissonPairStats;
	p->sumE = this->sumE;
	p->sumlgam = this->sumlgam;
	if (b2 != NULL) {p->sumE += ((PoissonPairStats*) b2)->sumE; p->sumlgam += ((PoissonPairStats*) b2)->sumlgam; }
	if (b3 != NULL) {p->sumE += ((PoissonPairStats*) b3)->sumE; p->sumlgam += ((PoissonPairStats*) b3)->sumlgam; }
	//p->sumlgam = this->sumlgam + ((PoissonPairStats*) b2)->sumlgam + ((PoissonPairStats*) b3)->sumlgam;
	return p;
}


PoissonPairStats* PoissonPairStats::Add2(ModelPairStatsBase* b2) {
	PoissonPairStats* p = new PoissonPairStats;
	p->sumE = this->sumE;
	p->sumlgam = this->sumlgam;
	if (b2 != NULL) {p->sumE += ((PoissonPairStats*) b2)->sumE; p->sumlgam += ((PoissonPairStats*) b2)->sumlgam; }
	return p;
}


float WPairStats::FBcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb) {
	float Hab, Eab, Eaa, Ebb, Haa, Hbb, ans;
	Eab = this->nE;
	if (aa==NULL) Eaa = 0;
	else Eaa = ((WPairStats* )aa)->nE;
	if (bb==NULL) Ebb = 0;
	else Ebb = ((WPairStats* )bb)->nE;
	Hab = ((WSelfStats*) sa)->degree * ((WSelfStats*) sb)->degree - Eab;
	Haa = ((WSelfStats*) sa)->selfMissing;
	Hbb = ((WSelfStats*) sb)->selfMissing;
	ans = lnBetaFunction(Eaa+Ebb+Eab+1,Haa+Hbb+Hab+1)
	    - lnBetaFunction(Eaa+1,Haa+1)
	    - lnBetaFunction(Ebb+1,Hbb+1)
	    - lnBetaFunction(Eab+1,Hab+1);
	return ans;
}

float WPairStats::MLcenterscore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, ModelPairStatsBase* aa, ModelPairStatsBase* bb) {
	float Tab, Eab, Eaa, Ebb, Haa, Hbb, Taa, Tbb, ans;
	Eab = this->nE;
	if (aa==NULL) Eaa = 0;
	else Eaa = ((WPairStats* )aa)->nE;
	if (bb==NULL) Ebb = 0;
	else Ebb = ((WPairStats* )bb)->nE;
	Haa = ((WSelfStats*) sa)->selfMissing;
	Hbb = ((WSelfStats*) sb)->selfMissing;
	Tab = ((WSelfStats*) sa)->degree * ((WSelfStats*) sb)->degree;
	Taa = Eaa + Haa;
	Tbb = Ebb + Hbb;
	float Etot = Eaa+Ebb+Eab;
	float Ttot = Taa+Tbb+Tab;
	float thetaTot = Etot/Ttot;
	float thetaA = Eaa/Taa;
	float thetaB = Ebb/Tbb;
	float thetaAB = Eab/Tab;
	ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
		- Eaa * mySafeLog(thetaA) - (Taa - Eaa) * mySafeLog(1-thetaA)
		- Ebb * mySafeLog(thetaB) - (Tbb - Ebb) * mySafeLog(1-thetaB)
		- Eab * mySafeLog(thetaAB) - (Tab - Eab) * mySafeLog(1-thetaAB);
	return ans;
}


float WPairStats::MLdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx) {
	float Eax, Ebx, Tax, Tbx, ans, da, db, dx;
	if (pax==NULL) Eax = 0;
	else Eax = ((WPairStats*) pax)->nE;
	if (pbx==NULL) Ebx = 0;
	else Ebx = ((WPairStats*) pbx)->nE;
	da = ((WSelfStats*) sa)->degree;
	db = ((WSelfStats*) sb)->degree;
	dx = ((WSelfStats*) sx)->degree;
	Tax = da * dx;
	Tbx = db * dx;
	float Etot = Eax+Ebx;
	float Ttot = Tax+Tbx;
	float thetaTot = Etot/Ttot;
	float thetaA = Eax/Tax;
	float thetaB = Ebx/Tbx;
	ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
		- Eax * mySafeLog(thetaA) - (Tax - Eax) * mySafeLog(1-thetaA)
		- Ebx * mySafeLog(thetaB) - (Tbx - Ebx) * mySafeLog(1-thetaB);
	return ans;
}


float WPairStats::FBdeltascore(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb,ModelSelfStatsBase* sx, ModelPairStatsBase* pax, ModelPairStatsBase* pbx) {
	float Eax, Ebx, Hax, Hbx, ans, da, db, dx;
	if (pax==NULL) Eax = 0;
	else Eax = ((WPairStats*) pax)->nE;
	if (pbx==NULL) Ebx = 0;
	else Ebx = ((WPairStats*) pbx)->nE;
	da = ((WSelfStats*) sa)->degree;
	db = ((WSelfStats*) sb)->degree;
	dx = ((WSelfStats*) sx)->degree;
	Hax = da * dx - Eax;
	Hbx = db * dx - Ebx;
	ans = lnBetaFunction(Eax+Ebx+1.0f,Hax+Hbx+1.0f)
              -lnBetaFunction(Eax+1.0f,Hax+1.0f)
              -lnBetaFunction(Ebx+1.0f,Hbx+1.0f);
	return ans;
}


void WPairStats::AddEdge(float x) {
	this->nE += x;
}

WPairStats* WPairStats::Add3(ModelPairStatsBase* b2, ModelPairStatsBase* b3) {
	WPairStats* p = new WPairStats;
	p->nE = this->nE;
	if (b2 != NULL) {p->nE += ((WPairStats*) b2)->nE; }
	if (b3 != NULL) {p->nE += ((WPairStats*) b3)->nE; }
	return p;
}


WPairStats* WPairStats::Add2(ModelPairStatsBase* b2) {
	WPairStats* p = new WPairStats;
	p->nE = this->nE;
	if (b2 != NULL) {p->nE += ((WPairStats*) b2)->nE; }
	return p;
}

#endif
