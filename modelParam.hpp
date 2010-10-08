#ifndef MODELPARAM_HPP
#define MODELPARAM_HPP
#include "modelStats.hpp"

class ModelParamBase {
	public:
		virtual void calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb)=0;
		virtual void collapse(ModelParamBase* pa, ModelParamBase* pb)=0;
		virtual void cleanup()=0;
		virtual void bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof)=0;
		virtual bool isZero()=0;
		virtual float predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb)=0;
		virtual float updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar)=0;
		virtual void init()=0;
		virtual std::string DerivedType()=0;
};

class BinomialParam: public ModelParamBase {
	public:
		float p;
		float numE;
		float numT;
		virtual void calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual void collapse(ModelParamBase* pa, ModelParamBase* pb);
		virtual void cleanup();
		virtual void bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof);
		virtual bool isZero();
		virtual float predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual float updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar);
		virtual void init();
		virtual std::string DerivedType() { return std::string("BinomialParam"); }
};

void BinomialParam::init() {
	p = 0;
	numE = 0;
	numT = 0;
}

float BinomialParam::updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar) {
	return 1 - (1 - oldsofar)*(1 - this->predict(sa,sb));
}

void BinomialParam::bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof) {
	this->p = std::max(0.0f,(((BinomialParam*) Ori)->p - ((BinomialParam*) Sof)->p)/(1 - ((BinomialParam*) Sof)->p));
}

void BinomialParam::collapse(ModelParamBase* pa, ModelParamBase* pb) {
	this->numE += ((BinomialParam*) pa)->numE + ((BinomialParam*) pb)->numE;
	this->numT += ((BinomialParam*) pa)->numT + ((BinomialParam*) pb)->numT;
}

float BinomialParam::predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	return ((BinomialSelfStats*) sa)->nV * this->p * ((BinomialSelfStats*) sb)->nV;
}

bool BinomialParam::isZero() {
	return (p==0);
}

void BinomialParam::calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	if (pab==NULL) this->numE = 0;
	else this->numE = ((BinomialPairStats*) pab)->nE;
	this->numT = ((BinomialSelfStats*) sa)->nV * ((BinomialSelfStats*) sb)->nV;
}

void BinomialParam::cleanup() {
	p = numE/numT;
	if (std::isnan(p)) p=0;
}
	

class PoissonParam: public ModelParamBase {
	public:
		float lambda;
		float numT;
		float sumE;
		float sumlgam;
		virtual void calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual void collapse(ModelParamBase* pa, ModelParamBase* pb);
		virtual void cleanup();
		virtual void bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof);
		virtual bool isZero();
		virtual float predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual float updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar);
		virtual void init();
		virtual std::string DerivedType() { return std::string("PoissonParam"); }
};

void PoissonParam::init() {
	lambda = 0;
	numT = 0;
	sumE = 0;
	sumlgam = 0;
}

float PoissonParam::updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar) {
	return oldsofar + this->predict(sa,sb);
}

void PoissonParam::bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof) {
	//this->p = std::max(0.0f,(((PoissonParam*) Ori)->p - ((PoissonParam*) Sof)->p)/(1 - ((PoissonParam*) Sof)->p));
	this->lambda = std::max(0.0f,(((PoissonParam*) Ori)->lambda - ((PoissonParam*) Sof)->lambda));
}

void PoissonParam::collapse(ModelParamBase* pa, ModelParamBase* pb) {
	assert(pa != NULL);
	assert(pb != NULL);
	this->sumE += ((PoissonParam*) pa)->sumE + ((PoissonParam*) pb)->sumE;
	this->sumlgam += ((PoissonParam*) pa)->sumlgam + ((PoissonParam*) pb)->sumlgam;
	this->numT += ((PoissonParam*) pa)->numT + ((PoissonParam*) pb)->numT;
}

float PoissonParam::predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	return ((PoissonSelfStats*) sa)->nV * this->lambda * ((PoissonSelfStats*) sb)->nV;
}

bool PoissonParam::isZero() {
	return (this->lambda ==0);
}

void PoissonParam::calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	if (pab==NULL) {
		this->sumE = 0;
		this->sumlgam = 0;
	} else {
		this->sumE = ((PoissonPairStats*) pab)->sumE;
		this->sumlgam = ((PoissonPairStats*) pab)->sumlgam;
	}
	this->numT = ((PoissonSelfStats*) sa)->nV * ((PoissonSelfStats*) sb)->nV;
}

void PoissonParam::cleanup() {
	this->lambda = this->sumE / this->numT;
	if (std::isnan(lambda)) lambda=0;
}


class WParam : public ModelParamBase {
	public:
		float p;
		float numerator;
		float denominator;
		// some other stuff to go here may be
		virtual void calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual void collapse(ModelParamBase* pa, ModelParamBase* pb);
		virtual void cleanup();
		virtual void bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof);
		virtual bool isZero();
		virtual float predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual float updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar);
		virtual void init();
		virtual std::string DerivedType() { return std::string("WParam"); }
};

void WParam::init() {
	p = 0;
	numerator = 0;
	denominator = 0;
}

float WParam::updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar) {
	return 1 - (1 - oldsofar)*(1 - this->predict(sa,sb));
}

void WParam::bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof) {
	this->p = std::max(0.0f,(((WParam*) Ori)->p - ((WParam*) Sof)->p)/(1 - ((WParam*) Sof)->p));
}

void WParam::collapse(ModelParamBase* pa, ModelParamBase* pb) {
	this->numerator += ((WParam*) pa)->numerator + ((WParam*) pb)->numerator;
	this->denominator += ((WParam*) pa)->denominator + ((WParam*) pb)->denominator;
}

float WParam::predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	return ((WSelfStats*) sa)->degree * this->p * ((WSelfStats*) sb)->degree;
}


void WParam::cleanup() {
	this->p = this->numerator / this->denominator;
	if (std::isnan(p)) p=0;
}

bool WParam::isZero() {
	return (this->p ==0);
}

void WParam::calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	if (pab==NULL) this->numerator = 0;
	else this->numerator = ((WPairStats*) pab)->nE;
	this->denominator = ((WSelfStats*) sa)->degree * ((WSelfStats*) sb)->degree;
}

class DcorrParam: public ModelParamBase {
	public:
		float lambda;
		float numT;
		float nE;
		virtual void calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual void collapse(ModelParamBase* pa, ModelParamBase* pb);
		virtual void cleanup();
		virtual void bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof);
		virtual bool isZero();
		virtual float predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb);
		virtual float updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar);
		virtual void init();
		virtual std::string DerivedType() { return std::string("DcorrParam"); }
};

void DcorrParam::calculate(ModelPairStatsBase* pab, ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	if (pab==NULL) this->nE = 0;
	else this->nE = ((DcorrPairStats*) pab)->nE;
	this->numT = ((DcorrSelfStats*) sa)->degree * ((DcorrSelfStats*) sb)->degree;
}

void DcorrParam::collapse(ModelParamBase* pa, ModelParamBase* pb) {
	this->nE += ((DcorrParam*) pa)->nE + ((DcorrParam*) pb)->nE;
	this->numT += ((DcorrParam*) pa)->numT + ((DcorrParam*) pb)->numT;
}


void DcorrParam::cleanup() {
	this->lambda = this->nE / this->numT;
	if (std::isnan(lambda)) lambda=0;
}

void DcorrParam::init() {
	lambda = 0;
	nE = 0;
	numT = 0;
}

float DcorrParam::updatedSoFar(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb, float oldsofar) {
	return oldsofar + predict(sa,sb);
}

void DcorrParam::bestfromSoFar(ModelParamBase* Ori, ModelParamBase* Sof) {
	this->lambda = std::max(0.0f, ((DcorrParam*) Ori)->lambda - ((DcorrParam*) Sof)->lambda);
}

bool DcorrParam::isZero() {
	return (this->lambda ==0);
}


float DcorrParam::predict(ModelSelfStatsBase* sa, ModelSelfStatsBase* sb) {
	return ((DcorrSelfStats*) sa)->degree * this->lambda * ((DcorrSelfStats*) sb)->degree;
}



#endif

