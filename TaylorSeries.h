#ifndef iRRAM_TaylorSeries_H
#define iRRAM_TaylorSeries_H

#include "iRRAM.h"
#include <vector>
#include "TaylorModel.h"



namespace iRRAM {
  
// For testing and debugging purposes we currently provide two variables:
unsigned int iRRAM_DEBUG_max_taylor_coeff=0; 
unsigned int iRRAM_DEBUG_last_taylor_coeff=0; 
// Their use is only safe for testing and debugging purposes!
// No decisions are allowed to be based on these variables.
// They might vanish in further version without any notice!!!


// compute the Taylor polynomial up to \pm tdelta, but excluding the first coefficient 
//(i.e. compute the difference to the initial value)
//Then take the maximum for all dimensions 
// If derive=n tells that we want the n-th derivative of the polynomial with max_index coefficients
// If include_a0==false, then the first coefficient is excluded from the sum

// Perhaps Euclidean norm instead of maximum norm would be better????
REAL poly_bound(const FUNCTION<std::vector<TM>,unsigned int> a, int max_index, const REAL &tdelta, int derive, bool include_a0) {
	c_int tc=c_int(INTERVAL(-tdelta,tdelta,true),INTERVAL(-tdelta,tdelta,true));
	
	REAL max_result(0);
	std::vector<TM> a0= a(0);
	INTEGER dfactor=1;
	for (int i=2; i<=derive; i++) dfactor *=i;
	for (unsigned nu=0; nu < a0.size(); nu++) {
	      c_int res = REAL(0);
	      INTEGER factor=dfactor;
	      if (include_a0) res = REAL(factor)*(a0[nu].to_real());
	      for (int k=1; k<=max_index; k++){
//cerr << " next k "<<k<< "\n"; 
		std::vector<TM> ak= a(k+derive);
//cerr << ak[nu]<<" k "<<k<< " maxindex "<< max_index<< "\n"; 
		REAL aknu = ak[nu].to_real(); 
//cerr << "chk " <<aknu << " "<< k <<" "<< nu<< " "<< mag(power(tc,k))<<  " ###\n";
		factor=(factor*(k+derive))/k;
		res+= REAL(factor)*(ak[nu].to_real())*power(tc,k);
//cerr << res._Ireal.low<< res._Ireal.upp << " +i* "<< res._Iimag.upp<< "####\n"; 
	      }
	      REAL mx = mag(res);
	      max_result = maximum(mx, max_result);
	}
	return max_result;
}


/****************************************************************************************/
/* summation of Taylor sequences							*/
/* The constructor needs a sequences of numbers with a radius of convergence 		*/
/* larger(!) than the given value 'radius'.						*/
/* 'bound_type' determines the meaning of 'bound':					*/
/* for 'bound_type=n', 'bound' has to be an upper bound 				*/
/* for the absolute value of the n-th derivative of the sum function 			*/
/* on a circle with the given 'radius' value 						*/
/****************************************************************************************/

class FUNCTIONAL_taylor_sum :public FUNCTIONAL_object<TM,REAL> {

/***********************************************************/
/* local data: coefficients (as a function),  radius and corresponding bound */

	FUNCTION<TM,unsigned int> _coeff;
	REAL _radius;
	REAL _bound;
	unsigned _bound_type;
public:
	FUNCTIONAL_taylor_sum(
			FUNCTION<TM,unsigned int> coeff,
			const REAL& radius,
			const REAL& bound,
			const unsigned bound_type
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
		_bound_type=bound_type;
	}

	~FUNCTIONAL_taylor_sum() {}

/****************************************************************************************/
/* Evaluation:										*/
/* Compute coefficients until the truncation error term is smaller than actual_prec	*/
/* or until the accumulated error is already of the order of the truncation error	*/
/* (whichever comes first)								*/
/****************************************************************************************/

	TM eval(const REAL &x)
	{
		single_valued code;

		TM sum=REAL(0);
		TM best=REAL(0);
		REAL factor=1;
		REAL error=_bound*power(_radius,_bound_type+1)/(_radius-abs(x));
		REAL errorfactor=abs(x)/_radius;
		iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum starting with precision "<<actual_stack().actual_prec 
		  << " at ratio "<< errorfactor.vsize.mantissa*pow(2,errorfactor.vsize.exponent)<<"\n";});

		sizetype sum_error,trunc_error,best_error,error_info;

		sizetype a,b;
		error.geterror(a);
		error.getsize(b);
		sizetype_add(best_error,a,b);
		int best_index=-1;

		for (unsigned int i=0;;i++){
			TM c = _coeff(i);
// 	cerr << c<<" c "<<i<<"\n";
// 	cerr << i<< " coeff: "<< c.to_real() <<"\n\n";
			sum= sum + c*factor;
// 	cerr << sum<<" sum "<<i<<"\n";
			factor=factor*x;

			error= error*errorfactor;

			  if ( i > _bound_type) {
			    REAL corrective = error;
			    for (unsigned j=1; j<=_bound_type;j++) corrective/=int(i+j);
			    sum.to_real().geterror(sum_error);
			    sum.c0.geterror(error_info);
			    sizetype a,b;
			    corrective.geterror(a);
			    corrective.getsize(b);
			    sizetype_add(trunc_error,a,b);

// test: how fast will things be is we assume the truncation error to be smaller than the current
// term?
//			(c*factor).to_real().getsize(trunc_error);
			
			    sizetype local_error;
			    sizetype_add(local_error,trunc_error,sum_error);
//			    if (sizetype_less(local_error,best_error)) {
			    if ( true ) {
				best=sum;
				best_error=local_error;
				best.c0.adderror(trunc_error);
				best_index=i;
			    }
			    sizetype_set(error_info,1,error_info.exponent);
			    if ( sizetype_less(trunc_error,error_info) 
//			      || trunc_error.exponent <= sum.c0.vsize.exponent-80+ACTUAL_STACK.actual_prec
			    ) {
				 iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum: stop at "
						<<i<< " with best at "<<best_index<<"\n";});
				if (i>iRRAM_DEBUG_max_taylor_coeff) iRRAM_DEBUG_max_taylor_coeff=i; // only statistical purpose!
				iRRAM_DEBUG_last_taylor_coeff = i;
//				cerr << sum<<" summe "<< i<< "\n";
//				cerr << best<<" best "<< best_index<< "\n";
				break;
			    }
			}
		}
		return best;
	}
};

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<TM,REAL> taylor_sum (
		FUNCTION<TM,unsigned int> coeff,
		const REAL& radius,
		const REAL& bound,
		const unsigned bound_type=0
) {
	return new FUNCTIONAL_taylor_sum(coeff,radius,bound,bound_type);
}

//********************************************************************************
// Class for summation of a vector of Taylor series
// 'bound' must be valid for 'radius' in each dimension
//********************************************************************************

class FUNCTIONAL_vector_taylor_sum :public FUNCTIONAL_object<std::vector<TM>,REAL > {

	REAL _radius;
	REAL _bound;
	unsigned _bound_type;
	FUNCTION<std::vector<TM>,unsigned int > _coeff; 
	std::vector<FUNCTION<TM,REAL >> _f_nu;

public:
	FUNCTIONAL_vector_taylor_sum(
			FUNCTION<std::vector<TM>,unsigned int > coeff,
			const REAL& radius,
			const REAL& bound,
			const int bound_type
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
		_bound_type=bound_type;
	}

	~FUNCTIONAL_vector_taylor_sum() {}

	std::vector<TM> eval(const REAL& t)
	{
		if (_f_nu.size()==0) { // on first call, initialize the projections
		  std::vector<TM> test=_coeff(0); // just to get the dimension
		  _f_nu.resize(test.size());
		  for (unsigned int nu=0; nu<_f_nu.size(); nu++){
			FUNCTION<TM,unsigned int>  coeff_nu=projection(_coeff,nu);
			_f_nu[nu]=taylor_sum(coeff_nu,_radius,_bound,_bound_type);
		  }
		}  
		std::vector<TM> result( _f_nu.size());
		for (unsigned int nu=0; nu<_f_nu.size(); nu++){
			result[nu]=_f_nu[nu](t);
		}
		return result;
	}
};

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<std::vector<TM>,REAL > taylor_sum (
		FUNCTION<std::vector<TM>,unsigned int > coeff,
		const REAL& radius,
		const REAL& maximum,
		const unsigned bound_type=0
) {
	return new FUNCTIONAL_vector_taylor_sum(coeff,radius,maximum,bound_type);
}


} // namespace iRRAM

#endif /* !  iRRAM_TaylorSeries_H */
