#ifndef iRRAM_TaylorSeries_H
#define iRRAM_TaylorSeries_H

#include <iRRAM/lib.h>
#include <vector>
#include "TaylorModel.h"
#include "ComplexIntervals.h"


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
template <typename T>
REAL poly_bound(const FUNCTION<std::vector<T>,unsigned int> a, int max_index, const REAL &tdelta, int derive, bool include_a0) {
	c_int tc=c_int(INTERVAL(-tdelta,tdelta,true),INTERVAL(-tdelta,tdelta,true));
	
	REAL max_result(0);
	std::vector<T> a0= a(0);
	INTEGER dfactor=1;
	for (int i=2; i<=derive; i++) dfactor *=i;
	for (unsigned nu=0; nu < a0.size(); nu++) {
	      c_int res = REAL(0);
	      INTEGER factor=dfactor;
	      if (include_a0) res = REAL(factor)*static_cast<REAL>(a0[nu]);
	      for (int k=1; k<=max_index; k++){
//cerr << " next k "<<k<< "\n"; 
		std::vector<T> ak= a(k+derive);
//cerr << ak[nu]<<" k "<<k<< " maxindex "<< max_index<< "\n"; 
		REAL aknu(ak[nu]);
//cerr << "chk " <<aknu << " "<< k <<" "<< nu<< " "<< mag(power(tc,k))<<  " ###\n";
		factor=(factor*(k+derive))/k;
		res += REAL(factor)*static_cast<REAL>(ak[nu])*power(tc,k);
//cerr << res._Ireal.low<< res._Ireal.upp << " +i* "<< res._Iimag.upp<< "####\n"; 
	      }
	      REAL mx = mag(res);
	      max_result = maximum(mx, max_result);
	}
	return max_result;
}

//********************************************************************************
// Class for summation of a vector of Taylor series
// 'bound' must be valid for 'radius' in each dimension
//********************************************************************************

template <typename T>
class FUNCTIONAL_vector_taylor_sum : public FUNCTIONAL_object<std::vector<T>,REAL> {

	std::vector<FUNCTION<T,REAL>> f;

public:
	FUNCTIONAL_vector_taylor_sum(
		FUNCTION<std::vector<T>,unsigned int> coeff,
		const REAL& radius,
		const REAL& bound,
		unsigned bound_type=0
	)
	: f(coeff(0).size()) // just to get the dimension
	{
		for (unsigned nu=0; nu<f.size(); nu++)
			f[nu] = taylor_sum(projection(coeff,nu), radius, bound, bound_type);
	}

	std::vector<T> eval(const REAL& t)
	{
		std::vector<T> result(f.size());
		for (unsigned int nu=0; nu<f.size(); nu++)
			result[nu] = f[nu](t);
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
	return new FUNCTIONAL_vector_taylor_sum<TM>(coeff,radius,maximum,bound_type);
}


//********************************************************************************
//********************************************************************************

/****************************************************************************************/
/* summation of Taylor sequences							*/
/* The constructor needs a sequences of numbers with a radius of convergence 		*/
/* larger(!) than the given value 'radius'.						*/
/* 'bound' has to be an upper bound for the absolute value of the sum function 		*/
/* on a circle with the given 'radius' value						*/
/****************************************************************************************/


unsigned int max_taylor_coeff=0; // only statistical purpose!
unsigned int last_taylor_coeff=0; // only statistical purpose!

template <typename T>
class FUNCTIONAL_taylor_sum : public FUNCTIONAL_object<T,REAL> {

/***********************************************************/
/* local data: coefficients (as a function),  radius and corresponding bound */

	FUNCTION<T,unsigned int> _coeff;
	REAL _radius;
	REAL _bound;
	unsigned _bound_type;

public:
	FUNCTIONAL_taylor_sum(
			FUNCTION<T,unsigned int> coeff,
			const REAL& radius,
			const REAL& bound,
			unsigned bound_type
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

	T eval(const REAL &x)
	{
		single_valued code;
		T sum(0);
		T best(0);
		REAL factor=1;
		REAL error=_bound*_radius/(_radius-abs(x));
		REAL errorfactor=abs(x)/_radius;
		iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum starting with precision "<<actual_stack().actual_prec
			<< " at ratio "<< errorfactor.vsize.mantissa*pow(2,errorfactor.vsize.exponent)<<"\n";});

		sizetype sum_error,trunc_error,best_error;

		sizetype_add(best_error,error.vsize,error.error);
		int best_index=-1;

		for (unsigned int i=0;;i++){
			sum= sum + _coeff(i)*factor;
			error= error*errorfactor;

			sum_error = geterror(sum);
			sizetype_add(trunc_error,error.vsize,error.error);

			sizetype local_error;
			sizetype_add(local_error,trunc_error,sum_error);
			if (sizetype_less(local_error,best_error)) {
				best=sum;
				best_error=local_error;
				adderror(best, local_error);
				best_index=i;
			}
			if (trunc_error.exponent < actual_stack().actual_prec ||
			     sizetype_less(trunc_error,sum_error)) {
				iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum: stop at "
						<<i<< " with best at "<<best_index<<"\n";});
				if (i>max_taylor_coeff) max_taylor_coeff=i; // only statistical purpose!
				last_taylor_coeff = i;
				break;
			}
			factor=factor*x;
		}
		return best;
	}
};

#if 0 /* this specialization is somehow broken */
/****************************************************************************************/
/* summation of Taylor sequences							*/
/* The constructor needs a sequences of numbers with a radius of convergence 		*/
/* larger(!) than the given value 'radius'.						*/
/* 'bound_type' determines the meaning of 'bound':					*/
/* for 'bound_type=n', 'bound' has to be an upper bound 				*/
/* for the absolute value of the n-th derivative of the sum function 			*/
/* on a circle with the given 'radius' value 						*/
/****************************************************************************************/

template <>
class FUNCTIONAL_taylor_sum<TM> :public FUNCTIONAL_object<TM,REAL> {

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

		sizetype_add(best_error,error.vsize,error.error);
		int best_index=-1;

		for (unsigned int i=0;;i++){
			TM c = _coeff(i);
// 	cerr << c<<" c "<<i<<"\n";
// 	cerr << i<< " coeff: "<< c.to_real() <<"\n\n";
			sum= sum + c*factor;
			error= error*errorfactor;
			factor=factor*x;
//	cerr << "sum    "<< sum << "\n";
//	cerr << "factor "<< factor << "\n";
			  if ( i > _bound_type) {
			    REAL corrective = error;
			    for (unsigned j=1; j<=_bound_type;j++) corrective/=int(i+j);
			    sum.to_real().geterror(sum_error);
//			    sum.geterror_info(error_info);
			    sum.c0.geterror(error_info);
			    sizetype_add(trunc_error,corrective.vsize,corrective.error);

// test: how fast will things be is we assume the truncation error to be smaller than the current
// term?
//			(c*factor).to_real().getsize(trunc_error);
			
			    sizetype local_error;
			    sizetype_add(local_error,trunc_error,sum_error);
			    if (sizetype_less(local_error,best_error)) {
				best=sum;
				best_error=local_error;
				best.c0.adderror(trunc_error);
				best_index=i;
			    }
sizetype_shift(trunc_error,trunc_error,20);
			    if (// trunc_error.exponent<ACTUAL_STACK.actual_prec ||
			      sizetype_less(trunc_error,error_info)) {
				 iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum: stop at "
						<<i<< " with best at "<<best_index<<"\n";});
				if (i>iRRAM_DEBUG_max_taylor_coeff) iRRAM_DEBUG_max_taylor_coeff=i; // only statistical purpose!
				iRRAM_DEBUG_last_taylor_coeff = i;
//				best.show();
				break;
			    }
			}
		}
		return best;
	}
};
#endif

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<TM,REAL> taylor_sum (
		FUNCTION<TM,unsigned int> coeff,
		const REAL& radius,
		const REAL& bound,
		unsigned bound_type=0
) {
	return new FUNCTIONAL_taylor_sum<TM>(coeff,radius,bound,bound_type);
}

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<REAL,REAL> taylor_sum (
		FUNCTION<REAL,unsigned int> coeff,
		const REAL& radius,
		const REAL& bound,
		unsigned bound_type=0
) {
	return new FUNCTIONAL_taylor_sum<REAL>(coeff,radius,bound,bound_type);
}

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<std::vector<REAL>,REAL > taylor_sum (
		FUNCTION<std::vector<REAL>,unsigned int > coeff,
		const REAL& radius,
		const REAL& maximum,
		unsigned bound_type=0
) {
	return new FUNCTIONAL_vector_taylor_sum<REAL>(coeff,radius,maximum,bound_type);
}


} // namespace iRRAM

#endif /* !  iRRAM_TaylorSeries_H */
