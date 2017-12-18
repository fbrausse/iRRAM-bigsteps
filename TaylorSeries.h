#ifndef iRRAM_TaylorSeries_H
#define iRRAM_TaylorSeries_H

#include "iRRAM.h"
#include <vector>
#include "TaylorModel.h"



namespace iRRAM {
  
/****************************************************************************************/
/* summation of Taylor sequences							*/
/* The constructor needs a sequences of numbers with a radius of convergence 		*/
/* larger(!) than the given value 'radius'.						*/
/* 'bound' has to be an upper bound for the absolute value of the sum function 		*/
/* on a circle with the given 'radius' value						*/
/****************************************************************************************/

class FUNCTIONAL_taylor_sum :public FUNCTIONAL_object<REAL,REAL> {

/***********************************************************/
/* local data: coefficients (as a function),  radius and corresponding bound */

	FUNCTION<REAL,unsigned int> _coeff;
	REAL _radius;
	REAL _bound;

public:
	FUNCTIONAL_taylor_sum(
			FUNCTION<REAL,unsigned int> coeff,
			const REAL& radius,
			const REAL& bound
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
	}

	~FUNCTIONAL_taylor_sum() {}

/****************************************************************************************/
/* Evaluation:										*/
/* Compute coefficients until the truncation error term is smaller than actual_prec	*/
/* or until the accumulated error is already of the order of the truncation error	*/
/* (whichever comes first)								*/
/****************************************************************************************/

	REAL eval(const REAL &x)
	{
		single_valued code;
		REAL sum=0;
		REAL best=0;
		REAL factor=1;
		REAL error=_bound*_radius/(_radius-abs(x));
		REAL errorfactor=abs(x)/_radius;
		iRRAM_DEBUG0(2,{cerr << "FUNCTIONAL_taylor_sum starting with precision "<<actual_stack().actual_prec
			<< " at ratio "<< errorfactor.vsize.mantissa*pow(2,errorfactor.vsize.exponent)<<"\n";});

		sizetype best_error;

		sizetype_add(best_error,error.vsize,error.error);
		int best_index=-1;

		for (unsigned int i=0;;i++){
			sizetype trunc_error;

			sum= sum + _coeff(i)*factor;
			error= error*errorfactor;

			sizetype sum_error = geterror(sum);
			sizetype_add(trunc_error,error.vsize,error.error);

			sizetype local_error;
			sizetype_add(local_error,trunc_error,sum_error);
			if (sizetype_less(local_error,best_error)) {
				best=sum;
				best_error=local_error;
				best.seterror(best_error);
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

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<REAL,REAL> taylor_sum (
		FUNCTION<REAL,unsigned int> coeff,
		const REAL& radius,
		const REAL& bound
) {
	return new FUNCTIONAL_taylor_sum(coeff,radius,bound);
}

//********************************************************************************
// Class for summation of a vector of Taylor series
// 'bound' must be valid for 'radius' in each dimension
//********************************************************************************

class FUNCTIONAL_vector_taylor_sum :public FUNCTIONAL_object<std::vector<REAL>,REAL > {

	REAL _radius;
	REAL _bound;
	FUNCTION<std::vector<REAL>,unsigned int > _coeff;

public:
	FUNCTIONAL_vector_taylor_sum(
			FUNCTION<std::vector<REAL>,unsigned int > coeff,
			const REAL& radius,
			const REAL& bound
	) {
		_coeff=coeff;
		_radius=radius;
		_bound=bound;
	}

	~FUNCTIONAL_vector_taylor_sum() {}

	std::vector<REAL> eval(const REAL& t)
	{
		std::vector<REAL> result=_coeff(0); // just to get the dimension
		for (unsigned int nu=0;nu<result.size();nu++){
			FUNCTION<REAL,unsigned int>  coeff_nu=projection(_coeff,nu);
			FUNCTION<REAL,REAL> f_nu=taylor_sum(coeff_nu,_radius,_bound);
			result[nu]=f_nu(t);
		}
		return result;
	}
};

/********************************************************/
/* corresponding constructor for FUNCTION		*/
/********************************************************/
inline FUNCTION<std::vector<REAL>,REAL > taylor_sum (
		FUNCTION<std::vector<REAL>,unsigned int > coeff,
		const REAL& radius,
		const REAL& maximum
) {
	return new FUNCTIONAL_vector_taylor_sum(coeff,radius,maximum);
}


} // namespace iRRAM

#endif /* !  iRRAM_TaylorSeries_H */
