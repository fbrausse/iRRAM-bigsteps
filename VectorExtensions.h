#include "iRRAM/core.h"

namespace iRRAM {

template <typename... P>
FUNCTION<std::vector<TM>,P...> lipschitzify(
	const FUNCTION<FUNCTION<std::vector<TM>,P...>,std::vector<TM>> &f,
	const FUNCTION<REAL,std::vector<TM>,P...> &lip_f,
	const FUNCTION<LAZY_BOOLEAN,std::vector<TM>> &on_domain,
	const std::vector<TM> &x
) {
	single_valued code_1; 

	if (!on_domain(x))
		iRRAM_REITERATE(0);

	std::vector<TM> x_new(x.size());
	sizetype arg_error;
	{
		sizetype_exact(arg_error);

		for (unsigned int i=0;i<x.size();i++) {
			sizetype x_error;
			DYADIC x_center;
			x[i].to_real().to_formal_ball(x_center,x_error);
			x_new[i]=TM(REAL(x_center));
			sizetype_max(arg_error,x_error,arg_error);
		}
	}

	cerr << "starting function lipschitz_maxnorm\n";

	FUNCTION<std::vector<TM>,P...> fl;
	{
		fl = f(x_new);
	}

	return from_algorithm(
		std::function<std::vector<TM>(const P &...)>(
			[fl,lip_f,arg_error,x](const P &... p)
			{
				single_valued code_1; 

				REAL lip_bound = lip_f(x,p...);
				cerr << "# lip_bound = "<< lip_bound<< "\n";

				std::vector<TM> lip_result;
				{
					lip_result = fl(p...);
				}

				sizetype lip_error,lip_size;
				lip_bound.getsize(lip_size);

				sizetype_mult(lip_error,lip_size,arg_error);

				for (unsigned int i=0;i<lip_result.size();i++) {
					lip_result[i].c0.adderror(lip_error);
				}

				return lip_result;
			}
		)
	);
}

std::vector<REAL> lipschitz_maxnorm (
    const FUNCTION<std::vector<REAL>,std::vector<REAL> >& f,
    const FUNCTION<REAL,std::vector<REAL>>& lip_f,
    const FUNCTION<LAZY_BOOLEAN,std::vector<REAL>>& on_domain,
    const std::vector<REAL>& x)
/* 
 * lipschitz_maxnorm(f,lip_f,on_domain,x) computes f(x) with an optimized algorithm,
 * 	with reduced error propagation
 * lip_f: Lipschitz-bound valid for f whereever on_domain(x)=true,
 * 	based on the maximum norm
 * on_domain: total, multi-valued with property: on_domain(x)=TRUE => x \in domain(f)
*/
{
  single_valued code_1; 

  if (!on_domain(x))
    iRRAM_REITERATE(0);
  REAL lip_bound=lip_f(x);
  cerr << "# lip_bound = "<< lip_bound<< "\n";

//  stiff code_2;
  
  std::vector<REAL> x_new(x.size());
  sizetype arg_error;
  sizetype_exact(arg_error);
  
  {
    sizetype x_error;
    DYADIC x_center;
    for (unsigned int i=0;i<x.size();i++){ 
      x[i].to_formal_ball(x_center,x_error);
      x_new[i]=x_center;
      sizetype_max(arg_error,x_error,arg_error);
     }
  }
  
  cerr << "starting function lipschitz_maxnorm\n";
  
  std::vector<REAL> lip_result = f(x_new);

  sizetype lip_error,lip_size;
  lip_bound.getsize(lip_size);

  sizetype_mult(lip_error,lip_size,arg_error);

  
  for (unsigned int i=0;i<lip_result.size();i++){ 
        lip_result[i].adderror(lip_error);
  }
//   { 
//     lip_result.geterror(lip_error);
//     fprintf(stderr,"end of lipschitz_new with error %d*2^(%d)\n",
//               lip_error.mantissa,lip_error.exponent);
//     fprintf(stderr,"  for argument with error %d*2^(%d)\n",
//               x_error.mantissa,x_error.exponent);
//   }

  return lip_result;
}


std::vector<REAL> lipschitz_maxnorm (
    const FUNCTION<std::vector<REAL>,std::vector<REAL> >& f,
    const REAL& lip_c,
    const FUNCTION<LAZY_BOOLEAN,std::vector<REAL>>& on_domain,
    const std::vector<REAL>& x){
  return lipschitz_maxnorm(f,from_value<REAL,std::vector<REAL>>(lip_c),on_domain,x);
  //return lipschitzify(f,from_value<REAL,std::vector<REAL>>(lip_c),on_domain,x)();
}

// FUNCTION<vector<REAL>, vector<REAL> > lipschitz_maxnorm (
//     const FUNCTION<vector<REAL>,vector<REAL> >& f,
//     const FUNCTION<vector<REAL>,REAL>& lip_f,
//     const FUNCTION<vector<REAL>,LAZY_BOOLEAN>& on_domain){
//   
// FUNCTION<vector<REAL>, vector<REAL> > g;   
//      return g;
// }
  
} /* ! namespace iRRAM */

