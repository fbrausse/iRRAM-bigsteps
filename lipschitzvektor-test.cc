#include "iRRAM.h"
#include <vector>
#include <list>
#include "iRRAM/limit_templates.h"


using namespace iRRAM;
using std::vector;
using std::pair;

#include "extension-vector-2.h"


vector<REAL> my_code(const vector<REAL> & x){vector<REAL> y=x; 
  for( unsigned int i=0; i< y.size(); i++){ y[i]=sqrt(y[i]);}
   return y;
}

LAZY_BOOLEAN my_test(const vector<REAL> & x){return true;}

REAL my_lipf(const vector<REAL> & x){return pi()/100000000;}

void compute () {

FUNCTION<vector<REAL>,vector<REAL> > f = my_code;
FUNCTION<vector<REAL>,LAZY_BOOLEAN> d = my_test;
FUNCTION<vector<REAL>,REAL> L = my_lipf;

vector<REAL> x(2);
x[0]=sqrt(REAL(2));
x[1]=sqrt(REAL(3));
  { 
    sizetype lip_error;
    x[0].geterror(lip_error);
    fprintf(stderr,"argument x[0] on interval has error %d*2^(%d)\n",
              lip_error.mantissa,lip_error.exponent);
  }


vector<REAL> y=f(x);
  { 
    sizetype lip_error;
    y[0].geterror(lip_error);
    fprintf(stderr,"f_original on interval has error %d*2^(%d)\n",
              lip_error.mantissa,lip_error.exponent);
  }

for (int i=0; i<= 10;i++)
  x=lipschitz_maxnorm(f,L,d,x);
  
  { 
    sizetype lip_error;
    x[0].geterror(lip_error);
    fprintf(stderr,"result has error %d*2^(%d)\n",
              lip_error.mantissa,lip_error.exponent);
  }

  cout << x[0]<< " , " << x[1]<< "\n";

}


