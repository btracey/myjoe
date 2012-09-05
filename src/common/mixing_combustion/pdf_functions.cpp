#include "pdf_functions.h"
double betai(double a, double b, double x) {
  double bt,betai;
  if ((x<0.0)||(x>1.0)) {cout << "bad argument X in BETAI. X = " << x << endl;
                         throw(-1);}

  if ((x == 0.0)||(x == 1.0)) { bt=0.0; }
  else { bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x)); }

  if (x < (a+1.0)/(a+b+2.0)) { betai = bt*betacf(a,b,x)/a; }
  else { betai = 1.0 - bt*betacf(b,a,1.0-x)/b; }

  return betai;
}

double betacf(double a, double b, double x) {
  double eps=3.0e-7,am,bm,az,qab,qap,qam,bz,aold,em,tem,d,ap,bp,app,bpp;
  int m,itmax=10000;

  am = 1.0;
  bm = 1.0;
  az = 1.0;
  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  bz = 1.0 - qab*x/qap;

  m=0;
  aold=0.0;
  while (m <= itmax)
  {
    m = m + 1;
    em = (double) m;
    tem = em + em;
    d = em*(b-(double) m)*x/((qam+tem)*(a+tem));
    ap = az + d*am;
    bp = bz + d*bm;

    d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem));
    app = ap + d*az;
    bpp = bp + d*bz;

    aold = az;
    am = ap/bpp;
    bm = bp/bpp;
    az = app/bpp;
    bz = 1.0;

    if (fabs(az-aold)<eps*fabs(az)) break;
  }

  if (m>itmax) { cout << "A or B too big, or ITMAX too small" << endl;
                 cout << "A =" << a <<endl;
                 cout << "B =" << b <<endl;
                 cout << "Number of iterations: "<< m <<endl;
                 cout << "Error: "<< fabs(az-aold)/fabs(az) << " vs tolerance: "<< eps<<endl;
                 throw(-1);}
  return az;
}

  double beta(double a, double b) {
    return  exp(lgamma(a)+lgamma(b)-lgamma(a+b));
  }

// Compute \int_0^x t^n beta(t,a,b) dt
double betai_monomial( double a, double b, double x, int n) {
  double result;
  double n_dble = (double) n;
  double I0 = betai(a+n_dble, b, x);
  if (n == 0) {
    result = I0;
  } else {
    //double factor = beta(a+n_dble,b)/beta(a,b);
    // This is bad behaved for a or b big
    // Instead use: beta(a+n_dble,b)/beta(a+n_dble-1.0,b) = (a+n_dble-1.0)/(a+b+n_dble-1.0)
    double factor = 1.0;
    double aux = n_dble;
    while (aux >= 1.0) {
      factor *= (a+aux-1.0)/(a+b+aux-1.0);
      aux -= 1.0;
    }
    result = factor*I0;
    if (result != result) {
      cout << "Nan in betai_monomial " << a <<" "<<b<<" "<<x<<" "<<n<<" "<<n_dble<<endl;
      cout << "  aux: "<< factor << " " << I0 <<endl;
      throw(-1);
    }
  }

  return result;
}

