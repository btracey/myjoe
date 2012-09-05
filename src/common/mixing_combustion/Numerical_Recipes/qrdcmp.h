#ifndef QRDCMP_H
#define QRDCMP_H

#include "nr3.h"

struct QRdcmp {
	Int n;
	MatDoub qt, r;
	Bool sing;
	QRdcmp(MatDoub_I &a);
	void solve(VecDoub_I &b, VecDoub_O &x);
	void qtmult(VecDoub_I &b, VecDoub_O &x);
	void rsolve(VecDoub_I &b, VecDoub_O &x);
	void update(VecDoub_I &u, VecDoub_I &v);
	void rotate(const Int i, const Doub a, const Doub b);
};
#endif

