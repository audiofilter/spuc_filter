/*
    Copyright (C) 2014 Tony Kirke

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <spuce/typedefs.h>
#include "other_freq.h"

namespace spuce {

void cic_freq(int rate, int order, int pts, double* w, double inc) {
	double db=0;
	double sum = 0;
	double wf = inc*M_PI/(double)pts;		
	for (int i=0;i<pts;i++) {
		if (i!=0) sum = (1.0/rate)*sin(0.5*wf*i*rate)/sin(0.5*wf*i);
    else sum = 0;
		db = 10.0*order*log(sum*sum)/log(10.0);
		w[i] = db;
	}
}
}
	
