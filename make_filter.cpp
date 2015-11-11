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
#include "make_filter.h"
#include "other_freq.h"
#include <spuce/filters/iir.h>
#include <spuce/filters/fir_coeff.h>
#include <spuce/filters/notch_allpass.h>
#include <spuce/filters/iir_coeff.h>
#include <spuce/filters/butterworth_iir.h>
#include <spuce/filters/chebyshev_iir.h>
#include <spuce/filters/elliptic_iir.h>
#include <spuce/filters/butterworth_fir.h>
#include <spuce/filters/gaussian_fir.h>
#include <spuce/filters/create_remez_lpfir.h>
#include <spuce/filters/raised_cosine.h>
#include <spuce/filters/root_raised_cosine.h>
#include <spuce/filters/shelf_allpass1.h>
namespace spuce {

  void make_filter::sel_filter(const char* c_sel) {
    std::string sel(c_sel);
    if (sel == "Chebyshev") change_filter(Chebyshev);
    else if (sel=="Maxflat Subband") change_filter(MaxflatHalfband);
    else if (sel=="Elliptic") change_filter(Elliptic);
    else if (sel=="Elliptic Subband") change_filter(EllipticHalfband);
    else if (sel=="Butterworth")	change_filter(Butterworth);
    else if (sel=="Maxflat FIR") change_filter(MaxflatFIR);
    else if (sel=="Gaussian FIR") change_filter(GaussianFIR);
    else if (sel=="Remez FIR") change_filter(RemezFIR);
    else if (sel=="Raised Cosine") change_filter(RaisedCosine);
    else if (sel=="Notch") change_filter(NotchIIR);
    else if (sel=="Cut/Boost") change_filter(CutBoost);
    else if (sel=="Shelf") change_filter(Shelf);
    else if (sel=="Root Raised Cosine") change_filter(RootRaisedCosine);
    else if (sel=="None") change_filter(None);
    else std::cout << "Invalid filter selection\n";
  }

  void make_filter::change_filter(fil_enum f) {
    last_shape = shape;
    shape = f;
  }

  void make_filter::set_filter_type(int t) {
    if (t==1) pass_type = high;
    else if (t==2) pass_type = band;
    else pass_type = low;
  }

  void make_filter::init(int points) { pts = points; }
  void make_filter::set_fs(double f) { fs = f*1000.0; }

  make_filter::~make_filter() {}
  make_filter::make_filter() : B_IIR(4), C_IIR(4), E_IIR(4),
                               CUT_B(0.125,0.9,0.001), 
                               B_Sub(0,3,2), E_Sub(0.3,3,2), 
                               S_IIR(0.1,0.4) {
    reset();
    fs = 44100;
    phase = 0;
  }


  void make_filter::clear_filters() {
    M_Fir.set_size(0);	
    R_Fir.set_size(0); 
    G_Fir.set_size(0);	
    RC_Fir.set_size(0);	
    RRC_Fir.set_size(0); 	
    E_IIR.clear(); 
    C_IIR.clear(); 
    B_IIR.clear(); 
  }

  void make_filter::reset() {
    nested_k = 0;
    hpf = false;
    pass_type = low;

    elliptic_fc = butterworth_fc = chebyshev_fc = 0.125;
    maxflat_fc = 0.16;
    gauss_fc = 0.06;
    notch_fc = 0.25;
    cut_fc = 0.125;
    rc_fc = rrc_fc = 0.5;
	
    remez_pass_edge = 0.4;
    remez_trans = 0.2;
    remez_stop_edge = remez_pass_edge+remez_trans;
    remez_stop_weight = 50;
	
    elliptic_pass_edge = 0.4;
    elliptic_trans = 0.2;
    elliptic_stop_edge = elliptic_pass_edge + elliptic_trans;
	
    elliptic_stop_db = 50;
    elliptic_ripple = 1.0;
	
    chebyshev_ripple = 1.0;
	
    rc_alpha = rrc_alpha = 0.25;
    notch_trans = cut_trans = notch50_trans = 0.9;
	
    elliptic_order =  4;
    butterworth_order = 4;
    chebyshev_order = 4;
    gauss_taps = 21;
    remez_taps = 33;
    maxflat_taps = 45;
    rc_taps = rrc_taps = 33;
	
    elliptic_halfband_order = maxflat_halfband_order = 3;
	
    elliptic_halfband_ripple = 0.3;
    variable_ripple = 0.3;
    elliptic_halfband_rate  = 2;
    maxflat_halfband_rate = 2;
    shape = Chebyshev;
    last_shape = shape;
    bits = 0;
    shelf_low = 9;
    shelf_high = 1;
    shelf_fc = 0.25;
    shelf_gain = 10.0; 
    low_shelf_gain = 10.0;
    high_shelf_gain = -10.0;
    phase = 0;
  }

  double make_filter::limit(double x, double mx, double mn) {
    if (x>mx) x = mx;
    else if (x<mn) x = mn;
    return x;
  }

  double make_filter::get_fc(int len,bool in_passband) {
    // Convert swipe to dB inc/decrease
    double fc=0.5;
    double gain = 	pow(2,0.002*len);

    //	std::cout << "(len) " << len << " gain = " << gain << "\n";
	
    switch (shape) {
		case EllipticHalfband: 
			fc = 1.0/elliptic_halfband_rate; 
			break;
		case Butterworth:      fc = limit(gain*butterworth_fc,0.95,0.001); 
			break;
		case Chebyshev:        fc = limit(gain*chebyshev_fc,0.95,0.001); 
			break;
		case Elliptic:         
			if (in_passband) {
				fc = limit(gain*elliptic_pass_edge,0.95-elliptic_trans,0.001); 
			} else {
				fc = elliptic_pass_edge;
			}
			break;
		case MaxflatHalfband: 
			fc = 1.0/maxflat_halfband_rate; 
			break;			
			
		case MaxflatFIR:  fc = 2.0*limit(gain*maxflat_fc,0.4,0.001);
			break;
		case GaussianFIR:  fc = 2.0*limit(gain*gauss_fc,0.4,0.001); 
			break;
		case RemezFIR: 
			if (in_passband) {
				fc = limit(gain*remez_pass_edge,1.0-remez_trans,0.001); 
			} else {
				fc = remez_pass_edge;
			}
			break;
		case NotchIIR: fc = 2.0*limit(gain*notch_fc,0.5,0.001); 
			break;
		case CutBoost: fc = 2.0*limit(gain*cut_fc,0.5,0.001);  
			break;
		case RaisedCosine:  
			if (in_passband) {
				fc = limit(gain*rc_fc,0.5,0.001); 
			} else {
				fc = rc_fc;
			}
			break;
		case RootRaisedCosine:
			break;
    case None:
      fc = 0;
      break;
    case Shelf:
      fc = 0;
      break;
    }
    return(fc);
  }

  int make_filter::get_order() {
    switch (shape) {
    case EllipticHalfband: return(elliptic_halfband_order); break;
    case Butterworth:      return(butterworth_order); break; 
    case Chebyshev:        return(chebyshev_order); break;
    case Elliptic:         return(elliptic_order); break;
    case MaxflatHalfband:  return(maxflat_halfband_order); break;
    case MaxflatFIR:       return(maxflat_taps); break;
    case GaussianFIR:      return(gauss_taps); break; 
    case RemezFIR:         return(remez_taps); break;
    case NotchIIR: return(2);
    case CutBoost: return(2);
    case RaisedCosine:     return(rc_taps); break; 
    case RootRaisedCosine: return(rrc_taps); break;
    case Shelf: return(0.0); break;
    case None: return(0);
    }
    return(0);
  }
  bool make_filter::is_bpf() {
    return(pass_type == band);
  }

  bool make_filter::is_fir() {
    switch (shape) {
		case MaxflatFIR:       return(true); break;
		case GaussianFIR:      return(true); break; 
		case RemezFIR:         return(true); break;
		case RaisedCosine:     return(true); break; 
		case RootRaisedCosine: return(true); break;
    default: break;
    }
    return(false);
  }
  double make_filter::ripple() {
    switch (shape) {
		case RemezFIR:         return(0.0); break; // for now
		case Elliptic: return(elliptic_ripple); break;
    case Chebyshev: return(chebyshev_ripple); break;
    default: return(0.0);
    }
    return(0.0);
  }

  double make_filter::stopdB() {
    if (shape == Elliptic) return(elliptic_stop_db); 
    else return(0.0);
  }

  double make_filter::change_center(double f) {
    return(0);
  }

  double make_filter::horiz_swipe(int len,bool in_passband) {
    // Convert swipe to dB inc/decrease
    const double min_fc = 0.0;
    
    double gain = 	pow(2,0.002*len);
    double inc;
	
    if (len < 0) inc = 2;
    else inc = 0.5;
    
    switch (shape) {
		case EllipticHalfband: 
			if (in_passband) {
				elliptic_halfband_rate = limit(inc*elliptic_halfband_rate,64,2);
        //				if ((elliptic_halfband_rate == 2) && (len > 0)) elliptic_halfband_hpf = true;
			} else {
				elliptic_halfband_ripple = limit(elliptic_halfband_ripple/gain,0.5,0.001);
			}
			break;
		case Butterworth:      butterworth_fc = limit(gain*butterworth_fc,0.95,min_fc); 
			break;
		case Chebyshev:        chebyshev_fc = limit(gain*chebyshev_fc,0.95,min_fc); 
			break;
		case Elliptic:         
			if (in_passband) {
				elliptic_pass_edge = limit(gain*elliptic_pass_edge,0.95-elliptic_trans,0.001); 
			} else {
        elliptic_stop_db = limit(gain*elliptic_stop_db,100,10.0);	
        elliptic_trans = limit(gain*elliptic_trans,0.95-elliptic_pass_edge,0.001); 
			}
			elliptic_stop_edge = elliptic_pass_edge+elliptic_trans;
			break;
		case MaxflatHalfband: 
			maxflat_halfband_rate = limit(inc*maxflat_halfband_rate,64,2);
			break;
		case MaxflatFIR:  maxflat_fc = limit(gain*maxflat_fc,0.4,0.001);
			break;
		case GaussianFIR:  gauss_fc = limit(gain*gauss_fc,0.4,min_fc); 
			break;
		case RemezFIR: 
			if (in_passband) {
				remez_pass_edge = limit(gain*remez_pass_edge,0.95-remez_trans,0.001); 
			} else {
				remez_trans = limit(gain*remez_trans,0.95-remez_pass_edge,0.001); 
			}
			remez_stop_edge = remez_pass_edge+remez_trans;
			break;
		case NotchIIR: notch_fc = limit(gain*notch_fc,0.5,min_fc); 
			break;
		case CutBoost: cut_fc = limit(gain*cut_fc,0.5,min_fc);  
			CUT_B.set_freq(cut_fc); 
			break;
		case RaisedCosine: 
			if (in_passband) {
				rc_fc = limit(gain*rc_fc,0.5,0.001);
			} else {
				rc_alpha = limit(gain*rc_alpha,1,0.01); 
			}
			break;
		case RootRaisedCosine: rrc_alpha = limit(gain*rrc_alpha,1,0.01); 
			break;
		case Shelf: shelf_fc	= limit(gain*shelf_fc,0.5,0.0); break;
    case None: break;
    }
    return(0.0);
	
  }
  void make_filter::vertical_swipe(int len, bool in_passband, bool above_stop) {
    const int MAX_IIR = 20;
    const int MIN_IIR = 1;
    const int MAX_FIR = 99;
    const int MIN_FIR = 15;
    int inc;
	
    if (len < 0) inc = 1;
    else inc = -1;
	
    // Convert swipe to dB inc/decrease
    double gain = 	pow(2,0.002*len);
    double ogain = 1.0/gain;
	
	
    switch (shape) {
		case Butterworth: butterworth_order = limit(butterworth_order+inc,MAX_IIR,MIN_IIR); break;
		case Chebyshev: 
			if (in_passband) {
				chebyshev_ripple = limit(ogain*chebyshev_ripple,20,0.0001);
			} else {
				chebyshev_order = limit(chebyshev_order+inc,MAX_IIR,MIN_IIR); 
			}
			break;
		case Elliptic: 
			if (in_passband) {
				elliptic_ripple = limit(ogain*elliptic_ripple,10,0.0001);
			} else {
				if (above_stop) {
					elliptic_stop_db = limit(ogain*elliptic_stop_db,100,10.0);					
				} else{
					elliptic_order = limit(elliptic_order+inc,MAX_IIR,MIN_IIR); 
				}
			}
			break;
		case EllipticHalfband: elliptic_halfband_order = limit(elliptic_halfband_order+inc,MAX_IIR,MIN_IIR); break;
			
		case MaxflatHalfband: maxflat_taps = limit(maxflat_taps+inc,MAX_FIR,MIN_FIR);break;
		case MaxflatFIR: maxflat_taps = limit(maxflat_taps+8*inc,MAX_FIR,MIN_FIR);break;
		case GaussianFIR: gauss_taps = limit(gauss_taps+inc,MAX_FIR,MIN_FIR);break;
			// FIX ME - should change to passband ripple while in passband
		case RemezFIR:	  
			if (in_passband) {
				remez_stop_weight = limit(ogain*remez_stop_weight,100,0.01);
			} else {
				if (above_stop) {
					remez_stop_weight = limit(ogain*remez_stop_weight,100,0.01);
				} else{
					remez_taps = limit(remez_taps+inc,MAX_FIR,MIN_FIR);break;
				}
			}			
		case RootRaisedCosine: rrc_taps = limit(rrc_taps+2*inc,MAX_FIR,MIN_FIR);break;
		case RaisedCosine: rc_taps = limit(rc_taps+2*inc,MAX_FIR,MIN_FIR);break;
			
		case NotchIIR: notch_trans = limit(gain*notch_trans,1,0.001); break;
		case CutBoost: cut_trans = limit(gain*cut_trans,10,0.0); CUT_B.set_depth(cut_trans); break;
		case Shelf: shelf_gain = limit(0.01*len+shelf_gain,100,-100);	break;
    default: break;
			
    }
  }
  double make_filter::update(double *w) {return(update(w,1.0)); }
  double make_filter::update(double *w, double inc) {
    double fc;
    double e_rate = elliptic_halfband_rate;
    double m_rate = maxflat_halfband_rate;

    switch (shape) {
    case None:
      for (int i=0;i<pts;i++) w[i] = 1.0;
      break;
		case MaxflatFIR:
			Maxflat_Fir.set_size(maxflat_taps);
      fc = maxflat_fc;
			butterworth_fir(Maxflat_Fir, fc);
			M_Fir.set_size(maxflat_taps);
			M_Fir.settaps(Maxflat_Fir);
			fir_freq(M_Fir,pts,w,inc);
      return(maxflat_fc);
			break;
		case RemezFIR:
      {
        fir_coeff<double> Remez_Fir(remez_taps);
        double pass_edge = remez_pass_edge;
        double stop_edge = remez_stop_edge;
        create_remez_lpfir(Remez_Fir,pass_edge,stop_edge,remez_stop_weight);
        R_Fir.set_size(remez_taps);
        R_Fir.settaps(Remez_Fir);
        fir_freq(R_Fir,pts,w,inc);
			
      }
      return(1.0-remez_pass_edge);
			break;
		case GaussianFIR:
      {
        fir_coeff<double> Gaussian_Fir(gauss_taps);
        fc = gauss_fc;
        gaussian_fir(Gaussian_Fir, fc);
        G_Fir.set_size(gauss_taps);
        G_Fir.settaps(Gaussian_Fir);
        fir_freq(G_Fir,pts,w,inc);

      }
      if (pass_type==high) return(1.0-gauss_fc);
      else return(gauss_fc);
      break;
		case RaisedCosine:
      {
        fir_coeff<double> RaisedCosine_Fir(rc_taps);
        fc = 1.0/rc_fc;
        raised_cosine(RaisedCosine_Fir, rc_alpha, fc);
        RC_Fir.set_size(rc_taps);
        RC_Fir.settaps(RaisedCosine_Fir);
        fir_freq(RC_Fir,pts,w,inc);
      }
			//std::cout << "RC alpha = " << rc_alpha << " taps = " << rc_taps << " ";
			//std::cout << " sum = " << RC_Fir.coeff_sum() << "\n";
      return(rc_fc);
			break;
		case RootRaisedCosine:
      {
        fir_coeff<double> RootRaisedCosine_Fir(rrc_taps);
        root_raised_cosine(RootRaisedCosine_Fir, rrc_alpha, 2);
        RRC_Fir.set_size(rrc_taps);
        RRC_Fir.settaps(RootRaisedCosine_Fir);
        fir_freq(RRC_Fir,pts,w,inc);
      }
			return(0.5);
			break;
			
			
      // Special Allpass based IIRs
    case MaxflatHalfband:
			B_Sub.reset();
			B_Sub.set_coeffs(0,maxflat_halfband_order,m_rate);
			other_freq(B_Sub,pts,w,inc);
			return(1.0/maxflat_halfband_rate);
			break;

		case EllipticHalfband:
			E_Sub.reset();			
			E_Sub.set_coeffs(elliptic_halfband_ripple,elliptic_halfband_order,
                       e_rate);
			other_freq(E_Sub,pts,w,inc);
			return(1.0/elliptic_halfband_order);
			break;
			
			
			////////// IIRs
			
		case Elliptic:
      {
        iir_coeff CF(elliptic_order);
        elliptic_iir(CF,elliptic_pass_edge,elliptic_ripple,elliptic_stop_db);
        E_IIR.realloc(CF);
        E_IIR.reset();
        iir_freq(CF,pass_type == high,pts,bits,w,inc);
      }
			return(elliptic_pass_edge);
			break;
		case Chebyshev:
      {
        iir_coeff CF(chebyshev_order);
        fc = chebyshev_fc;
        chebyshev_iir(CF,fc,chebyshev_ripple);
        C_IIR.realloc(CF);
        C_IIR.reset();
        iir_freq(CF,pass_type == high,pts,bits,w,inc);
      }
			return(chebyshev_fc);
			break;
		case Butterworth:
      {
        iir_coeff CF(butterworth_order);
        fc = butterworth_fc;
        butterworth_iir(CF,fc,3.0);
        B_IIR.realloc(CF);
        B_IIR.reset();
        iir_freq(CF,pass_type == high,pts,bits,w,inc);
      }
			return(butterworth_fc);
			break;
			
      // Special filters
		case NotchIIR:
			NOTCH.set_coeffs(notch_fc,notch_trans);
			NOTCH.reset();
			other_freq(NOTCH,pts,w,inc);
			return(2.0*notch_fc);
			break;
		case  CutBoost:
			CUT_B.reset();
			other_freq(CUT_B,pts,w,inc);
			return(2.0*cut_fc);
			break;
		case Shelf: 
			Z1_IIR.set_coeffs(shelf_fc,shelf_gain, -shelf_gain);
			other_freq(Z1_IIR, pts, w,inc);
			break;	
			
    }
    return(0);
  }

} // namespace SPUC
