#pragma once
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
#include <spuce/filters/fir_coeff.h>
#include <spuce/filters/iir.h>
#include <spuce/filters/fir.h>
#include <spuce/filters/cutboost.h>
#include <spuce/filters/notch_allpass.h>
//#include <spuce/filters/iir_allpass1.h>
//#include <spuce/filters/nested_iir_allpass_2.h>
//#include <spuce/filters/nested_shelf_allpass_2.h>
#include <spuce/filters/iir_allpass1_sections.h>
#include <spuce/filters/iir_allpass1_sections_variable_delay.h>
#include <spuce/filters/iir_shelf.h>
#include <spuce/filters/shelf_allpass1.h>
#include <spuce/filters/cascaded_cic.h>
#include <spuce/filters/iir_comb.h>
#include <spuce/filters/notch_comb.h>
namespace spuce {

#include "des_filter.h"

enum fil_enum {None,MaxflatHalfband, EllipticHalfband, 
			   Butterworth, Chebyshev, Elliptic,  
			   MaxflatFIR, GaussianFIR, RemezFIR, 
			   RaisedCosine, 
			   NotchIIR, CutBoost, Shelf, RootRaisedCosine};

class make_filter {
	
 enum fil_type {low, high, band};

 public:
  double elliptic_fc;
  double butterworth_fc;
  double chebyshev_fc;
  double notch_fc;
  double remez_pass_edge;
  double remez_stop_edge;
  double remez_stop_weight;
  double cut_stop_db;
  double cut_fc;
  double rc_fc;
  double rrc_fc;
  double maxflat_fc;
  double gauss_fc;
  double elliptic_halfband_ripple;
  double variable_ripple;
  double elliptic_pass_edge;
  double elliptic_stop_edge;
  double elliptic_trans;
  double remez_trans;
  double elliptic_stop_db;
  double elliptic_ripple;
  double chebyshev_ripple;
  double rc_alpha;
  double rrc_alpha;
  double notch_trans;
  double notch50_trans;
  double cut_trans;
  double shelf_low;
  double shelf_high;
  double shelf_fc;
  double shelf_gain;
  double comb_gain;
  double notch_comb_gain;
  
  double low_shelf_gain;
  double high_shelf_gain;
  double nested_k;
  bool   hpf;
	
  //  bool elliptic_halfband_hpf;

  int elliptic_halfband_rate;
  int maxflat_halfband_rate;
  int elliptic_order;
  int butterworth_order;
  int chebyshev_order;

  int gauss_taps;
  int remez_taps;
  int maxflat_taps;
  int rc_taps;
  int rrc_taps;
  int elliptic_halfband_order;
  int maxflat_halfband_order;
  int comb_rate;
  int notch_comb_rate;
  
  int pts;
  int bits;
  double fs;
  
  double bpf_freq;
  double phase;
	
  typedef double audio_data_type;
	
  fir_coeff<double> Remez_Fir;
  fir_coeff<double> Maxflat_Fir;
  fir_coeff<double> Gaussian_Fir;
  fir_coeff<double> RaisedCosine_Fir;
  fir_coeff<double> RootRaisedCosine_Fir;
	
  fir<audio_data_type > R_Fir;
  fir<audio_data_type > M_Fir;
  fir<audio_data_type > G_Fir;
  fir<audio_data_type > RC_Fir;
  fir<audio_data_type > RRC_Fir;
	
  iir<audio_data_type > B_IIR;
  iir<audio_data_type > C_IIR;
  iir<audio_data_type > E_IIR;
	
  cutboost<audio_data_type, double > CUT_B;
  notch_allpass<audio_data_type,double> NOTCH;
  notch_allpass<audio_data_type,double> N50;

  iir_allpass1_sections_variable_delay<audio_data_type,double> B_Sub;
  iir_allpass1_sections_variable_delay<audio_data_type,double> E_Sub;

  iir_shelf<audio_data_type,double> S_IIR;
  shelf_allpass1<audio_data_type,double> Z1_IIR;
	
  fil_enum shape;
  fil_enum last_shape;
  fil_type pass_type;
 
  double change_center(double f);

  double horiz_swipe(int len, bool in_passband);
  double get_fc(int len,bool in_passband);
  int get_order();
  bool is_fir();
  bool is_bpf();
  double ripple();
  double stopdB();
  void vertical_swipe(int len, bool in_passband,bool above_stop);
  double update(double *w);
  double update(double *w, double w_inc);

  void sel_filter(const char* sel);
  void change_filter(fil_enum f);
  double limit(double x, double mx, double min);
	
  make_filter();
  ~make_filter();
  void init(int points); 
  void reset();
  void clear_filters();
  void set_filter_type(int h);
  void set_fs(double f);

};
} // namespace spuce
