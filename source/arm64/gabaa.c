/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__GABAa_S
#define _nrn_initial _nrn_initial__GABAa_S
#define nrn_cur _nrn_cur__GABAa_S
#define _nrn_current _nrn_current__GABAa_S
#define nrn_jacob _nrn_jacob__GABAa_S
#define nrn_state _nrn_state__GABAa_S
#define _net_receive _net_receive__GABAa_S 
#define release release__GABAa_S 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gmax _p[0]
#define i _p[1]
#define g _p[2]
#define Ron _p[3]
#define Roff _p[4]
#define synon _p[5]
#define DRon _p[6]
#define DRoff _p[7]
#define _g _p[8]
#define _tsav _p[9]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
#define Alpha Alpha_GABAa_S
 double Alpha = 20;
#define Beta Beta_GABAa_S
 double Beta = 0.162;
#define Cdur Cdur_GABAa_S
 double Cdur = 0.3;
#define Cmax Cmax_GABAa_S
 double Cmax = 0.5;
#define Erev Erev_GABAa_S
 double Erev = -85;
#define Rtau Rtau_GABAa_S
 double Rtau = 0;
#define Rinf Rinf_GABAa_S
 double Rinf = 0;
#define deadtime deadtime_GABAa_S
 double deadtime = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "deadtime_GABAa_S", "ms",
 "Cmax_GABAa_S", "mM",
 "Cdur_GABAa_S", "ms",
 "Alpha_GABAa_S", "/ms",
 "Beta_GABAa_S", "/ms",
 "Erev_GABAa_S", "mV",
 "Rtau_GABAa_S", "ms",
 "gmax", "umho",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double Roff0 = 0;
 static double Ron0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "deadtime_GABAa_S", &deadtime_GABAa_S,
 "Cmax_GABAa_S", &Cmax_GABAa_S,
 "Cdur_GABAa_S", &Cdur_GABAa_S,
 "Alpha_GABAa_S", &Alpha_GABAa_S,
 "Beta_GABAa_S", &Beta_GABAa_S,
 "Erev_GABAa_S", &Erev_GABAa_S,
 "Rinf_GABAa_S", &Rinf_GABAa_S,
 "Rtau_GABAa_S", &Rtau_GABAa_S,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABAa_S",
 "gmax",
 0,
 "i",
 "g",
 0,
 "Ron",
 "Roff",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 10;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _gabaa_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 5;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABAa_S /Users/christopherjewell/Desktop/USC/Research/J4 Morphology/hoc/gabaa.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "simple GABAa receptors";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int release(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DRon = ( synon * Rinf - Ron ) / Rtau ;
   DRoff = - Beta * Roff ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DRon = DRon  / (1. - dt*( ( ( ( - 1.0 ) ) ) / Rtau )) ;
 DRoff = DRoff  / (1. - dt*( ( - Beta )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int release () {_reset=0;
 {
    Ron = Ron + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / Rtau)))*(- ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - Ron) ;
    Roff = Roff + (1. - exp(dt*(( - Beta )*( 1.0 ))))*(- ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - Roff) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 0.0 ) {
     if (  ! _args[1] ) {
       synon = synon + _args[0] ;
       _args[4] = _args[2] * exp ( - Beta * ( t - dt - _args[3] ) ) ;
       _args[2] = _args[2] * exp ( - Beta * ( t - _args[3] ) ) ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron + _args[4]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / Rtau ) ) )*( - ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron + _args[4] ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff - _args[2]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - Beta )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff - _args[2] ;
         }
 _args[3] = t ;
       _args[1] = 1.0 ;
       net_send ( _tqitem, _args, _pnt, t +  Cdur , 1.0 ) ;
       }
     }
   if ( _lflag  == 1.0 ) {
     synon = synon - _args[0] ;
     _args[4] = _args[0] * Rinf + ( _args[2] - _args[0] * Rinf ) * exp ( - ( t - dt - _args[3] ) / Rtau ) ;
     _args[2] = _args[0] * Rinf + ( _args[2] - _args[0] * Rinf ) * exp ( - ( t - _args[3] ) / Rtau ) ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron - _args[2]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / Rtau ) ) )*( - ( ( ( ( synon )*( Rinf ) ) ) / Rtau ) / ( ( ( ( - 1.0 ) ) ) / Rtau ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron - _args[2] ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff + _args[4]) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - Beta )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - Beta )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff + _args[4] ;
       }
 _args[3] = t ;
     net_send ( _tqitem, _args, _pnt, t +  deadtime , 2.0 ) ;
     }
   if ( _lflag  == 2.0 ) {
     _args[1] = 0.0 ;
     }
   } }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  Roff = Roff0;
  Ron = Ron0;
 {
   synon = 0.0 ;
   Rinf = Cmax * Alpha / ( Cmax * Alpha + Beta ) ;
   Rtau = 1.0 / ( ( Alpha * Cmax ) + Beta ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = ( Ron + Roff ) ;
   i = g * ( v - Erev ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  release();
 if(error){fprintf(stderr,"at line 107 in file gabaa.mod:\n	SOLVE release METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(Ron) - _p;  _dlist1[0] = &(DRon) - _p;
 _slist1[1] = &(Roff) - _p;  _dlist1[1] = &(DRoff) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/christopherjewell/Desktop/USC/Research/J4 Morphology/hoc/gabaa.mod";
static const char* nmodl_file_text = 
  "TITLE simple GABAa receptors\n"
  "\n"
  "COMMENT\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "	Simple model for GABAa receptors\n"
  "	================================\n"
  "\n"
  "  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS\n"
  "\n"
  "    Whole-cell recorded GABA-A postsynaptic currents (Otis et al, J. Physiol. \n"
  "    463: 391-407, 1993) were used to estimate the parameters of the present\n"
  "    model; the fit was performed using a simplex algorithm (see Destexhe et\n"
  "    al., J. Neurophysiol. 72: 803-818, 1994).\n"
  "\n"
  "  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)\n"
  "\n"
  "    The simplified model was obtained from a detailed synaptic model that \n"
  "    included the release of transmitter in adjacent terminals, its lateral \n"
  "    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe\n"
  "    and Sejnowski, 1995).  Short pulses of transmitter with first-order\n"
  "    kinetics were found to be the best fast alternative to represent the more\n"
  "    detailed models.\n"
  "\n"
  "  - ANALYTIC EXPRESSION\n"
  "\n"
  "    The first-order model can be solved analytically, leading to a very fast\n"
  "    mechanism for simulating synapses, since no differential equation must be\n"
  "    solved (see references below).\n"
  "\n"
  "\n"
  "\n"
  "References\n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for\n"
  "   computing synaptic conductances based on a kinetic model of receptor binding\n"
  "   Neural Computation 6: 10-14, 1994.  \n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for\n"
  "   excitable membranes, synaptic transmission and neuromodulation using a \n"
  "   common kinetic formalism, Journal of Computational Neuroscience 1: \n"
  "   195-230, 1994.\n"
  "\n"
  "See also:\n"
  "\n"
  "   http://cns.iaf.cnrs-gif.fr\n"
  "\n"
  "Written by A. Destexhe, 1995\n"
  "27-11-2002: the pulse is implemented using a counter, which is more\n"
  "	stable numerically (thanks to Yann LeFranc)\n"
  "\n"
  "\n"
  "Modified by Andrew Knox 2014\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS GABAa_S\n"
  "	RANGE g, gmax, R\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Rinf, Rtau\n"
  "}\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	dt		(ms)\n"
  "	deadtime = 1 (ms)\n"
  "	Cmax	= 0.5	(mM)		: max transmitter concentration\n"
  "	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)\n"
  "	Alpha	= 20	(/ms mM)	: forward (binding) rate\n"
  "	Beta	= 0.162	(/ms)		: backward (unbinding) rate\n"
  "	Erev	= -85	(mV)		: reversal potential\n"
  "	gmax		(umho)		: maximum conductance\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)		: postsynaptic voltage\n"
  "	i 		(nA)		: current = g*(v - Erev)\n"
  "	g 		(umho)		: conductance\n"
  "	Rinf				: steady state channels open\n"
  "	Rtau		(ms)		: time constant of channel binding\n"
  "	synon 			: sum of weights of all synapses in the \"onset\" state\n"
  "}\n"
  "\n"
  "STATE { Ron Roff }	: initialized to 0 by default\n"
  "				: total conductances of all synapses\n"
  "				: in the \"pulse on\" and \"pulse off\" states\n"
  "\n"
  "INITIAL {\n"
  "	synon = 0\n"
  "	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)\n"
  "	Rtau = 1 / ((Alpha * Cmax) + Beta)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE release METHOD cnexp\n"
  "	g = (Ron + Roff)\n"
  "	i = g*(v - Erev)\n"
  "}\n"
  "\n"
  "DERIVATIVE release { \n"
  "	Ron' = (synon*Rinf - Ron)/Rtau\n"
  "	Roff' = -Beta*Roff\n"
  "}\n"
  "\n"
  "\n"
  "NET_RECEIVE(weight, on, r0, t0 (ms), tmp)  {\n"
  "	if (flag == 0) {\n"
  "		:spike arrived, turn on\n"
  "		if (!on) {\n"
  "			: add to synapses in onset state\n"
  "			synon = synon + weight\n"
  "			tmp = r0*exp(-Beta*(t-dt-t0)) : matches old destexhe synapses better\n"
  "			r0 = r0*exp(-Beta*(t-t0))\n"
  "			Ron = Ron + tmp\n"
  "			Roff = Roff - r0\n"
  "			t0 = t 	\n"
  "			on = 1\n"
  "			net_send(Cdur,1)\n"
  "		}\n"
  "		:otherwise ignore new events \n"
  "	}\n"
  "	if (flag == 1) {\n"
  "		:turn off synapse\n"
  "		synon = synon - weight\n"
  "		: r0 at start of offset state\n"
  "		tmp = weight*Rinf + (r0-weight*Rinf)*exp(-(t-dt-t0)/Rtau)	: matches old destexhe synapses better\n"
  "		r0 = weight*Rinf + (r0-weight*Rinf)*exp(-(t-t0)/Rtau)\n"
  "		Ron = Ron - r0\n"
  "		Roff = Roff + tmp\n"
  "		t0 = t		\n"
  "		net_send(deadtime,2)	:flag = 2\n"
  "	}\n"
  "	if (flag == 2) {\n"
  "		on = 0 :now that dead time is passed, allow activity\n"
  "	}\n"
  "}\n"
  ;
#endif
