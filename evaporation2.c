#include "udf.h"

#include "sg_mphase.h"

#define LAT_HT 2257000
real T_SAT=300;



DEFINE_SOURCE(liq_src, cell, thread, dS, eqn)

{

Thread *pri_th, *sec_th,*t;

real m_dot_l,p_vap;
face_t f;

real vapor_p,p_operating;
Thread *tc = THREAD_SUPER_THREAD(thread); /*obtain mixture thread */
Domain *domain=Get_Domain(1);
pri_th =THREAD_SUB_THREAD(tc,0);

sec_th = THREAD_SUB_THREAD(tc,1);

p_operating = (RP_Get_Real ("operating-pressure")+C_P(cell,thread));
p_vap= (p_operating)*C_YI(cell,pri_th,0);

	if(p_vap>800&&p_vap<101325)
T_SAT=94.5992435905437-(-24.0095048284635)*log(p_vap+1211.31526043544);
else if (p_vap>1&&p_vap<=800)
	T_SAT=209.45062276709-(-8.6557493062369)*log(p_vap+0.0207409349306649);
else
	T_SAT=313.15;

if(C_T(cell, sec_th)>=T_SAT){

m_dot_l =0;

dS[eqn] = 0;

}
else {

m_dot_l =  15.1*C_VOF(cell, pri_th)*C_R(cell, pri_th)*fabs(C_T(cell, pri_th) - T_SAT)*C_YI(cell,pri_th,0)/T_SAT;

dS[eqn] = 0.;

}
 	dS[eqn]=0;
C_UDMI(cell,thread,0)=m_dot_l;
C_UDMI(cell,thread,1)=T_SAT;

return m_dot_l;

}

DEFINE_SOURCE(vap_src, cell, thread, dS, eqn)

{

Thread * sec_th, *pri_th,*t;
face_t f;
real m_dot_v;
	Thread *tc = THREAD_SUPER_THREAD(thread); /*obtain mixture thread */
Domain *domain=Get_Domain(1);
pri_th = THREAD_SUB_THREAD(tc,0);

sec_th = THREAD_SUB_THREAD(tc,1);

if(C_T(cell, pri_th)>=T_SAT){

m_dot_v = 0;

dS[eqn] = 0.;

}
else {

m_dot_v = -15.1*C_VOF(cell, pri_th)*C_R(cell, pri_th)*fabs(C_T(cell, pri_th) - T_SAT)*C_YI(cell,pri_th,0)/T_SAT;

dS[eqn] = -15.1*C_R(cell, pri_th)*fabs(C_T(cell, pri_th) - T_SAT)/T_SAT;

}
return m_dot_v;
}

DEFINE_SOURCE(enrg_src, cell, thread, dS, eqn)

{

  Thread *pri_th, *sec_th,*t;
  face_t f;
  real m_dot;
  Domain *domain=Get_Domain(1);
	Thread *tc = THREAD_SUPER_THREAD(thread); /*obtain mixture thread */
  pri_th = THREAD_SUB_THREAD(thread, 0);

  sec_th = THREAD_SUB_THREAD(thread, 1);

if(C_T(cell, thread)>=T_SAT){

m_dot = 0;

dS[eqn] =0;

}

else {

m_dot = 15.1*C_VOF(cell, pri_th)*C_R(cell, pri_th)*fabs(C_T(cell, thread) - T_SAT)/T_SAT;

dS[eqn] = 15.1*C_VOF(cell, pri_th)*C_R(cell, pri_th)/T_SAT;
}
return LAT_HT*m_dot;
}
