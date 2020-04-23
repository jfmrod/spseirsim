#include <eutils/emain.h>
#include <eutils/ernd.h>

#include <gsl/gsl_cdf.h>

typedef ebasicarray<uint32_t> euintarray;


void multinomial(int count,const edoublearray& w,euintarray& res,ernd& r)
{
  gsl_ran_multinomial(r.grng, w.size(), count, &w[0], &res[0]);
}

int binomial(int n,double p,const ernd& r)
{
  if (p <= 0.0) return 0.0; 
  return gsl_ran_binomial(r.grng, p, n);
}


struct sevent {
 public:
  uint8_t hsize;
  uint8_t hid;
  uint8_t nextevent;
};

class evqueuelist
{
  int ip;
  earray<ebasicarray<sevent> > eventarr;
  euintarray tmpres;

 public:
  evqueuelist();

  ebasicarray<sevent>& step();
  void add(int count,const edoublearray& evdist,ernd& r);
  void resize(int newsize);
};



class evqueue
{
  int _total;
  int ip;
  eintarray eventarr;
  euintarray tmpres;

 public:
  evqueue();

  inline int total(){ return(_total); }

  int step();
  void add(int count,const edoublearray& evdist,ernd& r);
  void resize(int newsize);
};

evqueue::evqueue(): _total(0),ip(0) {}

int evqueue::step(){
  if (eventarr.size()==0) return(0.0);

  ip=(ip+1)%eventarr.size();
  int tmpr=eventarr[ip];
  _total-=tmpr;
  eventarr[ip]=0.0;
  return(tmpr);
}

void evqueue::resize(int newsize){
  eventarr.reserve(newsize);
  int i=0;
  for (; i<ip && eventarr.size()<newsize; ++i){
    eventarr.add(eventarr[ip]);
    eventarr[ip]=0;
  }
  for (int ti=i; i<ip; ++i){
    eventarr[i-ti]=eventarr[i];
    eventarr[i]=0;
  }
  while (eventarr.size()<newsize)
    eventarr.add(0);
}

void evqueue::add(int count,const edoublearray& evdist,ernd& r){
  if (eventarr.size()<evdist.size())
    resize(evdist.size());
  _total+=count;

  tmpres.init(evdist.size(),0);
  multinomial(count,evdist,tmpres,r);
  
  for (int i=0; i<tmpres.size(); ++i)
    eventarr[(ip+i)%eventarr.size()]+=tmpres[i];
}

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

edoublearray delay_gamma(double mu,double shape,double tmax,double tstep)
{
  double scale=mu/shape;
  edoublearray tmp;
  tmp.reserve(int(tmax/tstep));
  double sumx=0.0;
  for (int i=0; i*tstep<tmax; ++i){
    tmp.add(gsl_cdf_gamma_P(i*tstep+tstep/2.0,shape,scale)-gsl_cdf_gamma_P(MAX(0.0,i*tstep-tstep/2.0),shape,scale));
    sumx+=tmp[i];
  }
  for (int i=0; i<tmp.size(); ++i)
    tmp[i]/=sumx;
  return(tmp);
}

struct shousehold {
 uint8_t size;
 uint8_t hi;
 uint8_t E;
 uint8_t Ia;
 uint8_t Ip;
 uint8_t Is;
};


int emain()
{
  eintarray agegroups;

  agegroups.init(18,500000);

  const double R0=2.7/5.0; // divided by number of days of infectiousness to obtain rate per day
  const double rSA=0.5; // ratio of Symptomatic/Asymptomatic
  const double iIa=0.5; // infectiosity of asymptomatic
  const double iIp=0.5; // infectiosity of presymptomatic

  const double tstep=0.25;

  int popsize=0;
  for (int i=0; i<agegroups.size(); ++i)
    popsize+=agegroups[i];

  eintarray S(agegroups); // susceptible
  eintarray F;
  F.init(agegroups.size(),0);

  eintarray R;
  R.init(agegroups.size(),0);

  eintarray H;
  H.init(agegroups.size(),0);


  eintarray households;
  households.init(5,0.0); // how many households with single individuals, two individuals, ... max 6
  households[0]=0; // household with 0 susceptible
  households[1]=0.20*popsize;
  households[2]=0.20*popsize/2.0;
  households[3]=0.30*popsize/3.0; // 3 individual households, e.g.: two adults one child
  households[4]=0.20*popsize/4.0;
  households[5]=0.10*popsize/5.0;
  int total_households=0;
  for (int i=0; i<households.size(); ++i)
    total_households+=households[i];

  earray<eintarray> householdsind;

  householdsind.init(households.size());

  householdstate.init(total_households);
  for (int i=0,hi=0; i<households.size(); ++i){
    householdsind.reserve(hoseholds[i]);
    for (int j=0; j<households[i]; ++j,++hi){
      householdsind.add(hi);
      householdstate[hi].size=i;
      householdstate[hi].hi=hi;
      householdstate[hi].Ia=0;
      householdstate[hi].Is=0;
      householdstate[hi].Ip=0;
      householdstate[hi].E=0;
    }
  }

  hhstate.init(households.size(),0);
  hhstate[1]=3;  // seed 3 infected individuals in 3 two person households

  int Iexposed=binomial(S[i],R0*(1.0-fE)*(0.5*allIa+0.5*allIp+allIs)*tstep/S[i],rnd);

  hhstate[0]
  hhstate[1] // 1H -> 1H1I -> 0H
  hhstate[2] // 2H -> 2H1I -> 1H, 2H2I -> 1H1I -> ...
  hhstate[2]*2 hhstate[3]
  hhstate[4]*3+hhstate[5]*2+hhstate[6]
  hhstate[7]*4+hhstate[8]*3+hhstate[9]*2+hhstate[10]

 


  // from Davies et al. 2020 covidm_params.R 
//  edoublearray dE(delay_gamma(4.0,4.0,60.0,0.25));  // Derived from Backer et al Eurosurveillance
//  edoublearray dIp(delay_gamma(2.4,4.0,60.0,0.25)); // Derived from Backer et al Eurosurveillance
//  edoublearray dIa(delay_gamma(7.0,4.0,60.0,0.25)); // Assumed 7 days subclinical shedding
//  edoublearray dIs(delay_gamma(3.2,3.7,60.0,0.25)); // Zhang et al 2020

  // from Davies et al. 2020 UK.R 
  edoublearray dE(delay_gamma(4.0,4.0,60.0,0.25));  // 4 days delay
  edoublearray dIp(delay_gamma(1.5,4.0,60.0,0.25)); // 1.5 days non-symptomatic
  edoublearray dIs(delay_gamma(3.5,4.0,60.0,0.25)); // 3.5 days as symptomatic
  edoublearray dIa(delay_gamma(5.0,4.0,60.0,0.25)); // 5 days infectious as asymptomatic

  edoublearray dHtoicu(delay_gamma(7.0,7.0,60.0,0.25)); // 7 days to go to ICU
  edoublearray dHtononicu(delay_gamma(7.0,7.0,60.0,0.25)); // 7 days to go to nonICU

  edoublearray dHicu(delay_gamma(10.0,10.0,60.0,0.25)); // 10 days in ICU
  edoublearray dHnonicu(delay_gamma(8.0,8.0,60.0,0.25)); // 8 days in nonICU

  edoublearray dIpD(delay_gamma(22.0,22.0,60.0,0.25)); // deaths occuring 22 days after symptoms

  ebasicarray<evqueue> E,Ip,Ia,Is,toicu,tononicu,icu,nonicu;
  E.init(agegroups.size());
  Ip.init(agegroups.size());
  Ia.init(agegroups.size());
  Is.init(agegroups.size());
  toicu.init(agegroups.size());
  tononicu.init(agegroups.size());
  icu.init(agegroups.size());
  nonicu.init(agegroups.size());


//  dH missing

//  eintarray Iexposed;
//  multinomial(Iexposed,
  double fE=0.0; // fraction of exposed (includes deaths, active infections and healed)


  int allE=0,allIp=0,allIs=0,allIa=0;

  allIp=3; // seed infection
  Ip[6].add(allIp,dIp,rnd);

  earray<edoublearray> p_icu_nonicu_noh;
  p_icu_nonicu_noh.init(agegroups.size()); // probabilities of icu hospitalization, nonicu hospitalization, and no hospitalization
  for (int i=0; i<p_icu_nonicu_noh.size(); ++i){
    p_icu_nonicu_noh[i].add(0.04*0.3);
    p_icu_nonicu_noh[i].add(0.04*(1.0-0.3));
    p_icu_nonicu_noh[i].add(1.0-0.04);
  }

  euintarray c_icu_nonicu_noh;

  double all_icu,all_nonicu;

  cout << "Time" << "\t" << "allE" << "\t" << "allIa" << "\t" << "allIp" << "\t" << "allIs" << "\t" << "allICU" << "\t" << "allNonICU" << endl;
  cout << "0.0" << "\t" << allE << "\t" << allIa << "\t" << allIp << "\t" << allIs << endl;
  for (int it=0; it<1000; ++it){
    all_icu=0.0;
    all_nonicu=0.0;

    fE=double(allE)/popsize;
    for (int i=0; i<agegroups.size(); ++i){
  //    double pi=U[i]*R0*(0.5*allIa+0.5*allIp+allIs);
  //    int Iexposed=(1.0-exp(-pi*tstep));

      // S -> E
      int Iexposed=binomial(S[i],R0*(1.0-fE)*(0.5*allIa+0.5*allIp+allIs)*tstep/S[i],rnd);
      E[i].add(Iexposed,dE,rnd);
      S[i]-=Iexposed;
      allE+=Iexposed;
  
      // E -> Ip, E -> Ia
      int outE=E[i].step();
      int newIp=binomial(outE,rSA,rnd);
      int newIa=outE-newIp;
      Ip[i].add(newIp,dIp,rnd);
      Ia[i].add(newIa,dIa,rnd);
      allIa+=newIa;
      allIp+=newIp;
  
      // Ip -> Is
      int outIp=Ip[i].step();
      Is[i].add(outIp,dIs,rnd);
      allIs+=outIp;
      allIp-=outIp;
  
      // Is -> toicu , Is -> tononicu , Is -> no hospitalization
      int outIs=Is[i].step();
      allIs-=outIs;
      c_icu_nonicu_noh.init(p_icu_nonicu_noh[i].size(),0);
      multinomial(outIs,p_icu_nonicu_noh[i],c_icu_nonicu_noh,rnd);
      toicu[i].add(c_icu_nonicu_noh[0],dHtoicu,rnd);
      tononicu[i].add(c_icu_nonicu_noh[1],dHtononicu,rnd);

      // toicu -> icu
      int out_toicu=toicu[i].step();
      icu[i].add(out_toicu,dHicu,rnd);
      
      // tononicu -> nonicu
      int out_tononicu=tononicu[i].step();
      nonicu[i].add(out_tononicu,dHnonicu,rnd);

      // icu ->
      icu[i].step();

      // nonicu ->
      nonicu[i].step();
  
      // Ia -> R
      int outIa=Ia[i].step();
      allIa-=outIa;
      R[i]+=outIa;

      all_icu+=icu[i].total();
      all_nonicu+=nonicu[i].total();
    }

    cout << it*tstep << "\t" << allE << "\t" << allIa << "\t" << allIp << "\t" << allIs << "\t" << all_icu << "\t" << all_nonicu << endl;
  }

/*
  for (int i=0; i<agegroups.size(); ++i){
    int Ip2Is=Ip[i].step();
    allIp-=Ip2Is;
    allIs+=Ip2Is;
    // missing hospitalizations and deaths
  }

  for (int i=0; i<agegroups.size(); ++i){
    int Ia2R=Ia[i].step(); // Missing asymptomatic deaths?
    allIa-=Ia2R;
    H[i]+=Ia2R;
  }
*/

  return(0);
}
