#include <eutils/emain.h>
#include <eutils/ernd.h>

#include <gsl/gsl_cdf.h>

#include <eutils/eparser.h>
#include <eutils/eparserinterpreter.h>
#include <eutils/etable.h>

#include <deque>

using namespace std;

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

edoublearray delay_gamma(double mu,double shape,double tmax,double tstep);


class rgamma
{
 public:
  double a;
  double b;
  int tmax;
  double tstep;

  rgamma();
  rgamma(double a,double b,int tmax,double tstep);
  int operator()(ernd& r);
};

rgamma::rgamma():a(1.0),b(1.0),tmax(1.0),tstep(1.0){}
rgamma::rgamma(double _a,double _b,int _tmax,double _tstep): a(_a),b(_b),tmax(_tmax),tstep(_tstep) {}

int rgamma::operator()(ernd& r)
{
  int tmpr=int(gsl_ran_gamma(r.grng,a,b)/tstep);
  if (tmpr>=tmax) tmpr=tmax-1;
  return(tmpr);
}


struct shousehold {
  uint8_t size;
  uint8_t hstind;
  uint8_t E;
  uint8_t Ia;
  uint8_t Ip;
  uint8_t Is;
  int ageind;
};

struct sevent {
  int hhid;
  uint8_t transition;
  uint8_t age;
};

class evqueuelist;
struct ssimstate;


class evqueuelist
{
  int ip;
  earray<deque<sevent> > eventarr;
  euintarray tmpres;

 public:
  evqueuelist();

  deque<sevent>& step();
  void add(ssimstate& st,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r);
  void add(ssimstate& st,sevent& ev,int newt);
  void resize(int newsize);
};

struct ssimstate {
  const double R0=2.7/5.0; // divided by number of days of infectiousness to obtain rate per day
  const double rSA=0.5; // ratio of Symptomatic/Asymptomatic
  const double iIa=0.5; // infectiosity of asymptomatic
  const double iIp=0.5; // infectiosity of presymptomatic

  const double tstep=0.25;

  edoublearray dE;  // 4 days delay
  edoublearray dIp; // 1.5 days non-symptomatic
  edoublearray dIs; // 3.5 days as symptomatic
  edoublearray dIa; // 5 days infectious as asymptomatic

  edoublearray dHtoicu; // 7 days to go to ICU
  edoublearray dHtononicu; // 7 days to go to nonICU

  edoublearray dHicu; // 10 days in ICU
  edoublearray dHnonicu; // 8 days in nonICU

  edoublearray dIpD; // deaths occuring 22 days after symptoms

  rgamma rIp;
  rgamma rIs;
  rgamma rIa;

  ebasicarray<uint8_t> pop_ages;
  ebasicarray<shousehold> households;

  earray<eintarray> hhLevels;
  eintarray hhIa,hhIp,hhIs;

  int allE=0,allIp=0,allIs=0,allIa=0;

  double all_icu,all_nonicu;

  evqueuelist evqueue;
};

void processEvents(ssimstate& st,deque<sevent>& evs,ernd& r)
{
  for (int i=0; i<evs.size(); ++i){
    sevent &e(evs[i]);
    switch(e.transition){
      case 0:{ // end of E state
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        if (rnd.uniform()<st.rSA){
          // E -> Ip
          ++st.allIp;
          ++st.hhIp[hhl];
          e.transition=2;
          st.evqueue.add(st,e,st.rIp(r));
        }else{
          // E -> Ia
          ++st.allIa;
          ++st.hhIa[hhl];
          e.transition=1;
          st.evqueue.add(st,e,st.rIa(r));
        }
       break;
      }
      case 1:{
        // Ia ->
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allIa;
        --st.hhIa[hhl];
       break;
      }
      case 2:{
        // Ip -> Is
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allIp;
        --st.hhIp[hhl];
        ++st.allIs;
        ++st.hhIs[hhl];
        e.transition=3;
        st.evqueue.add(st,e,st.rIs(r));
       break;
      }
      case 3:{
        // Is ->
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allIs;
        --st.hhIs[hhl];
       break;
      }
    }
  }
  evs.clear();
}


evqueuelist::evqueuelist(): ip(0) {}

void evqueuelist::add(ssimstate &st,sevent& ev,int newt)
{
//  tmpres.init(evdist.size(),0);
//  multinomial(1,evdist,tmpres,r);
//  for (int i=0; i<tmpres.size()-1; ++i){
//    if (tmpres[i]){
//      eventarr[(ip+i)%eventarr.size()].push_back(ev);
//      break;
//    }
//  }

  ldieif(newt>eventarr.size(),"should not happen");
  eventarr[(ip+newt)%eventarr.size()].push_back(ev);
}

void evqueuelist::add(ssimstate &st,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r)
{
  sevent ev;
  uint8_t hl=hhsize*int(hhsize+1)/2+(hhsize-hhexp);
//  uint8_t hhstnewlevel=hhsize*(hhsize+1)/2+(hhsize-(hhexp+1)); // increasing exposed count in household

  ev.transition=0; // end of E state
 
  for (int ic=1; ic<count.size(); ++ic){
    uint8_t hlnew=hl-ic; // transition level depends on number of new infections in household
    tmpres.init(evdist.size(),0);
    multinomial(count[ic],evdist,tmpres,r);
    for (int i=0; i<tmpres.size(); ++i){
      for (int j=0; j<tmpres[i]; ++j){
        ldieif(st.hhLevels[hl].size()==0,"hhLevels empty!");
        int ri=r.uniform()*st.hhLevels[hl].size();
        ldieif(ri>=st.hhLevels[hl].size(),"ri > array");

        ev.hhid=st.hhLevels[hl][ri];
        if (ri!=st.hhLevels[hl].size()-1)
          st.hhLevels[hl].swap(ri,st.hhLevels[hl].size()-1);
        st.households[st.hhLevels[hl][ri]].hstind=ri;

        st.hhLevels[hl].erase(st.hhLevels[hl].size()-1);
       
        shousehold &hst(st.households[ev.hhid]);
        if (hst.size != hhsize || hst.E != hhexp){
          printf("wrong transition: hst.size: %hhi hhsize: %hhi hst.E: %hhi hhexp: %hhi hhstlevel: %hhi\n",hst.size,hhsize,hst.E,hhexp,hl);
          exit(0);
        }

        if (hst.E>hst.size){ 
          printf("household exposed larger than household size: hst.E: %hhi hst.size: %hhi hhsize: %hhi hhexp: %hhi ic: %i count.size(): %i\n",hst.E,hst.size,hhsize,hhexp,ic,count.size());
          exit(0);
        }
        hst.hstind=st.hhLevels[hlnew].size();
        st.hhLevels[hlnew].add(ev.hhid);
        for (int l=0; l<ic; ++l){ // queue one event per person to track age of exposed individual
          ev.age=st.pop_ages[hst.ageind+hst.E+l];
          eventarr[(ip+i)%eventarr.size()].push_back(ev);
        }
        hst.E+=ic;
      }
    }
  }
}

void evqueuelist::resize(int newsize)
{
  eventarr.init(newsize);
}

deque<sevent>& evqueuelist::step()
{
  deque<sevent>& tmpevarr(eventarr[ip]);
  ip=(ip+1)%eventarr.size();
  return(tmpevarr);
}


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

int emain()
{
  double finter=0.5;
  double fintra=1.5;
  epregister(finter);
  epregister(fintra);
  eparseArgs();

  eintarray agegroups;

  agegroups.init(18,500000);


  int popsize=0;
  for (int i=0; i<agegroups.size(); ++i)
    popsize+=agegroups[i];

  cout << "# total population: " << popsize << endl;

  eintarray S(agegroups); // susceptible
  eintarray F;
  F.init(agegroups.size(),0);

  eintarray R;
  R.init(agegroups.size(),0);

  eintarray H;
  H.init(agegroups.size(),0);

  




  // Ferguson et al 2006 has household distribution sizes in supplementary graphs which the values below approximate
  eintarray householdSizeDist;
  householdSizeDist.init(6,0.0); // how many households with single individuals, two individuals, ... max 6
  householdSizeDist[0]=0; // household with 0 susceptible
  householdSizeDist[1]=0.31*popsize;      
  householdSizeDist[2]=0.34*popsize/2.0;
  householdSizeDist[3]=0.16*popsize/3.0; // 3 individual households, e.g.: two adults one child
  householdSizeDist[4]=0.13*popsize/4.0;
  householdSizeDist[5]=0.06*popsize/5.0;

  int total_households=0;
  for (int i=0; i<householdSizeDist.size(); ++i)
    total_households+=householdSizeDist[i];

  ssimstate st;

  evarhash options;
  options.add("header",1);

  etable agegroup_infparams(etableLoad("data/agegroup.infparams",options));

  cout <<  mularr(agegroup_infparams["Prop_symp_hospitalised"],agegroup_infparams["Prop_hospitalised_critical"]) << endl;
  cout <<  mularr(agegroup_infparams["Prop_symp_hospitalised"]),addarr(mularr(agegroup_infparams["Prop_hospitalised_critical"],-1.0),1.0)) << endl;
  cout <<  edoublearray(agegroup_infparams["Prop_critical_fatal"]) << endl;
  cout <<  edoublearray(agegroup_infparams["Prop_noncritical_fatal"]) << endl;
//  cout << agegroup_infparams << endl;
  exit(0);


  // from Davies et al. 2020 covidm_params.R 
//  edoublearray dE(delay_gamma(4.0,4.0,60.0,0.25));  // Derived from Backer et al Eurosurveillance
//  edoublearray dIp(delay_gamma(2.4,4.0,60.0,0.25)); // Derived from Backer et al Eurosurveillance
//  edoublearray dIa(delay_gamma(7.0,4.0,60.0,0.25)); // Assumed 7 days subclinical shedding
//  edoublearray dIs(delay_gamma(3.2,3.7,60.0,0.25)); // Zhang et al 2020

  // from Davies et al. 2020 UK.R 
  st.dE=delay_gamma(4.0,4.0,60.0,0.25);  // 4 days delay
  st.dIp=delay_gamma(1.5,4.0,60.0,0.25); // 1.5 days non-symptomatic
  st.dIs=delay_gamma(3.5,4.0,60.0,0.25); // 3.5 days as symptomatic
  st.dIa=delay_gamma(5.0,4.0,60.0,0.25); // 5 days infectious as asymptomatic

  st.rIp=rgamma(1.5,4.0,60.0,0.25); // 1.5 days non-symptomatic
  st.rIs=rgamma(3.5,4.0,60.0,0.25); // 3.5 days as symptomatic
  st.rIa=rgamma(5.0,4.0,60.0,0.25); // 5 days infectious as asymptomatic


  st.dHtoicu=delay_gamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.dHtononicu=delay_gamma(7.0,7.0,60.0,0.25); // 7 days to go to nonICU

  st.dHicu=delay_gamma(10.0,10.0,60.0,0.25); // 10 days in ICU
  st.dHnonicu=delay_gamma(8.0,8.0,60.0,0.25); // 8 days in nonICU

  st.dIpD=delay_gamma(22.0,22.0,60.0,0.25); // deaths occuring 22 days after symptoms

  int nlevels=householdSizeDist.size()*(householdSizeDist.size()+1)/2;
  st.hhLevels.init(nlevels);
  st.hhIa.init(nlevels,0);
  st.hhIp.init(nlevels,0);
  st.hhIs.init(nlevels,0);
//  uint8_t hhstlevel=hhsize*(hhsize+1)/2+(hhsize-hhexp);

  cout << "# total households: " << total_households << endl;
  st.households.init(total_households);
  int hi=0,ai=0;
  for (int i=0; i<householdSizeDist.size(); ++i){
    uint8_t hl=i*(i+1)/2+i;
    st.hhLevels[hl].reserve(householdSizeDist[i]);
    for (int j=0; j<householdSizeDist[i]; ++j,++hi){
      st.hhLevels[hl].add(hi);
      st.households[hi].size=i;
      st.households[hi].hstind=j;
      st.households[hi].Ia=0;
      st.households[hi].Is=0;
      st.households[hi].Ip=0;
      st.households[hi].E=0;
      st.households[hi].ageind=ai;
      ai+=i; // add size
    }
  }
  ldieif(hi!=total_households,"hi != total_households: "+estr(hi)+ " " + estr(total_households));



  st.pop_ages.init(popsize,0);

  
  int assigned=0;
  const int agegroupband=5; // years
  eintarray tmpag(agegroups);
  for (int i=householdSizeDist.size()-1; i>=1; --i){
    uint8_t hl=i*(i+1)/2+i;
    for (int j=0; j<st.hhLevels[hl].size(); ++j){
      shousehold &sh(st.households[st.hhLevels[hl][j]]);
      do {
        st.pop_ages[sh.ageind]=5+(i-2)+rnd.uniform()*(tmpag.size()-5-3*(i-2));
      } while (tmpag[st.pop_ages[sh.ageind]]==0);
      --tmpag[st.pop_ages[sh.ageind]];
      ++assigned;
      if (i>1){
        do {
          st.pop_ages[sh.ageind+1]=st.pop_ages[sh.ageind]+int(rnd.uniform()*3)-2;
        } while (tmpag[st.pop_ages[sh.ageind+1]]==0);
        --tmpag[st.pop_ages[sh.ageind+1]];
        ++assigned;
      }   
      // Children
      for (int l=2; l<i; ++l){
        do {
          st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)-2);
        } while (tmpag[st.pop_ages[sh.ageind+l]]==0);
        --tmpag[st.pop_ages[sh.ageind+l]];
        ++assigned;
      }
    }
  }


/*
  for (int i=0; i<100; ++i){
    int ri=rnd.uniform()*st.households.size();
    shousehold &hh(st.households[ri]);
    cout << estr().sprintf("%hi",hh.size);
    for (int j=0; j<hh.size; ++j)
      cout << "\t" << pop_ages[hh.ageind+j]*5;
    cout << endl;
  }
*/

  cout << "# done populating households: " << assigned << endl;

  st.evqueue.resize(st.dE.size()+1); // keep 1 extra for using as buffer for next step when adding new events

  edoublearray mp;
  euintarray counts;

  double fE=0.0; // fraction of exposed (includes deaths, active infections and healed)

  st.allE=3;
  counts.init(3,0);
  counts[1]=3;
  st.evqueue.add(st,2,0,counts,st.dE,rnd);

  cout << "Time" << "\t" << "fE" << "\t" << "allE" << "\t" << "allIa" << "\t" << "allIp" << "\t" << "allIs" << "\t" << "allICU" << "\t" << "allNonICU" << endl;
  for (int it=0; it<1000; ++it){

    st.all_icu=0.0;
    st.all_nonicu=0.0;

    fE=double(st.allE)/popsize;

    for (int hhsize=1; hhsize<householdSizeDist.size(); ++hhsize){
      for (int hhexp=0; hhexp<hhsize; ++hhexp){
        uint8_t hhstlevel=hhsize*(hhsize+1)/2+(hhsize-hhexp);
        double interhhrate=finter*(1.0-fE)*(0.5*st.allIa+0.5*st.allIp+st.allIs);
        // The intra household rate below is an approximation needed to avoid computing the exact probability of new infection per household which depends on the specific combination of Ia,Is and Ip
        double intrahhrate=fintra*(0.5*st.hhIa[hhstlevel]+0.5*st.hhIp[hhstlevel]+st.hhIs[hhstlevel]);

        int hS=hhsize-hhexp;
        int nS=st.hhLevels[hhstlevel].size()*hS; // total number of susceptible people in these households
        double p=st.R0*(interhhrate + intrahhrate)*st.tstep/nS; // probability of infection per susceptible person in these households

        mp.clear();
        mp.add(0.0);
        for (int ti=0; ti<hS; ++ti){
          mp.add(p);
          mp[0]+=p;
          p=p*p; // probability of more than 1 infection occuring
        }
        mp[0]=1.0-mp[0]; // probability of no infection
        counts.init(mp.size(),0);
        multinomial(st.hhLevels[hhstlevel].size(),mp,counts,rnd);
       
//        int newE=binomial(nS,,rnd)/(hhsize-hhexp);
//        if (newE>st.hhLevels[hhstlevel].size()) newE=st.hhLevels[hhstlevel].size(); // cap maximum number of new household exposures
        for (int ti=1; ti<counts.size(); ++ti)
          st.allE+=counts[ti]*ti;

        st.evqueue.add(st,hhsize,hhexp,counts,st.dE,rnd);
      }
      deque<sevent> &nextEvents(st.evqueue.step());
      processEvents(st,nextEvents,rnd);
//      st.allE+=nextEvents.size();
    }
    cout << it*st.tstep << "\t" << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.allIa << "\t" << st.allIp << "\t" << st.allIs << "\t" << st.all_icu << "\t" << st.all_nonicu << endl;
//   sleep(2);
  }


/*
  ebasicarray<evqueue> E,Ip,Ia,Is,toicu,tononicu,icu,nonicu;
  E.init(agegroups.size());
  Ip.init(agegroups.size());
  Ia.init(agegroups.size());
  Is.init(agegroups.size());
  toicu.init(agegroups.size());
  tononicu.init(agegroups.size());
  icu.init(agegroups.size());
  nonicu.init(agegroups.size());


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
*/

  return(0);
}
