#include <eutils/emain.h>
#include <eutils/ernd.h>

#include <gsl/gsl_cdf.h>

#include <eutils/eparser.h>
#include <eutils/eparserinterpreter.h>
#include <eutils/etable.h>

#include <deque>

using namespace std;

typedef ebasicarray<uint32_t> euintarray;

// needed when using random number generator in the arguments!
inline int maxint(int a,int b) { return(a>b?a:b); }
inline int minint(int a,int b) { return(a<b?a:b); }


void multinomial(int count,const edoublearray& w,euintarray& res,ernd& r)
{
  gsl_ran_multinomial(r.grng, w.size(), count, &w[0], &res[0]);
}

int binomial(int n,double p,const ernd& r)
{
  if (p <= 0.0) return 0.0; 
  return gsl_ran_binomial(r.grng, p, n);
}

void permute(ebasicarray<uint8_t>& pop_ages,int ind,int size,ernd& r){
  for (int i=0; i<size-1; ++i){
    int ir=r.uniform()*size;
    if (ir==i) continue;
    swap(pop_ages[ind+i],pop_ages[ind+ir]);
  }
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
rgamma::rgamma(double mu,double shape,int _tmax,double _tstep): tmax(_tmax),tstep(_tstep) {
  double scale=mu/shape;
  a=shape; b=scale;
}

int rgamma::operator()(ernd& r)
{
  int tmpr=int(gsl_ran_gamma(r.grng,a,b)/tstep+tstep/2.0);
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
  uint8_t hage;
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
  void add(ssimstate& st,uint8_t hage,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r,int ntransition=0);
  void add(ssimstate& st,sevent& ev,int newt);
  void resize(int newsize);
};

struct ssimstate {
  double R0=2.7; // divided by number of days of infectiousness to obtain rate per day
  const double rSA=0.66; // ratio of Symptomatic/Asymptomatic
//  const double iIa=0.5; // infectiosity of asymptomatic
//  const double iIp=1.0; // infectiosity of presymptomatic
//  const double iIs=1.0; // infectiosity of presymptomatic

  const double tstep=0.25;

  edoublearray dE;  // 4 days delay
  edoublearray dIp; // 1.5 days non-symptomatic
  edoublearray dIs; // 3.5 days as symptomatic
  edoublearray dIa; // 5 days infectious as asymptomatic

  edoublearray dHtoicu; // 7 days to go toicu state
  edoublearray dHtononicu; // 7 days to go tononicu state

  edoublearray dHicu; // 10 days in ICU
  edoublearray dHnonicu; // 8 days in nonICU

  edoublearray dIpD; // deaths occuring 22 days after symptoms

  rgamma rIp;
  rgamma rIs;
  rgamma rIa;

  rgamma rToicu;
  rgamma rTononicu;
  rgamma rInicu;
  rgamma rInnonicu;
  rgamma rIpDeath;

  edoublearray icu_symp;
  edoublearray nonicu_symp;
  edoublearray death_icu;
  edoublearray death_nonicu;

  ebasicarray<uint8_t> pop_ages;
  ebasicarray<shousehold> households;

  earray<earray<deque<int> > > hhLevels;
  earray<eintarray> hhIa,hhIp,hhIs;

  int allE=0,allIp=0,allIs=0,allIa=0;

  double allICU,allNonICU,allD;

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
          ++st.hhIp[hs.hage][hhl];
          ++hs.Ip;
          e.transition=2;
          st.evqueue.add(st,e,st.rIp(r));
        }else{
          // E -> Ia
          ++st.allIa;
          ++st.hhIa[hs.hage][hhl];
          ++hs.Ia;
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
        --st.hhIa[hs.hage][hhl];
        --hs.Ia;
       break;
      }
      case 2:{
        // Ip -> Is , Ip -> toicu , Ip -> tononicu , Ip -> death
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allIp;
        --st.hhIp[hs.hage][hhl];
        --hs.Ip;
        ++st.allIs;
        ++st.hhIs[hs.hage][hhl];
        ++hs.Is;
        e.transition=3;
        st.evqueue.add(st,e,st.rIs(r));
        double rf=rnd.uniform();
        rf-=st.icu_symp[e.age];
        if (rf<0.0){
          e.transition=4;
          st.evqueue.add(st,e,st.rToicu(r));
          if (rnd.uniform()<st.death_icu[e.age]){
            e.transition=6;
            st.evqueue.add(st,e,st.rIpDeath(r));
          }
        }else{
          rf-=st.nonicu_symp[e.age];
          if (rf<0.0){
            e.transition=5;
            st.evqueue.add(st,e,st.rTononicu(r));
            if (rnd.uniform()<st.death_nonicu[e.age]){
              e.transition=6;
              st.evqueue.add(st,e,st.rIpDeath(r));
            }
          }
        }
       break;
      }
      case 3:{
        // Ip -> 
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allIs;
        --st.hhIs[hs.hage][hhl];
        --hs.Is;
       break;
      }
      case 4:{
        // toicu -> icu
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        ++st.allICU;
        e.transition=7;
        st.evqueue.add(st,e,st.rInicu(r));
       break;
      }

      case 5:{
        // tononicu -> nonicu
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        ++st.allNonICU;
        e.transition=8;
        st.evqueue.add(st,e,st.rInnonicu(r));
       break;
      }
      case 6:{
        // death
        ++st.allD;
       break;
      }
      case 7:{
        // icu -> 
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allICU;
       break;
      }
      case 8:{
        // nonicu -> 
        shousehold &hs(st.households[e.hhid]);
        uint8_t hhl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
        --st.allNonICU;
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

void evqueuelist::add(ssimstate &st,uint8_t hage,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r,int ntransition)
{
  sevent ev;
  uint8_t hl=hhsize*int(hhsize+1)/2+(hhsize-hhexp);
//  uint8_t hhstnewlevel=hhsize*(hhsize+1)/2+(hhsize-(hhexp+1)); // increasing exposed count in household

  ev.transition=ntransition; // default: 0 = end of E state
 
  for (int ic=1; ic<count.size(); ++ic){
    if (count[ic]==0) continue;

    uint8_t hlnew=hl-ic; // transition level depends on number of new infections in household
    tmpres.init(evdist.size(),0);
    multinomial(count[ic],evdist,tmpres,r);
    for (int i=0; i<tmpres.size(); ++i){
      for (int j=0; j<tmpres[i]; ++j){
        ldieif(st.hhLevels[hage][hl].size()==0,"hhLevels empty!");

        // choose a random household
        int ri=r.uniform()*st.hhLevels[hage][hl].size();
        ldieif(ri>=st.hhLevels[hage][hl].size(),"ri > array");

        ev.hhid=st.hhLevels[hage][hl][ri];
        if (ri!=st.hhLevels[hage][hl].size()-1)
          swap(st.hhLevels[hage][hl][ri],st.hhLevels[hage][hl][st.hhLevels[hage][hl].size()-1]);
        st.households[st.hhLevels[hage][hl][ri]].hstind=ri;

//        st.hhLevels[hl].erase(st.hhLevels[hl].size()-1);
        st.hhLevels[hage][hl].pop_back();;
       
        shousehold &hst(st.households[ev.hhid]);
        ldieif(hst.hage!=hage,"mismatching hages!");

        // remove Ia, Ip, and Is counts from previous level
        st.hhIa[hage][hl]-=hst.Ia;
        st.hhIp[hage][hl]-=hst.Ip;
        st.hhIs[hage][hl]-=hst.Is;

        if (hst.size != hhsize || hst.E != hhexp){
          printf("wrong transition: hst.size: %hhi hhsize: %hhi hst.E: %hhi hhexp: %hhi hhstlevel: %hhi\n",hst.size,hhsize,hst.E,hhexp,hl);
          exit(0);
        }

        if (hst.E>hst.size){ 
          printf("household exposed larger than household size: hst.E: %hhi hst.size: %hhi hhsize: %hhi hhexp: %hhi ic: %i count.size(): %i\n",hst.E,hst.size,hhsize,hhexp,ic,count.size());
          exit(0);
        }
        hst.hstind=st.hhLevels[hage][hlnew].size();
        st.hhLevels[hage][hlnew].push_back(ev.hhid);

        // add Ia, Ip, and Is counts to new level
        st.hhIa[hage][hlnew]+=hst.Ia;
        st.hhIp[hage][hlnew]+=hst.Ip;
        st.hhIs[hage][hlnew]+=hst.Is;

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
  ssimstate st;
  epregister2(st.R0,"r0");

  double finter=1.0;
  double fintra=1.0;

  double foldr=1.0; // no difference in isolation for households with older members
  int tostart=-1; // start isolation
  int toend=-1; // end isolation
  epregister(foldr);
  epregister(tostart);
  epregister(toend);

  double tmax=300.0;
  epregister(tmax);
  epregister(finter);
  epregister(fintra);
  estr fpop="data/switzerland.agegroups";
  epregister(fpop);
  int fseed=3;
  epregister(fseed);

  eparseArgs();

  cout << "# R0: " << st.R0 << endl;

  eintarray agegroups;

//  agegroups.init(16,500000);


  evarhash options;
  options.add("header",1);
  options.add("sep","\t");
  etable agegroup_dem(etableLoad(fpop,options));
  agegroups=mularr(sumarr(agegroup_dem["\"f\""],agegroup_dem["\"m\""]),1000.0);
  while (agegroups.size()>16){
    agegroups[15]+=agegroups[agegroups.size()-1];
    agegroups.erase(agegroups.size()-1);
  }
  // limit to 16 age groups
//  cout << agegroups << endl;
//  exit(0);

  int popsize=0;
  for (int i=0; i<agegroups.size(); ++i)
    popsize+=agegroups[i];

  cout << "# total population: " << popsize << endl;

  



/*
  edoublearray householdSizeDist(0.0,0.31,0.34,0.16,0.13,0.06);
  euintarray householdCounts;
  householdCounts.init(6,0); // how many households with single individuals, two individuals, ... max 6
  multinomial(popsize,householdSizeDist,householdCounts,rnd);
*/

  // Ferguson et al 2006 has household distribution sizes in supplementary graphs which the values below approximate
  eintarray householdSizeDist;
  householdSizeDist.init(7,0.0); // how many households with single individuals, two individuals, ... max 6
  householdSizeDist[0]=0; // household with 0 susceptible
  householdSizeDist[1]=0.25*popsize;      
  householdSizeDist[2]=0.28*popsize/2.0;
  householdSizeDist[3]=0.20*popsize/3.0; // 3 individual households, e.g.: two adults one child
  householdSizeDist[4]=0.17*popsize/4.0;
  householdSizeDist[5]=0.08*popsize/5.0;
  householdSizeDist[6]=0.02*popsize/6.0;

  int total_households=0,total_hh_individuals=0;
  for (int i=0; i<householdSizeDist.size(); ++i){
    total_households+=householdSizeDist[i];
    total_hh_individuals+=householdSizeDist[i]*i;
  }
  householdSizeDist[1]+=popsize-total_hh_individuals; // Add missing individuals as single house holds
  total_households += popsize-total_hh_individuals;

  cout << "# hh individuals: " << total_hh_individuals << "    " << popsize << endl;
  

  etable agegroup_infparams(etableLoad("data/agegroup.infparams",options));

  st.icu_symp=mularr(agegroup_infparams["Prop_symp_hospitalised"],agegroup_infparams["Prop_hospitalised_critical"]);
  st.nonicu_symp=mularr(agegroup_infparams["Prop_symp_hospitalised"],sumarr(mularr(agegroup_infparams["Prop_hospitalised_critical"],-1.0),1.0));
  st.death_icu=edoublearray(agegroup_infparams["Prop_critical_fatal"]);
  st.death_nonicu=edoublearray(agegroup_infparams["Prop_noncritical_fatal"]);

  
//  cout << agegroup_infparams << endl;
//  exit(0);


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

  st.rToicu=rgamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.rTononicu=rgamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.rInicu=rgamma(10.0,10.0,60.0,0.25); // 10 days in ICU
  st.rInnonicu=rgamma(8.0,8.0,60.0,0.25); // 8 days in nonICU
  st.rIpDeath=rgamma(22.0,22.0,60.0,0.25); // deaths occuring 22 days after symptoms

  double R0day=st.R0/5.0; // per day rate given that individuals are infectious for 5 days



/*
  st.dHtoicu=delay_gamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.dHtononicu=delay_gamma(7.0,7.0,60.0,0.25); // 7 days to go to nonICU
  st.dHicu=delay_gamma(10.0,10.0,60.0,0.25); // 10 days in ICU
  st.dHnonicu=delay_gamma(8.0,8.0,60.0,0.25); // 8 days in nonICU
  st.dIpD=delay_gamma(22.0,22.0,60.0,0.25); // deaths occuring 22 days after symptoms
*/

  st.hhLevels.init(agegroups.size());
  st.hhIa.init(agegroups.size());
  st.hhIp.init(agegroups.size());
  st.hhIs.init(agegroups.size());

  int nlevels=householdSizeDist.size()*(householdSizeDist.size()+1)/2;
  for (int i=0; i<agegroups.size(); ++i){
    st.hhLevels[i].init(nlevels);
    st.hhIa[i].init(nlevels,0);
    st.hhIp[i].init(nlevels,0);
    st.hhIs[i].init(nlevels,0);
  }
//  uint8_t hhstlevel=hhsize*(hhsize+1)/2+(hhsize-hhexp);

  cout << "# total households: " << total_households << endl;
  st.households.init(total_households);


  earray<eintarray> hhlist;
  hhlist.init(householdSizeDist.size());

  int hi=0,ai=0;
  for (int i=1; i<householdSizeDist.size(); ++i){
//    uint8_t hl=i*(i+1)/2+i;
//    st.hhLevels[hl].reserve(householdSizeDist[i]);
    for (int j=0; j<householdSizeDist[i]; ++j,++hi){
      hhlist[i].add(hi);
//      st.hhLevels[hl].push_back(hi);
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
  ldieif(ai!=popsize,"ai != popsize "+estr(ai)+ " " + estr(popsize));



  st.pop_ages.init(popsize,0);

  
  cout << householdSizeDist << endl;

  int assigned=0;
  const int agegroupband=5; // years
  eintarray tmpag(agegroups);
  cout << tmpag << endl;

  // TODO: make gaussian distributed ages, this would take care of imperfect matches
  // TODO: make this step a setup step requiring a different option to run that generates the household structure file.
  for (int i=householdSizeDist.size()-1; i>=1; --i){
    cout << "# seeding individuals in households size: " << i << "   " << assigned << "   " << householdSizeDist[i] << "    " << (i==1?0:i-2)*householdSizeDist[i] << "    "  <<  tmpag[0]+tmpag[1]+tmpag[2]+tmpag[3] << endl;
    if (i==1)
      cout << tmpag << endl;
    for (int j=0; j<hhlist[i].size(); ++j){
      shousehold &sh(st.households[hhlist[i][j]]);
      do {
        ldieif(sh.ageind<0 || sh.ageind>=st.pop_ages.size(),"out of bounds: "+estr(sh.ageind)+" "+st.pop_ages.size());
        st.pop_ages[sh.ageind]=int(5+(i-2)*0.5)+rnd.uniform()*(tmpag.size()-int(5+(i-2)*0.5)-MAX(0,2*(i-2)));  // maxint is needed when the expression contains randomly generated numbers. Using the MAX macros causes the number to be generated twice!!
        ldieif(st.pop_ages[sh.ageind]<0 || st.pop_ages[sh.ageind]>=tmpag.size(),"out of bounds: "+estr().sprintf("%hhi",st.pop_ages[sh.ageind])+" "+tmpag.size());
//        st.pop_ages[sh.ageind]=5+(i-1)+rnd.uniform()*(tmpag.size()-5-(i-1)-MAX(0,2*(i-4)));
      } while (tmpag[st.pop_ages[sh.ageind]]==0);
      --tmpag[st.pop_ages[sh.ageind]];
      ++assigned;

//      // use single agegroup for hhLevels
//      sh.hage=0; //st.pop_ages[sh.ageind];

      // age group of "oldest" person in house
      sh.hage=st.pop_ages[sh.ageind];

      if (i>1){
        do {
          ldieif(sh.ageind+1<0 || sh.ageind+1>=st.pop_ages.size(),"out of bounds: "+estr(sh.ageind+1)+" "+st.pop_ages.size());
          st.pop_ages[sh.ageind+1]=minint(st.pop_ages[sh.ageind]+int(rnd.uniform()*3)-2,tmpag.size()-1);
          ldieif(st.pop_ages[sh.ageind+1]<0 || st.pop_ages[sh.ageind+1]>=tmpag.size(),"out of bounds: "+estr().sprintf("%hhi",st.pop_ages[sh.ageind+1])+" "+tmpag.size());
        } while (tmpag[st.pop_ages[sh.ageind+1]]==0);
        --tmpag[st.pop_ages[sh.ageind+1]];
        sh.hage=(sh.hage+st.pop_ages[sh.ageind+1])/2;
        ++assigned;
      }   
      // Children
      for (int l=2; l<i; ++l){
        do {
          if (tmpag[MAX(0,st.pop_ages[sh.ageind]-9)]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+1]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+2]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+3]==0)
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-9) + 4 + int(rnd.uniform()*3);
          else if (tmpag[MAX(0,st.pop_ages[sh.ageind]-9)]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+1]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+2]==0)
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-9) + 3;
          else
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-9) + int(rnd.uniform()*3);
/*
          if (tmpag[0]+tmpag[1]+tmpag[2]+tmpag[3]+tmpag[4]==0)
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)+1);
          else
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)-2);
*/
        } while (tmpag[st.pop_ages[sh.ageind+l]]==0);
        --tmpag[st.pop_ages[sh.ageind+l]];
        ++assigned;
        sh.hage=(sh.hage*l+st.pop_ages[sh.ageind+l])/(l+1);
      }
    }
  }

  int hhOlder=0;
  int hhOlderN=0;
  for (int i=1; i<hhlist.size(); ++i){
    uint8_t hl=i*(i+1)/2+i;
    for (int j=0; j<hhlist[i].size(); ++j){
      int hage=st.households[hhlist[i][j]].hage;
      st.hhLevels[hage][hl].push_back(hhlist[i][j]);
      if (hage>=10){ ++hhOlder; hhOlderN+=i; }
    }
  }

  cout << "# hhOlder: " << hhOlder << " (" << double(hhOlder)/total_households << ")" << endl;
  cout << "# hhOlderN: " << hhOlderN << " (" << double(hhOlderN)/popsize << ")" << endl;

  // randomize ages in households since infection always progresses sequentially in the age array
  for (int i=2; i<hhlist.size(); ++i){
    uint8_t hl=i*(i+1)/2+i;
    for (int j=0; j<hhlist[i].size(); ++j){
      shousehold &hh(st.households[hhlist[i][j]]);
      permute(st.pop_ages,hh.ageind,hh.size,rnd);
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

  // seed infections
  counts.init(3,0);
  counts[1]=fseed;

//  st.allIp=fseed; // directly as presymptomatic
//  st.evqueue.add(st,6,2,0,counts,st.dIp,rnd,2); 

  st.allE=fseed;
  st.evqueue.add(st,6,2,0,counts,st.dE,rnd);  // add Exposed

  st.allICU=0.0;
  st.allNonICU=0.0;
  st.allD=0.0;

  double tpeak=0.0;
  double peakIs=0.0;
  double peakICU=0.0;
  double peakNonICU=0.0;

  double soldr=1.0;

  cout << "Time" << "\t" << "fE" << "\t" << "allE" << "\t" << "allIa" << "\t" << "allIp" << "\t" << "allIs" << "\t" << "allICU" << "\t" << "allNonICU" << "\t" << "Deaths" << endl;
  for (int it=0; it*st.tstep<tmax; ++it){

    if (tostart>=0 && it*st.tstep > tostart){
      soldr=foldr;
      tostart=-1;
    }

    if (toend>0 && it*st.tstep > toend){
      soldr=1.0;
      toend=-1;
    }

//    double interhhrate=(hage>=10?folder:1.0)*finter*st.R0*(1.0-fE)*(0.5*st.allIa+st.allIp+st.allIs)/(popsize-st.allE);
    if (peakIs < st.allIs){
      peakIs=st.allIs;
      tpeak=it*st.tstep;
    }
    if (peakICU < st.allICU)
      peakICU=st.allICU;
    if (peakNonICU < st.allNonICU)
      peakNonICU=st.allNonICU;

    double interIprob=finter*R0day*(0.5*st.allIa+st.allIp+st.allIs)/popsize;

    fE=double(st.allE)/popsize;
    for (int hage=0; hage<agegroups.size(); ++hage){
      for (int hhsize=1; hhsize<householdSizeDist.size(); ++hhsize){
        for (int hhexp=hhsize-1; hhexp>=0; --hhexp){ // have to do decreasing here to avoid double infecting the same houses in the same step

          uint8_t hl=hhsize*(hhsize+1)/2+(hhsize-hhexp);
          if (st.hhLevels[hage][hl].size()==0) continue;

          int hS=hhsize-hhexp; // number of susceptible individuals per household
          int nS=st.hhLevels[hage][hl].size()*hS; // total number of susceptible individuals in these households

          // The intra household rate below is an approximation needed to avoid computing the exact probability of new infection per household which depends on the specific combination of Ia,Is and Ip
          double intraIprob=fintra*R0day*(0.5*st.hhIa[hage][hl]+st.hhIp[hage][hl]+st.hhIs[hage][hl])/nS;
  
          double p=((hage>=10?soldr:1.0)*interIprob + intraIprob)*st.tstep; // probability of infection per susceptible person in these households
  
          mp.clear();
          mp.add(0.0);
          mp.add(p*hS);
          mp[0]+=p*hS;
//          for (int ti=0; ti<hS; ++ti){
//            mp.add(p*(hS-ti));
//            mp[0]+=p*(hS-ti); //TODO: need to add a term for combinations of ti into hS
//            p=p*p; // probability of more than 1 infection occuring
//          }
          mp[0]=1.0-mp[0]; // probability of household having no infection
          counts.init(mp.size(),0);
          multinomial(st.hhLevels[hage][hl].size(),mp,counts,rnd);
         
  //        int newE=binomial(nS,,rnd)/(hhsize-hhexp);
  //        if (newE>st.hhLevels[hhstlevel].size()) newE=st.hhLevels[hhstlevel].size(); // cap maximum number of new household exposures
          for (int ti=1; ti<counts.size(); ++ti)
            st.allE+=counts[ti]*ti;
  
          st.evqueue.add(st,hage,hhsize,hhexp,counts,st.dE,rnd);
        }
      }
    }
    deque<sevent> &nextEvents(st.evqueue.step());
    processEvents(st,nextEvents,rnd);
    cout << it*st.tstep << "\t" << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.allIa << "\t" << st.allIp << "\t" << st.allIs << "\t" << st.allICU << "\t" << st.allNonICU << "\t" << st.allD << endl;
  }

  cout << "# fE: " << double(st.allE)/popsize << " tpeak: " << tpeak << " peakIs: " << peakIs << " peakICU: " << peakICU << " peakNonICU: " << peakNonICU << " deaths: " << st.allD << endl;


  return(0);
}
