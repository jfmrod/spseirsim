#include <eutils/emain.h>
#include <eutils/ernd.h>

#include <gsl/gsl_cdf.h>

#include <eutils/eparser.h>
#include <eutils/eparserinterpreter.h>
#include <eutils/etable.h>
#include <eutils/vector3.h>

#include <eutils/ethread.h>
#include <eutils/edaemon.h>
#include <eutils/esystem.h>

#include <deque>

#include "videoenc.h"


#include <geotiff.h>
#include <xtiffio.h>
#include <tiffio.h>
#include <geo_normalize.h>
#include <geo_simpletags.h>
#include <geovalues.h>

#include <shapefil.h>

#include <cairo.h>

cairo_surface_t *surface=0x00;
cairo_t *cr=0x00;


using namespace std;

const int vwidth=1024;
const int vheight=768;

struct scounts {
 int allCases=0;
 int allICU=0;
 int allNonICU=0;
 int allN=0;
 int allD=0;
};


double xlonMin=180.0,xlonMax=-180.0;
double ylatMin=90.0,ylatMax=-90.0;


typedef ebasicarray<uint32_t> euintarray;

// cannot use MIN or MAX macros when using random number generator in the arguments! 
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

template<class T>
void permute(T& arr,int s,int l,ernd& r){
  for (int i=0; i<l-1; ++i){
    int ir=i+r.uniformint(l-i);
    if (ir==i) continue;
    swap(arr[s+i],arr[s+ir]);
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
  uint8_t group;
  int gridpos;
  int ageind;
};

struct sevent {
  int hhid;
  uint8_t transition;
  uint8_t hi;
};

struct sindiv {
  uint8_t age=0;
  uint8_t state=0;
};

class evqueuelist;
struct ssimstate;
struct sthreadState;

/*
class elevel
{
 eintarray *inds;
 unsigned int s;
 unsigned int e;

 eintarray li;
 public:
   elevel(eintarray *_inds,unsigned int _s, unsigned int _e): inds(_inds),s(_s),e(_e) { }
   unsigned int swapleft(unsigned int p,unsigned int l) {
     swap(&inds[s+li[l]],&inds[s+p]); // swap element with position of first element in the "level"
     p=li[l];                         // update position of element
     ++li[l];                         // increase level start to exclude the element we just move there
   }
   unsigned int swap(unsigned int p,unsigned int ol,unsigned int nl) {
     for (;ol>nl; --ol)
       p=swapleft(p,ol);
     return(p);
   }
};
*/

class evqueuelist
{
  int ip;
  earray<deque<sevent> > eventarr;
  euintarray tmpres;

 public:
  evqueuelist();

  deque<sevent>& step();
  void add(ssimstate& st,sthreadState& ths,int hpos,uint8_t hgroup,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r,int ntransition=0);
  void add(ssimstate& st,sevent& ev,int newt);
  void resize(int newsize);
};

const int ngroups=2;
const int maxhousesize=6+1; // 
const int nlevels=maxhousesize*(maxhousesize+1)/2+maxhousesize;

struct sgrid
{
  int E;
  int N[ngroups];
  int Ia[ngroups];
  int Ip[ngroups];
  int Is[ngroups];
  int hhIa[ngroups][nlevels];
  int hhIp[ngroups][nlevels];
  int hhIs[ngroups][nlevels];
};

void addCounts(scounts& dst,scounts& src)
{
  dst.allCases+=src.allCases;
  dst.allICU+=src.allICU;
  dst.allNonICU+=src.allNonICU;
  dst.allD+=src.allD;
}

void addReset(ebasicarray<scounts>& dest,ebasicarray<scounts>& src)
{
  ldieif(dest.size()!=src.size(),"uninitialized src or dest?");
  for (int i=0; i<dest.size(); ++i){
    addCounts(dest[i],src[i]);
    src[i].allCases=0;
    src[i].allICU=0;
    src[i].allNonICU=0;
    src[i].allD=0;
  }
}


typedef void (*threadFunc_t)(ssimstate&,sthreadState&,int,int);

struct sthreadState {
  ernd r;
//  int tI=0;
//  int tN=0;

  evqueuelist evqueue;

  int allE=0;
  int allCases=0;

  ebasicarray<scounts> shapeCounts;
  ebasicarray<scounts> ageCounts;

  int Ip[ngroups]={0,0};
  int Is[ngroups]={0,0};
  int Ia[ngroups]={0,0};

  int allICU=0;
  int allNonICU=0;
  int allD=0;
};

struct ssimstate {
  double R0=2.7; // divided by number of days of infectiousness to obtain rate per day
  double rSA=0.66; // ratio of Symptomatic/Asymptomatic
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

//  const int ngroups=2; // defined as global now

  ebasicarray<SHPObject*> shapes;
  ebasicarray<uint8_t> gridMask;

  int spGridW=0;
  int spGridH=0;
  int spGridSize=0;
  ebasicarray<sgrid> spGrid;

  ebasicarray<sindiv> popindiv;
  ebasicarray<shousehold> households;

  eintarray hhLevels2;
  eintarray hhLevelsBegin;
  int grow=0;
  int arow=0;
//  earray<earray<deque<int> > > hhLevels;
//  earray<eintarray> hhIa,hhIp,hhIs;

  int allE=0;
  int allCases=0;

  int N[ngroups]={0,0};
  int Ip[ngroups]={0,0};
  int Is[ngroups]={0,0};
  int Ia[ngroups]={0,0};

  int allICU=0;
  int allNonICU=0;
  int allD=0;

  ebasicarray<scounts> shapeCounts;
  ebasicarray<scounts> ageCounts;


  edoublearray localIprob;
  edoublearray travelKernel;
  int travelKernelSize=0;

//  evqueuelist evqueue;
  double soldr=1.0;
  double smr=1.0;
  double finter=1.0;
  double fintra=1.0;
  double fglobal=0.001;
  double flocal=1.0;

  double R0day=R0/5.0;

  int householdSize=7;

  earrayof<int,int> iseedArr;

//  int fseed=10;
//  int seedGrid=0;

  emutex mutex;
  econdsig stateSignal;
  econdsig doneSignal;

  int nthreads=4;

  threadFunc_t threadFunc=0x00;
  earray<sthreadState> threadStates;
  int threadI=0;
  int threadDone=0;

  eintarray gridShuffle;  // needed to load balance threads
};



void processEvents(ssimstate& st,sthreadState& ths,deque<sevent>& evs,ernd& r)
{
  for (int i=0; i<evs.size(); ++i){
    sevent &e(evs[i]);
    shousehold &hs(st.households[e.hhid]);
    sindiv &hindiv(st.popindiv[hs.ageind+e.hi]);
    uint8_t hl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
    sgrid &g(st.spGrid[hs.gridpos]);

    switch(e.transition){
      case 0:{ // end of E state
        if (r.uniform()<st.rSA){
          // E -> Ip
          ++ths.Ip[hs.group];
          ++g.Ip[hs.group];
          ++g.hhIp[hs.group][hl];
          ++hs.Ip;
          hindiv.state=2; // presymptomatic
          e.transition=2;
          ths.evqueue.add(st,e,st.rIp(r));
        }else{
          // E -> Ia
          ++ths.Ia[hs.group];
          ++g.Ia[hs.group];
          ++g.hhIa[hs.group][hl];
          ++hs.Ia;
          hindiv.state=1; // asymptomatic
          e.transition=1;
          ths.evqueue.add(st,e,st.rIa(r));
        }
       break;
      }
      case 1:{
        // Ia ->
        hindiv.state=10; // recovered
        --ths.Ia[hs.group];
        --g.Ia[hs.group];
        --g.hhIa[hs.group][hl];
        --hs.Ia;
       break;
      }
      case 2:{
        // Ip -> Is , Ip -> toicu , Ip -> tononicu , Ip -> death
        --ths.Ip[hs.group];
        --g.Ip[hs.group];
        --g.hhIp[hs.group][hl];
        --hs.Ip;

        ++ths.allCases;
        ++ths.ageCounts[hindiv.age].allCases;
        ++ths.shapeCounts[st.gridMask[hs.gridpos]].allCases;

        ++ths.Is[hs.group];
        ++g.Is[hs.group];
        ++g.hhIs[hs.group][hl];
        ++hs.Is;

        hindiv.state=3;  // symptomatic
        double rf=r.uniform();
        rf-=st.icu_symp[hindiv.age];
        if (rf<0.0){ // going to icu
          if (r.uniform()<st.death_icu[hindiv.age]){
            hindiv.state=20; // will be fatal
            e.transition=6;  // death
            ths.evqueue.add(st,e,st.rIpDeath(r));
          }
          e.transition=4; //toicu
          ths.evqueue.add(st,e,st.rToicu(r));
        }else{
          rf-=st.nonicu_symp[hindiv.age];
          if (rf<0.0){ // going to non icu
            if (r.uniform()<st.death_nonicu[hindiv.age]){
              hindiv.state=20; // will be fatal
              e.transition=6; // death
              ths.evqueue.add(st,e,st.rIpDeath(r));
            }
            e.transition=5; //tononicu
            ths.evqueue.add(st,e,st.rTononicu(r));
          }else{
            e.transition=3; // recover from symptomatic
            ths.evqueue.add(st,e,st.rIs(r));
          }
        }
       break;
      }
      case 3:{
        // Ip -> 
        hindiv.state=10; // recovered
        --ths.Is[hs.group];
        --g.Is[hs.group];
        --g.hhIs[hs.group][hl];
        --hs.Is;
       break;
      }
      case 4:{
        // TODO: remove individual from population (stop being infective) once he joins ICU or nonICU, requires keeping track of individual to make sure he is not removed twice, once from Is -> and another time from the ICU, nonICU or death cases, maybe use age value as flag when -1 means he already was removed

        if (hindiv.state==30) // death occurred before NonICU
          break;

        // toicu -> icu

        --ths.Is[hs.group];
        --g.Is[hs.group];
        --g.hhIs[hs.group][hl];
        --hs.Is;

        ++ths.allICU;
        ++ths.ageCounts[hindiv.age].allICU;
        ++ths.shapeCounts[st.gridMask[hs.gridpos]].allICU;

        if (hindiv.state!=20){ // if state == 20 then individual will die and will not leave ICU
          e.transition=7;  // recover from ICU
          ths.evqueue.add(st,e,st.rInicu(r));
        }
        hindiv.state=7; // inicu

       break;
      }

      case 5:{
        // tononicu -> nonicu
        if (hindiv.state==30) // death occurred before NonICU
          break;

        --ths.Is[hs.group];
        --g.Is[hs.group];
        --g.hhIs[hs.group][hl];
        --hs.Is;

        ++ths.allNonICU;
        ++ths.ageCounts[hindiv.age].allNonICU;
        ++ths.shapeCounts[st.gridMask[hs.gridpos]].allNonICU;

        if (hindiv.state!=20){ // if state == 20 then individual will die and will not leave ICU
          e.transition=8; // recover from nonICU
          ths.evqueue.add(st,e,st.rInnonicu(r));
        }
        hindiv.state=8; // in nonicu
       break;
      }
      case 6:{
        // death
        ++ths.allD;
        ++ths.ageCounts[hindiv.age].allD;
        ++ths.shapeCounts[st.gridMask[hs.gridpos]].allD;
        switch (hindiv.state){
          case 7: 
            --ths.allICU;
            --ths.shapeCounts[st.gridMask[hs.gridpos]].allICU;
           break;
          case 8: 
            --ths.allNonICU;
            --ths.shapeCounts[st.gridMask[hs.gridpos]].allNonICU;
           break;
          case 20:
            // death occurred before ICU or NonICU
            --ths.Is[hs.group];
            --g.Is[hs.group];
            --g.hhIs[hs.group][hl];
            --hs.Is;
           break;
         default:
           ldie("state not expected");
          break;
        }
        hindiv.state=30;
       break;
      }
      case 7:{
        // icu -> 
        hindiv.state=10; // recovered
        --ths.allICU;
        --ths.shapeCounts[st.gridMask[hs.gridpos]].allICU;
       break;
      }
      case 8:{
        // nonicu -> 
        hindiv.state=10; // recovered
        --ths.allNonICU;
        --ths.shapeCounts[st.gridMask[hs.gridpos]].allNonICU;
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

void evqueuelist::add(ssimstate &st,sthreadState& ths,int hpos,uint8_t hgroup,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r,int ntransition)
{
  sevent ev;
  uint8_t hl=hhsize*int(hhsize+1)/2+(hhsize-hhexp);
//  uint8_t hhstnewlevel=hhsize*(hhsize+1)/2+(hhsize-(hhexp+1)); // increasing exposed count in household

  for (int ti=1; ti<count.size(); ++ti)
    ths.allE+=count[ti]*ti;

  ev.transition=ntransition; // default: 0 = end of E state

  sgrid &g(st.spGrid[hpos]);
 
  for (int ic=1; ic<count.size(); ++ic){
    if (count[ic]==0) continue;

    uint8_t hlnew=hl-ic; // transition level depends on number of new infections in household
    tmpres.init(evdist.size(),0);
    multinomial(count[ic],evdist,tmpres,r);
    for (int i=0; i<tmpres.size(); ++i){
      for (int j=0; j<tmpres[i]; ++j){
//        ldieif(st.hhLevels[hage][hl].size()==0,"hhLevels empty!");
        unsigned int hlb=st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl];
        unsigned int hlcount=st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl+1]-st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl];
        lddieif(hlcount==0,"hhLevels empty!");

        // choose a random household
        int ri=r.uniformint(hlcount);

        unsigned int hhid=st.hhLevels2[hlb+ri];
        ev.hhid=hhid;
        shousehold &hst(st.households[hhid]);

        if (hst.gridpos != hpos || hst.size != hhsize || hst.E != hhexp || hst.group != hgroup){
          printf("wrong transition: hst.gridpos: %i hpos: %i hst.size: %hhi hhsize: %hhi hst.E: %hhi hhexp: %hhi hhstlevel: %hhi\n",hst.gridpos,hpos,hst.size,hhsize,hst.E,hhexp,hl);
          exit(0);
        }

        if (hst.E>hst.size){ 
          printf("household exposed larger than household size: hst.E: %hhi hst.size: %hhi hhsize: %hhi hhexp: %hhi ic: %i count.size(): %i\n",hst.E,hst.size,hhsize,hhexp,ic,count.size());
          exit(0);
        }

        // remove Ia, Ip, and Is counts from previous level
        g.hhIa[hgroup][hl]-=hst.Ia;
        g.hhIp[hgroup][hl]-=hst.Ip;
        g.hhIs[hgroup][hl]-=hst.Is;

        // this code is implemented like this for efficiency reasons
        // households indices in the different levels are all stored in a single large array grouped by levels
        // the start of each subarray level is stored in the hhLevelsBegin array
        // to move a household to another level, I just swap the household to the beginning of the level subarray
        // and then increment the beginning of that subarray level, effectively removing it from the current level to the previous one
        // each level corresponds to a different number of exposed individuals in the households
        for (int l=0; l<ic; ++l){
          // move household chosen to beginning of level
          hlb=st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl-l];
          if (ri>0)
            swap(st.hhLevels2[hlb],st.hhLevels2[hlb+ri]); 
          st.households[st.hhLevels2[hlb+ri]].hstind=ri; // update position of not chosen but swapped household
          ++st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl-l];
          ri=st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl-l]-st.hhLevelsBegin[hpos*st.grow + hgroup*st.arow + hl-l-1]-1;
          hst.hstind=ri; // update position of moved household, to end of lower level
        }

        // add Ia, Ip, and Is counts to new level
        g.hhIa[hgroup][hlnew]+=hst.Ia;
        g.hhIp[hgroup][hlnew]+=hst.Ip;
        g.hhIs[hgroup][hlnew]+=hst.Is;

        for (int l=0; l<ic; ++l){ // queue one event per person to track age of exposed individual
          ev.hi=hst.E+l;
          eventarr[(ip+i)%eventarr.size()].push_back(ev);
        }
        hst.E+=ic;
        g.E+=ic;
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

int samplehh(eintarray& arr,ebasicarray<shousehold>& hhs,int s,ernd& r){
  int scount;
  int i;
  int tries=0;
  do{
    scount=0;
    for (i=0; i<arr.size() && scount < s; ++i){
      int ir=r.uniformint(arr.size()-i);
      scount+=hhs[arr[ir]].size;
      if (ir!=arr.size()-1-i)
        swap(arr[ir],arr[arr.size()-1-i]);
    }
    ++tries;
  } while (scount!=s && tries<50);
  return(i);
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

void colorRect(uint8_t *frameRaw,int xb,int yb,int xe,int ye,uint32_t color,int prow){
  if (xb>xe) swap(xb,xe);
  if (yb>ye) swap(yb,ye);
  for (int ty=yb; ty<ye; ++ty){
    for (int tx=xb; tx<xe; ++tx){
      int p=ty*prow*4+tx*4;
      frameRaw[p]=0xFF&color;
      frameRaw[p+1]=0xFF&(color>>8);
      frameRaw[p+2]=0xFF&(color>>16);
    }
  }
}

inline double clamp(double vmin,double vmax,double v){ v=(v<vmin?vmin:v); return(v>vmax?vmax:v); }


//SHPObject *psShape=0x00;


static int GTIFReportACorner( GTIF *gtif, GTIFDefn *defn, FILE * fp_out, const char * corner_name, double x, double y)
{
//  FILE *fp_out=stdout;
  double x_saved, y_saved;

  /* Try to transform the coordinate into PCS space */
  if ( !GTIFImageToPCS( gtif, &x, &y ) )
    return -1;
    
  x_saved = x;
  y_saved = y;

  fprintf( fp_out, "%-13s ", corner_name );

  if( defn->Model == ModelTypeGeographic )
  {
    fprintf( fp_out, "dec (%.7f,", x );
    fprintf( fp_out, "%.7f)\n", y );
    fprintf( fp_out, "(%s,", GTIFDecToDMS( x, "Long", 2 ) );
    fprintf( fp_out, "%s)\n", GTIFDecToDMS( y, "Lat", 2 ) );
  }
  else
  {
    fprintf( fp_out, " dec (%12.3f,%12.3f)", x, y );

    if ( GTIFProj4ToLatLong( defn, 1, &x, &y ) ) {
      fprintf( fp_out, "  (%.7f,", x );
      fprintf( fp_out, "%.7f)", y );
    } 
    else 
    {
      const char* pszLong = GTIFDecToDMS( x, "Long", 2 );
      if ( pszLong[0] == 0 ){
        fprintf( fp_out, "  (invalid)" );
      } else {
        fprintf( fp_out, "  (%s,", pszLong );
        fprintf( fp_out, "%s)", GTIFDecToDMS( y, "Lat", 2 ) );
      }
    }
  }

  fprintf( fp_out, "\n" );

//    if( inv_flag && GTIFPCSToImage( gtif, &x_saved, &y_saved ) ){
//        fprintf( fp_out, "      inverse (%11.3f,%11.3f)\n", x_saved, y_saved );
//    }  
  return 0;
}


GTIF *gtif=0x00;
float nodataval=0.0;
int rmaxpop=0.0;
int rpopsize=0;
int rwidth=0;
int rheight=0;
int imaxpop=-1;

eintarray popCounts;

inline int readPopDensR(uint32_t ix,uint32_t iy){
//  return(raster[iy*rwidth + ix]==nodataval?0.0:raster[iy*rwidth + ix]); // nodata values are already zeroed
  return(popCounts[iy*rwidth + ix]); 
}


uint32_t v3_col(const evector3& v){ return(0xFFu&uint32_t(v.x*0xFFu) | 0xFF00u&uint32_t(v.y*0xFF00u) | 0xFF0000u&uint32_t(v.z*0xFF0000)); }

uint32_t tricolorint(const evector3& col1,const evector3& col2,const evector3& col3,double p1,double p2){
  evector3 col12=col1*(1.0-p1)+col2*p1;
  return(v3_col((col1*(1.0-p1)+col2*p1)*(1.0-p2)+col3*p2));
}

evector3 tricolor(const evector3& col1,const evector3& col2,const evector3& col3,double p1,double p2){
  evector3 col12=col1*(1.0-p1)+col2*p1;
  return((col1*(1.0-p1)+col2*p1)*(1.0-p2)+col3*p2);
}


double plonMin=0.0;
double plonMax=0.0;
double lonRef=0.0;

double proj(double lon,double lat,double reflon){
  return(cos(M_PI*lat/180.0)*(lon-reflon));
}

double scale=1.0;
double xpos=0.0,ypos=0.0;

double revproj(double lon,double lat,double reflon){
  return(lon/cos(M_PI*lat/180.0)+reflon + plonMin);
}


void cairoDrawShape(cairo_t *cr,SHPObject *shp,float x,float y,float s,bool projection,float height,bool inverted){
  for (int j=0, iPart=1; j<shp->nVertices; ++j){
    float sx=s*(shp->padfX[j]-xlonMin)+x;
    if (projection)
      sx=s*(proj(shp->padfX[j],shp->padfY[j],lonRef)-plonMin)+x;
    
    float sy=s*(shp->padfY[j]-ylatMin)+y;
    if (inverted)
      sy=height-sy;
    if (j == 0 && shp->nParts > 0 )
      cairo_move_to(cr,sx,sy);
//       pszPartType = SHPPartTypeName( psShape->panPartType[0] );
            
    if (iPart < shp->nParts && shp->panPartStart[iPart]==j){
      cairo_close_path(cr);
      cairo_move_to(cr,sx,sy);
//       pszPartType = SHPPartTypeName( psShape->panPartType[0] );
//      pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
      iPart++;
    }else
      cairo_line_to(cr,sx,sy);
  }
}

/*
// multithreaded rendering
void threadRenderFrame(ssimstate& st,int i,int n){
  for (int i=st.spGrid.size()*i/n; i<st.spGrid.size()*(i+1)/n && i<st.spGrid.size(); ++i){
    int xi=i%st.spGridW;
    int yi=i/st.spGridH;

    double xlon=xi*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
    double ylat=(st.spGridH-yi-1)*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

    double xlonn=(xi+1)*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
    double ylatn=(st.spGridH-yi)*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

    double sx=scale*(proj(xlon,ylat,lonRef)-lonMin)+xpos;
    double sy=height-scale*(ylat-psShape->dfYMin)-ypos;

    double sxn=scale*(proj(xlonn,ylatn,lonRef)-lonMin)+xpos+1.0; // 1.0 is added to avoid tears
    double syn=height-scale*(ylatn-psShape->dfYMin)-ypos-1.0;

    sgrid &g(st.spGrid[i]);
    double tmpI=double(g.Ia[0]+g.Ia[1]+g.Ip[0]+g.Ip[1]+g.Is[0]+g.Is[1])/(g.N[0]+g.N[1]);
    double tmpE=(g.N[0]+g.N[1]==0?0.0:double(g.E)/(g.N[0]+g.N[1]));
    float pdens=clamp(0.0,1.0,(log(readPopDensR(i,j)+0.5)-log(0.5))/(log(rmaxpop+0.5)-log(0.5))*0.6+0.4);
    uint32_t color=tricolorint(evector3(1.0,1.0,1.0)*pdens,evector3(0.3,0.8,0.3),evector3(1.0,0.0,0.0),tmpE,clamp(0.0,1.0,(log(tmpI+0.001)-log(0.001))/(-log(0.001))));
    if (st.gridMask[i]==0)
      color=0xFFFFFF;
    colorRect(frameRaw,sx,sy,sxn,syn,color,width);
  }
}
*/

void renderFrame(char *daystr,float mitigation,int newCases,uint8_t *frameRaw,ssimstate& st,int width,int height){
//  memset(frameRaw, 0, 1920*1080*4);
//  videoPushFrame(frameRaw);
//  return;
  if ((plonMax-plonMin)/(ylatMax-ylatMin) > double(width)/height){
    scale=width/(plonMax-plonMin);
    ypos=(height-scale*(ylatMax-ylatMin))/2.0;
  }else{
    scale=height/(ylatMax-ylatMin);
    xpos=(width-scale*(plonMax-plonMin))/2.0;
  }
  scale=0.9*scale;
  ypos=(height-scale*(ylatMax-ylatMin))/2.0;
  xpos=(width-scale*(plonMax-plonMin))/2.0;

  double maxE=0.0;
  double maxI=0.0;
  for (int i=0; i<st.spGridSize; ++i){
    sgrid &g(st.spGrid[i]);
    double tmpI=double(g.Ia[0]+g.Ia[1]+g.Ip[0]+g.Ip[1]+g.Is[0]+g.Is[1])/(g.N[0]+g.N[1]);
    double tmpE=double(g.E)/(g.N[0]+g.N[1]);
    if (maxI<tmpI) maxI=tmpI;
    if (maxE<tmpE) maxE=tmpE;
  }

  for (int i=0; i<vwidth*vheight*4; i+=4){
    frameRaw[i]=0xff;
    frameRaw[i+1]=0xff;
    frameRaw[i+2]=0xff;
    frameRaw[i+3]=0x00;
  }

  for (int i=0; i<st.spGridW; ++i){
    for (int j=0; j<st.spGridH; ++j){
      if (st.gridMask[j*rwidth+i]==0) continue;
      double xlon=i*(xlonMax-xlonMin)/st.spGridW+xlonMin;
      double ylat=(st.spGridH-j-1)*(ylatMax-ylatMin)/st.spGridH+ylatMin;

      double xlonn=(i+1)*(xlonMax-xlonMin)/st.spGridW+xlonMin;
      double ylatn=(st.spGridH-j)*(ylatMax-ylatMin)/st.spGridH+ylatMin;

      double sx=scale*(proj(xlon,ylat,lonRef)-plonMin)+xpos;
      double sy=height-scale*(ylat-ylatMin)-ypos;

      double sxn=scale*(proj(xlonn,ylatn,lonRef)-plonMin)+xpos+1.0; // 1.0 is added to avoid tears
      double syn=height-scale*(ylatn-ylatMin)-ypos-1.0;

      sgrid &g(st.spGrid[j*st.spGridW+i]);
      double tmpI=double(g.Ia[0]+g.Ia[1]+g.Ip[0]+g.Ip[1]+g.Is[0]+g.Is[1])/(g.N[0]+g.N[1]);
      double tmpE=(g.N[0]+g.N[1]==0?0.0:double(g.E)/(g.N[0]+g.N[1]));
      float pdens=clamp(0.0,1.0,(log(readPopDensR(i,j)+0.5)-log(0.5))/(log(rmaxpop+0.5)-log(0.5))*0.6+0.4);
      uint32_t color=tricolorint(evector3(1.0,1.0,1.0)*pdens,evector3(0.3,0.8,0.3),evector3(1.0,0.0,0.0),tmpE,clamp(0.0,1.0,(log(tmpI+0.001)-log(0.001))/(-log(0.001))));
      if (st.gridMask[j*st.spGridW+i]==0)
        color=0xFFFFFF;
      colorRect(frameRaw,sx,sy,sxn,syn,color,width);
    }
  }

/*
  cairo_rectangle(cr,0.0,0.0,vwidth,vheight);
  cairo_set_source_rgb(cr, 1.0,1.0,1.0);
  cairo_fill(cr);

  evector3 color;
  for (int i=0; i<st.spGridW; ++i){
    for (int j=0; j<st.spGridH; ++j){
//      double sxlon=i*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
//      double sylat=j*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

      double xlon=i*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
      double ylat=(st.spGridH-j-1)*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

      double xlonn=(i+1)*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
      double ylatn=(st.spGridH-j)*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

      double sx=scale*(proj(xlon,ylat,lonRef)-lonMin)+xpos;
      double sy=height-scale*(ylat-psShape->dfYMin)-ypos;

      double sxn=scale*(proj(xlonn,ylatn,lonRef)-lonMin)+xpos+1.0; // 1.0 is added to avoid tears
      double syn=height-scale*(ylatn-psShape->dfYMin)-ypos-1.0;

      sgrid &g(st.spGrid[j*st.spGridW+i]);
      double tmpI=double(g.Ia[0]+g.Ia[1]+g.Ip[0]+g.Ip[1]+g.Is[0]+g.Is[1])/(g.N[0]+g.N[1]);
//      double tmpE=double(g.E)/(g.N[0]+g.N[1]);
      double tmpE=(g.N[0]+g.N[1]==0?0.0:double(g.E)/(g.N[0]+g.N[1]));
      float pdens=clamp(0.0,1.0,(log(readPopDensR(i,j)+0.5)-log(0.5))/(log(rmaxpop+0.5)-log(0.5))*0.6+0.4);
//      float pdens=clamp(0.0,1.0,readPopDensR(i+rxmin,rymax-j)/maxpdens*0.6+0.4);
//      float pdens=clamp(0.0,1.0,readPopDens(i*(psShape->dfXMax-psShape->dfXMin)/100.0+psShape->dfXMin,j*(psShape->dfYMax-psShape->dfYMin)/100.0+psShape->dfYMin)/maxpdens);
//      uint32_t color=tricolor(evector3(1.0,1.0,1.0)*pdens,evector3(0.3,0.8,0.3),evector3(1.0,0.0,0.0),tmpE,clamp(0.0,1.0,(log(tmpI+0.001)-log(0.001))/(-log(0.001))));
//      colorRect(frameRaw,sx,sy,sxn,syn,color,width);
      color=tricolor(evector3(1.0,1.0,1.0)*pdens,evector3(0.3,0.8,0.3),evector3(0.0,0.0,1.0),tmpE,clamp(0.0,1.0,(log(tmpI+0.001)-log(0.001))/(-log(0.001))));
      if (gridMask[j*st.spGridW+i]==0)
        color=evector3(1.0,1.0,1.0);
      cairo_move_to(cr,sx,sy);
      cairo_line_to(cr,sxn,sy);
      cairo_line_to(cr,sxn,syn);
      cairo_line_to(cr,sx,syn);
      cairo_close_path(cr);
      cairo_set_source_rgb(cr, color.x, color.y, color.z);
      cairo_fill(cr);
    }
  }
*/

  // tell cairo we have changed the image contents
  cairo_surface_mark_dirty(surface);

  for (int i=0; i<st.shapes.size(); ++i)
    cairoDrawShape(cr,st.shapes[i],xpos,ypos,scale,true,height,true);
  cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
  cairo_set_line_width(cr, 1);
  cairo_stroke (cr);

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_move_to (cr, 10.0, 50.0);
  cairo_show_text (cr, daystr);
  cairo_move_to (cr, 10.0, 100.0);
  cairo_show_text (cr, estr().sprintf("All cases: %.0lf",double(st.allCases))._str);
  cairo_move_to (cr, 10.0, 130.0);
  cairo_show_text (cr, estr().sprintf("New cases: %.0lf",double(newCases))._str);
  cairo_move_to (cr, 10.0, 160.0);
  cairo_show_text (cr, estr().sprintf("in ICU: %.0lf",double(st.allICU))._str);
  cairo_move_to (cr, 10.0, 190.0);
  cairo_show_text (cr, estr().sprintf("in nonICU: %.0lf",double(st.allNonICU))._str);
  cairo_move_to (cr, 10.0, 220.0);
  cairo_show_text (cr, estr().sprintf("Fatalities: %.0lf",double(st.allD))._str);
  cairo_move_to (cr, 10.0, vheight-32.0);
  cairo_show_text (cr, estr().sprintf("Exposed: %.1lf%%",double(st.allE)*100.0/(st.N[0]+st.N[1]))._str);
  if (mitigation<1.0){
    cairo_move_to (cr, 10.0, vheight-64.0);
    cairo_show_text (cr, estr().sprintf("Mitigation: %.2lf",mitigation)._str);
  }

/*
  cairo_set_source_rgb (cr, 1.0, 0.0, 1.0);
  cairo_rectangle (cr, 80.0, 60.0, 120.0, 80.0);
  cairo_fill (cr);
*/

  // finish any cairo drawing operations
  cairo_surface_flush(surface);

  videoPushFrame(frameRaw);
}

void initCairo(uint8_t *frameRaw,int width,int height)
{
//  cairo_surface_t *surface=cairo_image_surface_create(CAIRO_FORMAT_RGB32, width, height);
  surface=cairo_image_surface_create_for_data(frameRaw,CAIRO_FORMAT_ARGB32,width,height,width*4);
  int cstatus=cairo_surface_status(surface);
  cerr << "# stride RGB: " << cairo_format_stride_for_width (CAIRO_FORMAT_RGB24,width) << endl;
  cerr << "# stride ARGB: " << cairo_format_stride_for_width (CAIRO_FORMAT_ARGB32,width) << endl;
  ldieif(cstatus!=CAIRO_STATUS_SUCCESS,"error creating surface: "+estr(cstatus));
  cr=cairo_create(surface); 
  cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
//  cairo_set_antialias(cr,CAIRO_ANTIALIAS_FAST);

  cairo_select_font_face (cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size (cr, 32.0);
  cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);
}

#include <byteswap.h>

#ifndef TIFFTAG_GDAL_NODATA
#define TIFFTAG_GDAL_NODATA 42113
#endif

static const TIFFFieldInfo xgdaltiffFieldInfo[] = {
  /* XXX Insert Your tags here */
    { TIFFTAG_GDAL_NODATA,	-1,-1, TIFF_ASCII,	FIELD_CUSTOM,
      1,	0,	"GDAL nodata" }
};

#define	TIFF_N(a)	(sizeof (a) / sizeof (a[0]))

static TIFFExtendProc _ParentExtender = NULL;

static void _XTIFFDefaultDirectory(TIFF *tif)
{
  /* Install the extended Tag field info */
  TIFFMergeFieldInfo(tif, xgdaltiffFieldInfo, TIFF_N(xgdaltiffFieldInfo));

  /* Since an XTIFF client module may have overridden
   * the default directory method, we call it now to
   * allow it to set up the rest of its own methods.
   */

  if (_ParentExtender) 
      (*_ParentExtender)(tif);
}

static void _XTIFFInitialize(void)
{
  static int first_time=1;

  if (! first_time) return; /* Been there. Done that. */
  first_time = 0;

  /* Grab the inherited method and install */
  _ParentExtender = TIFFSetTagExtender(_XTIFFDefaultDirectory);
}



int loadPopDens(const estr& fname){
  _XTIFFInitialize(); // needed to setup GDAL novalue TAG info for TIFF, should just use the GDAL library for this

  TIFF *tif=0x00;
  tif=XTIFFOpen(fname._str,"r");
  if (!tif) return(-1);

  /* Install the extended Tag field info */
  TIFFMergeFieldInfo(tif, xgdaltiffFieldInfo, TIFF_N(xgdaltiffFieldInfo));
	
  gtif = GTIFNew(tif);
  if (!gtif) {
    fprintf(stderr,"failed in GTIFNew\n");
    return(-1);
  }

//  GTIFPrint(gtif,0,0);

  uint32_t ysize=0; //xsize declared above
  uint32_t xsize=0;
  TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &xsize );
  TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &ysize );
  uint32_t spp=0,sf=0,depth=0,bps=0,rowsperstrip=0;
  TIFFGetField( tif, TIFFTAG_BITSPERSAMPLE, &bps );
  TIFFGetField( tif, TIFFTAG_IMAGEDEPTH, &depth );
  TIFFGetField( tif, TIFFTAG_SAMPLESPERPIXEL, &spp );
  TIFFGetField( tif, TIFFTAG_SAMPLEFORMAT, &sf );
  TIFFGetField( tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

  char *tnodataval=0x00;
  TIFFGetField( tif, TIFFTAG_GDAL_NODATA, &tnodataval);
  nodataval=estr(tnodataval).f();
  
  cerr << "xsize: " << xsize << " ysize: " << ysize << " spp: " << spp << " sf: " << sf << " depth: " << depth << " bps: " << bps << " rowsperstrip: " << rowsperstrip << endl;

  GTIFDefn defn;
  if (GTIFGetDefn(gtif, &defn)){
//    printf( "\n" );
//    GTIFPrintDefnEx( gtif, &defn, stdout );

//    printf( "\n" );
//    printf( "PROJ.4 Definition: %s\n", GTIFGetProj4Defn(&defn));
  }

//  TIFFSetErrorHandler(TIFFError);

//	if (!TIFFReadRGBAImage(tif, xsize, ysize, raster, 0)) { lerror("reading TIFF data"); exit(-1); return(-1); }


//  printf( "\nCorner Coordinates:\n" );

  unsigned short raster_type = RasterPixelIsArea;
  GTIFKeyGetSHORT(gtif, GTRasterTypeGeoKey, &raster_type, 0, 1);

  double xmin = (raster_type == RasterPixelIsArea) ? 0.0 : -0.5;
  double ymin = xmin;
  double ymax = ymin + ysize;
  double xmax = xmin + xsize;

  ldieif(defn.Model != ModelTypeGeographic,"model used in GeoTIFF not implemented");
//  cout << "defn->model is ModelTypeGeographic? " << (defn.Model == ModelTypeGeographic) << endl; 

  double c_ul_x=xmin,c_ul_y=ymin;
  double c_bl_x=xmin,c_bl_y=ymax;
  double c_ur_x=xmax,c_ur_y=ymin;
  double c_br_x=xmax,c_br_y=ymax;
  if ( !GTIFImageToPCS( gtif, &c_ul_x, &c_ul_y ) )
    return -1;
  GTIFImageToPCS( gtif, &c_bl_x, &c_bl_y );
  GTIFImageToPCS( gtif, &c_ur_x, &c_ur_y );
  GTIFImageToPCS( gtif, &c_br_x, &c_br_y );

  cerr << "upper left: " << c_ul_x << " " << c_ul_y << endl;
  cerr << "upper right: " << c_ur_x << " " << c_ur_y << endl;
  cerr << "bottom left: " << c_bl_x << " " << c_bl_y << endl;
  cerr << "bottom right: " << c_br_x << " " << c_br_y << endl;

  uint32_t rxmin=0;
  uint32_t rxmax=0;
  uint32_t rymin=0;
  uint32_t rymax=0;

  double sX=xlonMin;
  double sY=ylatMin;
  cerr << "top left of shape: " << sX << " " << sY << endl;
  GTIFPCSToImage( gtif, &sX,&sY );
  cerr << "top left of shape: " << sX << " " << sY << endl;
//  cout << "value: " << raster[uint32_t(sY)*xsize + uint32_t(sX)] << " nodataval: " << nodataval << endl;
  rxmin=sX;
  rymin=sY;

  sX=xlonMax;
  sY=ylatMax;
  cerr << "bottom right of shape: " << sX << " " << sY << endl;
  GTIFPCSToImage( gtif, &sX,&sY );
  cerr << "bottom right of shape: " << sX << " " << sY << endl;
//  cout << "value: " << raster[uint32_t(sY)*xsize + uint32_t(sX)] << " nodataval: " << nodataval << endl;
  rxmax=sX;
  rymax=sY;

  if (rxmin>rxmax) swap(rxmin,rxmax);
  if (rymin>rymax) swap(rymin,rymax);
  cerr << "sxmin: " << rxmin << " sxmax: " << rxmax <<endl;
  cerr << "symin: " << rymin << " symax: " << rymax <<endl;
  rwidth=rxmax-rxmin;
  rheight=rymax-rymin;
  cerr << "rwidth: " << rwidth << " rheight: " << rheight <<endl;
 
	int npixels = rwidth * rheight;
	float *raster = (float*) _TIFFmalloc(npixels * sizeof (float));
	if (raster == NULL) return(-1);

  uint32_t config=0u;
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
  ldieif(config!=PLANARCONFIG_CONTIG,"separate planes in TIFF file not implemented");

  tstrip_t strip;
  int bytesread;
  ldieif(rowsperstrip>1,"GeoTIFF: more than one row per strip not implemented");

  uint32_t stripsize = TIFFStripSize(tif);
	float *sbuf = (float*) _TIFFmalloc(stripsize);
  uint32_t rp=0;
//  buf = _TIFFmalloc(stripsize);
  for (uint32_t strip=rymin, row=rymin; row<ysize && row<rymax ; ++strip,row+=rowsperstrip) {
//    if (rp+stripsize/sizeof(uint32_t) > npixels) ldie("TIFF data larger than expected: "+estr(row)+" "+estr(stripsize)+" "+estr(rp)+" "+estr(strip)+" "+estr(npixels)); 
    bytesread=TIFFReadEncodedStrip(tif, strip, sbuf, stripsize);
    ldieif(bytesread==-1,"Error reading GeoTIFF");
//    if (row>=rymin)
    memcpy(&raster[(row-rymin)*rwidth],&sbuf[rxmin],rwidth*4);
  }
  _TIFFfree(sbuf);

/*
  uint32_t nbytes;
  uint32_t rp=0u;
  uint32_t scansize=TIFFScanlineSize(tif);
  cout << "Scansize: " << scansize << endl;
 	for (uint32_t row=0; row<ysize; ++row){
 	  if (TIFFReadScanline(tif,&raster[rp],row) == -1){ lerror("reading TIFF file"); exit(-1); }
    rp+=scansize/sizeof(float);
// 	buf = _TIFFmalloc(TIFFScanlineSize(tif));
  }
*/
 
  popCounts.init(rheight*rwidth);
  for (int i=0; i<rheight*rwidth; ++i){
    if (raster[i]==nodataval) raster[i]=0.0;
    popCounts[i]=raster[i];
    ldieif(popCounts[i]<0,"negative population counts on grid?");
  }

  return(0);
}

int getGrid(double lat,double lon)
{
  if (ylatMin>lat || ylatMax<lat || xlonMin>lon || xlonMax<lon) return(-1);
  return((rheight-uint32_t((lat-ylatMin)*rheight/(ylatMax-ylatMin)))*rwidth+uint32_t((lon-xlonMin)*rwidth/(xlonMax-xlonMin)));
}

void adjustPopCounts(ssimstate& st,int popsize)
{
  cairo_surface_t *surfacemask=0x00;
  cairo_t *crmask=0x00;

  uint8_t *shapeMask=0x00;

  ldieif(popCounts.size()==0,"popCounts not initialized yet?");

  int shapeRow=cairo_format_stride_for_width(CAIRO_FORMAT_A8,rwidth);
  shapeMask=new uint8_t[rheight*shapeRow];
  memset(shapeMask,0,rheight*shapeRow);

  surfacemask=cairo_image_surface_create_for_data(shapeMask,CAIRO_FORMAT_A8,rwidth,rheight,shapeRow);
  int cstatus=cairo_surface_status(surfacemask);
  cerr << "# stride A8 width: " << rwidth << " stride: " << cairo_format_stride_for_width (CAIRO_FORMAT_A8,rwidth) << endl;
  ldieif(cstatus!=CAIRO_STATUS_SUCCESS,"error creating surface: "+estr(cstatus));
  crmask=cairo_create(surfacemask); 

  cairo_surface_t *colorSurf=0x00;
  int colorRow=cairo_format_stride_for_width(CAIRO_FORMAT_A8,1);
  uint8_t *colorData=new uint8_t[1*colorRow];
  colorData[0]=1u;
  colorSurf=cairo_image_surface_create_for_data(colorData,CAIRO_FORMAT_A8,1,1,colorRow);
//  cairo_pattern_t *colorPattern=cairo_pattern_create_for_surface(colorSurf);

  for (int i=0; i<st.shapes.size(); ++i){
    // color code each region
    colorData[0]=i+1;
    cairo_surface_mark_dirty(colorSurf);
    
    cairo_set_source_surface(crmask,colorSurf,0.0,0.0);
    cairo_pattern_set_extend(cairo_get_source(crmask),CAIRO_EXTEND_REPEAT); // repeat the pattern
//    cairo_set_source_rgb(crmask,(i+1+0.5)/255.0,(i+1+0.5)/255.0,(i+1+0.5)/255.0);
    cairoDrawShape(crmask,st.shapes[i],0.0,0.0,float(rwidth)/(xlonMax-xlonMin),false,rheight,true);
    cairo_fill(crmask);
  }

  cairo_surface_flush(surfacemask);
  cairo_destroy(crmask);
  cairo_surface_destroy(surfacemask);
//  cairo_pattern_destroy(colorPattern);
  cairo_surface_destroy(colorSurf);
  delete[] colorData;
 
  
  st.gridMask.init(rwidth*rheight);
  for (int i=0; i<rheight; ++i){
    for (int j=0; j<rwidth; ++j)
      st.gridMask[i*rwidth+j]=shapeMask[i*shapeRow+j];
  }
  delete[] shapeMask;

  rpopsize=0;
  rmaxpop=0;
  imaxpop=0;
  int totalpopsize=0;
  for (int i=0; i<popCounts.size(); ++i){
    totalpopsize+=popCounts[i];
    if (st.gridMask[i]==0) { popCounts[i]=0; continue; }
    if (popCounts[i]>rmaxpop) { imaxpop=i;  rmaxpop=popCounts[i]; }
    rpopsize+=popCounts[i];
  }
  cerr << "rpopsize: " << rpopsize << " totalpopsize: " << totalpopsize << endl;
  cerr << "rmaxpop: " << popCounts[imaxpop] << endl;


  

  cerr << "Adjusting raster population counts to: " << popsize << endl;
  while (fabs(rpopsize-popsize)>0.01*popCounts[imaxpop]){
    // adjust population counts to match current country population census
    long numfrac=0;
    for (int i=0; i<popCounts.size(); ++i){
      long tmpf=popCounts[i];
      tmpf=(tmpf*popsize+numfrac)/rpopsize;
      numfrac=(long(popCounts[i])*popsize+numfrac)%rpopsize;
      popCounts[i]=tmpf;
//    frac=(tmpf*double(popsize)/rpopsize+frac-floor(tmpf*double(popsize)/rpopsize+frac));
    }
//  rmaxpop=rmaxpop*double(popsize)/rpopsize;
    cerr << "rpopsize: " << rpopsize << " " << double(rpopsize)/popsize << endl;

    rpopsize=0;
    for (int i=0; i<popCounts.size(); ++i)
      rpopsize+=popCounts[i];
  }
  cerr << "rpopsize after initial adjustment: " << rpopsize << endl;
  cerr << "rmaxpop after initial adjustment: " << popCounts[imaxpop] << endl;
  cerr << "Final adjustment to grid with highest count:" << endl;
  popCounts[imaxpop]+=popsize-rpopsize;
  rmaxpop=popCounts[imaxpop];
  cerr << "adjusted rmaxpop: " << rmaxpop << endl;

  for (int i=0; i<st.gridMask.size(); ++i){
    int gm=st.gridMask[i];
    ldieif(gm<0 || gm>=st.shapeCounts.size(),"gridMask value out of bounds: "+estr(gm)+" "+st.shapeCounts.size());
    st.shapeCounts[gm].allN+=popCounts[i];
  }

  cerr << "# shapeCounts" << endl;
  for (int i=0; i<st.shapeCounts.size(); ++i){
    cerr << i << " " << st.shapeCounts[i].allN << endl;
  }



//  exit(0);
}

void loadDBF(const estr& fname,const estrarray& selection,estrarrayof<int>& shplist)
{
  DBFHandle hDBF = DBFOpen( fname._str, "rb" );
//  for (int i=0; i<selection.size(); ++i)   // this way does not allow for multiple matching
//    shplist.add(selection.values(i),-1);

  if (hDBF == NULL ){
    printf( "DBFOpen(%s,\"r\") failed.\n", fname._str );
    exit( 2 );
  }
    
  if ( DBFGetFieldCount(hDBF) == 0 ){
    printf( "There are no fields in this table!\n" );
    exit( 3 );
  }

  char szTitle[255],szFormat[255];

  estrarray fieldNames;

  for (int i = 0; i < DBFGetFieldCount(hDBF); ++i){
    DBFFieldType	eType;
    const char	 	*pszTypeName;
    char chNativeType;

    chNativeType = DBFGetNativeFieldType( hDBF, i );
    int nWidth,nDecimals;

    eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );
    if( eType == FTString )
        pszTypeName = "String";
    else if( eType == FTInteger )
        pszTypeName = "Integer";
    else if( eType == FTDouble )
        pszTypeName = "Double";
    else if( eType == FTInvalid )
        pszTypeName = "Invalid";
    fieldNames.add(szTitle);
//    printf( "Field %d: Type=%c/%s, Title=`%s', Width=%d, Decimals=%d\n",
//            i, chNativeType, pszTypeName, szTitle, nWidth, nDecimals );
  }
  cerr << "# shape dbf file fields: " << fieldNames << endl;

/*
  for( int i = 0; i < DBFGetFieldCount(hDBF); i++ ){
    DBFFieldType	eType;
    int nWidth,nDecimals;

    eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );
  }
  printf( "\n" );
*/

  estrarrayof<evar> valarr;

  bool is_hit;
  for ( int iRecord = 0; iRecord < DBFGetRecordCount(hDBF); iRecord++ ) {
    valarr.clear();
    for ( int i = 0; i < DBFGetFieldCount(hDBF); i++ ){
      DBFFieldType eType;
      int nWidth,nDecimals;
            
      eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );

      if( DBFIsAttributeNULL( hDBF, iRecord, i ) ) {
        valarr.add(szTitle,evar());
      } else {
        switch( eType ){
          case FTString:
            valarr.add(szTitle,DBFReadStringAttribute( hDBF, iRecord, i ));
            break;
          case FTInteger:
            valarr.add(szTitle,DBFReadIntegerAttribute( hDBF, iRecord, i ));
            break;
          case FTDouble:
            valarr.add(szTitle,DBFReadDoubleAttribute( hDBF, iRecord, i ));
            break;
          default:
            break;
        }
      }
/*
      if (selection.findkey(szTitle)==-1) continue; // skip fields not in selection list
      if( !DBFIsAttributeNULL( hDBF, iRecord, i ) && eType==FTString) {
        const char *fvalue=DBFReadStringAttribute( hDBF, iRecord, i );
        int iv=selection.find(fvalue);
        if (iv!=-1){
          is_hit=true;
          shplist.add(fvalue,iRecord);
//          sprintf( szFormat, "%i %%-%ds\n", iRecord, nWidth );
//          printf( szFormat, fvalue );
        }
      }
*/
    }
    bool is_hit=false,is_excluded=false;
    for (int j=0; j<selection.size(); ++j){
      estr fname=selection.keys(j);
      bool exclude=false;
      if (selection.keys(j)[0]=='-'){
        exclude=true;
        fname=selection.keys(j).substr(1);
      }
      if (valarr[fname]==selection.values(j)){
        if (!exclude) is_hit=true;
        else is_excluded=true;
      }
    }
    if (is_hit && !is_excluded){
      shplist.add(iRecord);
      cerr << valarr << endl;
    }
  }
/*
      if( DBFIsAttributeNULL( hDBF, iRecord, i ) ) {
        if( eType == FTString )
          sprintf( szFormat, "%%-%ds", nWidth );
        else
          sprintf( szFormat, "%%%ds", nWidth );

        printf( szFormat, "(NULL)" );
      } else {
        switch( eType ){
          case FTString:
            sprintf( szFormat, "%%-%ds", nWidth );
            printf( szFormat, 
                    DBFReadStringAttribute( hDBF, iRecord, i ) );
            break;
            
          case FTInteger:
            sprintf( szFormat, "%%%dd", nWidth );
            printf( szFormat, 
                    DBFReadIntegerAttribute( hDBF, iRecord, i ) );
            break;
            
          case FTDouble:
            sprintf( szFormat, "%%%d.%dlf", nWidth, nDecimals );
            printf( szFormat, 
                    DBFReadDoubleAttribute( hDBF, iRecord, i ) );
            break;
            
          default:
            break;
        }
      }
      printf("  ");
*/

  DBFClose( hDBF );
}

void loadShape(ssimstate& st,const estr& fname,const estrarray& selection,double lonLimit=-1000.0){

  estrarrayof<int> shplist;
  loadDBF(fname.substr(0,-5)+".dbf",selection,shplist);

  SHPHandle hSHP;
  int nShapeType, nEntities;
  double adfMinBound[4],adfMaxBound[4];

  hSHP = SHPOpen( fname._str, "rb" );
  SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

  int nPrecision = 15;
//  printf( "Shapefile Type: %s   # of Shapes: %d\n\n",
//            SHPTypeName( nShapeType ), nEntities );
    
/*
  printf( "File Bounds: (%.*g,%.*g,%.*g,%.*g)\n"
            "         to  (%.*g,%.*g,%.*g,%.*g)\n",
            nPrecision, adfMinBound[0], 
            nPrecision, adfMinBound[1], 
            nPrecision, adfMinBound[2], 
            nPrecision, adfMinBound[3], 
            nPrecision, adfMaxBound[0], 
            nPrecision, adfMaxBound[1], 
            nPrecision, adfMaxBound[2], 
            nPrecision, adfMaxBound[3] ); 
*/

  st.shapeCounts.add(scounts()); // outside the shapes

  for (int i=0; i<nEntities; ++i){
    if (nEntities>1 && shplist.find(i)==-1) continue; // if there is only one shape, take it, otherwise use the list

    int j;
    SHPObject *tmpShape;

    tmpShape = SHPReadObject(hSHP, i);

    ldieif (tmpShape == NULL,"Unable to read shape "+estr(i)+" in file: "+fname);
    xlonMin=MIN(xlonMin,tmpShape->dfXMin);
    xlonMax=MAX(xlonMax,tmpShape->dfXMax);
    ylatMin=MIN(ylatMin,tmpShape->dfYMin);
    ylatMax=MAX(ylatMax,tmpShape->dfYMax);
    st.shapes.add(tmpShape);
    st.shapeCounts.add(scounts());

/*
      printf( "\nShape:%d (%s)  nVertices=%d, nParts=%d\n"
                  "  Bounds:(%.*g,%.*g, %.*g)\n"
                  "      to (%.*g,%.*g, %.*g)\n",
                    0, SHPTypeName(psShape->nSHPType),
                    psShape->nVertices, psShape->nParts,
                    nPrecision, psShape->dfXMin,
                    nPrecision, psShape->dfYMin,
                    nPrecision, psShape->dfZMin,
                    nPrecision, psShape->dfXMax,
                    nPrecision, psShape->dfYMax,
                    nPrecision, psShape->dfZMax ); 
*/
    
//    SHPDestroyObject(psShape);
  }

/*
  int tj=0;
  int tiPart=1;

  eintarray newParts;
  newParts.add(0);
  newParts.add(0);

  // cut shape out of limits
  for (int j=0,iPart=1; lonLimit>-1000.0 && j<psShape->nVertices; ++j){
    if (iPart < psShape->nParts && psShape->panPartStart[iPart]==j){
      ++iPart;
      newParts.add(tj);
    }
    if (psShape->padfX[j]<lonLimit) continue;
    if (j!=tj){
      psShape->padfX[tj]=psShape->padfX[j];
      psShape->padfY[tj]=psShape->padfY[j];
    }
    if (tj==0 || psShape->padfX[tj]<psShape->dfXMin) psShape->dfXMin=psShape->padfX[tj];
    if (tj==0 || psShape->padfX[tj]>psShape->dfXMax) psShape->dfXMax=psShape->padfX[tj];
    if (tj==0 || psShape->padfY[tj]<psShape->dfYMin) psShape->dfYMin=psShape->padfY[tj];
    if (tj==0 || psShape->padfY[tj]>psShape->dfYMax) psShape->dfYMax=psShape->padfY[tj];
    ++tj;
    newParts[newParts.size()-1]=tj;
  }
  if (lonLimit>-1000.0){
//    for (int i=0; i<psShape->nParts; ++i)
//      cout << "i: "<< i << " " << psShape->panPartStart[i] << endl;
    psShape->nParts=1;
    for (int i=1,ti=1; i<newParts.size(); ++i){
      if (newParts[i]>newParts[i-1]){
        psShape->panPartStart[ti]=newParts[i];
        ++ti;
        ++psShape->nParts;
      }
    }
//    for (int i=0; i<psShape->nParts; ++i)
//      cerr << "i: "<< i << " " << psShape->panPartStart[i] << endl;
    cerr << "Reduced shape with lonLimit: " << psShape->nVertices << " to " << tj << endl;
    psShape->nVertices=tj;
*/
/*
      printf( "\nNew Shape:%d (%s)  nVertices=%d, nParts=%d\n"
                  "  Bounds:(%.*g,%.*g, %.*g)\n"
                  "      to (%.*g,%.*g, %.*g)\n",
                    0, SHPTypeName(psShape->nSHPType),
                    psShape->nVertices, psShape->nParts,
                    nPrecision, psShape->dfXMin,
                    nPrecision, psShape->dfYMin,
                    nPrecision, psShape->dfZMin,
                    nPrecision, psShape->dfXMax,
                    nPrecision, psShape->dfYMax,
                    nPrecision, psShape->dfZMax ); 
  }
*/
    


  // find projected min and max longitudes
  lonRef=(xlonMax+xlonMin)/2.0;
  for (int i=0; i<st.shapes.size(); ++i){
    for (int j=0; j<st.shapes[i]->nVertices; ++j){
      double tmpLon=proj(st.shapes[i]->padfX[j],st.shapes[i]->padfY[j],lonRef);
      if (i==0 && j==0 || tmpLon<plonMin) plonMin=tmpLon;
      if (i==0 && j==0 || tmpLon>plonMax) plonMax=tmpLon;
    }
  }

/*
  double cx=(psShape->dfXMax-psShape->dfXMin)/2.0;
  for (int j=0; j<psShape->nVertices; ++j)
    psShape->padfX[j] = cos(M_PI*psShape->padfY[j]/180.0)*(psShape->padfX[j]-cx);

  psShape->dfXMax=psShape->padfX[0];
  psShape->dfXMin=psShape->padfX[0];
  psShape->dfYMax=psShape->padfY[0];
  psShape->dfYMin=psShape->padfY[0];
  for (int j=1; j<psShape->nVertices; ++j){
    psShape->dfXMax = MAX(psShape->dfXMax,psShape->padfX[j]);
    psShape->dfXMin = MIN(psShape->dfXMin,psShape->padfX[j]);
    psShape->dfYMax = MAX(psShape->dfYMax,psShape->padfY[j]);
    psShape->dfYMin = MIN(psShape->dfYMin,psShape->padfY[j]);
  }
*/

  SHPClose(hSHP);
}

void blurKernel(int i,int j,edoublearray& barr,edoublearray& arr,int w,edoublearray& k,int ks){
  for (int ty=0; ty<ks; ++ty){
    for (int tx=0; tx<ks; ++tx)
      barr[(j+ty-ks/2)*w+i+tx-ks/2]+=k[ty*ks+tx]*arr[j*w+i];
  }
}

double travelKernel(ssimstate& st,int gx,int gy)
{
  int tymin=MAX(0,gy-st.travelKernelSize/2)-gy+st.travelKernelSize/2;
  int tymax=MIN(st.spGridH,gy+st.travelKernelSize/2+1)-gy+st.travelKernelSize/2;
  int txmin=MAX(0,gx-st.travelKernelSize/2)-gx+st.travelKernelSize/2;
  int txmax=MIN(st.spGridW,gx+st.travelKernelSize/2+1)-gx+st.travelKernelSize/2;

  double prob=0.0;
  double norm=0.0;
  for (int ty=tymin; ty<tymax; ++ty){
    for (int tx=txmin; tx<txmax; ++tx){
      norm+=st.travelKernel[ty*st.travelKernelSize+tx];
      prob+=st.travelKernel[ty*st.travelKernelSize+tx]*st.localIprob[(gy-st.travelKernelSize/2+ty)*st.spGridW+gx-st.travelKernelSize/2+tx];
    }
  }
  return(prob/norm);
}


void blur(edoublearray& barr,edoublearray& arr,int w,int h,edoublearray& k,int ks)
{
  for (int i=0; i<barr.size(); ++i)
    barr[i]=0.0;
  for (int j=ks/2; j<h-ks/2; ++j)
    for (int i=ks/2; i<w-ks/2; ++i)
      blurKernel(i,j,barr,arr,w,k,ks);
}


void threadLocalInfection(ssimstate& st,sthreadState& ths,int i,int n)
{
  for (int ti=st.spGrid.size()*i/n; ti<st.spGrid.size()*(i+1)/n && ti<st.spGrid.size(); ++ti){
    int gi=st.gridShuffle[ti];
    sgrid& g(st.spGrid[gi]);
    if (g.N[0]+g.N[1]==0) { st.localIprob[gi]=0.0; continue; }
    st.localIprob[gi]=st.finter*st.R0day*(0.5*g.Ia[0]+g.Ip[0]+g.Is[0] + (MIN(st.smr,st.soldr)/st.smr)*(0.5*g.Ia[1]+g.Ip[1]+g.Is[1]) )/(g.N[0]+(MIN(st.smr,st.soldr)/st.smr)*g.N[1]);
  }
}

void threadUpdateGrid(ssimstate& st,sthreadState& ths,int i,int n)
{
  // todo: include this in thread state
  edoublearray mp;
  euintarray counts;

  double globalIprob=st.finter*st.fglobal*st.R0day*(0.5*st.Ia[0]+st.Ip[0]+st.Is[0] + (MIN(st.smr,st.soldr)/st.smr)*(0.5*st.Ia[1]+st.Ip[1]+st.Is[1]) )/(st.N[0]+(MIN(st.smr,st.soldr)/st.smr)*st.N[1]);

  for (int ti=st.spGrid.size()*i/n; ti<st.spGrid.size()*(i+1)/n && ti<st.spGrid.size(); ++ti){
    int gi=st.gridShuffle[ti];
    sgrid& g(st.spGrid[gi]);
    if (g.N[0]+g.N[1]==0) { continue; } // skip empty grid

    int gx=gi%st.spGridW;
    int gy=gi/st.spGridW;

    int giup=((st.spGridH+gy-1)%st.spGridH)*st.spGridW+gx;
    int gidn=((gy+1)%st.spGridH)*st.spGridW+gx;
    int gilt=gy*st.spGridW+(st.spGridW+gx-1)%st.spGridW;
    int girt=gy*st.spGridW+(gx+1)%st.spGridW;

    // TODO: handle adjacent empty grids or big differences in population sizes. How should they affect the probability of infection of an adjacent grid?
//    double tmpLocalIprob=localIprobGaussian[gi]*flocal; // + localIprob[gilt]*(1.0-flocal);
//    double tmpLocalIprob=localIprob[gi]*flocal; // + localIprob[gilt]*(1.0-flocal);
//    double tmpLocalIprob=st.localIprob[gi]*0.6 + st.localIprob[gilt]*0.1 + st.localIprob[girt]*0.1 + st.localIprob[giup]*0.1 + st.localIprob[gidn]*0.1;
    double tmpLocalIprob=travelKernel(st,gx,gy);

    for (int hg=0; hg<ngroups; ++hg){
      for (int hhsize=1; hhsize<st.householdSize; ++hhsize){ // TODO: should move levels to outter loop for efficiency, inner loop should be longest sequential array of data
        for (int hhexp=hhsize-1; hhexp>=0; --hhexp){ // NOTE: have to do decreasing here to avoid infecting twice the same house in the same step

          uint8_t hl=hhsize*(hhsize+1)/2+(hhsize-hhexp);
          unsigned int hlcount=st.hhLevelsBegin[gi*st.grow + hg*st.arow + hl+1]-st.hhLevelsBegin[gi*st.grow + hg*st.arow +hl];
          if (hlcount==0) continue; // no households in level

          int hS=hhsize-hhexp; // number of susceptible individuals per household
          int nS=hlcount*hS; // total number of susceptible individuals in these households

          // The intra household rate below is an approximation needed to avoid computing the exact probability of new infection per household which depends on the specific combination of Ia,Is and Ip
          double intraIprob=st.fintra*st.R0day*(0.5*g.hhIa[hg][hl]+g.hhIp[hg][hl]+g.hhIs[hg][hl])/nS;
  
          double p=(MIN(st.smr,(hg>=1?st.soldr:1.0))*(globalIprob + tmpLocalIprob) + intraIprob)*st.tstep; // probability of infection per susceptible person in these households
  
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
          multinomial(hlcount,mp,counts,ths.r);
         
  //        int newE=binomial(nS,,rnd)/(hhsize-hhexp);
  //        if (newE>st.hhLevels[hhstlevel].size()) newE=st.hhLevels[hhstlevel].size(); // cap maximum number of new household exposures
 
          // each thread has its own event queue since every grid is independent of each other this ensures there are no race conditions over data access
          ths.evqueue.add(st,ths,gi,hg,hhsize,hhexp,counts,st.dE,ths.r);
//          st.evqueue.add(st,gi,hg,hhsize,hhexp,counts,st.dE,rnd);
        }
      }
    }
  }
  deque<sevent> &nextEvents(ths.evqueue.step());
  processEvents(st,ths,nextEvents,ths.r);
}

/*
void threadUpdateEvents(ssimstate& st,sthreadState& ths,int i,int n,ernd& r)
{
  deque<sevent> &nextEvents(ths.evqueue.step());
  processEvents(st,ths,nextEvents,r);
}
*/

void threadSeedInfections(ssimstate& st,sthreadState& ths,int i,int n)
{
//  if (st.seedGrid < st.spGrid.size()*i/n || st.seedGrid >= st.spGrid.size()*(i+1)/n) return;
  st.mutex.lock();
//  cerr << "i: " << i << " n: " << n << " " << st.spGrid.size()*i/n << " " << st.spGrid.size()*(i+1)/n << " " << st.spGrid.size() << endl;
  st.mutex.unlock();


  eintarray seedlist;
  for (int ti=st.spGrid.size()*i/n; ti<st.spGrid.size()*(i+1)/n && ti<st.spGrid.size(); ++ti){
    int gi=st.gridShuffle[ti];
    int ik=st.iseedArr.findkey(gi);
    if (ik!=-1) seedlist.add(ik);
  }
  if (seedlist.size()==0) return;

//  cerr << "# seeding initial infected: " << st.fseed << endl;

  // seed infections
  edoublearray mp;
  euintarray counts;

  for (int i=0; i<seedlist.size(); ++i){
    counts.init(3,0);
    counts[1]=st.iseedArr.values(seedlist[i]);
    ths.evqueue.add(st,ths,st.iseedArr.keys(seedlist[i]),0,2,0,counts,st.dE,ths.r);
  }
}

void threadCopyResetState(ssimstate &st,sthreadState& ths)
{
  st.allCases+=ths.allCases;
  ths.allCases=0;
  st.allE+=ths.allE;
  ths.allE=0;
  st.allICU+=ths.allICU; ths.allICU=0;
  st.allNonICU+=ths.allNonICU; ths.allNonICU=0;
  st.allD+=ths.allD; ths.allD=0;
  for (int i=0; i<ngroups; ++i){
    st.Ip[i]+=ths.Ip[i];
    st.Is[i]+=ths.Is[i];
    st.Ia[i]+=ths.Ia[i];
    ths.Ip[i]=0; ths.Is[i]=0; ths.Ia[i]=0;
  }
  addReset(st.shapeCounts,ths.shapeCounts);
  addReset(st.ageCounts,ths.ageCounts);
}

void threadRun(ssimstate &st)
{
//  sthreadState ths;
  int threadI=-1;
  threadFunc_t threadFunc=0x00;

  do {
    st.mutex.lock();
    
    --st.threadDone;
    // assign thread ID, needs to be fixed such that each thread handles always the same grid positions
//    if (threadI==-1){
//      threadI=st.threadDone;
//      cout << st.threadDone << endl;
//    }

    // copy thead state
    if (threadI>=0){
      sthreadState &ths(st.threadStates[threadI]);
      threadCopyResetState(st,ths);
    }

    st.doneSignal.signal();
    while (st.threadI==0) st.stateSignal.wait(st.mutex);

    threadFunc=st.threadFunc;
    --st.threadI;
    threadI=st.threadI;
//    cout << "Thread i: " << threadI << endl;
    st.mutex.unlock();

    if (threadFunc)
      (*threadFunc)(st,st.threadStates[threadI],threadI,st.nthreads);
  } while(threadFunc);
}

void actionDaemon()
{
}

estr sfile;

void bgThreadRun()
{
  startDaemon(sfile);
  getSystem().run();
}

int findPopGrid(ssimstate& st,double lat,double lon,int seedcount){
//  int seedGrid=getGrid(46.0,8.95);
  int seedGrid=getGrid(lat,lon);
//  if (st.seedGrid==-1) { cerr << "# coordinate not found chosing most populated area" << endl; st.seedGrid=imaxpop; }
  if (st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2]<seedcount){
//    cerr << "igrid: " << st.seedGrid << " levels empty, looking for another grid pos" << endl;
    int ix=seedGrid%st.spGridW;
    int ih=seedGrid/st.spGridW;
    for (int ti=0;ix<st.spGridW && ti<10 && st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2]<seedcount; ++ix,++seedGrid,++ti);
    cerr << "seedGrid: " << seedGrid << " " << st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2] << endl;
//    if (ix==st.spGridW) { cerr << "did not find any populated areas close to choice" << endl; seedGrid=imaxpop; }
    if (ix==st.spGridW || st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2]<seedcount) { cerr << "did not find any populated areas close to choice" << endl; seedGrid=-1; }
  }
  return(seedGrid);
}


int emain()
{
  bool printregions=false;
  epregister(printregions);

  int rseed=-1;
  epregister(rseed);

  ssimstate st;
  epregister2(st.R0,"r0");

  double cf=1.0;
  epregister2(cf,"cf"); // correction for ratio of Symptomatic/Asymptomatic, can be as high as 0.37 according to preliminary results from Streeck

  double df=1.0;
//  epregister2(df,"df"); // correction for unclassified covid deaths, preliminary results from increased number of fatalities compared to historical data show that covid deaths are being underreported by about 20%

  int agethres=50;
  epregister(agethres);

//  int avgage=0; // take most age of first individual (usually the oldest), or the avg age of the household for the household age group
//  epregister(avgage);

//  double finter=1.0;
//  double fintra=1.0;

/*
  double foldr=1.0; // isolation for households with older members
  int tostart=-1; // start isolation
  int toend=-1; // end isolation
  epregister(foldr);
  epregister(tostart);
  epregister(toend);

  double fmr=1.0; // mitigation factor
  int tmstart=-1; // start mitigation
  int tmend=-1; // end mitigation
  epregister(fmr);
  epregister(tmstart);
  epregister(tmend);
*/

  estr ovideo;
  epregister(ovideo);

  double tmax=300.0;
  epregister(tmax);
  epregister2(st.finter,"finter");
  epregister2(st.fintra,"fintra");
//  double fglobal=0.001;
  epregister2(st.fglobal,"fglobal");
//  double flocal=1.0;
  epregister2(st.flocal,"flocal");
  estr fpop="data/switzerland.agegroups";
  epregister(fpop);
//  estr fshape="data/popdensmaps/gadm36_CHE/gadm36_CHE_0.shp";
  estr fshape="data/popdensmaps/ref-nuts/NUTS_RG_01M_2016_4326_LEVL_2.shp";
 
  epregister(fshape);
//  int fseed=10;
//  epregister2(st.fseed,"fseed");

  estrarray iseed;
  epregister(iseed);


//  int nthreads=4;
  epregister2(st.nthreads,"nthreads");

  double ftravel=0.2;
  epregister(ftravel);

  estrarray events;
  epregister(events);

  estrarray oevents;
  epregister(oevents);

  double lonLimit=-1000.0;
  epregister(lonLimit);

  double fsinglehh=0.0;
  epregister(fsinglehh);
  double flargehh=0.0;
  epregister(flargehh);

//  bool icucontrol=false;
//  epregister(icucontrol);

  int casestrigger=-1;
  epregister(casestrigger);

  int icutrigger=-1;
  epregister(icutrigger);

  int oisotrigger=-1;
  epregister(oisotrigger);

  double oisorelease=-1;
  epregister(oisorelease);

  double foiso=0.05;
  epregister(foiso);


  double rthres=0.99;
  double rthres2=0.95;
  epregister(rthres);
  epregister(rthres2);


  int iculowlimit=-1;
  epregister(iculowlimit);
  int iculowlimit2=-1;
  epregister(iculowlimit2);
  double ficu=0.3;
  epregister(ficu);
  double ficu2=0.4;
  epregister(ficu2);

  int triggerdays=5;
  epregister(triggerdays);
  int triggercount=3;
  epregister(triggercount);

  double icuthres=0.0;
  epregister(icuthres);

  double fstep=0.1;
  epregister(fstep);


/*
  double ftrigger=1.0;
  epregister(ftrigger);
*/

//  daemonArgs(actionDaemon);
  epregister(sfile);
  int it;
  epregister(it);

  estrarray shpsel;
  shpsel.add("CNTR_CODE","CH");
  epregister(shpsel);

  eparseArgs();

  if (rseed!=-1){
    lwarn("setting seed. note: results are only reproducible when using the same number of threads");
    rnd.setSeed(rseed);
  }
  cout << "# random seed: " << rnd.seed << endl;
  cerr << "# random seed: " << rnd.seed << endl;



  st.rSA=st.rSA*cf; // correction for ratio of symptomatic/asymptomatic. Preliminary results from Hendrik Streeck show that this might be overestimated by 1.0/0.37 times

  ethreadFunc bgThread;
  if (sfile.len())
    bgThread.run(bgThreadRun);


  estrarray tmparr;
  earrayof<double,int> eventsArr;
  for (int i=0; i<events.size(); ++i){
    tmparr=events[i].explode(":");
    ldieif(tmparr.size()<2,"wrong format for event: <day>:<isolationfactor>,<day2>:<isolationfactor2>,...");
    eventsArr.add(tmparr[0].i(),tmparr[1].f());
  }

  earrayof<double,int> oeventsArr;
  for (int i=0; i<oevents.size(); ++i){
    tmparr=oevents[i].explode(":");
    ldieif(tmparr.size()<2,"wrong format for event: <day>:<isolationfactor>,<day2>:<isolationfactor2>,...");
    oeventsArr.add(tmparr[0].i(),tmparr[1].f());
  }

  cerr << "# R0: " << st.R0 << endl;

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

  cerr << "# total population: " << popsize << endl;

  cerr << "# shpsel: " << shpsel << endl;

  loadShape(st,fshape,shpsel,lonLimit);
//  loadPopDens("data/popdensmaps/gpw_v4_population_density_rev11_2020_2pt5_min.tif");

  loadPopDens("data/popdensmaps/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif"); // loads population counts raster data

  //TODO: cleanup the dependecies to make them explicit between the initializations and loading of data
  uint8_t *frameraw = new uint8_t[vwidth*vheight*4];
  initCairo(frameraw,vwidth,vheight); // initCairo after loadShape because it initializes a shapemask and after loading population raster data, to have a grid defined already

  adjustPopCounts(st,popsize);

  



/*
  edoublearray householdSizeDist(0.0,0.31,0.34,0.16,0.13,0.06);
  euintarray householdCounts;
  householdCounts.init(6,0); // how many households with single individuals, two individuals, ... max 6
  multinomial(popsize,householdSizeDist,householdCounts,rnd);
*/

  // Ferguson et al 2006 has household distribution sizes in supplementary graphs which the values below approximate
  // UK statistics have updated information on households
  st.householdSize=7;
  eintarray householdSizeDist;
  householdSizeDist.init(st.householdSize,0.0); // how many households with single individuals, two individuals, ... max 6
  householdSizeDist[0]=0; // household with 0 susceptible
  householdSizeDist[1]=(0.25+fsinglehh-flargehh*0.5)*popsize; 
  householdSizeDist[2]=(0.28-fsinglehh*0.2-flargehh*0.5)*popsize/2.0;
  householdSizeDist[3]=(0.20-fsinglehh*0.4+flargehh*0.2)*popsize/3.0; // 3 individual households, e.g.: two adults one child
  householdSizeDist[4]=(0.17-fsinglehh*0.3+flargehh*0.4)*popsize/4.0;
  householdSizeDist[5]=(0.08-fsinglehh*0.1+flargehh*0.2)*popsize/5.0;
  householdSizeDist[6]=(0.02+flargehh*0.2)*popsize/6.0;

  int total_households=0,total_hh_individuals=0;
  for (int i=0; i<householdSizeDist.size(); ++i){
    total_households+=householdSizeDist[i];
    total_hh_individuals+=householdSizeDist[i]*i;
  }
  householdSizeDist[1]+=popsize-total_hh_individuals; // Add missing individuals as single house holds
  total_households += popsize-total_hh_individuals;

  cerr << "# hh individuals: " << total_hh_individuals << "    " << popsize << endl;
  

  etable agegroup_infparams(etableLoad("data/agegroup.infparams",options));

  st.icu_symp=mularr(agegroup_infparams["Prop_symp_hospitalised"],agegroup_infparams["Prop_hospitalised_critical"]);
  st.nonicu_symp=mularr(agegroup_infparams["Prop_symp_hospitalised"],sumarr(mularr(agegroup_infparams["Prop_hospitalised_critical"],-1.0),1.0));
  st.death_icu=edoublearray(agegroup_infparams["Prop_critical_fatal"]);
  st.death_nonicu=edoublearray(agegroup_infparams["Prop_noncritical_fatal"]);
  // need to add a separate death from covid outside of hospital care to account for underreporting of covid deaths (around 20%)

  double totcritical=0.0,totdeaths=0.0;
  for (int i=0; i<agegroups.size(); ++i){
    totcritical+=agegroups[i]*st.rSA*st.icu_symp[i];
    totdeaths+=agegroups[i]*st.rSA*(st.icu_symp[i]*st.death_icu[i]+st.nonicu_symp[i]*st.death_nonicu[i])*df;
  }
  double tmpcritical=0.0,tmpdeaths=0.0;
  for (int i=0; i<agegroups.size(); ++i){
    tmpcritical+=agegroups[i]*st.rSA*st.icu_symp[i];
    tmpdeaths+=agegroups[i]*st.rSA*(st.icu_symp[i]*st.death_icu[i]+st.nonicu_symp[i]*st.death_nonicu[i])*df;
    cerr << "# age: " << 5*i << " critical: " << tmpcritical << " (" << tmpcritical/totcritical << ")" << " deaths: " << tmpdeaths << " (" << tmpdeaths/totdeaths << ")" <<endl;
  }


  
//  cout << agegroup_infparams << endl;
//  exit(0);


  // from Davies et al. 2020 covidm_params.R 
//  edoublearray dE(delay_gamma(4.0,4.0,60.0,0.25));  // Derived from Backer et al Eurosurveillance
//  edoublearray dIp(delay_gamma(2.4,4.0,60.0,0.25)); // Derived from Backer et al Eurosurveillance
//  edoublearray dIa(delay_gamma(7.0,4.0,60.0,0.25)); // Assumed 7 days subclinical shedding
//  edoublearray dIs(delay_gamma(3.2,3.7,60.0,0.25)); // Zhang et al 2020

  // from Davies et al. 2020 UK.R 
  st.dE=delay_gamma(4.0,4.0,60.0,0.25);  // 4 days delay
  st.dIp=delay_gamma(1.5,4.0,60.0,0.25); // 1.5 days presymptomatic
  st.dIs=delay_gamma(3.5,4.0,60.0,0.25); // 3.5 days symptomatic
  st.dIa=delay_gamma(5.0,4.0,60.0,0.25); // 5 days infectious as asymptomatic

  st.rIp=rgamma(1.5,4.0,60.0,0.25); // 1.5 days non-symptomatic
  st.rIs=rgamma(3.5,4.0,60.0,0.25); // 3.5 days as symptomatic
  st.rIa=rgamma(5.0,4.0,60.0,0.25); // 5 days infectious as asymptomatic

//  st.rTononicu=rgamma(4.0,4.0,60.0,0.25); // 4 days to go to NonICU, changed from Davis model which had 7 days
//  st.rToicu=rgamma(3.0,3.0,60.0,0.25); // 3 days to go to ICU after being in NonICU, this will now be after being in NonICU

  st.rToicu=rgamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.rTononicu=rgamma(7.0,7.0,60.0,0.25); // 7 days to go to NonICU

  st.rInicu=rgamma(10.0,10.0,60.0,0.25); // 10 days in ICU
  st.rInnonicu=rgamma(8.0,8.0,60.0,0.25); // 8 days in nonICU
  st.rIpDeath=rgamma(22.0,22.0,60.0,0.25); // deaths occuring 22 days after symptoms

  st.R0day=st.R0/5.0; // per day rate given that individuals are infectious for 5 days


  st.hhLevels2.init(total_households);
  ldieif(householdSizeDist.size()!=maxhousesize,"household size dist must be the same as maxhouse size");

  cerr << "# total households: " << total_households << endl;
  st.households.init(total_households);


  earray<eintarray> hhlist;
  hhlist.init(householdSizeDist.size());


  // initialize households
  int hi=0,ai=0;
  for (int i=1; i<householdSizeDist.size(); ++i){
//    uint8_t hl=i*(i+1)/2+i;
//    st.hhLevels[hl].reserve(householdSizeDist[i]);
    for (int j=0; j<householdSizeDist[i]; ++j,++hi){
      hhlist[i].add(hi);
//      st.hhLevels[hl].push_back(hi);
      st.households[hi].hage=0;
      st.households[hi].group=0;
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



  st.popindiv.init(popsize);

  
  cerr << householdSizeDist << endl;

  int assigned=0;
  const int agegroupband=5; // years
  int  tries=0;
  bool incomplete;

  do {
    ++tries;
    ldieif (tries>=4,"did not manage to seed households, try increasing number of larger households with -fflargehh 0.01");
    incomplete=false;
    eintarray tmpag(agegroups);
    cerr << tmpag << endl;
  
    // TODO: make gaussian distributed ages, this would take care of imperfect matches
    // TODO: make this step a setup step requiring a different option to run that generates the household structure file.
    for (int i=householdSizeDist.size()-1; i>=1; --i){
      cerr << "# seeding individuals in households size: " << i << "   " << assigned << "   " << householdSizeDist[i] << "    " << (i==1?0:i-2)*householdSizeDist[i] << "    "  <<  tmpag[0]+tmpag[1]+tmpag[2]+tmpag[3] << endl;
      if (i==1 && tmpag[0]+tmpag[1]+tmpag[2]+tmpag[3]>0){
        incomplete=true;
        break;
      }
      if (i==1)
        cerr << tmpag << endl;
      for (int j=0; j<hhlist[i].size(); ++j){
        shousehold &sh(st.households[hhlist[i][j]]);
        do {
          ldieif(sh.ageind<0 || sh.ageind>=st.popindiv.size(),"out of bounds: "+estr(sh.ageind)+" "+st.popindiv.size());
          st.popindiv[sh.ageind].age=int(5+(i-2)*0.5)+rnd.uniformint((tmpag.size()-int(5+(i-2)*0.5)-MAX(0,2*(i-2))));  // maxint is needed when the expression contains randomly generated numbers. Using the MAX macros causes the number to be generated twice!!
          ldieif(st.popindiv[sh.ageind].age<0 || st.popindiv[sh.ageind].age>=tmpag.size(),"out of bounds: "+estr().sprintf("%hhi",st.popindiv[sh.ageind].age)+" "+tmpag.size());
  //        st.pop_ages[sh.ageind]=5+(i-1)+rnd.uniform()*(tmpag.size()-5-(i-1)-MAX(0,2*(i-4)));
        } while (tmpag[st.popindiv[sh.ageind].age]==0 || i>=2 && tmpag[st.popindiv[sh.ageind].age]+tmpag[st.popindiv[sh.ageind].age-1]+tmpag[st.popindiv[sh.ageind].age-2]<=1); // if no individual of this age exists, or if the household is larger than 1 and only one individual exists that is younger or the same age (this will prevent the seeding from finishing in the next step)
        --tmpag[st.popindiv[sh.ageind].age];
        ++assigned;
  
  //      // use single agegroup for hhLevels
  //      sh.hage=0; //st.pop_ages[sh.ageind];
  
        // age group of "oldest" person in house
        sh.hage=st.popindiv[sh.ageind].age;
        sh.group=(sh.hage>=agethres/5?1:0); // set to 1 or 0 depending on age
  
        if (i>1){
          do {
            ldieif(sh.ageind+1<0 || sh.ageind+1>=st.popindiv.size(),"out of bounds: "+estr(sh.ageind+1)+" "+st.popindiv.size());
            st.popindiv[sh.ageind+1].age=minint(st.popindiv[sh.ageind].age+int(rnd.uniformint(3))-2,tmpag.size()-1);
            ldieif(st.popindiv[sh.ageind+1].age<0 || st.popindiv[sh.ageind+1].age>=tmpag.size(),"out of bounds: "+estr().sprintf("%hhi",st.popindiv[sh.ageind+1].age)+" "+tmpag.size());
          } while (tmpag[st.popindiv[sh.ageind+1].age]==0);
          --tmpag[st.popindiv[sh.ageind+1].age];
  //        if (avgage)
  //          sh.hage=(sh.hage+st.pop_ages[sh.ageind+1])/2;
          ++assigned;
        }   
        // Children
        for (int l=2; l<i; ++l){
          do {
            if (tmpag[MAX(0,st.popindiv[sh.ageind].age-9)]+tmpag[MAX(0,st.popindiv[sh.ageind].age-9)+1]+tmpag[MAX(0,st.popindiv[sh.ageind].age-9)+2]+tmpag[MAX(0,st.popindiv[sh.ageind].age-9)+3]==0)
              st.popindiv[sh.ageind+l].age=MAX(0,st.popindiv[sh.ageind].age-9) + 4 + int(rnd.uniformint(3));
            else if (tmpag[MAX(0,st.popindiv[sh.ageind].age-9)]+tmpag[MAX(0,st.popindiv[sh.ageind].age-9)+1]+tmpag[MAX(0,st.popindiv[sh.ageind].age-9)+2]==0)
              st.popindiv[sh.ageind+l].age=MAX(0,st.popindiv[sh.ageind].age-9) + 3;
            else
              st.popindiv[sh.ageind+l].age=MAX(0,st.popindiv[sh.ageind].age-9) + int(rnd.uniformint(3));
  /*
            if (tmpag[0]+tmpag[1]+tmpag[2]+tmpag[3]+tmpag[4]==0)
              st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)+1);
            else
              st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)-2);
  */
          } while (tmpag[st.popindiv[sh.ageind+l].age]==0);
          --tmpag[st.popindiv[sh.ageind+l].age];
          ++assigned;
  //        if (avgage)
  //          sh.hage=(sh.hage*l+st.pop_ages[sh.ageind+l])/(l+1);
        }
      }
    }
  }while(incomplete);
  cerr << "# done populating households: " << assigned << endl;


  st.spGridW=rwidth;
  st.spGridH=rheight;
  st.spGridSize=st.spGridW*st.spGridH;
  cerr << "# initializing spatial grid: " << st.spGridW << "x" << st.spGridH << endl;
  st.spGrid.init(st.spGridSize); // init 100 x 100 grid
  for (int i=0; i<st.spGrid.size(); ++i){
    for (int pg=0; pg<ngroups; ++pg){
      st.spGrid[i].E=0;
      st.spGrid[i].N[pg]=0;
      st.spGrid[i].Ip[pg]=0;
      st.spGrid[i].Is[pg]=0;
      st.spGrid[i].Ia[pg]=0;
      for (int l=0; l<nlevels; ++l){
        st.spGrid[i].hhIp[pg][l]=0;
        st.spGrid[i].hhIs[pg][l]=0;
        st.spGrid[i].hhIa[pg][l]=0;
      }
    }
  }

  // initialize number of people in group 0 (younger population) and 1 (older maybe isolated from 0)
  for (int i=0; i<st.households.size(); ++i){
    shousehold &hh(st.households[i]);
    st.N[hh.group]+=hh.size;
  }

  eintarray hhids;
  hhids.init(st.households.size());
  for (int i=0; i<hhids.size(); ++i)
    hhids[i]=i;
//  permute(hhids,0,hhids.size(),rnd);


  int gassign=0;
  int gpopAll=0;
  for (int i=0; i<st.spGrid.size(); ++i){
    sgrid &g(st.spGrid[i]);
    int gpop=popCounts[i];
    if (gpop==0) continue;

    if (popsize-gpopAll<10000) { cerr << " pop left to assign: " << popsize-gpopAll << " households left: " << hhids.size() << " assigned households: " << gassign << " gridposition: " << i << "  gpop: " << gpop << endl; }
    if (gpop > popsize-gpopAll) { cerr << "more individuals needed than available: " << gpop << " " << popsize-gpopAll << " gpos: " << i << " gsize: " << st.spGrid.size() << endl; /*exit(-1);*/ break; }
    int hhcount=samplehh(hhids,st.households,gpop,rnd); // find a random set of hh that sum to popsize in the grid and put it at the end of hhids
    for (int l=0; l<hhcount; ++l){
      shousehold &hh(st.households[hhids[hhids.size()-1]]);
      hh.gridpos=i;
      g.N[hh.group]+=hh.size;
      gpopAll+=hh.size;
      gpop-=hh.size;
      ++gassign;
      hhids.erase(hhids.size()-1);
    }
    ldieif(popsize-gpopAll>2000 && gpop!=0,"gpop not zero? gpop: "+estr(gpop)+" i: "+estr(i));
  }
  ldieif(popsize-gpopAll!=0 || hhids.size()>0,"seeding households on grid failed, population size: " + estr(popsize) + " individuals in households on grid: " +estr(gpopAll)+" households left: "+hhids.size());

/*
//  gpopAllCounts
  eintarray gposleft;

  int gassign=0;
  int gpopAll=0;
  hi=0;
  for (int i=0; i<rwidth*rheight; ++i){
    int gpop=raster[i];
    for (;hi<hhids.size() && st.households[hhids[hi]].size <= gpop; ++hi){
      shousehold &hh(st.households[hhids[hi]]);
      gpop-=hh.size;
      hh.gridpos=i;
      gpopAll+=hh.size;
      ++gassign;
      st.spGrid[hh.gridpos].N[hh.group]+=hh.size;
    }
    if (gpop>0)
      gposleft.add(i);
  }
  for ( ;hi<hhids.size(); ++hi){
    int ri;
    shousehold &hh(st.households[hhids[hi]]);
    cout << total_households-hi << " " << popsize-gpopAll << " " << hh.size << endl;
    do{
      ri=rnd.uniformint(gposleft.size());
    }while (hh.size>raster[gposleft[ri]]-st.spGrid[gposleft[ri]].N[0]-st.spGrid[gposleft[ri]].N[1]);
    sgrid &g(st.spGrid[gposleft[ri]]);
    hh.gridpos=gposleft[ri];
    gpopAll+=hh.size;
    ++gassign;
    g.N[hh.group]+=hh.size;
    if (raster[gposleft[ri]]-g.N[0]-g.N[1]==0){
      if (ri!=gposleft.size()-1)
        swap(gposleft[ri],gposleft[gposleft.size()-1]);
      gposleft.erase(gposleft.size()-1);
    }
  }
  cout << "gassign: " << gassign << endl;
  cout << "gridposleft: " << gposleft.size() << endl;
  cout << "gpopleft: " << popsize-gpopAll << endl;
*/

/*
  int gpopAllCounts=0;
  eintarray gpopCounts;
  gpopCounts.init(rwidth*rheight,0);
  for (int i=0; i<gpopCounts.size(); ++i){
    gpopCounts[i]=raster[i];
    gpopAllCounts+=raster[i];
  }

  // randomly place households across grid
  for (int i=0; i<st.households.size(); ++i){
    shousehold &hh(st.households[hhids[i]]);

    do{
      hh.gridpos=rnd.uniformint(st.spGrid.size());
    }while (gpopCounts[hh.gridpos]<hh.size);
    gpopCounts[hh.gridpops]-=hh.size;
    st.spGrid[hh.gridpos].N[hh.group]+=hh.size;
    gpopAllCounts-=hh.size;
    if (gpopAllCounts<100)
      cout << gpopAllCounts << endl;
  }
*/

  cerr << "# initializing levels" << endl;


  int hhOlder=0;
  int hhOlderN=0;
  int hhOlderNO=0;
  int hhOlderNY=0;
  int hhYoungNO=0;
  int hhYoungNY=0;

  eintarray hhcounts;

  // household ages and gridposition need to be initialized already
  hhcounts.init(st.spGrid.size()*ngroups*householdSizeDist.size(),0);
  for (int i=0; i<st.households.size(); ++i){
    shousehold& hh(st.households[i]);
    ++hhcounts[hh.gridpos*ngroups*householdSizeDist.size() + hh.group*householdSizeDist.size() + hh.size];
  }


//  cout << "# households by age: " << endl;
//  cout << hhagecount << endl;


  st.hhLevelsBegin.init(st.spGrid.size()*ngroups*householdSizeDist.size()*nlevels+1);

  st.grow=ngroups*householdSizeDist.size()*nlevels;
  st.arow=householdSizeDist.size()*nlevels;

  cerr << "# initializing levels indices" << endl;
  // initializes the start of each hhLevel subarray inside the large single array using the number of households per segment
  for (int gi=0,ipos=0; gi<st.spGrid.size(); ++gi){
    for (int hg=0; hg<ngroups; ++hg){
      for (int hsize=0; hsize<householdSizeDist.size(); ++hsize){
        for (int l=0; l<hsize; ++l) // set all levels to beginning of hh age subarray
          st.hhLevelsBegin[gi*st.grow + hg*st.arow + hsize*(hsize+1)/2+hsize-l]=ipos;
        ipos+=hhcounts[gi*ngroups*householdSizeDist.size() + hg*householdSizeDist.size() + hsize];
      }
    }
  }

  // end of level overlaps with beginning of next, therefore we have to reset last elements
  // once the levels are populated below, the elements will have the correct value
  for (int gi=0,ipos=0; gi<st.spGrid.size(); ++gi){
    for (int hg=0; hg<ngroups; ++hg){
      for (int hsize=0; hsize<householdSizeDist.size(); ++hsize)
        st.hhLevelsBegin[gi*st.grow + hg*st.arow + hsize*(hsize+1)/2+hsize+1]=st.hhLevelsBegin[gi*st.grow + hg*st.arow + hsize*(hsize+1)/2+hsize]; // this value will be incremented when populated below
    }
  }

  cerr << "# populating levels indices" << endl;
  // populate hhLevels
  for (int i=0; i<st.households.size(); ++i){
    shousehold &sh(st.households[i]);
    int hg=sh.group;
    uint8_t hl=sh.size*(sh.size+1)/2+sh.size;

    // store household indice in levels array
    st.hhLevels2[st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl+1]]=i;
    sh.hstind=st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl+1]-st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl];
    ++st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl+1]; // end of subarray is beginning of the next level

    sindiv *agearr=&st.popindiv[sh.ageind];
    if (hg>=1){
      ++hhOlder;
      hhOlderN+=sh.size;
      for (int l=0; l<sh.size; ++l){
        if (agearr[l].age>=agethres/5)
          ++hhOlderNO;
        else
          ++hhOlderNY;
      }
    }else{
      for (int l=0; l<sh.size; ++l){
        if (agearr[l].age>=agethres/5)
          ++hhYoungNO;
        else
          ++hhYoungNY;
      }
    }
  }


  cerr << "# hhOlder: " << hhOlder << " (" << double(hhOlder)/total_households << ")" << endl;
  cerr << "# hhOlderN: " << hhOlderN << " (" << double(hhOlderN)/popsize << ")" << endl;
  cerr << "# hhOlderNO: " << hhOlderNO << " (" << double(hhOlderNO)/popsize << ")" << endl;
  cerr << "# hhOlderNY: " << hhOlderNY << " (" << double(hhOlderNY)/popsize << ")" << endl;
  cerr << "# hhYoungNO: " << hhYoungNO << " (" << double(hhYoungNO)/popsize << ")" << endl;
  cerr << "# hhYoungNY: " << hhYoungNY << " (" << double(hhYoungNY)/popsize << ")" << endl;



  // randomize ages in households since infection always progresses sequentially in the age array
  for (int i=0; i<st.households.size(); ++i){
    shousehold &hh(st.households[i]);
    if (hh.size==1) continue; // cannot randomize single households
    uint8_t hl=hh.size*(hh.size+1)/2+hh.size;
    permute(st.popindiv,hh.ageind,hh.size,rnd);
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

//  st.evqueue.resize(st.dE.size()+1); // keep 1 extra for using as buffer for next step when adding new events

  double fE=0.0; // fraction of exposed (includes deaths, active infections and healed)


//  st.allIp=fseed; // directly as presymptomatic
//  st.evqueue.add(st,6,2,0,counts,st.dIp,rnd,2); 
//  st.allE=fseed;
//  cout << "hhLevelsBegin: " << st.hhLevelsBegin[1010*st.grow + 0*st.arow + 2*(2+1)/2+2+1] << " " << st.hhLevelsBegin[1010*st.grow + 0*st.arow + 2*(2+1)/2+2] << endl;

/*
  st.seedGrid=getGrid(46.0,8.95);
  if (st.seedGrid==-1) { cerr << "# coordinate not found chosing most populated area" << endl; st.seedGrid=imaxpop; }
  else if (st.hhLevelsBegin[st.seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[st.seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2]<st.fseed){
    cerr << "igrid: " << st.seedGrid << " levels empty, looking for another grid pos" << endl;
    int ix=st.seedGrid%st.spGridW;
    int ih=st.seedGrid/st.spGridW;
    for ( ;ix<st.spGridW && st.hhLevelsBegin[st.seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[st.seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2]<st.fseed; ++ix,++st.seedGrid);
    cerr << "seedGrid: " << st.seedGrid << " " << st.hhLevelsBegin[st.seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+3]-st.hhLevelsBegin[st.seedGrid*st.grow + 0*st.arow + 2*(2+1)/2+2] << endl;
    if (ix==st.spGridW) { cerr << "did not find any populated areas close to choice" << endl; st.seedGrid=imaxpop; }
  }
*/
//  st.evqueue.add(st,igrid,0,2,0,counts,st.dE,rnd);  // add Exposed to grid position with highest population count

  st.gridShuffle.init(st.spGrid.size());
  for (int i=0; i<st.gridShuffle.size(); ++i)
    st.gridShuffle[i]=i;
  permute(st.gridShuffle,0,st.gridShuffle.size(),rnd);


  cerr << "# starting simulation" << endl;

  double tpeak=0.0;
  double peakIs=0.0;
  double peakICU=0.0;
  double peakNonICU=0.0;

//  ldieif(videoOpen(1024,768,5,1000)!=0,"error creating video file");
  if (ovideo.len())
    ldieif(videoOpen(ovideo,1024,768,3,2000)!=0,"error creating video file");

//  edoublearray localIprob;
//  edoublearray localIprobGaussian;

/*
  edoublearray gaussianKernel;
  gaussianKernel.init(11*11,0.0);
  gaussianKernel[11*5+5]=0.6;
  gaussianKernel[11*6+5]=0.1;
  gaussianKernel[11*4+5]=0.1;
  gaussianKernel[11*5+4]=0.1;
  gaussianKernel[11*5+6]=0.1;
*/

  st.travelKernel.init(11*11,0.0);
  st.travelKernelSize=11;

/*
  st.travelKernel[11*5+5]=0.6;
  st.travelKernel[11*6+5]=0.1;
  st.travelKernel[11*4+5]=0.1;
  st.travelKernel[11*5+4]=0.1;
  st.travelKernel[11*5+6]=0.1;
*/

  double sump=0.0;
  for (int i=0; i<11; ++i){
    for (int j=0; j<11; ++j){
      int xf=double(i-5)*cos(M_PI*(ylatMin*0.5+ylatMax*0.5)/180.0),yf=double(j-5); // normalize x-axis distance by average latitude correction due to angle to distance correction
      st.travelKernel[j*11+i]=gsl_ran_gaussian_pdf(sqrt(xf*xf+yf*yf),ftravel);
      sump+=st.travelKernel[j*11+i];
    }
  }

//  cout << gaussianKernel << endl;
  for (int i=0; i<st.travelKernel.size(); ++i)
    st.travelKernel[i]/=sump;


  
  st.localIprob.init(st.spGrid.size());


  cout << "Time" << "\t" << "Isolation" << "\t" << "OlderIsolation" << "\t" << "Reff" << "\t" << "fE" << "\t" << "allE" << "\t" << "allCases" << "\t" << "newCases" << "\t" << "allIa" << "\t" << "allIp" << "\t" << "allIs" << "\t" << "allICU" << "\t" << "allNonICU" << "\t" << "Deaths";
  for (int i=1; i<printregions && st.shapes.size(); ++i)
    cout << "\tallCases_" << i << "\tallICU_" << i << "\tallNonICU_" << i << "\tallD_" << i;
  cout << endl;



  // setup seed array
  for (int i=0; i<iseed.size(); ++i){
    tmparr=iseed[i].explode(":");
    ldieif(tmparr.size()!=1 && tmparr.size()!=3,"wrong format for iseed: <lat>:<lon>:<infected>,<lat2>:<lon2>:<infected2>,...  or <infected>");
    if (tmparr.size()==3){
      int count=tmparr[2].i();
      int igrid=findPopGrid(st,tmparr[0].f(),tmparr[1].f(),count);
      if (igrid>=0){
        cerr << "# seeding position: " << iseed[i] << endl;
        st.iseedArr.add(igrid,count);
      }
    }else{
      cerr << "# seeding most populated area with: " << iseed[i] << endl;
      st.iseedArr.add(imaxpop,tmparr[0].i());
    }
  }


  st.ageCounts.init(agegroups.size());

  ethreads t;

  t.setThreads(st.nthreads);
  st.threadDone=st.nthreads;
  for (int i=0; i<st.nthreads; ++i){
    st.threadStates.addref(new sthreadState);
    st.threadStates[i].r.setSeed(rnd.uniformint(65565)); // initialize thread random seeds depending on main rng
    st.threadStates[i].evqueue.resize(st.dE.size()+1); // keep 1 extra for using as buffer for next step when adding new events
    st.threadStates[i].ageCounts.init(agegroups.size());
    st.threadStates[i].shapeCounts.init(st.shapeCounts.size());
  }

  t.run(threadRun,evararray(evarRef(st)));

  st.mutex.lock();
  while (st.threadDone>0) st.doneSignal.wait(st.mutex);
  st.mutex.unlock();
  cerr << "#starting simulation" << endl;





  st.mutex.lock();
  st.threadFunc=threadSeedInfections; // update events 
  st.threadI=st.nthreads;
  st.threadDone=st.nthreads;
  st.stateSignal.broadcast();
  while (st.threadDone>0) st.doneSignal.wait(st.mutex);
  st.mutex.unlock();



  time_t rawtime;
  time(&rawtime);
  tm *timeinfo = localtime(&rawtime);
  timeinfo->tm_year=2020;
  timeinfo->tm_mon=0;
  timeinfo->tm_mday=20;
  mktime(timeinfo);
  

  int lastCases=0,lastCases2=0;
  int newCases=0;
  int newICU=0,lastICU=0,lastICU2=0;
  char tmpsz[255];
  double Reff=1.0;
  double avgDeltaICU=0.0;
  int triggerCount=0;

  double tmpicusmr=-1.0;

  double lastTrigger=-10.0;

  double isoStart=-1.0,isoEnd=-1.0;

  for (it=0; it*st.tstep<tmax; ++it){
    fE=double(st.allE)/popsize;

    if (peakIs < st.Is[0]+st.Is[1]){
      peakIs=st.Is[0]+st.Is[1];
      tpeak=it*st.tstep;
    }
    if (peakICU < st.allICU)
      peakICU=st.allICU;
    if (peakNonICU < st.allNonICU)
      peakNonICU=st.allNonICU;

    while (eventsArr.size()>0 && eventsArr.keys(0) <= it*st.tstep){
      st.smr=eventsArr.values(0);
      eventsArr.erase(0);
    }
    while (oeventsArr.size()>0 && oeventsArr.keys(0) <= it*st.tstep){
      st.soldr=oeventsArr.values(0);
      oeventsArr.erase(0);
    }

    if (it%4==0){
      strftime(tmpsz,255,"%a %d %b",timeinfo);
      newCases=st.allCases-lastCases;
      newICU=st.allICU-lastICU;
//      newDeltaICU=newICU/lastNewICU;
      
//      Reff=(lastCases==lastCases2?1.0:5.0*double((st.allCases-lastCases))/(lastCases-lastCases2)); 
      Reff=Reff*2.0/3.0+((lastCases==lastCases2?1.0:pow(double(st.allCases-lastCases)/(lastCases-lastCases2),5.0))/3.0); 
//      ReffICU=ReffICU*2.0/3.0+((lastICU==lastICU2?1.0:pow(double(st.allICU-lastICU)/(lastICU-lastICU2),5.0))/3.0); 
      avgDeltaICU=avgDeltaICU*2.0/3.0+(st.allICU-lastICU); 
//      cout << "st.allCases: " << st.allCases << " lastCases: " << lastCases << " lastCases2: " << lastCases2 << endl;
//      Reff=(st.allCases-lastCases)/(lastCases-lastCases2); 
      if (ovideo.len())
        renderFrame(tmpsz,st.smr,newCases,frameraw,st,vwidth,vheight);
      lastCases2=lastCases;
      lastCases=st.allCases;
      lastICU2=lastICU;
      lastICU=st.allICU;
//      rawtime+=60*60*24; // number of seconds per day
      ++timeinfo->tm_mday;
      mktime(timeinfo);
//      if (tmpicusmr<0.0 && st.allICU+18*newICU+18*(18+1)*newDeltaICU/2>iculimit){

      if (oisotrigger>=0 && newCases>oisotrigger && isoStart<0.0){
        st.soldr=foiso;
        isoStart=floor(it*st.tstep);
      }
//      if (oisorelease>=0.0 && fE>=oisorelease && newCases<oisotrigger){
      if (oisorelease>=0.0 && fE>=0.5 && newCases<oisorelease && isoEnd<0.0){
        st.soldr=0.6;
        isoEnd=floor(it*st.tstep);
      }

//      if (iculowlimit>=0 && Reff<rthres && st.allICU<iculowlimit && tmpicusmr>=0.0 && it*st.tstep-lastTrigger>=triggerdays && ReffICU<=0.9){
      if (iculowlimit>=0 && Reff<rthres && st.allICU<iculowlimit && tmpicusmr>=0.0 && it*st.tstep-lastTrigger>=triggerdays && avgDeltaICU<icuthres){
        --triggercount;
        ficu=MIN(fstep+ficu,1.0);
//        if (triggercount==0 || st.allCases>0.08*popsize) ficu=1.0;
        if (triggercount==0) ficu=1.0;
        st.smr=ficu; //0.40*ftrigger;
        lastTrigger=it*st.tstep;
        cerr << "# Reducing isolation: " << Reff << " " << st.smr << endl;
        if (st.smr==1.0) tmpicusmr=-1.0;
      } else if (icutrigger>=0.0 && newICU>icutrigger && triggerCount<1){
        if (tmpicusmr<0.0){
          lastTrigger=it*st.tstep;
          tmpicusmr=st.smr;
//          st.smr=clamp(0.0,1.0,0.8*1.0/st.R0);
          st.smr=ficu; //0.3*ftrigger;
          cerr << "# ICU capacity triggered icu limit: " << newICU << " " << st.smr << endl;
          ++triggerCount;
        }
      }
/*
      if (iculowlimit>=0 && Reff<rthres2 && st.allICU<iculowlimit2 && tmpicusmr>=0.0){
        st.smr=tmpicusmr;
        tmpicusmr=-1.0;
        cerr << "# low ICU remove icu limit: " << " " << st.smr << endl;
      } else if (iculowlimit>=0 && Reff<rthres && st.allICU<iculowlimit && tmpicusmr>=0.0){
        st.smr=ficu2; //0.40*ftrigger;
        cerr << "# Reducing isolation: " << " " << st.smr << endl;
      } else if (icutrigger>=0.0 && newICU>icutrigger && triggerCount<1){
        if (tmpicusmr<0.0){
          tmpicusmr=st.smr;
//          st.smr=clamp(0.0,1.0,0.8*1.0/st.R0);
          st.smr=ficu; //0.3*ftrigger;
          cerr << "# ICU capacity triggered icu limit: " << st.allICU*exp(5*log(Reff)) << " " << st.smr << endl;
          ++triggerCount;
        }
      }
*/
/*
      if (iculimit>=0 && Reff<rthres2 && st.allICU<0.8*iculimit && tmpicusmr>=0.0){
        st.smr=tmpicusmr;
        tmpicusmr=-1.0;
        cerr << "# low ICU remove icu limit: " << " " << st.smr << endl;
      } else if (iculimit>=0 && Reff<rthres && st.allICU<0.9*iculimit && tmpicusmr>=0.0){
        st.smr=0.5*ftrigger;
        cerr << "# Reducing isolation: " << " " << st.smr << endl;
      } else if (casestrigger>=0.0 && newCases>casestrigger && triggerCount<1){
        if (tmpicusmr<0.0){
          tmpicusmr=st.smr;
//          st.smr=clamp(0.0,1.0,0.8*1.0/st.R0);
          st.smr=0.3*ftrigger;
          cerr << "# ICU capacity triggered icu limit: " << st.allICU*exp(5*log(Reff)) << " " << st.smr << endl;
          ++triggerCount;
        }
      }
*/
    }

//    double interhhrate=(hage>=10?folder:1.0)*finter*st.R0*(1.0-fE)*(0.5*st.allIa+st.allIp+st.allIs)/(popsize-st.allE);

    // The hhOlderN term is needed to prevent R0 from decreasing when older population is isolated from younger population
    // Whithout including this term, the simulation is equivalent to an older population considered immune but still interacting, thus reducing incorrectly the R0

    // Global probability of infection from random person in country (x0.01)
//    double globalIprob=0.01 * finter*R0day*(0.5*st.allIa+st.allIp+st.allIs - (1.0-MIN(smr,soldr)/smr)*(0.5*st.oallIa+st.oallIp+st.oallIs) )/(popsize-(1.0-MIN(smr,soldr)/smr)*hhOlderN);
//    double globalIprob=st.finter*fglobal*R0day*(0.5*st.Ia[0]+st.Ip[0]+st.Is[0] + (MIN(smr,soldr)/smr)*(0.5*st.Ia[1]+st.Ip[1]+st.Is[1]) )/(st.N[0]+(MIN(smr,soldr)/smr)*st.N[1]);

    // removed the cap because its more conservative and because there are situations where there should be no cap. For a more accurate representation of reality need to model how transmission occurs
    // through contacts, considering how time spent and number of contacts influence the transmission probability, and what is the effect of the isolation measures on either of these two variables

//    if (globalIprob>1.0) globalIprob=1.0; // The reason to cap the probability is to prevent the unrealistic situation where "isolation" has no effect on the probability of infection. Need to consider more carefully  what capping this value implies in the model and what the wider implications are
    // The cap means that the maximum probability of getting infected per day cannot be larger than 1, and because of this a mitigation factor of 50% will
    // guarantee a reduction in 50% probability of getting infected. This however implies there is a saturation of the probability of infection, 
    // To put this more mechanistically: if we assume isolation reduces the time a person spends with other contacts or number of contacts (with constant time per contact), and we assume the probability of infection is proportional to the time exposed to an infected individual. Then isolation will always have an effect on probability of infection.
    // There are cases when this is not true, if the probability of infection is so high that it is enough to have a single contact even of very short duration to have a probability of 100% of getting infected and the proportion of infected is high enough to guarantee there is one infected person in almost every group of contacts, then reducing contact (either by decreasing the time or by decreasing the number of contacts) should correctly result in no reduced probability of infection, up to a certain point.

    st.mutex.lock();
    st.threadFunc=threadLocalInfection; // update events 
    st.threadI=st.nthreads;
    st.threadDone=st.nthreads;
    st.stateSignal.broadcast();
    while (st.threadDone>0) st.doneSignal.wait(st.mutex);
    st.mutex.unlock();

/*
    for (int gi=0; gi<st.spGrid.size(); ++gi){
      sgrid& g(st.spGrid[gi]);
      if (g.N[0]+g.N[1]==0) { localIprob[gi]=0.0; continue; }
      localIprob[gi]=finter*R0day*(0.5*g.Ia[0]+g.Ip[0]+g.Is[0] + (MIN(smr,soldr)/smr)*(0.5*g.Ia[1]+g.Ip[1]+g.Is[1]) )/(g.N[0]+(MIN(smr,soldr)/smr)*g.N[1]);
    }
*/

//    blur(localIprobGaussian,localIprob,st.spGridW,st.spGridH,gaussianKernel,11);


    st.mutex.lock();
    st.threadFunc=threadUpdateGrid; // update events 
    st.threadI=st.nthreads;
    st.threadDone=st.nthreads;
    st.stateSignal.broadcast();
    while (st.threadDone>0) st.doneSignal.wait(st.mutex);
    st.mutex.unlock();

/*
    for (int gi=0; gi<st.spGrid.size(); ++gi){
      sgrid& g(st.spGrid[gi]);
      if (g.N[0]+g.N[1]==0) { continue; } // skip empty grid

      int gx=gi%st.spGridW;
      int gy=gi/st.spGridW;

      int giup=((st.spGridH+gy-1)%st.spGridH)*st.spGridW+gx;
      int gidn=((gy+1)%st.spGridH)*st.spGridW+gx;
      int gilt=gy*st.spGridW+(st.spGridW+gx-1)%st.spGridW;
      int girt=gy*st.spGridW+(gx+1)%st.spGridW;

      // TODO: handle adjacent empty grids or big differences in population sizes. How should they affect the probability of infection of an adjacent grid?
//      double tmpLocalIprob=localIprobGaussian[gi]*flocal; // + localIprob[gilt]*(1.0-flocal);
//      double tmpLocalIprob=localIprob[gi]*flocal; // + localIprob[gilt]*(1.0-flocal);
      double tmpLocalIprob=localIprob[gi]*0.6 + localIprob[gilt]*0.1 + localIprob[girt]*0.1 + localIprob[giup]*0.1 + localIprob[gidn]*0.1;

      for (int hg=0; hg<ngroups; ++hg){
        for (int hhsize=1; hhsize<householdSizeDist.size(); ++hhsize){ // TODO: should move levels to outter loop for efficiency, inner loop should be longest sequential array of data
          for (int hhexp=hhsize-1; hhexp>=0; --hhexp){ // NOTE: have to do decreasing here to avoid infecting twice the same house in the same step
  
            uint8_t hl=hhsize*(hhsize+1)/2+(hhsize-hhexp);
            unsigned int hlcount=st.hhLevelsBegin[gi*st.grow + hg*st.arow + hl+1]-st.hhLevelsBegin[gi*st.grow + hg*st.arow +hl];
            if (hlcount==0) continue; // no households in level
  
            int hS=hhsize-hhexp; // number of susceptible individuals per household
            int nS=hlcount*hS; // total number of susceptible individuals in these households
  
            // The intra household rate below is an approximation needed to avoid computing the exact probability of new infection per household which depends on the specific combination of Ia,Is and Ip
            double intraIprob=fintra*R0day*(0.5*g.hhIa[hg][hl]+g.hhIp[hg][hl]+g.hhIs[hg][hl])/nS;
    
            double p=(MIN(smr,(hg>=1?soldr:1.0))*(globalIprob + tmpLocalIprob) + intraIprob)*st.tstep; // probability of infection per susceptible person in these households
    
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
            multinomial(hlcount,mp,counts,rnd);
           
    //        int newE=binomial(nS,,rnd)/(hhsize-hhexp);
    //        if (newE>st.hhLevels[hhstlevel].size()) newE=st.hhLevels[hhstlevel].size(); // cap maximum number of new household exposures
   
            st.evqueue.add(st,gi,hg,hhsize,hhexp,counts,st.dE,rnd);
          }
        }
      }
    }
*/

/*
    st.mutex.lock();
    st.threadFunc=threadUpdateEvents; // update events 
    st.threadI=nthreads;
    st.threadDone=nthreads;
    st.stateSignal.signal();
    while (st.threadDone>0) st.doneSignal.wait(st.mutex);
    st.mutex.unlock();
*/

/*
    st.nextEvents=&nextEvents(st.evqueue.step());
    processEvents(st,nextEvents,rnd);
*/

//    cout << it*st.tstep << "\t" << (1.0-st.smr) << "\t" << (1.0-st.soldr) << "\t" << Reff << "\t" << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.allCases << "\t" << newCases << "\t" << st.Ia[0]+st.Ia[1] << "\t" << st.Ip[0]+st.Ip[1] << "\t" << st.Is[0]+st.Is[1] << "\t" << st.allICU << "\t" << st.allNonICU << "\t" << st.allD << endl;
    cout << it*st.tstep << "\t" << (1.0-st.smr) << "\t" << (1.0-st.soldr) << "\t" << Reff << "\t" << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.allCases << "\t" << newCases << "\t" << st.Ia[0]+st.Ia[1] << "\t" << st.Ip[0]+st.Ip[1] << "\t" << st.Is[0]+st.Is[1] << "\t" << st.allICU << "\t" << st.allNonICU << "\t" << st.allD;
    for (int i=1; i<printregions && st.shapes.size(); ++i)
      cout << "\t" << st.shapeCounts[i].allCases << "\t" << st.shapeCounts[i].allICU << "\t" << st.shapeCounts[i].allNonICU << "\t" << st.shapeCounts[i].allD;
    cout << endl;
//    if (st.Ia[0]+st.Ia[1]+st.Is[0]+st.Is[1]+st.Ip[0]+st.Ip[1]==0) break;
  }

  if (ovideo.len())
    videoClose();
//  cout << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.allCases << "\t" << newCases << "\t" << st.Ia[0]+st.Ia[1] << "\t" << st.Ip[0]+st.Ip[1] << "\t" << st.Is[0]+st.Is[1] << "\t" << st.allICU << "\t" << st.allNonICU << "\t" << st.allD << endl;

  cerr << "# rseed: " << rnd.seed << " fglobal: " << st.fglobal << " ftravel: " << ftravel << " cf: " << cf << " R0: " << st.R0 << " fE: " << double(st.allE)/popsize << " tpeak: " << tpeak << " peakIs: " << peakIs << " peakICU: " << peakICU << " peakNonICU: " << peakNonICU << " deaths: " << st.allD << " isotime: " << (isoEnd-isoStart)/30.0 << endl;

  st.mutex.lock();
  st.threadI=st.nthreads;
  st.threadFunc=0x00;
  st.stateSignal.broadcast();
  st.mutex.unlock();

  t.wait();

  return(0);
}
