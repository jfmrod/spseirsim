#include <eutils/emain.h>
#include <eutils/ernd.h>

#include <gsl/gsl_cdf.h>

#include <eutils/eparser.h>
#include <eutils/eparserinterpreter.h>
#include <eutils/etable.h>
#include <eutils/vector3.h>

#include <deque>

#include "videoenc.h"

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
    int ir=r.uniformint(size);
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
  uint8_t group;
  int gridpos;
  int ageind;
};

struct sevent {
  int hhid;
  uint8_t transition;
  uint8_t age;
};

class evqueuelist;
struct ssimstate;

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
  void add(ssimstate& st,int hpos,uint8_t hgroup,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r,int ntransition=0);
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

//  const int ngroups=2; // defined as global now

  int spGridW=0;
  int spGridH=0;
  int spGridSize=0;
  ebasicarray<sgrid> spGrid;

  ebasicarray<uint8_t> pop_ages;
  ebasicarray<shousehold> households;

  eintarray hhLevels2;
  eintarray hhLevelsBegin;
  int grow=0;
  int arow=0;
//  earray<earray<deque<int> > > hhLevels;
//  earray<eintarray> hhIa,hhIp,hhIs;

  int allE=0;

  int N[ngroups]={0,0};
  int Ip[ngroups]={0,0};
  int Is[ngroups]={0,0};
  int Ia[ngroups]={0,0};

  double allICU,allNonICU,allD;

  evqueuelist evqueue;
};

void processEvents(ssimstate& st,deque<sevent>& evs,ernd& r)
{
  for (int i=0; i<evs.size(); ++i){
    sevent &e(evs[i]);
    shousehold &hs(st.households[e.hhid]);
    uint8_t hl=hs.size*int(hs.size+1)/2+(hs.size-hs.E);
    sgrid &g(st.spGrid[hs.gridpos]);

    switch(e.transition){
      case 0:{ // end of E state
        if (rnd.uniform()<st.rSA){
          // E -> Ip
          ++st.Ip[hs.group];
          ++g.Ip[hs.group];
          ++g.hhIp[hs.group][hl];
          ++hs.Ip;
          e.transition=2;
          st.evqueue.add(st,e,st.rIp(r));
        }else{
          // E -> Ia
          ++st.Ia[hs.group];
          ++g.Ia[hs.group];
          ++g.hhIa[hs.group][hl];
          ++hs.Ia;
          e.transition=1;
          st.evqueue.add(st,e,st.rIa(r));
        }
       break;
      }
      case 1:{
        // Ia ->
        --st.Ia[hs.group];
        --g.Ia[hs.group];
        --g.hhIa[hs.group][hl];
        --hs.Ia;
       break;
      }
      case 2:{
        // Ip -> Is , Ip -> toicu , Ip -> tononicu , Ip -> death
        --st.Ip[hs.group];
        --g.Ip[hs.group];
        --g.hhIp[hs.group][hl];
        --hs.Ip;
        ++st.Is[hs.group];
        ++g.Is[hs.group];
        ++g.hhIs[hs.group][hl];
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
        --st.Is[hs.group];
        --g.Is[hs.group];
        --g.hhIs[hs.group][hl];
        --hs.Is;
       break;
      }
      case 4:{
        // TODO: remove individual from population (stop being infective) once he joins ICU or nonICU, requires keeping track of individual to make sure he is not removed twice, once from Is -> and another time from the ICU, nonICU or death cases, maybe use age value as flag when -1 means he already was removed

        // toicu -> icu
        ++st.allICU;
        e.transition=7;
        st.evqueue.add(st,e,st.rInicu(r));
       break;
      }

      case 5:{
        // tononicu -> nonicu
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
        --st.allICU;
       break;
      }
      case 8:{
        // nonicu -> 
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

void evqueuelist::add(ssimstate &st,int hpos,uint8_t hgroup,uint8_t hhsize,uint8_t hhexp,euintarray& count,const edoublearray& evdist,ernd& r,int ntransition)
{
  sevent ev;
  uint8_t hl=hhsize*int(hhsize+1)/2+(hhsize-hhexp);
//  uint8_t hhstnewlevel=hhsize*(hhsize+1)/2+(hhsize-(hhexp+1)); // increasing exposed count in household

  for (int ti=1; ti<count.size(); ++ti)
    st.allE+=count[ti]*ti;

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
        // all households indices in the different levels are all stored in a single large array grouped by levels
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
          ev.age=st.pop_ages[hst.ageind+hst.E+l];
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


SHPObject *psShape=0x00;

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
float *raster=0x00;
int rwidth=0;
int rheight=0;
float rmaxpop=0.0;
int rpopsize=0;

/*
inline float readPopDens(double lon,double lat){
  GTIFPCSToImage( gtif, &lon,&lat );
  return(raster[uint32_t(lat)*rwidth + uint32_t(lon)]==nodataval?0.0:raster[uint32_t(lat)*rwidth + uint32_t(lon)]);
}
*/

inline float readPopDensR(uint32_t ix,uint32_t iy){
  return(raster[iy*rwidth + ix]==nodataval?0.0:raster[iy*rwidth + ix]);
}



uint32_t v3_col(const evector3& v){ return(0xFFu&uint32_t(v.x*0xFFu) | 0xFF00u&uint32_t(v.y*0xFF00u) | 0xFF0000u&uint32_t(v.z*0xFF0000)); }

uint32_t tricolor(const evector3& col1,const evector3& col2,const evector3& col3,double p1,double p2){
  evector3 col12=col1*(1.0-p1)+col2*p1;
  return(v3_col((col1*(1.0-p1)+col2*p1)*(1.0-p2)+col3*p2));
}


double lonMin=0.0;
double lonMax=0.0;
double lonRef=0.0;

double proj(double lon,double lat,double reflon){
  return(cos(M_PI*lat/180.0)*(lon-reflon));
}

double scale=1.0;
double xpos=0.0,ypos=0.0;

double revproj(double lon,double lat,double reflon){
  return(lon/cos(M_PI*lat/180.0)+reflon + lonMin);
}

void cairoDrawShape(cairo_t *cr,SHPObject *shp){
  for (int j=0, iPart=1; j<shp->nVertices; ++j){
    if (j == 0 && shp->nParts > 0 )
      cairo_move_to(cr,scale*(proj(shp->padfX[j],shp->padfY[j],lonRef)-lonMin)+xpos,1080.0-scale*(shp->padfY[j]-shp->dfYMin)-ypos);
//       pszPartType = SHPPartTypeName( psShape->panPartType[0] );
            
    if (iPart < shp->nParts && shp->panPartStart[iPart]==j){
      cairo_close_path(cr);
      cairo_move_to(cr,scale*(proj(shp->padfX[j],shp->padfY[j],lonRef)-lonMin)+xpos,1080.0-scale*(shp->padfY[j]-shp->dfYMin)-ypos);
//       pszPartType = SHPPartTypeName( psShape->panPartType[0] );
//      pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
      iPart++;
    }else
      cairo_line_to(cr,scale*(proj(shp->padfX[j],shp->padfY[j],lonRef)-lonMin)+xpos,1080.0-scale*(shp->padfY[j]-shp->dfYMin)-ypos);
  }
}


void renderFrame(int day,float mitigation,uint8_t *frameRaw,ssimstate& st,int width,int height){
//  memset(frameRaw, 0, 1920*1080*4);
  for (int i=0; i<1920*1080*4; i+=4){
    frameRaw[i]=0xff;
    frameRaw[i+1]=0xff;
    frameRaw[i+2]=0xff;
    frameRaw[i+3]=0x00;
  }

  if ((lonMax-lonMin)/(psShape->dfYMax-psShape->dfYMin) > double(width)/height){
    scale=width/(lonMax-lonMin);
    ypos=(height-scale*(psShape->dfYMax-psShape->dfYMin))/2.0;
  }else{
    scale=height/(psShape->dfYMax-psShape->dfYMin);
    xpos=(width-scale*(lonMax-lonMin))/2.0;
  }
  scale=0.9*scale;
  ypos=(height-scale*(psShape->dfYMax-psShape->dfYMin))/2.0;
  xpos=(width-scale*(lonMax-lonMin))/2.0;

  double maxE=0.0;
  double maxI=0.0;
  for (int i=0; i<st.spGridSize; ++i){
    sgrid &g(st.spGrid[i]);
    double tmpI=double(g.Ia[0]+g.Ia[1]+g.Ip[0]+g.Ip[1]+g.Is[0]+g.Is[1])/(g.N[0]+g.N[1]);
    double tmpE=double(g.E)/(g.N[0]+g.N[1]);
    if (maxI<tmpI) maxI=tmpI;
    if (maxE<tmpE) maxE=tmpE;
  }

  for (int i=0; i<st.spGridW; ++i){
    for (int j=0; j<st.spGridH; ++j){
      double sxlon=i*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
      double sylat=j*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

      double xlon=(i-0.5)*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
      double ylat=(j-0.5)*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

      double xlonn=(i+1-0.5)*(psShape->dfXMax-psShape->dfXMin)/st.spGridW+psShape->dfXMin;
      double ylatn=(j+1-0.5)*(psShape->dfYMax-psShape->dfYMin)/st.spGridH+psShape->dfYMin;

      double sx=scale*(proj(xlon,ylat,lonRef)-lonMin)+xpos;
      double sy=1080.0-scale*(ylat-psShape->dfYMin)-ypos;
//      cout << "i: " << i << " j: " << j << " xlonylat: " << xlon << " " << ylat << " sxsy: " << sx << " " << sy <<endl;

      double sxn=scale*(proj(xlonn,ylatn,lonRef)-lonMin)+xpos;
      double syn=1080.0-scale*(ylatn-psShape->dfYMin)-ypos;

      sgrid &g(st.spGrid[j*st.spGridW+i]);
      double tmpI=double(g.Ia[0]+g.Ia[1]+g.Ip[0]+g.Ip[1]+g.Is[0]+g.Is[1])/(g.N[0]+g.N[1]);
      double tmpE=double(g.E)/(g.N[0]+g.N[1]);
      float pdens=clamp(0.0,1.0,(log(readPopDensR(i,rheight-j)+0.001)-log(0.001))/(log(rmaxpop+0.001)-log(0.001))*0.6+0.4);
//      float pdens=clamp(0.0,1.0,readPopDensR(i+rxmin,rymax-j)/maxpdens*0.6+0.4);
//      float pdens=clamp(0.0,1.0,readPopDens(i*(psShape->dfXMax-psShape->dfXMin)/100.0+psShape->dfXMin,j*(psShape->dfYMax-psShape->dfYMin)/100.0+psShape->dfYMin)/maxpdens);
      uint32_t color=tricolor(evector3(1.0,1.0,1.0)*pdens,evector3(0.3,0.8,0.3),evector3(1.0,0.0,0.0),tmpE,clamp(0.0,1.0,(log(tmpI+0.001)-log(0.001))/(-log(0.001))));
//      uint32_t color=0xFFFFFFu;
      colorRect(frameRaw,sx,sy,sxn,syn,color,width);
    }
  }


/*
  for (int i=(1920-1080)*4/2; i<width-(1920-1080)*4/2; ++i){
    frameRaw[0*4*width + (i*4)]=0x0u;
    frameRaw[0*4*width + (i*4)+1]=0x0u;
    frameRaw[0*4*width + (i*4)+2]=0x0u;
    frameRaw[1079*4*width + (i*4)]=0x0u;
    frameRaw[1079*4*width + (i*4)+1]=0x0u;
    frameRaw[1079*4*width + (i*4)+2]=0x0u;
  }
  for (int j=0; j<height; ++j){
    frameRaw[j*4*width + ((1920-1080)*4/2 + 0)]=0x0u;
    frameRaw[j*4*width + ((1920-1080)*4/2 + 0)+1]=0x0u;
    frameRaw[j*4*width + ((1920-1080)*4/2 + 0)+2]=0x0u;
    frameRaw[j*4*width + ((1920-1080)*4/2 + 1079*4)]=0x0u;
    frameRaw[j*4*width + ((1920-1080)*4/2 + 1079*4)+1]=0x0u;
    frameRaw[j*4*width + ((1920-1080)*4/2 + 1079*4)+2]=0x0u;
  }
*/

  // tell cairo we have changed the image contents
  cairo_surface_mark_dirty(surface);

  cairoDrawShape(cr,psShape);
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1);
  cairo_stroke (cr);

  cairo_move_to (cr, 10.0, 50.0);
  cairo_show_text (cr, (estr("Day: ")+day)._str);
  cairo_move_to (cr, 10.0, 1080.0-50.0);
  cairo_show_text (cr, estr().sprintf("Fraction exposed: %.3lf",double(st.allE)/(st.N[0]+st.N[1]))._str);
  if (mitigation<1.0){
    cairo_move_to (cr, 10.0, 190.0);
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

uint8_t *shapeMask=0x00;

void initCairo(uint8_t *frameRaw,int width,int height)
{
//  cairo_surface_t *surface=cairo_image_surface_create(CAIRO_FORMAT_RGB32, width, height);
  surface=cairo_image_surface_create_for_data(frameRaw,CAIRO_FORMAT_ARGB32,width,height,width*4);
  int cstatus=cairo_surface_status(surface);
  cout << "# stride RGB: " << cairo_format_stride_for_width (CAIRO_FORMAT_RGB24,width) << endl;
  cout << "# stride ARGB: " << cairo_format_stride_for_width (CAIRO_FORMAT_ARGB32,width) << endl;
  ldieif(cstatus!=CAIRO_STATUS_SUCCESS,"error creating surface: "+estr(cstatus));
  cr=cairo_create(surface); 

  cairo_select_font_face (cr, "serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size (cr, 32.0);
  cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);

  cairo_surface_t *surfacemask=0x00;
  cairo_t *crmask=0x00;


  shapeMask=new uint8_t[rwidth*rheight];
  memset(shapeMask,0,rwidth*rheight);
  surfacemask=cairo_image_surface_create_for_data(shapeMask,CAIRO_FORMAT_A8,rwidth,rheight,rwidth);
  cstatus=cairo_surface_status(surfacemask);
  ldieif(cstatus!=CAIRO_STATUS_SUCCESS,"error creating surface: "+estr(cstatus));
  crmask=cairo_create(surfacemask); 
  cairoDrawShape(crmask,psShape);
  cairo_set_source_rgb(crmask,1.0,1.0,1.0);
  cairo_fill(crmask);
  cairo_surface_flush(surfacemask);
  cairo_destroy(crmask);
  cairo_surface_destroy(surfacemask);
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



int loadPopDens(const estr& fname,int popsize){
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

  GTIFPrint(gtif,0,0);

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
  

  cout << "xsize: " << xsize << " ysize: " << ysize << " spp: " << spp << " sf: " << sf << " depth: " << depth << " bps: " << bps << " rowsperstrip: " << rowsperstrip << endl;

  GTIFDefn defn;
  if (GTIFGetDefn(gtif, &defn)){
    printf( "\n" );
    GTIFPrintDefnEx( gtif, &defn, stdout );

    printf( "\n" );
    printf( "PROJ.4 Definition: %s\n", GTIFGetProj4Defn(&defn));
  }

//  TIFFSetErrorHandler(TIFFError);

//	if (!TIFFReadRGBAImage(tif, xsize, ysize, raster, 0)) { lerror("reading TIFF data"); exit(-1); return(-1); }


  printf( "\nCorner Coordinates:\n" );

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

  cout << "upper left: " << c_ul_x << " " << c_ul_y << endl;
  cout << "upper right: " << c_ur_x << " " << c_ur_y << endl;
  cout << "bottom left: " << c_bl_x << " " << c_bl_y << endl;
  cout << "bottom right: " << c_br_x << " " << c_br_y << endl;

  uint32_t rxmin=0;
  uint32_t rxmax=0;
  uint32_t rymin=0;
  uint32_t rymax=0;

  double sX=psShape->dfXMin;
  double sY=psShape->dfYMin;
  cout << "top left of shape: " << sX << " " << sY << endl;
  GTIFPCSToImage( gtif, &sX,&sY );
  cout << "top left of shape: " << sX << " " << sY << endl;
//  cout << "value: " << raster[uint32_t(sY)*xsize + uint32_t(sX)] << " nodataval: " << nodataval << endl;
  rxmin=sX;
  rymin=sY;

  sX=psShape->dfXMax;
  sY=psShape->dfYMax;
  cout << "bottom right of shape: " << sX << " " << sY << endl;
  GTIFPCSToImage( gtif, &sX,&sY );
  cout << "bottom right of shape: " << sX << " " << sY << endl;
//  cout << "value: " << raster[uint32_t(sY)*xsize + uint32_t(sX)] << " nodataval: " << nodataval << endl;
  rxmax=sX;
  rymax=sY;

  if (rxmin>rxmax) swap(rxmin,rxmax);
  if (rymin>rymax) swap(rymin,rymax);
  cout << "sxmin: " << rxmin << " sxmax: " << rxmax <<endl;
  cout << "symin: " << rymin << " symax: " << rymax <<endl;
  rwidth=rxmax-rxmin;
  rheight=rymax-rymin;
  cout << "rwidth: " << rwidth << " rheight: " << rheight <<endl;
 
	int npixels = rwidth * rheight;
	raster = (float*) _TIFFmalloc(npixels * sizeof (float));
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
  for (uint32_t strip=0, row=0; row<ysize && row<rymax ; ++strip,row+=rowsperstrip) {
//    if (rp+stripsize/sizeof(uint32_t) > npixels) ldie("TIFF data larger than expected: "+estr(row)+" "+estr(stripsize)+" "+estr(rp)+" "+estr(strip)+" "+estr(npixels)); 
    bytesread=TIFFReadEncodedStrip(tif, strip, sbuf, stripsize);
    ldieif(bytesread==-1,"Error reading GeoTIFF");
    if (row>=rymin)
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
 
  rpopsize=0;
  int totalpopsize=0;
  rmaxpop=(raster[0]==nodataval?0.0:raster[0]);
  for (int i=0; i<rheight*rwidth; ++i){
    float tmppden=(raster[i]==nodataval?0.0:raster[i]);
    if (tmppden>rmaxpop) rmaxpop=tmppden;
    if (shapeMask[i]!=0)
      rpopsize+=tmppden;
    totalpopsize+=tmppden;
  }
  cout << "rpopsize: " << rpopsize << " totalpopsize: " << totalpopsize << endl;

  cout << "Adjusting raster population counts to: " << popsize << endl;
  // adjust population counts to match current country population census
  for (int i=0; i<rheight*rwidth; ++i){
    float tmpf=(raster[i]==nodataval?0.0:raster[i]);
    raster[i]=tmpf*double(rpopsize)/popsize;
  }
  rmaxpop=rmaxpop*double(rpopsize)/popsize;
  cout << "rmaxpop: " << rmaxpop << endl;

//  exit(0);

  return(0);
}


void loadShape(const estr& fname){
  SHPHandle hSHP;
  int nShapeType, nEntities;
  double adfMinBound[4],adfMaxBound[4];

  hSHP = SHPOpen( fname._str, "rb" );
  SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

  int nPrecision = 15;
  printf( "Shapefile Type: %s   # of Shapes: %d\n\n",
            SHPTypeName( nShapeType ), nEntities );
    
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

  ldieif(nEntities>1,"more than one shape entity found");

//  for (int i=0; i<nEntities; ++i){
    int j;
//    SHPObject *psShape;

    psShape = SHPReadObject(hSHP, 0);

    if (psShape == NULL) {
      fprintf( stderr,
                   "Unable to read shape %d, terminating object reading.\n",
                    0 );
//      break;
    }else{

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
    }
//    SHPDestroyObject(psShape);
//  }

  lonRef=(psShape->dfXMax+psShape->dfXMin)/2.0;
  for (int j=0; j<psShape->nVertices; ++j){
    double tmpLon=proj(psShape->padfX[j],psShape->padfY[j],lonRef);
    if (j==0 || tmpLon<lonMin) lonMin=tmpLon;
    if (j==0 || tmpLon>lonMax) lonMax=tmpLon;
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


int emain()
{
  ssimstate st;
  epregister2(st.R0,"r0");

  int avgage=0; // take most age of first individual (usually the oldest), or the avg age of the household for the household age group
  epregister(avgage);

  double finter=1.0;
  double fintra=1.0;

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

  double tmax=300.0;
  epregister(tmax);
  epregister(finter);
  epregister(fintra);
  double fglobal=0.01;
  epregister(fglobal);
  estr fpop="data/switzerland.agegroups";
  epregister(fpop);
  int fseed=10;
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

  loadShape("data/popdensmaps/gadm36_CHE/gadm36_CHE_0.shp");
//  loadPopDens("data/popdensmaps/gpw_v4_population_density_rev11_2020_2pt5_min.tif");

  //TODO: cleanup the dependecies to make them explicit between the initializations and loading of data
  uint8_t *frameraw = new uint8_t[1920*1080*4];
  initCairo(frameraw,1920,1080); // loadPopDens has to be after initCairo, initCairo after loadShape because it initializes a shapemask

  loadPopDens("data/popdensmaps/gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif",popsize); // loads population counts and adjusts using more recent total census data

  



/*
  edoublearray householdSizeDist(0.0,0.31,0.34,0.16,0.13,0.06);
  euintarray householdCounts;
  householdCounts.init(6,0); // how many households with single individuals, two individuals, ... max 6
  multinomial(popsize,householdSizeDist,householdCounts,rnd);
*/

  // Ferguson et al 2006 has household distribution sizes in supplementary graphs which the values below approximate
  // UK statistics have updated information on households
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

  double totcritical=0.0,totdeaths=0.0;
  for (int i=0; i<agegroups.size(); ++i){
    totcritical+=agegroups[i]*st.rSA*st.icu_symp[i];
    totdeaths+=agegroups[i]*st.rSA*(st.icu_symp[i]*st.death_icu[i]+st.nonicu_symp[i]*st.death_nonicu[i]);
  }
  double tmpcritical=0.0,tmpdeaths=0.0;
  for (int i=0; i<agegroups.size(); ++i){
    tmpcritical+=agegroups[i]*st.rSA*st.icu_symp[i];
    tmpdeaths+=agegroups[i]*st.rSA*(st.icu_symp[i]*st.death_icu[i]+st.nonicu_symp[i]*st.death_nonicu[i]);
    cout << "# age: " << 5*i << " critical: " << tmpcritical << " (" << tmpcritical/totcritical << ")" << " deaths: " << tmpdeaths << " (" << tmpdeaths/totdeaths << ")" <<endl;
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

  st.rToicu=rgamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.rTononicu=rgamma(7.0,7.0,60.0,0.25); // 7 days to go to ICU
  st.rInicu=rgamma(10.0,10.0,60.0,0.25); // 10 days in ICU
  st.rInnonicu=rgamma(8.0,8.0,60.0,0.25); // 8 days in nonICU
  st.rIpDeath=rgamma(22.0,22.0,60.0,0.25); // deaths occuring 22 days after symptoms

  double R0day=st.R0/5.0; // per day rate given that individuals are infectious for 5 days


  st.hhLevels2.init(total_households);
  ldieif(householdSizeDist.size()!=maxhousesize,"household size dist must be the same as maxhouse size");

  cout << "# total households: " << total_households << endl;
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
        st.pop_ages[sh.ageind]=int(5+(i-2)*0.5)+rnd.uniformint((tmpag.size()-int(5+(i-2)*0.5)-MAX(0,2*(i-2))));  // maxint is needed when the expression contains randomly generated numbers. Using the MAX macros causes the number to be generated twice!!
        ldieif(st.pop_ages[sh.ageind]<0 || st.pop_ages[sh.ageind]>=tmpag.size(),"out of bounds: "+estr().sprintf("%hhi",st.pop_ages[sh.ageind])+" "+tmpag.size());
//        st.pop_ages[sh.ageind]=5+(i-1)+rnd.uniform()*(tmpag.size()-5-(i-1)-MAX(0,2*(i-4)));
      } while (tmpag[st.pop_ages[sh.ageind]]==0);
      --tmpag[st.pop_ages[sh.ageind]];
      ++assigned;

//      // use single agegroup for hhLevels
//      sh.hage=0; //st.pop_ages[sh.ageind];

      // age group of "oldest" person in house
      sh.hage=st.pop_ages[sh.ageind];
      sh.group=(sh.hage>=10?1:0); // set to 1 or 0 depending on age

      if (i>1){
        do {
          ldieif(sh.ageind+1<0 || sh.ageind+1>=st.pop_ages.size(),"out of bounds: "+estr(sh.ageind+1)+" "+st.pop_ages.size());
          st.pop_ages[sh.ageind+1]=minint(st.pop_ages[sh.ageind]+int(rnd.uniformint(3))-2,tmpag.size()-1);
          ldieif(st.pop_ages[sh.ageind+1]<0 || st.pop_ages[sh.ageind+1]>=tmpag.size(),"out of bounds: "+estr().sprintf("%hhi",st.pop_ages[sh.ageind+1])+" "+tmpag.size());
        } while (tmpag[st.pop_ages[sh.ageind+1]]==0);
        --tmpag[st.pop_ages[sh.ageind+1]];
        if (avgage)
          sh.hage=(sh.hage+st.pop_ages[sh.ageind+1])/2;
        ++assigned;
      }   
      // Children
      for (int l=2; l<i; ++l){
        do {
          if (tmpag[MAX(0,st.pop_ages[sh.ageind]-9)]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+1]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+2]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+3]==0)
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-9) + 4 + int(rnd.uniformint(3));
          else if (tmpag[MAX(0,st.pop_ages[sh.ageind]-9)]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+1]+tmpag[MAX(0,st.pop_ages[sh.ageind]-9)+2]==0)
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-9) + 3;
          else
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-9) + int(rnd.uniformint(3));
/*
          if (tmpag[0]+tmpag[1]+tmpag[2]+tmpag[3]+tmpag[4]==0)
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)+1);
          else
            st.pop_ages[sh.ageind+l]=MAX(0,st.pop_ages[sh.ageind]-6+int(rnd.uniform()*3)-2);
*/
        } while (tmpag[st.pop_ages[sh.ageind+l]]==0);
        --tmpag[st.pop_ages[sh.ageind+l]];
        ++assigned;
        if (avgage)
          sh.hage=(sh.hage*l+st.pop_ages[sh.ageind+l])/(l+1);
      }
    }
  }
  cout << "# done populating households: " << assigned << endl;


  st.spGridW=rwidth;
  st.spGridH=rheight;
  st.spGridSize=st.spGridW*st.spGridH;
  cout << "# initializing spatial grid: " << st.spGridW << "x" << st.spGridH << endl;
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

  // randomly place households across grid
  for (int i=0; i<st.households.size(); ++i){
    shousehold &hh(st.households[i]);

    hh.gridpos=rnd.uniformint(st.spGrid.size());
    st.spGrid[hh.gridpos].N[hh.group]+=hh.size;
  }

  cout << "# initializing levels" << endl;


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

  cout << "# initializing levels indices" << endl;
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

  cout << "# populating levels indices" << endl;
  // populate hhLevels
  for (int i=0; i<st.households.size(); ++i){
    shousehold &sh(st.households[i]);
    int hg=sh.group;
    uint8_t hl=sh.size*(sh.size+1)/2+sh.size;

    // store household indice in levels array
    st.hhLevels2[st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl+1]]=i;
    sh.hstind=st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl+1]-st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl];
    ++st.hhLevelsBegin[sh.gridpos*st.grow + hg*st.arow + hl+1]; // end of subarray is beginning of the next level

    uint8_t *agearr=&st.pop_ages[sh.ageind];
    if (hg>=1){
      ++hhOlder;
      hhOlderN+=sh.size;
      for (int l=0; l<sh.size; ++l){
        if (agearr[l]>=10)
          ++hhOlderNO;
        else
          ++hhOlderNY;
      }
    }else{
      for (int l=0; l<sh.size; ++l){
        if (agearr[l]>=10)
          ++hhYoungNO;
        else
          ++hhYoungNY;
      }
    }
  }


  cout << "# hhOlder: " << hhOlder << " (" << double(hhOlder)/total_households << ")" << endl;
  cout << "# hhOlderN: " << hhOlderN << " (" << double(hhOlderN)/popsize << ")" << endl;
  cout << "# hhOlderNO: " << hhOlderNO << " (" << double(hhOlderNO)/popsize << ")" << endl;
  cout << "# hhOlderNY: " << hhOlderNY << " (" << double(hhOlderNY)/popsize << ")" << endl;
  cout << "# hhYoungNO: " << hhYoungNO << " (" << double(hhYoungNO)/popsize << ")" << endl;
  cout << "# hhYoungNY: " << hhYoungNY << " (" << double(hhYoungNY)/popsize << ")" << endl;



  // randomize ages in households since infection always progresses sequentially in the age array
  for (int i=0; i<st.households.size(); ++i){
    shousehold &hh(st.households[i]);
    if (hh.size==1) continue; // cannot randomize single households
    uint8_t hl=hh.size*(hh.size+1)/2+hh.size;
    permute(st.pop_ages,hh.ageind,hh.size,rnd);
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

  st.evqueue.resize(st.dE.size()+1); // keep 1 extra for using as buffer for next step when adding new events

  edoublearray mp;
  euintarray counts;

  double fE=0.0; // fraction of exposed (includes deaths, active infections and healed)

  cout << "# seeding initial infected: " << fseed << endl;
  // seed infections
  counts.init(3,0);
  counts[1]=fseed;

//  st.allIp=fseed; // directly as presymptomatic
//  st.evqueue.add(st,6,2,0,counts,st.dIp,rnd,2); 

  st.allE=fseed;
  cout << "hhLevelsBegin: " << st.hhLevelsBegin[1010*st.grow + 0*st.arow + 2*(2+1)/2+2+1] << " " << st.hhLevelsBegin[1010*st.grow + 0*st.arow + 2*(2+1)/2+2] << endl;
//  st.evqueue.add(st,1010,0,2,0,counts,st.dE,rnd);  // add Exposed grid position 10,10

  cout << "# starting simulation" << endl;

  st.allICU=0.0;
  st.allNonICU=0.0;
  st.allD=0.0;

  double tpeak=0.0;
  double peakIs=0.0;
  double peakICU=0.0;
  double peakNonICU=0.0;

  double soldr=1.0;
  double smr=1.0;

  ldieif(videoOpen()!=0,"error creating video file");

  edoublearray localIprob;
  localIprob.init(st.spGrid.size());

  cout << "Time" << "\t" << "fE" << "\t" << "allE" << "\t" << "allIa" << "\t" << "allIp" << "\t" << "allIs" << "\t" << "allICU" << "\t" << "allNonICU" << "\t" << "Deaths" << endl;
  for (int it=0; it*st.tstep<tmax; ++it){
    renderFrame(it*st.tstep,smr,frameraw,st,1920,1080);

    if (tmstart>=0 && it*st.tstep >= tmstart){
      smr=fmr;
      tmstart=-1;
    }

    if (tmend>0 && it*st.tstep >= tmend){
      smr=1.0;
      tmend=-1;
    }

    if (tostart>=0 && it*st.tstep >= tostart){
      soldr=foldr;
      tostart=-1;
    }

    if (toend>0 && it*st.tstep >= toend){
      soldr=1.0;
      toend=-1;
    }

//    double interhhrate=(hage>=10?folder:1.0)*finter*st.R0*(1.0-fE)*(0.5*st.allIa+st.allIp+st.allIs)/(popsize-st.allE);
    if (peakIs < st.Is[0]+st.Is[1]){
      peakIs=st.Is[0]+st.Is[1];
      tpeak=it*st.tstep;
    }
    if (peakICU < st.allICU)
      peakICU=st.allICU;
    if (peakNonICU < st.allNonICU)
      peakNonICU=st.allNonICU;

    // The hhOlderN term is needed to prevent R0 from decreasing when older population is isolated from younger population
    // Whithout including this term, the simulation is equivalent to an older population considered immune but still interacting, thus reducing incorrectly the R0

    // Global probability of infection from random person in country (x0.01)
//    double globalIprob=0.01 * finter*R0day*(0.5*st.allIa+st.allIp+st.allIs - (1.0-MIN(smr,soldr)/smr)*(0.5*st.oallIa+st.oallIp+st.oallIs) )/(popsize-(1.0-MIN(smr,soldr)/smr)*hhOlderN);
    double globalIprob=finter*fglobal*R0day*(0.5*st.Ia[0]+st.Ip[0]+st.Is[0] + (MIN(smr,soldr)/smr)*(0.5*st.Ia[1]+st.Ip[1]+st.Is[1]) )/(st.N[0]+(MIN(smr,soldr)/smr)*st.N[1]);

    // removed the cap because its more conservative and because there are situations where there should be no cap. For a more accurate representation of reality need to model how transmission occurs
    // through contacts, considering how time spent and number of contacts influence the transmission probability, and what is the effect of the isolation measures on either of these two variables

//    if (globalIprob>1.0) globalIprob=1.0; // The reason to cap the probability is to prevent the unrealistic situation where "isolation" has no effect on the probability of infection. Need to consider more carefully  what capping this value implies in the model and what the wider implications are
    // The cap means that the maximum probability of getting infected per day cannot be larger than 1, and because of this a mitigation factor of 50% will
    // guarantee a reduction in 50% probability of getting infected. This however implies there is a saturation of the probability of infection, 
    // To put this more mechanistically: if we assume isolation reduces the time a person spends with other contacts or number of contacts (with constant time per contact), and we assume the probability of infection is proportional to the time exposed to an infected individual. Then isolation will always have an effect on probability of infection.
    // There are cases when this is not true, if the probability of infection is so high that it is enough to have a single contact even of very short duration to have a probability of 100% of getting infected and the proportion of infected is high enough to guarantee there is one infected person in almost every group of contacts, then reducing contact (either by decreasing the time or by decreasing the number of contacts) should correctly result in no reduced probability of infection, up to a certain point.

    fE=double(st.allE)/popsize;
    for (int gi=0; gi<st.spGrid.size(); ++gi){
      sgrid& g(st.spGrid[gi]);
      // TODO: precompute all infection probabilities before updating data!
      localIprob[gi]=finter*R0day*(0.5*g.Ia[0]+g.Ip[0]+g.Is[0] + (MIN(smr,soldr)/smr)*(0.5*g.Ia[1]+g.Ip[1]+g.Is[1]) )/(g.N[0]+(MIN(smr,soldr)/smr)*g.N[1]);
    }

    for (int gi=0; gi<st.spGrid.size(); ++gi){
      sgrid& g(st.spGrid[gi]);
      int gx=gi%st.spGridW;
      int gy=gi/st.spGridW;

      int giup=((st.spGridH+gy-1)%st.spGridH)*st.spGridW+gx;
      int gidn=((gy+1)%st.spGridH)*st.spGridW+gx;
      int gilt=gy*st.spGridW+(st.spGridW+gx-1)%st.spGridW;
      int girt=gy*st.spGridW+(gx+1)%st.spGridW;

      double tmpLocalIprob=localIprob[gi]*0.6 + localIprob[gilt]*0.1 + localIprob[girt]*0.1 + localIprob[giup]*0.1 + localIprob[gidn]*0.1;

      // TODO: precompute all infection probabilities before updating data!
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
    deque<sevent> &nextEvents(st.evqueue.step());
    processEvents(st,nextEvents,rnd);
    cout << it*st.tstep << "\t" << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.Ia[0]+st.Ia[1] << "\t" << st.Ip[0]+st.Ip[1] << "\t" << st.Is[0]+st.Is[1] << "\t" << st.allICU << "\t" << st.allNonICU << "\t" << st.allD << endl;
//    if (st.Ia[0]+st.Ia[1]+st.Is[0]+st.Is[1]+st.Ip[0]+st.Ip[1]==0) break;
  }
  videoClose();
  cout << double(st.allE)/popsize << "\t" << st.allE << "\t" << st.Ia[0]+st.Ia[1] << "\t" << st.Ip[0]+st.Ip[1] << "\t" << st.Is[0]+st.Is[1] << "\t" << st.allICU << "\t" << st.allNonICU << "\t" << st.allD << endl;

  cout << "# fE: " << double(st.allE)/popsize << " tpeak: " << tpeak << " peakIs: " << peakIs << " peakICU: " << peakICU << " peakNonICU: " << peakNonICU << " deaths: " << st.allD << endl;

  return(0);
}
