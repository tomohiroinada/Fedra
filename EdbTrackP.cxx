//-- Author :  Valeri Tioukov   19.05.2002
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbPattern                                                           //
//                                                                      //
// Segments pattern                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TIndexCell.h"
#include "TPolyLine.h"
#include "EdbAffine.h"
#include "EdbPattern.h"
#include "EdbVertex.h"
#include "EdbPhys.h"
#include "EdbLog.h"
#include "vt++/CMatrix.hh"
#include "vt++/VtVector.hh"

ClassImp(EdbSegmentsBox)
ClassImp(EdbTrackP)
ClassImp(EdbPattern)
ClassImp(EdbPatternsVolume)

using namespace MATRIX;
using namespace TMath;

//______________________________________________________________________________
EdbSegmentsBox::EdbSegmentsBox(int nseg)
{
  if(nseg>0) eSegments = new TClonesArray("EdbSegP",nseg);
  else       eSegments = new TClonesArray("EdbSegP");
  Set0();
}
 
//______________________________________________________________________________
EdbSegmentsBox::EdbSegmentsBox(float x0, float y0, float z0, int nseg)
{
  if(nseg>0) eSegments = new TClonesArray("EdbSegP",nseg);
  else       eSegments = new TClonesArray("EdbSegP");
  eX=x0;  eY=y0;  eZ=z0;
  eDZkeep=0;
}
 
//______________________________________________________________________________
EdbSegmentsBox::~EdbSegmentsBox( )
{
  if(eSegments) {
    eSegments->Delete();
    delete eSegments;
    eSegments=0;
  }
}
 
//______________________________________________________________________________
void EdbSegmentsBox::Set0()
{
  eX=0;  eY=0;  eZ=0;
  eDZkeep=0;
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment(int id, float x, float y, float tx, float ty, 
			    float w, int flag)
{
  return new((*eSegments)[N()])  EdbSegP( id,x,y,tx,ty,w,flag );
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment( int i, EdbSegP &s )
{
  return new((*eSegments)[i])  EdbSegP( s );
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment( EdbSegP &s )
{
  return new((*eSegments)[N()])  EdbSegP( s );
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment( EdbSegP &s1, EdbSegP &s2 )
{
  EdbSegP *s = new((*eSegments)[N()])  EdbSegP( s1 );
  s->MergeTo(s2);
  return s;
}
 
//______________________________________________________________________________
void EdbSegmentsBox::Reset()
{
  Set0();
  if(eSegments) eSegments->Clear();
}
 
//______________________________________________________________________________
float EdbSegmentsBox::DiffAff(EdbAffine2D *aff)
{
  EdbSegmentsBox p,p1;
  p.AddSegment(0,Xmin(),Ymin(),0.,0.);
  p.AddSegment(0,Xmax(),Ymin(),0.,0.);
  p.AddSegment(0,Xmin(),Ymax(),0.,0.);
  p.AddSegment(0,Xmax(),Ymax(),0.,0.);

  p1.AddSegment(0,Xmin(),Ymin(),0.,0.);
  p1.AddSegment(0,Xmax(),Ymin(),0.,0.);
  p1.AddSegment(0,Xmin(),Ymax(),0.,0.);
  p1.AddSegment(0,Xmax(),Ymax(),0.,0.);

  p1.Transform(aff);
  return p.Diff(p1);
}
 
//______________________________________________________________________________
float EdbSegmentsBox::Diff(EdbSegmentsBox &p)
{
  // return the mean difference beteween pattern elements
  int nseg = TMath::Min( N(), p.N() );

  if(nseg<1) return 0;

  EdbSegP *s1=0, *s2=0;

  float dx=0, dy=0;
  double sdx=0;
  for(int i=0; i<nseg; i++ ) {
    s1 =   GetSegment(i);
    s2 = p.GetSegment(i);
    dx = s2->X() - s1->X();
    dy = s2->Y() - s1->Y();
    sdx += TMath::Sqrt( dx*dx + dy*dy);
  }
  return sdx/nseg;
}

//______________________________________________________________________________
float EdbSegmentsBox::GetSize(Int_t XorY)
{
  // return the size of the segments box in X//Y
  // simple approach using max/min values of segments
  // works only if 2 or more segemts are in it.
  int nseg = N();
  if(nseg<2) return 0;
  EdbSegP *s1=0;

  float minX=0; float maxX=0;
  float minY=0; float maxY=0;
  for(int i=0; i<nseg; i++ ) {
    s1 =   GetSegment(i);
    if (i==0) { 
      minX=s1->X(); maxX=s1->X(); 
      minY=s1->Y(); maxY=s1->Y();
    }
    maxX=TMath::Max(maxX,s1->X());
    minX=TMath::Min(minX,s1->X());
    maxY=TMath::Max(maxY,s1->Y());
    minY=TMath::Min(minY,s1->Y());
  }
  float sizeX=TMath::Abs(maxX-minX);
  float sizeY=TMath::Abs(maxY-minY);
  
  if (XorY==1) { return sizeY; }
  else { return sizeX; }
}

//______________________________________________________________________________
float EdbSegmentsBox::GetSizeXY()
{
  return (GetSize(0)*GetSize(1));
}
 
//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsPlate(int plate)
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetPlate( plate );
  }
}
 
//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsZ()
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetZ( Z() );
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsDZ(float dz)
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetDZ( dz );
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::ProjectTo(const float dz)
{
  eZ += dz;  eDZkeep += dz;
  int nseg = N();
  for(int i=0; i<nseg; i++ ) GetSegment(i)->PropagateToDZ(dz);
/*  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);
    p->SetX( p->X() + p->TX()*dz );
    p->SetY( p->Y() + p->TY()*dz );
    p->SetZ( Z() );
  }*/
}
 
//______________________________________________________________________________
int EdbSegmentsBox::CalculateXY( EdbSegmentsBox *pat, EdbAffine2D *aff )
{
  int n=N();
  if( n>pat->N() )  n=pat->N();
  if(n<2) return 0;

  aff->Calculate(this,pat);
  return 1;
}

//______________________________________________________________________________
int EdbSegmentsBox::CalculateAXAY( EdbSegmentsBox *pat, EdbAffine2D *aff )
{
  int n=N();
  if( n>pat->N() )  n=pat->N();
  if(n<2) return 0;

  float *ax1 = new float[n];
  float *ay1 = new float[n];
  float *ax2 = new float[n];
  float *ay2 = new float[n];

  EdbSegP *p1,*p2;
  for(int i=0; i<n; i++ ) {
    p1 = GetSegment(i);
    p2 = pat->GetSegment(i);
    ax1[i] = p1->TX();
    ay1[i] = p1->TY();
    ax2[i] = p2->TX();
    ay2[i] = p2->TY();
  }
  aff->Calculate(n,ax1,ay1,ax2,ay2);

  delete ax1;
  delete ay1;
  delete ax2;
  delete ay2;

  return 1;
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformA( const EdbAffine2D *aff )
{
  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = aff->A11()*p->TX() + aff->A12()*p->TY() + aff->B1();
    ty = aff->A21()*p->TX() + aff->A22()*p->TY() + aff->B2();
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformShr( float shr )
{
  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = p->TX()/shr;
    ty = p->TY()/shr;
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformARot( const EdbAffine2D *aff )
{
  // apply to the angles only rotation members of transformation

  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = aff->A11()*p->TX() + aff->A12()*p->TY();
    ty = aff->A21()*p->TX() + aff->A22()*p->TY();
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::Print(Option_t *opt) const
{
  int nseg=GetN();
  printf("EdbSegmentsBox: %d segments\n", nseg );
  for(int i=0; i<nseg; i++) GetSegment(i)->Print();
} 

//______________________________________________________________________________
//______________________________________________________________________________
EdbTrackP::EdbTrackP(int nseg)
{
  Set0();
  if(nseg>0) eS  = new TSortedList();
  if(nseg>0) { eSF = new TSortedList();    eSF->SetOwner(); }
}

//______________________________________________________________________________
EdbTrackP::EdbTrackP(EdbSegP &seg) : EdbSegP(seg)
{
  //create empty track with segment parameters
  eS=0;
  eSF=0;
  eM=0;
  eDE=0;
  ePDG=-999;
  eVTAS = 0;
  eVTAE = 0;
  eNpl=0;
  eN0=0;
}
 
//______________________________________________________________________________
EdbTrackP::EdbTrackP(EdbSegP *seg, float m) : EdbSegP( *seg )
{
  eS=0;
  eSF=0;
  eM=0;
  eDE=0;
  ePDG=-999;
  eVTAS = 0;
  eVTAE = 0;
  eNpl=0;
  eN0=0;
  AddSegment(seg);
  SetM(m);
}
 
//______________________________________________________________________________
EdbTrackP::~EdbTrackP()
{
  if(eS)    { eS->Clear();  delete eS;  eS=0;  }
  if(eSF)   { eSF->Clear(); delete eSF; eSF=0; }
}

//______________________________________________________________________________
void EdbTrackP::Set0()
{
  ((EdbSegP*)this)->Set0();
  eS=0;
  eSF=0;
  eM=0;
  eDE=0;
  ePDG=-999;
  eVTAS = 0;
  eVTAE = 0;
  eNpl=0;
  eN0=0;
  ePerrUp=0;
  ePerrDown=0;
  return;
}

//______________________________________________________________________________
void EdbTrackP::AddVTA(EdbVTA *vta)
{
  if(vta) {
    if(vta->Zpos()==0)      eVTAE=vta;
    else if(vta->Zpos()==1) eVTAS=vta;
  }
}

//______________________________________________________________________________
void EdbTrackP::ClearVTA()
{
  eVTAE=0;
  eVTAS=0;
}

//______________________________________________________________________________
void EdbTrackP::ClearVTA(EdbVTA *vta)
{
  if(vta->Zpos()==0)      eVTAE=0;
  else if(vta->Zpos()==1) eVTAS=0;
}

//______________________________________________________________________________
void EdbTrackP::Copy(const EdbTrackP &tr)
{
  // do the physical copy of segments
  Clear();
  ((EdbSegP*)(this))->Copy( *((EdbSegP*)(&tr)) );
  SetM(tr.M());
  SetPDG(tr.PDG());
  SetNpl(tr.Npl());  
  SetN0(tr.N0());
  SetDE(tr.DE());
  AddVTA(tr.VTAS());
  AddVTA(tr.VTAE());

  int nseg=tr.N();
  for(int i=0; i<nseg; i++)
    AddSegment(new EdbSegP(*tr.GetSegment(i)));
  nseg=tr.NF();
  for(int i=0; i<nseg; i++)
    AddSegmentF(new EdbSegP(*tr.GetSegmentF(i)));
  if (eS) eS->SetOwner();
  if (eSF) eSF->SetOwner();
}

//______________________________________________________________________________
void EdbTrackP::Transform(const EdbAffine2D &aff)
{
  // apply transformation to all track elements
  ((EdbSegP*)(this))->Transform(&aff);
  for(int i=0; i<N(); i++)  GetSegment(i)->Transform(&aff);
  for(int i=0; i<NF(); i++) GetSegmentF(i)->Transform(&aff);
}

//______________________________________________________________________________
int EdbTrackP::CheckAliasSegments()
{
  int nalias=0;
  for(int i=0; i<N(); i++) 
    if( GetSegment(i)->Track() != ID()) nalias++;
  return nalias;
}

//______________________________________________________________________________
int EdbTrackP::RemoveAliasSegments()
{
  int nalias=0;
  EdbSegP *s=0;
  for(int i=0; i<N(); i++) {
    s = GetSegment(i);
    if( s->Track() != ID()) { 
      nalias++;
      RemoveSegment(s);
    }
  }
  return nalias;
}

//______________________________________________________________________________
int EdbTrackP::CheckMaxGap()
{
  int ngap=0;
  if(N()<2) return 0;
  EdbSegP *s1=0, *s2=0;
  int gap = 0;
  for(int i=0; i<N()-1; i++) {
    s1 = GetSegment(i);
    s2 = GetSegment(i+1);
    gap = TMath::Abs(s2->PID()-s1->PID());
    if( gap > ngap) ngap = gap;
  }
  return ngap;
}
//______________________________________________________________________________
int EdbTrackP::GetSegmentsMCTrack( int& nsegf ) const
{
  // return: ID of the MC track contributed to this reconstructed track 
  //         with the maximal number of segments
  // nsegf - return number of segments originated from the selected MC track

  int nseg = N();
  if (nseg < 2) return -1;
  if (nseg > 300) nseg = 300;
  int count[300];
  int i, f = 0;
  for(i=0; i<nseg; i++)
  {
    count[i] = 0;
    f = GetSegment(i)->MCTrack();
    if (f < 0) continue;
    for (int j=0; j<nseg; j++)
    {
	if ( f == GetSegment(j)->MCTrack()) count[i]++;
    }
  }
  int cmax = 0;
  int mctrackmax = -1;
  for(i=0; i<nseg; i++)
  {
	if (count[i] > cmax)
	{
	    cmax = count[i];
	    mctrackmax = GetSegment(i)->MCTrack();
	}
  }
  nsegf = cmax;
  return mctrackmax;
}

//______________________________________________________________________________
Float_t EdbTrackP::Wmean() const
{
  int n = N();
  double wtot=0;
  for(int i=0; i<n; i++) wtot+=GetSegment(i)->W();
  wtot/=n;
  return (Float_t)wtot;
}

//______________________________________________________________________________
int EdbTrackP::GetSegmentsFlag( int& nsegf ) const
{
  printf(" EdbTrackP::GetSegmentsFlag: obsolete, to be deleted!  now replaced by GetSegmentsMCTrack\n");

  int nseg = N();
  if (nseg < 2) return -1;
  if (nseg > 300) nseg = 300;
  int count[300];
  int i, f = 0;
  for(i=0; i<nseg; i++)
  {
    count[i] = 0;
    f = GetSegment(i)->Flag();
    if (f < 0) continue;
    for (int j=0; j<nseg; j++)
    {
	if ( f == GetSegment(j)->Flag()) count[i]++;
    }
  }
  int cmax = 0;
  int flagmax = -1;
  for(i=0; i<nseg; i++)
  {
	if (count[i] > cmax)
	{
	    cmax = count[i];
	    flagmax = GetSegment(i)->Flag();
	}
  }
  nsegf = cmax;
  return flagmax;
}
//______________________________________________________________________________
int EdbTrackP::GetSegmentsAid( int& nsegf ) const
{
  int nseg = N();
  if (nseg < 2) return -1;
  if (nseg > 300) nseg = 300;
  int count[300];
  int aid[300], ai;
  int i, f = 0;
  for(i=0; i<nseg; i++)
  {
    count[i] = 0;
    aid[i] = -1;
    ai = GetSegment(i)->Aid(1);
    f = ai + (GetSegment(i)->Aid(0))*1000;
    if (f < 0) continue;
    aid[i] = ai;
    for (int j=0; j<nseg; j++)
    {
	if ( f == (GetSegment(j)->Aid(1) + (GetSegment(i)->Aid(0))*1000)) count[i]++;
    }
  }
  int cmax = 0;
  int aidmax = -1;
  for(i=0; i<nseg; i++)
  {
	if (count[i] > cmax)
	{
	    cmax = count[i];
	    aidmax = aid[i];
	}
  }
  nsegf = cmax;
  return aidmax;
}
//______________________________________________________________________________
float EdbTrackP::Wgrains() const
{
  float w=0.;
  int nseg=N();
  for(int i=0; i<nseg; i++)    w+=GetSegment(i)->W();
  return w;
}

//______________________________________________________________________________
void EdbTrackP::AddTrack(const EdbTrackP &tr)
{
  int nseg=tr.N();
  int nsegf=tr.NF();
  for(int i=0; i<nseg; i++)
    AddSegment(tr.GetSegment(i));
  for(int i=0; i<nsegf; i++)
    AddSegmentF(new EdbSegP(*(tr.GetSegmentF(i))));
}

//______________________________________________________________________________
EdbSegP *EdbTrackP::GetSegmentWithClosestZ( float z, float dzMax )
{
  float dzmin=dzMax;
  EdbSegP *sbest=0;
  int nseg=N();
  for(int i=0; i<nseg; i++) {
     EdbSegP *s = GetSegment(i);
     float dz = Abs(s->eZ-z);
     if( dz < dzmin ) {
       dzmin=dz;
       sbest=s;
     }
   }
  return sbest;
}

//______________________________________________________________________________
void EdbTrackP::FitTrack()
{
  // track fit by averaging of segments parameters

  int nseg=N();
  float x=0,y=0,z=0,tx=0,ty=0,w=0;
  EdbSegP *seg=0;
  for(int i=0; i<nseg; i++) {
    seg = GetSegment(i);
    x  += seg->X();
    y  += seg->Y();
    z  += seg->Z();
    tx += seg->TX();
    ty += seg->TY();
    w  += seg->W();
  }
  x  /= nseg;
  y  /= nseg;
  z  /= nseg;
  tx /= nseg;
  ty /= nseg;
  Set(ID(),x,y,tx,ty,w,Flag());
  SetZ(z);
}

//______________________________________________________________________________
int  EdbTrackP::FitTrackKFS( bool zmax, float X0, int design )
{
  // if (zmax==true)  track parameters are calculated at segment with max Z
  // if (zmax==false) track parameters are calculated at segment with min Z
  // TODO: track parameters??
  //
  // Note: momentum for track must be defined with SetP() prior to call this routine!
  // Note: momentum dispersion for track should be defined with SetErrorP()
  //	   prior to call this routine - it is necessary for vertex fit
  //
  // on the output: Chi2: the full chi-square (not divided by NDF); NDF = 4
  //                Prob: is Chi2 probability (area of the tail of Chi2-distribution)
  //                      If we accept events with Prob >= ProbMin then ProbMin is the 
  //                      probability to reject the good event


  int nseg=N();
  if(nseg<=0) return 0;
  if(NF()) {
    eSF->Delete();
  }

  if(SP()<0.000001) SetErrorP(1.);    // TODO: razobratsia s etimi impul'sami!
  //printf("%d segments to fit\n",N());

  //TODO - eliminate constants!!!

  float dPb = 0., p = 0., m = 0.13957, e = 0.13957, de = 0., pa = 0., pn = 0.;
  float e0 = e;
  float eTPb = 1000./1300.;
  float pcut = 0.050;
  double teta0sq;
  double dz, ptx, pty;
  int step;
  int istart, iend;

  VtVector *par[260], *parpred[260], *pars[260], *meas[260];
  VtSqMatrix *pred[260];
  VtSymMatrix *cov[260], *covpred[260], *covpredinv[260], *covs[260], *dmeas[260];
 
  int i=0;
  if(nseg == 1) {
    EdbSegP *s = GetSegment(0);
    EdbSegP segf(*s);

    segf.Set(s->ID(),s->X(),s->Y(),s->TX(),s->TY(),1.,s->Flag());
    segf.SetZ(s->Z());
    segf.SetErrorP ( SP() );
    segf.SetChi2(0.);
    segf.SetProb(1.);
    segf.SetP( P() );
    segf.SetPID( s->PID() );
    segf.SetDZ(s->DZ());

    AddSegmentF(new EdbSegP(segf));
    return 0;
  }

  if(nseg>259)   return -1;

  EdbSegP *seg0=0;
  EdbSegP *seg=0;

  if(  GetSegment(N()-1)->Z() <  GetSegment(0)->Z() )
  {
    if (zmax)
    {
	step=-1;
	istart=nseg-1;
	iend=0;
    }
    else
    {
	step=1;
	istart=0;
	iend=nseg-1;
    }
  }
  else
  {
    if (!zmax)
    {
	step=-1;
	istart=nseg-1;
	iend=0;
    }
    else
    {
	step=1;
	istart=0;
	iend=nseg-1;
    }

  }

  seg0 = GetSegment(istart);

  par[istart] = new VtVector(  (double)(seg0->X()), 
			       (double)(seg0->Y()),  
			       (double)(seg0->TX()), 
			       (double)(seg0->TY()) );
  meas[istart] = new VtVector(*par[istart]);
  pred[iend] = new VtSqMatrix(4);
//  (*pred[istart]).clear();
//  (*pred[istart])(0,0) = 1.;
//  (*pred[istart])(1,1) = 1.;
//  (*pred[istart])(2,2) = 1.;
//  (*pred[istart])(3,3) = 1.;
  cov[istart] = new VtSymMatrix(4);             // covariance matrix for seg0
  for(int k=0; k<4; k++) 
    for(int l=0; l<4; l++) (*cov[istart])(k,l) = (seg0->COV())(k,l);
  dmeas[istart] = new VtSymMatrix(*cov[istart]);             // covariance matrix for seg0

  Double_t chi2=0.; 

  i = istart;
  p = P();
  m = M();
  if (p < pcut) p = pcut; 
  e = TMath::Sqrt((double)p*(double)p + (double)m*(double)m);
  e0 = e;
  while( (i+=step) != iend+step ) {

    seg = GetSegment(i);
 
    VtSymMatrix dms(4);   // multiple scattering matrix
    dms.clear();

    dz = seg->Z()-seg0->Z();
    ptx = (*par[i-step])(2);                        // previous tx
    pty = (*par[i-step])(3);                        // previous ty
    dPb = dz*TMath::Sqrt(1.+ptx*ptx+pty*pty); // thickness of the Pb+emulsion cell in microns
    if ((design != 0) && (p > pcut))
    {
	de = EdbPhysics::DeAveragePb(p, m, TMath::Abs(eTPb*dPb));
	if (de < 0.) de = 0.;
	if (design < 0) de = -de;
	if (dz >= 0.)
	    e = e - de;
	else
	    e = e + de;
	if (e < m) e = m;
	pn = TMath::Sqrt((double)e*(double)e - (double)m*(double)m);
	if (pn <= pcut) pn = pcut;
	pa = 0.5*(p + pn);
	p  = pn;
    }
    else
    {
	pa = p;
    }
    teta0sq = EdbPhysics::ThetaMS2( pa, m, dPb, X0 );
    dms(0,0) = teta0sq*dz*dz/3.;
    dms(1,1) = dms(0,0);
    dms(2,2) = teta0sq;
    dms(3,3) = dms(2,2);
//    dms(2,0) = teta0sq*TMath::Abs(dz)/2.;
    dms(2,0) = teta0sq*dz/2.;
    dms(3,1) = dms(2,0);
    dms(0,2) = dms(2,0);
    dms(1,3) = dms(2,0);

    pred[i-step] = new VtSqMatrix(4);        //propagation matrix for track parameters (x,y,tx,ty)
    pred[i-step]->clear();

    (*pred[i-step])(0,0) = 1.;
    (*pred[i-step])(1,1) = 1.;
    (*pred[i-step])(2,2) = 1.;
    (*pred[i-step])(3,3) = 1.;
    (*pred[i-step])(0,2) = dz;
    (*pred[i-step])(1,3) = dz;

    parpred[i] = new VtVector(4);            // prediction from seg0 to seg
    *parpred[i] = (*pred[i-step])*(*par[i-step]);

    covpred[i] = new VtSymMatrix(4);         // covariation matrix for prediction
    *covpred[i] = (*pred[i-step])*((*cov[i-step])*((*pred[i-step]).T()))+dms;

    dmeas[i] = new VtSymMatrix(4);           // original covariation  matrix for seg
    for(int k=0; k<4; k++) 
      for(int l=0; l<4; l++) (*dmeas[i])(k,l) = (seg->COV())(k,l);

    covpredinv[i] = new VtSymMatrix(4);
    (*covpredinv[i]) = (*covpred[i]).dsinv();
    VtSymMatrix dmeasinv(4);
    dmeasinv  = (*dmeas[i]).dsinv();
    cov[i] = new VtSymMatrix(4);
    (*cov[i]) = (*covpredinv[i]) + dmeasinv;
    (*cov[i]) = (*cov[i]).dsinv();

    meas[i] = new VtVector( (double)(seg->X()), 
			    (double)(seg->Y()),  
			    (double)(seg->TX()), 
			    (double)(seg->TY()) );

    par[i] = new VtVector(4);
    (*par[i]) = (*cov[i])*((*covpredinv[i])*(*parpred[i]) + dmeasinv*(*meas[i]));   // new parameters for seg

//    chi2 += ((*par[i])-(*parpred[i]))*((*covpredinv[i])*((*par[i])-(*parpred[i]))) + 
//	    ((*par[i])-(*meas[i]))*(dmeasinv*((*par[i])-(*meas[i])));

    VtSymMatrix dresid(4);
    dresid = (*dmeas[i]) - (*cov[i]);
    dresid = dresid.dsinv();

    chi2 += ((*par[i])-(*meas[i]))*(dresid*((*par[i])-(*meas[i])));

    seg0 = seg;
  }

  Set(ID(),(float)(*par[iend])(0),(float)(*par[iend])(1),
	   (float)(*par[iend])(2),(float)(*par[iend])(3),1.,Flag());
  SetZ(GetSegment(iend)->Z());
  SetPID(GetSegment(iend)->PID());
  SetCOV( (*cov[iend]).array(), 4 );

  //SetChi2((float)chi2);
  //SetProb( (float)TMath::Prob(chi2,nseg*4));

// Smoothing

  double chi2p=0;

  pars[iend] = new VtVector(*par[iend]);
  covs[iend] = new VtSymMatrix(*cov[iend]);
  VtSymMatrix dresid(4);
  dresid = (*dmeas[iend]) - (*covs[iend]);
  dresid = dresid.dsinv();
  chi2p = ((*pars[iend])-(*meas[iend]))*(dresid*((*pars[iend])-(*meas[iend])));

  EdbSegP segf;
  segf.Set(ID(),(float)(*pars[iend])(0),(float)(*pars[iend])(1),
	   (float)(*pars[iend])(2),(float)(*pars[iend])(3),1.,Flag());
  segf.SetZ(GetSegment(iend)->Z());
  segf.SetCOV( (*covs[iend]).array(), 4 );
  segf.SetErrorP ( SP() );
  segf.SetChi2((float)chi2p);
  segf.SetProb( (float)TMath::Prob(chi2p,4));
  segf.SetW( (float)nseg );
  segf.SetP( P() );
  segf.SetPID( GetSegment(iend)->PID() );
  segf.SetDZ(seg->DZ());
  
  AddSegmentF(new EdbSegP(segf));

  i=iend; 
  double DE=0.;
  EdbPhysics::DeAveragePbFastSet(P(), M());
   
  while( (i-=step) != istart-step ) {
	VtSqMatrix BackTr(4);
	BackTr = (*cov[i])*(((*pred[i]).T())*(*covpredinv[i+step]));
	pars[i] = new VtVector(4);
	covs[i] = new VtSymMatrix(4);
	(*pars[i]) = (*par[i]) + BackTr*((*pars[i+step])-(*parpred[i+step]));
	(*covs[i]) = (*cov[i]) + BackTr*(((*covs[i+step])-(*covpred[i+step]))*BackTr.T());
	dresid = (*dmeas[i]) - (*covs[i]);
	dresid = dresid.dsinv();
	chi2p = ((*pars[i])-(*meas[i]))*(dresid*((*pars[i])-(*meas[i])));
//	chi2 += chi2p;

	segf.Set(ID(),(float)(*pars[i])(0),(float)(*pars[i])(1),
		 (float)(*pars[i])(2),(float)(*pars[i])(3),1.,Flag());
	segf.SetZ(GetSegment(i)->Z());
	segf.SetCOV( (*covs[i]).array(), 4 );
	segf.SetErrorP ( SP() );
	segf.SetChi2((float)chi2p);
	segf.SetProb( (float)TMath::Prob(chi2p,4));
	segf.SetW( (float)nseg );
	segf.SetP( P() );
	segf.SetPID( GetSegment(i)->PID() );
	segf.SetDZ(seg->DZ());
	AddSegmentF(new EdbSegP(segf));
	dz = eTPb*TMath::Abs(GetSegment(i)->Z() - GetSegment(i+step)->Z()); 
        dPb = dz*TMath::Sqrt(1.+(*pars[i])(2)*(*pars[i])(2)+(*pars[i])(3)*(*pars[i])(3)); // thickness of the Pb+emulsion cell in microns
	DE += EdbPhysics::DeAveragePbFast(P(),M(),TMath::Abs(dPb));
  }
  SetChi2((float)chi2);
  SetProb((float)TMath::Prob(chi2,4*(nseg-1)));
  SetW( (float)nseg );
  SetDE( (float)DE );

//  DEBUG
//  printf(" e0 - e = %f, de = %f\n", e0-e, DE);

// Delete matrixes and vectors

  delete par[istart];
  par[istart] = 0;
  delete cov[istart];
  cov[istart] = 0;
  delete meas[istart];
  meas[istart] = 0;
  delete dmeas[istart];
  dmeas[istart] = 0;
  delete pred[istart];
  pred[istart] = 0;
  delete pars[istart];
  pars[istart] = 0;
  delete covs[istart];
  covs[istart] = 0;
  i=istart; 
  while( (i+=step) != iend+step ) {
    delete pred[i];
    pred[i] = 0;
    delete parpred[i];
    parpred[i] = 0;
    delete covpred[i];
    covpred[i] = 0;
    delete covpredinv[i];
    covpredinv[i] = 0;
    delete par[i];
    par[i] = 0;
    delete cov[i];
    cov[i] = 0;
    delete meas[i];
    meas[i] = 0;
    delete dmeas[i];
    dmeas[i] = 0;
    delete pars[i];
    pars[i] = 0;
    delete covs[i];
    covs[i] = 0;
  }
  return 0;
}

///______________________________________________________________________________
float   EdbTrackP::MakePredictionTo( Float_t z, EdbSegP &ss )
{
  float dz = Zmax()-Zmin();
  const EdbSegP *tr=0;
  if     ( z <= Zmin() )      tr = TrackZmin();   //TODO: fitted - not fitted
  else if( z >= Zmax() )      tr = TrackZmax();
  else {
    for(int i=0; i<N(); i++)       {
      EdbSegP *s = GetSegment(i);
      if( Abs(s->Z()-z) < dz ) { dz = Abs(s->Z()-z); tr=s; }     // select nearest segment (TODO: correct interpolation)
    }
  }
  if(!tr) // no segments in this track: use track body
    tr = this;
  dz = Abs(tr->Z()-z);
  ss.Copy(*tr);
  ss.PropagateTo(z);
  return dz;
}

///______________________________________________________________________________
int EdbTrackP::MakeSelector( EdbSegP &ss, bool followZ )
{
  if(N()<2) return 0;
  ss.SetCOV( GetSegment(0)->COV() );             // TODO ?
  const EdbSegP *tr;
  if (NF())
  {
    if( followZ ) tr = TrackZmax();
    else  tr = TrackZmin();
  }
  else
  {
    if( followZ ) tr = GetSegmentLast();
    else  tr = GetSegmentFirst();
  }
  ss.SetTX(tr->TX());
  ss.SetTY(tr->TY());
  ss.SetX(tr->X());
  ss.SetY(tr->Y());
  ss.SetZ(tr->Z());
  ss.SetPID(tr->PID());

  if( tr->PID() > GetSegment(0)->PID() )     return 1;
  if( tr->PID() > GetSegment(N()-1)->PID() ) return 1;
  return -1;
}

//______________________________________________________________________________
float EdbTrackP::CHI2()
{
  double dtx=0,dty=0,chi2=0;
  EdbSegP *seg=0;
  int    nseg=N();
  for(int i=0; i<nseg; i++) {
    seg = GetSegment(i);
    dtx = seg->TX()-TX();
    dty = seg->TY()-TY();
    chi2 += TMath::Sqrt( dtx*dtx/seg->STX() + 
			 dty*dty/seg->STY() );
  }
  chi2  /= nseg;
  return chi2;
}

//______________________________________________________________________________
float EdbTrackP::CHI2F()
{
  double dtx=0,dty=0,chi2=0;
  EdbSegP *s=0, *sf=0;
  int    nseg=N();
  for(int i=0; i<nseg; i++) {
    s  = GetSegment(i);
    sf = GetSegmentF(i);
    dtx = s->TX() - sf->TX();
    dty = s->TY() - sf->TY();
    chi2 += TMath::Sqrt( dtx*dtx/s->STX() + 
			 dty*dty/s->STY() );
  }
  chi2  /= nseg;
  return chi2;
}

//______________________________________________________________________________
void EdbTrackP::Print() 
{
  int nseg=0, nsegf=0;
  nseg = N();
  nsegf = NF();

  printf("EdbTrackP with %d segments and %d fitted segments\n", nseg, nsegf );
  printf("particle mass = %f\n", M() );
  ((EdbSegP*)this)->Print(); 

  if(nseg) 
    for(int i=0; i<nseg; i++)
      ((EdbSegP*)(eS->At(i)))->Print(); 

//    if(nsegf) 
//      for(int i=0; i<nsegf; i++) 
//        ((EdbSegP*)(eSF->At(i)))->Print();
}

//______________________________________________________________________________
void EdbTrackP::PrintNice() 
{
  int nseg=0, nsegf=0;
  nseg = N();
  nsegf = NF();

  printf("EdbTrackP with %d segments and %d fitted segments:\n", nseg, nsegf );
  printf("  N  PID     ID          X             Y              Z        TX       TY     W      P      Flag     MC     Track    Chi2    Prob     Mass\n");
  printf("    %3d %8d %13.2f %13.2f %13.2f %7.4f  %7.4f %5.1f  %7.2f %7d %7d %7d  %7.4f  %7.4f  %5.3f\n",
	   PID(), ID(),X(),Y(),Z(),     TX(),   TY(),  W(),  P(), Flag(), MCEvt(),   Track(), Chi2(), Prob(),   M());

  EdbSegP *s=0;
  if(nseg) 
    for(int i=0; i<nseg; i++) {
      s = GetSegment(i);
      printf("%3d %3d %8d %13.2f %13.2f %13.2f %7.4f  %7.4f %5.1f  %7.2f %7d %7d %7d  %7.4f  %7.4f\n",
	     i, s->PID(), s->ID(),s->X(),s->Y(),s->Z(),  s->TX(),s->TY(),  s->W(),  s->P(), s->Flag(), s->MCEvt(), s->Track(),s->Chi2(),s->Prob());
    }
}


//______________________________________________________________________________
EdbVertex  *EdbTrackP::VertexS()
{
    if(eVTAS) return eVTAS->GetVertex();
    return 0;
}

//______________________________________________________________________________
EdbVertex  *EdbTrackP::VertexE()
{
    if(eVTAE) return eVTAE->GetVertex();
    return 0;
}

//______________________________________________________________________________
Float_t EdbTrackP::Zmin() const
{
  Float_t zmin = Z(); 
  if (N() && GetSegmentFirst()->Z() < zmin) zmin = GetSegmentFirst()->Z();
  return zmin;
}

//______________________________________________________________________________
Float_t EdbTrackP::Zmax() const
{
  Float_t zmax = Z(); 
  if (N() && GetSegmentLast()->Z() > zmax) zmax = GetSegmentLast()->Z();
  return zmax;
}
//______________________________________________________________________________
void EdbTrackP::SetPerrUp(Float_t perrUp)
{
  ePerrUp=perrUp;
}
//______________________________________________________________________________
void EdbTrackP::SetPerrDown(Float_t perrDown)
{
  ePerrDown=perrDown;
}
//______________________________________________________________________________
void EdbTrackP::SetPerr(Float_t perrDown, Float_t perrUp)
{
  ePerrUp=perrUp;
  ePerrDown=perrDown;
}

//______________________________________________________________________________
const EdbSegP  *EdbTrackP::TrackStart() const 
{ 
  if(!N())    return (EdbSegP*)this;
  if(Dir()<0) return GetSegmentLast();
  if(Dir()>0) return GetSegmentFirst();
  return (EdbSegP*)this;
}

//______________________________________________________________________________
const EdbSegP  *EdbTrackP::TrackEnd() const 
{ 
  if(!N())    return (EdbSegP*)this;
  if(Dir()>0) return GetSegmentLast();
  if(Dir()<0) return GetSegmentFirst();
  return (EdbSegP*)this;
}

//==============================================================================
EdbPattern::EdbPattern() : EdbSegmentsBox()
{
  eCell     = new TIndexCell();
  Set0();
}

//______________________________________________________________________________
EdbPattern::EdbPattern(float x0, float y0, float z0, int n) : EdbSegmentsBox(x0,y0,z0,n) 
{
  Set0();
  eCell     = new TIndexCell();
}

//______________________________________________________________________________
EdbPattern::EdbPattern(EdbPattern &p):EdbSegmentsBox(p)
{
  eID = p.eID;
  ePID = p.ePID;
  eSide = p.eSide;
  eScanID = p.eScanID;
  eFlag = p.eFlag;
  eNAff = p.eNAff;
  eStepX = p.eStepX;  eStepY = p.eStepY;
  eStepTX = p.eStepTX;  eStepTY = p.eStepTY;
  for(int i=0; i<4; i++) eSigma[i] = p.eSigma[i];
  if(p.eCell) eCell = new TIndexCell(*(p.eCell));
}
 
//______________________________________________________________________________
EdbPattern::~EdbPattern()
{
  if(eCell)  { delete eCell;  eCell=0; }
}

//______________________________________________________________________________
void EdbPattern::Set0()
{
  eID   = 0;
  ePID  = 0;
  eFlag = 0;
  eNAff = 0;
  eStepX=eStepY=eStepTX=eStepTY=0;
  for(int i=0; i<4; i++) eSigma[i]=0;
  eSide=0;
}

//______________________________________________________________________________
EdbSegP *EdbPattern::FindSegment(int id)
{
  for(int i=0; i<N(); i++) if(GetSegment(i)->ID()==id) return GetSegment(i);
  return 0;
}

//______________________________________________________________________________
float EdbPattern::SummaryPath()
{
  // calculate the microscope path of the eventual prediction scan of this pattern
  float sum=0,dx,dy;
  if(N()<2) return sum;
  for(int i=0; i<N()-1; i++ ) {
    dx = GetSegment(i+1)->X()-GetSegment(i)->X();
    dy = GetSegment(i+1)->Y()-GetSegment(i)->Y();
    sum += Sqrt(dx*dx+dy*dy);
  }
  return sum;
}

//______________________________________________________________________________
void EdbPattern::FillCell( float stepx, float stepy, float steptx, float stepty )
{
  // fill cells with fixed size at z=zPat

  TIndexCell *cell = Cell();
  if(cell) cell->Drop();
  SetStep(stepx,stepy,steptx,stepty);

  float x,y,tx,ty,dz;
  Long_t  val[5];  // x,y,ax,ay,i
  EdbSegP *p;
  int npat = N();
  for(int i=0; i<npat; i++ ) {
    p = GetSegment(i);
    dz = Z() - p->Z();
    tx = p->TX();
    ty = p->TY();
    x  = p->X() + tx*dz;
    y  = p->Y() + ty*dz;
    val[0]= (Long_t)( x / stepx  );
    val[1]= (Long_t)( y / stepy  );
    val[2]= (Long_t)( tx/ steptx );
    val[3]= (Long_t)( ty/ stepty );
    val[4]= (Long_t)(i);
    cell->Add(5,val);
  }
  cell->Sort();

}

//______________________________________________________________________________
int EdbPattern::FindCompliments(EdbSegP &s, TObjArray &arr, float nsigx, float nsigt)
{
  // return the array of segments compatible with the
  // prediction segment s with the accuracy of nsig (in sigmas)

  float dz = Z()-s.Z();
  float x = s.X() + s.TX()*dz;
  float y = s.Y() + s.TY()*dz;
  float sx = TMath::Sqrt( s.SX() + s.STX()*dz*dz );
  float sy = TMath::Sqrt( s.SY() + s.STY()*dz*dz );

  float stx = s.STX();
  if (stx <= 0.) stx = 1.;
  float sty = s.STY();
  if (sty <= 0.) sty = 1.;

  //printf("Step: %f %f %f %f\n",StepX(),StepY(),StepTX(),StepTY());
  //printf("Sigm: %f %f %f %f\n",sx,sy,TMath::Sqrt(stx),TMath::Sqrt(sty));


  long vcent[4] = { (long)(x/StepX()),
		    (long)(y/StepY()),
		    (long)(s.TX()/StepTX()),
		    (long)(s.TY()/StepTY())  };
  long vdiff[4] = { (long)(sx*nsigx/StepX()+1),
		    (long)(sy*nsigx/StepY()+1),
		    (long)(TMath::Sqrt(stx)*nsigt/StepTX()+1),
		    (long)(TMath::Sqrt(sty)*nsigt/StepTY()+1) };

  sy*=sy;
  sx*=sx;

  float dx,dy,dtx,dty;

  int nseg=0;
  EdbSegP *seg=0;


  //printf("find compliments in pattern %d\n",ID());
  //s.Print();

  long vmin[4],vmax[4];
  for(int i=0; i<4; i++) {
    vmin[i] = vcent[i]-vdiff[i];
    vmax[i] = vcent[i]+vdiff[i];
    //printf("%f +- %f \t",(float)vcent[i],(float)vdiff[i]);
  }
  //printf("\n");

  ///TODO: move this cycle into iterator

  TIndexCell *c1=0,*c2=0,*c3=0,*c4=0;
  for(vcent[0]=vmin[0]; vcent[0]<=vmax[0]; vcent[0]++) {
    c1 = eCell->Find(vcent[0]);
    if(!c1) continue;
    for(vcent[1]=vmin[1]; vcent[1]<=vmax[1]; vcent[1]++) {
      c2 = c1->Find(vcent[1]);
      if(!c2) continue;
      for(vcent[2]=vmin[2]; vcent[2]<=vmax[2]; vcent[2]++) {
	c3 = c2->Find(vcent[2]);
	if(!c3) continue;
	for(vcent[3]=vmin[3]; vcent[3]<=vmax[3]; vcent[3]++) {
	  c4 = c3->Find(vcent[3]);
	  if(!c4) continue;

	  for(int i=0; i<c4->N(); i++) {
	    seg = GetSegment(c4->At(i)->Value());
	    dtx=s.TX()-seg->TX();
	    if( dtx*dtx > stx*nsigt*nsigt )    continue;
	    dty=s.TY()-seg->TY();
	    if( dty*dty > sty*nsigt*nsigt )    continue;
	    dz = seg->Z()-s.Z();
	    dx=s.X()+s.TX()*dz-seg->X();
	    if( dx*dx > sx*nsigx*nsigx )           continue;
	    dy=s.Y()+s.TY()*dz-seg->Y();
	    if( dy*dy > sy*nsigx*nsigx )           continue;
	    //if(!s.IsCompatible(*seg,nsigx,nsigt)) continue;
	    arr.Add(seg);
	    nseg++;
	  }

	}
      }
    }
  }


//    TIndexCellIterV itr( eCell, 4, vcent, vdiff );
//    const TIndexCell *c=0;
//    while( (c=itr.Next()) ) 
//      for(int i=0; i<c->N(); i++) {
//        seg = GetSegment(c->At(i)->Value());
//        if(!s.IsCompatible(*seg,nsigx,nsigt)) continue;
//        arr.Add(seg);
//        nseg++;
//      }

  return nseg;
}

//______________________________________________________________________________
void EdbPattern::SetSegmentsPID()
{
  int nseg = N();
  for(int i=0; i<nseg; i++) {
    GetSegment(i)->SetPID(ID());  //PID of the segment must be ID of the pattern!
  }
}

//______________________________________________________________________________
void EdbPattern::SetSegmentsScanID(EdbID id)
{
  for(int i=0; i<N(); i++)     GetSegment(i)->SetScanID(id);
}

//______________________________________________________________________________
EdbPattern *EdbPattern::ExtractSubPattern(float min[5], float max[5], int MCEvt)
{
  //
  // return the copy of selected segments
  //
  EdbPattern *pat = new EdbPattern( X(), Y(), Z() );
  EdbSegP *s;
  int nseg = N();

  for(int i=0; i<nseg; i++) {
    s = GetSegment(i);
		
    if (s->MCEvt()!=MCEvt && s->MCEvt()>0 && MCEvt>=0) continue;
    // Do not continue not in case MCEvt was not specified at all.
    // This allows backward compability.
    
    if(s->X()  < min[0])   continue;
    if(s->Y()  < min[1])   continue;
    if(s->TX() < min[2])   continue;
    if(s->TY() < min[3])   continue;
    if(s->W()  < min[4])   continue;

    if(s->X()  > max[0])   continue;
    if(s->Y()  > max[1])   continue;
    if(s->TX() > max[2])   continue;
    if(s->TY() > max[3])   continue;
    if(s->W()  > max[4])   continue;

    pat->AddSegment(*s);
  }

  return pat;
}

//______________________________________________________________________________
void EdbPattern::Reset()
{
  ((EdbSegmentsBox *)this)->Reset();
  Set0();
  if(eCell)     eCell->Drop();
}

//______________________________________________________________________________
float EdbPattern::Xmean()
{
  Double_t s=0;
  int n=N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += GetSegment(i)->eX;
  s /= n;
  return (float)s;
}

//______________________________________________________________________________
float EdbPattern::Ymean()
{
  Double_t s=0;
  int n=N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += GetSegment(i)->eY;
  s /= n;
  return (float)s;
}

 
////////////////////////////////////////////////////////////////////////////////
//
//   EdbPatternsVolume
//
////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
EdbPatternsVolume::EdbPatternsVolume()
{
  ePatterns   = new TObjArray();
  eTracksCell = 0;
  ePatternsCell = 0;
  Set0();
}

//______________________________________________________________________________
EdbPatternsVolume::EdbPatternsVolume(EdbPatternsVolume &pvol)
{
  ePatterns   = new TObjArray();
  eTracksCell = 0;
  ePatternsCell = 0;
  Set0();

  pvol.PassProperties(*this);
  int npat,nseg;
  npat = Npatterns();
  for(int j=0; j<npat; j++) {
    nseg = GetPattern(j)->N();
    for(int i=0; i<nseg; i++ ) {
      pvol.GetPattern(j)->AddSegment( *(GetPattern(j)->GetSegment(i)) );
    }
  }
}

//______________________________________________________________________________
EdbPatternsVolume::~EdbPatternsVolume()
{
  if(ePatterns) {
    ePatterns->Delete();
    delete ePatterns;
    ePatterns=0;
  }
  if(eTracksCell)   { delete eTracksCell;   eTracksCell=0; }
  if(ePatternsCell) { delete ePatternsCell; ePatternsCell=0; }
}

//______________________________________________________________________________
void EdbPatternsVolume::Set0()
{
  eX=0; eY=0; eZ=0;
  eDescendingZ=0;
}
 
//______________________________________________________________________________
int EdbPatternsVolume::DropCouples()
{
  int count=0;
  int npat=Npatterns();
  for(int i=0; i<npat; i++ )
    count += GetPattern(i)->Cell()->DropCouples(4);
  if(count) printf("%d couples are dropped in volume cells\n",count);
  return count;
}
 
//______________________________________________________________________________
void EdbPatternsVolume::SetPatternsID()
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ )
    GetPattern(i)->SetID(i);
}
 
//______________________________________________________________________________
void EdbPatternsVolume::Transform( const EdbAffine2D *aff )
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ )  {
    GetPattern(i)->Transform(aff);
    GetPattern(i)->TransformARot(aff);
  }
}
 
//______________________________________________________________________________
void EdbPatternsVolume::Centralize()
{
  // find geometrical center (XY) of all patterns and set it as the center of 
  // coordinates  to simplify transformations
  // To be used before any operations on patterns

  float xc=0;
  float yc=0;
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    xc += GetPattern(i)->Xmax() + GetPattern(i)->Xmin();
    yc += GetPattern(i)->Ymax() + GetPattern(i)->Ymin();
  }
  xc = xc/Npatterns()/2;
  yc = yc/Npatterns()/2;

  Centralize(xc,yc);
}

//______________________________________________________________________________
void EdbPatternsVolume::Centralize( float xc, float yc )
{
  eX = xc;  eY=yc;
  Shift(-xc,-yc);
  float npat = Npatterns();
  for(int i=0; i<npat; i++ ) 
    GetPattern(i)->SetKeep(1,0,0,1,0,0);
}

//______________________________________________________________________________
void EdbPatternsVolume::PrintAff() const
{
  EdbAffine2D a;
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    GetPattern(i)->GetKeep(a);
    printf(" %3d  (%5d) Z=%13.3f ",i,GetPattern(i)->NAff(), GetPattern(i)->Z() ); a.Print();
  }
}

//______________________________________________________________________________
void EdbPatternsVolume::PrintStat( Option_t *opt) const
{
  int npat = Npatterns();
  printf("\nVolume statistics for %d patterns\n",npat);

  float dx,dy;
  EdbPattern *pat=0;
  printf("pat# \t segments \t dX \t\tdY \t meanDist \n");
  int i;
  for(i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    dx = pat->Xmax() - pat->Xmin();
    dy = pat->Ymax()- pat->Ymin();
    printf(" %d\t %d\t %10.2f \t %10.2f \t %10.4f \n", 
	   i, pat->GetN(),dx,dy, TMath::Sqrt(dx*dy/pat->GetN()) );
  }

  npat=Npatterns();
  for(i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    pat->Cell()->PrintStat();
  }
}
 
//______________________________________________________________________________
void EdbPatternsVolume::PrintStat(EdbPattern &pat) const
{
} 

//______________________________________________________________________________
void EdbPatternsVolume::Print() const
{
  int npat = Npatterns();
  printf("\nEdbPatternsVolume with %d patterns\n",npat);
  EdbPattern *pat=0;
  for(int i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    printf(" id=%3d pid=%3d   x:y:z =  %12.3f %12.3f %12.3f  n= %8d   ScanID = %s \tside=%d\n", 
	   pat->ID(), pat->PID(), pat->X(),pat->Y(),pat->Z(),pat->N(), pat->ScanID().AsString(), pat->Side() );
  }
  printf("\n");
} 

//______________________________________________________________________________
void EdbPatternsVolume::AddPattern( EdbPattern *pat )
{
  ePatterns->Add(pat);
}

//______________________________________________________________________________
void EdbPatternsVolume::AddPatternAt( EdbPattern *pat,int id )
{
  if(ePatterns->GetSize()<id+1) ePatterns->Expand(id+1);
  pat->SetID(id);
  ePatterns->AddAt(pat,id);
}

//______________________________________________________________________________
EdbPattern *EdbPatternsVolume::GetPattern( int id ) const
{
  if(Npatterns()>id) return (EdbPattern*)ePatterns->At(id);
  else return 0;
}


//______________________________________________________________________________
void EdbPatternsVolume::PassProperties( EdbPatternsVolume &pvol )
{
  pvol.SetXYZ(eX,eY,eZ);
  EdbAffine2D a;
  EdbPattern *p=0;

  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    p = GetPattern(i);
    p->GetKeep(a);
    EdbPattern *psel = new EdbPattern( p->X(),p->Y(),p->Z() );
    psel->SetKeep( a.A11(), a.A12(), a.A21(), a.A22(), a.B1(),a.B2() );
    pvol.AddPattern(psel);
  }
  pvol.SetPatternsID();
}


//______________________________________________________________________________
void EdbPatternsVolume::Shift( float x, float y )
{
  EdbAffine2D aff;
  aff.ShiftX(x);
  aff.ShiftY(y);

  int npat = Npatterns();
  for(int i=0; i<npat; i++)
    GetPattern(i)->Transform(&aff);
}

//______________________________________________________________________________
void EdbPatternsVolume::DropCell()
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) GetPattern(i)->Cell()->Drop();
}

//______________________________________________________________________________
EdbPattern  *EdbPatternsVolume::GetPatternByPlate(int plate, int side)
{
  EdbPattern *p=0;
  for(int i=0; i<Npatterns(); i++) {
    p = GetPattern(i);
    if(p) if( p->Plate()==plate && p->Side()==side ) return p;
  }
  return 0;
}

///______________________________________________________________________________
EdbPattern  *EdbPatternsVolume::InsertPattern(EdbPattern *pat, Bool_t descendingZ)
{
  // insert new pattern using it's Z as a main criteria for the positioning
  EdbPattern *p  = GetPatternByPlate( pat->Plate(), pat->Side() );
  if(p)  {
    Log(1,"EdbPatternsVolume::InsertPattern",
	"ERROR: Attempt to insert new pattern for already existing plate/side z = %f znew = %f", 
	pat->Plate(), pat->Side(), p->Z(), pat->Z());
    return p;
  }
  AddPattern(pat);
  SortPatternsByZ(descendingZ);
  return pat;
}

//______________________________________________________________________________
void EdbPatternsVolume::SortPatternsByZ(Bool_t descendingZ)
{
  int n = Npatterns();
  int   *ind = new int[n];
  float *z   = new float[n];
  for(int i=0; i<Npatterns(); i++) {
    ind[i] = i; 
    z[i] = GetPattern(i)->Z();
  }
  TMath::Sort(n,z,ind,descendingZ);
  TObjArray *pnew = new TObjArray(ePatterns->GetSize());
  for(int i=0; i<n; i++)     pnew->Add(GetPattern(ind[i]));
  ePatterns->Clear();
  delete ePatterns;
  ePatterns = pnew;
  delete ind;
  delete z;
}



//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternZLowestHighest(Bool_t lowestZ) const
{
  // By default it returns the pattern with the lowest Z position. This is
  // NOT necessarily the first pattern of the patterns volume.
  
  if (Npatterns()==0) {
    Log(1,"EdbPatternsVolume::GetPatternZLowestHighest","ERROR: Attempt to get Pattern, but NO patterns are in EdbPatternsVolume!");
    return 0;
  }
  if (gEDBDEBUGLEVEL>3) cout << "EdbPatternsVolume::GetPatternZLowestHighest  lowestZ = " << lowestZ << endl;
  
  EdbPattern *pat=0;
  pat=GetPattern(0);
  if (Npatterns()==1) return pat;
  float zpat0,zpat1;
  zpat0=GetPattern(0)->Z();
  zpat1=GetPattern(Npatterns()-1)->Z();
  if (gEDBDEBUGLEVEL>2) {
  cout << "EdbPatternsVolume  Pattern At Position 0    has Z = " << zpat0 << endl;
  cout << "EdbPatternsVolume  Pattern At Position last has Z = " << zpat1 << endl;
  }
  if (zpat0>zpat1) pat=GetPattern(Npatterns()-1);
  if (kTRUE==lowestZ) return pat;
  if (kFALSE==lowestZ && zpat0<zpat1) return GetPattern(Npatterns()-1);
  return GetPattern(0);
}


//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternSucceding(EdbPattern* pat) const
{
	return GetPatternNext(pat->Z(),1);
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternPreceding(EdbPattern* pat) const
{
	return GetPatternNext(pat->Z(),-1);
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternNext(float z, int dir) const
{
  // Return next pattern in either downstream (dir=+1) (i.e. Z gets bigger)
  //                  or in either   upstream (dir=-1) (i.e. Z gets smaller)
  // direction.
  // If the last/first pattern is already given and one ask therefor 
  // for a pattern that does not exists, then return a NULL pattern.
  
  EdbPattern *pat=0;
  float dz=dir*99999999.;
  float zpat;
  for(int i=0; i<Npatterns(); i++) {
    zpat=GetPattern(i)->Z();
    if(dir>0) if(zpat>z) if(zpat-z<dz) {dz=zpat-z; pat=GetPattern(i);}
    if(dir<0) if(zpat<z) if(zpat-z>dz) {dz=zpat-z; pat=GetPattern(i);}
  }
  if (gEDBDEBUGLEVEL>2) 
		if (pat == 0) cout <<  "EdbPattern* EdbPatternsVolume::GetPatternNext()   WARNING: pattern is NULL ! Check your input!" << endl;
  return pat;
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternByPID(int pid) const
{
  // Return pattern having PID()  == pid.
  // If there is no such pattern, return a NULL pattern.
  //
  // This is important because some codings (libShower for example).
  // rely on getting Patterns by the PIDs.
  
  EdbPattern *pat=0;
	for(int i=0; i<Npatterns(); i++) {
		pat = GetPattern(i);
    if (GetPattern(i)->PID() == pid ) return pat;
	}
	
  if (gEDBDEBUGLEVEL>2) 
		if (pat == 0) cout <<  "EdbPattern* EdbPatternsVolume::GetPatternNext()   WARNING: pattern is NULL ! Check your input!" << endl;
  return NULL;
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternByZ(float z) const
{
  // Return pattern having Z()  == z +- 5(microns) due to roundings.
  // If there is no such pattern, return a NULL pattern.
  //
  // This is important because some codings (libShower for example).
  // rely on getting Patterns by the Z.
  
  EdbPattern *pat=0;
	for(int i=0; i<Npatterns(); i++) {
		pat = GetPattern(i);
    if ( TMath::Abs(GetPattern(i)->Z() - z)<5 ) return pat;
	}
	
  if (gEDBDEBUGLEVEL>2) 
		if (pat == 0) cout <<  "EdbPattern* EdbPatternsVolume::GetPatternNext()   WARNING: pattern is NULL ! Check your input!" << endl;
  return NULL;
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::NextPattern(float z, int dir) const
{
  // Return next pattern in either downstream (dir=+1)
  //                  or in either   upstream (dir=-1)
  // direction.
  
  EdbPattern *pat=0;
  float dz=dir*99999999.;
  float zpat;
  for(int i=0; i<Npatterns(); i++) {
    zpat=GetPattern(i)->Z();
    if(dir>0) if(zpat>z) if(zpat-z<dz) {dz=zpat-z; pat=GetPattern(i);}
    if(dir<0) if(zpat<z) if(zpat-z>dz) {dz=zpat-z; pat=GetPattern(i);}
  }
  return pat;
}

//______________________________________________________________________________
int EdbPatternsVolume::FindComplimentsVol(EdbSegP &ss, TObjArray &arr, float nsig, float nsigt, int Dpat)
{
  EdbPattern *pat  = 0;
  int npat = Npatterns();
  bool p_inverse_z = (GetPattern(1)->Z() - GetPattern(0)->Z()) > 0. ? false : true ;
  int ii = 0;
  int dpat = 1;
  int pid = 0;
  if (p_inverse_z) dpat = -1;
  for (int i = 0; i<npat; i++)
  {
    ii = i;
    if (p_inverse_z) ii = npat-1-i;
    pat = GetPattern(ii);  
    if (ss.Z() < pat->Z())
    {
	pid = ii - dpat;
	if ( pid >= npat )  ss.SetPID(npat-1);
	if ( pid <  0 )     ss.SetPID(0);
	else	            ss.SetPID(pid);
	break;
    }
    else if (ss.Z() == pat->Z())
    {
	ss.SetPID(ii);
	break;
    }
  }

  int p0 = ss.PID();
  int pend   = p0 + Dpat;
  if (pend >= npat) pend = npat - 1;
  int pstart   = p0 - Dpat;
  if (pstart < 0) pstart = 0;

  //arr.Clear();
  int nseg = 0;
  int n = 0;
  for(int i=pstart; i <=pend; i++ ) {
    pat = GetPattern(i);
    if(!pat)                   continue;
    n = pat->FindCompliments(ss,arr,nsig,nsigt);
    nseg += n;
  }

  return nseg;
}

//______________________________________________________________________________
Float_t EdbPatternsVolume::Xmean()
{
	if (Npatterns()==0) return 0;
	EdbPattern* pat = GetPattern(0);
	Double_t s=0;
  int n=pat->N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += pat->GetSegment(i)->eX;
  s /= n;
  return (float)s;
}

//______________________________________________________________________________
Float_t EdbPatternsVolume::Ymean()
{
	if (Npatterns()==0) return 0;
	EdbPattern* pat = GetPattern(0);
	Double_t s=0;
  int n=pat->N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += pat->GetSegment(i)->eY;
  s /= n;
  return (float)s;
}
 EdbPattern.cxx:1
 EdbPattern.cxx:2
 EdbPattern.cxx:3
 EdbPattern.cxx:4
 EdbPattern.cxx:5
 EdbPattern.cxx:6
 EdbPattern.cxx:7
 EdbPattern.cxx:8
 EdbPattern.cxx:9
 EdbPattern.cxx:10
 EdbPattern.cxx:11
 EdbPattern.cxx:12
 EdbPattern.cxx:13
 EdbPattern.cxx:14
 EdbPattern.cxx:15
 EdbPattern.cxx:16
 EdbPattern.cxx:17
 EdbPattern.cxx:18
 EdbPattern.cxx:19
 EdbPattern.cxx:20
 EdbPattern.cxx:21
 EdbPattern.cxx:22
 EdbPattern.cxx:23
 EdbPattern.cxx:24
 EdbPattern.cxx:25
 EdbPattern.cxx:26
 EdbPattern.cxx:27
 EdbPattern.cxx:28
 EdbPattern.cxx:29
 EdbPattern.cxx:30
 EdbPattern.cxx:31
 EdbPattern.cxx:32
 EdbPattern.cxx:33
 EdbPattern.cxx:34
 EdbPattern.cxx:35
 EdbPattern.cxx:36
 EdbPattern.cxx:37
 EdbPattern.cxx:38
 EdbPattern.cxx:39
 EdbPattern.cxx:40
 EdbPattern.cxx:41
 EdbPattern.cxx:42
 EdbPattern.cxx:43
 EdbPattern.cxx:44
 EdbPattern.cxx:45
 EdbPattern.cxx:46
 EdbPattern.cxx:47
 EdbPattern.cxx:48
 EdbPattern.cxx:49
 EdbPattern.cxx:50
 EdbPattern.cxx:51
 EdbPattern.cxx:52
 EdbPattern.cxx:53
 EdbPattern.cxx:54
 EdbPattern.cxx:55
 EdbPattern.cxx:56
 EdbPattern.cxx:57
 EdbPattern.cxx:58
 EdbPattern.cxx:59
 EdbPattern.cxx:60
 EdbPattern.cxx:61
 EdbPattern.cxx:62
 EdbPattern.cxx:63
 EdbPattern.cxx:64
 EdbPattern.cxx:65
 EdbPattern.cxx:66
 EdbPattern.cxx:67
 EdbPattern.cxx:68
 EdbPattern.cxx:69
 EdbPattern.cxx:70
 EdbPattern.cxx:71
 EdbPattern.cxx:72
 EdbPattern.cxx:73
 EdbPattern.cxx:74
 EdbPattern.cxx:75
 EdbPattern.cxx:76
 EdbPattern.cxx:77
 EdbPattern.cxx:78
 EdbPattern.cxx:79
 EdbPattern.cxx:80
 EdbPattern.cxx:81
 EdbPattern.cxx:82
 EdbPattern.cxx:83
 EdbPattern.cxx:84
 EdbPattern.cxx:85
 EdbPattern.cxx:86
 EdbPattern.cxx:87
 EdbPattern.cxx:88
 EdbPattern.cxx:89
 EdbPattern.cxx:90
 EdbPattern.cxx:91
 EdbPattern.cxx:92
 EdbPattern.cxx:93
 EdbPattern.cxx:94
 EdbPattern.cxx:95
 EdbPattern.cxx:96
 EdbPattern.cxx:97
 EdbPattern.cxx:98
 EdbPattern.cxx:99
 EdbPattern.cxx:100
 EdbPattern.cxx:101
 EdbPattern.cxx:102
 EdbPattern.cxx:103
 EdbPattern.cxx:104
 EdbPattern.cxx:105
 EdbPattern.cxx:106
 EdbPattern.cxx:107
 EdbPattern.cxx:108
 EdbPattern.cxx:109
 EdbPattern.cxx:110
 EdbPattern.cxx:111
 EdbPattern.cxx:112
 EdbPattern.cxx:113
 EdbPattern.cxx:114
 EdbPattern.cxx:115
 EdbPattern.cxx:116
 EdbPattern.cxx:117
 EdbPattern.cxx:118
 EdbPattern.cxx:119
 EdbPattern.cxx:120
 EdbPattern.cxx:121
 EdbPattern.cxx:122
 EdbPattern.cxx:123
 EdbPattern.cxx:124
 EdbPattern.cxx:125
 EdbPattern.cxx:126
 EdbPattern.cxx:127
 EdbPattern.cxx:128
 EdbPattern.cxx:129
 EdbPattern.cxx:130
 EdbPattern.cxx:131
 EdbPattern.cxx:132
 EdbPattern.cxx:133
 EdbPattern.cxx:134
 EdbPattern.cxx:135
 EdbPattern.cxx:136
 EdbPattern.cxx:137
 EdbPattern.cxx:138
 EdbPattern.cxx:139
 EdbPattern.cxx:140
 EdbPattern.cxx:141
 EdbPattern.cxx:142
 EdbPattern.cxx:143
 EdbPattern.cxx:144
 EdbPattern.cxx:145
 EdbPattern.cxx:146
 EdbPattern.cxx:147
 EdbPattern.cxx:148
 EdbPattern.cxx:149
 EdbPattern.cxx:150
 EdbPattern.cxx:151
 EdbPattern.cxx:152
 EdbPattern.cxx:153
 EdbPattern.cxx:154
 EdbPattern.cxx:155
 EdbPattern.cxx:156
 EdbPattern.cxx:157
 EdbPattern.cxx:158
 EdbPattern.cxx:159
 EdbPattern.cxx:160
 EdbPattern.cxx:161
 EdbPattern.cxx:162
 EdbPattern.cxx:163
 EdbPattern.cxx:164
 EdbPattern.cxx:165
 EdbPattern.cxx:166
 EdbPattern.cxx:167
 EdbPattern.cxx:168
 EdbPattern.cxx:169
 EdbPattern.cxx:170
 EdbPattern.cxx:171
 EdbPattern.cxx:172
 EdbPattern.cxx:173
 EdbPattern.cxx:174
 EdbPattern.cxx:175
 EdbPattern.cxx:176
 EdbPattern.cxx:177
 EdbPattern.cxx:178
 EdbPattern.cxx:179
 EdbPattern.cxx:180
 EdbPattern.cxx:181
 EdbPattern.cxx:182
 EdbPattern.cxx:183
 EdbPattern.cxx:184
 EdbPattern.cxx:185
 EdbPattern.cxx:186
 EdbPattern.cxx:187
 EdbPattern.cxx:188
 EdbPattern.cxx:189
 EdbPattern.cxx:190
 EdbPattern.cxx:191
 EdbPattern.cxx:192
 EdbPattern.cxx:193
 EdbPattern.cxx:194
 EdbPattern.cxx:195
 EdbPattern.cxx:196
 EdbPattern.cxx:197
 EdbPattern.cxx:198
 EdbPattern.cxx:199
 EdbPattern.cxx:200
 EdbPattern.cxx:201
 EdbPattern.cxx:202
 EdbPattern.cxx:203
 EdbPattern.cxx:204
 EdbPattern.cxx:205
 EdbPattern.cxx:206
 EdbPattern.cxx:207
 EdbPattern.cxx:208
 EdbPattern.cxx:209
 EdbPattern.cxx:210
 EdbPattern.cxx:211
 EdbPattern.cxx:212
 EdbPattern.cxx:213
 EdbPattern.cxx:214
 EdbPattern.cxx:215
 EdbPattern.cxx:216
 EdbPattern.cxx:217
 EdbPattern.cxx:218
 EdbPattern.cxx:219
 EdbPattern.cxx:220
 EdbPattern.cxx:221
 EdbPattern.cxx:222
 EdbPattern.cxx:223
 EdbPattern.cxx:224
 EdbPattern.cxx:225
 EdbPattern.cxx:226
 EdbPattern.cxx:227
 EdbPattern.cxx:228
 EdbPattern.cxx:229
 EdbPattern.cxx:230
 EdbPattern.cxx:231
 EdbPattern.cxx:232
 EdbPattern.cxx:233
 EdbPattern.cxx:234
 EdbPattern.cxx:235
 EdbPattern.cxx:236
 EdbPattern.cxx:237
 EdbPattern.cxx:238
 EdbPattern.cxx:239
 EdbPattern.cxx:240
 EdbPattern.cxx:241
 EdbPattern.cxx:242
 EdbPattern.cxx:243
 EdbPattern.cxx:244
 EdbPattern.cxx:245
 EdbPattern.cxx:246
 EdbPattern.cxx:247
 EdbPattern.cxx:248
 EdbPattern.cxx:249
 EdbPattern.cxx:250
 EdbPattern.cxx:251
 EdbPattern.cxx:252
 EdbPattern.cxx:253
 EdbPattern.cxx:254
 EdbPattern.cxx:255
 EdbPattern.cxx:256
 EdbPattern.cxx:257
 EdbPattern.cxx:258
 EdbPattern.cxx:259
 EdbPattern.cxx:260
 EdbPattern.cxx:261
 EdbPattern.cxx:262
 EdbPattern.cxx:263
 EdbPattern.cxx:264
 EdbPattern.cxx:265
 EdbPattern.cxx:266
 EdbPattern.cxx:267
 EdbPattern.cxx:268
 EdbPattern.cxx:269
 EdbPattern.cxx:270
 EdbPattern.cxx:271
 EdbPattern.cxx:272
 EdbPattern.cxx:273
 EdbPattern.cxx:274
 EdbPattern.cxx:275
 EdbPattern.cxx:276
 EdbPattern.cxx:277
 EdbPattern.cxx:278
 EdbPattern.cxx:279
 EdbPattern.cxx:280
 EdbPattern.cxx:281
 EdbPattern.cxx:282
 EdbPattern.cxx:283
 EdbPattern.cxx:284
 EdbPattern.cxx:285
 EdbPattern.cxx:286
 EdbPattern.cxx:287
 EdbPattern.cxx:288
 EdbPattern.cxx:289
 EdbPattern.cxx:290
 EdbPattern.cxx:291
 EdbPattern.cxx:292
 EdbPattern.cxx:293
 EdbPattern.cxx:294
 EdbPattern.cxx:295
 EdbPattern.cxx:296
 EdbPattern.cxx:297
 EdbPattern.cxx:298
 EdbPattern.cxx:299
 EdbPattern.cxx:300
 EdbPattern.cxx:301
 EdbPattern.cxx:302
 EdbPattern.cxx:303
 EdbPattern.cxx:304
 EdbPattern.cxx:305
 EdbPattern.cxx:306
 EdbPattern.cxx:307
 EdbPattern.cxx:308
 EdbPattern.cxx:309
 EdbPattern.cxx:310
 EdbPattern.cxx:311
 EdbPattern.cxx:312
 EdbPattern.cxx:313
 EdbPattern.cxx:314
 EdbPattern.cxx:315
 EdbPattern.cxx:316
 EdbPattern.cxx:317
 EdbPattern.cxx:318
 EdbPattern.cxx:319
 EdbPattern.cxx:320
 EdbPattern.cxx:321
 EdbPattern.cxx:322
 EdbPattern.cxx:323
 EdbPattern.cxx:324
 EdbPattern.cxx:325
 EdbPattern.cxx:326
 EdbPattern.cxx:327
 EdbPattern.cxx:328
 EdbPattern.cxx:329
 EdbPattern.cxx:330
 EdbPattern.cxx:331
 EdbPattern.cxx:332
 EdbPattern.cxx:333
 EdbPattern.cxx:334
 EdbPattern.cxx:335
 EdbPattern.cxx:336
 EdbPattern.cxx:337
 EdbPattern.cxx:338
 EdbPattern.cxx:339
 EdbPattern.cxx:340
 EdbPattern.cxx:341
 EdbPattern.cxx:342
 EdbPattern.cxx:343
 EdbPattern.cxx:344
 EdbPattern.cxx:345
 EdbPattern.cxx:346
 EdbPattern.cxx:347
 EdbPattern.cxx:348
 EdbPattern.cxx:349
 EdbPattern.cxx:350
 EdbPattern.cxx:351
 EdbPattern.cxx:352
 EdbPattern.cxx:353
 EdbPattern.cxx:354
 EdbPattern.cxx:355
 EdbPattern.cxx:356
 EdbPattern.cxx:357
 EdbPattern.cxx:358
 EdbPattern.cxx:359
 EdbPattern.cxx:360
 EdbPattern.cxx:361
 EdbPattern.cxx:362
 EdbPattern.cxx:363
 EdbPattern.cxx:364
 EdbPattern.cxx:365
 EdbPattern.cxx:366
 EdbPattern.cxx:367
 EdbPattern.cxx:368
 EdbPattern.cxx:369
 EdbPattern.cxx:370
 EdbPattern.cxx:371
 EdbPattern.cxx:372
 EdbPattern.cxx:373
 EdbPattern.cxx:374
 EdbPattern.cxx:375
 EdbPattern.cxx:376
 EdbPattern.cxx:377
 EdbPattern.cxx:378
 EdbPattern.cxx:379
 EdbPattern.cxx:380
 EdbPattern.cxx:381
 EdbPattern.cxx:382
 EdbPattern.cxx:383
 EdbPattern.cxx:384
 EdbPattern.cxx:385
 EdbPattern.cxx:386
 EdbPattern.cxx:387
 EdbPattern.cxx:388
 EdbPattern.cxx:389
 EdbPattern.cxx:390
 EdbPattern.cxx:391
 EdbPattern.cxx:392
 EdbPattern.cxx:393
 EdbPattern.cxx:394
 EdbPattern.cxx:395
 EdbPattern.cxx:396
 EdbPattern.cxx:397
 EdbPattern.cxx:398
 EdbPattern.cxx:399
 EdbPattern.cxx:400
 EdbPattern.cxx:401
 EdbPattern.cxx:402
 EdbPattern.cxx:403
 EdbPattern.cxx:404
 EdbPattern.cxx:405
 EdbPattern.cxx:406
 EdbPattern.cxx:407
 EdbPattern.cxx:408
 EdbPattern.cxx:409
 EdbPattern.cxx:410
 EdbPattern.cxx:411
 EdbPattern.cxx:412
 EdbPattern.cxx:413
 EdbPattern.cxx:414
 EdbPattern.cxx:415
 EdbPattern.cxx:416
 EdbPattern.cxx:417
 EdbPattern.cxx:418
 EdbPattern.cxx:419
 EdbPattern.cxx:420
 EdbPattern.cxx:421
 EdbPattern.cxx:422
 EdbPattern.cxx:423
 EdbPattern.cxx:424
 EdbPattern.cxx:425
 EdbPattern.cxx:426
 EdbPattern.cxx:427
 EdbPattern.cxx:428
 EdbPattern.cxx:429
 EdbPattern.cxx:430
 EdbPattern.cxx:431
 EdbPattern.cxx:432
 EdbPattern.cxx:433
 EdbPattern.cxx:434
 EdbPattern.cxx:435
 EdbPattern.cxx:436
 EdbPattern.cxx:437
 EdbPattern.cxx:438
 EdbPattern.cxx:439
 EdbPattern.cxx:440
 EdbPattern.cxx:441
 EdbPattern.cxx:442
 EdbPattern.cxx:443
 EdbPattern.cxx:444
 EdbPattern.cxx:445
 EdbPattern.cxx:446
 EdbPattern.cxx:447
 EdbPattern.cxx:448
 EdbPattern.cxx:449
 EdbPattern.cxx:450
 EdbPattern.cxx:451
 EdbPattern.cxx:452
 EdbPattern.cxx:453
 EdbPattern.cxx:454
 EdbPattern.cxx:455
 EdbPattern.cxx:456
 EdbPattern.cxx:457
 EdbPattern.cxx:458
 EdbPattern.cxx:459
 EdbPattern.cxx:460
 EdbPattern.cxx:461
 EdbPattern.cxx:462
 EdbPattern.cxx:463
 EdbPattern.cxx:464
 EdbPattern.cxx:465
 EdbPattern.cxx:466
 EdbPattern.cxx:467
 EdbPattern.cxx:468
 EdbPattern.cxx:469
 EdbPattern.cxx:470
 EdbPattern.cxx:471
 EdbPattern.cxx:472
 EdbPattern.cxx:473
 EdbPattern.cxx:474
 EdbPattern.cxx:475
 EdbPattern.cxx:476
 EdbPattern.cxx:477
 EdbPattern.cxx:478
 EdbPattern.cxx:479
 EdbPattern.cxx:480
 EdbPattern.cxx:481
 EdbPattern.cxx:482
 EdbPattern.cxx:483
 EdbPattern.cxx:484
 EdbPattern.cxx:485
 EdbPattern.cxx:486
 EdbPattern.cxx:487
 EdbPattern.cxx:488
 EdbPattern.cxx:489
 EdbPattern.cxx:490
 EdbPattern.cxx:491
 EdbPattern.cxx:492
 EdbPattern.cxx:493
 EdbPattern.cxx:494
 EdbPattern.cxx:495
 EdbPattern.cxx:496
 EdbPattern.cxx:497
 EdbPattern.cxx:498
 EdbPattern.cxx:499
 EdbPattern.cxx:500
 EdbPattern.cxx:501
 EdbPattern.cxx:502
 EdbPattern.cxx:503
 EdbPattern.cxx:504
 EdbPattern.cxx:505
 EdbPattern.cxx:506
 EdbPattern.cxx:507
 EdbPattern.cxx:508
 EdbPattern.cxx:509
 EdbPattern.cxx:510
 EdbPattern.cxx:511
 EdbPattern.cxx:512
 EdbPattern.cxx:513
 EdbPattern.cxx:514
 EdbPattern.cxx:515
 EdbPattern.cxx:516
 EdbPattern.cxx:517
 EdbPattern.cxx:518
 EdbPattern.cxx:519
 EdbPattern.cxx:520
 EdbPattern.cxx:521
 EdbPattern.cxx:522
 EdbPattern.cxx:523
 EdbPattern.cxx:524
 EdbPattern.cxx:525
 EdbPattern.cxx:526
 EdbPattern.cxx:527
 EdbPattern.cxx:528
 EdbPattern.cxx:529
 EdbPattern.cxx:530
 EdbPattern.cxx:531
 EdbPattern.cxx:532
 EdbPattern.cxx:533
 EdbPattern.cxx:534
 EdbPattern.cxx:535
 EdbPattern.cxx:536
 EdbPattern.cxx:537
 EdbPattern.cxx:538
 EdbPattern.cxx:539
 EdbPattern.cxx:540
 EdbPattern.cxx:541
 EdbPattern.cxx:542
 EdbPattern.cxx:543
 EdbPattern.cxx:544
 EdbPattern.cxx:545
 EdbPattern.cxx:546
 EdbPattern.cxx:547
 EdbPattern.cxx:548
 EdbPattern.cxx:549
 EdbPattern.cxx:550
 EdbPattern.cxx:551
 EdbPattern.cxx:552
 EdbPattern.cxx:553
 EdbPattern.cxx:554
 EdbPattern.cxx:555
 EdbPattern.cxx:556
 EdbPattern.cxx:557
 EdbPattern.cxx:558
 EdbPattern.cxx:559
 EdbPattern.cxx:560
 EdbPattern.cxx:561
 EdbPattern.cxx:562
 EdbPattern.cxx:563
 EdbPattern.cxx:564
 EdbPattern.cxx:565
 EdbPattern.cxx:566
 EdbPattern.cxx:567
 EdbPattern.cxx:568
 EdbPattern.cxx:569
 EdbPattern.cxx:570
 EdbPattern.cxx:571
 EdbPattern.cxx:572
 EdbPattern.cxx:573
 EdbPattern.cxx:574
 EdbPattern.cxx:575
 EdbPattern.cxx:576
 EdbPattern.cxx:577
 EdbPattern.cxx:578
 EdbPattern.cxx:579
 EdbPattern.cxx:580
 EdbPattern.cxx:581
 EdbPattern.cxx:582
 EdbPattern.cxx:583
 EdbPattern.cxx:584
 EdbPattern.cxx:585
 EdbPattern.cxx:586
 EdbPattern.cxx:587
 EdbPattern.cxx:588
 EdbPattern.cxx:589
 EdbPattern.cxx:590
 EdbPattern.cxx:591
 EdbPattern.cxx:592
 EdbPattern.cxx:593
 EdbPattern.cxx:594
 EdbPattern.cxx:595
 EdbPattern.cxx:596
 EdbPattern.cxx:597
 EdbPattern.cxx:598
 EdbPattern.cxx:599
 EdbPattern.cxx:600
 EdbPattern.cxx:601
 EdbPattern.cxx:602
 EdbPattern.cxx:603
 EdbPattern.cxx:604
 EdbPattern.cxx:605
 EdbPattern.cxx:606
 EdbPattern.cxx:607
 EdbPattern.cxx:608
 EdbPattern.cxx:609
 EdbPattern.cxx:610
 EdbPattern.cxx:611
 EdbPattern.cxx:612
 EdbPattern.cxx:613
 EdbPattern.cxx:614
 EdbPattern.cxx:615
 EdbPattern.cxx:616
 EdbPattern.cxx:617
 EdbPattern.cxx:618
 EdbPattern.cxx:619
 EdbPattern.cxx:620
 EdbPattern.cxx:621
 EdbPattern.cxx:622
 EdbPattern.cxx:623
 EdbPattern.cxx:624
 EdbPattern.cxx:625
 EdbPattern.cxx:626
 EdbPattern.cxx:627
 EdbPattern.cxx:628
 EdbPattern.cxx:629
 EdbPattern.cxx:630
 EdbPattern.cxx:631
 EdbPattern.cxx:632
 EdbPattern.cxx:633
 EdbPattern.cxx:634
 EdbPattern.cxx:635
 EdbPattern.cxx:636
 EdbPattern.cxx:637
 EdbPattern.cxx:638
 EdbPattern.cxx:639
 EdbPattern.cxx:640
 EdbPattern.cxx:641
 EdbPattern.cxx:642
 EdbPattern.cxx:643
 EdbPattern.cxx:644
 EdbPattern.cxx:645
 EdbPattern.cxx:646
 EdbPattern.cxx:647
 EdbPattern.cxx:648
 EdbPattern.cxx:649
 EdbPattern.cxx:650
 EdbPattern.cxx:651
 EdbPattern.cxx:652
 EdbPattern.cxx:653
 EdbPattern.cxx:654
 EdbPattern.cxx:655
 EdbPattern.cxx:656
 EdbPattern.cxx:657
 EdbPattern.cxx:658
 EdbPattern.cxx:659
 EdbPattern.cxx:660
 EdbPattern.cxx:661
 EdbPattern.cxx:662
 EdbPattern.cxx:663
 EdbPattern.cxx:664
 EdbPattern.cxx:665
 EdbPattern.cxx:666
 EdbPattern.cxx:667
 EdbPattern.cxx:668
 EdbPattern.cxx:669
 EdbPattern.cxx:670
 EdbPattern.cxx:671
 EdbPattern.cxx:672
 EdbPattern.cxx:673
 EdbPattern.cxx:674
 EdbPattern.cxx:675
 EdbPattern.cxx:676
 EdbPattern.cxx:677
 EdbPattern.cxx:678
 EdbPattern.cxx:679
 EdbPattern.cxx:680
 EdbPattern.cxx:681
 EdbPattern.cxx:682
 EdbPattern.cxx:683
 EdbPattern.cxx:684
 EdbPattern.cxx:685
 EdbPattern.cxx:686
 EdbPattern.cxx:687
 EdbPattern.cxx:688
 EdbPattern.cxx:689
 EdbPattern.cxx:690
 EdbPattern.cxx:691
 EdbPattern.cxx:692
 EdbPattern.cxx:693
 EdbPattern.cxx:694
 EdbPattern.cxx:695
 EdbPattern.cxx:696
 EdbPattern.cxx:697
 EdbPattern.cxx:698
 EdbPattern.cxx:699
 EdbPattern.cxx:700
 EdbPattern.cxx:701
 EdbPattern.cxx:702
 EdbPattern.cxx:703
 EdbPattern.cxx:704
 EdbPattern.cxx:705
 EdbPattern.cxx:706
 EdbPattern.cxx:707
 EdbPattern.cxx:708
 EdbPattern.cxx:709
 EdbPattern.cxx:710
 EdbPattern.cxx:711
 EdbPattern.cxx:712
 EdbPattern.cxx:713
 EdbPattern.cxx:714
 EdbPattern.cxx:715
 EdbPattern.cxx:716
 EdbPattern.cxx:717
 EdbPattern.cxx:718
 EdbPattern.cxx:719
 EdbPattern.cxx:720
 EdbPattern.cxx:721
 EdbPattern.cxx:722
 EdbPattern.cxx:723
 EdbPattern.cxx:724
 EdbPattern.cxx:725
 EdbPattern.cxx:726
 EdbPattern.cxx:727
 EdbPattern.cxx:728
 EdbPattern.cxx:729
 EdbPattern.cxx:730
 EdbPattern.cxx:731
 EdbPattern.cxx:732
 EdbPattern.cxx:733
 EdbPattern.cxx:734
 EdbPattern.cxx:735
 EdbPattern.cxx:736
 EdbPattern.cxx:737
 EdbPattern.cxx:738
 EdbPattern.cxx:739
 EdbPattern.cxx:740
 EdbPattern.cxx:741
 EdbPattern.cxx:742
 EdbPattern.cxx:743
 EdbPattern.cxx:744
 EdbPattern.cxx:745
 EdbPattern.cxx:746
 EdbPattern.cxx:747
 EdbPattern.cxx:748
 EdbPattern.cxx:749
 EdbPattern.cxx:750
 EdbPattern.cxx:751
 EdbPattern.cxx:752
 EdbPattern.cxx:753
 EdbPattern.cxx:754
 EdbPattern.cxx:755
 EdbPattern.cxx:756
 EdbPattern.cxx:757
 EdbPattern.cxx:758
 EdbPattern.cxx:759
 EdbPattern.cxx:760
 EdbPattern.cxx:761
 EdbPattern.cxx:762
 EdbPattern.cxx:763
 EdbPattern.cxx:764
 EdbPattern.cxx:765
 EdbPattern.cxx:766
 EdbPattern.cxx:767
 EdbPattern.cxx:768
 EdbPattern.cxx:769
 EdbPattern.cxx:770
 EdbPattern.cxx:771
 EdbPattern.cxx:772
 EdbPattern.cxx:773
 EdbPattern.cxx:774
 EdbPattern.cxx:775
 EdbPattern.cxx:776
 EdbPattern.cxx:777
 EdbPattern.cxx:778
 EdbPattern.cxx:779
 EdbPattern.cxx:780
 EdbPattern.cxx:781
 EdbPattern.cxx:782
 EdbPattern.cxx:783
 EdbPattern.cxx:784
 EdbPattern.cxx:785
 EdbPattern.cxx:786
 EdbPattern.cxx:787
 EdbPattern.cxx:788
 EdbPattern.cxx:789
 EdbPattern.cxx:790
 EdbPattern.cxx:791
 EdbPattern.cxx:792
 EdbPattern.cxx:793
 EdbPattern.cxx:794
 EdbPattern.cxx:795
 EdbPattern.cxx:796
 EdbPattern.cxx:797
 EdbPattern.cxx:798
 EdbPattern.cxx:799
 EdbPattern.cxx:800
 EdbPattern.cxx:801
 EdbPattern.cxx:802
 EdbPattern.cxx:803
 EdbPattern.cxx:804
 EdbPattern.cxx:805
 EdbPattern.cxx:806
 EdbPattern.cxx:807
 EdbPattern.cxx:808
 EdbPattern.cxx:809
 EdbPattern.cxx:810
 EdbPattern.cxx:811
 EdbPattern.cxx:812
 EdbPattern.cxx:813
 EdbPattern.cxx:814
 EdbPattern.cxx:815
 EdbPattern.cxx:816
 EdbPattern.cxx:817
 EdbPattern.cxx:818
 EdbPattern.cxx:819
 EdbPattern.cxx:820
 EdbPattern.cxx:821
 EdbPattern.cxx:822
 EdbPattern.cxx:823
 EdbPattern.cxx:824
 EdbPattern.cxx:825
 EdbPattern.cxx:826
 EdbPattern.cxx:827
 EdbPattern.cxx:828
 EdbPattern.cxx:829
 EdbPattern.cxx:830
 EdbPattern.cxx:831
 EdbPattern.cxx:832
 EdbPattern.cxx:833
 EdbPattern.cxx:834
 EdbPattern.cxx:835
 EdbPattern.cxx:836
 EdbPattern.cxx:837
 EdbPattern.cxx:838
 EdbPattern.cxx:839
 EdbPattern.cxx:840
 EdbPattern.cxx:841
 EdbPattern.cxx:842
 EdbPattern.cxx:843
 EdbPattern.cxx:844
 EdbPattern.cxx:845
 EdbPattern.cxx:846
 EdbPattern.cxx:847
 EdbPattern.cxx:848
 EdbPattern.cxx:849
 EdbPattern.cxx:850
 EdbPattern.cxx:851
 EdbPattern.cxx:852
 EdbPattern.cxx:853
 EdbPattern.cxx:854
 EdbPattern.cxx:855
 EdbPattern.cxx:856
 EdbPattern.cxx:857
 EdbPattern.cxx:858
 EdbPattern.cxx:859
 EdbPattern.cxx:860
 EdbPattern.cxx:861
 EdbPattern.cxx:862
 EdbPattern.cxx:863
 EdbPattern.cxx:864
 EdbPattern.cxx:865
 EdbPattern.cxx:866
 EdbPattern.cxx:867
 EdbPattern.cxx:868
 EdbPattern.cxx:869
 EdbPattern.cxx:870
 EdbPattern.cxx:871
 EdbPattern.cxx:872
 EdbPattern.cxx:873
 EdbPattern.cxx:874
 EdbPattern.cxx:875
 EdbPattern.cxx:876
 EdbPattern.cxx:877
 EdbPattern.cxx:878
 EdbPattern.cxx:879
 EdbPattern.cxx:880
 EdbPattern.cxx:881
 EdbPattern.cxx:882
 EdbPattern.cxx:883
 EdbPattern.cxx:884
 EdbPattern.cxx:885
 EdbPattern.cxx:886
 EdbPattern.cxx:887
 EdbPattern.cxx:888
 EdbPattern.cxx:889
 EdbPattern.cxx:890
 EdbPattern.cxx:891
 EdbPattern.cxx:892
 EdbPattern.cxx:893
 EdbPattern.cxx:894
 EdbPattern.cxx:895
 EdbPattern.cxx:896
 EdbPattern.cxx:897
 EdbPattern.cxx:898
 EdbPattern.cxx:899
 EdbPattern.cxx:900
 EdbPattern.cxx:901
 EdbPattern.cxx:902
 EdbPattern.cxx:903
 EdbPattern.cxx:904
 EdbPattern.cxx:905
 EdbPattern.cxx:906
 EdbPattern.cxx:907
 EdbPattern.cxx:908
 EdbPattern.cxx:909
 EdbPattern.cxx:910
 EdbPattern.cxx:911
 EdbPattern.cxx:912
 EdbPattern.cxx:913
 EdbPattern.cxx:914
 EdbPattern.cxx:915
 EdbPattern.cxx:916
 EdbPattern.cxx:917
 EdbPattern.cxx:918
 EdbPattern.cxx:919
 EdbPattern.cxx:920
 EdbPattern.cxx:921
 EdbPattern.cxx:922
 EdbPattern.cxx:923
 EdbPattern.cxx:924
 EdbPattern.cxx:925
 EdbPattern.cxx:926
 EdbPattern.cxx:927
 EdbPattern.cxx:928
 EdbPattern.cxx:929
 EdbPattern.cxx:930
 EdbPattern.cxx:931
 EdbPattern.cxx:932
 EdbPattern.cxx:933
 EdbPattern.cxx:934
 EdbPattern.cxx:935
 EdbPattern.cxx:936
 EdbPattern.cxx:937
 EdbPattern.cxx:938
 EdbPattern.cxx:939
 EdbPattern.cxx:940
 EdbPattern.cxx:941
 EdbPattern.cxx:942
 EdbPattern.cxx:943
 EdbPattern.cxx:944
 EdbPattern.cxx:945
 EdbPattern.cxx:946
 EdbPattern.cxx:947
 EdbPattern.cxx:948
 EdbPattern.cxx:949
 EdbPattern.cxx:950
 EdbPattern.cxx:951
 EdbPattern.cxx:952
 EdbPattern.cxx:953
 EdbPattern.cxx:954
 EdbPattern.cxx:955
 EdbPattern.cxx:956
 EdbPattern.cxx:957
 EdbPattern.cxx:958
 EdbPattern.cxx:959
 EdbPattern.cxx:960
 EdbPattern.cxx:961
 EdbPattern.cxx:962
 EdbPattern.cxx:963
 EdbPattern.cxx:964
 EdbPattern.cxx:965
 EdbPattern.cxx:966
 EdbPattern.cxx:967
 EdbPattern.cxx:968
 EdbPattern.cxx:969
 EdbPattern.cxx:970
 EdbPattern.cxx:971
 EdbPattern.cxx:972
 EdbPattern.cxx:973
 EdbPattern.cxx:974
 EdbPattern.cxx:975
 EdbPattern.cxx:976
 EdbPattern.cxx:977
 EdbPattern.cxx:978
 EdbPattern.cxx:979
 EdbPattern.cxx:980
 EdbPattern.cxx:981
 EdbPattern.cxx:982
 EdbPattern.cxx:983
 EdbPattern.cxx:984
 EdbPattern.cxx:985
 EdbPattern.cxx:986
 EdbPattern.cxx:987
 EdbPattern.cxx:988
 EdbPattern.cxx:989
 EdbPattern.cxx:990
 EdbPattern.cxx:991
 EdbPattern.cxx:992
 EdbPattern.cxx:993
 EdbPattern.cxx:994
 EdbPattern.cxx:995
 EdbPattern.cxx:996
 EdbPattern.cxx:997
 EdbPattern.cxx:998
 EdbPattern.cxx:999
 EdbPattern.cxx:1000
 EdbPattern.cxx:1001
 EdbPattern.cxx:1002
 EdbPattern.cxx:1003
 EdbPattern.cxx:1004
 EdbPattern.cxx:1005
 EdbPattern.cxx:1006
 EdbPattern.cxx:1007
 EdbPattern.cxx:1008
 EdbPattern.cxx:1009
 EdbPattern.cxx:1010
 EdbPattern.cxx:1011
 EdbPattern.cxx:1012
 EdbPattern.cxx:1013
 EdbPattern.cxx:1014
 EdbPattern.cxx:1015
 EdbPattern.cxx:1016
 EdbPattern.cxx:1017
 EdbPattern.cxx:1018
 EdbPattern.cxx:1019
 EdbPattern.cxx:1020
 EdbPattern.cxx:1021
 EdbPattern.cxx:1022
 EdbPattern.cxx:1023
 EdbPattern.cxx:1024
 EdbPattern.cxx:1025
 EdbPattern.cxx:1026
 EdbPattern.cxx:1027
 EdbPattern.cxx:1028
 EdbPattern.cxx:1029
 EdbPattern.cxx:1030
 EdbPattern.cxx:1031
 EdbPattern.cxx:1032
 EdbPattern.cxx:1033
 EdbPattern.cxx:1034
 EdbPattern.cxx:1035
 EdbPattern.cxx:1036
 EdbPattern.cxx:1037
 EdbPattern.cxx:1038
 EdbPattern.cxx:1039
 EdbPattern.cxx:1040
 EdbPattern.cxx:1041
 EdbPattern.cxx:1042
 EdbPattern.cxx:1043
 EdbPattern.cxx:1044
 EdbPattern.cxx:1045
 EdbPattern.cxx:1046
 EdbPattern.cxx:1047
 EdbPattern.cxx:1048
 EdbPattern.cxx:1049
 EdbPattern.cxx:1050
 EdbPattern.cxx:1051
 EdbPattern.cxx:1052
 EdbPattern.cxx:1053
 EdbPattern.cxx:1054
 EdbPattern.cxx:1055
 EdbPattern.cxx:1056
 EdbPattern.cxx:1057
 EdbPattern.cxx:1058
 EdbPattern.cxx:1059
 EdbPattern.cxx:1060
 EdbPattern.cxx:1061
 EdbPattern.cxx:1062
 EdbPattern.cxx:1063
 EdbPattern.cxx:1064
 EdbPattern.cxx:1065
 EdbPattern.cxx:1066
 EdbPattern.cxx:1067
 EdbPattern.cxx:1068
 EdbPattern.cxx:1069
 EdbPattern.cxx:1070
 EdbPattern.cxx:1071
 EdbPattern.cxx:1072
 EdbPattern.cxx:1073
 EdbPattern.cxx:1074
 EdbPattern.cxx:1075
 EdbPattern.cxx:1076
 EdbPattern.cxx:1077
 EdbPattern.cxx:1078
 EdbPattern.cxx:1079
 EdbPattern.cxx:1080
 EdbPattern.cxx:1081
 EdbPattern.cxx:1082
 EdbPattern.cxx:1083
 EdbPattern.cxx:1084
 EdbPattern.cxx:1085
 EdbPattern.cxx:1086
 EdbPattern.cxx:1087
 EdbPattern.cxx:1088
 EdbPattern.cxx:1089
 EdbPattern.cxx:1090
 EdbPattern.cxx:1091
 EdbPattern.cxx:1092
 EdbPattern.cxx:1093
 EdbPattern.cxx:1094
 EdbPattern.cxx:1095
 EdbPattern.cxx:1096
 EdbPattern.cxx:1097
 EdbPattern.cxx:1098
 EdbPattern.cxx:1099
 EdbPattern.cxx:1100
 EdbPattern.cxx:1101
 EdbPattern.cxx:1102
 EdbPattern.cxx:1103
 EdbPattern.cxx:1104
 EdbPattern.cxx:1105
 EdbPattern.cxx:1106
 EdbPattern.cxx:1107
 EdbPattern.cxx:1108
 EdbPattern.cxx:1109
 EdbPattern.cxx:1110
 EdbPattern.cxx:1111
 EdbPattern.cxx:1112
 EdbPattern.cxx:1113
 EdbPattern.cxx:1114
 EdbPattern.cxx:1115
 EdbPattern.cxx:1116
 EdbPattern.cxx:1117
 EdbPattern.cxx:1118
 EdbPattern.cxx:1119
 EdbPattern.cxx:1120
 EdbPattern.cxx:1121
 EdbPattern.cxx:1122
 EdbPattern.cxx:1123
 EdbPattern.cxx:1124
 EdbPattern.cxx:1125
 EdbPattern.cxx:1126
 EdbPattern.cxx:1127
 EdbPattern.cxx:1128
 EdbPattern.cxx:1129
 EdbPattern.cxx:1130
 EdbPattern.cxx:1131
 EdbPattern.cxx:1132
 EdbPattern.cxx:1133
 EdbPattern.cxx:1134
 EdbPattern.cxx:1135
 EdbPattern.cxx:1136
 EdbPattern.cxx:1137
 EdbPattern.cxx:1138
 EdbPattern.cxx:1139
 EdbPattern.cxx:1140
 EdbPattern.cxx:1141
 EdbPattern.cxx:1142
 EdbPattern.cxx:1143
 EdbPattern.cxx:1144
 EdbPattern.cxx:1145
 EdbPattern.cxx:1146
 EdbPattern.cxx:1147
 EdbPattern.cxx:1148
 EdbPattern.cxx:1149
 EdbPattern.cxx:1150
 EdbPattern.cxx:1151
 EdbPattern.cxx:1152
 EdbPattern.cxx:1153
 EdbPattern.cxx:1154
 EdbPattern.cxx:1155
 EdbPattern.cxx:1156
 EdbPattern.cxx:1157
 EdbPattern.cxx:1158
 EdbPattern.cxx:1159
 EdbPattern.cxx:1160
 EdbPattern.cxx:1161
 EdbPattern.cxx:1162
 EdbPattern.cxx:1163
 EdbPattern.cxx:1164
 EdbPattern.cxx:1165
 EdbPattern.cxx:1166
 EdbPattern.cxx:1167
 EdbPattern.cxx:1168
 EdbPattern.cxx:1169
 EdbPattern.cxx:1170
 EdbPattern.cxx:1171
 EdbPattern.cxx:1172
 EdbPattern.cxx:1173
 EdbPattern.cxx:1174
 EdbPattern.cxx:1175
 EdbPattern.cxx:1176
 EdbPattern.cxx:1177
 EdbPattern.cxx:1178
 EdbPattern.cxx:1179
 EdbPattern.cxx:1180
 EdbPattern.cxx:1181
 EdbPattern.cxx:1182
 EdbPattern.cxx:1183
 EdbPattern.cxx:1184
 EdbPattern.cxx:1185
 EdbPattern.cxx:1186
 EdbPattern.cxx:1187
 EdbPattern.cxx:1188
 EdbPattern.cxx:1189
 EdbPattern.cxx:1190
 EdbPattern.cxx:1191
 EdbPattern.cxx:1192
 EdbPattern.cxx:1193
 EdbPattern.cxx:1194
 EdbPattern.cxx:1195
 EdbPattern.cxx:1196
 EdbPattern.cxx:1197
 EdbPattern.cxx:1198
 EdbPattern.cxx:1199
 EdbPattern.cxx:1200
 EdbPattern.cxx:1201
 EdbPattern.cxx:1202
 EdbPattern.cxx:1203
 EdbPattern.cxx:1204
 EdbPattern.cxx:1205
 EdbPattern.cxx:1206
 EdbPattern.cxx:1207
 EdbPattern.cxx:1208
 EdbPattern.cxx:1209
 EdbPattern.cxx:1210
 EdbPattern.cxx:1211
 EdbPattern.cxx:1212
 EdbPattern.cxx:1213
 EdbPattern.cxx:1214
 EdbPattern.cxx:1215
 EdbPattern.cxx:1216
 EdbPattern.cxx:1217
 EdbPattern.cxx:1218
 EdbPattern.cxx:1219
 EdbPattern.cxx:1220
 EdbPattern.cxx:1221
 EdbPattern.cxx:1222
 EdbPattern.cxx:1223
 EdbPattern.cxx:1224
 EdbPattern.cxx:1225
 EdbPattern.cxx:1226
 EdbPattern.cxx:1227
 EdbPattern.cxx:1228
 EdbPattern.cxx:1229
 EdbPattern.cxx:1230
 EdbPattern.cxx:1231
 EdbPattern.cxx:1232
 EdbPattern.cxx:1233
 EdbPattern.cxx:1234
 EdbPattern.cxx:1235
 EdbPattern.cxx:1236
 EdbPattern.cxx:1237
 EdbPattern.cxx:1238
 EdbPattern.cxx:1239
 EdbPattern.cxx:1240
 EdbPattern.cxx:1241
 EdbPattern.cxx:1242
 EdbPattern.cxx:1243
 EdbPattern.cxx:1244
 EdbPattern.cxx:1245
 EdbPattern.cxx:1246
 EdbPattern.cxx:1247
 EdbPattern.cxx:1248
 EdbPattern.cxx:1249
 EdbPattern.cxx:1250
 EdbPattern.cxx:1251
 EdbPattern.cxx:1252
 EdbPattern.cxx:1253
 EdbPattern.cxx:1254
 EdbPattern.cxx:1255
 EdbPattern.cxx:1256
 EdbPattern.cxx:1257
 EdbPattern.cxx:1258
 EdbPattern.cxx:1259
 EdbPattern.cxx:1260
 EdbPattern.cxx:1261
 EdbPattern.cxx:1262
 EdbPattern.cxx:1263
 EdbPattern.cxx:1264
 EdbPattern.cxx:1265
 EdbPattern.cxx:1266
 EdbPattern.cxx:1267
 EdbPattern.cxx:1268
 EdbPattern.cxx:1269
 EdbPattern.cxx:1270
 EdbPattern.cxx:1271
 EdbPattern.cxx:1272
 EdbPattern.cxx:1273
 EdbPattern.cxx:1274
 EdbPattern.cxx:1275
 EdbPattern.cxx:1276
 EdbPattern.cxx:1277
 EdbPattern.cxx:1278
 EdbPattern.cxx:1279
 EdbPattern.cxx:1280
 EdbPattern.cxx:1281
 EdbPattern.cxx:1282
 EdbPattern.cxx:1283
 EdbPattern.cxx:1284
 EdbPattern.cxx:1285
 EdbPattern.cxx:1286
 EdbPattern.cxx:1287
 EdbPattern.cxx:1288
 EdbPattern.cxx:1289
 EdbPattern.cxx:1290
 EdbPattern.cxx:1291
 EdbPattern.cxx:1292
 EdbPattern.cxx:1293
 EdbPattern.cxx:1294
 EdbPattern.cxx:1295
 EdbPattern.cxx:1296
 EdbPattern.cxx:1297
 EdbPattern.cxx:1298
 EdbPattern.cxx:1299
 EdbPattern.cxx:1300
 EdbPattern.cxx:1301
 EdbPattern.cxx:1302
 EdbPattern.cxx:1303
 EdbPattern.cxx:1304
 EdbPattern.cxx:1305
 EdbPattern.cxx:1306
 EdbPattern.cxx:1307
 EdbPattern.cxx:1308
 EdbPattern.cxx:1309
 EdbPattern.cxx:1310
 EdbPattern.cxx:1311
 EdbPattern.cxx:1312
 EdbPattern.cxx:1313
 EdbPattern.cxx:1314
 EdbPattern.cxx:1315
 EdbPattern.cxx:1316
 EdbPattern.cxx:1317
 EdbPattern.cxx:1318
 EdbPattern.cxx:1319
 EdbPattern.cxx:1320
 EdbPattern.cxx:1321
 EdbPattern.cxx:1322
 EdbPattern.cxx:1323
 EdbPattern.cxx:1324
 EdbPattern.cxx:1325
 EdbPattern.cxx:1326
 EdbPattern.cxx:1327
 EdbPattern.cxx:1328
 EdbPattern.cxx:1329
 EdbPattern.cxx:1330
 EdbPattern.cxx:1331
 EdbPattern.cxx:1332
 EdbPattern.cxx:1333
 EdbPattern.cxx:1334
 EdbPattern.cxx:1335
 EdbPattern.cxx:1336
 EdbPattern.cxx:1337
 EdbPattern.cxx:1338
 EdbPattern.cxx:1339
 EdbPattern.cxx:1340
 EdbPattern.cxx:1341
 EdbPattern.cxx:1342
 EdbPattern.cxx:1343
 EdbPattern.cxx:1344
 EdbPattern.cxx:1345
 EdbPattern.cxx:1346
 EdbPattern.cxx:1347
 EdbPattern.cxx:1348
 EdbPattern.cxx:1349
 EdbPattern.cxx:1350
 EdbPattern.cxx:1351
 EdbPattern.cxx:1352
 EdbPattern.cxx:1353
 EdbPattern.cxx:1354
 EdbPattern.cxx:1355
 EdbPattern.cxx:1356
 EdbPattern.cxx:1357
 EdbPattern.cxx:1358
 EdbPattern.cxx:1359
 EdbPattern.cxx:1360
 EdbPattern.cxx:1361
 EdbPattern.cxx:1362
 EdbPattern.cxx:1363
 EdbPattern.cxx:1364
 EdbPattern.cxx:1365
 EdbPattern.cxx:1366
 EdbPattern.cxx:1367
 EdbPattern.cxx:1368
 EdbPattern.cxx:1369
 EdbPattern.cxx:1370
 EdbPattern.cxx:1371
 EdbPattern.cxx:1372
 EdbPattern.cxx:1373
 EdbPattern.cxx:1374
 EdbPattern.cxx:1375
 EdbPattern.cxx:1376
 EdbPattern.cxx:1377
 EdbPattern.cxx:1378
 EdbPattern.cxx:1379
 EdbPattern.cxx:1380
 EdbPattern.cxx:1381
 EdbPattern.cxx:1382
 EdbPattern.cxx:1383
 EdbPattern.cxx:1384
 EdbPattern.cxx:1385
 EdbPattern.cxx:1386
 EdbPattern.cxx:1387
 EdbPattern.cxx:1388
 EdbPattern.cxx:1389
 EdbPattern.cxx:1390
 EdbPattern.cxx:1391
 EdbPattern.cxx:1392
 EdbPattern.cxx:1393
 EdbPattern.cxx:1394
 EdbPattern.cxx:1395
 EdbPattern.cxx:1396
 EdbPattern.cxx:1397
 EdbPattern.cxx:1398
 EdbPattern.cxx:1399
 EdbPattern.cxx:1400
 EdbPattern.cxx:1401
 EdbPattern.cxx:1402
 EdbPattern.cxx:1403
 EdbPattern.cxx:1404
 EdbPattern.cxx:1405
 EdbPattern.cxx:1406
 EdbPattern.cxx:1407
 EdbPattern.cxx:1408
 EdbPattern.cxx:1409
 EdbPattern.cxx:1410
 EdbPattern.cxx:1411
 EdbPattern.cxx:1412
 EdbPattern.cxx:1413
 EdbPattern.cxx:1414
 EdbPattern.cxx:1415
 EdbPattern.cxx:1416
 EdbPattern.cxx:1417
 EdbPattern.cxx:1418
 EdbPattern.cxx:1419
 EdbPattern.cxx:1420
 EdbPattern.cxx:1421
 EdbPattern.cxx:1422
 EdbPattern.cxx:1423
 EdbPattern.cxx:1424
 EdbPattern.cxx:1425
 EdbPattern.cxx:1426
 EdbPattern.cxx:1427
 EdbPattern.cxx:1428
 EdbPattern.cxx:1429
 EdbPattern.cxx:1430
 EdbPattern.cxx:1431
 EdbPattern.cxx:1432
 EdbPattern.cxx:1433
 EdbPattern.cxx:1434
 EdbPattern.cxx:1435
 EdbPattern.cxx:1436
 EdbPattern.cxx:1437
 EdbPattern.cxx:1438
 EdbPattern.cxx:1439
 EdbPattern.cxx:1440
 EdbPattern.cxx:1441
 EdbPattern.cxx:1442
 EdbPattern.cxx:1443
 EdbPattern.cxx:1444
 EdbPattern.cxx:1445
 EdbPattern.cxx:1446
 EdbPattern.cxx:1447
 EdbPattern.cxx:1448
 EdbPattern.cxx:1449
 EdbPattern.cxx:1450
 EdbPattern.cxx:1451
 EdbPattern.cxx:1452
 EdbPattern.cxx:1453
 EdbPattern.cxx:1454
 EdbPattern.cxx:1455
 EdbPattern.cxx:1456
 EdbPattern.cxx:1457
 EdbPattern.cxx:1458
 EdbPattern.cxx:1459
 EdbPattern.cxx:1460
 EdbPattern.cxx:1461
 EdbPattern.cxx:1462
 EdbPattern.cxx:1463
 EdbPattern.cxx:1464
 EdbPattern.cxx:1465
 EdbPattern.cxx:1466
 EdbPattern.cxx:1467
 EdbPattern.cxx:1468
 EdbPattern.cxx:1469
 EdbPattern.cxx:1470
 EdbPattern.cxx:1471
 EdbPattern.cxx:1472
 EdbPattern.cxx:1473
 EdbPattern.cxx:1474
 EdbPattern.cxx:1475
 EdbPattern.cxx:1476
 EdbPattern.cxx:1477
 EdbPattern.cxx:1478
 EdbPattern.cxx:1479
 EdbPattern.cxx:1480
 EdbPattern.cxx:1481
 EdbPattern.cxx:1482
 EdbPattern.cxx:1483
 EdbPattern.cxx:1484
 EdbPattern.cxx:1485
 EdbPattern.cxx:1486
 EdbPattern.cxx:1487
 EdbPattern.cxx:1488
 EdbPattern.cxx:1489
 EdbPattern.cxx:1490
 EdbPattern.cxx:1491
 EdbPattern.cxx:1492
 EdbPattern.cxx:1493
 EdbPattern.cxx:1494
 EdbPattern.cxx:1495
 EdbPattern.cxx:1496
 EdbPattern.cxx:1497
 EdbPattern.cxx:1498
 EdbPattern.cxx:1499
 EdbPattern.cxx:1500
 EdbPattern.cxx:1501
 EdbPattern.cxx:1502
 EdbPattern.cxx:1503
 EdbPattern.cxx:1504
 EdbPattern.cxx:1505
 EdbPattern.cxx:1506
 EdbPattern.cxx:1507
 EdbPattern.cxx:1508
 EdbPattern.cxx:1509
 EdbPattern.cxx:1510
 EdbPattern.cxx:1511
 EdbPattern.cxx:1512
 EdbPattern.cxx:1513
 EdbPattern.cxx:1514
 EdbPattern.cxx:1515
 EdbPattern.cxx:1516
 EdbPattern.cxx:1517
 EdbPattern.cxx:1518
 EdbPattern.cxx:1519
 EdbPattern.cxx:1520
 EdbPattern.cxx:1521
 EdbPattern.cxx:1522
 EdbPattern.cxx:1523
 EdbPattern.cxx:1524
 EdbPattern.cxx:1525
 EdbPattern.cxx:1526
 EdbPattern.cxx:1527
 EdbPattern.cxx:1528
 EdbPattern.cxx:1529
 EdbPattern.cxx:1530
 EdbPattern.cxx:1531
 EdbPattern.cxx:1532
 EdbPattern.cxx:1533
 EdbPattern.cxx:1534
 EdbPattern.cxx:1535
 EdbPattern.cxx:1536
 EdbPattern.cxx:1537
 EdbPattern.cxx:1538
 EdbPattern.cxx:1539
 EdbPattern.cxx:1540
 EdbPattern.cxx:1541
 EdbPattern.cxx:1542
 EdbPattern.cxx:1543
 EdbPattern.cxx:1544
 EdbPattern.cxx:1545
 EdbPattern.cxx:1546
 EdbPattern.cxx:1547
 EdbPattern.cxx:1548
 EdbPattern.cxx:1549
 EdbPattern.cxx:1550
 EdbPattern.cxx:1551
 EdbPattern.cxx:1552
 EdbPattern.cxx:1553
 EdbPattern.cxx:1554
 EdbPattern.cxx:1555
 EdbPattern.cxx:1556
 EdbPattern.cxx:1557
 EdbPattern.cxx:1558
 EdbPattern.cxx:1559
 EdbPattern.cxx:1560
 EdbPattern.cxx:1561
 EdbPattern.cxx:1562
 EdbPattern.cxx:1563
 EdbPattern.cxx:1564
 EdbPattern.cxx:1565
 EdbPattern.cxx:1566
 EdbPattern.cxx:1567
 EdbPattern.cxx:1568
 EdbPattern.cxx:1569
 EdbPattern.cxx:1570
 EdbPattern.cxx:1571
 EdbPattern.cxx:1572
 EdbPattern.cxx:1573
 EdbPattern.cxx:1574
 EdbPattern.cxx:1575
 EdbPattern.cxx:1576
 EdbPattern.cxx:1577
 EdbPattern.cxx:1578
 EdbPattern.cxx:1579
 EdbPattern.cxx:1580
 EdbPattern.cxx:1581
 EdbPattern.cxx:1582
 EdbPattern.cxx:1583
 EdbPattern.cxx:1584
 EdbPattern.cxx:1585
 EdbPattern.cxx:1586
 EdbPattern.cxx:1587
 EdbPattern.cxx:1588
 EdbPattern.cxx:1589
 EdbPattern.cxx:1590
 EdbPattern.cxx:1591
 EdbPattern.cxx:1592
 EdbPattern.cxx:1593
 EdbPattern.cxx:1594
 EdbPattern.cxx:1595
 EdbPattern.cxx:1596
 EdbPattern.cxx:1597
 EdbPattern.cxx:1598
 EdbPattern.cxx:1599
 EdbPattern.cxx:1600
 EdbPattern.cxx:1601
 EdbPattern.cxx:1602
 EdbPattern.cxx:1603
 EdbPattern.cxx:1604
 EdbPattern.cxx:1605
 EdbPattern.cxx:1606
 EdbPattern.cxx:1607
 EdbPattern.cxx:1608
 EdbPattern.cxx:1609
 EdbPattern.cxx:1610
 EdbPattern.cxx:1611
 EdbPattern.cxx:1612
 EdbPattern.cxx:1613
 EdbPattern.cxx:1614
 EdbPattern.cxx:1615
 EdbPattern.cxx:1616
 EdbPattern.cxx:1617
 EdbPattern.cxx:1618
 EdbPattern.cxx:1619
 EdbPattern.cxx:1620
 EdbPattern.cxx:1621
 EdbPattern.cxx:1622
 EdbPattern.cxx:1623
 EdbPattern.cxx:1624
 EdbPattern.cxx:1625
 EdbPattern.cxx:1626
 EdbPattern.cxx:1627
 EdbPattern.cxx:1628
 EdbPattern.cxx:1629
 EdbPattern.cxx:1630
 EdbPattern.cxx:1631
 EdbPattern.cxx:1632
 EdbPattern.cxx:1633
 EdbPattern.cxx:1634
 EdbPattern.cxx:1635
 EdbPattern.cxx:1636
 EdbPattern.cxx:1637
 EdbPattern.cxx:1638
 EdbPattern.cxx:1639
 EdbPattern.cxx:1640
 EdbPattern.cxx:1641
 EdbPattern.cxx:1642
 EdbPattern.cxx:1643
 EdbPattern.cxx:1644
 EdbPattern.cxx:1645
 EdbPattern.cxx:1646
 EdbPattern.cxx:1647
 EdbPattern.cxx:1648
 EdbPattern.cxx:1649
 EdbPattern.cxx:1650
 EdbPattern.cxx:1651
 EdbPattern.cxx:1652
 EdbPattern.cxx:1653
 EdbPattern.cxx:1654
 EdbPattern.cxx:1655
 EdbPattern.cxx:1656
 EdbPattern.cxx:1657
 EdbPattern.cxx:1658
 EdbPattern.cxx:1659
 EdbPattern.cxx:1660
 EdbPattern.cxx:1661
 EdbPattern.cxx:1662
 EdbPattern.cxx:1663
 EdbPattern.cxx:1664
 EdbPattern.cxx:1665
 EdbPattern.cxx:1666
 EdbPattern.cxx:1667
 EdbPattern.cxx:1668
 EdbPattern.cxx:1669
 EdbPattern.cxx:1670
 EdbPattern.cxx:1671
 EdbPattern.cxx:1672
 EdbPattern.cxx:1673
 EdbPattern.cxx:1674
 EdbPattern.cxx:1675
 EdbPattern.cxx:1676
 EdbPattern.cxx:1677
 EdbPattern.cxx:1678
 EdbPattern.cxx:1679
 EdbPattern.cxx:1680
 EdbPattern.cxx:1681
 EdbPattern.cxx:1682
 EdbPattern.cxx:1683
 EdbPattern.cxx:1684
 EdbPattern.cxx:1685
 EdbPattern.cxx:1686
 EdbPattern.cxx:1687
 EdbPattern.cxx:1688
 EdbPattern.cxx:1689
 EdbPattern.cxx:1690
 EdbPattern.cxx:1691
 EdbPattern.cxx:1692
 EdbPattern.cxx:1693
 EdbPattern.cxx:1694
 EdbPattern.cxx:1695
 EdbPattern.cxx:1696
 EdbPattern.cxx:1697
 EdbPattern.cxx:1698
 EdbPattern.cxx:1699
 EdbPattern.cxx:1700
 EdbPattern.cxx:1701
 EdbPattern.cxx:1702
 EdbPattern.cxx:1703
 EdbPattern.cxx:1704
 EdbPattern.cxx:1705
 EdbPattern.cxx:1706
 EdbPattern.cxx:1707
 EdbPattern.cxx:1708
 EdbPattern.cxx:1709
 EdbPattern.cxx:1710
 EdbPattern.cxx:1711
 EdbPattern.cxx:1712
 EdbPattern.cxx:1713
 EdbPattern.cxx:1714
 EdbPattern.cxx:1715
 EdbPattern.cxx:1716
 EdbPattern.cxx:1717
 EdbPattern.cxx:1718
 EdbPattern.cxx:1719
 EdbPattern.cxx:1720
 EdbPattern.cxx:1721
 EdbPattern.cxx:1722
 EdbPattern.cxx:1723
 EdbPattern.cxx:1724
 EdbPattern.cxx:1725
 EdbPattern.cxx:1726
 EdbPattern.cxx:1727
 EdbPattern.cxx:1728
 EdbPattern.cxx:1729
 EdbPattern.cxx:1730
 EdbPattern.cxx:1731
 EdbPattern.cxx:1732
 EdbPattern.cxx:1733
 EdbPattern.cxx:1734
 EdbPattern.cxx:1735
 EdbPattern.cxx:1736
 EdbPattern.cxx:1737
 EdbPattern.cxx:1738
 EdbPattern.cxx:1739
 EdbPattern.cxx:1740
 EdbPattern.cxx:1741
 EdbPattern.cxx:1742
 EdbPattern.cxx:1743
 EdbPattern.cxx:1744
 EdbPattern.cxx:1745
 EdbPattern.cxx:1746
 EdbPattern.cxx:1747
 EdbPattern.cxx:1748
 EdbPattern.cxx:1749
 EdbPattern.cxx:1750
 EdbPattern.cxx:1751
 EdbPattern.cxx:1752
 EdbPattern.cxx:1753
 EdbPattern.cxx:1754
 EdbPattern.cxx:1755
 EdbPattern.cxx:1756
 EdbPattern.cxx:1757
 EdbPattern.cxx:1758
 EdbPattern.cxx:1759
 EdbPattern.cxx:1760
 EdbPattern.cxx:1761
 EdbPattern.cxx:1762
 EdbPattern.cxx:1763
 EdbPattern.cxx:1764
 EdbPattern.cxx:1765
 EdbPattern.cxx:1766
 EdbPattern.cxx:1767
 EdbPattern.cxx:1768
 EdbPattern.cxx:1769
 EdbPattern.cxx:1770
 EdbPattern.cxx:1771
 EdbPattern.cxx:1772
 EdbPattern.cxx:1773
 EdbPattern.cxx:1774
 EdbPattern.cxx:1775
 EdbPattern.cxx:1776
 EdbPattern.cxx:1777
 EdbPattern.cxx:1778
 EdbPattern.cxx:1779
 EdbPattern.cxx:1780
 EdbPattern.cxx:1781
 EdbPattern.cxx:1782
 EdbPattern.cxx:1783
 EdbPattern.cxx:1784
 EdbPattern.cxx:1785
 EdbPattern.cxx:1786
 EdbPattern.cxx:1787
 EdbPattern.cxx:1788
 EdbPattern.cxx:1789
 EdbPattern.cxx:1790
 EdbPattern.cxx:1791
 EdbPattern.cxx:1792
 EdbPattern.cxx:1793
 EdbPattern.cxx:1794
 EdbPattern.cxx:1795
 EdbPattern.cxx:1796
 EdbPattern.cxx:1797
 EdbPattern.cxx:1798
 EdbPattern.cxx:1799
 EdbPattern.cxx:1800
 EdbPattern.cxx:1801
 EdbPattern.cxx:1802
 EdbPattern.cxx:1803
 EdbPattern.cxx:1804
 EdbPattern.cxx:1805
 EdbPattern.cxx:1806
 EdbPattern.cxx:1807
 EdbPattern.cxx:1808
 EdbPattern.cxx:1809
 EdbPattern.cxx:1810
 EdbPattern.cxx:1811
 EdbPattern.cxx:1812
 EdbPattern.cxx:1813
 EdbPattern.cxx:1814
 EdbPattern.cxx:1815
 EdbPattern.cxx:1816
 EdbPattern.cxx:1817
 EdbPattern.cxx:1818
 EdbPattern.cxx:1819
 EdbPattern.cxx:1820
 EdbPattern.cxx:1821
 EdbPattern.cxx:1822
 EdbPattern.cxx:1823
 EdbPattern.cxx:1824
 EdbPattern.cxx:1825
 EdbPattern.cxx:1826
 EdbPattern.cxx:1827
 EdbPattern.cxx:1828
 EdbPattern.cxx:1829
 EdbPattern.cxx:1830
 EdbPattern.cxx:1831
 EdbPattern.cxx:1832
 EdbPattern.cxx:1833
 EdbPattern.cxx:1834
 EdbPattern.cxx:1835
 EdbPattern.cxx:1836
 EdbPattern.cxx:1837
 EdbPattern.cxx:1838
 EdbPattern.cxx:1839
 EdbPattern.cxx:1840
 EdbPattern.cxx:1841
 EdbPattern.cxx:1842
 EdbPattern.cxx:1843
 EdbPattern.cxx:1844
 EdbPattern.cxx:1845
 EdbPattern.cxx:1846
 EdbPattern.cxx:1847
 EdbPattern.cxx:1848
 EdbPattern.cxx:1849
 EdbPattern.cxx:1850
 EdbPattern.cxx:1851
 EdbPattern.cxx:1852
 EdbPattern.cxx:1853
 EdbPattern.cxx:1854
 EdbPattern.cxx:1855
 EdbPattern.cxx:1856
 EdbPattern.cxx:1857
 EdbPattern.cxx:1858
 EdbPattern.cxx:1859
 EdbPattern.cxx:1860
 EdbPattern.cxx:1861
 EdbPattern.cxx:1862
 EdbPattern.cxx:1863
 EdbPattern.cxx:1864
 EdbPattern.cxx:1865
 EdbPattern.cxx:1866
 EdbPattern.cxx:1867
 EdbPattern.cxx:1868
 EdbPattern.cxx:1869
 EdbPattern.cxx:1870
 EdbPattern.cxx:1871
 EdbPattern.cxx:1872
 EdbPattern.cxx:1873
 EdbPattern.cxx:1874
 EdbPattern.cxx:1875
 EdbPattern.cxx:1876
 EdbPattern.cxx:1877
 EdbPattern.cxx:1878
 EdbPattern.cxx:1879
 EdbPattern.cxx:1880
 EdbPattern.cxx:1881
 EdbPattern.cxx:1882
 EdbPattern.cxx:1883
 EdbPattern.cxx:1884
 EdbPattern.cxx:1885
 EdbPattern.cxx:1886
 EdbPattern.cxx:1887
 EdbPattern.cxx:1888
 EdbPattern.cxx:1889
 EdbPattern.cxx:1890
 EdbPattern.cxx:1891
 EdbPattern.cxx:1892
 EdbPattern.cxx:1893
 EdbPattern.cxx:1894
 EdbPattern.cxx:1895
 EdbPattern.cxx:1896
 EdbPattern.cxx:1897
 EdbPattern.cxx:1898
 EdbPattern.cxx:1899
 EdbPattern.cxx:1900
 EdbPattern.cxx:1901
 EdbPattern.cxx:1902
 EdbPattern.cxx:1903
 EdbPattern.cxx:1904
 EdbPattern.cxx:1905
 EdbPattern.cxx:1906
 EdbPattern.cxx:1907
 EdbPattern.cxx:1908
 EdbPattern.cxx:1909
 EdbPattern.cxx:1910
 EdbPattern.cxx:1911
 EdbPattern.cxx:1912
 EdbPattern.cxx:1913
 EdbPattern.cxx:1914
 EdbPattern.cxx:1915
 EdbPattern.cxx:1916
 EdbPattern.cxx:1917
 EdbPattern.cxx:1918
 EdbPattern.cxx:1919
 EdbPattern.cxx:1920

» Last changed: 2017-07-03 09:59 » Last generated: 2017-07-03 09:59
This page has been automatically generated. For comments or suggestions regarding the documentation or ROOT in general please send a mail to ROOT support.
