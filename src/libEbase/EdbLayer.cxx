#include "EdbLog.h"
#include "EdbLayer.h"
#include "EdbPhys.h"

ClassImp(EdbLayer)
ClassImp(EdbCorrectionMap)

using namespace TMath;

///==============================================================================
EdbLayer::EdbLayer()
{
  Set0();
}

///==============================================================================
EdbLayer::EdbLayer(const EdbLayer &l)
{
  Set0();
  Copy(l);
}

///==============================================================================
void EdbLayer::Set0()
{
  eID=0;
  eZ=0; eZmin=0; eZmax=0;
  eX=0; eY=0; eTX=0; eTY=0;
  eDX=eDY=kMaxInt;
  eShr = 1.;
  eZcorr=0;
}

///______________________________________________________________________________
void EdbLayer::Copy(const EdbLayer &l)
{
  eID = l.ID();
  eZ = l.Z(); 
  eZmin=l.Zmin(); 
  eZmax=l.Zmax();
  eX=l.X(); 
  eY=l.Y(); 
  eTX=l.TX(); 
  eTY=l.TY();
  eDX=l.DX(); 
  eDY=l.DY();
  CopyCorr(l);
  eNcp = l.eNcp;
  eMap.Copy(l.eMap);
}

///______________________________________________________________________________
void EdbLayer::CopyCorr(const EdbLayer &l)
{
  ResetCorr();
  ApplyCorrections( l.Shr(), l.Zcorr(), *(l.AffineXY()),  *(l.AffineTXTY()));
}

///______________________________________________________________________________
void EdbLayer::Invert()
{
  eAffXY.Invert();
  eAffTXTY.Invert();
  eShr = 1./eShr;
  //TODO: invert also zcorr?
}

///______________________________________________________________________________
void EdbLayer::CorrectSeg( EdbSegP &s)
{
  s.Set( s.ID(), X(s), Y(s), TX(s), TY(s), s.W(), s.Flag() );
}

///______________________________________________________________________________
void EdbLayer::CorrectSegLocal( EdbSegP &s)
{
  eMap.CorrectSeg(s);
}

///______________________________________________________________________________
void EdbLayer::ApplyCorrections(EdbLayer &la)
{
  ApplyCorrections( la.Shr(), la.Zcorr(), *(la.GetAffineXY()), *(la.GetAffineTXTY()) );
  SetNcp(la.Ncp());
}

///______________________________________________________________________________
void EdbLayer::ApplyCorrections(float shr, float zcorr, const EdbAffine2D &affxy, const EdbAffine2D &afftxty)
{
  eShr *= shr;
  eZcorr += zcorr;
  Log(3,"EdbLayer::ApplyCorrections","eShr: %f eZcorr: %f ", eShr, eZcorr );
  Log(3,"EdbLayer::ApplyCorrections","aff before: %s", eAffXY.AsString());
  Log(3,"EdbLayer::ApplyCorrections","corr      : %s", affxy.AsString());
  eAffXY.Transform(affxy);
  Log(3,"EdbLayer::ApplyCorrections","aff after: %s", eAffXY.AsString());
  Log(3,"EdbLayer::ApplyCorrections","afftxty  : %s", afftxty.AsString());
  eAffTXTY.Transform(afftxty);
}

///______________________________________________________________________________
void EdbLayer::SubstructCorrections(EdbLayer &l)
{
  // apply inverse corrections of layer l
  EdbAffine2D aff( *(l.GetAffineXY()) ), afftxty( *(l.GetAffineTXTY()) );
  aff.Invert();
  afftxty.Invert();
  ApplyCorrections( 1./l.Shr(), -l.Zcorr() , aff, afftxty );
}

///______________________________________________________________________________
void EdbLayer::ResetCorr()
{
  eShr = 1.;
  eZcorr=0;
  eAffXY.Reset();
  eAffTXTY.Reset();
}

///______________________________________________________________________________
void EdbLayer::ShiftZ(float dz)
{
  eZ += dz;  eZmin +=dz;  eZmax +=dz;
}

///______________________________________________________________________________
bool EdbLayer::IsInside(float x, float y, float z)
{
  if( x<Xmin() )  return false;
  if( x>Xmax() )  return false;
  if( y<Ymin() )  return false;
  if( y>Ymax() )  return false;
  if( z<Zmin() )  return false;
  if( z>Zmax() )  return false;
  return true;
}

///______________________________________________________________________________
bool EdbLayer::IsInside(float x, float y)
{
  if( x<Xmin() )  return false;
  if( x>Xmax() )  return false;
  if( y<Ymin() )  return false;
  if( y>Ymax() )  return false;
  return true;
}

///______________________________________________________________________________
void EdbLayer::Print()
{
  printf("Layer:\t%d\n",eID);
  printf("Ncp:\t%d\n", eNcp );
  printf("ZLAYER\t%f %f %f\n",eZ,eZmin,eZmax);
  printf("SHR\t%f\n",eShr);
  printf("Zcorr\t%f\n",eZcorr);
  printf("AFFXY\t");        eAffXY.Print();
  printf("AFFTXTY\t");      eAffTXTY.Print();
  printf("Local Corrections Table: %d x %d \n", eMap.NX(), eMap.NY() );
}

//=====================================================================================================
EdbCorrectionMap::~EdbCorrectionMap()
{
}

//----------------------------------------------------------------------------------------------------
void EdbCorrectionMap::Init( EdbCell2 &c )
{
  Init(c.NX(),c.Xmin(),c.Xmax(),c.NY(),c.Ymin(),c.Ymax());
}

//-----------------------------------------------------------------------------------------------------
void EdbCorrectionMap::Init(  int nx, float minx, float maxx, int ny, float miny, float maxy )
{
  EdbCell2::InitCell( nx, minx, maxx, ny, miny, maxy, 1);
  for(int i=0; i<Ncell(); i++) AddObject(i, new EdbLayer() );
}

//-----------------------------------------------------------------------------------------------------
void EdbCorrectionMap::CorrectSeg( EdbSegP &s)
{
  EdbLayer *loc = GetLayer(s.X(), s.Y());
  if(loc) loc->CorrectSeg(s);
}

//-----------------------------------------------------------------------------------------------------
void EdbCorrectionMap::PrintDZ()
{
  int n=Ncell();
  for(int i=0; i<n; i++) printf("%d  dz=%7.2f\n",i,GetLayer(i)->Zcorr());
}

//-----------------------------------------------------------------------------------------------------
void EdbCorrectionMap::ApplyCorrections(EdbCorrectionMap &map)
{
  if( !map.Ncell() ) return;
  if( NX() != map.NX() || 
      NY() != map.NY() ) {
      Log(1,"EdbCorrectionMap::ApplyCorrections","Warning: incompatible maps: %dx%d  vs  %dx%d   Reset map to new one!", NX(),NY(), map.NX(), map.NY() );
      Init( map );
  }
  for(int i=0; i<Ncell(); i++) {
    EdbLayer *loc  = GetLayer(i);
    EdbLayer *locD = map.GetLayer(i);
    loc->ApplyCorrections( *locD );
  }
}

//-----------------------------------------------------------------------------------------------------
int EdbCorrectionMap::Ncp()
{
  int ncp=0;
  for(int i=0; i<Ncell(); i++) {
    EdbLayer *loc  = GetLayer(i);
    if(loc) ncp+=loc->Ncp();
  }
  return ncp;
}

/*
///______________________________________________________________________________
EdbSegP *EdbCorrectionMap::CorrLoc(int j)
{
  EdbSegP *s = new EdbSegP();
  s->Set(j,Xj(j), Yj(j), 0, 0, 0, 0);
  GetLayer(j)->CorrectSeg(*s);
  s->SetX( s->X() - Xj(j) );
  s->SetY( s->Y() - Yj(j) );
  s->SetZ(GetLayer(j)->Zcorr());
  return s;
}
*/

///______________________________________________________________________________
EdbSegCorr EdbCorrectionMap::CorrLoc(int j)
{
  EdbSegP s;
  float x0=Xj(j), y0=Yj(j);
  s.Set(j, x0,y0, 0, 0, 0, 0);
  EdbLayer *la = GetLayer(j);
  la->CorrectSeg(s);
  s.SetX( s.X() - x0 );
  s.SetY( s.Y() - y0 );
  s.SetZ(la->Zcorr());
  EdbSegCorr corr;
  corr.SetV( 0, s.X());
  corr.SetV( 1, s.Y());
  corr.SetV( 2, s.Z());
  corr.SetV( 3, s.TX());
  corr.SetV( 4, s.TY());
  corr.SetV( 5, la->Shr());
  return corr;
}

///______________________________________________________________________________
EdbSegCorr EdbCorrectionMap::CorrLoc(float x, float y)
{
  EdbSegCorr corr;
//TODO
  return corr;
}
