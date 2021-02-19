//-- Author :  Valeri Tioukov   18/12/2011
// read tracks tree, analyse it and produce the report

#include <string.h>
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TEnv.h>
#include <TMath.h>
#include <TRint.h>
#include <TH1F.h>
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbAlignmentV.h"
#include "EdbTrackFitter.h"
#include "EdbCorrectionMapper.h"
#include "EdbCouplesTree.h"
#include "EdbMomentumEstimator.h"

using namespace std;
using namespace TMath;

bool MakeCorrectionMap(EdbPVRec &ali, TEnv &cenv);
void DoGlobalCorr(EdbPVRec &ali, TEnv &cenv);
void DoGlobalCorrV2(EdbPVRec &ali, TEnv &cenv);
void GlobalDiff(EdbPVRec &ali, const char *suff);
void GlobalEff(EdbPVRec &ali, const char *suff);
void CheckGap(EdbPVRec &ali, const char *suff);
void CheckMom(EdbPVRec &ali, TEnv &cenv);
bool ReadFunParameters(const char* key);
void MakeCorrectionMap_test(EdbPVRec &ali, TEnv &cenv);

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout<< "\nUsage: \n\t  emtrackan -set=ID [-corraff -divide=NXxNY -o=DATA_DIRECTORY -v=DEBUG] \n";
  cout<< "\t  emtrackan -file=File [-o=DATA_DIRECTORY -v=DEBUG] \n\n";
  cout<< "\t  Analyse the tracks tree and produce the report file: ID.trk.an.root\n";
  cout<< "\t\t -corraff    -  save corrections to set.root file\n";
  cout<< "\t\t -global     -  global corrections using long tracks to unbend the full set\n";
  cout<< "\t\t -globV2     -  global corrections using long tracks to unbend the full set NEW!\n";
  cout<< "\t\t -momentum   -  tracks momentum estimation\n";
  cout<< "\t\t -divide=NxM -  divide set into into NxM zones and calculate local correction in each zone\n";
  cout<< "\t\t                the following tracking may use this correction if the parameter fedra.track.do_local_corr in track.rootrc is set to 1\n";
  cout<< "\t\t                normally is required several iterations (>3) like emtra.../emtrackan... to converge the local corrections. \n";
  cout<< "\t\t                Map granularity NxM must be the same for all iterations\n";
  cout<< "\n\t Example: emtrackan -set=700000.0.1.0 -divide=10x10 -corraff -v=2 \n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default(TEnv &cenv)
{
  // default parameters for the new linking
  
  cenv.SetValue("trackan.outdir"          , "..");
  cenv.SetValue("trackan.read_cut"        , "1");
  cenv.SetValue("trackan.NCPmin"          , "50");
  cenv.SetValue("trackan.EdbDebugLevel"   ,   1 );
  cenv.SetValue("trackan.DoRefitLine"            , 1 );
  cenv.SetValue("trackan.global.doXYcorr"        , 1);
  cenv.SetValue("trackan.global.doTXTYcorr"      , 1);
  cenv.SetValue("trackan.global.doZcorr"         , 1);
  cenv.SetValue("trackan.MomEst.Alg",0);
  cenv.SetValue("trackan.MomEst.DT", "0.0021 0.0054 0");
  cenv.SetValue("trackan.MomEst.DTx","0.0021 0.0093 0");
  cenv.SetValue("trackan.MomEst.DTy","0.0021 0 0");
  cenv.SetValue("trackan.MomEst.X0", 5600);
  cenv.SetValue("trackan.MomEst.M", 0.13957);
  cenv.SetValue("trackan.global.trFlag", 0);
  cenv.SetValue("trackan.global.noScale", 0);
}

bool      do_set      = false;
bool      do_file     = false;
bool      do_corraff  = false;
bool      do_momentum = false;
EdbID idset;
EdbScanProc sproc;

//----------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("trackanenv");
  set_default(cenv);
  gEDBDEBUGLEVEL        = cenv.GetValue("trackan.EdbDebugLevel" ,  1  );
  const char *outdir    = cenv.GetValue("trackan.outdir"        , "..");

  const char *name=0;
  
  bool do_global=false;
  bool do_globV2=false;
  bool do_cmtest=false;
  
  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];
    if(!strncmp(key,"-set=",5))
    {
      if(strlen(key)>5)	if(idset.Set(key+5))   do_set=true;
    }
    if(!strncmp(key,"-file=",6))
    {
      if(strlen(key)>6) { name=(key+6);   do_file=true;}
    }
    if(!strncmp(key,"-divide=",8))
    {
      if(strlen(key)>8) {
        int nx=0,ny=0; 
        if(2 == sscanf(key+8,"%dx%d",&nx,&ny) ) {
          cenv.SetValue("trackan.NX" , nx);
          cenv.SetValue("trackan.NY" , ny);
        }
      }
    }
    if(!strncmp(key,"-corraff",8))
    {
      do_corraff=true;
    }
    else if(!strncmp(key,"-momentum",9))
    {
      do_momentum=true;
    }
    else if(!strncmp(key,"-global",7))
    {
      do_global=true;
    }
    else if(!strncmp(key,"-globV2",7))
    {
      do_globV2=true;
    }
    else if(!strncmp(key,"-cmtest",7))
    {
      do_cmtest=true;
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
  }
  
  sproc.eProcDirClient=outdir;
  cenv.ReadFile( "trackan.rootrc" ,kEnvLocal);
  //cenv.SetValue("emtra.outdir"               , outdir);
  cenv.WriteFile("trackan.save.rootrc");
  
  const char *cut = cenv.GetValue("trackan.read_cut"        , "1");

  EdbPVRec ali;
  if(do_set)       sproc.ReadTracksTree( idset,ali, cut );
  else if(do_file) sproc.ReadTracksTree( name ,ali, cut );
  //ali.Print();
  
  if(do_momentum)     CheckMom(ali,cenv);
  else if(do_global)  DoGlobalCorr(ali,cenv);
  else if(do_globV2)  DoGlobalCorrV2(ali,cenv);
  else if(do_cmtest)  MakeCorrectionMap_test(ali,cenv);
  else                MakeCorrectionMap(ali,cenv);

  //EdbDataProc::MakeTracksTree( *(ali.eTracks), 0., 0., "al.trk.root" );

  return 0;
}

//------------------------------------------------------------------------------
bool MakeCorrectionMap(EdbPVRec &ali, TEnv &cenv)
{
  float xmin = ali.Xmin();
  float xmax = ali.Xmax();
  float ymin = ali.Ymin();
  float ymax = ali.Ymax();
  
  EdbPattern *patMin = ali.GetPatternZLowestHighest(1);   // get pattern with lowest z
  EdbPattern *patMax = ali.GetPatternZLowestHighest(0);   // get pattern with highest z
  
  bool direction = 0;    // 1-increaseZ, 2-decreaseZ
  int idstart,idend, step;
  if(direction) {
    idstart = patMin->ID();
    idend   = patMax->ID();
  }
  else {
    idstart = patMax->ID();
    idend   = patMin->ID(); 
  }
  step    = idstart<idend?1:-1;
  
  int ref_pl =  ali.GetPattern(idend)->ScanID().ePlate;

  int npat = ali.Npatterns();
  if(npat<2)  return 0; 
  if(Abs(idend-idstart)+1 != npat) return 0;
  int idmap[npat]; int j=idstart;
  for(int i=0; i<npat; i++) { 
    idmap[i]=j; j+=step; 
  }
  
  EdbScanSet *ss =0;
  if(do_set) ss=sproc.ReadScanSet(idset);

  int NX = cenv.GetValue("trackan.NX" , 0);
  int NY = cenv.GetValue("trackan.NY" , 0);
  
  int NCPMIN = cenv.GetValue("trackan.NCPmin" , 50);

  EdbCorrectionMapper cmtmp[npat-1];
  EdbCorrectionMapper cm[npat-1];

  for(int i=0; i<npat-1; i++) {
    EdbPattern *p1 = ali.GetPattern(idmap[i]);
    EdbPattern *p2 = ali.GetPattern(idmap[i+1]);
    cm[i].eID1 = p1->ScanID();
    cm[i].eID2 = p2->ScanID();
    cm[i].eZ1  = p1->Z();
    cm[i].eZ2  = p2->Z();
    cm[i].eNcpMin = NCPMIN;
    
    if(NX*NY) cm[i].InitMap( NX,xmin, xmax, NY, ymin,ymax);
    if(NX*NY) cmtmp[i].InitMap( NX,xmin, xmax, NY, ymin,ymax);
    
    cm[i].eHdty_ty.InitH2( 50,-2,2, 50,-0.1,0.1 );
    cm[i].eHshr.InitH1( 100, 0.8, 1.2 );

    cm[i].eHdxy.InitH2( 61,-60.5,60.5, 61,-60.5,60.5 );
    cm[i].eHdtxy.InitH2( 41,-0.105,0.105, 41,-0.105,0.105 );
    cm[i].eHdz.InitH1( 50,-100,100 );
    
    cm[i].eHxy1.InitH2(100,xmin, xmax, 100, ymin,ymax);
    cm[i].eHtxty1.InitH2(100,-2, 2, 100, -2,2);
    cm[i].eHxy2.InitH2(100,xmin, xmax, 100, ymin,ymax);
    cm[i].eHtxty2.InitH2(100,-2, 2, 100, -2,2);
  }


  int ntr  = ali.Ntracks();

  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = ali.GetTrack(i);
    int nseg = t->N();
    if(nseg<2) continue;
    
    for(int j=0; j<nseg-1; j++) {
      EdbSegP *s1,*s2;
      if(direction) {
        s1 = t->GetSegment(j);
        s2 = t->GetSegment(j+1);
      }
      else {
        s2 = t->GetSegment(j);
        s1 = t->GetSegment(j+1);
      }
      if( Abs(s2->PID()-s1->PID())!=1 )  continue;

      cm[idmap[s1->PID()]].Fill(*s1,*s2);                // fill general hists
      cmtmp[idmap[s1->PID()]].AddSegCouple(s1,s2);
    }
  }
  
  //EnrichPoorBins:
  for(int i=0; i<npat-1; i++) {
    int nbin = cmtmp[i].eMapAl.Ncell();
    for(int j=0; j<nbin; j++) {
      EdbCorrectionBin *bin  = cm[i].GetBin(j);
      EdbCorrectionBin *bint = cmtmp[i].GetBin(j);

      cm[i].AddBin(*bint); bin->AddBin(*bint);             // TODO enrichment algorithm is disabled now
/*      if(bint->eAl.Ncp() >= NCPMIN)      { cm[i].AddBin(*bint); bin->AddBin(*bint); }
      else      {
        cm[i].AddBin(*bint); 
        TObjArray arr;
        int ne = cmtmp[i].eMapAl.SelectObjectsCJ(j,1,arr);
        for(int ie=0; ie<ne; ie++) {
          EdbCorrectionBin *bine = (EdbCorrectionBin *)(arr.At(ie));
          bin->AddBin(*bine);
        }
      }*/
    }
    if(nbin<1) cm[i].AddBin( cmtmp[i] );
  }

  for(int i=0; i<npat-1; i++) {
    cm[i].SetCorrDZ();
    cm[i].CalculateCorrections();
    cm[i].MakeCorrectionsTable();
  }
  

  TString name;
  sproc.MakeFileName(name,idset,"trk.an.root",false);
  Log(2,"MakeCorrectionMap","%s",name.Data());
  TFile f(name,"RECREATE");
  
  if(ss) {
    for(int i=0; i<npat-1; i++) {
	      ss->eReferencePlate=ref_pl;
	
	      if(NX*NY) cm[i].UpdateLayerWithLocalCorr( cm[i].eLayer );
        
        cm[i].eLayer.Print();
        
        ss->UpdateBrickWithP2P( cm[i].eLayer, cm[i].eID1.ePlate, cm[i].eID2.ePlate );
        
        gROOT->SetBatch();
        TCanvas *csum = cm[i].DrawSum(Form("%d",i));
        csum->Write();
        TCanvas *cdiff = cm[i].DrawMap(cm[i].eMap, Form("diff_%d",i));
        cdiff->Write();
        TCanvas *cabs = cm[i].DrawMap( ss->GetPlate(cm[i].eID1.ePlate)->Map(), Form("%d",i));
        cabs->Write();

    }
  }
  f.Close();

  printf("********************write set\n");
  
  if(do_corraff) if(ss) sproc.WriteScanSet(idset,*ss);

  return 1;
}

//------------------------------------------------------------------------------
void MakeCorrectionMap_test(EdbPVRec &ali, TEnv &cenv)
{
  int NCPMIN = cenv.GetValue("trackan.NCPmin" , 50);
  int npat = ali.Npatterns();
  int ntr  = ali.Ntracks();
  int   NX = cenv.GetValue("trackan.NX" , 0);
  int   NY = cenv.GetValue("trackan.NY" , 0);
  Log(2,"MakeCorrectionMap_test","with %d tracks and %d patterns map: %d x %d",ntr, npat, NX,NY);
  if(npat<2)     return;
  if(ntr<NCPMIN) return;
  if(!NX*NY)     return;
  
  EdbScanSet *ss = sproc.ReadScanSet(idset);
  struct PatPlat{
    EdbPattern *pat;
    EdbPlateP  *pl;
    EdbCorrectionMapper cm;
  };
  PatPlat sequence[npat];
  for( int ipat=0; ipat<npat; ipat++ )  //TODO: select direction here
  {
    sequence[ipat].pat = ali.GetPattern(ipat); sequence[ipat].pat->SetPID(ipat);
    sequence[ipat].pl  = ss->GetPlate( sequence[ipat].pat->ScanID().ePlate );
  }

  float xmin = ali.Xmin();
  float xmax = ali.Xmax();
  float ymin = ali.Ymin();
  float ymax = ali.Ymax();
  bool trFlag     = cenv.GetValue("trackan.global.trFlag"          , 0);

  for(int i=0; i<npat-1; i++) {
    EdbPattern *p1          = sequence[i].pat;
    EdbPattern *p2          = sequence[i+1].pat;
    EdbCorrectionMapper &cm = sequence[i].cm;
    cm.eID1 = p1->ScanID();
    cm.eID2 = p2->ScanID();
    cm.eZ1  = p1->Z();
    cm.eZ2  = p2->Z();
    cm.eNcpMin = NCPMIN;
    
    cm.InitMap( NX,xmin, xmax, NY, ymin,ymax);
    /*
    cm[i].eHdty_ty.InitH2( 50,-2,2, 50,-0.1,0.1 );
    cm[i].eHshr.InitH1( 100, 0.8, 1.2 );
    cm[i].eHdxy.InitH2( 61,-60.5,60.5, 61,-60.5,60.5 );
    cm[i].eHdtxy.InitH2( 41,-0.105,0.105, 41,-0.105,0.105 );
    cm[i].eHdz.InitH1( 50,-100,100 );
    cm[i].eHxy1.InitH2(100,xmin, xmax, 100, ymin,ymax);
    cm[i].eHtxty1.InitH2(100,-2, 2, 100, -2,2);
    cm[i].eHxy2.InitH2(100,xmin, xmax, 100, ymin,ymax);
    cm[i].eHtxty2.InitH2(100,-2, 2, 100, -2,2);
    */
    for(int itr=0; itr<ntr; itr++) {
      EdbTrackP *t = ali.GetTrack(itr);
      if(trFlag!=0) if( trFlag!=t->Flag() )  continue;
      EdbSegP *s1=0, *s2=0;
      for(int j=0; j<t->N(); j++) {
        EdbSegP *s  = t->GetSegment(j);
        if(s->PID()==p1->PID()) s1=s;
        if(s->PID()==p2->PID()) s2=s;
      }
      if(s1&&s2) { cm.AddSegCouple(s1,s2); }
    }
  }
  sequence[npat-1].cm.InitMap( NX,xmin, xmax, NY, ymin,ymax);

  for(int i=0; i<npat-1; i++) 
  {
    int nbin = sequence[i].cm.eMapAl.Ncell();
    for(int j=0;j<nbin; j++) 
    {
      printf("%s -> %s dz: %f \n", 
             sequence[i].cm.eID1.AsString(), 
             sequence[i].cm.eID2.AsString(), 
             sequence[i].cm.eZ2-sequence[i].cm.eZ1
            );
      EdbAlignmentV  &al  = sequence[i].cm.GetBin(j)->eAl;
      EdbAlignmentV  &aln = sequence[i+1].cm.GetBin(j)->eAl;
      if(al.Ncp()>=NCPMIN)
      {
        al.CalcAffFull();   // update coord of side0 and ang of side1
        aln.eCorrL[0].GetAffineTXTY()->Transform( al.eCorrL[1].GetAffineTXTY() );
      }
    }
  }

  
  for(int j=0;j<NX*NY; j++)
  {
    for( int i=0; i<npat; i++ ) {
      EdbAlignmentV  &al = sequence[i].cm.GetBin(j)->eAl;
      EdbLayer *l = sequence[i].cm.eMap.GetLayer(j);
      
      float dzcorr= al.eCorrL[0].Zcorr()- (sequence[i].cm.eZ2-sequence[i].cm.eZ1);
      l->GetAffineTXTY()->Transform(  al.eCorrL[0].GetAffineTXTY() );
      l->SetNcp(  al.Ncp() );

      for( int ii=0; ii<=i; ii++ ) {
        EdbLayer *ll = sequence[ii].cm.eMap.GetLayer(j);
        ll->GetAffineXY()->Transform(  al.eCorrL[0].GetAffineXY() );
        ll->SetZcorr( ll->Zcorr() + dzcorr  );
      }
    }
  }

  for( int i=0; i<npat; i++ ) {
    printf( Form("\n************** %s -> %s\n", 
            sequence[i].cm.eID1.AsString(), 
            sequence[i].cm.eID2.AsString()) );
    sequence[i].pl->ApplyCorrectionsLocal(sequence[i].cm.eMap);
    //sequence[i].pl->Print();
  }

  if(do_corraff) if(ss) sproc.WriteScanSet(idset,*ss);

  //ss->Print();
  /*
  TString name;
  sproc.MakeFileName(name,idset,"trk.an.root",false);
  Log(2,"MakeCorrectionMap","%s",name.Data());
  TFile f(name,"RECREATE");
  
  if(ss) {
    for(int i=0; i<npat-1; i++) {
      ss->eReferencePlate=ref_pl;
  
      if(NX*NY) cm[i].UpdateLayerWithLocalCorr( cm[i].eLayer );
        
      cm[i].eLayer.Print();
        
      ss->UpdateBrickWithP2P( cm[i].eLayer, cm[i].eID1.ePlate, cm[i].eID2.ePlate );
        
      gROOT->SetBatch();
      TCanvas *csum = cm[i].DrawSum(Form("%d",i));
      csum->Write();
      TCanvas *cdiff = cm[i].DrawMap(cm[i].eMap, Form("diff_%d",i));
      cdiff->Write();
      TCanvas *cabs = cm[i].DrawMap( ss->GetPlate(cm[i].eID1.ePlate)->Map(), Form("%d",i));
      cabs->Write();
    }
  }
  f.Close();

  printf("********************write set\n");
  
  if(do_corraff) if(ss) sproc.WriteScanSet(idset,*ss);
  */
}

//------------------------------------------------------------------------------
void DoGlobalCorr(EdbPVRec &ali, TEnv &cenv)
{
  int ntr  = ali.Ntracks();
  int npat = ali.Npatterns();
  Log(2,"DoGlobalCorr","with %d tracks and %d patterns",ntr, npat);
  
  cenv.Print();
  bool doRefitLine= cenv.GetValue("trackan.DoRefitLine"            , 1);
  bool doXYcorr   = cenv.GetValue("trackan.global.doXYcorr"        , 1);
  bool doTXTYcorr = cenv.GetValue("trackan.global.doTXTYcorr"      , 1);
  bool doZcorr    = cenv.GetValue("trackan.global.doZcorr"         , 1);
  bool trFlag     = cenv.GetValue("trackan.global.trFlag"          , 0);
  bool noScale    = cenv.GetValue("trackan.global.noScale"         , 0);


  GlobalDiff(ali,"before_fit");
  CheckGap(ali,"before_fit");
  
  if(doRefitLine) {
    EdbTrackFitter fitter;
    for(int i=0; i<ntr; i++) {
      EdbTrackP *t = ali.GetTrack(i);
      fitter.FitTrackLine(*t);
    }
  }
  
  GlobalDiff(ali,"before_corr");

  EdbCorrectionMapper cm[npat];

  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = ali.GetTrack(i);
    if(trFlag!=0) if( trFlag!=t->Flag() )  continue;
    int nseg = t->N();
    for(int j=0; j<nseg; j++) {
      EdbSegP *s  = t->GetSegment(j);
      EdbSegP *sf = t->GetSegmentF(j);
      cm[s->PID()].FillS(s,sf);
    }
  }
  
  if(doXYcorr) {
    for(int i=0; i<npat; i++) {
      int n=0;
      if(noScale) n = cm[i].eAl.CalculateAffXYTurn( *(cm[i].eAl.eCorrL[0].GetAffineXY()) );
      else        n = cm[i].eAl.CalculateAffXY( *(cm[i].eAl.eCorrL[0].GetAffineXY()) );
      Log(2,"CalculateAffXY","%6d  %s", n, cm[i].eAl.eCorrL[0].GetAffineXY()->AsString() );
    }
  
    for(int i=0; i<npat; i++) {
      ali.GetPattern(i)->Transform( cm[i].eAl.eCorrL[0].GetAffineXY() );
    }
    GlobalDiff(ali,"afterXY");
  }

  if(doTXTYcorr) {
    for(int i=0; i<npat; i++) {
      int n = cm[i].eAl.CalculateAffTXTY( *(cm[i].eAl.eCorrL[0].GetAffineTXTY()) );
      Log(2,"CalculateAffTXTY","%6d  %s", n, cm[i].eAl.eCorrL[0].GetAffineTXTY()->AsString() );
    }
  
    for(int i=0; i<npat; i++) {
      ali.GetPattern(i)->TransformA( cm[i].eAl.eCorrL[0].GetAffineTXTY() );
    }
  
    GlobalDiff(ali,"afterTXTY");
  }
  
  if(doZcorr) {
    GlobalDiff(ali,"beforeDZ");
    for(int i=0; i<npat; i++)
    {
      float dz = cm[i].eAl.FineCorrZ();
      cm[i].eAl.eCorrL[0].SetZcorr(dz);
      printf("%d oldz = %f  dz = %f \n",i, ali.GetPattern(i)->Z(), dz );
      ali.GetPattern(i)->ProjectTo(dz);
    }
    GlobalDiff(ali,"afterDZ");
  }
  
  int plate_pid[npat];
  for(int i=0; i<npat; i++) {
    plate_pid[i] = ali.GetPattern(i)->ScanID().ePlate;
  }
  
  EdbScanSet *ss =0;
  if(do_set) ss=sproc.ReadScanSet(idset);
  
  for(int i=0; i<npat; i++) {
    EdbPlateP *plate = ss->GetPlate( plate_pid[i] );
    plate->ApplyCorrections( cm[i].eAl.eCorrL[0] );
  }
  
  GlobalEff(ali,"a");

  if(do_corraff) if(ss) sproc.WriteScanSet(idset,*ss);

}


//------------------------------------------------------------------------------
void DoGlobalCorrV2(EdbPVRec &ali, TEnv &cenv)
{
  int ntr  = ali.Ntracks();
  int npat = ali.Npatterns();
  Log(2,"DoGlobalCorrV2","with %d tracks and %d patterns",ntr, npat);
  if(npat<2) return;
  
  cenv.Print();
  bool doXYcorr   = cenv.GetValue("trackan.global.doXYcorr"        , 1);
  bool doTXTYcorr = cenv.GetValue("trackan.global.doTXTYcorr"      , 1);
  bool doZcorr    = cenv.GetValue("trackan.global.doZcorr"         , 1);
  bool trFlag     = cenv.GetValue("trackan.global.trFlag"          , 0);
  bool noScale    = cenv.GetValue("trackan.global.noScale"         , 0);

  EdbScanSet *ss = sproc.ReadScanSet(idset);
  
  int niterset=1;
  for( int iii=0; iii<2*niterset; iii++ )
  {
    printf("\n ******** iii = %d *********\n",iii);
    bool go_forward = 1;
    //bool go_forward = iii%2;
    for( int ipat=1; ipat<npat; ipat++ )
    {
      int id1,id2;
      if(go_forward) { id1 = ipat-1;    id2 = ipat; }
      else           { id1 = npat-ipat; id2 = npat-ipat-1; }
      
      EdbPattern *p1 = ali.GetPattern(id1);  p1->SetPID(id1);
      EdbPattern *p2 = ali.GetPattern(id2);  p2->SetPID(id2);
      printf("plate: %d->%d pid: %d->%d  z: %.1f ->%.1f = %.1f \n", p1->ScanID().ePlate, p2->ScanID().ePlate,p1->PID(), p2->PID(), p1->Z(), p2->Z(), p2->Z()-p1->Z() );
      TObjArray p1corr;
      TObjArray p2corr;
      int ncp=0;
      for(int itr=0; itr<ntr; itr++) {
        EdbTrackP *t = ali.GetTrack(itr);
        if(trFlag!=0) if( trFlag!=t->Flag() )  continue;
        EdbSegP *s1=0, *s2=0;
        for(int j=0; j<t->N(); j++) {
          EdbSegP *s  = t->GetSegment(j);
          if(s->PID()==p1->PID()) s1=s;
          if(s->PID()==p2->PID()) s2=s;
        }
        if(s1&&s2) { p1corr.Add(s1); p2corr.Add(s2); ncp++; }
      }
    
      EdbAlignmentV al;
      for(int i=0; i<ncp; i++)
      {
        al.eS[0].Add(p1corr.At(i));
        al.eS[1].Add(p2corr.At(i));
      }
    
      EdbAffine2D aff1to2XYcum;
      EdbAffine2D aff2to1TXTYcum;
      float dz01corr=0;
      float dz01set=p2->Z()-p1->Z();
      
      for(int iter=0; iter<3; iter++) {
        printf("***** iter %d\n",iter);
        
        float dz01 = dz01set;
        if(doZcorr) 
        {
          dz01 = al.FineCorrZ( al.eS[0], al.eS[1] );
          printf("dz(0->1) = %f\n",dz01);
        }
    
        p1->ProjectTo( dz01 );
        EdbAffine2D aff1to2XY;
        al.CalculateAffXY( al.eS[0], al.eS[1], aff1to2XY);
        p1->Transform(&aff1to2XY);
        aff1to2XYcum.Transform(&aff1to2XY);
        aff1to2XY.Print();
       
        EdbAffine2D aff2to1TXTY;
        al.CalculateAffTXTY( al.eS[1], al.eS[0], aff2to1TXTY);
        p2->TransformA(&aff2to1TXTY);
        aff2to1TXTY.Print();
        aff2to1TXTYcum.Transform(&aff2to1TXTY);
    
        p1->ProjectTo( -dz01 );

        float dz10 = al.FineCorrZ( al.eS[1], al.eS[0] );
        printf("dz(1->0) = %f\n",dz10);
        
        dz01corr = dz01-dz01set;
        printf("dzcorr = %.1f\n",dz01corr);
      }
    
      EdbPlateP *plate2 = ss->GetPlate( p2->ScanID().ePlate );
      plate2->GetAffineTXTY()->Transform(&aff2to1TXTYcum);

      if(go_forward) {
        for(int ipatu=0; ipatu<ipat; ipatu++)
        {
          EdbPlateP *plate = ss->GetPlate( ali.GetPattern(ipatu)->ScanID().ePlate );
          plate->GetAffineXY()->Transform(&aff1to2XYcum);
          ali.GetPattern(ipatu)->SetZ( ali.GetPattern(ipatu)->Z() - dz01corr);
          plate->SetZlayer( plate->Z() - dz01corr, plate->Zmin(), plate->Zmax() );
        }
      }
      else 
      {
        for(int ipatu=npat-1; ipatu>=ipat; ipatu--)
        {
          EdbPlateP *plate = ss->GetPlate( ali.GetPattern(ipatu)->ScanID().ePlate );
          plate->GetAffineXY()->Transform(&aff1to2XYcum);
        }
      }
    
    }
  }
  if(do_corraff) if(ss) sproc.WriteScanSet(idset,*ss);
}

//---------------------------------------------------------------------------
void GlobalDiff(EdbPVRec &ali, const char *suff="")
{
  int   ntr  = ali.Ntracks();
  int   npat = ali.Npatterns();
  float xmin = ali.Xmin();
  float xmax = ali.Xmax();
  float ymin = ali.Ymin();
  float ymax = ali.Ymax();
  
  int pmin = 512, pmax=0;
  for(int i=0; i<npat; i++) {
    int p = ali.GetPattern(i)->ScanID().ePlate;
    if(p<pmin) pmin=p;
    if(p>pmax) pmax=p;
  }

  TH2F hDXR("dxr","dxr vs x",    (pmax-pmin+3)*100, pmin-1, pmax+2, 200, -200, 200 );
  TH2F hDYR("dyr","dyr vs y",    (pmax-pmin+3)*100, pmin-1, pmax+2, 200, -200, 200 );
  TH2F hDTXR("dtxr","dtxr vs x", (pmax-pmin+3)*100, pmin-1, pmax+2, 140, -0.07, 0.07 );
  TH2F hDTYR("dtyr","dtyr vs y", (pmax-pmin+3)*100, pmin-1, pmax+2, 140, -0.07, 0.07 );
  
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = ali.GetTrack(i);
    int nseg = t->N();
    for(int j=0; j<nseg; j++) {
      EdbSegP *s  = t->GetSegment(j);
      EdbSegP *sf = t->GetSegmentF(j);
      
      int plateid = s->ScanID().ePlate;
      hDXR.Fill(  (s->X()-xmin)/(xmax-xmin+1000) + plateid, sf->X()-s->X() );
      hDYR.Fill(  (s->Y()-ymin)/(ymax-ymin+1000) + plateid, sf->Y()-s->Y() );
      hDTXR.Fill( (s->X()-xmin)/(xmax-xmin+1000) + plateid, sf->TX()-s->TX() );
      hDTYR.Fill( (s->Y()-ymin)/(ymax-ymin+1000) + plateid, sf->TY()-s->TY() );
    }
  }
  
  TString name;
  sproc.MakeFileName(name,idset,"trk.an.root",false);
  Log(2,"MakeCorrectionMap","%s",name.Data());
  TFile f(name,"UPDATE");
  
  TCanvas *c = new TCanvas(Form("GDiff_%s",suff), Form("global_diff_%s",suff), 800,600);
  c->Divide(2,2);
  c->cd(1);   hDXR.Draw(); 
  c->cd(2);   hDYR.Draw();
  c->cd(3);   hDTXR.Draw();
  c->cd(4);   hDTYR.Draw();
  c->Write();
  f.Close();
}

//---------------------------------------------------------------------------
void GlobalEff(EdbPVRec &ali, const char *suff="")
{
  gROOT->SetStyle("Pub");
  
  int   ntr  = ali.Ntracks();
  int   npat = ali.Npatterns();
  int nseg0 = npat>7 ? 7 : npat;
  
  int pmin = 57, pmax=0;
  for(int i=0; i<npat; i++) {
    int p = ali.GetPattern(i)->ScanID().ePlate;
    if(p<pmin) pmin=p;
    if(p>pmax) pmax=p;
  }

  TH1F hAll("all","StartEnd vs plateID",    (pmax-pmin+3), pmin-1, pmax+2 );
  TH1F hStart("start","start",    (pmax-pmin+3), pmin-1, pmax+2 );
  TH1F hEnd("end","end",    (pmax-pmin+3), pmin-1, pmax+2 );
  
  TH1F hNSeg("nseg","nseg vs plateID",    (pmax-pmin+3), pmin-1, pmax+2 );
  
  TH1F hHoles("holes","holes vs plateID", (pmax-pmin+3), pmin-1, pmax+2 );
  
  TH1F hFillAbs("FillAbs","FillAbs", (pmax-pmin+3), pmin-1, pmax+2 );
  TH1F hFillReal("FillReal","plates fill factor", (pmax-pmin+3), pmin-1, pmax+2 );
  
  TH2F hEff("eff","eff vs theta", 50, 0, 2,  55, 0, 1.1 );
  TH2F hChi( "Chi","Chi vs plateID",   (pmax-pmin+3), pmin-1, pmax+2, 100, 0,3 );
  TH2F hW(     "W",  "W vs plateID",   (pmax-pmin+3), pmin-1, pmax+2, 40, 0,40 );
  
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = ali.GetTrack(i);
    int nseg = t->N();
    for(int j=0; j<nseg; j++) {
      EdbSegP *s  = t->GetSegment(j);
      int plateid = s->ScanID().ePlate;
      
      hNSeg.Fill(plateid);
      hChi.Fill(plateid,s->Chi2());
      hW.Fill(plateid,s->W());
      
      if(j>0&&j<nseg-1) hFillReal.Fill(plateid);
    
    }
    int pstart = t->GetSegmentFirst()->ScanID().ePlate;
    int pend   = t->GetSegmentLast()->ScanID().ePlate;
    hStart.Fill(pstart);
    hEnd.Fill(pend);
    
    for(int i=pstart; i<=pend; i++) hHoles.Fill(i);
    for(int i=pstart+1; i<pend; i++) hFillAbs.Fill(i);
    
    if(nseg>=nseg0) { hEff.Fill( t->Theta(), 1.*(t->N()-2)/(t->Npl()-2) ); }
  }
  
  hAll.Add( &hStart );  hAll.Add( &hEnd );
  hHoles.Add( &hNSeg, -1);
  
  TH1F *ratio = (TH1F*)(hFillReal.Clone("ratio"));
  ratio->Divide(&hFillAbs);
  
  TString name;
  sproc.MakeFileName(name,idset,"trk.an.root",false);
  TFile f(name,"UPDATE");
  
  TCanvas *c = new TCanvas(Form("GEff_%s",suff), Form("global_eff_%s",suff), 900,900);
  c->Divide(3,3);
  c->cd(1);   hNSeg.Draw();
  c->cd(4);   hHoles.Draw();
  c->cd(7);   { hAll.SetLineColor(kBlack); hAll.Draw(); hStart.SetLineColor(kRed); hStart.Draw("same"); hEnd.SetLineColor(kBlue); hEnd.Draw("same"); }
  c->cd(2);   ratio->Draw();
  c->cd(5);   hChi.ProfileX()->Draw();
  c->cd(8);   hW.ProfileX()->Draw();
  c->cd(3);   hEff.ProfileX()->Draw();
  c->cd(6);   hEff.ProjectionY()->Draw();
  c->Write();
  
  TCanvas *cc = new TCanvas(Form("GTr_%s",suff), Form("global_tr_%s",suff), 900,900);
  cc->Divide(1,3);
  cc->cd(1)->SetGrid();   hNSeg.Draw();
  cc->cd(2)->SetGrid();   { hAll.SetLineColor(kBlack); hAll.Draw(); hStart.SetLineColor(kRed); hStart.Draw("same"); hEnd.SetLineColor(kBlue); hEnd.Draw("same"); }
  cc->cd(3)->SetGrid();   ratio->Draw();
  cc->Write();
  
  
  f.Close();
}


int ReadFunParameters(const char* line, TF1& fun){
  const int Npars=10;
  
  double p[Npars];
  int NparsRead=0;
  std::string lineString(line);
  size_t pos=0;
  
  for (int i = 0; i < Npars; ++i){
    p[i]=0;
    if(sscanf(line+pos,"%lg",p+i)<=0)break;
    NparsRead++;
    // printf("Read parameter#%d=%g\n",i,p[i]);
    pos=lineString.find_first_of(' ',pos+1);
  }
  fun.SetParameters(p);
  return NparsRead;
}
//---------------------------------------------------------------------------
void CheckMom(EdbPVRec &ali, TEnv &cenv)
{
  int   ntr  = ali.Ntracks();
  Log(2,"Check momentum","for %d tracks",ntr);
  
  EdbMomentumEstimator mes;

  //get fitter parameters from TEnv
  TString line;
  line=cenv.GetValue("trackan.MomEst.DT","0.0021 0.0054 0");
  if(ReadFunParameters(line.Data(),mes.eDTsErrorFun)<=0) throw Form("Wrong line for \"trackan.MomEst.DT\":\'%s\'",line.Data());
  line=cenv.GetValue("trackan.MomEst.DTx","0.0021 0.0093 0");
  if(ReadFunParameters(line.Data(),mes.eDTxErrorFun)<=0) throw Form("Wrong line for \"trackan.MomEst.DTx\":\'%s\'",line.Data());
  line=cenv.GetValue("trackan.MomEst.DTy","0.0021 0 0");
  if(ReadFunParameters(line.Data(),mes.eDTyErrorFun)<=0) throw Form("Wrong line for \"trackan.MomEst.DTy\":\'%s\'",line.Data());
  mes.eAlg=cenv.GetValue("trackan.MomEst.Alg",0);
  mes.eX0=cenv.GetValue("trackan.MomEst.X0",5600);
  mes.eM=cenv.GetValue("trackan.MomEst.M",0.13957);
  mes.Print();
  //mes.eDTx0=0.001;  mes.eDTx1=0; mes.eDTx2=0;
  //mes.eDTy0=0.001;  mes.eDTy1=0; mes.eDTy2=0;

  TObjArray fitted_tracks;
  
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = ali.GetTrack(i);
    //int nseg = t->N();
//    Float_t P = mes.PMSang(*t);
    Float_t P = mes.PMS(*t);
    t->SetP(P);
    fitted_tracks.Add(t);
  }
  
  EdbDataProc::MakeTracksTree( fitted_tracks, 0., 0., Form("b%s.trk.mom.root", idset.AsString()) );
}


//---------------------------------------------------------------------------
void CheckGap(EdbPVRec &ali, const char *suff="")
{
  int   ntr  = ali.Ntracks();
  Log(2,"CheckGap","for %d tracks",ntr);
  //TH1F hGapN("gapN","GapN",   20, 0, 20 );
  //TH1F hGapZ("gapZ","GapZ",   200, 0, 20000 );
  //TH2F hDXfwd("DXfwd","DXfwd vs gap",  20, 0,20, 100, -200.,200. );
  //TH2F hDYfwd("DYfwd","DYfwd vs gap",  20, 0,20, 100, -200.,200. );
  //TH2F hDXbwd("DXbwd","DXbwd vs gap",  20, 0,20, 100, -200.,200. );
  //TH2F hDYbwd("DYbwd","DYbwd vs gap",  20, 0,20, 100, -200.,200. );
  
  TString name;
  sproc.MakeFileName(name,idset,"trk.an.root",false);
  TFile f(name,"UPDATE");
  
  EdbCouplesTree cpt;
  cpt.InitCouplesTree("gaps", 0, "WRITE");
  
  int gapcnt=0;
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = ali.GetTrack(i);
    int nseg = t->N();
    for(int j=0; j<nseg-1; j++) {
      EdbSegP *s  = t->GetSegment(j);
      EdbSegP *sn  = t->GetSegment(j+1);
      //int   ngap    = sn->ScanID().ePlate - s->ScanID().ePlate;
      //float dz      = sn->Z() - s->Z();
      //float dxfwd   = sn->X() - (s->X() + s->TX()*dz);
      //float dyfwd   = sn->Y() - (s->Y() + s->TY()*dz);
      //float dxbwd   = s->X() - (sn->X() - sn->TX()*dz);
      //float dybwd   = s->Y() - (sn->Y() - sn->TY()*dz);
      
      cpt.Fill(s,sn,t);
      gapcnt++;
      //hGapN.Fill(Abs(ngap));
      //hGapZ.Fill(Abs(dz));
      //hDXfwd.Fill( Abs(ngap) , dxfwd );
      //hDYfwd.Fill( Abs(ngap) , dyfwd );
      //hDXbwd.Fill( Abs(ngap) , dxbwd );
      //hDYbwd.Fill( Abs(ngap) , dybwd );
    }
  }
  /*
  TCanvas *c = new TCanvas(Form("Gap_%s",suff), Form("Gap_%s",suff), 900,900);
  c->Divide(2,3);
  c->cd(1)->SetLogy();   hGapN.Draw();
  c->cd(2)->SetLogy();   hGapZ.Draw();
  c->cd(3)->SetLogz();   hDXfwd.Draw("colz");
  c->cd(4)->SetLogz();   hDYfwd.Draw("colz");
  c->cd(5)->SetLogz();   hDXbwd.Draw("colz");
  c->cd(6)->SetLogz();   hDYbwd.Draw("colz");
  
  c->Write();
  */
  cpt.Close();
  Log(2,"CheckGap"," %d gaps saved",gapcnt);
  //f.Close();
}
