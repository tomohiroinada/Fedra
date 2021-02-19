//-- Author :  Luca Scotto & Valeri Tioukov   15.08.2007

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbRunTracking                                                       //
//                                                                      //
// Prediction search in the raw data                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "EdbRunTracking.h"
#include "EdbBrick.h"
#include "EdbLog.h"
#include "EdbTrackFitter.h"

using namespace TMath;

ClassImp(EdbRunTracking)

///_________________________________________________________________________
EdbRunTracking::~EdbRunTracking()
{

}

//_________________________________________________________________________
void EdbRunTracking::Set0()
{
  eCondMT.SetDefault();
  eCondBT.SetDefault();

  eDeltaRview           = 400;
  eDeltaTheta           = 0.15;
  ePreliminaryPulsMinMT = 5;
  ePreliminaryChi2MaxMT = 5;

  eDeltaR      = 20;
  ePulsMinBT   = 18;
  ePulsMinDegradBT = 0;
  eChi2MaxBT   = 2.5;

  eChi2MaxMT   = 1.6;
  ePulsMinDegradMT = 0;
  ePulsMinMT   = 10.;

  eDegradPos   = 0;
  eDegradSlope = 0;
  ePredictionScan = false;
}

//----------------------------------------------------------------------------------------
void EdbRunTracking::SetPred(const EdbSegP &pred)
{
  ePred.Copy(pred);
  eNext.Set0();
  eS.Set0(); eS1.Set0(); eS2.Set0();
  eS1cnd.Reset();
  eS2cnd.Reset();
  eScnd.Reset();
  eS1pre.Reset();
  eS2pre.Reset();
  eSpre.Reset();
}

//----------------------------------------------------------------------------------------
void EdbRunTracking::Print()
{
  printf("EdbRunTracking selection criteria:\n"); 
  printf("Preliminary mt selection:\n");
  printf("                         eDeltaRview           = %5.3f\n", eDeltaRview);
  printf("                         eDeltaTheta           = %5.3f\n", eDeltaTheta);
  printf("                         ePreliminaryPulsMinMT = %3.1f\n", ePreliminaryPulsMinMT);
  printf("                         ePreliminaryChi2MaxMT = %5.3f\n", ePreliminaryChi2MaxMT);
  printf("Preliminary bt selection:\n");
  printf("                         eDeltaR               = %3.1f\n", eDeltaR);
  printf("                         eChi2MaxBT            = %5.3f\n", eChi2MaxBT);
  printf("Final       bt selection:\n");
  printf("                         ePulsMinBT            = %3.1f\n", ePulsMinBT);
  printf("                         ePulsMinDegradBT      = %3.1f\n", ePulsMinDegradBT);
  printf("Final       mt selection:\n");
  printf("                         ePulsMinMT            = %3.1f\n", ePulsMinMT);
  printf("                         ePulsMinDegradMT      = %3.1f\n", ePulsMinDegradMT);
  printf("                         eChi2MaxMT            = %5.3f\n", eChi2MaxMT);
  printf("\n");
}

//----------------------------------------------------------------------------------------
int EdbRunTracking::UpdateFlag(int flag, int status)
{
  // status: -1 -found nothing, 0-bt, 1-mt1, 2-mt2

  int bth = flag/10000;
  int mth = (flag/100)%100;
  int tb  = flag%10;

  switch (status) {
  case -1: bth++; mth++;       break;
  case  0: bth=0; mth=0; tb=0; break;
  case  1: bth++; mth=0; tb=1; break;
  case  2: bth++; mth=0; tb=2; break;
  }

  return( bth*10000+mth*100+tb );
}

//----------------------------------------------------------------------------------------
int EdbRunTracking::ExtrapolateCond(EdbScanCond &inputcond, int flag, EdbScanCond &outputcond )
{
  // TODO: tuning the dependency of sigma
  // NOTED by Artem: when we do jumping we have no "holes" stored in the flag, so this function do not 
  //                 extrapolate errors correctly so one had to encrease the acceptance manually, 
  //                 possible solution can be to add total plates number before the last found 
  //                 track into the flag and use it only for extrapolation

  int bth = GetBTHoles(flag);
  int mth = GetMTHoles(flag);
  outputcond = inputcond;
  outputcond.SetSigma0( inputcond.SigmaX(0)  + eDegradPos   * mth,
			inputcond.SigmaY(0)  + eDegradPos   * mth,
			inputcond.SigmaTX(0) + eDegradSlope * bth,
			inputcond.SigmaTY(0) + eDegradSlope * bth );
  return(0);
}

//----------------------------------------------------------------------------------------
int EdbRunTracking::FindCompliments( EdbSegP &s, EdbPattern &pat, TObjArray &found, float chi2max, TArrayF &chiarr )
{
  // return found sorted by increasing chi2

  int nfound=0;
  int maxcand=chiarr.GetSize();
  TArrayF   chi2arr(maxcand);
  TObjArray arr(maxcand);
  TArrayI   ind(maxcand);

  int nseg = pat.FindCompliments(s,arr,30,200);  // acceptance (prelim): s.SX()*30; s.STX*200
  //  printf("\nnseg = %d\n",nseg);
  if(nseg>maxcand)  {
    printf("Warning!: Too many segments %d, accept only the first %d \n", nseg, maxcand);
    nseg = maxcand;
  }
  if(nseg<=0) return 0;

  EdbSegP *s2=0;
  for(int j=0; j<nseg; j++) {
    s2 = (EdbSegP *)arr.At(j);
    EdbSegP s3;
    s3.Copy(s);
    chi2arr[j] = EdbTrackFitter::Chi2Seg(&s3, s2);
    //chi2arr[j] = EdbTrackFitter::Chi2SegM(s, *s2,s3,eCondBT,eCondMT);
  }
  TMath::Sort(nseg,chi2arr.GetArray(),ind.GetArray(),0);
  //  printf("pred = %f %f %f %f\n",s.X(),s.Y(),s.TX(),s.TY());
  for(int j=0; j<nseg; j++) {
    s2 = (EdbSegP *)arr.At(ind[j]);
    //    printf("j = %d, ind = %d, chi2 = %f    %f\n",j,ind[j],chi2arr[ind[j]],chi2max);
    //    printf("Tx = %f %f %f %f\n",s2->X(),s2->Y(),s2->TX(),s2->TY());
    if(chi2arr[ind[j]] > chi2max ) break;
    chiarr[j] = chi2arr[ind[j]];
    s2->SetMC(s.MCEvt(),s.MCTrack());
    found.Add(s2);
    nfound++;
  }

  //  printf("nfound = %d\n",nfound);
  return nfound;
}


//----------------------------------------------------------------------------------------
int EdbRunTracking::FindCandidates( EdbSegP &spred, EdbPattern &fndbt, EdbPattern &fnds1, EdbPattern &fnds2 )
{
  // Find microtracks and basetracks for the prediction segment "spred"
  // Selection criteria:
  //  1) select all microtracks with (puls >= ePreliminaryPulsMinMT) (6) and (chi2 < ePreliminaryChi2MaxMT) for both side
  //  2) select all basetracks by using the selected microtracks, with (chi2_bt < eChi2MaxBT)
  // Microtracks and basetracks in output are sorted according to chi2 (ascending order = the first one is the best one)
  //
  // Input:
  //   spred - track prediction
  // Output: 
  //   fnds1 - microtracks having puls>=ePreliminaryPulsMinMT and chi2<ePreliminaryChi2MaxMT 
  //           in the top side of the emulsion plate
  //   fnds2 - microtracks having puls>=ePreliminaryPulsMinMT and chi2<ePreliminaryChi2MaxMT 
  //           in the bottom side of the emulsion plate
  //   fndbt - basetracks built from the microtracks and (chi2_bt < eChi2MaxBT)

  ClearCuts();

  EdbScanCond condBT;
  ExtrapolateCond(eCondBT,spred.Flag(),condBT);
  if (gEDBDEBUGLEVEL>=3) condBT.Print();

  spred.SetZ(107.);                                      // TODO!
  spred.SetErrors();
  condBT.FillErrorsCov( spred.TX(), spred.TY(), spred.COV() );

  //EdbPVRec aview;  //with 2 patterns of preselected microtracks
  //aview.AddPattern( new EdbPattern(0,0,214));  // TODO! sequence??
  //aview.AddPattern( new EdbPattern(0,0,0)  );

  float xmin[5]={-500, -500, spred.TX()-eDeltaTheta, spred.TY()-eDeltaTheta, ePreliminaryPulsMinMT };         //TODO!!
  float xmax[5]={ 500,  500, spred.TX()+eDeltaTheta, spred.TY()+eDeltaTheta, 50 };

  for(int side=1; side<=2; side++) {
    EdbPattern pat;
    AddSegmentCut(side,1,xmin,xmax);
    if(ePredictionScan) GetPatternDataForPrediction( spred.ID(), side, pat );
    else                GetPatternXY( spred, side,  pat, eDeltaRview);
    Log(2,"FindCandidates","%d microtracks in side %d",pat.N(),side);

    for(int i=0; i<pat.N(); i++) {
      EdbSegP *s = pat.GetSegment(i);
      s->SetErrors();
      eCondMT.FillErrorsCov( spred.TX(), spred.TY(), s->COV() );
    }

    pat.FillCell(10,10,0.01,0.01);  //divide view on this cells

    TArrayF chi2arr(10000);  //TODO!
    TObjArray found;
    FindCompliments( spred, pat, found, ePreliminaryChi2MaxMT, chi2arr );

    for(int j=0; j<found.GetEntries(); j++) {
      EdbSegP *s = (EdbSegP *)(found.At(j));
      s->SetChi2(chi2arr[j]);
      if     (side==1) fnds1.AddSegment(*s);
      else if(side==2) fnds2.AddSegment(*s);
    }

  }

  // filling fndbt
  
  EdbPattern bt;

  for(int is1=0; is1<fnds1.N(); is1++) {
    for(int is2=0; is2<fnds2.N(); is2++) {

      EdbSegP *s1 = fnds1.GetSegment(is1);
      EdbSegP *s2 = fnds2.GetSegment(is2);

      float dx1=s1->X()-(spred.X()+spred.TX()*(s1->Z()-spred.Z()));
      float dy1=s1->Y()-(spred.Y()+spred.TY()*(s1->Z()-spred.Z()));
      float dx2=s2->X()-(spred.X()+spred.TX()*(s2->Z()-spred.Z()));
      float dy2=s2->Y()-(spred.Y()+spred.TY()*(s2->Z()-spred.Z()));
      float r = Sqrt( (dx1-dx2)*(dx1-dx2) + (dy1-dy2)*(dy1-dy2) ); 

      if(r<eDeltaR) {  // has good BT
	EdbSegP s3;
	s3.Copy(spred);
	s3.SetX( 0.5*(s1->X() + s2->X()) );
	s3.SetY( 0.5*(s1->Y() + s2->Y()) );
	s3.SetZ( 0.5*(s1->Z() + s2->Z()) );
	s3.SetTX( (s2->X() - s1->X()) / (s2->Z() - s1->Z()) );
	s3.SetTY( (s2->Y() - s1->Y()) / (s2->Z() - s1->Z()) );
	s3.SetW(s1->W()+s2->W());

	s3.SetFlag(is2*10000+is1);

	EdbSegP s4(spred);
	float chi = EdbTrackFitter::Chi2Seg(&s4, &s3);  // depends on the sequence (best defined should be first)
	//float chi = EdbTrackFitter::Chi2SegM(spred, s3, s4, eCondBT, eCondBT);
	if(chi<eChi2MaxBT) {
	  s3.SetChi2(chi);
	  bt.AddSegment(s3);
	}
      }
    }
  }

  TArrayF chi2arr(bt.N());
  for (int i=0;i<bt.N();i++) chi2arr[i]=bt.GetSegment(i)->Chi2();
  TArrayI   ind(bt.N());
  TMath::Sort(bt.N(),chi2arr.GetArray(),ind.GetArray(),false);
  for(int i=0; i<bt.N(); i++) {
    EdbSegP *s = bt.GetSegment(ind[i]);
    fndbt.AddSegment(*s);
  }

  Log(2,"FindCandidates","Found %d basetrack candidate in %d+%d preselected microtracks",fndbt.N(),fnds1.N(),fnds2.N());

  return 1; //TODO!
}


//----------------------------------------------------------------------------------------
int EdbRunTracking::FindCandidateMT( EdbPattern &fnds1, EdbPattern &fnds2, EdbSegP &fnd )
{
  EdbSegP s1,s2;
  int n1=FindBestCandidate( fnds1, s1, eS1cnd, ePulsMinMT, ePulsMinDegradMT, eChi2MaxMT);
  int n2=FindBestCandidate( fnds2, s2, eS2cnd, ePulsMinMT, ePulsMinDegradMT, eChi2MaxMT);
  Log(2,"FindCandidateMT","Found %d+%d microtrack candidates after cuts",n1,n2);

  if( n1==0&&n2==0 )  return 0;

  if( s1.Chi2() <= s2.Chi2() ) {
    fnd.Copy(s1);
    return 1;
  }
  else {
    fnd.Copy(s2);
    return 2;
  }
  return 0;
}

//----------------------------------------------------------------------------------------
int EdbRunTracking::FindBestCandidate(EdbPattern &cand, EdbSegP &fnd, EdbPattern &passed, float wmin, float wmindegrad, float chi2max )
{
  int n=0;
  fnd.Set0();
  fnd.SetChi2(10000.+chi2max);
  for (int i=0; i<cand.N(); i++) {
    EdbSegP *s = cand.GetSegment(i);
    if ( s->W()<wmin+s->Chi2()*wmindegrad )  continue;
    if ( s->Chi2()>chi2max )                 continue;
    n++;
    passed.AddSegment(*s);
    if (s->Chi2()<fnd.Chi2())      fnd.Copy(*s);
  }
  return n;
}

//----------------------------------------------------------------------------------------
int EdbRunTracking::FindPrediction( EdbSegP &spred, EdbSegP &fndbt, EdbSegP &fnds1, EdbSegP &fnds2, EdbSegP &snewpred )
{
  // Select the best (micro or base) track matching with the prediction
  // and prepare for a new search.
  //
  // Selection criteria:
  //  1) Call FindCandidates having the list of basetrack and microtrack candidates
  //  2) Call FindCandidateBT which looks for the best basetrack, if any
  //  3) If no basetrack is found, call FindCandidateMT which looks fot the best microtrack, if any.
  //     Microtracks accepted shold satisfy the following cut: (puls >= ePulsMinMT) (10) and (chi2 < eChi2MaxMT)
  //
  // Input:
  //   spred - track prediction
  // Output:
  //   - if a basetrack is found (status 0):
  //       fndbt - basetrack found
  //       fnds1 - microtrack top contained in the found basetrack
  //       fnds2 - microtrack bottom contained in the found basetrack
  //       snewpred - fndbt with flag equal to 0
  //
  //   - if a microtrack top is found (status 1):
  //       fndbt - dummy
  //       fnds1 - microtrack top found
  //       fnds2 - dummy
  //       snewpred - a track with slopes from the prediction spred and positions 
  //                  from an extrapolation of fnds1. The flag is updated by UpdateFlag
  //
  //   - if a microtrack bottom is found (status 2):
  //       fndbt - dummy
  //       fnds1 - dummy
  //       fnds2 - microtrack bottom found
  //       snewpred - a track with slopes from the prediction spred and positions 
  //                  from an extrapolation of fnds2. The flag is updated by UpdateFlag
  //
  //   - if nothing is found (status -1):
  //       fndbt - dummy
  //       fnds1 - dummy
  //       fnds2 - dummy
  //       snewpred - the prediction spred. The flag is updated by UpdateFlag
  //
  // Return:
  //   -1: no track found
  //    0: basetrack found
  //    1: microtrack top found
  //    2: microtrack bottom found

  //EdbPattern vfndbt,vfnds1,vfnds2;

  eStatus = -1;
  SetPred(spred);
  FindCandidates( spred, eSpre, eS1pre, eS2pre );

  EdbSegP fnd;
  int nbt = FindBestCandidate(eSpre,fnd, eScnd, ePulsMinBT, ePulsMinDegradBT, eChi2MaxBT);
  if ( nbt > 0 ) {
    eS.Copy(fnd);
    eS1.Copy(*(eS1pre.GetSegment(fnd.Flag()%10000)));
    eS2.Copy(*(eS2pre.GetSegment(fnd.Flag()/10000)));
    eS.SetFlag(0);
    eNext.Copy(fnd);
    eNext.SetFlag(UpdateFlag(spred.Flag(),0));          // if bt found : bth=0, mth=0, tb=0
    eStatus = 0;
   }

  int if_mt = FindCandidateMT(eS1pre,eS2pre,fnd);
  if(eStatus!=-1) goto RESUME;

  switch(if_mt) {
  case 0:                       // find nothing
    eNext.Copy(spred);
    eNext.SetFlag(UpdateFlag(spred.Flag(),-1));      // hole: if not found: bth++, mth++, tb= keep last value
    eNext.SetW(0);
    eStatus = -1;    goto RESUME;
  case 1:         // best microtrack is on the 1-st side
    eS1.Copy(fnd);
    eNext.Copy(spred);
    eNext.SetX( fnd.X() + spred.TX()*(spred.Z()-fnd.Z()) );
    eNext.SetY( fnd.Y() + spred.TY()*(spred.Z()-fnd.Z()) );
    eNext.SetZ(spred.Z());
    eNext.SetFlag(UpdateFlag(spred.Flag(),1));          // if mt found : bth++, mth=0, tb=1
    eNext.SetW(fnd.W());
    eStatus = 1;    goto RESUME;
  case 2:         // best microtrack is on the 2-d side
    eS2.Copy(fnd);
    eNext.Copy(spred);
    eNext.SetX( fnd.X() + spred.TX()*(spred.Z()-fnd.Z()) );
    eNext.SetY( fnd.Y() + spred.TY()*(spred.Z()-fnd.Z()) );
    eNext.SetZ(spred.Z());
    eNext.SetFlag(UpdateFlag(spred.Flag(),2));         // if mt found : bth++, mth=0, tb=2
    eNext.SetW(fnd.W());
    eStatus = 2;    goto RESUME;
  }

 RESUME:
  
  Log(2,"FindPrediction","status = %d, good candidates [s:s1:s2] %d:%d:%d ;  preliminary [s:s1:s2] %d:%d:%d",
      eStatus, 
      eScnd.N(),eS1cnd.N(),eS2cnd.N(),
      eSpre.N(),eS1pre.N(),eS2pre.N()
      );
  snewpred.Copy(eNext);
  fndbt.Copy(eS);
  fnds1.Copy(eS1);
  fnds2.Copy(eS2);
  return eStatus;
}


//----------------------------------------------------------------------------------------
TTree *EdbRunTracking::InitSBtree(const char *file_name, const char *mode)
{
  const char *tree_name="sbt";
  TTree *tree=0;
  if (!tree) {
    TFile *f = new TFile(file_name,mode);
    if (f)  tree = (TTree*)f->Get(tree_name);
    if(!tree) {

      f->cd();
      tree = new TTree(tree_name,tree_name);
      tree->SetMaxTreeSize(15000000000LL);   //set 15 Gb file size limit)

      Int_t idp[4], idf[4], stat;
      EdbSegP *s_pred = 0, *s_bt = 0, *s_mt1=0, *s_mt2=0, *s_next=0;
      TClonesArray *scnd   = new TClonesArray("EdbSegP");
      TClonesArray *s1cnd  = new TClonesArray("EdbSegP");
      TClonesArray *s2cnd  = new TClonesArray("EdbSegP");
      TClonesArray *spre   = new TClonesArray("EdbSegP");
      TClonesArray *s1pre  = new TClonesArray("EdbSegP");
      TClonesArray *s2pre  = new TClonesArray("EdbSegP");
     
      tree->Branch("idpred",idp,"idp[4]/I");
      tree->Branch("idfound",idf,"idf[4]/I");
      tree->Branch("stat",&stat,"stat/I");
      tree->Branch("pred.", "EdbSegP", &s_pred);
      tree->Branch("s.", "EdbSegP", &s_bt);
      tree->Branch("s1.", "EdbSegP", &s_mt1);
      tree->Branch("s2.", "EdbSegP", &s_mt2);
      tree->Branch("next.", "EdbSegP", &s_next);
      tree->Branch("scnd",&scnd);
      tree->Branch("s1cnd",&s1cnd);
      tree->Branch("s2cnd",&s2cnd);
      tree->Branch("spre",&spre);
      tree->Branch("s1pre",&s1pre);
      tree->Branch("s2pre",&s2pre);
      tree->Write();
    }
  }

  if(!tree) Log(1,"InitSBtree","ERROR!!! can't initialize tree at %s as %s\n",file_name,mode);
  return tree;
}

//----------------------------------------------------------------------------------------
bool  EdbRunTracking::UpdateSBtree( TTree &tsbt, int idp[4], int idf[4])
  // int stat, 
  //			    EdbSegP &ps, EdbSegP &fndbt, EdbSegP &fnds1, EdbSegP &fnds2, EdbSegP &nextpred, 
  //			    EdbPattern *scnd, EdbPattern *s1cnd, EdbPattern *s2cnd )
{
  EdbSegP *s_pred = &ePred, *s_bt = &eS, *s_mt1=&eS1, *s_mt2=&eS2, *s_next=&eNext;

  tsbt.SetBranchAddress("idpred",idp);
  tsbt.SetBranchAddress("idfound",idf);
  tsbt.SetBranchAddress("stat",&eStatus);
  tsbt.SetBranchAddress("pred.", &s_pred);
  tsbt.SetBranchAddress("s.", &s_bt);
  tsbt.SetBranchAddress("s1.", &s_mt1);
  tsbt.SetBranchAddress("s2.", &s_mt2);
  tsbt.SetBranchAddress("next.", &s_next);

  TClonesArray *s_cnd = eScnd.GetSegments();
  tsbt.SetBranchAddress("scnd", &s_cnd);
  TClonesArray *s1_cnd = eS1cnd.GetSegments();
  tsbt.SetBranchAddress("s1cnd", &s1_cnd);
  TClonesArray *s2_cnd = eS2cnd.GetSegments();
  tsbt.SetBranchAddress("s2cnd", &s2_cnd);

  TClonesArray *s_pre = eSpre.GetSegments();
  tsbt.SetBranchAddress("spre", &s_pre);
  TClonesArray *s1_pre = eS1pre.GetSegments();
  tsbt.SetBranchAddress("s1pre", &s1_pre);
  TClonesArray *s2_pre = eS2pre.GetSegments();
  tsbt.SetBranchAddress("s2pre", &s2_pre);

  tsbt.Fill();
  return true;
}


//----------------------------------------------------------------------------------------
bool  EdbRunTracking::GetSBtreeEntry( int entry, TTree &tsbt)
{
  EdbSegP *s_pred = &ePred, *s_bt = &eS, *s_mt1=&eS1, *s_mt2=&eS2, *s_next=&eNext;

  tsbt.SetBranchAddress("idpred",eIdp);
  tsbt.SetBranchAddress("idfound",eIdf);
  tsbt.SetBranchAddress("stat",&eStatus);
  tsbt.SetBranchAddress("pred.", &s_pred);
  tsbt.SetBranchAddress("s.", &s_bt);
  tsbt.SetBranchAddress("s1.", &s_mt1);
  tsbt.SetBranchAddress("s2.", &s_mt2);
  tsbt.SetBranchAddress("next.", &s_next);

  TClonesArray *s_cnd = eScnd.GetSegments();
  tsbt.SetBranchAddress("scnd", &s_cnd);
  TClonesArray *s1_cnd = eS1cnd.GetSegments();
  tsbt.SetBranchAddress("s1cnd", &s1_cnd);
  TClonesArray *s2_cnd = eS2cnd.GetSegments();
  tsbt.SetBranchAddress("s2cnd", &s2_cnd);

  TClonesArray *s_pre = eSpre.GetSegments();
  tsbt.SetBranchAddress("spre", &s_pre);
  TClonesArray *s1_pre = eS1pre.GetSegments();
  tsbt.SetBranchAddress("s1pre", &s1_pre);
  TClonesArray *s2_pre = eS2pre.GetSegments();
  tsbt.SetBranchAddress("s2pre", &s2_pre);

  tsbt.GetEntry(entry);
  return true;
}

//______________________________________________________________________________
void EdbRunTracking::CloseSBtree(TTree *tree)
{
  tree->AutoSave();
  TFile *f=0;
  f = tree->GetCurrentFile();
  if(f) {
    f->Purge();
    f->Close();
  }
  tree=0;
}

/*
//______________________________________________________________________________
float EdbRunTracking::TrackExtrapolationToZ(EdbTrackP &t, float z, EdbSegP &ps)
{
  if     (z<=t.Zmin()) ps.Copy(*t.());
  else if(z>=t.Zmax()) ps.Copy(*t.GetSegmentFLast());
  else         // z is inside the track
    { Log(2, "EdbRunTracking::TrackExtrapolationToZ", " TODO!\n"); }
  float dz = Abs(ps.Z() - z);
  ps.PropagateTo( z );
  return dz;
}
*/

//______________________________________________________________________________
int EdbRunTracking::FindTrack(EdbTrackP &pred, EdbTrackP &found, EdbPlateP &plate)
{
  // look for tracks in this plate
  // track - input track in brick RS - will be updated on output
  // plate  - all plate parameters including affine transformation plate-to-brick

  int status=-100;
  float DZmax = 1350*100;

  EdbAffine2D p2b(*(plate.GetAffineXY()));   // from plate to brick
  EdbAffine2D b2p(p2b); b2p.Invert();        // from brick to plate
  EdbSegP ps;
  
  float dz = pred.MakePredictionTo(plate.Z(), ps);
  printf("EdbRunTracking::FindTrack: dz: %f - %f = %f\n", plate.Z(), pred.Z(), dz );
  
  if(Abs(dz)>DZmax)               return status;
  if(GetBTHoles(pred.Flag())>5)   return status;
  if(GetMTHoles(pred.Flag())>3)   return status;

  ps.SetFlag(pred.Flag());
  ps.PrintNice();
  ps.Transform(&b2p);                     // plate.Transoform(seg) ???
  ps.PrintNice();
 
  EdbSegP fndbt, fnds1, fnds2, snewpred;
  status = FindPrediction( ps, fndbt, fnds1, fnds2, snewpred ); // -1: not found; 0-bt, 1-bot, 2-top
  found.SetFlag(snewpred.Flag());
  fnds1.PrintNice(); 
  fnds2.PrintNice();

  TransformFromPlateRS(plate);      // transform all components into brick RS

  if(status>=0) {
    found.AddSegment(  new EdbSegP(eNext) );
    found.AddSegmentF( new EdbSegP(ePred) );   // add prediction as "fitted segment" because it is an extrapolation
    found.SetSegmentsTrack();
  }

  Log(2,"EdbRunTracking::FindTracks","status = %d",status);
  return status;
}

//______________________________________________________________________________
void EdbRunTracking::TransformFromPlateRS(EdbPlateP &plate)
{
  EdbAffine2D p2b(*(plate.GetAffineXY()));   // from plate to brick

  ePred.Transform(&p2b);
  eNext.Transform(&p2b);
  eS.Transform(&p2b);
  eS1.Transform(&p2b);
  eS2.Transform(&p2b);
  eScnd.Transform(&p2b);
  eS1cnd.Transform(&p2b);
  eS2cnd.Transform(&p2b);
  eSpre.Transform(&p2b);
  eS1pre.Transform(&p2b);
  eS2pre.Transform(&p2b);
  
  ePred.SetPID( plate.ID() );
  eNext.SetPID( plate.ID() );
  eS.SetPID(    plate.ID() );
  eS1.SetPID(   plate.ID() );
  eS2.SetPID(   plate.ID() );
  for(int i=0; i<eScnd.N(); i++)  eScnd.GetSegment(i)->SetPID( plate.ID() );
  for(int i=0; i<eS1cnd.N(); i++) eS1cnd.GetSegment(i)->SetPID( plate.ID() );
  for(int i=0; i<eS2cnd.N(); i++) eS2cnd.GetSegment(i)->SetPID( plate.ID() );
  for(int i=0; i<eSpre.N(); i++)  eSpre.GetSegment(i)->SetPID( plate.ID() );
  for(int i=0; i<eS1pre.N(); i++) eS1pre.GetSegment(i)->SetPID( plate.ID() );
  for(int i=0; i<eS2pre.N(); i++) eS2pre.GetSegment(i)->SetPID( plate.ID() );

  ePred.SetZ( plate.Z() );
  eNext.SetZ( plate.Z() );
  eS.SetZ(    plate.Z() );
  eS1.SetZ(   plate.GetLayer(1)->Z() + plate.Z() );
  eS2.SetZ(   plate.GetLayer(2)->Z() + plate.Z() );
  for(int i=0; i<eScnd.N(); i++)  eScnd.GetSegment(i)->SetZ( plate.Z() );
  for(int i=0; i<eS1cnd.N(); i++) eS1cnd.GetSegment(i)->SetZ( plate.GetLayer(1)->Z() + plate.Z() );
  for(int i=0; i<eS2cnd.N(); i++) eS2cnd.GetSegment(i)->SetZ( plate.GetLayer(2)->Z() + plate.Z() );
  for(int i=0; i<eSpre.N(); i++)  eSpre.GetSegment(i)->SetZ( plate.Z() );
  for(int i=0; i<eS1pre.N(); i++) eS1pre.GetSegment(i)->SetZ( plate.GetLayer(1)->Z() + plate.Z() );
  for(int i=0; i<eS2pre.N(); i++) eS2pre.GetSegment(i)->SetZ( plate.GetLayer(2)->Z() + plate.Z() );
  
  eNext.SetDZ( plate.GetLayer(0)->DZ() );
  eS.SetDZ(    plate.GetLayer(0)->DZ() );
  eS1.SetDZ(   plate.GetLayer(1)->DZ());
  eS2.SetDZ(   plate.GetLayer(2)->DZ());
  for(int i=0; i<eScnd.N(); i++)  eScnd.GetSegment(i)->SetDZ( plate.GetLayer(0)->DZ() );
  for(int i=0; i<eS1cnd.N(); i++) eS1cnd.GetSegment(i)->SetDZ( plate.GetLayer(1)->DZ() );
  for(int i=0; i<eS2cnd.N(); i++) eS2cnd.GetSegment(i)->SetDZ( plate.GetLayer(2)->DZ() );
  for(int i=0; i<eSpre.N(); i++)  eSpre.GetSegment(i)->SetDZ( plate.GetLayer(0)->DZ() );
  for(int i=0; i<eS1pre.N(); i++) eS1pre.GetSegment(i)->SetDZ( plate.GetLayer(1)->DZ() );
  for(int i=0; i<eS2pre.N(); i++) eS2pre.GetSegment(i)->SetDZ( plate.GetLayer(2)->DZ() );
}

//______________________________________________________________________________
void EdbRunTracking::CheckZ( float &z1, float &z2 )
{
  EdbSegP *s1 = (EdbSegP*)(eS1pre.At(0));
  EdbSegP *s2 = (EdbSegP*)(eS2pre.At(0));
  if(s1&&s2)
  {
    z1=s1->Z();
    z2=s2->Z();
  }
  else
  {
    z1 = 214;
    z2 = 0;
  }
}

//______________________________________________________________________________
int EdbRunTracking::GetSegmentsForDB( EdbSegP &s, EdbSegP &s1, EdbSegP &s2 )
{
  // return s,s1,s2 if status >=0, create fake microtracks if necessary
  // for fakes:        seg->GetPuls(), seg->GetVolume(), seg->GetSigmaX()            //0, 0, -1
  float z1=0,z2=0;
  CheckZ(z1,z2);
  s.Copy(eNext);  //has z in the center of base (107)
  if (eStatus==0) // both microtracks and basetrack found
  {
    s1.Copy(eS1);
    s2.Copy(eS2);
  }
  else if (eStatus==1) // s1 found
  {
    s1.Copy(eS1);
    s2.Copy(eNext);
    s2.PropagateTo(z1);
    s2.SetW(0); s2.SetChi2(0);
   }
  else if (eStatus==2) // s1 found
  {
    s2.Copy(eS2);
    s1.Copy(eNext);
    s1.PropagateTo(z2);
    s1.SetW(0); s1.SetChi2(0);
  }
  else 
  {
    s1.Copy(eNext);
    s2.Copy(eNext);
  }
  Log(2,"EdbRunTracking::GetSegmentsForDB","status = %d, z1=%f  z2= %f",eStatus,s1.Z(),s2.Z());
  return eStatus;
}
