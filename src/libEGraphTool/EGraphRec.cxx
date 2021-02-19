//----------------------------------------------------------------------------
// Program: Emulsion Graphical Reconstraction Toolkit - library
//
// Author:  Artem Chukanov (chukanov at nusun.jinr.ru)
//          31.01.2008
//
//----------------------------------------------------------------------------

#include "EGraphRec.h"
#include "EdbLog.h"
#include <TEnv.h>
#include <TGTab.h>
#include <TGLabel.h>
#include <TGTextView.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGFileDialog.h>
#include <TBranchClones.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TGLSAViewer.h>
#include <TView.h>
#include <TThread.h>
#include <TMath.h>
#include <iostream>

#include "EdbScanProc.h"
#include "EdbDisplay.h"

using namespace std;
using namespace TMath;

ClassImp(EGraphRec)


//----------------------------------------------------------------------------
EGraphRec::EGraphRec()
{
  InitVariables();
  InitDrawVariables();
  ReadCmdConfig();      // reading config file
}


//----------------------------------------------------------------------------
EGraphRec::~EGraphRec() 
{
  SafeDelete(fThSBProcess);
  SafeDelete(fThSBCheckProcess);
  SafeDelete(fPredTracks);
  SafeDelete(fFoundTracks);
  SafeDelete(fRecProc);
}


//----------------------------------------------------------------------------
void EGraphRec::ProcessEvent()
{
  // DrawEvent(nentries);
  // InitScanSet();
}


//----------------------------------------------------------------------------
void EGraphRec::ResetProcess()
{
  // Reset process

//   cout << "Reset Process" << endl;
  
//   if (fThProcessEvent) {
//     TThread::Delete(fThProcessEvent);
//     SafeDelete(fThProcessEvent);
//     fThProcessEvent      = new TThread("ThProcessEvent", 
//  				       ThProcessEvent, (void*) this);
//   }
}


//----------------------------------------------------------------------------
void EGraphRec::StartScanBack()
{
  // start scan back

  // Initialization of brick to be processed

  fBrickToProc.brickId    = fEntryProcBrickId->GetIntNumber();
  fBrickToProc.ver        = fEntryProcVer->GetIntNumber();
  fBrickToProc.firstPlate = fEntrySBFirstPlate->GetIntNumber();
  fBrickToProc.lastPlate  = fEntrySBLastPlate->GetIntNumber();
  fBrickToProc.step       = fEntrySBStep->GetIntNumber();

  InitScanSet();

  SafeDelete(fFoundTracks);
  fFoundTracks = new EdbPVRec();
  fThSBProcess = new TThread("ThSBProcess", ThSBProcess, (void*) this);

  fThSBProcess->Run();                // Run scan back Process in Thread mode
  fThSBCheckProcess->Run();           // Check the end of job
  fButtonSBStart->SetEnabled(kFALSE); // Disable Execution button
}


//----------------------------------------------------------------------------
void EGraphRec::StartVertexRec()
{
  // start vertex reconstruction

  // Initialization of brick to be processed

  fBrickToProc.brickId    = fEntryProcBrickId->GetIntNumber();
  fBrickToProc.ver        = fEntryProcVer->GetIntNumber();
  fBrickToProc.firstPlate = fEntrySBFirstPlate->GetIntNumber();
  fBrickToProc.lastPlate  = fEntrySBLastPlate->GetIntNumber();
  fBrickToProc.step       = fEntrySBStep->GetIntNumber();

  fVertexRecOpt.QualityMode = fEntryQualityMode->GetIntNumber();
  fVertexRecOpt.UseMom      = fCheckUseMom->IsDown();
  fVertexRecOpt.UseSegPar   = fCheckUseSegPar->IsDown();
  fVertexRecOpt.DZmax       = fEntryDZmax->GetNumber();
  fVertexRecOpt.ProbMinV    = fEntryProbMinV->GetNumber();
  fVertexRecOpt.ImpMax      = fEntryImpMax->GetNumber();

  fScanProc->eProcDirClient = fDataDir;     // brick directory initialization

  fRecProc->SetScanProc(fScanProc);         // Scan Proc
  fRecProc->SetBrickToProc(fBrickToProc);   // Brick to process
  fRecProc->SetProcId(fProcId);             // Process Ids
  fRecProc->SetVertexRecOpt(fVertexRecOpt); // vertex reconstruction options

  fPVRec = fRecProc->VertexRec();

  DrawEvent();

  // delete recProc;
}


//----------------------------------------------------------------------------
void EGraphRec::DrawEvent()
{
  fDisplay->fCanvas->Clear("D");
  fDisplay->fCanvas->cd();

  fDisplay->SetDrawTracks(4);
  // fDisplay->SetDrawVertex(1);

  if (fFoundTracks) fDisplay->SetArrTr(fFoundTracks->eTracks);
  if (fPVRec) fDisplay->SelectVertexTracks(fPVRec->eVTX);

  fDisplay->Draw();
  fDisplay->fCanvas->Update();

}


//----------------------------------------------------------------------------
void EGraphRec::ReconstructTracks()
{

}


//----------------------------------------------------------------------------
void EGraphRec::SetTree(TTree *tree)
{
  fEvent = new EdbView();
  fEvtTree = tree;

  TBranch       *headerBranch   = fEvtTree->GetBranch("headers");
  TBranchClones *clustersBranch = (TBranchClones*)(fEvtTree->GetBranch("clusters"));
  TBranchClones *segmentsBranch = (TBranchClones*)(fEvtTree->GetBranch("segments"));
  TBranchClones *tracksBranch   = (TBranchClones*)(fEvtTree->GetBranch("tracks"));
  TBranchClones *framesBranch   = (TBranchClones*)(fEvtTree->GetBranch("frames"));

  clustersBranch->ResetReadEntry();
  segmentsBranch->ResetReadEntry();
  tracksBranch->ResetReadEntry();
  framesBranch->ResetReadEntry();

  headerBranch->SetAddress(fEvent->GetHeaderAddr());
  clustersBranch->SetAddress(fEvent->GetClustersAddr());
  segmentsBranch->SetAddress(fEvent->GetSegmentsAddr());
  tracksBranch->SetAddress(fEvent->GetTracksAddr());
  framesBranch->SetAddress(fEvent->GetFramesAddr());
}


//----------------------------------------------------------------------------
void EGraphRec::AddProcBrickFrame(TGVerticalFrame *workframe)
{
  // TODO: add the scanning options window

  TGLabel *label;
  TGHorizontalFrame *frame;

  TGVButtonGroup *GroupProcBrick = new TGVButtonGroup(workframe,
						      "Processing brick");

  // BrickId

  frame = new TGHorizontalFrame(GroupProcBrick);
  label = new TGLabel(frame, "Brick Id");
  fEntryProcBrickId = new TGNumberEntry(frame, fProcBrick.brickId, 7, 0, 
					TGNumberEntry::kNESInteger,
					TGNumberEntry::kNEANonNegative,
					TGNumberEntry::kNELLimitMin, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntryProcBrickId, fLayoutRightExpY);
  GroupProcBrick->AddFrame(frame, fLayout1);

  // Version

  frame = new TGHorizontalFrame(GroupProcBrick);
  label = new TGLabel(frame, "Base version");
  fEntryProcVer = new TGNumberEntry(frame, fProcBrick.ver, 3, 3, 
				    TGNumberEntry::kNESInteger,
				    TGNumberEntry::kNEANonNegative,
				    TGNumberEntry::kNELLimitMin, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntryProcVer, fLayoutRightExpY);
  GroupProcBrick->AddFrame(frame, fLayout1);

  workframe->AddFrame(GroupProcBrick, fLayout1);
}


//----------------------------------------------------------------------------
void EGraphRec::AddProcListFrame(TGVerticalFrame *workframe)
{
  TGVButtonGroup *GroupProcList = new TGVButtonGroup(workframe,"Process list");
  fCheckProcScan = new TGCheckButton(GroupProcList, "Scanning");
  fCheckProcLink = new TGCheckButton(GroupProcList, "Linking tracks");
  fCheckProcAlgn = new TGCheckButton(GroupProcList, "Alignment plates");
  fCheckProcTrks = new TGCheckButton(GroupProcList, "Reconstruct tracks");
  fCheckProcVrtx = new TGCheckButton(GroupProcList, "Reconstruct vertex");

//   fCheckProcLink->SetState(kButtonDown);
//   fCheckProcAlgn->SetState(kButtonDown);
  fCheckProcScan->SetEnabled(kFALSE);
  fCheckProcTrks->SetEnabled(kFALSE);
  fCheckProcVrtx->SetEnabled(kFALSE);

  workframe->AddFrame(GroupProcList, fLayout1);

  // Break process

  TGTextButton *endprocess = new TGTextButton(workframe, "Reset process");
  endprocess->Connect("Clicked()", "EGraphRec", this, "ResetProcess()");
  // endprocess->Associate(this);
  workframe->AddFrame(endprocess, fLayout1);

  // Running process

  fTextProcEvent = new TGTextButton(workframe, "Execute event");
  fTextProcEvent->Connect("Clicked()", "EGraphRec", this, "ProcessEvent()");

  workframe->AddFrame(fTextProcEvent, fLayout1);
}


//----------------------------------------------------------------------------
void EGraphRec::AddRecOptFrame(TGTab *worktab)
{
  // Create a tab with reconstruction options buttons

  TGCompositeFrame *tf = worktab->AddTab("Reconstruction options");
  TGLabel *label;
  TGHorizontalFrame *frame;

  // fRecParams = fGR->GetRecParameters();

  // Momentum

  frame = new TGHorizontalFrame(tf);
  label = new TGLabel(frame, "momentum");
  fRecOptEntry[0] = new TGNumberEntry(frame, 1, //fRecParams->p, 
				      7, 0, 
				      TGNumberEntry::kNESRealTwo,
				      TGNumberEntry::kNEANonNegative,
				      TGNumberEntry::kNELLimitMin, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fRecOptEntry[0], fLayoutRightExpY);
  tf->AddFrame(frame, fLayout1);

  // Prob min

  frame = new TGHorizontalFrame(tf);
  label  = new TGLabel(frame, "probmin");
  fRecOptEntry[1] = new TGNumberEntry(frame, 1, // fRecParams->probmin, 
				      7, 1,
				      TGNumberEntry::kNESRealTwo,
				      TGNumberEntry::kNEANonNegative,
				      TGNumberEntry::kNELLimitMinMax, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fRecOptEntry[1], fLayoutRightExpY);
  tf->AddFrame(frame, fLayout1);

  // Nseg min

  frame = new TGHorizontalFrame(tf);
  label  = new TGLabel(frame, "nsegmin");
  fRecOptEntry[2] = new TGNumberEntry(frame, 1, // fRecParams->nsegmin, 
				      7, 2,
				      TGNumberEntry::kNESInteger,
				      TGNumberEntry::kNEAPositive,
				      TGNumberEntry::kNELLimitMin, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fRecOptEntry[2], fLayoutRightExpY);
  tf->AddFrame(frame, fLayout1); 

  // Max gap

  frame = new TGHorizontalFrame(tf);
  label  = new TGLabel(frame, "maxgap");
  fRecOptEntry[3] = new TGNumberEntry(frame, 1, // fRecParams->maxgap, 
				      7, 3,
				      TGNumberEntry::kNESInteger,
				      TGNumberEntry::kNEAPositive,
				      TGNumberEntry::kNELLimitMin, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fRecOptEntry[3], fLayoutRightExpY);
  tf->AddFrame(frame, fLayout1);

  // Scanning conditions

  label = new TGLabel(tf,"Scanning conditions:");
  tf->AddFrame(label, fLayout1);

  // Coords sigma

  frame = new TGHorizontalFrame(tf);
  label  = new TGLabel(frame, "Coords Sigma");
  fRecOptEntry[4] = new TGNumberEntry(frame, 1, //fGR->GetScanCond()->SigmaX(0), 
				      7, 4,
				      TGNumberEntry::kNESRealTwo,
				      TGNumberEntry::kNEANonNegative,
				      TGNumberEntry::kNELLimitMinMax, 0, 10);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fRecOptEntry[4], fLayoutRightExpY);
  tf->AddFrame(frame, fLayout1);

  // Tangents sigma

  frame = new TGHorizontalFrame(tf);
  label  = new TGLabel(frame, "Tangents Sigma");
  fRecOptEntry[5] = new TGNumberEntry(frame,1, // fGR->GetScanCond()->SigmaTX(0), 
				      7, 5,
				      TGNumberEntry::kNESRealFour,
				      TGNumberEntry::kNEANonNegative,
				      TGNumberEntry::kNELLimitMinMax, 0, 1);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fRecOptEntry[5], fLayoutRightExpY);
  tf->AddFrame(frame, fLayout1);

  // buttons
    
  TGHorizontalFrame *Frame1 = new TGHorizontalFrame(tf);

  TGTextButton *Apply = new TGTextButton(Frame1, "Apply", 2);
  // Apply->Connect("Clicked()", "EGraphRec", this, "DoApplyRec()");
  Frame1->AddFrame(Apply, fLayout1);

  TGTextButton* Reset = new TGTextButton(Frame1, "Reset", 3);
  // Reset->Connect("Clicked()", "EGraphRec", this, "DoResetRec()");
  Frame1->AddFrame(Reset, fLayout1);

  tf->AddFrame(Frame1, fLayout1);
}


//----------------------------------------------------------------------------
void EGraphRec::AddScanBackFrame(TGTab *worktab)
{
  TGLabel           *label;
  TGHorizontalFrame *frame;

  // Create a tab for scan back procedure

  TGCompositeFrame *sb_tab = worktab->AddTab("Scan Back");

  // Parameters

  TGVButtonGroup *GroupSBPar = new TGVButtonGroup(sb_tab, "Parameters");

  // First plate to process

  frame = new TGHorizontalFrame(GroupSBPar);
  label = new TGLabel(frame, "First plate");
  fEntrySBFirstPlate = new TGNumberEntry(frame, fProcBrick.firstPlate, 3, 1, 
					 TGNumberEntry::kNESInteger,
					 TGNumberEntry::kNEANonNegative,
					 TGNumberEntry::kNELLimitMinMax,0,58);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntrySBFirstPlate, fLayoutRightExpY);
  GroupSBPar->AddFrame(frame, fLayout1);

  // Last plate to process

  frame = new TGHorizontalFrame(GroupSBPar);
  label = new TGLabel(frame, "Last plate");
  fEntrySBLastPlate = new TGNumberEntry(frame, fProcBrick.lastPlate, 3, 2, 
					TGNumberEntry::kNESInteger,
					TGNumberEntry::kNEANonNegative,
					TGNumberEntry::kNELLimitMinMax,0,58);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntrySBLastPlate, fLayoutRightExpY);
  GroupSBPar->AddFrame(frame, fLayout1);

  // Step

  frame = new TGHorizontalFrame(GroupSBPar);
  label = new TGLabel(frame, "Step");
  fEntrySBStep = new TGNumberEntry(frame, 1, 3, 3, 
				   TGNumberEntry::kNESInteger,
				   TGNumberEntry::kNEANonNegative,
				   TGNumberEntry::kNELLimitMinMax,0,58);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntrySBStep, fLayoutRightExpY);
  GroupSBPar->AddFrame(frame, fLayout1);

  sb_tab->AddFrame(GroupSBPar, fLayout1);

  // Scan back process list

  TGVButtonGroup *GroupSBProcList = new TGVButtonGroup(sb_tab, "Process list");
  fCheckSBScan = new TGCheckButton(GroupSBProcList, "Scanning");
  fCheckSBLink = new TGCheckButton(GroupSBProcList, "Linking");
  fCheckSBAlgn = new TGCheckButton(GroupSBProcList, "Alignment plates");
  fCheckSBTrks = new TGCheckButton(GroupSBProcList, "Reconstruct tracks");

//   fCheckSBLink->SetState(kButtonDown);
//   fCheckSBAlgn->SetState(kButtonDown);
  fCheckSBScan->SetEnabled(kFALSE);
  fCheckSBTrks->SetEnabled(kFALSE);

  sb_tab->AddFrame(GroupSBProcList, fLayout1);

  // Read prediction

  TGTextButton *readPred = new TGTextButton(sb_tab, "Read prediction");
  readPred->Connect("Clicked()", "EGraphRec", this, "ReadSBPred()");
  sb_tab->AddFrame(readPred, fLayout1);

  // start process

  fButtonSBStart = new TGTextButton(sb_tab, "Start Scan Back");
  fButtonSBStart->Connect("Clicked()","EGraphRec",this,"StartScanBack()");
  sb_tab->AddFrame(fButtonSBStart, fLayout1);
}


//----------------------------------------------------------------------------
void EGraphRec::AddVertexRecFrame(TGTab *worktab)
{
  // create a tab for vertex reconstruction

  TGLabel *label;
  TGHorizontalFrame *frame;
  TGCompositeFrame *vr_tab = worktab->AddTab("Vertex Rec");

  // use or not track momentum for vertex calculations

  fCheckUseMom = new TGCheckButton(vr_tab, "Use track momentum");
  if (fVertexRecOpt.UseMom) fCheckUseMom->SetState(kButtonDown, kTRUE);
  else fCheckUseMom->SetState(kButtonUp, kTRUE);
  vr_tab->AddFrame(fCheckUseMom, fLayout1);

  // use only the nearest measured segments for vertex fit

  fCheckUseSegPar = new TGCheckButton(vr_tab, "Use nearest measured segments");
  if (fVertexRecOpt.UseSegPar) fCheckUseSegPar->SetState(kButtonDown, kTRUE);
  else fCheckUseSegPar->SetState(kButtonUp, kTRUE);
  vr_tab->AddFrame(fCheckUseSegPar, fLayout1);

  // vertex quality estimation method

  frame = new TGHorizontalFrame(vr_tab);
  label = new TGLabel(frame, "Quality Mode");
  fEntryQualityMode = new TGNumberEntry(frame, fVertexRecOpt.QualityMode,
					5, 0, 
					TGNumberEntry::kNESInteger,
					TGNumberEntry::kNEANonNegative,
					TGNumberEntry::kNELLimitMin);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntryQualityMode, fLayoutRightExpY);
  vr_tab->AddFrame(frame, fLayout1);

  // maximum z-gap in the track-vertex group

  frame = new TGHorizontalFrame(vr_tab);
  label = new TGLabel(frame, "Maximum Z gap");
  fEntryDZmax = new TGNumberEntry(frame, fVertexRecOpt.DZmax,
				  7, 0, 
				  TGNumberEntry::kNESRealOne,
				  TGNumberEntry::kNEANonNegative,
				  TGNumberEntry::kNELLimitMin);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntryDZmax, fLayoutRightExpY);
  vr_tab->AddFrame(frame, fLayout1);

  // minimum acceptable probability for chi2-distance between tracks

  frame = new TGHorizontalFrame(vr_tab);
  label = new TGLabel(frame, "Minimal vertex prob");
  fEntryProbMinV = new TGNumberEntry(frame, fVertexRecOpt.ProbMinV,
				     7, 0, 
				     TGNumberEntry::kNESReal,
				     TGNumberEntry::kNEANonNegative,
				     TGNumberEntry::kNELLimitMin);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntryProbMinV, fLayoutRightExpY);
  vr_tab->AddFrame(frame, fLayout1);


  // maximal acceptable impact parameter

  frame = new TGHorizontalFrame(vr_tab);
  label = new TGLabel(frame, "Maximal impact param");
  fEntryImpMax = new TGNumberEntry(frame, fVertexRecOpt.ImpMax,
				   7, 0, 
				   TGNumberEntry::kNESRealTwo,
				   TGNumberEntry::kNEANonNegative,
				   TGNumberEntry::kNELLimitMin);
  frame->AddFrame(label, fLayoutLeftExpY);
  frame->AddFrame(fEntryImpMax, fLayoutRightExpY);
  vr_tab->AddFrame(frame, fLayout1);

  // start process

  TGTextButton *buttonVRStart = new TGTextButton(vr_tab, "Start Vertex Rec");
  buttonVRStart->Connect("Clicked()","EGraphRec",this,"StartVertexRec()");
  vr_tab->AddFrame(buttonVRStart, fLayout1);
}


//----------------------------------------------------------------------------
void EGraphRec::AddCanvasFrame(TGTab *worktab)
{
  TGCompositeFrame *tf;

  tf = worktab->AddTab("cFedra");
  fDisplayFedra = new TRootEmbeddedCanvas("FEDRA Viewer", tf, 900, 650);
  tf->AddFrame(fDisplayFedra, fLayout3);

  TCanvas *CanvasFedra = fDisplayFedra->GetCanvas();

  fDisplay = new EdbDisplay("FEDRA Viewer (orig)", -50000., 50000., 
			    -50000., 50000., -150., 37850., CanvasFedra);

//   tf = worktab->AddTab("GL Viewer");
//   fDisplayHitsGL = new TRootEmbeddedCanvas("Reconstruction GL", tf, 900, 650);
//   fGLViewer = new TGLSAViewer(tf, fDisplayHitsGL->GetCanvas());

  // tf->AddFrame(fDisplayHitsGL, fLayout3);
}


//----------------------------------------------------------------------------
void EGraphRec::AddInfoFrame(TGVerticalFrame *workframe)
{
  fTextInfo = new TGTextView(workframe, 900, 200, kFixedHeight);
  const TGFont *font = gClient->GetFont("-*-courier-bold-r-*-*-18-*-*-*-*-*-*-*");
  fTextInfo->SetFont(font->GetFontStruct());
  workframe->AddFrame(fTextInfo, fLayout1);

  WriteInfo();
}


//----------------------------------------------------------------------------
void EGraphRec::ReadSBPred()
{
  // const char *filetypes[] = {"ROOT files", "*.root", 0, 0};

  static TString dir(".");
  TGFileInfo fi;
  // fi.fFileTypes = filetypes;
  fi.fIniDir    = StrDup(dir.Data());
  new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, &fi);
  dir = fi.fIniDir;

  fPredTracks->Reset();

  // open file with predictions

  if (fi.fFilename) fScanProc->ReadPatTXT(fi.fFilename, *fPredTracks);

  // Create directory and write the predictions to the root file

}


//----------------------------------------------------------------------------
void EGraphRec::Set3DViewer()
{
  fGLViewer = (TGLSAViewer*)fDisplay->fPad->GetViewer3D("ogl");
  fGLViewer->SetClearColor(38);
}


void EGraphRec::ZoomIn()
{
//   fDisplayHits->GetCanvas()->cd();
//   fDisplayHits->GetCanvas()->GetView()->Zoom();
//   fDisplayHits->GetCanvas()->Modified();
//   fDisplayHits->GetCanvas()->Update();
}


//----------------------------------------------------------------------------
void EGraphRec::ZoomOut()
{
//   fDisplayHits->GetCanvas()->cd();
//   fDisplayHits->GetCanvas()->GetView()->UnZoom();
//   fDisplayHits->GetCanvas()->Modified();
//   fDisplayHits->GetCanvas()->Update();
}


//----------------------------------------------------------------------------
void EGraphRec::WriteInfo()
{
  fTextInfo->Clear();
  fTextInfo->AddLine("");
  fTextInfo->AddLine(" Brick data directory: " + fDataDir);
}


//----------------------------------------------------------------------------
void EGraphRec::InitScanSet()
{
  fScanProc->eProcDirClient = fDataDir;  // brick directory initialization

  // Scan set initialization

//   if (fScanSet && fScanSet->Brick().ID() != fBrickToProc.brickId)
//     SafeDelete(fScanSet);

//   if (!fScanSet) {
//     fScanSet = new EdbScanSet();
//     fScanSet->Brick().SetID(fBrickToProc.brickId);
//   }

  SafeDelete(fScanSet);
  fScanSet = new EdbScanSet();
  fScanSet->Brick().SetID(fBrickToProc.brickId);
}


//----------------------------------------------------------------------------
void EGraphRec::InitVariables()
{
  fEvent            = NULL;
  fFoundTracks      = NULL;
  fPVRec            = NULL;
  fThSBProcess      = NULL;
  fThSBCheckProcess = NULL;
  fScanSet          = NULL;
  fThSBProcess      = NULL;
  fDataDir          = "";

  fProcBrick.brickId = fProcBrick.firstPlate = fProcBrick.lastPlate = 
    fProcBrick.ver = -1;
  fProcId.interCalib = fProcId.volumeScan = fProcId.predScan = 
    fProcId.scanForth = -1;

  // init default vertex reconstruction options

  fVertexRecOpt.QualityMode = 0;
  fVertexRecOpt.UseMom      = kTRUE;
  fVertexRecOpt.UseSegPar   = kFALSE;
  fVertexRecOpt.DZmax       = 3000.;
  fVertexRecOpt.ProbMinV    = 0.001;
  fVertexRecOpt.ImpMax      = 20.;

  fScanProc     = new EdbScanProc();
  fPredTracks   = new EdbPattern();
  fRecProc      = new EGraphRecProc();

  // TThread functions

  // fThSBProcess      = new TThread("ThSBProcess", ThSBProcess, (void*) this);
  fThSBCheckProcess = new TThread("ThSBCheckProcess",
				  ThSBCheckProcess, (void*) this);

  gEDBDEBUGLEVEL = 2;
}

//----------------------------------------------------------------------------
void EGraphRec::InitDrawVariables()
{
  fLayout1 = new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2);
  fLayout2 = new TGLayoutHints(kLHintsExpandY, 2, 2, 2, 2);
  fLayout3 = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);
  fLayoutLeftExpY  = new TGLayoutHints(kLHintsLeft|kLHintsExpandY, 2, 4, 0, 0);
  fLayoutRightExpY = new TGLayoutHints(kLHintsRight|kLHintsExpandY,4, 2, 0, 0);
}


//----------------------------------------------------------------------------
void EGraphRec::ReadCmdConfig()
{
  TString configDir  = (TString)getenv("FEDRA_ROOT") + "/config/";
  TString configFile = configDir + "EGraphRec";

  if (gEnv->ReadFile(configFile + ".cfg", kEnvChange) == -1) 
    if (gEnv->ReadFile(configFile + "_default.cfg", kEnvChange) == -1) return;

  if (fDataDir == "")
    fDataDir = gEnv->GetValue("EGraphRec.DataDir", "");

  // brick Id

  if (fProcBrick.brickId < 0)
    fProcBrick.brickId = gEnv->GetValue("EGraphRec.ProcBrick.brickId", 0);
  if (fProcBrick.firstPlate < 0)
    fProcBrick.firstPlate = gEnv->GetValue("EGraphRec.ProcBrick.firstPlate", 0);
  if (fProcBrick.lastPlate < 0)
    fProcBrick.lastPlate = gEnv->GetValue("EGraphRec.ProcBrick.lastPlate", 0);
  if (fProcBrick.ver < 0)
    fProcBrick.ver = gEnv->GetValue("EGraphRec.ProcBrick.ver", 0);

  // processes Ids

  if (fProcId.interCalib < 0)
    fProcId.interCalib = gEnv->GetValue("EGraphRec.ProcId.interCalib", 0);
  if (fProcId.volumeScan < 0)
    fProcId.volumeScan = gEnv->GetValue("EGraphRec.ProcId.volumeScan", 0);
  if (fProcId.predScan < 0)
    fProcId.predScan = gEnv->GetValue("EGraphRec.ProcId.predScan", 0);
  if (fProcId.scanForth < 0)
    fProcId.scanForth = gEnv->GetValue("EGraphRec.ProcId.scanForth", 0);

  // vertex reconstruction options

  fVertexRecOpt.QualityMode = gEnv->GetValue("VertexRecOpt.QualityMode", 
					     fVertexRecOpt.QualityMode);
  fVertexRecOpt.UseMom      = gEnv->GetValue("VertexRecOpt.UseMom", 
					     fVertexRecOpt.UseMom);
  fVertexRecOpt.UseSegPar   = gEnv->GetValue("VertexRecOpt.UseSegPar", 
					     fVertexRecOpt.UseSegPar);
  fVertexRecOpt.DZmax       = gEnv->GetValue("VertexRecOpt.DZmax", 
					     fVertexRecOpt.DZmax);
  fVertexRecOpt.ProbMinV    = gEnv->GetValue("VertexRecOpt.ProbMinV", 
					     fVertexRecOpt.ProbMinV);
  fVertexRecOpt.ImpMax      = gEnv->GetValue("VertexRecOpt.ImpMax", 
					     fVertexRecOpt.ImpMax);
}


//----------------------------------------------------------------------------
void EGraphRec::ClearEvent()
{
  SafeDelete(fEvent);
}
