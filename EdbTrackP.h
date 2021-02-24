#ifndef ROOT_EdbPattern
#define ROOT_EdbPattern
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbPattern                                                           //
//                                                                      //
// pattern of segments                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TSortedList.h"
#include "TArrayL.h"
#include "TClonesArray.h"
#include "EdbSegP.h"

class EdbAffine2D;
class TIndexCell;
class EdbVTA;
class EdbVertex;


//______________________________________________________________________________
class EdbSegmentsBox : public TObject, public EdbPointsBox2D {
 
 private:
 
  Float_t eX,eY,eZ;         // central point

  TClonesArray *eSegments;  // collection of segments (EdbSegP)

  Float_t eDZkeep;          // eDZkeep = eZ - Zoriginal (before any projections)

  //  void (*PF_COUPLE)(EdbSegP *seg1,EdbSegP *seg2); //pointer to the function

  //protected:
 public:
  void       SetZ(float z) { eZ=z; }

 public:
  EdbSegmentsBox(int nseg=0);
  EdbSegmentsBox(EdbSegmentsBox &box)
    {
      eX=box.X();
      eY=box.Y();
      eZ=box.Z();
      int n=box.N();
      eSegments = new TClonesArray("EdbSegP",n);
      for(int i=0; i<n; i++) 
	AddSegment( *(box.GetSegment(i)) );
    }

  EdbSegmentsBox(float x0, float y0, float z0, int nseg=0);
  virtual ~EdbSegmentsBox();
 
  EdbSegP  *AddSegment(int i, EdbSegP &s);
  EdbSegP  *AddSegment(EdbSegP &s);
  EdbSegP  *AddSegment(EdbSegP &s1, EdbSegP &s2); 
  EdbSegP  *AddSegment(int id, float x, float y, float tx, float ty, 
		      float w=0, int flag=0);
  Int_t     GetN() const {return eSegments ? eSegments->GetEntriesFast() : 0;}
  EdbSegP  *GetSegment(int i) const {return (EdbSegP*)eSegments->At(i); }
  EdbSegP  *GetSegmentLast()  const {return GetN() ? (EdbSegP*)eSegments->Last() : 0;}

  TClonesArray *GetSegments()     const { return eSegments; }
  void         *GetSegmentsAddr()       { return &eSegments; }

  void      SetSegmentsZ();
  void      SetSegmentsDZ(float dz);

  void      SetSegmentsPlate(int plate);

  // mandatory virtual functions:

  void      SetX(float x) { eX=x; }
  void      SetY(float y) { eY=y; }

  Float_t   X()       const {return eX;}
  Float_t   Y()       const {return eY;}
  Float_t   Z()       const {return eZ;}
  Float_t   DZ()      const {return eDZkeep;}
  Int_t     N()       const {return GetN();}
  EdbPoint *At(int i) const {return (EdbPoint*)GetSegment(i);}
 
  // other functions

  void    Print( Option_t *opt="") const;
  void    Set0();
  void    Reset();
  void    ProjectTo(const float dz);
  void    TransformA(    const EdbAffine2D *affA );
  void    TransformARot( const EdbAffine2D *affA );
  void    TransformShr(  const float shr );
  Int_t   CalculateXY(  EdbSegmentsBox *p , EdbAffine2D *aff  );
  Int_t   CalculateAXAY(EdbSegmentsBox *p , EdbAffine2D *affA );
  Float_t DiffAff( EdbAffine2D *aff );
  Float_t Diff( EdbSegmentsBox &p );
  Float_t GetSize(Int_t XorY);
  Float_t GetSizeX() { return GetSize(0); }
  Float_t GetSizeY() { return GetSize(1); }
  Float_t GetSizeXY();

  ClassDef(EdbSegmentsBox,1)  // collection of segments
};

//______________________________________________________________________________
class EdbTrackP : public EdbSegP {
 
 private:

  TSortedList *eS;     // array of segments
  TSortedList *eSF;    // array of fitted segments
  Int_t        eNpl;   // number of plates passed throw
  Int_t        eN0;    // number of holes (if any)
  Float_t      eM;     // invariant mass of the particle
  Float_t      eDE;    // total energy loss of the particle between first and last segments
  Int_t        ePDG;   // particle ID from PDG
  Float_t      ePerrUp;    // error of P() in upper direction, obtained by MCS,or shower-algorithm
  Float_t      ePerrDown;    // error of P() in lower direction, obtained by MCS,or shower-algorithm
  EdbVTA      *eVTAS;  //! vertex track start is attached to
  EdbVTA      *eVTAE;  //! vertex track end is attached to

 public:

  EdbTrackP(int nseg=0);
  EdbTrackP(EdbSegP &seg);
  EdbTrackP(EdbSegP *seg, float m=0.12);
  EdbTrackP(EdbTrackP &track) { Set0(); Copy(track); }

  virtual ~EdbTrackP();

  void       Set0();
  void       SetOwner() { if(eS) eS->SetOwner(true); }
  void       AddVTA(EdbVTA *vta);
  void       ClearVTA();
  void       ClearVTA(EdbVTA *vta);
  EdbVTA    *VTAS() const {return eVTAS;}
  EdbVTA    *VTAE() const {return eVTAE;}
  EdbVertex *VertexS();
  EdbVertex *VertexE();
  EdbVertex *Vertex(int zpos) {return zpos? VertexS(): VertexE();}

  void    SetPDG( int pdg )  { ePDG=pdg; }
  Int_t   PDG()      const {return ePDG;}
  
  Float_t Wmean()      const;

  void    SetM( float m )  { eM=m; }
  Float_t M()      const {return eM;}
  

  //void    SetErrorsCOV(EdbScanCond &cond);
  void    SetCounters() { SetNpl(); SetN0(); }

  void    SetN0( int n0 )  { eN0 = n0; }
  void    SetN0()          { eN0 = eNpl-N(); }
  Int_t   N0()      const  {return eN0;}

  void    SetDE( float de )  { eDE=de; }
  Float_t DE()      const  {return eDE;}

  void    SetNpl( int npl )  { eNpl=npl; }
  void    SetNpl()  
    { if(eS) eNpl = 1+TMath::Abs(GetSegment(0)->PID() - GetSegment(N()-1)->PID()); }
  Int_t   Npl()      const   {return eNpl;}

  //  TList *S()  const { return eS; }
  //  TList *SF() const { return eSF; }

  //int      N() const  { if(eS)  return eS->GetEntries(); else return 0; } //TODO fast

  Int_t    N()  const  { return (eS) ?  eS->GetSize()  : 0; }
  Int_t    NF() const  { return (eSF)?  eSF->GetSize() : 0; }

  float    Wgrains() const;
  int	   GetSegmentsFlag( int &nseg ) const;
  int	   GetSegmentsAid( int &nseg ) const;
  int	   GetSegmentsMCTrack( int &nseg ) const;

  EdbSegP *GetSegmentWithClosestZ( float z, float dz);
  
  EdbSegP *GetSegmentFirst()  const {return (eS) ? (EdbSegP*)(eS->First()) : 0;}
  EdbSegP *GetSegmentLast()   const {return (eS) ? (EdbSegP*)(eS->Last()) : 0;}

  EdbSegP *GetSegmentFFirst() const {return (eSF) ? (EdbSegP*)(eSF->First()) : 0;}
  EdbSegP *GetSegmentFLast()  const {return (eSF) ? (EdbSegP*)(eSF->Last()) : 0;}

  EdbSegP *GetSegment(int i)  const {return (eS) ? (EdbSegP*)(eS->At(i)) : 0; }
  EdbSegP *GetSegmentF(int i) const {return (eSF) ? (EdbSegP*)(eSF->At(i)) : 0;}

  EdbSegP *TrackZmin(bool usesegpar=false) const { if(usesegpar || (!eSF)) return GetSegmentFirst(); else return GetSegmentFFirst(); }
  EdbSegP *TrackZmax(bool usesegpar=false) const { if(usesegpar || (!eSF)) return GetSegmentLast();  else return GetSegmentFLast(); }

  EdbSegP *TrackExtremity(bool zpos, bool usesegpar=false) 
    const { return zpos? TrackZmin(usesegpar) : TrackZmax(usesegpar); } // 0-end, 1-strart (as in vertex class)

  const EdbSegP *TrackStart() const;
  const EdbSegP *TrackEnd()   const;

  Int_t   Dir()    const {return (DZ()<0) ? -1 : 1;}
  Float_t Zmax()   const;
  Float_t Zmin()   const;
  Float_t Zstart() const {return TrackStart()->Z();}
  Float_t Zend()   const {return TrackEnd()->Z();}
  void    AddTrack(const EdbTrackP &tr);

  void  AddSegment(EdbSegP *s)
  { 
    if(!eS) eS = new TSortedList();
    eS->Add(s);
  }
  void  RemoveSegment(EdbSegP *s)
  { 
    if(!eS) return;
    eS->Remove(s);
    s->SetTrack(-1);
    SetCounters();
  }
  void  SubstituteSegment( EdbSegP *sold, EdbSegP *snew )
  { 
    if(!eS) return;
    eS->Remove(sold);
    eS->Add(snew);
    sold->SetTrack(-1);
  }
  void  AddSegmentF(EdbSegP *s)  
  { 
    if(!eSF) 
      {
        eSF = new TSortedList();
	eSF->SetOwner();
      }
    eSF->Add(s);
  }

  int   RemoveAliasSegments();
  int   CheckMaxGap();
  int   CheckAliasSegments();
  int   SetSegmentsTrack(int id) {for(int i=0; i<N(); i++) GetSegment(i)->SetTrack(id); return N();}
  int   SetSegmentsTrack() {return SetSegmentsTrack(ID());}
  int   FitTrackKFS( bool zmax=false, float X0=5810., int design = 0 );
  int   MakeSelector( EdbSegP &ss, bool followZ=true );
  float MakePredictionTo( Float_t z, EdbSegP &ss );
  float CHI2();
  float CHI2F();
  void  FitTrack();
  void  Copy(const EdbTrackP &tr);
  void  Transform(const EdbAffine2D &tr);
  
  void SetPerrUp(Float_t perrUp);
  void SetPerrDown(Float_t perrDown);
  void SetPerr(Float_t perrDown, Float_t perrUp);
  Float_t PerrUp()      const {return ePerrUp;}
  Float_t PerrDown()      const {return ePerrDown;}
  

  void  Clear()  { if(eS)  eS->Clear(); if(eSF) eSF->Clear(); }
  void  ClearF() { if(eSF) eSF->Clear(); }
  void  Print();
  void  PrintNice();

  ClassDef(EdbTrackP,7)  // track consists of segments
};

//______________________________________________________________________________
class EdbPattern : public EdbSegmentsBox {
 
 private:

  Int_t       eID;             // pattern id in the volume
  Int_t       ePID;            // correspond to the piece ID
  TIndexCell *eCell;           //! associated with eSegments
  Float_t     eStepX ,eStepY;  // bin size for the index cell 
  Float_t     eStepTX,eStepTY; // bin size for the index cell 
  Float_t     eSigma[4];       // accuracy in comparison to neibour pattern
  Int_t       eFlag;           // pattern flag
  Int_t       eNAff;           // number of segments selected for affine calculation
  EdbID       eScanID;         //! main scanning ID for this pattern, defined when possible
  Int_t       eSide;           //! emulsion side 0/1/2, defined when possible

 public:

  EdbPattern();
  EdbPattern( EdbPattern &p );
  EdbPattern(float x0, float y0, float z0, int n=0 );
  virtual ~EdbPattern();
 
  void  Set0();
  void  SetScanID(EdbID id) {eScanID=id;}
  void  SetSigma(float sx,float sy,float stx,float sty) 
                {eSigma[0]=sx; eSigma[1]=sy; eSigma[2]=stx; eSigma[3]=sty;}
  void  SetStep(float stepx, float stepy, float steptx, float stepty ) 
                {eStepX=stepx; eStepY=stepy; eStepTX=steptx; eStepTY=stepty;}
  void  Reset();
  float SummaryPath();
  void  FillCell( float stepx, float stepy, float steptx, float stepty );
  int   FindCompliments(EdbSegP &s, TObjArray &arr, float nsig, float nsigt);
  void  SetSegmentsPID();
  void  SetSegmentsScanID(EdbID id);
  EdbPattern *ExtractSubPattern(float min[5], float max[5], int MCevt=-1);

  void  SetID(int id)   {eID=id;}
  void  SetPID(int pid) {ePID=pid;}
  void  SetNAff(int n)  {eNAff=n;}
  void  SetSide(int side)  {eSide=side;}

  Int_t       NAff()   const {return eNAff;}
  float       StepX()  const {return eStepX;}
  float       StepY()  const {return eStepY;}
  float       StepTX() const {return eStepTX;}
  float       StepTY() const {return eStepTY;}
  int         ID()     const {return eID;}
  int         PID()    const {return ePID;}
  TIndexCell *Cell()   const {return eCell;}
  EdbSegP    *FindSegment(int id);
  float      Xmean();
  float      Ymean();

  EdbID      ScanID() const {return eScanID;}
  Int_t      Plate() const {return eScanID.ePlate;}
  Int_t      Side()  const {return eSide;}

  ClassDef(EdbPattern,2)  // pattern of segments
};

//______________________________________________________________________________
class EdbPatternsVolume : public TObject {
 
 private:
 
  Float_t    eX,eY,eZ;       // central point

 public:

  TObjArray  *ePatterns;      // collection of patterns

  TIndexCell *eTracksCell;    //! "vidt:vids" - connected segments cell
  TIndexCell *ePatternsCell;  //! "pid:id1:chi2:id2" - all found couples

  Bool_t      eDescendingZ;    // if =0 - z increase in the pattrens array; if =1 - decrease 

 public:
  EdbPatternsVolume();
  EdbPatternsVolume(EdbPatternsVolume &pvol);
  virtual ~EdbPatternsVolume();

  void Set0();
  void SetPatternsID();

  void Transform( const EdbAffine2D *aff );
  void Shift(float x, float y);
  void Centralize();
  void Centralize(float xc, float yc);
  void SetXYZ( float x, float y, float z ) { eX=x; eY=y; eZ=z; }

  Float_t X() const {return eX;}
  Float_t Y() const {return eY;}
  Float_t Z() const {return eZ;}
  Int_t   Npatterns() const { if(ePatterns) 
                        return ePatterns->GetEntriesFast();
                        else return 0; }
  Float_t      Xmean();
  Float_t      Ymean();

  void PassProperties(EdbPatternsVolume &pvol);

  void         AddPattern( EdbPattern *pat );
  void         AddPatternAt( EdbPattern *pat, int id );
  EdbPattern  *GetPattern( int id ) const;
  EdbPattern  *GetPatternZLowestHighest(Bool_t lowestZ=kTRUE) const;
	EdbPattern* GetPatternPreceding(EdbPattern* pat) const;
	EdbPattern* GetPatternSucceding(EdbPattern* pat) const;
	EdbPattern* NextPattern(float z, int dir) const;
  EdbPattern* GetPatternNext(float z, int dir) const;
	EdbPattern* GetPatternByPID(int pid) const;
	EdbPattern* GetPatternByZ(float z) const;

  void   DropCell();

  void   PrintStat( Option_t *opt="") const;
  void   PrintStat( EdbPattern &pat ) const;
  void   PrintAff() const;
  int    DropCouples();

  Long_t Vid(int pid, int sid) const { return pid*1000000+sid; }
  Int_t  Pid(Long_t vid)       const { return vid/1000000; }
  Int_t  Sid(Long_t vid)       const { return vid%1000000; }

  EdbSegP *GetSegment(Long_t vid) const 
    {return GetPattern(Pid(vid))->GetSegment( Sid(vid) );}

  int FindComplimentsVol(EdbSegP &s, TObjArray &arr, float nsig, float nsigt, int dpat);

  EdbPattern  *GetPatternByPlate(int plate, int side);
  EdbPattern  *InsertPattern(EdbPattern *pat, Bool_t descendingZ=0);
  void         SortPatternsByZ(Bool_t descendingZ=0);

  void   Print() const;

  ClassDef(EdbPatternsVolume,1)  // patterns nostri
};

#endif /* ROOT_EdbPattern */
 EdbPattern.h:1
 EdbPattern.h:2
 EdbPattern.h:3
 EdbPattern.h:4
 EdbPattern.h:5
 EdbPattern.h:6
 EdbPattern.h:7
 EdbPattern.h:8
 EdbPattern.h:9
 EdbPattern.h:10
 EdbPattern.h:11
 EdbPattern.h:12
 EdbPattern.h:13
 EdbPattern.h:14
 EdbPattern.h:15
 EdbPattern.h:16
 EdbPattern.h:17
 EdbPattern.h:18
 EdbPattern.h:19
 EdbPattern.h:20
 EdbPattern.h:21
 EdbPattern.h:22
 EdbPattern.h:23
 EdbPattern.h:24
 EdbPattern.h:25
 EdbPattern.h:26
 EdbPattern.h:27
 EdbPattern.h:28
 EdbPattern.h:29
 EdbPattern.h:30
 EdbPattern.h:31
 EdbPattern.h:32
 EdbPattern.h:33
 EdbPattern.h:34
 EdbPattern.h:35
 EdbPattern.h:36
 EdbPattern.h:37
 EdbPattern.h:38
 EdbPattern.h:39
 EdbPattern.h:40
 EdbPattern.h:41
 EdbPattern.h:42
 EdbPattern.h:43
 EdbPattern.h:44
 EdbPattern.h:45
 EdbPattern.h:46
 EdbPattern.h:47
 EdbPattern.h:48
 EdbPattern.h:49
 EdbPattern.h:50
 EdbPattern.h:51
 EdbPattern.h:52
 EdbPattern.h:53
 EdbPattern.h:54
 EdbPattern.h:55
 EdbPattern.h:56
 EdbPattern.h:57
 EdbPattern.h:58
 EdbPattern.h:59
 EdbPattern.h:60
 EdbPattern.h:61
 EdbPattern.h:62
 EdbPattern.h:63
 EdbPattern.h:64
 EdbPattern.h:65
 EdbPattern.h:66
 EdbPattern.h:67
 EdbPattern.h:68
 EdbPattern.h:69
 EdbPattern.h:70
 EdbPattern.h:71
 EdbPattern.h:72
 EdbPattern.h:73
 EdbPattern.h:74
 EdbPattern.h:75
 EdbPattern.h:76
 EdbPattern.h:77
 EdbPattern.h:78
 EdbPattern.h:79
 EdbPattern.h:80
 EdbPattern.h:81
 EdbPattern.h:82
 EdbPattern.h:83
 EdbPattern.h:84
 EdbPattern.h:85
 EdbPattern.h:86
 EdbPattern.h:87
 EdbPattern.h:88
 EdbPattern.h:89
 EdbPattern.h:90
 EdbPattern.h:91
 EdbPattern.h:92
 EdbPattern.h:93
 EdbPattern.h:94
 EdbPattern.h:95
 EdbPattern.h:96
 EdbPattern.h:97
 EdbPattern.h:98
 EdbPattern.h:99
 EdbPattern.h:100
 EdbPattern.h:101
 EdbPattern.h:102
 EdbPattern.h:103
 EdbPattern.h:104
 EdbPattern.h:105
 EdbPattern.h:106
 EdbPattern.h:107
 EdbPattern.h:108
 EdbPattern.h:109
 EdbPattern.h:110
 EdbPattern.h:111
 EdbPattern.h:112
 EdbPattern.h:113
 EdbPattern.h:114
 EdbPattern.h:115
 EdbPattern.h:116
 EdbPattern.h:117
 EdbPattern.h:118
 EdbPattern.h:119
 EdbPattern.h:120
 EdbPattern.h:121
 EdbPattern.h:122
 EdbPattern.h:123
 EdbPattern.h:124
 EdbPattern.h:125
 EdbPattern.h:126
 EdbPattern.h:127
 EdbPattern.h:128
 EdbPattern.h:129
 EdbPattern.h:130
 EdbPattern.h:131
 EdbPattern.h:132
 EdbPattern.h:133
 EdbPattern.h:134
 EdbPattern.h:135
 EdbPattern.h:136
 EdbPattern.h:137
 EdbPattern.h:138
 EdbPattern.h:139
 EdbPattern.h:140
 EdbPattern.h:141
 EdbPattern.h:142
 EdbPattern.h:143
 EdbPattern.h:144
 EdbPattern.h:145
 EdbPattern.h:146
 EdbPattern.h:147
 EdbPattern.h:148
 EdbPattern.h:149
 EdbPattern.h:150
 EdbPattern.h:151
 EdbPattern.h:152
 EdbPattern.h:153
 EdbPattern.h:154
 EdbPattern.h:155
 EdbPattern.h:156
 EdbPattern.h:157
 EdbPattern.h:158
 EdbPattern.h:159
 EdbPattern.h:160
 EdbPattern.h:161
 EdbPattern.h:162
 EdbPattern.h:163
 EdbPattern.h:164
 EdbPattern.h:165
 EdbPattern.h:166
 EdbPattern.h:167
 EdbPattern.h:168
 EdbPattern.h:169
 EdbPattern.h:170
 EdbPattern.h:171
 EdbPattern.h:172
 EdbPattern.h:173
 EdbPattern.h:174
 EdbPattern.h:175
 EdbPattern.h:176
 EdbPattern.h:177
 EdbPattern.h:178
 EdbPattern.h:179
 EdbPattern.h:180
 EdbPattern.h:181
 EdbPattern.h:182
 EdbPattern.h:183
 EdbPattern.h:184
 EdbPattern.h:185
 EdbPattern.h:186
 EdbPattern.h:187
 EdbPattern.h:188
 EdbPattern.h:189
 EdbPattern.h:190
 EdbPattern.h:191
 EdbPattern.h:192
 EdbPattern.h:193
 EdbPattern.h:194
 EdbPattern.h:195
 EdbPattern.h:196
 EdbPattern.h:197
 EdbPattern.h:198
 EdbPattern.h:199
 EdbPattern.h:200
 EdbPattern.h:201
 EdbPattern.h:202
 EdbPattern.h:203
 EdbPattern.h:204
 EdbPattern.h:205
 EdbPattern.h:206
 EdbPattern.h:207
 EdbPattern.h:208
 EdbPattern.h:209
 EdbPattern.h:210
 EdbPattern.h:211
 EdbPattern.h:212
 EdbPattern.h:213
 EdbPattern.h:214
 EdbPattern.h:215
 EdbPattern.h:216
 EdbPattern.h:217
 EdbPattern.h:218
 EdbPattern.h:219
 EdbPattern.h:220
 EdbPattern.h:221
 EdbPattern.h:222
 EdbPattern.h:223
 EdbPattern.h:224
 EdbPattern.h:225
 EdbPattern.h:226
 EdbPattern.h:227
 EdbPattern.h:228
 EdbPattern.h:229
 EdbPattern.h:230
 EdbPattern.h:231
 EdbPattern.h:232
 EdbPattern.h:233
 EdbPattern.h:234
 EdbPattern.h:235
 EdbPattern.h:236
 EdbPattern.h:237
 EdbPattern.h:238
 EdbPattern.h:239
 EdbPattern.h:240
 EdbPattern.h:241
 EdbPattern.h:242
 EdbPattern.h:243
 EdbPattern.h:244
 EdbPattern.h:245
 EdbPattern.h:246
 EdbPattern.h:247
 EdbPattern.h:248
 EdbPattern.h:249
 EdbPattern.h:250
 EdbPattern.h:251
 EdbPattern.h:252
 EdbPattern.h:253
 EdbPattern.h:254
 EdbPattern.h:255
 EdbPattern.h:256
 EdbPattern.h:257
 EdbPattern.h:258
 EdbPattern.h:259
 EdbPattern.h:260
 EdbPattern.h:261
 EdbPattern.h:262
 EdbPattern.h:263
 EdbPattern.h:264
 EdbPattern.h:265
 EdbPattern.h:266
 EdbPattern.h:267
 EdbPattern.h:268
 EdbPattern.h:269
 EdbPattern.h:270
 EdbPattern.h:271
 EdbPattern.h:272
 EdbPattern.h:273
 EdbPattern.h:274
 EdbPattern.h:275
 EdbPattern.h:276
 EdbPattern.h:277
 EdbPattern.h:278
 EdbPattern.h:279
 EdbPattern.h:280
 EdbPattern.h:281
 EdbPattern.h:282
 EdbPattern.h:283
 EdbPattern.h:284
 EdbPattern.h:285
 EdbPattern.h:286
 EdbPattern.h:287
 EdbPattern.h:288
 EdbPattern.h:289
 EdbPattern.h:290
 EdbPattern.h:291
 EdbPattern.h:292
 EdbPattern.h:293
 EdbPattern.h:294
 EdbPattern.h:295
 EdbPattern.h:296
 EdbPattern.h:297
 EdbPattern.h:298
 EdbPattern.h:299
 EdbPattern.h:300
 EdbPattern.h:301
 EdbPattern.h:302
 EdbPattern.h:303
 EdbPattern.h:304
 EdbPattern.h:305
 EdbPattern.h:306
 EdbPattern.h:307
 EdbPattern.h:308
 EdbPattern.h:309
 EdbPattern.h:310
 EdbPattern.h:311
 EdbPattern.h:312
 EdbPattern.h:313
 EdbPattern.h:314
 EdbPattern.h:315
 EdbPattern.h:316
 EdbPattern.h:317
 EdbPattern.h:318
 EdbPattern.h:319
 EdbPattern.h:320
 EdbPattern.h:321
 EdbPattern.h:322
 EdbPattern.h:323
 EdbPattern.h:324
 EdbPattern.h:325
 EdbPattern.h:326
 EdbPattern.h:327
 EdbPattern.h:328
 EdbPattern.h:329
 EdbPattern.h:330
 EdbPattern.h:331
 EdbPattern.h:332
 EdbPattern.h:333
 EdbPattern.h:334
 EdbPattern.h:335
 EdbPattern.h:336
 EdbPattern.h:337
 EdbPattern.h:338
 EdbPattern.h:339
 EdbPattern.h:340
 EdbPattern.h:341
 EdbPattern.h:342
 EdbPattern.h:343
 EdbPattern.h:344
 EdbPattern.h:345
 EdbPattern.h:346
 EdbPattern.h:347
 EdbPattern.h:348
 EdbPattern.h:349
 EdbPattern.h:350
 EdbPattern.h:351
 EdbPattern.h:352
 EdbPattern.h:353
 EdbPattern.h:354
 EdbPattern.h:355
 EdbPattern.h:356
 EdbPattern.h:357
 EdbPattern.h:358
 EdbPattern.h:359
 EdbPattern.h:360
 EdbPattern.h:361
 EdbPattern.h:362
 EdbPattern.h:363
 EdbPattern.h:364
 EdbPattern.h:365
 EdbPattern.h:366
 EdbPattern.h:367
 EdbPattern.h:368
 EdbPattern.h:369
 EdbPattern.h:370
 EdbPattern.h:371
 EdbPattern.h:372
 EdbPattern.h:373
 EdbPattern.h:374
 EdbPattern.h:375
 EdbPattern.h:376
 EdbPattern.h:377
 EdbPattern.h:378
 EdbPattern.h:379
 EdbPattern.h:380
 EdbPattern.h:381
 EdbPattern.h:382
 EdbPattern.h:383
 EdbPattern.h:384
 EdbPattern.h:385
 EdbPattern.h:386
 EdbPattern.h:387
 EdbPattern.h:388
 EdbPattern.h:389
 EdbPattern.h:390
 EdbPattern.h:391
 EdbPattern.h:392
 EdbPattern.h:393
 EdbPattern.h:394
 EdbPattern.h:395
 EdbPattern.h:396
 EdbPattern.h:397
 EdbPattern.h:398
 EdbPattern.h:399
 EdbPattern.h:400
 EdbPattern.h:401
 EdbPattern.h:402
 EdbPattern.h:403
 EdbPattern.h:404
 EdbPattern.h:405
 EdbPattern.h:406

» Last changed: 2017-07-03 09:59 » Last generated: 2017-07-03 09:59
This page has been automatically generated. For comments or suggestions regarding the documentation or ROOT in general please send a mail to ROOT support.