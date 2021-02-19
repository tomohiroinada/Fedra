#ifndef ROOT_TOracleServerE2W
#define ROOT_TOracleServerE2W
 
#include "TOracleServerE2.h"
#include "EdbView.h"
#include "EdbBrick.h"

class TTree;
class EdbPattern;
class EdbPatternsVolume;

class TOracleServerE2W : public TOracleServerE2 {

 private:

  Int_t eNviewsPerArea;

 public:
  TOracleServerE2W(const char *db, const char *uid, const char *pw):
    TOracleServerE2(db, uid, pw){eNviewsPerArea=0;}
    ~TOracleServerE2W(){}

    Int_t       AddEventBricks(char *databrick);
    Int_t       AddBrick_Set(char *id, char *idrange_min, char *idrange_max, char *id_partition);
    Int_t       AddBrick_Space(char *id_brick, char *id_set);
    Int_t       AddPlate(char *id_eventbrick, char *dataplate);
    Int_t       AddPlateCalibration(char *id_eventbrick, char *id_process_operation, char *datacalibration);
    Int_t       AddZone(char *data);
    Int_t       AddProcessOperation(char *id_machine, char *id_programsettings, char *id_requester, 
				    char *id_parent_operation, char *id_eventbrick, char *id_plate, char *driverlevel,
				    char *id_calibration, char *starttime, char *finishtime, char *success, 
				    char *notes);
    Int_t       AddView(char *dataView);
    Int_t       AddBaseTracks(EdbPattern &pat, char *id_eventbrick, char *id_zone);
    Int_t       AddBaseTracks(TTree *tree, char *id_eventbrick, char *id_zone, bool usebuffer=true);
    Int_t       AddMicroTrack(char *datamicro);
    Int_t       AddBaseTrack(char *database);
    Int_t       AddScanbackPath(char *datapath);
    Int_t       AddScanbackPath(char *id_eventbrick, char *id_header_operation, int id_path, int id_start_plate, int skipCSconnection=0);
    Int_t       AddScanbackPrediction(char *dataprediciton);
    Int_t       AddTemplateMarkSets(char *datamarks);
    
    Int_t       AddPlateCalibration(char *id_eventbrick, char *id_process_operation, EdbPlateP *plate);
    Int_t       AddView(EdbView *view, int id_view, char *id_eventbrick, char *id_zone, bool usebuffer=true);
    Int_t       AddViews(EdbRun *run, char *id_eventbrick, char *id_zone, bool usebuffer=true);

    Int_t       AddVolume(char *id_eventbrick, char *id_process_operation, int ivolume);
    Int_t       AddVolumeSlice(char *datavolumeslice);
    Int_t       AddBSBpathsVolumes(char *databsbpathsvolumes);
    
    Int_t       DeleteBrick(char *id_eventbrick);
    Int_t       DeleteBrickSpace(char *id_brick);
    Int_t       DeleteOperation(char *id_brick, char *id_process_operation);
    Int_t       DeletePlateOperation(char *id_brick, char *id_process_operation, char *id_plate);

    Int_t       NviewsPerArea() {return eNviewsPerArea;};

    ClassDef(TOracleServerE2W,1)  // Write enabled access to the OPERA db 
};

#endif
