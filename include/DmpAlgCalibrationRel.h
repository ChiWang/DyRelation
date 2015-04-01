/*
 *  $Id: DmpAlgCalibrationRel.h, 2015-03-03 23:17:03 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

#ifndef DmpAlgCalibrationRel_H
#define DmpAlgCalibrationRel_H

#include <map>
#include <fstream>
#include "DmpVAlg.h"
#include "DmpEvtPsdRaw.h"

class TH2D;

class DmpAlgCalibrationRel : public DmpVAlg{
/*
 *  DmpAlgCalibrationRel
 *
 */
public:
  DmpAlgCalibrationRel();
  ~DmpAlgCalibrationRel();

  bool Initialize();
  bool ProcessThisEvent();
  bool Finalize();

private:
  DmpEvtBgoRaw          *fEvtBgo;   // without pedetal, pure signal
  DmpEvtPsdRaw          *fEvtPsd;

  std::map<short, std::map<short,std::map<short,std::map<short,TH2D*> > > > fBgoRelHist;//[DmpParameterBgo::kPlaneNo*2][DmpParameterBgo::kBarNo][2][2];        // 14: layer, 22: bar, 2: side, 2: dy relation
  std::ofstream         o_RelData_Bgo;      //

  std::map<short,std::map<short,std::map<short,TH2D*> > > fPsdRelHist;//[DmpParameterPsd::kPlaneNo*2][DmpParameterPsd::kStripNo][2];
  std::ofstream         o_RelData_Psd;      //

};

#endif

