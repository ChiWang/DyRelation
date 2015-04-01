/*
 *  $Id: DmpAlgCalRelBinding.cc, 2015-03-03 18:36:32 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/09/2014
*/

#include <boost/python.hpp>
#include "DmpAlgCalibrationRel.h"

BOOST_PYTHON_MODULE(libDmpCalRel){
  using namespace boost::python;

  class_<DmpAlgCalibrationRel,boost::noncopyable,bases<DmpVAlg> >("DmpAlgCalibrationRel",init<>())
    ;
}

