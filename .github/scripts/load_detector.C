

#ifdef __CLING__
R__LOAD_LIBRARY(libDDCore)
#endif



void load_detector(std::string compact){
  gSystem->Load("libDDCore");
  using namespace dd4hep;
  Detector& description = Detector::getInstance();
  description.fromCompact(compact);
}
