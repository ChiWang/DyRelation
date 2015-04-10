
void Run()
{
  FileStat_t x;
  TString libName0 ="libDCRel"; // DAMPE Calibration Relation
  if(gSystem->GetPathInfo(libName0,x)){
    gSystem->Load("$DMPSWSYS/lib/libDmpBase.so");
    gInterpreter->AddIncludePath("$DMPSWSYS/include");
    gSystem->CompileMacro("./RelationStatus.C","k",libName0);
  }else{
    gSystem->Load(libName0);
  }
}

