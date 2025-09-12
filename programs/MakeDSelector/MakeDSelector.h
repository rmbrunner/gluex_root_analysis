#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <string>

#include <TFile.h>
#include <TTree.h>

#include "../../libraries/DSelector/DTreeInterface.h"
#include "../../libraries/DSelector/particleType.h"

void Print_Usage(void);
void Print_HeaderFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap);
void Print_SourceFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap,
                      bool extraDefaults);
void Print_HeaderFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface);
void Print_SourceFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                            bool extraDefaults);

