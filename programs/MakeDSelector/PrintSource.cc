#include "MakeDSelector.h"
#include <string>

inline bool IsRelevant(const string &name)
{
    return (name != "ComboBeam" && name != "Target" && name.find("Decaying") == string::npos);
}
static std::vector<std::string> tokenize(const std::string &s, char sep = '_')
{
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, sep))
    {
        if (!item.empty())
        {
            tokens.push_back(item);
        }
    }
    return tokens;
}

void Print_SourceFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap,
                      bool extraDefaults)
{
    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locSourceName = locSelectorName + string(".C");
    ofstream locSourceStream, csvOut("branches.csv");
    locSourceStream.open(locSourceName.c_str());

    locSourceStream << R"(#include ")" << locSelectorName << R"(.h"
#include <TLorentzRotation.h>

void )" << locSelectorName
                    << R"CODE(::Init(TTree *locTree)
{
   // USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

   // The Init() function is called when the selector needs to initialize a new tree or chain.
   // Typically here the branch addresses and branch pointers of the tree will be set.
   // Init() will be called many times when running on PROOF (once per file to be processed).

   //USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
   dOutputFileName = ""; 
   dOutputTreeFileName = ""; //"" for none
   dFlatTreeFileName = ")"
                    << locSelectorBaseName
                    << R"(.root"; //output flat tree (one combo per tree entry), "" for none
   dFlatTreeName = ""; //if blank, default name will be chosen
   //dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.
   //dSaveTLorentzVectorsAsFundamentaFlatTree = true; // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.

   //Because this function gets called for each TTree in the TChain, we must be careful:
   	//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
   bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
   DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
   //gDirectory now points to the output file with name dOutputFileName (if any)
   if(locInitializedPriorFlag)
   	return; //have already created histograms, etc. below: exit

   Get_ComboWrappers();
   dPreviousRunNumber = 0;

   /*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

   // EXAMPLE: Create deque for histogramming particle masses:
   // For histogramming the phi mass in phi -> K+ K-
   // Be sure to change this and dAnalyzeCutActions to match reaction
   std::deque<Particle_t> MyPhi;
   MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

   //ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
   //false/true below: use measured/kinfit data

   //PID
   dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
   //below: value: +/- N ns, UnknownParticle: All PIDs, SYS_NULL: all timing systems
   //dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

   //PIDFOM (for charged tracks)
   dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
   //dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
   //dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));

   //MASSES
   //dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
   //dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

   //KINFIT RESULTS
   dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

   //CUT MISSING MASS
   //dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

   //CUT ON SHOWER QUALITY
   //dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

   //BEAM ENERGY
   dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
   //dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999

   //KINEMATICS
   dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

   // ANALYZE CUT ACTIONS
   // // Change MyPhi to match reaction
   dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

   //INITIALIZE ACTIONS
   //If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
   Initialize_Actions();
   dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

   /******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

   //EXAMPLE MANUAL HISTOGRAMS:
   dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
   dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
   dHist_BeamEnergy_BestChiSq = new TH1I("BeamEnergy_BestChiSq", ";Beam Energy (GeV)", 600, 0.0, 12.0);

}
)CODE" << endl;

    set<string> branches;
    string finalState;
    if (extraDefaults)
    {
        locSourceStream << "    // == Extra default branches ==\n";
        // individual masses
        for (auto &step : locComboInfoMap)
        {
            for (auto &p : step.second)
            {
                string name = p.second.second;
                if (!IsRelevant(name))
                {
                    continue;
                }
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"mass_"
                                << name << "\");\n";
                // locSourceStream << "    "
                //                    "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"mass_"
                //                 << name << "_measured\");\n";
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                                   "(\"costh_lab_"
                                << name << "\");\n";
                // locSourceStream << "    "
                //                    "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                //                    "(\"costh_lab_"
                //                 << name << "_measured\");\n";
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<"
                                   "Double_t>(\"phi_lab_"
                                << name << "\");\n";
                // locSourceStream << "    "
                //                    "dFlatTreeInterface->Create_Branch_Fundamental<"
                //                    "Double_t>(\"phi_lab_"
                //                 << name << "_measured\");\n";
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<"
                                   "Double_t>(\"costh_GJ_"
                                << name << "\");\n";
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<"
                                   "Double_t>(\"phi_GJ_"
                                << name << "\");\n";
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                                   "(\"costh_H_"
                                << name << "\");\n";
                locSourceStream << "    "
                                   "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                                   "(\"phi_H_"
                                << name << "\");\n";
                csvOut << name << endl;
            }
        }

        // sorted combinations
        vector<string> raw;
        for (auto &step : locComboInfoMap)
        {
            for (auto &p : step.second)
            {

                const string &n = p.second.second;
                if (IsRelevant(n))
                {
                    raw.push_back(n);
                }
            }
        }

        // define desired order of particle types
        vector<string> order = {"Pi0",    "PiPlus", "PiMinus", "Photon",  "KLong",
                                "KShort", "KPlus",  "KMinus",  "Neutron", "Proton"};

        auto typeIndex = [&](const string &s) {
            for (size_t i = 0; i < order.size(); ++i)
            {
                if (s.rfind(order[i], 0) == 0) // starts with
                {
                    return int(i);
                }
            }
            return int(order.size());
        };

        auto numericSuffix = [&](const string &s) {
            size_t pos = s.find_last_not_of("0123456789");
            if (pos + 1 < s.size())
            {
                return stoi(s.substr(pos + 1));
            }
            return 0;
        };

        sort(raw.begin(), raw.end(), [&](const string &a, const string &b) {
            int ai = typeIndex(a), bi = typeIndex(b);
            if (ai != bi)
            {
                return ai < bi;
            }
            return numericSuffix(a) < numericSuffix(b);
        });

        // now generate combinations in that sorted order
        vector<string> names(raw.begin(), raw.end());
        for (size_t r = 2; r <= names.size(); ++r)
        {
            vector<int> idx(r);
            function<void(int, int)> comb = [&](int start, unsigned long int depth) {
                if (depth == r)
                {
                    string combo;
                    for (unsigned long int i = 0; i < r; ++i)
                    {
                        combo += names[idx[i]];
                        if (i + 1 < r)
                        {
                            combo += "_";
                        }
                    }
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<"
                                       "Double_t>(\"mass_"
                                    << combo << "\");\n";
                    //    "    "
                    //    "dFlatTreeInterface->Create_Branch_Fundamental<"
                    //    "Double_t>(\"mass_"
                    // << combo << "_measured\");\n";
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                                       "(\"costh_lab_"
                                    << combo << "\");\n";
                    // locSourceStream << "    "
                    //                    "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                    //                    "(\"costh_lab_"
                    //                 << combo << "_measured\");\n";
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<"
                                       "Double_t>(\"phi_lab_"
                                    << combo << "\");\n";
                    // locSourceStream << "    "
                    //                    "dFlatTreeInterface->Create_Branch_Fundamental<"
                    //                    "Double_t>(\"phi_lab_"
                    //                 << combo << "_measured\");\n";
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<"
                                       "Double_t>(\"costh_GJ_"
                                    << combo << "\");\n";
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<"
                                       "Double_t>(\"phi_GJ_"
                                    << combo << "\");\n";
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                                       "(\"costh_H_"
                                    << combo << "\");\n";
                    locSourceStream << "    "
                                       "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                                       "(\"phi_H_"
                                    << combo << "\");\n";

                    csvOut << combo << endl;
                    branches.insert(combo);

                    finalState = combo;

                    return;
                }
                for (int i = start; i < int(names.size()); ++i)
                {
                    idx[depth] = i;
                    comb(i + 1, depth + 1);
                }
            };
            comb(0, 0);
        }

        // fit stats
        // locSourceStream << "    "
        //                   "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
        //                   "(\"chisq\");\n";
        // locSourceStream << "    "
        //                   "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
        //                   "(\"clevel\");\n";
        // dE/dx for any charged
        // for (auto &step : locComboInfoMap) {
        //  for (auto &p : step.second) {
        //    if (ParticleCharge(p.second.first) != 0) {
        //      string name = p.second.second;
        //      if (!IsRelevant(name)) {
        //        continue;
        //      }
        //      locSourceStream << "    "
        //                         "dFlatTreeInterface->Create_Branch_Fundamental<"
        //                         "Double_t>(\"dEdx_cdc_"
        //                      << name << "\");\n";
        //      locSourceStream << "    "
        //                         "dFlatTreeInterface->Create_Branch_Fundamental<"
        //                         "Double_t>(\"dEdx_fdc_"
        //                      << name << "\");\n";
        //    }
        //  }
        locSourceStream << "    // == End extra defaults == \n";
    }
    // locSourceStream << "}\n\n";

    locSourceStream
        << R"(   /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

  //EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
  //The type for the branch must be included in the brackets
  //1st function argument is the name of the branch
  //2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
  /*
  dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
  dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array","my_int");
  dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
  dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
  dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
  */

  /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

  // RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
  // dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");

  //EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
  //The type for the branch must be included in the brackets
  //1st function argument is the name of the branch
  //2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
  /*
  dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
  dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
  dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
  dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
  */

  /************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

  //TO SAVE PROCESSING TIME
     //If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
     //By default, for each event, the data is retrieved for all branches
     //If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
     //Do this by doing something similar to the commented code below

  //dTreeInterface->Clear_GetEntryBranches(); //now get none
  //dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

  /************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

  dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

)"
        << "Bool_t " << locSelectorName << R"(::Process(Long64_t locEntry)
{
  // The Process() function is called for each entry in the tree. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  // Use fStatus to set the return value of TTree::Process().
  // The return value is currently not used.

  //CALL THIS FIRST
  DSelector::Process(locEntry); //Gets the data from the tree for the entry
  //cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
  //TLorentzVector locProductionX4 = Get_X4_Production();

  /******************************************** GET POLARIZATION ORIENTATION ******************************************/

  //Only if the run number changes
  //RCDB environment must be setup in order for this to work! (Will return false otherwise)
  UInt_t locRunNumber = Get_RunNumber();
  if(locRunNumber != dPreviousRunNumber)
  {
     dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
     dPreviousRunNumber = locRunNumber;
  }

  /********************************************* SETUP UNIQUENESS TRACKING ********************************************/

  //ANALYSIS ACTIONS: Reset uniqueness tracking for each action
  //For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
  Reset_Actions_NewEvent();
  dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

  //PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
     //Sometimes, some content is the exact same between one combo and the next
     	//e.g. maybe two combos have different beam particles, but the same data for the final-state
     //When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
     //So, for each quantity you histogram, keep track of what particles you used (for a given combo)
     //Then for each combo, just compare to what you used before, and make sure it's unique

  //EXAMPLE 0: Event-specific info:
  Bool_t locUsedSoFar_Event = false; // Flag used to mark if the best chi-squared combo is filled in the histogram

  //EXAMPLE 1: Particle-specific info:
  set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search. This container is used for the "hybrid" method dealing with combinatorics.

  //EXAMPLE 2: Combo-specific info:
     //In general: Could have multiple particles with the same PID: Use a set of Int_t's
     //In general: Multiple PIDs, so multiple sets: Contain within a map
     //Multiple combos: Contain maps within a set (easier, faster to search)
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

  //INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

  /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

  /*
  Int_t locMyInt = 7;
  dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

  TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
  dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

  for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
     dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
  */

  /************************************************* LOOP OVER COMBOS *************************************************/

  // Vector to store combo information
  std::vector<std::pair<UInt_t, Double_t>> loc_combos;

  // Pre-loop to gather kinfit ComboIndex-chiSq pairing and sort by chiSq value ascendingly
  for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
  {
     dComboWrapper->Set_ComboIndex(loc_i);
     Double_t locChiSq = dComboWrapper->Get_ChiSq_KinFit("");
     loc_combos.push_back(std::make_pair(loc_i, locChiSq));
  }
  // Sort the combos by ChiSq
  std::sort(loc_combos.begin(), loc_combos.end(), [](const std::pair<UInt_t, Double_t>& a, const std::pair<UInt_t, Double_t>& b) {
     return a.second < b.second;
  });

  //Loop over combos
  for(const auto& loc_combo : loc_combos)
  {
     UInt_t loc_i = loc_combo.first;
     //Set branch array indices for combo and all combo particles
     dComboWrapper->Set_ComboIndex(loc_i);

     // Is used to indicate when combos have been cut
     if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
     	continue; // Combo has been cut previously

     /********************************************** GET PARTICLE INDICES *********************************************/)"
        << endl
        << endl;

    // print particle indices
    map<int, map<int, pair<Particle_t, string>>>::iterator locStepIterator =
        locComboInfoMap.begin();
    for (; locStepIterator != locComboInfoMap.end(); ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locSourceStream << "   	//Step " << locStepIndex << endl;

        map<int, pair<Particle_t, string>> &locStepInfoMap = locStepIterator->second;
        map<int, pair<Particle_t, string>>::iterator locParticleIterator = locStepInfoMap.begin();
        for (; locParticleIterator != locStepInfoMap.end(); ++locParticleIterator)
        {
            Particle_t locPID = locParticleIterator->second.first;
            string locParticleName = locParticleIterator->second.second;

            if (locPID == UnknownParticle)
            {
                continue;
            }
            else if (locParticleName == "ComboBeam")
            {
                locSourceStream << "   	Int_t locBeamID = "
                                   "dComboBeamWrapper->Get_BeamID();"
                                << endl;
            }
            else if (locParticleName.substr(0, 6) == "Target")
            {
                continue;
            }
            else if (locParticleName.substr(0, 8) == "Decaying")
            {
                continue;
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                continue;
            }
            else if (ParticleCharge(locPID) != 0)
            {
                locSourceStream << "   	Int_t loc" << locParticleName << "TrackID = d"
                                << locParticleName << "Wrapper->Get_TrackID();" << endl;
            }
            else
            {

                locSourceStream << "   	Int_t loc" << locParticleName << "NeutralID = d"
                                << locParticleName

                                << "Wrapper->Get_NeutralID();" << endl;
            }
        }
        locSourceStream << endl;
    }

    locSourceStream
        << R"(/*********************************************** GET FOUR-MOMENTUM **********************************************/)"
        << endl
        << endl;

    // get p4's
    locSourceStream << "   	// Get P4\'s: //is kinfit if kinfit "
                       "performed, else is measured"
                    << endl;
    locSourceStream << "   	//dTargetP4 is target p4" << endl;

    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locSourceStream << "   	//Step " << locStepIndex << endl;

        map<int, pair<Particle_t, string>> &locStepInfoMap = locStepIterator->second;
        map<int, pair<Particle_t, string>>::iterator locParticleIterator = locStepInfoMap.begin();

        for (; locParticleIterator != locStepInfoMap.end(); ++locParticleIterator)
        {
            int locParticleIndex = locParticleIterator->first;
            Particle_t locPID = locParticleIterator->second.first;
            string locParticleName = locParticleIterator->second.second;

            if (locPID == UnknownParticle)
            {
                continue;
            }
            else if (locParticleName == "ComboBeam")
            {
                locSourceStream << "   	TLorentzVector locBeamP4 = "
                                   "dComboBeamWrapper->Get_P4();"
                                << endl;
            }
            else if (locParticleName.substr(0, 6) == "Target")
            {

                continue;
            }
            else if (locParticleName.substr(0, 8) == "Decaying")
            {
                string locBranchName = locParticleName + string("__P4_KinFit");

                if ((locTreeInterface->Get_Branch(locBranchName) != NULL) &&
                    (locParticleIndex < 0)) // else not reconstructed
                {
                    locSourceStream << "   	TLorentzVector loc" << locParticleName << "P4 = d"
                                    << locParticleName << "Wrapper->Get_P4();" << endl;
                }
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                string locBranchName = locParticleName + string("__P4_KinFit");

                if (locTreeInterface->Get_Branch(locBranchName) != NULL) // else not reconstructed
                {
                    locSourceStream << "   	TLorentzVector loc" << locParticleName << "P4 = d"
                                    << locParticleName << "Wrapper->Get_P4();" << endl;
                }
            }
            else // detected
            {
                locSourceStream << "   	TLorentzVector loc" << locParticleName << "P4 = d"
                                << locParticleName << "Wrapper->Get_P4();" << endl;
            }
        }
    }
    locSourceStream << endl;

    // get measured p4's
    locSourceStream << "   	// Get Measured P4\'s:" << endl;
    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locSourceStream << "   	//Step " << locStepIndex << endl;

        map<int, pair<Particle_t, string>> &locStepInfoMap = locStepIterator->second;
        map<int, pair<Particle_t, string>>::iterator locParticleIterator = locStepInfoMap.begin();
        for (; locParticleIterator != locStepInfoMap.end(); ++locParticleIterator)
        {
            Particle_t locPID = locParticleIterator->second.first;
            string locParticleName = locParticleIterator->second.second;

            if (locPID == UnknownParticle)
            {
                continue;
            }
            else if (locParticleName == "ComboBeam")
            {
                locSourceStream << "   	TLorentzVector locBeamP4_Measured = "
                                   "dComboBeamWrapper->Get_P4_Measured();"
                                << endl;
            }
            else if (locParticleName.substr(0, 6) == "Target")
            {
                continue;
            }
            else if (locParticleName.substr(0, 8) == "Decaying")
            {
                continue;
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                continue;
            }
            else // detected

            {
                locSourceStream << "   	TLorentzVector loc" << locParticleName << "P4_Measured = d"
                                << locParticleName << "Wrapper->Get_P4_Measured();" << endl;
            }
        }
    }

    locSourceStream
        << endl
        << R"(/********************************************* GET COMBO RF TIMING INFO *****************************************/

    TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
    // Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
    // Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
    // Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
    // Int_t locNumOutOfTimeBunchesInTree = XXX; //YOU need to specify this number
       //Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches)

    // Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
    // Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree;

    // Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
    // Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
    // Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
    // if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
    //    dComboWrapper->Set_IsComboCut(true);
    //    continue;
    // }

    /********************************************* COMBINE FOUR-MOMENTUM ********************************************/

    // DO YOUR STUFF HERE

    // Combine 4-vectors
    TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
    locMissingP4_Measured -= )"
        << endl;

    // calc missing p4
    bool locFirstFlag = true;
    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        map<int, pair<Particle_t, string>> &locStepInfoMap = locStepIterator->second;
        map<int, pair<Particle_t, string>>::iterator locParticleIterator = locStepInfoMap.begin();
        for (; locParticleIterator != locStepInfoMap.end(); ++locParticleIterator)
        {
            Particle_t locPID = locParticleIterator->second.first;
            string locParticleName = locParticleIterator->second.second;

            if (locPID == UnknownParticle)
            {
                continue;
            }
            else if (locParticleName == "ComboBeam")
            {
                continue;
            }
            else if (locParticleName.substr(0, 6) == "Target")
            {
                continue;
            }
            else if (locParticleName.substr(0, 8) == "Decaying")
            {
                continue;
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                continue;
            }

            // detected
            if (!locFirstFlag)
            {
                locSourceStream << " + ";
            }
            locFirstFlag = false;
            locSourceStream << "loc" << locParticleName << "P4_Measured";
        }
    }
    locSourceStream << ";" << endl << endl;

    if (extraDefaults)
    {
        locSourceStream << "TLorentzVector CM_P4 = locBeamP4 + dTargetP4;" << endl
                        << "TLorentzRotation CM_Boost(-CM_P4.BoostVector());" << endl;
        locSourceStream << "TLorentzVector beamCM = CM_Boost * locBeamP4;" << endl
                        << "TLorentzVector targetCM = CM_Boost * dTargetP4;" << endl;
        for (auto &branch : branches)
        {
            vector<string> particles = tokenize(branch);
            for (UInt_t i = 0; i < particles.size(); i++)
            {
                locSourceStream << endl
                                << "TLorentzVector particleXCM = CM_Boost * particleXP4_" << branch
                                << ";" << endl;
                locSourceStream << endl
                                << "TLorentzRotation restFrameXBoost_" << branch << "(-particleXCM_"
                                << branch << ".BoostVector());" << endl;
                locSourceStream << "TLorentzVector particleXP4_" << branch << " = restFrameXBoost_"
                                << branch << " * (";
                locSourceStream << "TLorentzVector particleXP4 = CM_Boost * (";
                if (i)
                {
                    locSourceStream << " + ";
                }
                locSourceStream << "loc" << particles[i] << "P4";
            }
            locSourceStream << ");" << endl;
        }

        for (auto &step : locComboInfoMap)
        {
            for (auto &p : step.second)
            {
                string name = p.second.second;
                if (!IsRelevant(name))
                {
                    continue;
                }
                locSourceStream << "TLorentzVector " << name << "CM = CM_Boost *"
                                << " loc" << name << "P4;" << endl;
            }
        }

        for (auto &branch : branches)
        {
            locSourceStream << "TLorentzVector referenceGJ_" << branch << " = restFrameXBoost_"
                            << branch << " * (";
            vector<string> particles = tokenize(branch);
            for (UInt_t i = 0; i < particles.size(); i++)
            {
                if (i)
                {
                    locSourceStream << " + ";
                }
                locSourceStream << particles[i] << "CM";
            }
            locSourceStream << ");" << endl;
            locSourceStream << endl
                            << "TLorentzVector beamGJ_" << branch << " = restFrameXBoost_" << branch
                            << " * beamCM;" << endl
                            << "TLorentzVector targetGJ_" << branch
                            << " = restFrameXBoost * targetCM;" << endl;
        }

        for (auto &step : locComboInfoMap)
        {
            for (auto &p : step.second)
            {
                string name = p.second.second;
                if (!IsRelevant(name))
                {
                    continue;
                }
                locSourceStream << "TLorentzVector " << name << "GJ = restFrameXBoost * " << name
                                << "CM;" << endl
                                << "TVector3 " << name << "P3 = " << name << "GJ.Vect();" << endl;
            }
        }

        for (auto &branch : branches)
        {
            string uv = branch;
            string lv = "";
            vector<string> fsparticles = tokenize(finalState);
            UInt_t counter = 0;
            for (auto &fsp : fsparticles)
            {

                if (uv.rfind(fsp) == string::npos)
                {
                    if (counter != 0)
                    {
                        lv.append(" + ");
                    }
                    counter++;
                    lv.append(fsp + "CM");
                }
            }
            if (lv != "")
            {
                locSourceStream << "TVector3 z_GJ_" << uv << " = (beamGJ.Vect()).Unit();" << endl
                                << "TVector3 y_GJ_" << uv << " = ((beamCM.Vect()).Cross(-(" << lv
                                << ").Vect())).Unit();" << endl
                                << "TVector3 x_GJ_" << uv << " = ((y_GJ_" << uv << ").Cross(z_GJ_"
                                << uv << ")).Unit();" << endl
                                << "double costh_GJ_" << uv << " = (referenceGJ_" << uv
                                << ".Vect()).Dot(z_GJ_" << uv << ") / (referenceGJ_" << uv
                                << ".Vect()).Mag();" << endl
                                << "double phi_GJ_" << uv << " = TMath::ATan2((referenceGJ_" << uv
                                << ".Vect()).Dot(y_GJ_" << uv << "), (referenceGJ_" << uv
                                << ".Vect()).Dot(x_GJ_" << uv << "));" << endl
                                << "TVector3 z_H_" << uv << " = particleXCM.Vect().Unit();" << endl
                                << "TVector3 y_H_" << uv << " = y_GJ_" << uv << ";" << endl
                                << "TVector3 x_H_" << uv << " = y_H_" << uv << ".Cross(z_H_" << uv
                                << ").Unit();" << endl
                                << "double costh_H_" << uv << " = referenceGJ_" << uv
                                << ".Vect().Dot(z_H_" << uv << ") / referenceGJ_" << uv
                                << ".Vect().Mag();" << endl
                                << "double phi_H_" << uv << " = TMath::ATan2(referenceGJ_" << uv
                                << ".Vect().Dot(y_H_" << uv << "), referenceGJ_" << uv
                                << ".Vect().Dot(x_H_" << uv << "));" << endl;
            }
        }
    }

    locSourceStream
        << R"(/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

    // Loop through the analysis actions, executing them in order for the active particle combo
    dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
    if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
       continue;

    //if you manually execute any actions, and it fails a cut, be sure to call:
       //dComboWrapper->Set_IsComboCut(true);

    /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

    /*
    TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
    //for arrays below: 2nd argument is value, 3rd is array index
    //NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
       //So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
    dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
    dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
    */

    /**************************************** EXAMPLE: BEST chi2 METHOD *****************************************/

    //Need to uncomment the section computing combo timing info before running this block of code
    //if(locUsedSoFar_Event == false)
    //{
       //Fill the histogram only when the beam bunch is in-time.
       //if(!locRelBeamBucket)
       //{
       //	dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());
       //	locUsedSoFar_Event = true;
       //}
    //}

    /**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

    //Histogram beam energy (if haven't already)
    if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
    {
       dHist_BeamEnergy->Fill(locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
       //dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); // Alternate version with accidental subtraction

       locUsedSoFar_BeamEnergy.insert(locBeamID);
    }

    /************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

    //Missing Mass Squared
    double locMissingMassSquared = locMissingP4_Measured.M2();

    //Uniqueness tracking: Build the map of particles used for the missing mass
       //For beam: Don't want to group with final-state photons. Instead use "UnknownParticle" PID (not ideal, but it's easy).
    map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
)" << endl;

    // insert uniqueness tracking
    // print particle indices
    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        map<int, pair<Particle_t, string>> &locStepInfoMap = locStepIterator->second;
        map<int, pair<Particle_t, string>>::iterator locParticleIterator = locStepInfoMap.begin();

        for (; locParticleIterator != locStepInfoMap.end(); ++locParticleIterator)
        {
            Particle_t locPID = locParticleIterator->second.first;
            string locParticleName = locParticleIterator->second.second;

            if (locPID == UnknownParticle)
            {
                continue;
            }
            else if (locParticleName == "ComboBeam")
            {
                locSourceStream << "   	"

                                   "locUsedThisCombo_MissingMass[UnknownParticle]."
                                   "insert(locBeamID); "
                                   "//beam"
                                << endl;
            }
            else if (locParticleName.substr(0, 6) == "Target")
            {
                continue;
            }
            else if (locParticleName.substr(0, 8) == "Decaying")
            {
                continue;
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                continue;
            }
            else if (ParticleCharge(locPID) != 0)
            {
                locSourceStream << "   	locUsedThisCombo_MissingMass[" << EnumString(locPID)
                                << "].insert(loc" << locParticleName << "TrackID);" << endl;
            }
            else
            {
                locSourceStream << "   	locUsedThisCombo_MissingMass[" << EnumString(locPID)
                                << "].insert(loc" << locParticleName << "NeutralID);" << endl;
            }
        }
    }
    locSourceStream << R"(

   	//compare to what's been used so far
   	if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
   	{
   		//unique missing mass combo: histogram it, and register this combo of particles
   		dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
   		//dHist_MissingMassSquared->Fill(locMissingMassSquared, locHistAccidWeightFactor); // Alternate version with accidental subtraction

   		locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
   	}

   	//E.g. Cut
   	//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
   	//{
   	//	dComboWrapper->Set_IsComboCut(true);
   	//	continue;
   	//}

   	/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

   	// RECOMMENDED: FILL ACCIDENTAL WEIGHT
   	// dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight", locHistAccidWeightFactor);)"
                    << endl;

    if (extraDefaults)
    {
        locSourceStream << "    // == Extra default fills ==\n";
        // individual masses
        for (auto &step : locComboInfoMap)
        {
            for (auto &p : step.second)
            {
                string name = p.second.second;
                if (!IsRelevant(name))
                {
                    continue;
                }
                locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                                << name << "\", loc" << name << "P4.M());\n";
                // locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                //                 << name << "_measured\", loc" << name << "P4_Measured.M());\n";
            }
        }
        // combinations
        //    // combinations (sorted by type and numeric suffix)
        vector<string> raw;
        for (auto &step : locComboInfoMap)
        {
            for (auto &p : step.second)
            {
                const string &n = p.second.second;
                if (IsRelevant(n))
                {
                    raw.push_back(n);
                }
            }
        }

        // desired particle order
        static const vector<string> order = {"Pi0",    "PiPlus", "PiMinus", "Photon",  "KLong",
                                             "KShort", "KPlus",  "KMinus",  "Neutron", "Proton"};
        auto typeIndex = [&](const string &s) {
            for (size_t i = 0; i < order.size(); ++i)
            {
                if (s.rfind(order[i], 0) == 0)
                {
                    return int(i);
                }
            }
            return int(order.size());
        };

        auto numericSuffix = [&](const string &s) {
            size_t pos = s.find_last_not_of("0123456789");
            if (pos + 1 < s.size())
            {
                return stoi(s.substr(pos + 1));
            }
            return 0;
        };

        sort(raw.begin(), raw.end(), [&](const string &a, const string &b) {
            int ia = typeIndex(a), ib = typeIndex(b);
            if (ia != ib)
            {
                return ia < ib;
            }
            return numericSuffix(a) < numericSuffix(b);
        });

        vector<string> names(raw.begin(), raw.end());

        // now generate all r‐body combinations
        for (size_t r = 2; r <= names.size(); ++r)
        {
            vector<int> idx(r);
            function<void(int, int)> comb = [&](int start, unsigned long int depth) {
                if (depth == r)
                {
                    // build the combo key, e.g. "PiPlus_PiMinus_Proton"
                    string combo;
                    for (unsigned long int i = 0; i < r; ++i)
                    {
                        combo += names[idx[i]];
                        if (i + 1 < r)
                        {
                            combo += "_";
                        }
                    }

                    // kinematic‐fit mass
                    locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                                    << combo << "\", (";

                    for (unsigned long int i = 0; i < r; ++i)
                    {
                        if (i)
                        {
                            locSourceStream << " + ";
                        }
                        locSourceStream << "loc" << names[idx[i]] << "P4";
                    }
                    locSourceStream << ").M());\n";

                    // // measured mass
                    // locSourceStream << " dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                    //                 << combo << "_measured\", (";
                    // for (unsigned long int i = 0; i < r; ++i)
                    // {
                    //     if (i)
                    //     {
                    //         locSourceStream << " + ";
                    //     }
                    //     locSourceStream << "loc" << names[idx[i]] << "P4_Measured";
                    // }
                    // locSourceStream << ").M());\n";

                    // costh_lab
                    locSourceStream
                        << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_lab_"
                        << combo << "\", (";

                    for (unsigned long int i = 0; i < r; ++i)
                    {
                        if (i)
                        {
                            locSourceStream << " + ";
                        }
                        locSourceStream << "loc" << names[idx[i]] << "P4";
                    }
                    locSourceStream << ").Vect().CosTheta());\n";

                    // measured costh_lab
                    // locSourceStream
                    //     << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_lab_"
                    //     << combo << "_measured\", (";
                    // for (unsigned long int i = 0; i < r; ++i)
                    // {
                    //     if (i)
                    //     {
                    //         locSourceStream << " + ";
                    //     }
                    //     locSourceStream << "loc" << names[idx[i]] << "P4_Measured";
                    // }
                    // locSourceStream << ").Vect().CosTheta());\n";

                    // phi_lab
                    locSourceStream
                        << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_lab_" << combo
                        << "\", (";

                    for (unsigned long int i = 0; i < r; ++i)
                    {
                        if (i)
                        {
                            locSourceStream << " + ";
                        }
                        locSourceStream << "loc" << names[idx[i]] << "P4";
                    }
                    locSourceStream << ").Vect().Phi());\n";

                    // measured phi_lab
                    // locSourceStream
                    //     << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_lab_" <<
                    //     combo
                    //     << "_measured\", (";
                    // for (unsigned long int i = 0; i < r; ++i)
                    // {
                    //     if (i)
                    //     {
                    //         locSourceStream << " + ";
                    //     }
                    //     locSourceStream << "loc" << names[idx[i]] << "P4_Measured";
                    // }
                    // locSourceStream << ").Vect().Phi());\n";
                    if (combo != finalState)
                    {
                        // costh_GJ_
                        locSourceStream
                            << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_GJ_"
                            << combo << "\", costh_GJ_" << combo << ");\n";

                        // phi_GJ
                        locSourceStream
                            << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_GJ_"
                            << combo << "\", phi_GJ_" << combo << ");\n";

                        // costh_H_
                        locSourceStream
                            << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_H_"
                            << combo << "\", costh_H_" << combo << ");\n";

                        // phi_H
                        locSourceStream
                            << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_H_"
                            << combo << "\", phi_H_" << combo << ");\n";
                    }
                    return;
                }
                for (int i = start; i < int(names.size()); ++i)
                {
                    idx[depth] = i;
                    comb(i + 1, depth + 1);
                }
            };
            comb(0, 0);
        }

        // fit stats
        // locSourceStream
        //    << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"chisq\", "
        //       "chisq/ndof);\n";
        // locSourceStream
        //    << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"clevel\", "
        //       "clevel);\n";
        // dE/dx for any charged
        // for (auto &step : locComboInfoMap) {
        //  for (auto &p : step.second) {
        //    if (ParticleCharge(p.second.first) != 0) {
        //      string name = p.second.second;
        //      if (!IsRelevant(name)) {
        //        continue;
        //      }
        //      locSourceStream
        //          << "    "
        //             "dFlatTreeInterface->Fill_Fundamental<Double_t>(\"dEdx_cdc_"
        //          << name << "\", dEdx_cdc_" << name << ");\n";
        //      locSourceStream
        //          << "    "
        //             "dFlatTreeInterface->Fill_Fundamental<Double_t>(\"dEdx_fdc_"
        //          << name << "\", dEdx_fdc_" << name << ");\n";
        //    }
        //  }
        //}
        // frames

        locSourceStream << "    // == End extra fills ==\n";
    }

    locSourceStream << R"(

     /*
     //FILL ANY CUSTOM BRANCHES FIRST!!
     Int_t locMyInt_Flat = 7;
     dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

     TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
     dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

     for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
     {
     	dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
     	TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
     	dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
     }
     */

     //FILL FLAT TREE
     Fill_FlatTree(); //for the active combo
  } // end of combo loop

  //FILL HISTOGRAMS: Num combos / events surviving actions
  Fill_NumCombosSurvivedHists();

  /******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
  /*
  //Thrown beam: just use directly
  if(dThrownBeam != NULL)
     double locEnergy = dThrownBeam->Get_P4().E();

  //Loop over throwns
  for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
  {
     //Set branch array indices corresponding to this particle
     dThrownWrapper->Set_ArrayIndex(loc_i);

     //Do stuff with the wrapper here ...
  }
  */

  /****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
  /*
  //Loop over beam particles (note, only those appearing in combos are present)
  for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
  {
     //Set branch array indices corresponding to this particle
     dBeamWrapper->Set_ArrayIndex(loc_i);

     //Do stuff with the wrapper here ...
  }

  //Loop over charged track hypotheses
  for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
  {
     //Set branch array indices corresponding to this particle
     dChargedHypoWrapper->Set_ArrayIndex(loc_i);

     //Do stuff with the wrapper here ...
  }

  //Loop over neutral particle hypotheses
  for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
  {
     //Set branch array indices corresponding to this particle
     dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

     //Do stuff with the wrapper here ...
  }
  */

  /************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
  /*
  Bool_t locIsEventCut = true;
  for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
     //Set branch array indices for combo and all combo particles
     dComboWrapper->Set_ComboIndex(loc_i);
     // Is used to indicate when combos have been cut
     if(dComboWrapper->Get_IsComboCut())
     	continue;
     locIsEventCut = false; // At least one combo succeeded
     break;
  }
  if(!locIsEventCut && dOutputTreeFileName != "")
     Fill_OutputTree();
  */

   return kTRUE;
}

void )" << locSelectorName
                    << R"(::Finalize(void)
{
   //Save anything to output here that you do not want to be in the default DSelector output ROOT file.

    //Otherwise, don't do anything else (especially if you are using PROOF).
   	//If you are using PROOF, this function is called on each thread,
   	//so anything you do will not have the combined information from the various threads.
   	//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

   //DO YOUR STUFF HERE

   //CALL THIS LAST
   DSelector::Finalize(); //Saves results to the output file
}
)" << endl;

    locSourceStream.close();
}
