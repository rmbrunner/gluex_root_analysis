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

    locSourceStream << "#include \"" << locSelectorName << ".h\"" << endl
                    << "#include <TLorentzRotation.h>" << endl
                    << endl
                    << "void " << locSelectorName << "::Init(TTree *locTree)" << endl
                    << "{" << endl
                    << "	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH "
                       "A \"USER\" OR "
                       "\"EXAMPLE\" LABEL. LEAVE THE REST ALONE."
                    << endl
                    << endl
                    << "	// The Init() function is called when the selector "
                       "needs to initialize "
                       "a new tree or chain."
                    << endl
                    << "	// Typically here the branch addresses and branch "
                       "pointers of the tree will be set."
                    << endl
                    << "	// Init() will be called many times when running on "
                       "PROOF (once per file to be processed)."
                    << endl
                    << endl
                    << "	//USERS: SET OUTPUT FILE NAME //can be overriden by "
                       "user in PROOF"
                    << endl
                    << "	dOutputFileName = \"\"; " << endl
                    << "	dOutputTreeFileName = \"\"; //\"\" for none" << endl
                    << "	dFlatTreeFileName = \"" << locSelectorBaseName
                    << ".root\"; //output flat tree (one combo per tree entry), \"\" for none"
                    << endl
                    << "	dFlatTreeName = \"\"; //if blank, default name will be chosen" << endl
                    << "	//dSaveDefaultFlatBranches = true; // False: don't "
                       "save default "
                       "branches, reduce disk footprint."
                    << endl
                    << "	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // "
                       "Default (or false): save "
                       "particles as TLorentzVector objects. True: save as four "
                       "doubles instead."
                    << endl
                    << endl
                    << "	//Because this function gets called for each TTree in "
                       "the TChain, we "
                       "must be careful:"
                    << endl
                    << "		//We need to re-initialize the tree interface "
                       "& branch wrappers, "
                       "but don't want to recreate histograms"
                    << endl
                    << "	bool locInitializedPriorFlag = dInitializedFlag; "
                       "//save whether have "
                       "been initialized previously"
                    << endl
                    << "	DSelector::Init(locTree); //This must be called to "
                       "initialize wrappers "
                       "for each new TTree"
                    << endl
                    << "	//gDirectory now points to the output file with name "
                       "dOutputFileName (if any)"
                    << endl
                    << "	if(locInitializedPriorFlag)" << endl
                    << "		return; //have already created histograms, "
                       "etc. below: exit"
                    << endl
                    << endl
                    << "	Get_ComboWrappers();" << endl
                    << "	dPreviousRunNumber = 0;" << endl
                    << endl
                    << "	/*********************************** EXAMPLE USER "
                       "INITIALIZATION: "
                       "ANALYSIS ACTIONS **********************************/"
                    << endl
                    << endl
                    << "	// EXAMPLE: Create deque for histogramming particle masses:" << endl
                    << "	// // For histogramming the phi mass in phi -> K+ K-" << endl
                    << "	// // Be sure to change this and dAnalyzeCutActions to "
                       "match reaction"
                    << endl
                    << "	std::deque<Particle_t> MyPhi;" << endl
                    << "	MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);" << endl
                    << endl
                    << "	//ANALYSIS ACTIONS: //Executed in order if added to "
                       "dAnalysisActions"
                    << endl
                    << "	//false/true below: use measured/kinfit data" << endl
                    << endl
                    << "	//PID" << endl
                    << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_ParticleID(dComboWrapper, false));"
                    << endl
                    << "	//below: value: +/- N ns, UnknownParticle: All PIDs, "
                       "SYS_NULL: all timing systems"
                    << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DCutAction_PIDDeltaT(dComboWrapper, "
                       "false, 0.5, KPlus, SYS_BCAL));"
                    << endl
                    << endl
                    << "	//PIDFOM (for charged tracks)" << endl
                    << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_PIDFOM(dComboWrapper));"
                    << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));"
                    << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DCutAction_EachPIDFOM(dComboWrapper, 0.1));"
                    << endl
                    << endl
                    << "	//MASSES" << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DHistogramAction_InvariantMass(dComboWrapper, "
                       "false, Lambda, 1000, 1.0, 1.2, \"Lambda\"));"
                    << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DHistogramAction_MissingMassSquared(dComboWrapper, "
                       "false, 1000, -0.1, 0.1));"
                    << endl
                    << endl
                    << "	//KINFIT RESULTS" << endl
                    << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_KinFitResults(dComboWrapper));"
                    << endl
                    << endl
                    << "	//CUT MISSING MASS" << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));"
                    << endl
                    << endl
                    << "	//CUT ON SHOWER QUALITY" << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));"
                    << endl
                    << endl
                    << "	//BEAM ENERGY" << endl
                    << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_BeamEnergy(dComboWrapper, false));"

                    << endl
                    << "	//dAnalysisActions.push_back(new "
                       "DCutAction_BeamEnergy(dComboWrapper, "
                       "false, 8.2, 8.8));  // Coherent peak for runs in the "
                       "range 30000-59999"
                    << endl
                    << endl
                    << "	//KINEMATICS" << endl
                    << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_ParticleComboKinematics(dComboWrapper, false));"
                    << endl
                    << endl
                    << "	// ANALYZE CUT ACTIONS" << endl
                    << "	// // Change MyPhi to match reaction" << endl
                    << "	dAnalyzeCutActions = new "
                       "DHistogramAction_AnalyzeCutActions( dAnalysisActions, "
                       "dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "
                       "\"CutActionEffect\" );"
                    << endl
                    << endl
                    << "	//INITIALIZE ACTIONS" << endl
                    << "	//If you create any actions that you want to run manually "
                       "(i.e. don't "
                       "add to dAnalysisActions), be sure to initialize them here as well"
                    << endl
                    << "	Initialize_Actions();" << endl
                    << "	dAnalyzeCutActions->Initialize(); // manual action, "
                       "must call Initialize()"
                    << endl
                    << endl
                    << "	/******************************** EXAMPLE USER INITIALIZATION: "
                       "STAND-ALONE HISTOGRAMS *******************************/"
                    << endl
                    << endl
                    << "	//EXAMPLE MANUAL HISTOGRAMS:" << endl
                    << "	dHist_MissingMassSquared = new "
                       "TH1I(\"MissingMassSquared\", \";Missing "
                       "Mass Squared (GeV/c^{2})^{2}\", 600, -0.06, 0.06);"
                    << endl
                    << "	dHist_BeamEnergy = new TH1I(\"BeamEnergy\", \";Beam "
                       "Energy (GeV)\", "
                       "600, 0.0, 12.0);"
                    << endl
                    << "	dHist_BeamEnergy_BestChiSq = new "
                       "TH1I(\"BeamEnergy_BestChiSq\", "
                       "\";Beam Energy (GeV)\", 600, 0.0, 12.0);"
                    << endl
                    << endl;

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

        // sorted combinations: collect, sort, and generate
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
        locSourceStream << "    // == End extra defaults ==\n";
    }
    // locSourceStream << "}\n\n";

    locSourceStream << "	/************************** EXAMPLE USER "
                       "INITIALIZATION: CUSTOM OUTPUT "
                       "BRANCHES - MAIN TREE *************************/"
                    << endl
                    << endl
                    << "	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT "
                       "ROOT FILE NAME MUST FIRST "
                       "BE GIVEN!!!! (ABOVE: TOP)):"
                    << endl
                    << "	//The type for the branch must be included in the brackets" << endl
                    << "	//1st function argument is the name of the branch" << endl
                    << "	//2nd function argument is the name of the "
                       "branch that contains the "
                       "size of the array (for fundamentals only)"
                    << endl
                    << "	/*" << endl
                    << "	dTreeInterface->Create_Branch_Fundamental<Int_t>(\"my_int\"); "
                       "//fundamental = char, int, float, double, etc."
                    << endl
                    << "	dTreeInterface->Create_Branch_FundamentalArray<"
                       "Int_t>(\"my_int_array\","
                       " \"my_int\");"
                    << endl
                    << "	dTreeInterface->Create_Branch_FundamentalArray<"
                       "Float_t>(\"my_combo_"
                       "array\", \"NumCombos\");"
                    << endl
                    << "	dTreeInterface->Create_Branch_NoSplitTObject<"
                       "TLorentzVector>(\"my_p4\");"
                    << endl
                    << "	dTreeInterface->Create_Branch_ClonesArray<"
                       "TLorentzVector>(\"my_p4_array\");"
                    << endl
                    << "	*/" << endl
                    << endl
                    << "	/************************** EXAMPLE USER "
                       "INITIALIZATION: CUSTOM OUTPUT "
                       "BRANCHES - FLAT TREE *************************/"
                    << endl
                    << endl
                    << "	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH"

                    << endl
                    << "	// "
                       "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                       "(\"accidweight\");"
                    << endl
                    << endl
                    << "	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT "

                       "ROOT FILE NAME MUST FIRST "
                       "BE GIVEN!!!! (ABOVE: TOP)):"
                    << endl
                    << "	//The type for the branch must be included in the brackets" << endl
                    << "	//1st function argument is the name of the branch" << endl
                    << "	//2nd function argument is the name of the "
                       "branch that contains the "
                       "size of the array (for fundamentals only)"
                    << endl
                    << "	/*" << endl
                    << "	dFlatTreeInterface->Create_Branch_Fundamental<"
                       "Int_t>(\"flat_my_int\"); "
                       "//fundamental = char, int, float, double, etc."
                    << endl
                    << "	dFlatTreeInterface->Create_Branch_"

                       "FundamentalArray<Int_t>(\"flat_my_"
                       "int_array\", \"flat_my_int\");"
                    << endl
                    << "	dFlatTreeInterface->Create_Branch_"
                       "NoSplitTObject<TLorentzVector>("
                       "\"flat_my_p4\");"
                    << endl
                    << "	dFlatTreeInterface->Create_Branch_ClonesArray<"
                       "TLorentzVector>(\"flat_"
                       "my_p4_array\");"
                    << endl
                    << "	*/" << endl

                    << endl
                    << "	/************************************* " << endl
                    << "ADVANCED EXAMPLE: CHOOSE "
                       "BRANCHES TO READ "
                    << endl
                    << "************************************/" << endl
                    << endl
                    << "	//TO SAVE PROCESSING TIME" << endl
                    << "		//If you know you don't need all of "
                       "the branches/data, but just a "
                       "subset of it, you can speed things up"
                    << endl
                    << "		//By default, for each event, the data "
                       "is retrieved for all branches"
                    << endl
                    << "		//If you know you only need data for "
                       "some branches, you can skip "
                       "grabbing data from the branches you don't need"
                    << endl
                    << "		//Do this by doing something similar "
                       "to the commented code below"
                    << endl
                    << endl
                    << "	//dTreeInterface->Clear_GetEntryBranches(); //now get none" << endl
                    << "	//dTreeInterface->Register_GetEntryBranch(\"Proton__P4\"); "
                       "//manually "
                       "set the branches you want"
                    << endl
                    << endl
                    << "	/**************************************" << endl
                    << " DETERMINE IF ANALYZING SIMULATED DATA " << endl
                    << "*************************************/" << endl
                    << endl
                    << "	dIsMC = (dTreeInterface->Get_Branch(\"MCWeight\") != NULL);" << endl
                    << endl
                    << "}" << endl
                    << endl
                    << "Bool_t " << locSelectorName << "::Process(Long64_t locEntry)" << endl
                    << "{" << endl
                    << "	// The Process() function is called for each "
                       "entry in the tree. "
                       "The entry argument"
                    << endl
                    << "	// specifies which entry in the currently "
                       "loaded tree is to be processed."
                    << endl
                    << "	//" << endl
                    << "	// This function should contain the \"body\" "
                       "of the analysis. It can contain"
                    << endl
                    << "	// simple or elaborate selection criteria, run "
                       "algorithms on the data"
                    << endl
                    << "	// of the event and typically fill histograms." << endl
                    << "	//" << endl
                    << "	// The processing can be stopped by calling Abort()." << endl
                    << "	// Use fStatus to set the return value of TTree::Process()." << endl
                    << "	// The return value is currently not used." << endl
                    << endl
                    << "	//CALL THIS FIRST" << endl
                    << "	DSelector::Process(locEntry); //Gets the data "
                       "from the tree for the entry"
                    << endl
                    << "	//cout << \"RUN \" << Get_RunNumber() << \", EVENT \" << "
                       "Get_EventNumber() << endl;"
                    << endl
                    << "	//TLorentzVector locProductionX4 = Get_X4_Production();" << endl
                    << endl
                    << "	/******************************************** " << endl
                    << "GET POLARIZATION ORIENTATION" << endl
                    << " ******************************************/" << endl
                    << endl
                    << "	//Only if the run number changes" << endl
                    << "	//RCDB environment must be setup in order for "
                       "this to work! (Will return false otherwise)"
                    << endl
                    << "	UInt_t locRunNumber = Get_RunNumber();" << endl
                    << "	if(locRunNumber != dPreviousRunNumber)" << endl
                    << "	{" << endl
                    << "		dIsPolarizedFlag = "
                       "dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);"
                    << endl
                    << "		dPreviousRunNumber = locRunNumber;" << endl
                    << "	}" << endl
                    << endl
                    << "	/********************************************* " << endl
                    << "SETUP UNIQUENESS TRACKING " << endl
                    << "********************************************/" << endl
                    << endl
                    << "	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action" << endl
                    << "	//For any actions that you are executing "
                       "manually, be sure to call "
                       "Reset_NewEvent() on them here"
                    << endl
                    << "	Reset_Actions_NewEvent();" << endl
                    << "	dAnalyzeCutActions->Reset_NewEvent(); // "
                       "manual action, must call "
                       "Reset_NewEvent()"
                    << endl
                    << endl
                    << "	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING" << endl
                    << "		//Sometimes, some content is the exact "
                       "same between one combo and the next"
                    << endl
                    << "			//e.g. maybe two combos have "
                       "different beam particles, but the "
                       "same data for the final-state"
                    << endl
                    << "		//When histogramming, you don\'t want to double-count "
                       "when this "
                       "happens: artificially inflates your signal (or background)"
                    << endl
                    << "		//So, for each quantity you histogram, "
                       "keep track of what "
                       "particles you used (for a given combo)"
                    << endl
                    << "		//Then for each combo, just compare to "
                       "what you used before, and "
                       "make sure it\'s unique"
                    << endl
                    << endl
                    << "	//EXAMPLE 0: Event-specific info:" << endl
                    << "	Bool_t locUsedSoFar_Event = false; // Flag "
                       "used to mark if the best "
                       "chi-squared combo is filled in the histogram"
                    << endl
                    << endl
                    << "	//EXAMPLE 1: Particle-specific info:" << endl
                    << "	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: "
                       "Unique ID for beam "
                       "particles. set: easy to use, fast to search. This "
                       "container is used for "
                       "the \"hybrid\" method dealing with combinatorics."
                    << endl
                    << endl
                    << "	//EXAMPLE 2: Combo-specific info:" << endl

                    << "		//In general: Could have multiple "
                       "particles with the same PID: Use "
                       "a set of Int_t\'s"
                    << endl
                    << "		//In general: Multiple PIDs, so "
                       "multiple sets: Contain within a map"

                    << endl
                    << "		//Multiple combos: Contain maps within "
                       "a set (easier, faster to search)"
                    << endl
                    << "	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;" << endl
                    << endl
                    << "	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE" << endl
                    << endl
                    << "	/**************************************** " << endl
                    << "EXAMPLE: FILL CUSTOM OUTPUT BRANCHES" << endl
                    << "**************************************/" << endl
                    << endl
                    << "	/*" << endl
                    << "	Int_t locMyInt = 7;" << endl
                    << "	dTreeInterface->Fill_Fundamental<Int_t>(\"my_int\", locMyInt);" << endl
                    << endl
                    << "	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);" << endl
                    << "	dTreeInterface->Fill_TObject<TLorentzVector>("
                       "\"my_p4\", locMyP4);"
                    << endl
                    << endl
                    << "	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)" << endl
                    << "		"
                       "dTreeInterface->Fill_Fundamental<Int_t>(\"my_int_array\", 3*loc_i, "
                       "loc_i); //2nd argument = value, 3rd = array index"
                    << endl
                    << "	*/" << endl
                    << endl
                    << "	/************************************************* " << endl
                    << " LOOP OVER COMBOS " << endl
                    << "*************************************************/" << endl
                    << endl
                    << "	// Vector to store combo information" << endl
                    << "	std::vector<std::pair<UInt_t, Double_t>> loc_combos;" << endl
                    << endl
                    << "	// Pre-loop to gather kinfit ComboIndex-chiSq "
                       "pairing and sort by "
                       "chiSq value ascendingly"
                    << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)" << endl
                    << "	{" << endl
                    << "		dComboWrapper->Set_ComboIndex(loc_i);" << endl
                    << "		Double_t locChiSq = "
                       "dComboWrapper->Get_ChiSq_KinFit(\"\");"
                    << endl
                    << "		loc_combos.push_back(std::make_pair(loc_i, locChiSq));" << endl
                    << "	}" << endl
                    << "	// Sort the combos by ChiSq" << endl

                    << "	std::sort(loc_combos.begin(), loc_combos.end(), [](const "
                       "std::pair<UInt_t, Double_t>& a, const std::pair<UInt_t, Double_t>& "
                       "b) {"
                    << endl

                    << "		return a.second < b.second;" << endl
                    << "	});" << endl
                    << endl
                    << "	//Loop over combos" << endl
                    << "	for(const auto& loc_combo : loc_combos)" << endl
                    << "	{" << endl
                    << "		UInt_t loc_i = loc_combo.first;" << endl
                    << "		//Set branch array indices for combo "
                       "and all combo particles"
                    << endl
                    << "		dComboWrapper->Set_ComboIndex(loc_i);" << endl
                    << endl
                    << "		// Is used to indicate when combos have been cut" << endl
                    << "		if(dComboWrapper->Get_IsComboCut()) // "
                       "Is false when tree "
                       "originally created"
                    << endl
                    << "			continue; // Combo has been cut previously" << endl
                    << endl
                    << "		/********************************************** " << endl
                    << "GET PARTICLE INDICES" << endl
                    << "*********************************************/" << endl
                    << endl
                    << "		//Used for tracking uniqueness when "
                       "filling histograms, and for "
                       "determining unused particles"
                    << endl
                    << endl;

    // print particle indices
    map<int, map<int, pair<Particle_t, string>>>::iterator locStepIterator =
        locComboInfoMap.begin();
    for (; locStepIterator != locComboInfoMap.end(); ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locSourceStream << "		//Step " << locStepIndex << endl;

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
                locSourceStream << "		Int_t locBeamID = "
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
                locSourceStream << "		Int_t loc" << locParticleName << "TrackID = d"
                                << locParticleName << "Wrapper->Get_TrackID();" << endl;
            }
            else
            {

                locSourceStream << "		Int_t loc" << locParticleName << "NeutralID = d"
                                << locParticleName

                                << "Wrapper->Get_NeutralID();" << endl;
            }
        }
        locSourceStream << endl;
    }
    locSourceStream << "		/*********************************************** GET "
                       "FOUR-MOMENTUM "
                       "**********************************************/"
                    << endl;
    locSourceStream << endl;

    // get p4's
    locSourceStream << "		// Get P4\'s: //is kinfit if kinfit "
                       "performed, else is measured"
                    << endl;
    locSourceStream << "		//dTargetP4 is target p4" << endl;
    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locSourceStream << "		//Step " << locStepIndex << endl;

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
                locSourceStream << "		TLorentzVector locBeamP4 = "
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
                    locSourceStream << "		TLorentzVector loc" << locParticleName << "P4 = d"
                                    << locParticleName << "Wrapper->Get_P4();" << endl;
                }
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                string locBranchName = locParticleName + string("__P4_KinFit");
                if (locTreeInterface->Get_Branch(locBranchName) != NULL) // else not reconstructed
                {
                    locSourceStream << "		TLorentzVector loc" << locParticleName << "P4 = d"
                                    << locParticleName << "Wrapper->Get_P4();" << endl;
                }
            }
            else // detected
            {
                locSourceStream << "		TLorentzVector loc" << locParticleName << "P4 = d"
                                << locParticleName << "Wrapper->Get_P4();" << endl;
            }
        }
    }
    locSourceStream << endl;

    // get measured p4's
    locSourceStream << "		// Get Measured P4\'s:" << endl;
    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locSourceStream << "		//Step " << locStepIndex << endl;

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
                locSourceStream << "		TLorentzVector locBeamP4_Measured = "
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
                locSourceStream << "		TLorentzVector loc" << locParticleName
                                << "P4_Measured = d" << locParticleName
                                << "Wrapper->Get_P4_Measured();" << endl;
            }
        }
    }

    locSourceStream << endl
                    << "		/********************************************* GET "
                       "COMBO RF TIMING "
                       "INFO *****************************************/"
                    << endl
                    << endl
                    << "		TLorentzVector locBeamX4_Measured = "
                       "dComboBeamWrapper->Get_X4_Measured();"
                    << endl
                    << "		// Double_t locBunchPeriod = "
                       "dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());"
                    << endl
                    << "		// Double_t locDeltaT_RF = "
                       "dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), "
                       "locBeamX4_Measured, dComboWrapper);"
                    << endl
                    << "		// Int_t locRelBeamBucket = "
                       "dAnalysisUtilities.Get_RelativeBeamBucket(Get_"
                       "RunNumber(), locBeamX4_Measured, "
                       "dComboWrapper); // 0 for in-time events, non-zero "
                       "integer for out-of-time photons"
                    << endl
                    << "		// Int_t locNumOutOfTimeBunchesInTree "
                       "= XXX; //YOU need to "
                       "specify this number"
                    << endl
                    << "			//Number of out-of-time beam "
                       "bunches in tree (on a single "
                       "side, so that total number out-of-time bunches "
                       "accepted is 2 times this "
                       "number for left + right bunches) "
                    << endl
                    << endl
                    << "		// Bool_t locSkipNearestOutOfTimeBunch = true; // "
                       "True: skip "
                       "events from nearest out-of-time bunch on either side (recommended)."
                    << endl
                    << "		// Int_t locNumOutOfTimeBunchesToUse = "
                       "locSkipNearestOutOfTimeBunch ? "
                       "locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; "

                    << endl
                    << "		// Double_t locAccidentalScalingFactor = "
                       "dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), "
                       "locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations "

                       "require "
                       "added factor, which is different for data and MC."
                    << endl
                    << "		// Double_t locAccidentalScalingFactorError = "
                       "dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber()"
                       ", "
                       "locBeamP4.E()); "
                       "// Ideal value would be 1, but deviations observed, need added "
                       "factor."
                    << endl
                    << "		// Double_t locHistAccidWeightFactor = "
                       "locRelBeamBucket==0 ? 1 : "
                       "-locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // "
                       "Weight by "
                       "1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time"
                    << endl
                    << "		// if(locSkipNearestOutOfTimeBunch && "
                       "abs(locRelBeamBucket)==1) { // Skip "
                       "nearest out-of-time bunch: tails of in-time "
                       "distribution also leak in"
                    << endl
                    << "		// 	dComboWrapper->Set_IsComboCut(true); " << endl
                    << "		// 	continue; " << endl
                    << "		// } " << endl

                    << endl
                    << "		/********************************************* COMBINE "
                       "FOUR-MOMENTUM ********************************************/"
                    << endl
                    << endl
                    << "		// DO YOUR STUFF HERE" << endl
                    << endl
                    << "		// Combine 4-vectors" << endl
                    << "		TLorentzVector locMissingP4_Measured = "
                       "locBeamP4_Measured + dTargetP4;"
                    << endl
                    << "		locMissingP4_Measured -= ";

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
    locSourceStream << ";" << endl;
    locSourceStream << endl;
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
                                << "double costh_GJ_" << uv
                                << " = "
                                   "(referenceGJ_"
                                << uv << ".Vect()).Dot(z_GJ_" << uv
                                << ") / "
                                   "(referenceGJ_"
                                << uv << ".Vect()).Mag();" << endl
                                << "double phi_GJ_" << uv
                                << " = "
                                   "TMath::ATan2((referenceGJ_"
                                << uv << ".Vect()).Dot(y_GJ_" << uv << "), (referenceGJ_" << uv
                                << ".Vect()).Dot(x_GJ_" << uv << "));" << endl
                                << "TVector3 z_H_" << uv << " = particleXCM.Vect().Unit();" << endl
                                << "TVector3 y_H_" << uv << " = y_GJ_" << uv << ";" << endl
                                << "TVector3 x_H_" << uv << " = y_H_" << uv << ".Cross(z_H_" << uv
                                << ").Unit();" << endl
                                << "double costh_H_" << uv << " = referenceGJ_" << uv
                                << ".Vect().Dot(z_H_" << uv
                                << ") / "
                                   "referenceGJ_"
                                << uv << ".Vect().Mag();" << endl
                                << "double phi_H_" << uv
                                << " = "
                                   "TMath::ATan2(referenceGJ_"
                                << uv << ".Vect().Dot(y_H_" << uv << "), referenceGJ_" << uv
                                << ".Vect().Dot(x_H_" << uv << "));" << endl;
            }
        }
    }
    locSourceStream << "		/******************************************** EXECUTE "
                       "ANALYSIS "
                       "ACTIONS *******************************************/"
                    << endl
                    << endl
                    << "		// Loop through the analysis actions, "
                       "executing them in order for "
                       "the active particle combo"
                    << endl
                    << "		dAnalyzeCutActions->Perform_Action(); "
                       "// Must be executed before "
                       "Execute_Actions()"
                    << endl
                    << "		if(!Execute_Actions()) //if the active "
                       "combo fails a cut, "
                       "IsComboCutFlag automatically set"
                    << endl
                    << "			continue;" << endl
                    << endl
                    << "		//if you manually execute any actions, "
                       "and it fails a cut, be sure to call:"
                    << endl
                    << "			//dComboWrapper->Set_IsComboCut(true);" << endl
                    << endl
                    << "		/**************************************** EXAMPLE: "
                       "FILL CUSTOM "
                       "OUTPUT BRANCHES **************************************/"
                    << endl
                    << endl
                    << "		/*" << endl
                    << "		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);" << endl
                    << "		//for arrays below: 2nd argument is "
                       "value, 3rd is array index"
                    << endl
                    << "		//NOTE: By filling here, AFTER the cuts above, some "
                       "indices won't "
                       "be updated (and will be whatever they were from the last event)"
                    << endl
                    << "			//So, when you draw the "
                       "branch, be sure to cut on "
                       "\"IsComboCut\" to avoid these."
                    << endl
                    << "		"
                       "dTreeInterface->Fill_Fundamental<Float_t>(\"my_combo_array\", "
                       "-2*loc_i, loc_i);"
                    << endl
                    << "		"
                       "dTreeInterface->Fill_TObject<TLorentzVector>(\"my_p4_array\", "
                       "locMyComboP4, loc_i);"
                    << endl
                    << "		*/" << endl
                    << endl
                    << "		/**************************************** EXAMPLE: "
                       "BEST chi2 "
                       "METHOD *****************************************/"
                    << endl
                    << endl
                    << "        // Need to uncomment the section computing "
                       "combo timing info "

                       "before running this block of code"
                    << endl
                    << "		//if(locUsedSoFar_Event == false)" << endl
                    << "		//{" << endl
                    << "			// Fill the histogram only "
                       "when the beam bunch is in-time. "
                    << endl
                    << "			//if(!locRelBeamBucket)" << endl
                    << "			//{" << endl
                    << "			//	"
                       "dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());"
                    << endl
                    << "			//	locUsedSoFar_Event = true;" << endl
                    << "			//}" << endl
                    << "		//}" << endl
                    << endl
                    << "		/**************************************** EXAMPLE: "
                       "HISTOGRAM BEAM "
                       "ENERGY *****************************************/"
                    << endl
                    << endl
                    << "		//Histogram beam energy (if haven\'t already)" << endl
                    << "		if(locUsedSoFar_BeamEnergy.find(locBeamID) == "
                       "locUsedSoFar_BeamEnergy.end())"
                    << endl
                    << "		{" << endl
                    << "			dHist_BeamEnergy->Fill(locBeamP4.E()); // "
                       "Fills in-time and "
                       "out-of-time beam photon combos"
                    << endl
                    << "			"
                       "//dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); "
                       "// "
                       "Alternate version with accidental subtraction"
                    << endl
                    << endl
                    << "			locUsedSoFar_BeamEnergy.insert(locBeamID);" << endl
                    << "		}" << endl
                    << endl
                    << "		/************************************ "
                       "EXAMPLE: HISTOGRAM MISSING "
                       "MASS SQUARED ************************************/"
                    << endl
                    << endl
                    << "		//Missing Mass Squared" << endl

                    << "		double locMissingMassSquared = "
                       "locMissingP4_Measured.M2();"
                    << endl
                    << endl
                    << "		//Uniqueness tracking: Build the map "
                       "of particles used for the missing mass"
                    << endl
                    << "			//For beam: Don\'t want to group with "
                       "final-state photons. "

                       "Instead use \"UnknownParticle\" PID (not ideal, but it\'s easy)."
                    << endl
                    << "		map<Particle_t, set<Int_t> > "
                       "locUsedThisCombo_MissingMass;"
                    << endl;

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
                locSourceStream << "		"

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
                locSourceStream << "		locUsedThisCombo_MissingMass[" << EnumString(locPID)
                                << "].insert(loc" << locParticleName << "TrackID);" << endl;
            }
            else
            {
                locSourceStream << "		locUsedThisCombo_MissingMass[" << EnumString(locPID)
                                << "].insert(loc" << locParticleName << "NeutralID);" << endl;
            }
        }
    }
    locSourceStream << endl

                    << "		//compare to what\'s been used so far" << endl
                    << "		"
                       "if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == "
                       "locUsedSoFar_MissingMass.end())"
                    << endl
                    << "		{" << endl
                    << "			//unique missing mass combo: "
                       "histogram it, and register this "
                       "combo of particles"
                    << endl
                    << "			"
                       "dHist_MissingMassSquared->Fill(locMissingMassSquared); // "
                       "Fills in-time and out-of-time beam photon combos"
                    << endl
                    << "			"
                       "//"
                       "dHist_MissingMassSquared->Fill(locMissingMassSquared,"
                       "locHistAccidWeightFactor); "
                       "// "
                       "Alternate version with accidental subtraction"
                    << endl
                    << endl
                    << "			"
                       "locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);"
                    << endl
                    << "		}" << endl
                    << endl
                    << "		//E.g. Cut" << endl
                    << "		//if((locMissingMassSquared < -0.04) "
                       "|| (locMissingMassSquared > 0.04))"
                    << endl
                    << "		//{" << endl
                    << "		//	dComboWrapper->Set_IsComboCut(true);" << endl
                    << "		//	continue;" << endl
                    << "		//}" << endl
                    << endl
                    << "		/****************************************** FILL FLAT "
                       "TREE (IF "
                       "DESIRED) ******************************************/"
                    << endl
                    << endl
                    << "		// RECOMMENDED: FILL ACCIDENTAL WEIGHT" << endl
                    << "		// "
                       "dFlatTreeInterface->Fill_Fundamental<Double_t>(\"accidweight\","
                       "locHistAccidWeightFactor);"
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

    locSourceStream << endl
                    << "		/*" << endl
                    << "		//FILL ANY CUSTOM BRANCHES FIRST!!" << endl

                    << "		Int_t locMyInt_Flat = 7;" << endl
                    << "		"
                       "dFlatTreeInterface->Fill_Fundamental<Int_t>(\"flat_my_int\", "
                       "locMyInt_Flat);"
                    << endl
                    << endl
                    << "		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);" << endl
                    << "		"
                       "dFlatTreeInterface->Fill_TObject<TLorentzVector>(\"flat_my_p4\", "
                       "locMyP4_Flat);"
                    << endl
                    << endl
                    << "		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)" << endl
                    << "		{" << endl
                    << "			"
                       "dFlatTreeInterface->Fill_Fundamental<Int_t>(\"flat_my_int_array\", "
                       "3*loc_j, loc_j); //2nd argument = value, 3rd = array index"
                    << endl
                    << "			TLorentzVector locMyComboP4_Flat(8.0, "
                       "7.0, 6.0, 5.0);"
                    << endl
                    << "			"
                       "dFlatTreeInterface->Fill_TObject<TLorentzVector>(\"flat_"
                       "my_p4_array\", "
                       "locMyComboP4_Flat, loc_j);"
                    << endl
                    << "		}" << endl
                    << "		*/" << endl
                    << endl
                    << "		//FILL FLAT TREE" << endl
                    << "		Fill_FlatTree(); //for the active combo" << endl
                    << "	} // end of combo loop" << endl
                    << endl
                    << "	//FILL HISTOGRAMS: Num combos / events surviving actions" << endl
                    << "	Fill_NumCombosSurvivedHists();" << endl
                    << endl
                    << "	/******************************************* LOOP OVER "
                       "THROWN DATA "
                       "(OPTIONAL) ***************************************/"
                    << endl
                    << "/*" << endl
                    << "	//Thrown beam: just use directly" << endl
                    << "	if(dThrownBeam != NULL)" << endl
                    << "		double locEnergy = dThrownBeam->Get_P4().E();" << endl
                    << endl
                    << "	//Loop over throwns" << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)" << endl
                    << "	{" << endl
                    << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl
                    << "		dThrownWrapper->Set_ArrayIndex(loc_i);" << endl
                    << endl
                    << "		//Do stuff with the wrapper here ..." << endl
                    << "	}" << endl
                    << "*/" << endl
                    << "	/****************************************** LOOP OVER "
                       "OTHER ARRAYS "
                       "(OPTIONAL) ***************************************/"
                    << endl
                    << "/*" << endl
                    << "	//Loop over beam particles (note, only those appearing in "
                       "combos are present)"
                    << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)" << endl
                    << "	{" << endl
                    << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl
                    << "		dBeamWrapper->Set_ArrayIndex(loc_i);" << endl
                    << endl
                    << "		//Do stuff with the wrapper here ..." << endl
                    << "	}" << endl
                    << endl
                    << "	//Loop over charged track hypotheses" << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)" << endl
                    << "	{" << endl
                    << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl
                    << "		dChargedHypoWrapper->Set_ArrayIndex(loc_i);" << endl
                    << endl
                    << "		//Do stuff with the wrapper here ..." << endl
                    << "	}" << endl
                    << endl
                    << "	//Loop over neutral particle hypotheses" << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)" << endl
                    << "	{" << endl
                    << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl
                    << "		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);" << endl
                    << endl
                    << "		//Do stuff with the wrapper here ..." << endl
                    << "	}" << endl
                    << "*/" << endl
                    << endl
                    << "	/************************************ EXAMPLE: FILL CLONE OF "
                       "TTREE "
                       "HERE WITH CUTS APPLIED ************************************/"
                    << endl
                    << "/*" << endl
                    << "	Bool_t locIsEventCut = true;" << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {" << endl
                    << "		//Set branch array indices for combo and all "
                       "combo particles"
                    << endl
                    << "		dComboWrapper->Set_ComboIndex(loc_i);" << endl
                    << "		// Is used to indicate when combos have been cut" << endl
                    << "		if(dComboWrapper->Get_IsComboCut())" << endl
                    << "			continue;" << endl
                    << "		locIsEventCut = false; // At least one combo succeeded" << endl
                    << "		break;" << endl
                    << "	}" << endl
                    << "	if(!locIsEventCut && dOutputTreeFileName != \"\")" << endl
                    << "		Fill_OutputTree();" << endl
                    << "*/" << endl
                    << endl
                    << "	return kTRUE;" << endl
                    << "}" << endl
                    << endl
                    << "void " << locSelectorName << "::Finalize(void)" << endl
                    << "{" << endl
                    << "	//Save anything to output here that you do not want to "
                       "be in the "
                       "default DSelector output ROOT file."
                    << endl
                    << endl
                    << "	//Otherwise, don\'t do anything else (especially if "
                       "you are using PROOF)."
                    << endl
                    << "		//If you are using PROOF, this function is "
                       "called on each thread,"
                    << endl
                    << "		//so anything you do will not have the "
                       "combined information from "
                       "the various threads."
                    << endl
                    << "		//Besides, it is best-practice to do "
                       "post-processing (e.g. "
                       "fitting) separately, in case there is a problem."
                    << endl
                    << endl
                    << "	//DO YOUR STUFF HERE" << endl
                    << endl
                    << "	//CALL THIS LAST" << endl
                    << "	DSelector::Finalize(); //Saves results to the output file" << endl
                    << "}" << endl;

    locSourceStream.close();
}
