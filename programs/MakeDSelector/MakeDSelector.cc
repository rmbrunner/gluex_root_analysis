#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <string>

#include <TFile.h>
#include <TTree.h>

#include "../../libraries/DSelector/DTreeInterface.h"
#include "../../libraries/DSelector/particleType.h"

// function declarations
void Print_Usage(void);
void Print_HeaderFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap);
void Print_SourceFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap,
                      bool extraDefaults);
void Print_HeaderFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface);
void Print_SourceFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                            bool extraDefaults);

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
int main(int argc, char *argv[])
{
    bool extraDefaults = false;
    vector<string> args;

    // parse --extra-defaults flag
    for (int i = 1; i < argc; ++i)
    {
        string s(argv[i]);
        if (s == "--extra-defaults")
        {
            extraDefaults = true;
        }
        else
        {
            args.push_back(s);
        }
    }
    if (args.size() != 3)
    {
        Print_Usage();
        return 0;
    }

    string locInputFileName = args[0];
    string locTreeName = args[1];
    string locSelectorBaseName = args[2];

    TFile *locInputFile = new TFile(locInputFileName.c_str(), "READ");
    TTree *locTree = (TTree *)locInputFile->Get(locTreeName.c_str());

    DTreeInterface *locTreeInterface = new DTreeInterface(locTree, true); // true: is input

    bool locIsMCGenTreeFlag = (locTreeInterface->Get_Branch("NumCombos") == NULL);
    if (locIsMCGenTreeFlag)
    {
        Print_HeaderFile_MCGen(locSelectorBaseName, locTreeInterface);
        Print_SourceFile_MCGen(locSelectorBaseName, locTreeInterface, extraDefaults);
        string locSelectorName = string("DSelector_") + locSelectorBaseName;
        cout << "Selector files " << locSelectorName << ".* generated." << endl;
        return 0;
    }

    // get combo info
    map<int, map<int, pair<Particle_t, string>>> locComboInfoMap;
    locTreeInterface->Get_ComboInfo(locComboInfoMap);

    Print_HeaderFile(locSelectorBaseName, locTreeInterface, locComboInfoMap);
    Print_SourceFile(locSelectorBaseName, locTreeInterface, locComboInfoMap, extraDefaults);

    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    cout << "Selector files " << locSelectorName << ".* generated." << endl;

    return 0;
}

void Print_Usage(void)
{
    cout << endl
         << "Makes a custom DSelector for the input TTree created by the DANA "
            "ANALYSIS library."
         << endl
         << "Usage: MakeDSelector [--extra-defaults] <input_root_file> "
            "<tree_name> "
            "<selector_base_name>"
         << endl
         << "    --extra-defaults   : auto-generate common flat-tree branches "
            "(1D mass plots, 2D correlation plots, Dalitz plots, angular distributions)"
         << endl
         << "  input_root_file     : ROOT file containing the TTree" << endl
         << "  tree_name           : name of TTree inside ROOT file" << endl
         << "  selector_base_name  : base name for generated DSelector files" << endl
         << endl;
}

void Print_HeaderFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap)
{
    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locHeaderName = locSelectorName + string(".h");
    ofstream locHeaderStream;
    locHeaderStream.open(locHeaderName.c_str());

    locHeaderStream << "#ifndef " << locSelectorName << "_h" << endl;
    locHeaderStream << "#define " << locSelectorName << "_h" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#include <iostream>" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#include \"DSelector/DSelector.h\"" << endl;
    locHeaderStream << "#include \"DSelector/DHistogramActions.h\"" << endl;
    locHeaderStream << "#include \"DSelector/DCutActions.h\"" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#include \"TH1I.h\"" << endl;
    locHeaderStream << "#include \"TH2I.h\"" << endl;
    locHeaderStream << endl;
    locHeaderStream << "class " << locSelectorName << " : public DSelector" << endl;
    locHeaderStream << "{" << endl;
    locHeaderStream << "	public:" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		" << locSelectorName
                    << "(TTree* locTree = NULL) : DSelector(locTree){}" << endl;
    locHeaderStream << "		virtual ~" << locSelectorName << "(){}" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		void Init(TTree *tree);" << endl;
    locHeaderStream << "		Bool_t Process(Long64_t entry);" << endl;
    locHeaderStream << endl;
    locHeaderStream << "	private:" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		void Get_ComboWrappers(void);" << endl;
    locHeaderStream << "		void Finalize(void);" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		// BEAM POLARIZATION INFORMATION" << endl;
    locHeaderStream << "		UInt_t dPreviousRunNumber;" << endl;
    locHeaderStream << "		bool dIsPolarizedFlag; //else is AMO" << endl;
    locHeaderStream << "		bool dIsPARAFlag; //else is PERP or AMO" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		bool dIsMC;" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		// ANALYZE CUT ACTIONS" << endl;
    locHeaderStream << "		// // Automatically makes mass histograms "
                       "where one cut is missing"
                    << endl;
    locHeaderStream << "		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS" << endl;
    locHeaderStream << endl;

    // print particle step, particle wrapper declarations
    map<int, map<int, pair<Particle_t, string>>>::iterator locStepIterator =
        locComboInfoMap.begin();
    for (; locStepIterator != locComboInfoMap.end(); ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locHeaderStream << "		//Step " << locStepIndex << endl;
        locHeaderStream << "		DParticleComboStep* dStep" << locStepIndex << "Wrapper;"
                        << endl;

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
                locHeaderStream << "		DBeamParticle* dComboBeamWrapper;" << endl;
            }
            else if (locParticleName.substr(0, 6) == "Target")
            {
                continue;
            }
            else if (locParticleName.substr(0, 8) == "Decaying")
            {
                string locBranchName = locParticleName + string("__P4_KinFit");
                if ((locTreeInterface->Get_Branch(locBranchName) != NULL) &&
                    (locParticleIndex < 0)) // else not reconstructed or in final state
                {
                    locHeaderStream << "		DKinematicData* d" << locParticleName << "Wrapper;"
                                    << endl;
                }
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                string locBranchName = locParticleName + string("__P4_KinFit");
                if (locTreeInterface->Get_Branch(locBranchName) != NULL) // else not reconstructed
                {
                    locHeaderStream << "		DKinematicData* d" << locParticleName << "Wrapper;"
                                    << endl;
                }
            }
            else if (ParticleCharge(locPID) != 0)
            {
                locHeaderStream << "		DChargedTrackHypothesis* d" << locParticleName
                                << "Wrapper;" << endl;
            }
            else
            {
                locHeaderStream << "		DNeutralParticleHypothesis* d" << locParticleName
                                << "Wrapper;" << endl;
            }
        }
        locHeaderStream << endl;
    }

    // resume
    locHeaderStream << "		// DEFINE YOUR HISTOGRAMS HERE" << endl;
    locHeaderStream << "		// EXAMPLES:" << endl;
    locHeaderStream << "		TH1I* dHist_MissingMassSquared;" << endl;
    locHeaderStream << "		TH1I* dHist_BeamEnergy;" << endl;
    locHeaderStream << "		TH1I* dHist_BeamEnergy_BestChiSq;" << endl;
    locHeaderStream << endl;
    locHeaderStream << "	ClassDef(" << locSelectorName << ", 0);" << endl;
    locHeaderStream << "};" << endl;
    locHeaderStream << endl;
    locHeaderStream << "void " << locSelectorName << "::Get_ComboWrappers(void)"

                    << endl;
    locHeaderStream << "{" << endl;

    // print particle step, particle wrapper assignments
    for (locStepIterator = locComboInfoMap.begin(); locStepIterator != locComboInfoMap.end();
         ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        if (locStepIndex != 0)
        {
            locHeaderStream << endl;
        }
        locHeaderStream << "	//Step " << locStepIndex << endl;
        locHeaderStream << "	dStep" << locStepIndex
                        << "Wrapper = dComboWrapper->Get_ParticleComboStep(" << locStepIndex << ");"
                        << endl;

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
                locHeaderStream << "	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep"
                                << locStepIndex << "Wrapper->Get_InitialParticle());" << endl;
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
                    locHeaderStream << "	d" << locParticleName << "Wrapper = dStep"
                                    << locStepIndex << "Wrapper->Get_InitialParticle();" << endl;
                }
            }
            else if (locParticleName.substr(0, 7) == "Missing")
            {
                string locBranchName = locParticleName + string("__P4_KinFit");
                if (locTreeInterface->Get_Branch(locBranchName) != NULL) // else not reconstructed
                {
                    locHeaderStream << "	d" << locParticleName << "Wrapper = dStep"
                                    << locStepIndex << "Wrapper->Get_FinalParticle("
                                    << locParticleIndex << ");" << endl;
                }
            }
            else if (ParticleCharge(locPID) != 0)
            {
                locHeaderStream << "	d" << locParticleName
                                << "Wrapper = static_cast<DChargedTrackHypothesis*>(dStep"
                                << locStepIndex << "Wrapper->Get_FinalParticle(" << locParticleIndex
                                << "));" << endl;
            }
            else
            {
                locHeaderStream << "	d" << locParticleName
                                << "Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep"
                                << locStepIndex << "Wrapper->Get_FinalParticle(" << locParticleIndex
                                << "));" << endl;
            }
        }
    }

    // resume
    locHeaderStream << "}" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#endif // " << locSelectorName << "_h" << endl;

    locHeaderStream.close();
}

void Print_SourceFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap,
                      bool extraDefaults)
{
    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locSourceName = locSelectorName + string(".C");
    ofstream locSourceStream, csvOut("branches.csv");
    locSourceStream.open(locSourceName.c_str());

    locSourceStream << "#include \"" << locSelectorName << ".h\"" << endl;
    locSourceStream << "#include <TLorentzRotation.h>" << endl;
    locSourceStream << endl;
    locSourceStream << "void " << locSelectorName << "::Init(TTree *locTree)" << endl;
    locSourceStream << "{" << endl;
    locSourceStream << "	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH "
                       "A \"USER\" OR "
                       "\"EXAMPLE\" LABEL. LEAVE THE REST ALONE."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	// The Init() function is called when the selector "
                       "needs to initialize "
                       "a new tree or chain."
                    << endl;
    locSourceStream << "	// Typically here the branch addresses and branch "
                       "pointers of the tree will be set."
                    << endl;
    locSourceStream << "	// Init() will be called many times when running on "
                       "PROOF (once per "
                       "file to be processed)."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//USERS: SET OUTPUT FILE NAME //can be overriden by "
                       "user in PROOF"
                    << endl;
    locSourceStream << "	dOutputFileName = \"\"; " << endl;
    locSourceStream << "	dOutputTreeFileName = \"\"; //\"\" for none" << endl;
    locSourceStream << "	dFlatTreeFileName = \"" << locSelectorBaseName
                    << ".root\"; //output flat tree (one "
                       "combo per tree "
                       "entry), \"\" for none"
                    << endl;
    locSourceStream << "	dFlatTreeName = \"\"; //if blank, default name will be chosen" << endl;
    locSourceStream << "	//dSaveDefaultFlatBranches = true; // False: don't "
                       "save default "
                       "branches, reduce disk footprint."
                    << endl;
    locSourceStream << "	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // "
                       "Default (or false): save "
                       "particles as TLorentzVector objects. True: save as four "
                       "doubles instead."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Because this function gets called for each TTree in "
                       "the TChain, we "
                       "must be careful:"
                    << endl;
    locSourceStream << "		//We need to re-initialize the tree interface "
                       "& branch wrappers, "
                       "but don't want to recreate histograms"
                    << endl;
    locSourceStream << "	bool locInitializedPriorFlag = dInitializedFlag; "
                       "//save whether have "
                       "been initialized previously"
                    << endl;
    locSourceStream << "	DSelector::Init(locTree); //This must be called to "
                       "initialize wrappers "
                       "for each new TTree"
                    << endl;
    locSourceStream << "	//gDirectory now points to the output file with name "
                       "dOutputFileName (if any)"
                    << endl;
    locSourceStream << "	if(locInitializedPriorFlag)" << endl;
    locSourceStream << "		return; //have already created histograms, "
                       "etc. below: exit"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	Get_ComboWrappers();" << endl;
    locSourceStream << "	dPreviousRunNumber = 0;" << endl;
    locSourceStream << endl;
    locSourceStream << "	/*********************************** EXAMPLE USER "
                       "INITIALIZATION: "
                       "ANALYSIS ACTIONS **********************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	// EXAMPLE: Create deque for histogramming particle masses:" << endl;
    locSourceStream << "	// // For histogramming the phi mass in phi -> K+ K-" << endl;
    locSourceStream << "	// // Be sure to change this and dAnalyzeCutActions to "
                       "match reaction"
                    << endl;
    locSourceStream << "	std::deque<Particle_t> MyPhi;" << endl;
    locSourceStream << "	MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);" << endl;
    locSourceStream << endl;
    locSourceStream << "	//ANALYSIS ACTIONS: //Executed in order if added to "
                       "dAnalysisActions"
                    << endl;
    locSourceStream << "	//false/true below: use measured/kinfit data" << endl;
    locSourceStream << endl;
    locSourceStream << "	//PID" << endl;
    locSourceStream << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_ParticleID(dComboWrapper, false));"
                    << endl;
    locSourceStream << "	//below: value: +/- N ns, UnknownParticle: All PIDs, "
                       "SYS_NULL: all timing systems"
                    << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DCutAction_PIDDeltaT(dComboWrapper, "
                       "false, 0.5, KPlus, SYS_BCAL));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//PIDFOM (for charged tracks)" << endl;
    locSourceStream << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_PIDFOM(dComboWrapper));"
                    << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));"
                    << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DCutAction_EachPIDFOM(dComboWrapper, 0.1));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//MASSES" << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DHistogramAction_InvariantMass(dComboWrapper, "
                       "false, Lambda, 1000, 1.0, 1.2, \"Lambda\"));"
                    << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DHistogramAction_MissingMassSquared(dComboWrapper, "
                       "false, 1000, -0.1, 0.1));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//KINFIT RESULTS" << endl;
    locSourceStream << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_KinFitResults(dComboWrapper));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//CUT MISSING MASS" << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//CUT ON SHOWER QUALITY" << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//BEAM ENERGY" << endl;
    locSourceStream << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_BeamEnergy(dComboWrapper, false));"

                    << endl;
    locSourceStream << "	//dAnalysisActions.push_back(new "
                       "DCutAction_BeamEnergy(dComboWrapper, "
                       "false, 8.2, 8.8));  // Coherent peak for runs in the "
                       "range 30000-59999"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//KINEMATICS" << endl;
    locSourceStream << "	dAnalysisActions.push_back(new "
                       "DHistogramAction_ParticleComboKinematics(dComboWrapper, false));"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	// ANALYZE CUT ACTIONS" << endl;
    locSourceStream << "	// // Change MyPhi to match reaction" << endl;
    locSourceStream << "	dAnalyzeCutActions = new "
                       "DHistogramAction_AnalyzeCutActions( dAnalysisActions, "
                       "dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "
                       "\"CutActionEffect\" );"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//INITIALIZE ACTIONS" << endl;
    locSourceStream << "	//If you create any actions that you want to run manually "
                       "(i.e. don't "
                       "add to dAnalysisActions), be sure to initialize them here as well"
                    << endl;
    locSourceStream << "	Initialize_Actions();" << endl;
    locSourceStream << "	dAnalyzeCutActions->Initialize(); // manual action, "
                       "must call Initialize()"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************** EXAMPLE USER INITIALIZATION: "
                       "STAND-ALONE HISTOGRAMS *******************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//EXAMPLE MANUAL HISTOGRAMS:" << endl;
    locSourceStream << "	dHist_MissingMassSquared = new "
                       "TH1I(\"MissingMassSquared\", \";Missing "
                       "Mass Squared (GeV/c^{2})^{2}\", 600, -0.06, 0.06);"
                    << endl;
    locSourceStream << "	dHist_BeamEnergy = new TH1I(\"BeamEnergy\", \";Beam "
                       "Energy (GeV)\", "
                       "600, 0.0, 12.0);"
                    << endl;
    locSourceStream << "	dHist_BeamEnergy_BestChiSq = new "
                       "TH1I(\"BeamEnergy_BestChiSq\", "
                       "\";Beam Energy (GeV)\", 600, 0.0, 12.0);"
                    << endl;
    locSourceStream << endl;

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
                    csvOut << combo << endl;
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
        string uv = "PiPlus_PiMinus1_Photon1_Photon2"; // TODO: currently must be user defined
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<"
                           "Double_t>(\"costh_GJ_"
                        << uv << "\");\n";
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<"
                           "Double_t>(\"phi_GJ_"
                        << uv << "\");\n";
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                           "(\"costh_H_"
                        << uv << "\");\n";
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                           "(\"phi_H_"
                        << uv << "\");\n";
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
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT "
                       "ROOT FILE NAME MUST FIRST "
                       "BE GIVEN!!!! (ABOVE: TOP)):"
                    << endl;
    locSourceStream << "	//The type for the branch must be included in the brackets" << endl;
    locSourceStream << "	//1st function argument is the name of the branch" << endl;
    locSourceStream << "	//2nd function argument is the name of the "
                       "branch that contains the "
                       "size of the array (for fundamentals only)"
                    << endl;
    locSourceStream << "	/*" << endl;
    locSourceStream << "	dTreeInterface->Create_Branch_Fundamental<Int_t>(\"my_int\"); "
                       "//fundamental = char, int, float, double, etc."
                    << endl;
    locSourceStream << "	dTreeInterface->Create_Branch_FundamentalArray<"
                       "Int_t>(\"my_int_array\","
                       " \"my_int\");"
                    << endl;
    locSourceStream << "	dTreeInterface->Create_Branch_FundamentalArray<"
                       "Float_t>(\"my_combo_"
                       "array\", \"NumCombos\");"
                    << endl;
    locSourceStream << "	dTreeInterface->Create_Branch_NoSplitTObject<"
                       "TLorentzVector>(\"my_p4\");"
                    << endl;
    locSourceStream << "	dTreeInterface->Create_Branch_ClonesArray<"
                       "TLorentzVector>(\"my_p4_array\");"
                    << endl;
    locSourceStream << "	*/" << endl;
    locSourceStream << endl;
    locSourceStream << "	/************************** EXAMPLE USER "
                       "INITIALIZATION: CUSTOM OUTPUT "
                       "BRANCHES - FLAT TREE *************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH"

                    << endl;
    locSourceStream << "	// "
                       "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                       "(\"accidweight\");"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT "

                       "ROOT FILE NAME MUST FIRST "
                       "BE GIVEN!!!! (ABOVE: TOP)):"
                    << endl;
    locSourceStream << "	//The type for the branch must be included in the brackets" << endl;
    locSourceStream << "	//1st function argument is the name of the branch" << endl;
    locSourceStream << "	//2nd function argument is the name of the "
                       "branch that contains the "
                       "size of the array (for fundamentals only)"
                    << endl;
    locSourceStream << "	/*" << endl;
    locSourceStream << "	dFlatTreeInterface->Create_Branch_Fundamental<"
                       "Int_t>(\"flat_my_int\"); "
                       "//fundamental = char, int, float, double, etc."
                    << endl;
    locSourceStream << "	dFlatTreeInterface->Create_Branch_"

                       "FundamentalArray<Int_t>(\"flat_my_"
                       "int_array\", \"flat_my_int\");"
                    << endl;
    locSourceStream << "	dFlatTreeInterface->Create_Branch_"
                       "NoSplitTObject<TLorentzVector>("
                       "\"flat_my_p4\");"
                    << endl;
    locSourceStream << "	dFlatTreeInterface->Create_Branch_ClonesArray<"
                       "TLorentzVector>(\"flat_"
                       "my_p4_array\");"
                    << endl;
    locSourceStream << "	*/" << endl;

    locSourceStream << endl;
    locSourceStream << "	/************************************* "
                       "ADVANCED EXAMPLE: CHOOSE "
                       "BRANCHES TO READ ************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//TO SAVE PROCESSING TIME" << endl;
    locSourceStream << "		//If you know you don't need all of "
                       "the branches/data, but just a "
                       "subset of it, you can speed things up"
                    << endl;
    locSourceStream << "		//By default, for each event, the data "
                       "is retrieved for all branches"
                    << endl;
    locSourceStream << "		//If you know you only need data for "
                       "some branches, you can skip "
                       "grabbing data from the branches you don't need"
                    << endl;
    locSourceStream << "		//Do this by doing something similar "
                       "to the commented code below"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//dTreeInterface->Clear_GetEntryBranches(); //now get none" << endl;
    locSourceStream << "	//dTreeInterface->Register_GetEntryBranch(\"Proton__P4\"); "
                       "//manually "
                       "set the branches you want"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	/************************************** DETERMINE IF ANALYZING "
                       "SIMULATED DATA *************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	dIsMC = (dTreeInterface->Get_Branch(\"MCWeight\") != NULL);" << endl;
    locSourceStream << endl;
    locSourceStream << "}" << endl;
    locSourceStream << endl;
    locSourceStream << "Bool_t " << locSelectorName << "::Process(Long64_t locEntry)" << endl;
    locSourceStream << "{" << endl;
    locSourceStream << "	// The Process() function is called for each "
                       "entry in the tree. "
                       "The entry argument"
                    << endl;
    locSourceStream << "	// specifies which entry in the currently "
                       "loaded tree is to be processed."
                    << endl;
    locSourceStream << "	//" << endl;
    locSourceStream << "	// This function should contain the \"body\" "
                       "of the analysis. It can contain"
                    << endl;
    locSourceStream << "	// simple or elaborate selection criteria, run "
                       "algorithms on the data"
                    << endl;
    locSourceStream << "	// of the event and typically fill histograms." << endl;
    locSourceStream << "	//" << endl;
    locSourceStream << "	// The processing can be stopped by calling Abort()." << endl;
    locSourceStream << "	// Use fStatus to set the return value of TTree::Process()." << endl;
    locSourceStream << "	// The return value is currently not used." << endl;
    locSourceStream << endl;
    locSourceStream << "	//CALL THIS FIRST" << endl;
    locSourceStream << "	DSelector::Process(locEntry); //Gets the data "
                       "from the tree for the entry"
                    << endl;
    locSourceStream << "	//cout << \"RUN \" << Get_RunNumber() << \", EVENT \" << "
                       "Get_EventNumber() << endl;"
                    << endl;
    locSourceStream << "	//TLorentzVector locProductionX4 = Get_X4_Production();" << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************************** GET POLARIZATION "
                       "ORIENTATION ******************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Only if the run number changes" << endl;
    locSourceStream << "	//RCDB environment must be setup in order for "
                       "this to work! (Will "
                       "return false otherwise)"
                    << endl;
    locSourceStream << "	UInt_t locRunNumber = Get_RunNumber();" << endl;
    locSourceStream << "	if(locRunNumber != dPreviousRunNumber)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		dIsPolarizedFlag = "
                       "dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);"
                    << endl;
    locSourceStream << "		dPreviousRunNumber = locRunNumber;" << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << endl;
    locSourceStream << "	/********************************************* "
                       "SETUP UNIQUENESS "
                       "TRACKING ********************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action" << endl;
    locSourceStream << "	//For any actions that you are executing "
                       "manually, be sure to call "
                       "Reset_NewEvent() on them here"
                    << endl;
    locSourceStream << "	Reset_Actions_NewEvent();" << endl;
    locSourceStream << "	dAnalyzeCutActions->Reset_NewEvent(); // "
                       "manual action, must call "
                       "Reset_NewEvent()"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING" << endl;
    locSourceStream << "		//Sometimes, some content is the exact "
                       "same between one combo and the next"
                    << endl;
    locSourceStream << "			//e.g. maybe two combos have "
                       "different beam particles, but the "
                       "same data for the final-state"
                    << endl;
    locSourceStream << "		//When histogramming, you don\'t want to double-count "
                       "when this "
                       "happens: artificially inflates your signal (or background)"
                    << endl;
    locSourceStream << "		//So, for each quantity you histogram, "
                       "keep track of what "
                       "particles you used (for a given combo)"
                    << endl;
    locSourceStream << "		//Then for each combo, just compare to "
                       "what you used before, and "
                       "make sure it\'s unique"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//EXAMPLE 0: Event-specific info:" << endl;
    locSourceStream << "	Bool_t locUsedSoFar_Event = false; // Flag "
                       "used to mark if the best "
                       "chi-squared combo is filled in the histogram"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//EXAMPLE 1: Particle-specific info:" << endl;
    locSourceStream << "	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: "
                       "Unique ID for beam "
                       "particles. set: easy to use, fast to search. This "
                       "container is used for "
                       "the \"hybrid\" method dealing with combinatorics."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//EXAMPLE 2: Combo-specific info:" << endl;

    locSourceStream << "		//In general: Could have multiple "
                       "particles with the same PID: Use "
                       "a set of Int_t\'s"
                    << endl;
    locSourceStream << "		//In general: Multiple PIDs, so "
                       "multiple sets: Contain within a map"

                    << endl;
    locSourceStream << "		//Multiple combos: Contain maps within "
                       "a set (easier, faster to search)"
                    << endl;
    locSourceStream << "	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;" << endl;
    locSourceStream << endl;
    locSourceStream << "	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE" << endl;
    locSourceStream << endl;
    locSourceStream << "	/**************************************** "
                       "EXAMPLE: FILL CUSTOM OUTPUT "
                       "BRANCHES **************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	/*" << endl;
    locSourceStream << "	Int_t locMyInt = 7;" << endl;
    locSourceStream << "	dTreeInterface->Fill_Fundamental<Int_t>(\"my_int\", locMyInt);" << endl;
    locSourceStream << endl;
    locSourceStream << "	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);" << endl;
    locSourceStream << "	dTreeInterface->Fill_TObject<TLorentzVector>("
                       "\"my_p4\", locMyP4);"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)" << endl;
    locSourceStream << "		"
                       "dTreeInterface->Fill_Fundamental<Int_t>(\"my_int_array\", 3*loc_i, "
                       "loc_i); //2nd argument = value, 3rd = array index"
                    << endl;
    locSourceStream << "	*/" << endl;
    locSourceStream << endl;
    locSourceStream << "	/************************************************* LOOP OVER "
                       "COMBOS "
                       "*************************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	// Vector to store combo information" << endl;
    locSourceStream << "	std::vector<std::pair<UInt_t, Double_t>> loc_combos;" << endl;
    locSourceStream << endl;
    locSourceStream << "	// Pre-loop to gather kinfit ComboIndex-chiSq "
                       "pairing and sort by "
                       "chiSq value ascendingly"
                    << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		dComboWrapper->Set_ComboIndex(loc_i);" << endl;
    locSourceStream << "		Double_t locChiSq = "
                       "dComboWrapper->Get_ChiSq_KinFit(\"\");"
                    << endl;
    locSourceStream << "		loc_combos.push_back(std::make_pair(loc_i, locChiSq));" << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << "	// Sort the combos by ChiSq" << endl;
    locSourceStream

        << "	std::sort(loc_combos.begin(), loc_combos.end(), [](const "
           "std::pair<UInt_t, Double_t>& a, const std::pair<UInt_t, Double_t>& "
           "b) {"
        << endl;

    locSourceStream << "		return a.second < b.second;" << endl;
    locSourceStream << "	});" << endl;
    locSourceStream << endl;
    locSourceStream << "	//Loop over combos" << endl;
    locSourceStream << "	for(const auto& loc_combo : loc_combos)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		UInt_t loc_i = loc_combo.first;" << endl;
    locSourceStream << "		//Set branch array indices for combo "
                       "and all combo particles"
                    << endl;
    locSourceStream << "		dComboWrapper->Set_ComboIndex(loc_i);" << endl;
    locSourceStream << endl;
    locSourceStream << "		// Is used to indicate when combos have been cut" << endl;
    locSourceStream << "		if(dComboWrapper->Get_IsComboCut()) // "
                       "Is false when tree "
                       "originally created"
                    << endl;
    locSourceStream << "			continue; // Combo has been cut previously" << endl;
    locSourceStream << endl;
    locSourceStream << "		/********************************************** GET "
                       "PARTICLE "
                       "INDICES *********************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		//Used for tracking uniqueness when "
                       "filling histograms, and for "
                       "determining unused particles"
                    << endl;
    locSourceStream << endl;

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

    locSourceStream << endl;
    locSourceStream << "		/********************************************* GET "
                       "COMBO RF TIMING "
                       "INFO *****************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		TLorentzVector locBeamX4_Measured = "
                       "dComboBeamWrapper->Get_X4_Measured();"
                    << endl;
    locSourceStream << "		// Double_t locBunchPeriod = "
                       "dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());"
                    << endl;
    locSourceStream << "		// Double_t locDeltaT_RF = "
                       "dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), "
                       "locBeamX4_Measured, dComboWrapper);"
                    << endl;
    locSourceStream << "		// Int_t locRelBeamBucket = "
                       "dAnalysisUtilities.Get_RelativeBeamBucket(Get_"
                       "RunNumber(), locBeamX4_Measured, "
                       "dComboWrapper); // 0 for in-time events, non-zero "
                       "integer for out-of-time photons"
                    << endl;
    locSourceStream << "		// Int_t locNumOutOfTimeBunchesInTree "
                       "= XXX; //YOU need to "
                       "specify this number"
                    << endl;
    locSourceStream << "			//Number of out-of-time beam "
                       "bunches in tree (on a single "
                       "side, so that total number out-of-time bunches "
                       "accepted is 2 times this "
                       "number for left + right bunches) "
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		// Bool_t locSkipNearestOutOfTimeBunch = true; // "
                       "True: skip "
                       "events from nearest out-of-time bunch on either side (recommended)."
                    << endl;
    locSourceStream << "		// Int_t locNumOutOfTimeBunchesToUse = "
                       "locSkipNearestOutOfTimeBunch ? "
                       "locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; "

                    << endl;
    locSourceStream << "		// Double_t locAccidentalScalingFactor = "
                       "dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), "
                       "locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations "

                       "require "
                       "added factor, which is different for data and MC."
                    << endl;
    locSourceStream << "		// Double_t locAccidentalScalingFactorError = "
                       "dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber()"
                       ", "
                       "locBeamP4.E()); "
                       "// Ideal value would be 1, but deviations observed, need added "
                       "factor."
                    << endl;
    locSourceStream << "		// Double_t locHistAccidWeightFactor = "
                       "locRelBeamBucket==0 ? 1 : "
                       "-locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // "
                       "Weight by "
                       "1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time"
                    << endl;
    locSourceStream << "		// if(locSkipNearestOutOfTimeBunch && "
                       "abs(locRelBeamBucket)==1) { // Skip "
                       "nearest out-of-time bunch: tails of in-time "
                       "distribution also leak in"
                    << endl;
    locSourceStream << "		// 	dComboWrapper->Set_IsComboCut(true); " << endl;
    locSourceStream << "		// 	continue; " << endl;
    locSourceStream << "		// } " << endl;

    locSourceStream << endl;
    locSourceStream << "		/********************************************* COMBINE "
                       "FOUR-MOMENTUM ********************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		// DO YOUR STUFF HERE" << endl;
    locSourceStream << endl;
    locSourceStream << "		// Combine 4-vectors" << endl;
    locSourceStream << "		TLorentzVector locMissingP4_Measured = "
                       "locBeamP4_Measured + dTargetP4;"
                    << endl;
    locSourceStream << "		locMissingP4_Measured -= ";

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
                        << "TLorentzRotation CM_Boost(-CM_P4.BoostVector());" << endl
                        << "TLorentzVector particleXP4 = CM_Boost * (locPiPlusP4 + locPiMinus1P4 "
                           "+ locPhoton1P4 + locPhoton1P4); // For now must be user defined!!!"
                        << endl;
        locSourceStream << "TLorentzVector beamCM = CM_Boost * locBeamP4;" << endl
                        << "TLorentzVector targetCM = CM_Boost * dTargetP4;" << endl;

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

        locSourceStream << endl
                        << "TLorentzVector particleXCM = CM_Boost * particleXP4;" << endl
                        << endl
                        << "TLorentzRotation restFrameXBoost(-particleXCM.BoostVector());" << endl;

        locSourceStream << "TLorentzVector referenceGJ = restFrameXBoost * (PiPlusCM + PiMinus1CM "
                           "+ Photon1CM + Photon2CM); // For now must "
                           "be user defined!!!"
                        << endl;
        locSourceStream << "TLorentzVector beamGJ = restFrameXBoost * beamCM;" << endl
                        << "TLorentzVector targetGJ = restFrameXBoost * targetCM;" << endl;
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

        locSourceStream
            << "TVector3 z_GJ = (beamGJ.Vect()).Unit();" << endl
            << "TVector3 y_GJ = ((beamCM.Vect()).Cross(-ProtonCM.Vect())).Unit();" << endl
            << "TVector3 x_GJ = ((y_GJ).Cross(z_GJ)).Unit();" << endl
            << "double costh_GJ_PiPlus_PiMinus1_Photon1_Photon2 = (referenceGJ.Vect()).Dot(z_GJ) / "
               "(referenceGJ.Vect()).Mag();"
            << endl
            << "double phi_GJ_PiPlus_PiMinus1_Photon1_Photon2 = "
               "TMath::ATan2((referenceGJ.Vect()).Dot(y_GJ), (referenceGJ.Vect()).Dot(x_GJ));"
            << endl
            << "TVector3 z_H = particleXCM.Vect().Unit();" << endl
            << "TVector3 y_H = y_GJ;" << endl
            << "TVector3 x_H = y_H.Cross(z_H).Unit();" << endl
            << "double costh_H_PiPlus_PiMinus1_Photon1_Photon2 = referenceGJ.Vect().Dot(z_H) / "
               "referenceGJ.Vect().Mag();"
            << endl
            << "double phi_H_PiPlus_PiMinus1_Photon1_Photon2 = "
               "TMath::ATan2(referenceGJ.Vect().Dot(y_H), referenceGJ.Vect().Dot(x_H));"
            << endl;
    }
    locSourceStream << "		/******************************************** EXECUTE "
                       "ANALYSIS "
                       "ACTIONS *******************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		// Loop through the analysis actions, "
                       "executing them in order for "
                       "the active particle combo"
                    << endl;
    locSourceStream << "		dAnalyzeCutActions->Perform_Action(); "
                       "// Must be executed before "
                       "Execute_Actions()"
                    << endl;
    locSourceStream << "		if(!Execute_Actions()) //if the active "
                       "combo fails a cut, "
                       "IsComboCutFlag automatically set"
                    << endl;
    locSourceStream << "			continue;" << endl;
    locSourceStream << endl;
    locSourceStream << "		//if you manually execute any actions, "
                       "and it fails a cut, be sure to call:"
                    << endl;
    locSourceStream << "			//dComboWrapper->Set_IsComboCut(true);" << endl;
    locSourceStream << endl;
    locSourceStream << "		/**************************************** EXAMPLE: "
                       "FILL CUSTOM "
                       "OUTPUT BRANCHES **************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		/*" << endl;
    locSourceStream << "		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);" << endl;
    locSourceStream << "		//for arrays below: 2nd argument is "
                       "value, 3rd is array index"
                    << endl;
    locSourceStream << "		//NOTE: By filling here, AFTER the cuts above, some "
                       "indices won't "
                       "be updated (and will be whatever they were from the last event)"
                    << endl;
    locSourceStream << "			//So, when you draw the "
                       "branch, be sure to cut on "
                       "\"IsComboCut\" to avoid these."
                    << endl;
    locSourceStream << "		"
                       "dTreeInterface->Fill_Fundamental<Float_t>(\"my_combo_array\", "
                       "-2*loc_i, loc_i);"
                    << endl;
    locSourceStream << "		"
                       "dTreeInterface->Fill_TObject<TLorentzVector>(\"my_p4_array\", "
                       "locMyComboP4, loc_i);"
                    << endl;
    locSourceStream << "		*/" << endl;
    locSourceStream << endl;
    locSourceStream << "		/**************************************** EXAMPLE: "
                       "BEST chi2 "
                       "METHOD *****************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "        // Need to uncomment the section computing "
                       "combo timing info "

                       "before running this block of code"
                    << endl;
    locSourceStream << "		//if(locUsedSoFar_Event == false)" << endl;
    locSourceStream << "		//{" << endl;
    locSourceStream << "			// Fill the histogram only "
                       "when the beam bunch is in-time. "
                    << endl;
    locSourceStream << "			//if(!locRelBeamBucket)" << endl;
    locSourceStream << "			//{" << endl;
    locSourceStream << "			//	"
                       "dHist_BeamEnergy_BestChiSq->Fill(locBeamP4.E());"
                    << endl;
    locSourceStream << "			//	locUsedSoFar_Event = true;" << endl;
    locSourceStream << "			//}" << endl;
    locSourceStream << "		//}" << endl;
    locSourceStream << endl;
    locSourceStream << "		/**************************************** EXAMPLE: "
                       "HISTOGRAM BEAM "
                       "ENERGY *****************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		//Histogram beam energy (if haven\'t already)" << endl;
    locSourceStream << "		if(locUsedSoFar_BeamEnergy.find(locBeamID) == "
                       "locUsedSoFar_BeamEnergy.end())"
                    << endl;
    locSourceStream << "		{" << endl;
    locSourceStream << "			dHist_BeamEnergy->Fill(locBeamP4.E()); // "
                       "Fills in-time and "
                       "out-of-time beam photon combos"
                    << endl;
    locSourceStream << "			"
                       "//dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); "
                       "// "
                       "Alternate version with accidental subtraction"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "			locUsedSoFar_BeamEnergy.insert(locBeamID);" << endl;
    locSourceStream << "		}" << endl;
    locSourceStream << endl;
    locSourceStream << "		/************************************ "
                       "EXAMPLE: HISTOGRAM MISSING "
                       "MASS SQUARED ************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		//Missing Mass Squared" << endl;

    locSourceStream << "		double locMissingMassSquared = "
                       "locMissingP4_Measured.M2();"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		//Uniqueness tracking: Build the map "
                       "of particles used for the missing mass"
                    << endl;
    locSourceStream << "			//For beam: Don\'t want to group with "
                       "final-state photons. "

                       "Instead use \"UnknownParticle\" PID (not ideal, but it\'s easy)."
                    << endl;
    locSourceStream << "		map<Particle_t, set<Int_t> > "
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
    locSourceStream << endl;

    locSourceStream << "		//compare to what\'s been used so far" << endl;
    locSourceStream << "		"
                       "if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == "
                       "locUsedSoFar_MissingMass.end())"
                    << endl;
    locSourceStream << "		{" << endl;
    locSourceStream << "			//unique missing mass combo: "
                       "histogram it, and register this "
                       "combo of particles"
                    << endl;
    locSourceStream << "			"
                       "dHist_MissingMassSquared->Fill(locMissingMassSquared); // "
                       "Fills in-time and out-of-time beam photon combos"
                    << endl;
    locSourceStream << "			"
                       "//"
                       "dHist_MissingMassSquared->Fill(locMissingMassSquared,"
                       "locHistAccidWeightFactor); "
                       "// "
                       "Alternate version with accidental subtraction"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "			"
                       "locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);"
                    << endl;
    locSourceStream << "		}" << endl;
    locSourceStream << endl;
    locSourceStream << "		//E.g. Cut" << endl;
    locSourceStream << "		//if((locMissingMassSquared < -0.04) "
                       "|| (locMissingMassSquared > 0.04))"
                    << endl;
    locSourceStream << "		//{" << endl;
    locSourceStream << "		//	dComboWrapper->Set_IsComboCut(true);" << endl;
    locSourceStream << "		//	continue;" << endl;
    locSourceStream << "		//}" << endl;
    locSourceStream << endl;
    locSourceStream << "		/****************************************** FILL FLAT "
                       "TREE (IF "
                       "DESIRED) ******************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		// RECOMMENDED: FILL ACCIDENTAL WEIGHT" << endl;
    locSourceStream << "		// "
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

        string uv = "PiPlus_PiMinus1_Photon1_Photon2"; // TODO: currently must be user defined
        // costh_GJ_
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_GJ_" << uv
                        << "\", costh_GJ_" << uv << ");\n";

        // phi_GJ
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_GJ_" << uv
                        << "\", phi_GJ_" << uv << ");\n";

        // costh_H_
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_H_" << uv
                        << "\", costh_H_" << uv << ");\n";

        // phi_H
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_H_" << uv
                        << "\", phi_H_" << uv << ");\n";
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

    locSourceStream << endl;
    locSourceStream << "		/*" << endl;
    locSourceStream << "		//FILL ANY CUSTOM BRANCHES FIRST!!" << endl;

    locSourceStream << "		Int_t locMyInt_Flat = 7;" << endl;
    locSourceStream << "		"
                       "dFlatTreeInterface->Fill_Fundamental<Int_t>(\"flat_my_int\", "
                       "locMyInt_Flat);"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);" << endl;
    locSourceStream << "		"
                       "dFlatTreeInterface->Fill_TObject<TLorentzVector>(\"flat_my_p4\", "
                       "locMyP4_Flat);"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)" << endl;
    locSourceStream << "		{" << endl;
    locSourceStream << "			"
                       "dFlatTreeInterface->Fill_Fundamental<Int_t>(\"flat_my_int_array\", "
                       "3*loc_j, loc_j); //2nd argument = value, 3rd = array index"
                    << endl;
    locSourceStream << "			TLorentzVector locMyComboP4_Flat(8.0, "
                       "7.0, 6.0, 5.0);"
                    << endl;
    locSourceStream << "			"
                       "dFlatTreeInterface->Fill_TObject<TLorentzVector>(\"flat_"
                       "my_p4_array\", "
                       "locMyComboP4_Flat, loc_j);"
                    << endl;
    locSourceStream << "		}" << endl;
    locSourceStream << "		*/" << endl;
    locSourceStream << endl;
    locSourceStream << "		//FILL FLAT TREE" << endl;
    locSourceStream << "		Fill_FlatTree(); //for the active combo" << endl;
    locSourceStream << "	} // end of combo loop" << endl;
    locSourceStream << endl;
    locSourceStream << "	//FILL HISTOGRAMS: Num combos / events surviving actions" << endl;
    locSourceStream << "	Fill_NumCombosSurvivedHists();" << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************************* LOOP OVER "
                       "THROWN DATA "
                       "(OPTIONAL) ***************************************/"
                    << endl;
    locSourceStream << "/*" << endl;
    locSourceStream << "	//Thrown beam: just use directly" << endl;
    locSourceStream << "	if(dThrownBeam != NULL)" << endl;
    locSourceStream << "		double locEnergy = dThrownBeam->Get_P4().E();" << endl;
    locSourceStream << endl;
    locSourceStream << "	//Loop over throwns" << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl;
    locSourceStream << "		dThrownWrapper->Set_ArrayIndex(loc_i);" << endl;
    locSourceStream << endl;
    locSourceStream << "		//Do stuff with the wrapper here ..." << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << "*/" << endl;
    locSourceStream << "	/****************************************** LOOP OVER "
                       "OTHER ARRAYS "
                       "(OPTIONAL) ***************************************/"
                    << endl;
    locSourceStream << "/*" << endl;
    locSourceStream << "	//Loop over beam particles (note, only those appearing in "
                       "combos are present)"
                    << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl;
    locSourceStream << "		dBeamWrapper->Set_ArrayIndex(loc_i);" << endl;
    locSourceStream << endl;
    locSourceStream << "		//Do stuff with the wrapper here ..." << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << endl;
    locSourceStream << "	//Loop over charged track hypotheses" << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl;
    locSourceStream << "		dChargedHypoWrapper->Set_ArrayIndex(loc_i);" << endl;
    locSourceStream << endl;
    locSourceStream << "		//Do stuff with the wrapper here ..." << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << endl;
    locSourceStream << "	//Loop over neutral particle hypotheses" << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl;
    locSourceStream << "		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);" << endl;
    locSourceStream << endl;
    locSourceStream << "		//Do stuff with the wrapper here ..." << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << "*/" << endl;
    locSourceStream << endl;
    locSourceStream << "	/************************************ EXAMPLE: FILL CLONE OF "
                       "TTREE "
                       "HERE WITH CUTS APPLIED ************************************/"
                    << endl;
    locSourceStream << "/*" << endl;
    locSourceStream << "	Bool_t locIsEventCut = true;" << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {" << endl;
    locSourceStream << "		//Set branch array indices for combo and all "
                       "combo particles"
                    << endl;
    locSourceStream << "		dComboWrapper->Set_ComboIndex(loc_i);" << endl;
    locSourceStream << "		// Is used to indicate when combos have been cut" << endl;
    locSourceStream << "		if(dComboWrapper->Get_IsComboCut())" << endl;
    locSourceStream << "			continue;" << endl;
    locSourceStream << "		locIsEventCut = false; // At least one combo succeeded" << endl;
    locSourceStream << "		break;" << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << "	if(!locIsEventCut && dOutputTreeFileName != \"\")" << endl;
    locSourceStream << "		Fill_OutputTree();" << endl;
    locSourceStream << "*/" << endl;
    locSourceStream << endl;
    locSourceStream << "	return kTRUE;" << endl;
    locSourceStream << "}" << endl;
    locSourceStream << endl;
    locSourceStream << "void " << locSelectorName << "::Finalize(void)" << endl;
    locSourceStream << "{" << endl;
    locSourceStream << "	//Save anything to output here that you do not want to "
                       "be in the "
                       "default DSelector output ROOT file."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Otherwise, don\'t do anything else (especially if "
                       "you are using PROOF)."
                    << endl;
    locSourceStream << "		//If you are using PROOF, this function is "
                       "called on each thread,"
                    << endl;
    locSourceStream << "		//so anything you do will not have the "
                       "combined information from "
                       "the various threads."
                    << endl;
    locSourceStream << "		//Besides, it is best-practice to do "
                       "post-processing (e.g. "
                       "fitting) separately, in case there is a problem."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//DO YOUR STUFF HERE" << endl;
    locSourceStream << endl;
    locSourceStream << "	//CALL THIS LAST" << endl;
    locSourceStream << "	DSelector::Finalize(); //Saves results to the output file" << endl;
    locSourceStream << "}" << endl;

    locSourceStream.close();
}

void Print_HeaderFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface)
{
    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locHeaderName = locSelectorName + string(".h");
    ofstream locHeaderStream;
    locHeaderStream.open(locHeaderName.c_str());

    locHeaderStream << "#ifndef " << locSelectorName << "_h" << endl;
    locHeaderStream << "#define " << locSelectorName << "_h" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#include <iostream>" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#include \"DSelector/DSelector.h\"" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#include \"TH1I.h\"" << endl;
    locHeaderStream << "#include \"TH2I.h\"" << endl;
    locHeaderStream << endl;
    locHeaderStream << "class " << locSelectorName << " : public DSelector" << endl;
    locHeaderStream << "{" << endl;
    locHeaderStream << "	public:" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		" << locSelectorName
                    << "(TTree* locTree = NULL) : DSelector(locTree){}" << endl;
    locHeaderStream << "		virtual ~" << locSelectorName << "(){}" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		void Init(TTree *tree);" << endl;
    locHeaderStream << "		Bool_t Process(Long64_t entry);" << endl;
    locHeaderStream << endl;
    locHeaderStream << "	private:" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		void Finalize(void);" << endl;
    locHeaderStream << endl;
    locHeaderStream << "		// BEAM POLARIZATION INFORMATION" << endl;
    locHeaderStream << "		UInt_t dPreviousRunNumber;" << endl;
    locHeaderStream << "		bool dIsPolarizedFlag; //else is AMO" << endl;
    locHeaderStream << "		bool dIsPARAFlag; //else is PERP or AMO" << endl;
    locHeaderStream << endl;
    locHeaderStream << "	ClassDef(" << locSelectorName << ", 0);" << endl;
    locHeaderStream << "};" << endl;
    locHeaderStream << endl;
    locHeaderStream << "#endif // " << locSelectorName << "_h" << endl;

    locHeaderStream.close();
}

void Print_SourceFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                            bool extraDefaults)
{
    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locSourceName = locSelectorName + string(".C");
    ofstream locSourceStream;
    locSourceStream.open(locSourceName.c_str());

    locSourceStream << "#include \"" << locSelectorName << ".h\"" << endl;
    locSourceStream << endl;
    locSourceStream << "void " << locSelectorName << "::Init(TTree *locTree)" << endl;
    locSourceStream << "{" << endl;
    locSourceStream << "	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH "
                       "A \"USER\" OR "
                       "\"EXAMPLE\" LABEL. LEAVE THE REST ALONE."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	// The Init() function is called when the selector "
                       "needs to initialize "
                       "a new tree or chain."
                    << endl;
    locSourceStream << "	// Typically here the branch addresses and branch "
                       "pointers of the "
                       "tree will be set."
                    << endl;
    locSourceStream << "	// Init() will be called many times when running on "
                       "PROOF (once per "
                       "file to be processed)."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//USERS: SET OUTPUT FILE NAME //can be overriden by "
                       "user in PROOF"
                    << endl;
    locSourceStream << "	dFlatTreeFileName = \"" << locSelectorBaseName
                    << ".root\"; //\"\" for none" << endl;
    locSourceStream << "	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning "
                       "into separate "
                       "files for AmpTools"
                    << endl;
    locSourceStream << "	//dOutputTreeFileNameMap[\"Bin1\"] = "
                       "\"mcgen_bin1.root\"; //key is "
                       "user-defined, value is output file name"
                    << endl;
    locSourceStream << "	//dOutputTreeFileNameMap[\"Bin2\"] = "
                       "\"mcgen_bin2.root\"; //key is "
                       "user-defined, value is output file name"
                    << endl;
    locSourceStream << "	//dOutputTreeFileNameMap[\"Bin3\"] = "
                       "\"mcgen_bin3.root\"; //key is "
                       "user-defined, value is output file name"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Because this function gets called for each TTree in "
                       "the TChain, we "
                       "must be careful:"
                    << endl;
    locSourceStream << "		//We need to re-initialize the tree interface "
                       "& branch wrappers, "
                       "but don't want to recreate histograms"
                    << endl;
    locSourceStream << "	bool locInitializedPriorFlag = dInitializedFlag; "
                       "//save whether have "
                       "been initialized previously"
                    << endl;
    locSourceStream << "	DSelector::Init(locTree); //This must be called to "
                       "initialize wrappers "
                       "for each new TTree"
                    << endl;
    locSourceStream << "	//gDirectory now points to the output file with name "
                       "dOutputFileName (if any)"
                    << endl;
    locSourceStream << "	if(locInitializedPriorFlag)" << endl;
    locSourceStream << "		return; //have already created histograms, "
                       "etc. below: exit"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	dPreviousRunNumber = 0;" << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************** EXAMPLE USER INITIALIZATION: "
                       "STAND-ALONE HISTOGRAMS *******************************/"
                    << endl;

    if (extraDefaults)
    {
        ifstream csvIn("branches.csv");
        locSourceStream << "    // == Extra default branches ==\n";
        vector<TString> sel;
        std::string line;
        // read in and sanitize combo names
        while (std::getline(csvIn, line))
        {
            if (!line.empty())
            {
                line.erase(line.find_last_not_of(" \n\r\t") + 1);
                sel.push_back(line);
            }
        }
        for (auto &name : sel)
        {

            locSourceStream << "    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"mass_"
                            << name << "\", loc" << name << "P4.M());\n";
            locSourceStream
                << "    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"costh_lab_"
                << name << "\", loc" << name << "P4.Vect().CosTheta();\n";
            locSourceStream
                << "    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"phi_lab_" << name
                << "\", loc" << name << "P4.Vect().Phi();\n";
        }

        string uv = "PiPlus_PiMinus1_Photon1_Photon2"; // TODO: currently must be user defined
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<"
                           "Double_t>(\"costh_GJ_"
                        << uv << "\");\n";
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<"
                           "Double_t>(\"phi_GJ_"
                        << uv << "\");\n";
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                           "(\"costh_H_"
                        << uv << "\");\n";
        locSourceStream << "    "
                           "dFlatTreeInterface->Create_Branch_Fundamental<Double_t>"
                           "(\"phi_H_"
                        << uv << "\");\n";
        // individual masses
        locSourceStream << "    // == End extra defaults ==\n";
    }
    locSourceStream << endl;
    locSourceStream << "	/************************************* ADVANCED "
                       "EXAMPLE: CHOOSE "
                       "BRANCHES TO READ ************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//TO SAVE PROCESSING TIME" << endl;
    locSourceStream << "		//If you know you don't need all of the "
                       "branches/data, but just a "
                       "subset of it, you can speed things up"
                    << endl;
    locSourceStream << "		//By default, for each event, the data is "
                       "retrieved for all branches"
                    << endl;
    locSourceStream << "		//If you know you only need data for some "
                       "branches, you can skip "
                       "grabbing data from the branches you don't need"
                    << endl;
    locSourceStream << "		//Do this by doing something similar to the "
                       "commented code below"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//dTreeInterface->Clear_GetEntryBranches(); //now get none" << endl;
    locSourceStream << "	//dTreeInterface->Register_GetEntryBranch(\"Proton__P4\"); "
                       "//manually "
                       "set the branches you want"
                    << endl;
    locSourceStream << "}" << endl;
    locSourceStream << endl;

    locSourceStream << "Bool_t " << locSelectorName << "::Process(Long64_t locEntry)" << endl;
    locSourceStream << "{" << endl;
    locSourceStream << "	// The Process() function is called for each entry in "
                       "the tree. "
                       "The entry argument"
                    << endl;
    locSourceStream << "	// specifies which entry in the currently loaded tree "
                       "is to be processed."
                    << endl;
    locSourceStream << "	//" << endl;
    locSourceStream << "	// This function should contain the \"body\" of the analysis. "
                       "It can contain"
                    << endl;
    locSourceStream << "	// simple or elaborate selection criteria, run "
                       "algorithms on the data"

                    << endl;
    locSourceStream << "	// of the event and typically fill histograms." << endl;
    locSourceStream << "	//" << endl;
    locSourceStream << "	// The processing can be stopped by calling Abort()." << endl;
    locSourceStream << "	// Use fStatus to set the return value of TTree::Process()." << endl;
    locSourceStream << "	// The return value is currently not used." << endl;
    locSourceStream << endl;
    locSourceStream << "	//CALL THIS FIRST" << endl;
    locSourceStream << "	DSelector::Process(locEntry); //Gets the data from the "
                       "tree for the entry"
                    << endl;
    locSourceStream << "	//cout << \"RUN \" << Get_RunNumber() << \", EVENT \" << "
                       "Get_EventNumber() << endl;"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************************** GET POLARIZATION "
                       "ORIENTATION ******************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Only if the run number changes" << endl;
    locSourceStream << "	//RCDB environment must be setup in order for this to "
                       "work! (Will "
                       "return false otherwise)"
                    << endl;
    locSourceStream << "	UInt_t locRunNumber = Get_RunNumber();" << endl;
    locSourceStream << "	if(locRunNumber != dPreviousRunNumber)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		dIsPolarizedFlag = "
                       "dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);"
                    << endl;
    locSourceStream << "		dPreviousRunNumber = locRunNumber;" << endl;
    locSourceStream << "	}" << endl;
    locSourceStream << endl;
    locSourceStream << "	/********************************************* SETUP "
                       "UNIQUENESS "
                       "TRACKING ********************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE" << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************************* LOOP OVER "
                       "THROWN DATA "
                       "***************************************/"
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Thrown beam: just use directly" << endl;
    locSourceStream << "	double locBeamEnergyUsedForBinning = 0.0;" << endl;
    locSourceStream << "	if(dThrownBeam != NULL)" << endl;
    locSourceStream << "		locBeamEnergyUsedForBinning = "
                       "dThrownBeam->Get_P4().E();"
                    << endl;
    if (extraDefaults)
    {
        ifstream csvIn("branches.csv");
        vector<string> sel;
        std::string line;
        // read in and sanitize combo names
        while (std::getline(csvIn, line))
        {
            if (!line.empty())
            {
                line.erase(line.find_last_not_of(" \n\r\t") + 1);
                sel.push_back(line);
            }
        }
        locSourceStream << "    TLorentzVector locZeroVec(0,0,0,0);";
        for (auto &name : sel)
        {
            if (name.rfind("_") == 0)
            {
                locSourceStream << "TLorentzVector loc" << name << "P4 = locZeroVec;" << endl;
                locSourceStream << "UInt_t" << name << "_number = 0;" << endl;
            }
        }
        locSourceStream << "TLorentzVector dTargetP4(0,0,0,0.938);" << endl
                        << "TLorentzVector locBeamP4 dThrownBeam->Get_P4();" << endl;
    }
    locSourceStream << endl;
    locSourceStream << "	//Loop over throwns" << endl;
    locSourceStream << "	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)" << endl;
    locSourceStream << "	{" << endl;
    locSourceStream << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl;
    locSourceStream << "		dThrownWrapper->Set_ArrayIndex(loc_i);" << endl;
    locSourceStream << endl;
    locSourceStream << "		//Do stuff with the wrapper here ..." << endl;
    locSourceStream << "		Particle_t locPID = dThrownWrapper->Get_PID();" << endl;
    locSourceStream << "		TLorentzVector locThrownP4 = dThrownWrapper->Get_P4();"

                    << endl;
    locSourceStream << "		//cout << \"Thrown \" << loc_i << \": \" << locPID << "
                       "\", \" << "
                       "locThrownP4.Px() << \", \" << locThrownP4.Py() << \", \" << "
                       "locThrownP4.Pz() << \", \" << locThrownP4.E() << endl;"
                    << endl;
    if (extraDefaults)
    {
        // FINISH CODE
    }
    locSourceStream << "	}" << endl;
    locSourceStream << endl;
    locSourceStream << "	//OR Manually:" << endl;
    locSourceStream << "	//BEWARE: Do not expect the particles to be at the "
                       "same array indices "
                       "from one event to the next!!!!"
                    << endl;
    locSourceStream << "	//Why? Because while your channel may be the same, the "
                       "pions/kaons/etc. will decay differently each event."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//BRANCHES: "
                       "https://halldweb.jlab.org/wiki/index.php/"
                       "Analysis_TTreeFormat#TTree_Format:_Simulated_Data"
                    << endl;
    locSourceStream << "	TClonesArray** locP4Array = "
                       "dTreeInterface->Get_PointerToPointerTo_TClonesArray(\"Thrown__P4\");"
                    << endl;
    locSourceStream << "	TBranch* locPIDBranch = "
                       "dTreeInterface->Get_Branch(\"Thrown__PID\");"
                    << endl;
    locSourceStream << "/*" << endl;
    locSourceStream << "	Particle_t locThrown1PID = "
                       "PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);"
                    << endl;
    locSourceStream << "	TLorentzVector locThrown1P4 = "
                       "*((TLorentzVector*)(*locP4Array)->At(0));"
                    << endl;
    locSourceStream << "	cout << \"Particle 1: \" << locThrown1PID << \", \" << "
                       "locThrown1P4.Px() << \", \" << locThrown1P4.Py() << \", \" << "
                       "locThrown1P4.Pz() << \", \" << locThrown1P4.E() << endl;"
                    << endl;
    locSourceStream << "	Particle_t locThrown2PID = "
                       "PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);"
                    << endl;
    locSourceStream << "	TLorentzVector locThrown2P4 = "
                       "*((TLorentzVector*)(*locP4Array)->At(1));"
                    << endl;
    locSourceStream << "	cout << \"Particle 2: \" << locThrown2PID << \", \" << "
                       "locThrown2P4.Px() << \", \" << locThrown2P4.Py() << \", \" << "
                       "locThrown2P4.Pz() << \", \" << locThrown2P4.E() << endl;"
                    << endl;
    locSourceStream << "*/" << endl;
    locSourceStream << endl;
    locSourceStream << endl;
    locSourceStream << "	/******************************************* BIN THROWN DATA "
                       "INTO "
                       "SEPARATE TREES FOR AMPTOOLS ***************************************/"
                    << endl;

    if (extraDefaults)
    {
        locSourceStream << "    // == Extra default fills ==\n";
        // individual masses
        ifstream csvIn("branches.csv");
        vector<string> sel;
        std::string line;
        // read in and sanitize combo names
        while (std::getline(csvIn, line))
        {
            if (!line.empty())
            {
                line.erase(line.find_last_not_of(" \n\r\t") + 1);
                sel.push_back(line);
            }
        }
        for (auto &name : sel)
        {
            if (name.rfind("_") == 0)
            {
                locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                                << name << "\", loc" << name << "P4.M());\n";
            }
            else
            {

                vector<string> particles = tokenize(name);

                // now generate all r‐body combinations
                for (size_t r = 2; r <= particles.size(); ++r)
                {
                    vector<int> idx(r);
                    function<void(int, unsigned long)> comb = [&](int start,
                                                                  unsigned long int depth) {
                        if (depth == r)
                        {
                            // kinematic‐fit mass
                            locSourceStream
                                << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                                << name << "\", (";
                            for (unsigned long int i = 0; i < r; ++i)
                            {
                                if (i)
                                {
                                    locSourceStream << " + ";
                                }
                                locSourceStream << "loc" << particles[idx[i]] << "P4";
                            }
                            locSourceStream << ").M());\n";

                            // costh_lab
                            locSourceStream
                                << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_lab_"
                                << name << "\", (";
                            for (unsigned long int i = 0; i < r; ++i)
                            {
                                if (i)
                                {
                                    locSourceStream << " + ";
                                }
                                locSourceStream << "loc" << particles[idx[i]] << "P4";
                            }
                            locSourceStream << ").Vect().CosTheta());\n";

                            // phi_lab
                            locSourceStream
                                << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_lab_"
                                << name << "\", (";
                            for (unsigned long int i = 0; i < r; ++i)
                            {
                                if (i)
                                {
                                    locSourceStream << " + ";
                                }
                                locSourceStream << "loc" << particles[idx[i]] << "P4";
                            }
                            locSourceStream << ").Vect().Phi());\n";

                            return;
                        }
                        for (int i = start; i < int(particles.size()); ++i)
                        {
                            idx[depth] = i;
                            comb(i + 1, depth + 1);
                        }
                    };
                    comb(0, 0);
                }
            }
        }
        string uv = "PiPlus_PiMinus1_Photon1_Photon2"; // TODO: currently must be user defined
        // costh_GJ_
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_GJ_" << uv
                        << "\", costh_GJ_" << uv << ");\n";

        // phi_GJ
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_GJ_" << uv
                        << "\", phi_GJ_" << uv << ");\n";

        // costh_H_
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_H_" << uv
                        << "\", costh_H_" << uv << ");\n";

        // phi_H
        locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_H_" << uv
                        << "\", phi_H_" << uv << ");\n";
    }
    locSourceStream << endl;
    locSourceStream << "/*" << endl;
    locSourceStream << "	//THESE KEYS MUST BE DEFINED IN THE INIT SECTION "
                       "(along with the "
                       "output file names)"
                    << endl;
    locSourceStream << "	if((locBeamEnergyUsedForBinning >= 8.0) && "
                       "(locBeamEnergyUsedForBinning < 9.0))"
                    << endl;
    locSourceStream << "		Fill_OutputTree(\"Bin1\"); //your user-defined key" << endl;
    locSourceStream << "	else if((locBeamEnergyUsedForBinning >= 9.0) && "
                       "(locBeamEnergyUsedForBinning < 10.0))"
                    << endl;
    locSourceStream << "		Fill_OutputTree(\"Bin2\"); //your user-defined key" << endl;
    locSourceStream << "	else if((locBeamEnergyUsedForBinning >= 10.0) && "
                       "(locBeamEnergyUsedForBinning < 11.0))"
                    << endl;
    locSourceStream << "		Fill_OutputTree(\"Bin3\"); //your user-defined key" << endl;
    locSourceStream << "*/" << endl;
    locSourceStream << endl;
    locSourceStream << "	return kTRUE;" << endl;
    locSourceStream << "}" << endl;
    locSourceStream << endl;
    locSourceStream << "void " << locSelectorName << "::Finalize(void)" << endl;
    locSourceStream << "{" << endl;
    locSourceStream << "	//Save anything to output here that you do not want to "
                       "be in the "
                       "default DSelector output ROOT file."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//Otherwise, don\'t do anything else (especially if "
                       "you are using PROOF)."
                    << endl;
    locSourceStream << "		//If you are using PROOF, this function is "
                       "called on each thread,"
                    << endl;
    locSourceStream << "		//so anything you do will not have the "
                       "combined information from "
                       "the various threads."
                    << endl;
    locSourceStream << "		//Besides, it is best-practice to do "
                       "post-processing (e.g. "
                       "fitting) separately, in case there is a problem."
                    << endl;
    locSourceStream << endl;
    locSourceStream << "	//DO YOUR STUFF HERE" << endl;
    locSourceStream << endl;
    locSourceStream << "	//CALL THIS LAST" << endl;
    locSourceStream << "	DSelector::Finalize(); //Saves results to the output file" << endl;
    locSourceStream << "}" << endl;

    locSourceStream.close();
}
