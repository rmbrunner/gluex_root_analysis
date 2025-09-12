#include "MakeDSelector.h"

static std::string NormalizeName(const std::string &name)
{
    auto remapSuffix = [&](const std::string &suffix,
                           const std::string &replacement) -> std::string {
        if (name.size() >= suffix.size() &&
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
        {
            return name.substr(0, name.size() - suffix.size()) + replacement;
        }
        return {};
    };
    string s = remapSuffix("Minus", "-");
    if (!s.empty())
    {
        return s;
    }
    s = remapSuffix("Plus", "+");
    if (!s.empty())
    {
        return s;
    }

    return name;
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
static std::size_t count_if_contains(const std::vector<std::string> &vec, const std::string &sub)
{
    return std::count_if(vec.begin(), vec.end(),
                         [&](const std::string &s) { return s.find(sub) != std::string::npos; });
}

void Print_SourceFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                            bool extraDefaults)
{

    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locSourceName = locSelectorName + string(".C");
    ofstream locSourceStream;
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
                       "pointers of the "
                       "tree will be set."
                    << endl
                    << "	// Init() will be called many times when running on "
                       "PROOF (once per "
                       "file to be processed)."
                    << endl
                    << endl
                    << "	//USERS: SET OUTPUT FILE NAME //can be overriden by "
                       "user in PROOF"
                    << endl
                    << "	dFlatTreeFileName = \"" << locSelectorBaseName
                    << ".root\"; //\"\" for none" << endl
                    << "	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning "
                       "into separate "
                       "files for AmpTools"
                    << endl
                    << "	//dOutputTreeFileNameMap[\"Bin1\"] = "
                       "\"mcgen_bin1.root\"; //key is "
                       "user-defined, value is output file name"
                    << endl
                    << "	//dOutputTreeFileNameMap[\"Bin2\"] = "
                       "\"mcgen_bin2.root\"; //key is "
                       "user-defined, value is output file name"
                    << endl
                    << "	//dOutputTreeFileNameMap[\"Bin3\"] = "
                       "\"mcgen_bin3.root\"; //key is "
                       "user-defined, value is output file name"
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
                    << "	dPreviousRunNumber = 0;" << endl
                    << endl
                    << "	/******************************** EXAMPLE USER INITIALIZATION: "
                       "STAND-ALONE HISTOGRAMS *******************************/"
                    << endl;

    string finalState;
    if (extraDefaults)
    {
        ifstream csvIn("branches.csv");
        locSourceStream << "    // == Extra default branches ==\n";
        vector<TString> sel;
        string line;
        // read in and sanitize combo names
        while (std::getline(csvIn, line))
        {
            if (!line.empty())
            {
                line.erase(line.find_last_not_of(" \n\r\t") + 1);
                sel.push_back(line);
                finalState = line;
            }
        }
        for (auto &name : sel)
        {

            locSourceStream << "    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"mass_"
                            << name << "\");\n";
            locSourceStream
                << "    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"costh_lab_"
                << name << "\");\n";
            locSourceStream
                << "    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(\"phi_lab_" << name
                << "\");\n";
            if (name != finalState)
            {
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
            }
        }

        // individual masses
        locSourceStream << "    // == End extra defaults ==\n";
        csvIn.close();
    }
    locSourceStream << endl
                    << "	/************************************* ADVANCED "
                       "EXAMPLE: CHOOSE "
                       "BRANCHES TO READ ************************************/"
                    << endl
                    << endl
                    << "	//TO SAVE PROCESSING TIME" << endl
                    << "		//If you know you don't need all of the "
                       "branches/data, but just a "
                       "subset of it, you can speed things up"
                    << endl
                    << "		//By default, for each event, the data is "
                       "retrieved for all branches"
                    << endl
                    << "		//If you know you only need data for some "
                       "branches, you can skip "
                       "grabbing data from the branches you don't need"
                    << endl
                    << "		//Do this by doing something similar to the "
                       "commented code below"
                    << endl
                    << endl
                    << "	//dTreeInterface->Clear_GetEntryBranches(); //now get none" << endl
                    << "	//dTreeInterface->Register_GetEntryBranch(\"Proton__P4\"); "
                       "//manually "
                       "set the branches you want"
                    << endl
                    << "}" << endl
                    << endl

                    << "Bool_t " << locSelectorName << "::Process(Long64_t locEntry)" << endl
                    << "{" << endl
                    << "	// The Process() function is called for each entry in "
                       "the tree. "
                       "The entry argument"
                    << endl
                    << "	// specifies which entry in the currently loaded tree "
                       "is to be processed."
                    << endl
                    << "	//" << endl
                    << "	// This function should contain the \"body\" of the analysis. "
                       "It can contain"
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
                    << "	DSelector::Process(locEntry); //Gets the data from the "
                       "tree for the entry"
                    << endl
                    << "	//cout << \"RUN \" << Get_RunNumber() << \", EVENT \" << "
                       "Get_EventNumber() << endl;"
                    << endl
                    << endl
                    << "	/******************************************** GET POLARIZATION "
                       "ORIENTATION ******************************************/"
                    << endl
                    << endl
                    << "	//Only if the run number changes" << endl
                    << "	//RCDB environment must be setup in order for this to "
                       "work! (Will "
                       "return false otherwise)"
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
                    << "	/********************************************* SETUP "
                       "UNIQUENESS "
                       "TRACKING ********************************************/"
                    << endl
                    << endl
                    << "	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE" << endl
                    << endl
                    << "	/******************************************* LOOP OVER "
                       "THROWN DATA "
                       "***************************************/"
                    << endl
                    << endl
                    << "	//Thrown beam: just use directly" << endl
                    << "	double locBeamEnergyUsedForBinning = 0.0;" << endl
                    << "	if(dThrownBeam != NULL)" << endl
                    << "		locBeamEnergyUsedForBinning = "
                       "dThrownBeam->Get_P4().E();"
                    << endl;
    if (extraDefaults)
    {
        ifstream csvIn("branches.csv");
        vector<string> sel;
        string line;
        // read in and sanitize combo names
        while (std::getline(csvIn, line))
        {
            if (!line.empty())
            {
                line.erase(line.find_last_not_of(" \n\r\t") + 1);
                sel.push_back(line);
            }
        }
        set<string> primitives;
        for (auto &name : sel)
        {
            if (name.rfind("_") == string::npos)
            {
                string primitive = name.substr(0, name.find_last_not_of("0123456789") + 1);
                primitives.insert(primitive);
            }
        }

        locSourceStream << "    TLorentzVector locZeroVec(0,0,0,0);";
        for (auto &name : sel)
        {
            if (name.rfind("_") == string::npos)
            {
                locSourceStream << "TLorentzVector loc" << name << "P4 = locZeroVec;" << endl;
            }
        }
        for (auto &name : primitives)
        {
            locSourceStream << "    UInt_t " << name << "_number = 0;" << endl;
        }
        locSourceStream << "TLorentzVector dTargetP4(0,0,0,0.938);" << endl
                        << "TLorentzVector locBeamP4 = dThrownBeam->Get_P4();" << endl;
        csvIn.close();
    }

    locSourceStream << endl
                    << "	//Loop over throwns" << endl
                    << "	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)" << endl
                    << "	{" << endl
                    << "		//Set branch array indices corresponding to "
                       "this particle"
                    << endl
                    << "		dThrownWrapper->Set_ArrayIndex(loc_i);" << endl
                    << endl
                    << "		//Do stuff with the wrapper here ..." << endl
                    << "		Particle_t locPID = dThrownWrapper->Get_PID();" << endl
                    << "		TLorentzVector locThrownP4 = dThrownWrapper->Get_P4();"

                    << endl
                    << "		//cout << \"Thrown \" << loc_i << \": \" << locPID << "
                       "\", \" << "
                       "locThrownP4.Px() << \", \" << locThrownP4.Py() << \", \" << "
                       "locThrownP4.Pz() << \", \" << locThrownP4.E() << endl;"
                    << endl;
    if (extraDefaults)
    {
        ifstream csvIn("branches.csv");
        vector<string> sel;
        std::string line;
        // read in and sanitize combo names
        while (std::getline(csvIn, line))
        {
            if (!line.empty() && line.rfind("_") == string::npos)
            {
                line.erase(line.find_last_not_of(" \n\r\t") + 1);
                sel.push_back(line);
            }
        }
        set<string> primitives;
        for (auto &name : sel)
        {
            string primitive = name.substr(0, name.find_last_not_of("0123456789") + 1);
            primitives.insert(primitive);
        }
        for (auto &primitive : primitives)
        {
            locSourceStream << "    if(locPID == " << ParticleEnum(NormalizeName(primitive).c_str())
                            << ")" << endl
                            << "    {" << endl;
            if (count_if_contains(sel, primitive) == 1)
            {
                locSourceStream << "        " << primitive << "_number++;" << endl
                                << "        loc" << primitive << "P4 = locThrownP4;" << endl;
            }
            else
            {
                locSourceStream << "        " << primitive << "_number++;" << endl;

                for (UInt_t i = 0; i < count_if_contains(sel, primitive); i++)
                {
                    locSourceStream << "            "
                                    << "if (" << primitive << "_number == " << i + 1 << ")" << endl
                                    << "            {"
                                    << "                loc" << primitive << (i + 1)
                                    << "P4 = locThrownP4;"
                                    << "            }" << endl
                                    << endl;
                }
            }
            locSourceStream << "    }" << endl;
        }
    }
    locSourceStream << "	}" << endl;

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
        locSourceStream
            << "    bool pimIsNull = " << endl
            << "        ((locPiMinus1P4 == locZeroVec) || (locPiMinus2P4 == locZeroVec));" << endl
            << "    if (!pimIsNull) {" << endl;

        locSourceStream << "TLorentzVector CM_P4 = locBeamP4 + dTargetP4;" << endl
                        << "TLorentzRotation CM_Boost(-CM_P4.BoostVector());" << endl
                        << "TLorentzVector particleXP4 = CM_Boost * (locPiPlusP4 + locPiMinus1P4 "
                           "+ locPhoton1P4 + locPhoton1P4); // For now must be user defined!!!"
                        << endl;
        locSourceStream << "TLorentzVector beamCM = CM_Boost * locBeamP4;" << endl
                        << "TLorentzVector targetCM = CM_Boost * dTargetP4;" << endl;

        for (auto &name : sel)
        {
            if (name.rfind("_") == string::npos)
            {
                locSourceStream << "TLorentzVector " << name << "CM = CM_Boost *"
                                << " loc" << name << "P4;" << endl;
            }
        }

        locSourceStream << endl
                        << "TLorentzVector particleXCM = CM_Boost * particleXP4;" << endl
                        << endl
                        << "TLorentzRotation restFrameXBoost(-particleXCM.BoostVector());" << endl;

        for (auto &branch : sel)
        {
            locSourceStream << "TLorentzVector referenceGJ_" << branch << " = restFrameXBoost * (";
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
        }
        locSourceStream << endl
                        << "TLorentzVector beamGJ = restFrameXBoost * beamCM;" << endl
                        << "TLorentzVector targetGJ = restFrameXBoost * targetCM;" << endl;

        for (auto &name : sel)
        {
            if (name.rfind("_") == string::npos)
            {
                locSourceStream << "TLorentzVector " << name << "GJ = restFrameXBoost * " << name
                                << "CM;" << endl
                                << "TVector3 " << name << "P3 = " << name << "GJ.Vect();" << endl;
            }
        }
        for (auto &branch : sel)
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
        csvIn.close();
    }

    locSourceStream << endl
                    << "	//OR Manually:" << endl
                    << "	//BEWARE: Do not expect the particles to be at the "
                       "same array indices "
                       "from one event to the next!!!!"
                    << endl
                    << "	//Why? Because while your channel may be the same, the "
                       "pions/kaons/etc. will decay differently each event."
                    << endl
                    << endl
                    << "	//BRANCHES: "
                       "https://halldweb.jlab.org/wiki/index.php/"
                       "Analysis_TTreeFormat#TTree_Format:_Simulated_Data"
                    << endl
                    << "	TClonesArray** locP4Array = "
                       "dTreeInterface->Get_PointerToPointerTo_TClonesArray(\"Thrown__P4\");"
                    << endl
                    << "	TBranch* locPIDBranch = "
                       "dTreeInterface->Get_Branch(\"Thrown__PID\");"
                    << endl
                    << "/*" << endl
                    << "	Particle_t locThrown1PID = "
                       "PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);"
                    << endl
                    << "	TLorentzVector locThrown1P4 = "
                       "*((TLorentzVector*)(*locP4Array)->At(0));"
                    << endl
                    << "	cout << \"Particle 1: \" << locThrown1PID << \", \" << "
                       "locThrown1P4.Px() << \", \" << locThrown1P4.Py() << \", \" << "
                       "locThrown1P4.Pz() << \", \" << locThrown1P4.E() << endl;"
                    << endl
                    << "	Particle_t locThrown2PID = "
                       "PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);"
                    << endl
                    << "	TLorentzVector locThrown2P4 = "
                       "*((TLorentzVector*)(*locP4Array)->At(1));"
                    << endl
                    << "	cout << \"Particle 2: \" << locThrown2PID << \", \" << "
                       "locThrown2P4.Px() << \", \" << locThrown2P4.Py() << \", \" << "
                       "locThrown2P4.Pz() << \", \" << locThrown2P4.E() << endl;"
                    << endl
                    << "*/" << endl
                    << endl
                    << endl
                    << "	/******************************************* BIN THROWN DATA "
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
            if (name.rfind("_") == string::npos)
            {
                locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                                << name << "\", loc" << name << "P4.M());\n";
            }
            else
            {

                vector<string> particles = tokenize(name);

                // kinematic‐fit mass
                locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"mass_"
                                << name << "\", (";
                for (unsigned long int i = 0; i < particles.size(); i++)
                {
                    if (i)
                    {
                        locSourceStream << " + ";
                    }
                    locSourceStream << "loc" << particles[i] << "P4";
                }
                locSourceStream << ").M());\n";

                // costh_lab
                locSourceStream << "    "
                                   "dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_lab_"
                                << name << "\", (";
                for (unsigned long int i = 0; i < particles.size(); i++)
                {
                    if (i)
                    {
                        locSourceStream << " + ";
                    }
                    locSourceStream << "loc" << particles[i] << "P4";
                }
                locSourceStream << ").Vect().CosTheta());\n";

                // phi_lab
                locSourceStream << "    "
                                   "dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_lab_"
                                << name << "\", (";
                for (unsigned long int i = 0; i < particles.size(); i++)
                {
                    if (i)
                    {
                        locSourceStream << " + ";
                    }
                    locSourceStream << "loc" << particles[i] << "P4";
                }
                locSourceStream << ").Vect().Phi());\n";

                if (name != finalState)
                {
                    // costh_GJ_
                    locSourceStream
                        << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_GJ_" << name
                        << "\", costh_GJ_" << name << ");\n";

                    // phi_GJ
                    locSourceStream
                        << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_GJ_" << name
                        << "\", phi_GJ_" << name << ");\n";

                    // costh_H_
                    locSourceStream
                        << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"costh_H_" << name
                        << "\", costh_H_" << name << ");\n";

                    // phi_H
                    locSourceStream << "    dFlatTreeInterface->Fill_Fundamental<Double_t>(\"phi_H_"
                                    << name << "\", phi_H_" << name << ");\n";
                }
            }
        }
        csvIn.close();
    }
    locSourceStream << endl
                    << "/*" << endl
                    << "	//THESE KEYS MUST BE DEFINED IN THE INIT SECTION "
                       "(along with the "
                       "output file names)"
                    << endl
                    << "	if((locBeamEnergyUsedForBinning >= 8.0) && "
                       "(locBeamEnergyUsedForBinning < 9.0))"
                    << endl
                    << "		Fill_OutputTree(\"Bin1\"); //your user-defined key" << endl
                    << "	else if((locBeamEnergyUsedForBinning >= 9.0) && "
                       "(locBeamEnergyUsedForBinning < 10.0))"
                    << endl
                    << "		Fill_OutputTree(\"Bin2\"); //your user-defined key" << endl
                    << "	else if((locBeamEnergyUsedForBinning >= 10.0) && "
                       "(locBeamEnergyUsedForBinning < 11.0))"
                    << endl
                    << "		Fill_OutputTree(\"Bin3\"); //your user-defined key" << endl
                    << "*/" << endl
                    << "Fill_FlatTree();" << endl
                    << "}" << endl
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
