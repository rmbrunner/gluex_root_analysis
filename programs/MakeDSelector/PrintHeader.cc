#include "MakeDSelector.h"

void Print_HeaderFile(string locSelectorBaseName, DTreeInterface *locTreeInterface,
                      map<int, map<int, pair<Particle_t, string>>> &locComboInfoMap)
{
    string locSelectorName = string("DSelector_") + locSelectorBaseName;
    string locHeaderName = locSelectorName + string(".h");
    ofstream locHeaderStream;
    locHeaderStream.open(locHeaderName.c_str());

    locHeaderStream << "#ifndef " << locSelectorName << "_h" << endl
                    << "#define " << locSelectorName << "_h" << endl
                    << endl
                    << "#include <iostream>" << endl
                    << endl
                    << "#include \"DSelector/DSelector.h\"" << endl
                    << "#include \"DSelector/DHistogramActions.h\"" << endl
                    << "#include \"DSelector/DCutActions.h\"" << endl
                    << endl
                    << "#include \"TH1I.h\"" << endl
                    << "#include \"TH2I.h\"" << endl
                    << endl
                    << "class " << locSelectorName << " : public DSelector" << endl
                    << "{" << endl
                    << "	public:" << endl
                    << endl
                    << "		" << locSelectorName
                    << "(TTree* locTree = NULL) : DSelector(locTree){}" << endl
                    << "		virtual ~" << locSelectorName << "(){}" << endl
                    << endl
                    << "		void Init(TTree *tree);" << endl
                    << "		Bool_t Process(Long64_t entry);" << endl
                    << endl
                    << "	private:" << endl
                    << endl
                    << "		void Get_ComboWrappers(void);" << endl
                    << "		void Finalize(void);" << endl
                    << endl
                    << "		// BEAM POLARIZATION INFORMATION" << endl
                    << "		UInt_t dPreviousRunNumber;" << endl
                    << "		bool dIsPolarizedFlag; //else is AMO" << endl
                    << "		bool dIsPARAFlag; //else is PERP or AMO" << endl
                    << endl
                    << "		bool dIsMC;" << endl
                    << endl
                    << "		// ANALYZE CUT ACTIONS" << endl
                    << "		// // Automatically makes mass histograms "
                       "where one cut is missing"
                    << endl
                    << "		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;" << endl
                    << endl
                    << "		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS" << endl
                    << endl;

    // print particle step, particle wrapper declarations
    map<int, map<int, pair<Particle_t, string>>>::iterator locStepIterator =
        locComboInfoMap.begin();
    for (; locStepIterator != locComboInfoMap.end(); ++locStepIterator)
    {
        int locStepIndex = locStepIterator->first;
        locHeaderStream << "		//Step " << locStepIndex << endl
                        << "		DParticleComboStep* dStep" << locStepIndex << "Wrapper;"
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
    locHeaderStream << "		// DEFINE YOUR HISTOGRAMS HERE" << endl
                    << "		// EXAMPLES:" << endl
                    << "		TH1I* dHist_MissingMassSquared;" << endl
                    << "		TH1I* dHist_BeamEnergy;" << endl
                    << "		TH1I* dHist_BeamEnergy_BestChiSq;" << endl
                    << endl
                    << "	ClassDef(" << locSelectorName << ", 0);" << endl
                    << "};" << endl
                    << endl
                    << "void " << locSelectorName << "::Get_ComboWrappers(void)"

                    << endl
                    << "{" << endl;

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
    locHeaderStream << "}" << endl << endl << "#endif // " << locSelectorName << "_h" << endl;

    locHeaderStream.close();
}
