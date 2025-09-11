#include "DSelector_OmegaKLambdaThrown.h"

void DSelector_OmegaKLambdaThrown::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "OmegaKLambdaThrown.root"; //"" for none
	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning into separate files for AmpTools
	//dOutputTreeFileNameMap["Bin1"] = "mcgen_bin1.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["Bin2"] = "mcgen_bin2.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["Bin3"] = "mcgen_bin3.root"; //key is user-defined, value is output file name

	dFlatTreeFileName = "OmegaKLambdaFlat.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.
	
        //Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit
        dSkipNoTriggerEvents = false;
	dPreviousRunNumber = 0;

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pim_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pim_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pip");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pi0");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_p");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_kplus");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_3pi_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_3pi_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_3pik_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_3pik_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_Lambda_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_Lambda_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimpip_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimpi0_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimk_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimp_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimpip_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimpi0_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimk_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pimp_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pi0pip");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("t_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("tprime_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("t_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("tprime_2");

  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_Kpippi0");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_Lambdapippi0_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_Lambdapippi0_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_KLambda_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_KLambda_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_Lambda_1_omega_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_Lambda_2_omega_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_Lambdapi0_1");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_Lambdapi0_2");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_Kpi0");

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_OmegaKLambdaThrown::Process(Long64_t locEntry)
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

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE



	/******************************************* LOOP OVER THROWN DATA ***************************************/

	//Thrown beam: just use directly
	double locBeamEnergyUsedForBinning = 0.0;
	if(dThrownBeam != NULL)
		locBeamEnergyUsedForBinning = dThrownBeam->Get_P4().E();
        
        TLorentzVector locZeroVec(0,0,0,0);
        TLorentzVector locPiPlusP4,locPiMinus1P4,locPiMinus2P4,locKPlusP4,locProtonP4,locPhoton1P4,locPhoton2P4=locZeroVec;
        TLorentzVector dTargetP4(0,0,0,0.938);
        TLorentzVector locBeamP4 = dThrownBeam->Get_P4();


        UInt_t photonNumber = 0;
        UInt_t pimNumber = 0;
	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);
                
		//Do stuff with the wrapper here ...
		Particle_t locPID = dThrownWrapper->Get_PID();
		TLorentzVector locThrownP4 = dThrownWrapper->Get_P4();
		//cout << "Thrown " << loc_i << ": " << locPID << ", " << locThrownP4.Px() << ", " << locThrownP4.Py() << ", " << locThrownP4.Pz() << ", " << locThrownP4.E() << endl;
	
                if (locPID == 8) locPiPlusP4 = locThrownP4;
                if (locPID == 9) {
                    //if(dThrownWrapper->Get_ParentIndex() < 0){ locPiMinus2P4 = locThrownP4;
                    //}else{locPiMinus1P4 = locThrownP4;}
                    pimNumber++;
                    if(pimNumber == 1) locPiMinus1P4 = locThrownP4;
                    if(pimNumber == 2) locPiMinus2P4 = locThrownP4;
                }
                if (locPID == 11) locKPlusP4 = locThrownP4;
                if (locPID == 14) locProtonP4 = locThrownP4;
                if (locPID == 1) {
                    photonNumber++;
                    if(photonNumber == 1) locPhoton1P4 = locThrownP4;
                    if(photonNumber == 2) locPhoton2P4 = locThrownP4;
                }


        }
        bool pimIsNull = ((locPiMinus1P4 == locZeroVec) || (locPiMinus2P4 == locZeroVec));
        if (!pimIsNull){

	TLorentzVector locMass_3Pi_1 = locPiPlusP4 + locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
	TLorentzVector locMass_3Pi_2 = locPiPlusP4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
	TLorentzVector locMass_Pi0 = locPhoton1P4 + locPhoton2P4;
	TLorentzVector locMass_Lambda_1 = locProtonP4 + locPiMinus1P4;
	TLorentzVector locMass_Lambda_2 = locProtonP4 + locPiMinus2P4;
	TLorentzVector locMass_3PiK_1 = locMass_3Pi_1 + locKPlusP4;
	TLorentzVector locMass_3PiK_2 = locMass_3Pi_2 + locKPlusP4;   


  TLorentzVector locMass_PimPip_1 = locPiPlusP4 + locPiMinus1P4;
  TLorentzVector locMass_PimPi0_1 = locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
  TLorentzVector locMass_Pi0Pip = locPiPlusP4 + locPhoton1P4 + locPhoton2P4;
  TLorentzVector locMass_PimKp_1 = locPiMinus1P4 + locKPlusP4;
  TLorentzVector locMass_PimP_1 = locPiMinus1P4 + locProtonP4;
  TLorentzVector locMass_PimPip_2 = locPiPlusP4 + locPiMinus2P4;
  TLorentzVector locMass_PimPi0_2 = locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
  TLorentzVector locMass_PimKp_2 = locPiMinus2P4 + locKPlusP4;
  TLorentzVector locMass_PimP_2 = locPiMinus2P4 + locProtonP4;

  TLorentzVector locMass_PipPi0K =
      locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locKPlusP4;
  TLorentzVector locMass_PipPi0Lambda1 =
      locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locMass_Lambda_1;
  TLorentzVector locMass_PipPi0Lambda2 =
      locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locMass_Lambda_2;
  TLorentzVector locMass_KpLambda1 = locKPlusP4 + locMass_Lambda_1;
  TLorentzVector locMass_KpLambda2 = locKPlusP4 + locMass_Lambda_2;
  TLorentzVector locMass_Lambda1omega2 = locMass_Lambda_1 + locMass_3Pi_2;
  TLorentzVector locMass_Lambda2omega1 = locMass_Lambda_2 + locMass_3Pi_1;
  TLorentzVector locMass_Pi0Lambda1 =
      locPhoton1P4 + locPhoton2P4 + locMass_Lambda_1;
  TLorentzVector locMass_Pi0Lambda2 =
      locPhoton1P4 + locPhoton2P4 + locMass_Lambda_2;
  TLorentzVector locMass_Pi0K = locPhoton1P4 + locPhoton2P4 + locKPlusP4;

	//Calculating t min for omega_1 k+ system
	TLorentzVector locMomentumTransfer_1 = locBeamP4  - locMass_3PiK_1;//[needs to be squared to be t [.M2()]
	TLorentzVector sMandelstam_1 = locBeamP4 + dTargetP4; //This is not s until .M2() [it needs to be squared]

	double E3CM_1 = ( sMandelstam_1.M2() + locMass_3PiK_1.M2() - locMass_Lambda_2.M2() ) / (2*sMandelstam_1.M());// equation 47.36 pdg book
	double p3CM_1 = sqrt( E3CM_1*E3CM_1 - locMass_3PiK_1.M2()); // equation 47.37
	double p1CM_1 = ( locBeamP4.Vect() ).Mag() * dTargetP4.M() / sMandelstam_1.M() ; // equation 47.37
	double t_1term_1 =  ( locBeamP4.M2() - locMass_3PiK_1.M2() - dTargetP4.M2() + locMass_Lambda_2.M2() )/ (2*sMandelstam_1.M()); //  First term in equation 47.35
	double tmin_1 = (t_1term_1)*(t_1term_1) - (p1CM_1 - p3CM_1)*(p1CM_1 - p3CM_1); // Eq 47.35
	double tprime_1 = (locMomentumTransfer_1.M2() - tmin_1);
	tprime_1 = -1 * tprime_1;

	//Calculating t min for omega_2 k+ system
	TLorentzVector locMomentumTransfer_2 = locBeamP4  - locMass_3PiK_2;//[needs to be squared to be t [.M2()]
	TLorentzVector sMandelstam_2 = locBeamP4 + dTargetP4; //This is not s until .M2() [it needs to be squared]

	double E3CM_2 = ( sMandelstam_2.M2() + locMass_3PiK_2.M2() - locMass_Lambda_1.M2() ) / (2*sMandelstam_2.M());// equation 47.36 pdg book
	double p3CM_2 = sqrt( E3CM_2*E3CM_2 - locMass_3PiK_2.M2()); // equation 47.37
	double p1CM_2 = ( locBeamP4.Vect() ).Mag() * dTargetP4.M() / sMandelstam_2.M() ; // equation 47.37
	double t_1term_2 =  ( locBeamP4.M2() - locMass_3PiK_2.M2() - dTargetP4.M2() + locMass_Lambda_1.M2() )/ (2*sMandelstam_2.M()); //  First term in equation 47.35
	double tmin_2 = (t_1term_2)*(t_1term_2) - (p1CM_2 - p3CM_2)*(p1CM_2 - p3CM_2); // Eq 47.35
	double tprime_2 = (locMomentumTransfer_2.M2() - tmin_2);
	tprime_2 = -1 * tprime_2;

	//OR Manually:
	//BEWARE: Do not expect the particles to be at the same array indices from one event to the next!!!!
	//Why? Because while your channel may be the same, the pions/kaons/etc. will decay differently each event.

	//BRANCHES: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat#TTree_Format:_Simulated_Data
	TClonesArray** locP4Array = dTreeInterface->Get_PointerToPointerTo_TClonesArray("Thrown__P4");
	TBranch* locPIDBranch = dTreeInterface->Get_Branch("Thrown__PID");
/*
	Particle_t locThrown1PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);
	TLorentzVector locThrown1P4 = *((TLorentzVector*)(*locP4Array)->At(0));
	cout << "Particle 1: " << locThrown1PID << ", " << locThrown1P4.Px() << ", " << locThrown1P4.Py() << ", " << locThrown1P4.Pz() << ", " << locThrown1P4.E() << endl;
	Particle_t locThrown2PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);
	TLorentzVector locThrown2P4 = *((TLorentzVector*)(*locP4Array)->At(1));
	cout << "Particle 2: " << locThrown2PID << ", " << locThrown2P4.Px() << ", " << locThrown2P4.Py() << ", " << locThrown2P4.Pz() << ", " << locThrown2P4.E() << endl;
*/


	/******************************************* BIN THROWN DATA INTO SEPARATE TREES FOR AMPTOOLS ***************************************/

/*
	//THESE KEYS MUST BE DEFINED IN THE INIT SECTION (along with the output file names)
	if((locBeamEnergyUsedForBinning >= 8.0) && (locBeamEnergyUsedForBinning < 9.0))
		Fill_OutputTree("Bin1"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 9.0) && (locBeamEnergyUsedForBinning < 10.0))
		Fill_OutputTree("Bin2"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 10.0) && (locBeamEnergyUsedForBinning < 11.0))
		Fill_OutputTree("Bin3"); //your user-defined key
*/


  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_3pi_1",
                                                 locMass_3Pi_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_3pi_2",
                                                 locMass_3Pi_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_3pik_1",
                                                 locMass_3PiK_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_3pik_2",
                                                 locMass_3PiK_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pim_1",
                                                 locPiMinus1P4.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pim_2",
                                                 locPiMinus2P4.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pip", locPiPlusP4.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pi0", locMass_Pi0.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambda_1",
                                                 locMass_Lambda_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambda_2",
                                                 locMass_Lambda_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_kplus", locKPlusP4.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_p", locProtonP4.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimpip_1",
                                                 locMass_PimPip_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimpi0_1",
                                                 locMass_PimPi0_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimk_1",
                                                 locMass_PimKp_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimp_1",
                                                 locMass_PimP_1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimpip_2",
                                                 locMass_PimPip_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimpi0_2",
                                                 locMass_PimPi0_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimk_2",
                                                 locMass_PimKp_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pimp_2",
                                                 locMass_PimP_2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("t_1",
                                                 -locMomentumTransfer_1.M2());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("tprime_1", tprime_1);
  dFlatTreeInterface->Fill_Fundamental<Double_t>("t_2",
                                                 -locMomentumTransfer_2.M2());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("tprime_2", tprime_2);

  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Kpippi0",
                                                 locMass_PipPi0K.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambdapippi0_1",
                                                 locMass_PipPi0Lambda1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambdapippi0_2",
                                                 locMass_PipPi0Lambda2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_KLambda_1",
                                                 locMass_KpLambda1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_KLambda_2",
                                                 locMass_KpLambda2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambda_1_omega_2",
                                                 locMass_Lambda1omega2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambda_2_omega_1",
                                                 locMass_Lambda2omega1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambdapi0_1",
                                                 locMass_Pi0Lambda1.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Lambdapi0_2",
                                                 locMass_Pi0Lambda2.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_Kpi0", locMass_Pi0K.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pi0pip",
                                                 locMass_Pi0Pip.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pi0",
                                                 locMass_Pi0.M());

        Fill_FlatTree();
}
	return kTRUE;
}

void DSelector_OmegaKLambdaThrown::Finalize(void)
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
