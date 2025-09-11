#include "DSelector_OmegaKLambdaData.h"

void DSelector_OmegaKLambdaData::Init(TTree *locTree) {
  // USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE"
  // LABEL. LEAVE THE REST ALONE.

  // The Init() function is called when the selector needs to initialize a new
  // tree or chain. Typically here the branch addresses and branch pointers of
  // the tree will be set. Init() will be called many times when running on
  // PROOF (once per file to be processed).

  // USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
  dOutputFileName = "";     //"" for none
  dOutputTreeFileName = ""; //"" for none
  dFlatTreeFileName =
      "OmegaKLambdaDataFlatNoZNew.root"; // output flat tree (one combo per tree
                                         // entry), "" for none
  dFlatTreeName = "";               // if blank, default name will be chosen
  dSaveDefaultFlatBranches = true; // False: don't save default branches,
  // reduce disk footprint. dSaveTLorentzVectorsAsFundamentaFlatTree = false; //
  // Default (or false): save particles as TLorentzVector objects. True: save as
  // four doubles instead.

  // Because this function gets called for each TTree in the TChain, we must be
  // careful: We need to re-initialize the tree interface & branch wrappers, but
  // don't want to recreate histograms
  bool locInitializedPriorFlag =
      dInitializedFlag; // save whether have been initialized previously
  DSelector::Init(
      locTree); // This must be called to initialize wrappers for each new TTree
  // gDirectory now points to the output file with name dOutputFileName (if any)
  if (locInitializedPriorFlag)
    return; // have already created histograms, etc. below: exit

  Get_ComboWrappers();
  dPreviousRunNumber = 0;

  /*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS
   * ACTIONS **********************************/

  // EXAMPLE: Create deque for histogramming particle masses:
  // // For histogramming the phi mass in phi -> K+ K-
  // // Be sure to change this and dAnalyzeCutActions to match reaction
  // std::deque<Particle_t> MyPhi;
  // MyPhi.push_back(KPlus);
  // MyPhi.push_back(KMinus);

  // // ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
  // // false/true below: use measured/kinfit data

  // // PID
  // dAnalysisActions.push_back(
  //     new DHistogramAction_ParticleID(dComboWrapper, false));
  // // below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
  // // dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper,
  // false,
  // // 0.5, KPlus, SYS_BCAL));

  // // PIDFOM (for charged tracks)
  // dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
  // dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus,
  // 0.1)); dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper,
  // 0.1));

  // MASSES
  // dAnalysisActions.push_back(new
  // DHistogramAction_InvariantMass(dComboWrapper, false, Lambda,
  // 1000, 1.0, 1.2, "Lambda")); dAnalysisActions.push_back(new
  // DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1,
  // 0.1));

  // KINFIT RESULTS
  // dAnalysisActions.push_back(new
  // DHistogramAction_KinFitResults(dComboWrapper));

  // CUT MISSING MASS
  // dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper,
  // false, -0.03, 0.02));

  // CUT ON SHOWER QUALITY
  // dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper,
  // SYS_FCAL, 0.5));

  // BEAM ENERGY
  // dAnalysisActions.push_back(
  //  new DHistogramAction_BeamEnergy(dComboWrapper, false));
  // dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper,
  // false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999
  // KINEMATICS
  // dAnalysisActions.push_back(
  //  new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
  // ANALYZE CUT ACTIONS
  // // Change MyPhi to match reaction
  // dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions(
  //  dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4,
  //"CutActionEffect");

  // INITIALIZE ACTIONS
  // If you create any actions that you want to run manually (i.e. don't add to
  // dAnalysisActions), be sure to initialize them here as well
  Initialize_Actions();
  // dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

  /******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE
   * HISTOGRAMS *******************************/

  // EXAMPLE MANUAL HISTOGRAMS:
  // dHist_MissingMassSquared =
  //  new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}",
  //         600, -0.06, 0.06);
  // dHist_BeamEnergy =
  //   new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

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
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mass_pi0_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_Lambda_1_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_Lambda_2_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_3pik_1_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_3pik_2_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_3pi_1_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "mass_3pi_2_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("missing_mass_sq");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("missing_energy");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "momentum_p_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "momentum_pip_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "momentum_kp_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "dEdx_cdc_p_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "dEdx_cdc_kp_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "dEdx_cdc_pip_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "dEdx_fdc_p_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "dEdx_fdc_kp_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>(
      "dEdx_fdc_pip_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("beta_p_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("beta_kp_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("beta_pip_measured");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("clevel");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("chisq");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("beam_time_delta");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("lambda_param");

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
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("detached_vertex");
  dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("path_length");
  dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("path_length_sigma");

  /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT
   * BRANCHES - MAIN TREE *************************/

  // EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE
  // GIVEN!!!! (ABOVE: TOP)): The type for the branch must be included in the
  // brackets 1st function argument is the name of the branch 2nd function
  // argument is the name of the branch that contains the size of the array (for
  // fundamentals only)
  /*
  dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental =
  char, int, float, double, etc.
  dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array",
  "my_int");
  dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array",
  "NumCombos");
  dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
  dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
  */

  /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT
   * BRANCHES - FLAT TREE *************************/

  // RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
  // dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");

  // EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE
  // GIVEN!!!! (ABOVE: TOP)): The type for the branch must be included in the
  // brackets 1st function argument is the name of the branch 2nd function
  // argument is the name of the branch that contains the size of the array (for
  // fundamentals only)
  /*
  dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int");
  //fundamental = char, int, float, double, etc.
  dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array",
  "flat_my_int");
  dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
  dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
  */

  /************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO
   * READ ************************************/

  // TO SAVE PROCESSING TIME
  // If you know you don't need all of the branches/data, but just a subset of
  // it, you can speed things up By default, for each event, the data is
  // retrieved for all branches If you know you only need data for some
  // branches, you can skip grabbing data from the branches you don't need Do
  // this by doing something similar to the commented code below

  // dTreeInterface->Clear_GetEntryBranches(); //now get none
  // dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the
  // branches you want

  /************************************** DETERMINE IF ANALYZING SIMULATED DATA
   * *************************************/

  dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);
}

Bool_t DSelector_OmegaKLambdaData::Process(Long64_t locEntry) {
  // The Process() function is called for each entry in the tree. The entry
  // argument specifies which entry in the currently loaded tree is to be
  // processed.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  // Use fStatus to set the return value of TTree::Process().
  // The return value is currently not used.

  // CALL THIS FIRST
  DSelector::Process(locEntry); // Gets the data from the tree for the entry
  // cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() <<
  // endl; TLorentzVector locProductionX4 = Get_X4_Production();

  /******************************************** GET POLARIZATION ORIENTATION
   * ******************************************/

  // Only if the run number changes
  // RCDB environment must be setup in order for this to work! (Will return
  // false otherwise)
  UInt_t locRunNumber = Get_RunNumber();
  if (locRunNumber != dPreviousRunNumber) {
    dIsPolarizedFlag =
        dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
    dPreviousRunNumber = locRunNumber;
  }

  /********************************************* SETUP UNIQUENESS TRACKING
   * ********************************************/

  // ANALYSIS ACTIONS: Reset uniqueness tracking for each action
  // For any actions that you are executing manually, be sure to call
  // Reset_NewEvent() on them here
  Reset_Actions_NewEvent();
  // dAnalyzeCutActions
  //   ->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

  // PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
  // Sometimes, some content is the exact same between one combo and the next
  // e.g. maybe two combos have different beam particles, but the same data for
  // the final-state
  // When histogramming, you don't want to double-count when this happens:
  // artificially inflates your signal (or background) So, for each quantity you
  // histogram, keep track of what particles you used (for a given combo) Then
  // for each combo, just compare to what you used before, and make sure it's
  // unique

  // EXAMPLE 1: Particle-specific info:
  set<Int_t> locUsedSoFar_BeamEnergy; // Int_t: Unique ID for beam particles.
                                      // set: easy to use, fast to search

  // EXAMPLE 2: Combo-specific info:
  // In general: Could have multiple particles with the same PID: Use a set of
  // Int_t's In general: Multiple PIDs, so multiple sets: Contain within a map
  // Multiple combos: Contain maps within a set (easier, faster to search)
  set<map<Particle_t, set<Int_t>>> locUsedSoFar_MissingMass;

  // INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

  /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES
   * **************************************/

  /*
  Int_t locMyInt = 7;
  dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

  TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
  dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

  for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
          dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i,
  loc_i); //2nd argument = value, 3rd = array index
  */

  /************************************************* LOOP OVER COMBOS
   * *************************************************/
  double best_chisq = 1000000000;
  int best_combo = -1;
  // Loop over combos
  for (UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
    // Set branch array indices for combo and all combo particles
    dComboWrapper->Set_ComboIndex(loc_i);

    // Is used to indicate when combos have been cut
    if (dComboWrapper
            ->Get_IsComboCut()) // Is false when tree originally created
      continue;                 // Combo has been cut previously

    /********************************************** GET PARTICLE INDICES
     * *********************************************/
    double locRFTime = dComboWrapper->Get_RFTime();
    TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();

    // Beam
    double locPropagatedRFTimeBeam =
        locRFTime + (locBeamX4_Measured.Z() - dTargetCenter.Z()) / 29.9792458;
    double locBeamDeltaT = locBeamX4_Measured.T() - locPropagatedRFTimeBeam;

    double chisq = dComboWrapper->Get_ChiSq_KinFit();
    double ndof = dComboWrapper->Get_NDF_KinFit();
    // Used for tracking uniqueness when filling histograms, and for determining
    // unused particles
    if (chisq / ndof < best_chisq) {
      best_chisq = chisq / ndof;
      best_combo = loc_i;
    }
  }

  if (best_combo == -1) {
    return kTRUE;
  }
  dComboWrapper->Set_ComboIndex(best_combo);
  // Step 0
  Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
  Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
  Int_t locPiMinus1TrackID = dPiMinus1Wrapper->Get_TrackID();
  Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();

  // Step 1
  Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
  Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

  // Step 2
  Int_t locPiMinus2TrackID = dPiMinus2Wrapper->Get_TrackID();
  Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

  /*********************************************** GET FOUR-MOMENTUM
   * **********************************************/

  // Get P4's: //is kinfit if kinfit performed, else is measured
  // dTargetP4 is target p4
  // Step 0
  TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
  TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
  TLorentzVector locPiMinus1P4 = dPiMinus1Wrapper->Get_P4();
  TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
  // Step 1
  TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
  TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
  TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
  // Step 2
  TLorentzVector locPiMinus2P4 = dPiMinus2Wrapper->Get_P4();
  TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();

  // Get Measured P4's:
  // Step 0
  TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
  TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
  TLorentzVector locPiMinus1P4_Measured = dPiMinus1Wrapper->Get_P4_Measured();
  TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
  // Step 1
  TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
  TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
  // Step 2
  TLorentzVector locPiMinus2P4_Measured = dPiMinus2Wrapper->Get_P4_Measured();
  TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();

  /********************************************* GET COMBO RF TIMING INFO
   * *****************************************/

  TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
  TLorentzVector locProtonX4_Measured = dProtonWrapper->Get_X4_Measured();
  // Double_t locBunchPeriod =
  // dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
  //  Double_t locDeltaT_RF =
  //  dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured,
  //  dComboWrapper); Int_t locRelBeamBucket =
  //  dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(),
  //  locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero
  //  integer for out-of-time photons Int_t locNumOutOfTimeBunchesInTree =
  //  XXX; //YOU need to specify this number
  // Number of out-of-time beam bunches in tree (on a single side, so that
  // total number out-of-time bunches accepted is 2 times this number for left
  // + right bunches)

  // Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from
  // nearest out-of-time bunch on either side (recommended). Int_t
  // locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ?
  // locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; Double_t
  // locAccidentalScalingFactor =
  // dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(),
  // locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require
  // added factor, which is different for data and MC. Double_t
  // locAccidentalScalingFactorError =
  // dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(),
  // locBeamP4.E()); // Ideal value would be 1, but deviations observed, need
  // added factor. Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1
  // : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight
  // by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
  // if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip
  // nearest out-of-time bunch: tails of in-time distribution also leak in
  // 	dComboWrapper->Set_IsComboCut(true);
  // 	continue;
  // }

  /********************************************* COMBINE FOUR-MOMENTUM
   * ********************************************/

  // DO YOUR STUFF HERE

  double locPathLength = -1;
  DParticleComboStep *locStep = dStep2Wrapper;
  Particle_t locInitialPID = locStep->Get_InitialPID();
  TLorentzVector locStepSpacetimeVertex = dStep2Wrapper->Get_X4();
  int locFromStepIndex = dStep2Wrapper->Get_InitDecayFromIndices().first;
  TLorentzVector locFromSpacetimeVertex = dStep0Wrapper->Get_X4();
  TLorentzVector locDeltaSpacetime =
      locStepSpacetimeVertex - locFromSpacetimeVertex;

  // dHistMap_DetachedPathLength[loc_i]->Fill(locPathLength);
  // dHistMap_DetachedLifetime[loc_i]->Fill(locDeltaSpacetime.T());
  locPathLength = locDeltaSpacetime.Vect().Mag();


  double dEdx_cdc_p = dProtonWrapper->Get_dEdx_CDC();
  double dEdx_cdc_pip = dPiPlusWrapper->Get_dEdx_CDC();
  double dEdx_cdc_kp = dKPlusWrapper->Get_dEdx_CDC();
  double dEdx_fdc_p = dProtonWrapper->Get_dEdx_FDC();
  double dEdx_fdc_pip = dPiPlusWrapper->Get_dEdx_FDC();
  double dEdx_fdc_kp = dKPlusWrapper->Get_dEdx_FDC();
  double beta_p = dProtonWrapper->Get_Beta_Timing_Measured();
  double beta_pip = dPiPlusWrapper->Get_Beta_Timing_Measured();
  double beta_kp = dKPlusWrapper->Get_Beta_Timing_Measured();
  double clevel = dComboWrapper->Get_ConfidenceLevel_KinFit("");
  double chisq = dComboWrapper->Get_ChiSq_KinFit("");
  double ndof = dComboWrapper->Get_NDF_KinFit("");

  double locRFTime = dComboWrapper->Get_RFTime();
  // Beam
  double locPropagatedRFTimeBeam =
      locRFTime + (locBeamX4_Measured.Z() - dTargetCenter.Z()) / 29.9792458;
  double locBeamDeltaT = locBeamX4_Measured.T() - locPropagatedRFTimeBeam;

  // Combine 4-vectors
  TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
  locMissingP4_Measured -= locPiPlusP4_Measured + locPiMinus1P4_Measured +
                           locKPlusP4_Measured + locPhoton1P4_Measured +
                           locPhoton2P4_Measured + locPiMinus2P4_Measured +
                           locProtonP4_Measured;

  TLorentzVector locMass_3Pi_1 =
      locPiPlusP4 + locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
  TLorentzVector locMass_3Pi_2 =
      locPiPlusP4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
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

  TLorentzVector locMass_3Pi_1_Measured =
      locPiPlusP4_Measured + locPiMinus1P4_Measured + locPhoton1P4_Measured +
      locPhoton2P4_Measured;
  TLorentzVector locMass_3Pi_2_Measured =
      locPiPlusP4_Measured + locPiMinus2P4_Measured + locPhoton1P4_Measured +
      locPhoton2P4_Measured;
  TLorentzVector locMass_Pi0_Measured =
      locPhoton1P4_Measured + locPhoton2P4_Measured;
  TLorentzVector locMass_Lambda_1_Measured =
      locProtonP4_Measured + locPiMinus1P4_Measured;
  TLorentzVector locMass_Lambda_2_Measured =
      locProtonP4_Measured + locPiMinus2P4_Measured;
  TLorentzVector locMass_3PiK_1_Measured =
      locMass_3Pi_1_Measured + locKPlusP4_Measured;
  TLorentzVector locMass_3PiK_2_Measured =
      locMass_3Pi_2_Measured + locKPlusP4_Measured;

  // Calculating t min for omega_1 k+ system
  TLorentzVector locMomentumTransfer_1 =
      locBeamP4 - locMass_3PiK_1; //[needs to be squared to be t [.M2()]
  TLorentzVector sMandelstam_1 =
      locBeamP4 +
      dTargetP4; // This is not s until .M2() [it needs to be squared]

  double E3CM_1 =
      (sMandelstam_1.M2() + locMass_3PiK_1.M2() - locMass_Lambda_2.M2()) /
      (2 * sMandelstam_1.M()); // equation 47.36 pdg book
  double p3CM_1 = sqrt(E3CM_1 * E3CM_1 - locMass_3PiK_1.M2()); // equation 47.37
  double p1CM_1 = (locBeamP4.Vect()).Mag() * dTargetP4.M() /
                  sMandelstam_1.M(); // equation 47.37
  double t_1term_1 = (locBeamP4.M2() - locMass_3PiK_1.M2() - dTargetP4.M2() +
                      locMass_Lambda_2.M2()) /
                     (2 * sMandelstam_1.M()); //  First term in equation 47.35
  double tmin_1 = (t_1term_1) * (t_1term_1) -
                  (p1CM_1 - p3CM_1) * (p1CM_1 - p3CM_1); // Eq 47.35
  double tprime_1 = (locMomentumTransfer_1.M2() - tmin_1);
  tprime_1 = -1 * tprime_1;

  // Calculating t min for omega_2 k+ system
  TLorentzVector locMomentumTransfer_2 =
      locBeamP4 - locMass_3PiK_2; //[needs to be squared to be t [.M2()]
  TLorentzVector sMandelstam_2 =
      locBeamP4 +
      dTargetP4; // This is not s until .M2() [it needs to be squared]

  double E3CM_2 =
      (sMandelstam_2.M2() + locMass_3PiK_2.M2() - locMass_Lambda_1.M2()) /
      (2 * sMandelstam_2.M()); // equation 47.36 pdg book
  double p3CM_2 = sqrt(E3CM_2 * E3CM_2 - locMass_3PiK_2.M2()); // equation 47.37
  double p1CM_2 = (locBeamP4.Vect()).Mag() * dTargetP4.M() /
                  sMandelstam_2.M(); // equation 47.37
  double t_1term_2 = (locBeamP4.M2() - locMass_3PiK_2.M2() - dTargetP4.M2() +
                      locMass_Lambda_1.M2()) /
                     (2 * sMandelstam_2.M()); //  First term in equation 47.35
  double tmin_2 = (t_1term_2) * (t_1term_2) -
                  (p1CM_2 - p3CM_2) * (p1CM_2 - p3CM_2); // Eq 47.35
  double tprime_2 = (locMomentumTransfer_2.M2() - tmin_2);
  tprime_2 = -1 * tprime_2;

  // Boost  FROM: LAB  ---> TO: CM Frame
  TLorentzVector centerOfMass = locBeamP4 + dTargetP4;
  TLorentzRotation centerofMassBoost(-centerOfMass.BoostVector());

  TLorentzVector beamCM = centerofMassBoost * locBeamP4;
  TLorentzVector targetCM = centerofMassBoost * dTargetP4;
  TLorentzVector piPlusCM = centerofMassBoost * locPiPlusP4;
  TLorentzVector piMinus1CM = centerofMassBoost * locPiMinus1P4;
  TLorentzVector piMinus2CM = centerofMassBoost * locPiMinus2P4;
  TLorentzVector KpCM = centerofMassBoost * locKPlusP4;
  TLorentzVector protonCM = centerofMassBoost * locProtonP4;
  TLorentzVector photon1CM = centerofMassBoost * locPhoton1P4;
  TLorentzVector photon2CM = centerofMassBoost * locPhoton2P4;
  TLorentzVector particleX1CM =
      KpCM + piPlusCM + piMinus1CM + photon1CM + photon2CM;
  // Boost!  FROM: CM Frame ---> TO: X Rest Frame (GJ)
  TLorentzRotation restFrameXBoost(-particleX1CM.BoostVector());

  TLorentzVector beamGJ = restFrameXBoost * beamCM;
  TLorentzVector targetGJ = restFrameXBoost * targetCM;
  TLorentzVector protonGJ = restFrameXBoost * protonCM;
  TLorentzVector KpGJ = restFrameXBoost * locKPlusP4;
  TLorentzVector piPlusGJ = restFrameXBoost * piPlusCM;
  TLorentzVector piMinus1GJ = restFrameXBoost * piMinus1CM;
  TLorentzVector piMinus2GJ = restFrameXBoost * piMinus2CM;
  TLorentzVector photon1GJ = restFrameXBoost * photon1CM;
  TLorentzVector photon2GJ = restFrameXBoost * photon2CM;
  TLorentzVector threePiGJ1 = piPlusGJ + piMinus1GJ + photon1GJ + photon2GJ;

  // Boost!  FROM: X Rest Frame ---> TO: Omega1 Rest Frame
  TLorentzRotation restFrameOGBoost(-threePiGJ1.BoostVector());
  TLorentzVector beamOG = restFrameOGBoost * beamGJ;
  TLorentzVector piPlusOG = restFrameOGBoost * piPlusGJ;
  TLorentzVector piMinus1OG = restFrameOGBoost * piMinus1GJ;
  TLorentzVector photon1OG = restFrameOGBoost * photon1GJ;
  TLorentzVector photon2OG = restFrameOGBoost * photon2GJ;

  // Calculating lambda param for omega1
  TVector3 piZeroOG3V = photon1OG.Vect() + photon2OG.Vect();
  TVector3 piPlusOG3V = piPlusOG.Vect();
  TVector3 piMinusOG13V = piMinus1OG.Vect();

  double piZeroOG3VMag = piZeroOG3V.Mag();

  TLorentzVector m3Pions = photon1OG + photon2OG + piPlusOG + piMinus1OG;

  double denominator_lambda = ((1. / 9.) * m3Pions.M2()) - (.140 * .140);
  double lambdaMax = (3. / 4.) * (denominator_lambda * denominator_lambda);

  TVector3 lambdaPions = piPlusOG3V.Cross(piMinusOG13V);

  double ratioLambda = lambdaPions.Mag2() / lambdaMax;
  Float_t PathLengthSigma = Get_Fundamental<Float_t>("DecayingLambda__PathLengthSigma", best_combo);

  /******************************************** EXECUTE ANALYSIS ACTIONS
   * *******************************************/

  // Loop through the analysis actions, executing them in order for the active
  // particle combo
  // dAnalyzeCutActions
  //  ->Perform_Action(); // Must be executed before Execute_Actions()
  if (!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag
                          // automatically set
    return kTRUE;

  // if you manually execute any actions, and it fails a cut, be sure to call:
  // dComboWrapper->Set_IsComboCut(true);

  bool confidenceLevel = dComboWrapper->Get_ConfidenceLevel_KinFit("") < 1e-5;

  // bool omegaWindow = (locMass_3Pi_1.M() > .66 && locMass_3Pi_1.M() < .99)
  // || (locMass_3Pi_2.M() > .66 && locMass_3Pi_2.M() < .99);
  //  suggested in the FCAL quality paper. This could be studied in detail if
  //  needed
  bool showerQuality_Photon1 = false;
  bool showerQuality_Photon2 = false;

  if (dPhoton1Wrapper->Get_Detector_System_Timing() == SYS_FCAL)
    showerQuality_Photon1 = dPhoton1Wrapper->Get_Shower_Quality() < .5;
  if (dPhoton2Wrapper->Get_Detector_System_Timing() == SYS_FCAL)
    showerQuality_Photon2 = dPhoton2Wrapper->Get_Shower_Quality() < .5;
  bool zVertex =
      locProtonX4_Measured.Z() < 52 ||
      locProtonX4_Measured.Z() > 78; // Used in other GlueX analysis, //NOTE:
                                     // improper to cut on proton zVertex
  // it might be the standard
  if (confidenceLevel || showerQuality_Photon1 ||
      showerQuality_Photon2 /*|| zVertex*/) {
    dComboWrapper->Set_IsComboCut(true);
    return kTRUE;
  }
  /**************************************** EXAMPLE: FILL CUSTOM OUTPUT
   * BRANCHES **************************************/

  /*
  TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
  //for arrays below: 2nd argument is value, 3rd is array index
  //NOTE: By filling here, AFTER the cuts above, some indices won't be updated
  (and will be whatever they were from the last event)
          //So, when you draw the branch, be sure to cut on "IsComboCut" to
  avoid these. dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array",
  -2*loc_i, loc_i);
  dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4,
  loc_i);
  */

  /**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY
   * *****************************************/

  // Histogram beam energy (if haven't already)
  /*if (locUsedSoFar_BeamEnergy.find(locBeamID) ==
      locUsedSoFar_BeamEnergy.end()) {
    dHist_BeamEnergy->Fill(
        locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
    // dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); //
    // Alternate version with accidental subtraction

    locUsedSoFar_BeamEnergy.insert(locBeamID);
  }*/

  /************************************ EXAMPLE: HISTOGRAM MISSING MASS
   * SQUARED ************************************/

  // Missing Mass Squared
  double locMissingMassSquared = locMissingP4_Measured.M2();

  // Uniqueness tracking: Build the map of particles used for the missing mass
  // For beam: Don't want to group with final-state photons. Instead use
  // "Unknown" PID (not ideal, but it's easy).
  map<Particle_t, set<Int_t>> locUsedThisCombo_MissingMass;
  locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); // beam
  locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
  locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus1TrackID);
  locUsedThisCombo_MissingMass[KPlus].insert(locKPlusTrackID);
  locUsedThisCombo_MissingMass[Gamma].insert(locPhoton1NeutralID);
  locUsedThisCombo_MissingMass[Gamma].insert(locPhoton2NeutralID);
  locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus2TrackID);
  locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);

  // compare to what's been used so far
  /*if (locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) ==
      locUsedSoFar_MissingMass.end()) {
    // unique missing mass combo: histogram it, and register this combo of
    // particles
    dHist_MissingMassSquared->Fill(
        locMissingMassSquared); // Fills in-time and out-of-time beam photon
                                // combos
    //
  dHist_MissingMassSquared->Fill(locMissingMassSquared,locHistAccidWeightFactor);
    // // Alternate version with accidental subtraction

    locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
  }*/

  // E.g. Cut
  // if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
  //{
  //	dComboWrapper->Set_IsComboCut(true);
  //	continue;
  // }

  /****************************************** FILL FLAT TREE (IF DESIRED)
   * ******************************************/

  // RECOMMENDED: FILL ACCIDENTAL WEIGHT
  // dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight",locHistAccidWeightFactor);

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
  dFlatTreeInterface->Fill_Fundamental<Double_t>("mass_pi0_measured",
                                                 locMass_Pi0_Measured.M());
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "mass_3pi_1_measured",
      locMass_3Pi_1_Measured.M()); // mass_Lambda_1_measured
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "mass_3pi_2_measured",
      locMass_3Pi_2_Measured.M()); // mass_Lambda_2_measured
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "mass_3pik_1_measured",
      locMass_3PiK_1_Measured.M()); // mass_3pik_1_measured")
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "mass_3pik_2_measured",
      locMass_3PiK_2_Measured.M()); // mass_3pik_1_measured")
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "mass_Lambda_1_measured",
      locMass_Lambda_1_Measured.M()); // mass_3pi_1_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "mass_Lambda_2_measured",
      locMass_Lambda_2_Measured.M()); // mass_3pi_2_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "missing_mass_sq", locMissingP4_Measured.M2()); // missing_mass_sq");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "missing_energy", locMissingP4_Measured.E()); // missing_energy");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "momentum_p_measured",
      locProtonP4_Measured.Vect().Mag()); // momentum_p_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "momentum_pip_measured",
      locPiPlusP4_Measured.Vect().Mag()); // momentum_pip_measured"
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "momentum_kp_measured",
      locKPlusP4_Measured.Vect().Mag()); // momentum_kp_measured")
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "dEdx_cdc_p_measured", dEdx_cdc_p); // dEdx_cdc_p_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "dEdx_cdc_kp_measured", dEdx_cdc_kp); // dEdx_cdc_kp_measured")
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "dEdx_cdc_pip_measured", dEdx_cdc_pip); // dEdx_cdc_pip_measured"
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "dEdx_fdc_p_measured", dEdx_fdc_p); // dEdx_fdc_p_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "dEdx_fdc_kp_measured", dEdx_fdc_kp); // dEdx_fdc_kp_measured")
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "dEdx_fdc_pip_measured", dEdx_fdc_pip); // dEdx_fdc_pip_measured"
  dFlatTreeInterface->Fill_Fundamental<Double_t>("beta_p_measured",
                                                 beta_p); // beta_p_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "beta_kp_measured", beta_kp); // beta_kp_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "beta_pip_measured", beta_pip); // beta_pip_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "clevel", clevel); // beta_pip_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "chisq", chisq / ndof); // beta_pip_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>(
      "beam_time_delta", locBeamDeltaT); // beta_pip_measured");
  dFlatTreeInterface->Fill_Fundamental<Double_t>("lambda_param", ratioLambda);
  dFlatTreeInterface->Fill_Fundamental<Double_t>("detached_vertex",
                                                 locStepSpacetimeVertex.Z());
  dFlatTreeInterface->Fill_Fundamental<Double_t>("path_length", locPathLength);
  dFlatTreeInterface->Fill_Fundamental<Float_t>("path_length_sigma", PathLengthSigma);
  /*
  //FILL ANY CUSTOM BRANCHES FIRST!!
  Int_t locMyInt_Flat = 7;
  dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

  TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
  dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4",
  locMyP4_Flat);

  for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
  {
          dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array",
  3*loc_j, loc_j); //2nd argument = value, 3rd = array index TLorentzVector
  locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
          dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array",
  locMyComboP4_Flat, loc_j);
  }
  */

  // FILL FLAT TREE
  Fill_FlatTree(); // for the active combo
                   // end of combo loop

  // FILL HISTOGRAMS: Num combos / events surviving actions
  // Fill_NumCombosSurvivedHists();

  /******************************************* LOOP OVER THROWN DATA (OPTIONAL)
   * ***************************************/
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
  /****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL)
   * ***************************************/
  /*
          //Loop over beam particles (note, only those appearing in combos are
     present) for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
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

  /************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH
   * CUTS APPLIED ************************************/
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

void DSelector_OmegaKLambdaData::Finalize(void) {
  // Save anything to output here that you do not want to be in the default
  // DSelector output ROOT file.

  // Otherwise, don't do anything else (especially if you are using PROOF).
  // If you are using PROOF, this function is called on each thread,
  // so anything you do will not have the combined information from the various
  // threads. Besides, it is best-practice to do post-processing (e.g. fitting)
  // separately, in case there is a problem.

  // DO YOUR STUFF HERE

  // CALL THIS LAST
  DSelector::Finalize(); // Saves results to the output file
}
