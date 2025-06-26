#include "MakeDSelector.h"
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
            "(-1D mass plots, 2D correlation plots, Dalitz plots, angular distributions)"
         << endl
         << "  input_root_file     : ROOT file containing the TTree" << endl
         << "  tree_name           : name of TTree inside ROOT file" << endl
         << "  selector_base_name  : base name for generated DSelector files" << endl
         << endl;
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


