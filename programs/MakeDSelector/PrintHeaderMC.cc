#include "MakeDSelector.h"


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
