#include "MakeDSelector.h"


void Print_HeaderFile_MCGen(string locSelectorBaseName, DTreeInterface *locTreeInterface)
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
                    << "		void Finalize(void);" << endl
                    << endl
                    << "		// BEAM POLARIZATION INFORMATION" << endl
                    << "		UInt_t dPreviousRunNumber;" << endl
                    << "		bool dIsPolarizedFlag; //else is AMO" << endl
                    << "		bool dIsPARAFlag; //else is PERP or AMO" << endl
                    << endl
                    << "	ClassDef(" << locSelectorName << ", 0);" << endl
                    << "};" << endl
                    << endl
                    << "#endif // " << locSelectorName << "_h" << endl;

    locHeaderStream.close();
}
