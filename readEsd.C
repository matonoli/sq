// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>

using namespace std;
using namespace RooFit;

// GLOBALS
TChain* mChain;
TFile* mFout;

bool makeChain(const Char_t *inputFile="test.list") {

	if (!mChain) mChain = new TChain("esdTree");
	TString inputFileStr(inputFile);
	flagMC = inputFileStr.Contains("MC");
	string const dirFile = inputFileStr.Data();
	if (dirFile.find(".lis") != string::npos)	{
		
		ifstream inputStream(dirFile.c_str());
		if (!inputStream)	{
			cout << "ERROR: Cannot open list file " << dirFile << endl;
			return false;	}

		int nFile = 0;
		string file;
		while (getline(inputStream, file))	{
	  		if (file.find(".root") != string::npos)	{
				TFile* ftmp = TFile::Open(file.c_str());
				if (ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())	{
		  			cout << " Read in V0 tree file " << file << endl;
		  			mChain->Add(file.c_str());
		  			++nFile;	}
				if (ftmp) ftmp->Close();
	  		}
		}

	cout << " Total " << nFile << " files have been read in. " << endl;
	return true;
	}

	else if (dirFile.find(".root") != string::npos)	{
		mChain->Add(dirFile.c_str());	
		return true;	}
  	else	{
		cout << " No good input file to read ... " << endl;
		return false;	}
}

void readEsd(Int_t nEvents=10, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();
	TList* lHist = gDirectory->GetList();
	int iHist = 0;

	if (!makeChain(inputFile)) printf("Couldn't create the chain! \n", );
	else printf("Chain created with %i entries \n", mChain->GetEntries());

	nEvents = (nEvents < mChain->GetEntries()) ? nEvents : mChain->GetEntries();

	AliESDEvent* ev = new AliESDEvent();
	ev->ReadFromTree(mChain);
	AliMCEvent* mcev = new AliMCEvent();
	mcev->ReadFromTree(mChain,"");
	
	mChain->Show(180);

	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		ev->Reset();
		mcev->Reset();
		if (iEv!=180) continue;
		
		//if (iEv%10000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		mChain->GetEntry(iEv);
		//if (!mEvent) continue;
		printf("MCev %i \n", mcev->GetRunNumber());
		AliMultSelection *ams = (AliMultSelection*)ev->FindListObject("MultSelection");

		Float_t vz = ev->GetPrimaryVertexTracks()->GetZ();
		Int_t refmult = AliESDtrackCuts::GetReferenceMultiplicity(ev, AliESDtrackCuts::kTrackletsITSTPC,0.8);
		Float_t v0mult = (ams) ? ams->GetMultiplicityPercentile("V0M") : -999;
		Int_t mcmulti = 0;


		printf("EVENT #%i , VZ is %9.9f , rmult is %i , v0mult is %f \n", iEv, vz, refmult, v0mult);

		
	}	// EVENT LOOP FINISHED


	// WRITING OBJECTS TO OUTPUT FILE
	if (outputFile!="")	mFout = new TFile(outputFile,"RECREATE");
	iHist = 0; while (lHist->At(iHist)) {			// should use an iterator...
		TString objName(lHist->At(iHist)->GetName());
		if (objName.BeginsWith("h")) lHist->At(iHist)->Write();
		iHist++;
	}	// can be replaced with embed->GetHistList()->Write(); ?
}