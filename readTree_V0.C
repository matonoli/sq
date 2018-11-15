// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
//#include <AliAnalysisPIDTrack.h>

#include <iostream>

using namespace std;

TChain* mChain;
AliAnalysisPIDEvent* mEvent;
TClonesArray* bTracks;
TClonesArray* bV0s;

bool makeChain(const Char_t *inputFile="test.list") {

	if (!mChain) mChain = new TChain("PIDTree");
	TString inputFileStr(inputFile);
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

bool SelectEvent(AliAnalysisPIDEvent* e) {

	//cout << "ae is " << e->AcceptEvent(1,1) << endl;
	if (!e->AcceptEvent(1,0)) 			return false;
	if (bV0s->GetEntriesFast() < 1) 	return false;
	return true; 

}

void readTree_V0(Int_t nEvents=10, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();

	if (!makeChain(inputFile)) printf("Couldn't create the chain! \n", );
	else printf("Chain created with %i entries \n", mChain->GetEntries());

	*bTracks 	= new TClonesArray("AliAnalysisPIDTrack");
	*bV0s 		= new TClonesArray("AliAnalysisPIDV0");
	mChain->SetBranchAddress("AnalysisTrack",&bTracks);
	mChain->SetBranchAddress("AnalysisV0Track",&bV0s);
	mChain->SetBranchAddress("AnalysisEvent",&mEvent);

	TH1F* hEventMonitor			= new TH1F("hEventMonitor","",10,-0.5,9.5);

	nEvents = (nEvents < mChain->GetEntries()) ? nEvents : mChain->GetEntries();

	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		hEventMonitor->Fill(0);
		bTracks->Clear();
		mChain->GetEntry(iEv);
		//printf("event vz is %f \n", mEvent->GetVertexZ());

		if (iEv==0) mEvent->PrintEventSelection();
		//precut histos
		if (!SelectEvent(mEvent)) continue;
		//aftercut histos
		hEventMonitor->Fill(1);

		Int_t nTracks = bTracks->GetEntriesFast();
		printf("nTracks is %i \n", nTracks);
		for (int iTr = 0; iTr < nTracks; ++iTr)	{
			AliAnalysisPIDTrack* track = (AliAnalysisPIDTrack*)bTracks->At(iTr);
			//printf("pt is %f \n", track->GetPt());

			//precut track histos
			//select track
			//aftercut track histos
		}

		//v0loop

		
	}

	printf(" WHAT IS UP \n", );
	hEventMonitor->Draw();
}