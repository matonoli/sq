// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
//#include <AliAnalysisPIDTrack.h>

#include <iostream>

using namespace std;

TChain* mChain;
AliAnalysisPIDEvent* mEvent;
TClonesArray* bTracks = 0;
TClonesArray* bV0s = 0;

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
	if (!e->AcceptEvent(1,0)) 				return false;
	if (bV0s->GetEntriesFast() < 1) 		return false;
	if (bTracks->GetEntriesFast() < 1) 		return false;

	return true; 
}

bool SelectV0(AliAnalysisPIDV0* v0) {

	if (fabs(v0->GetEta()) > 0.8) 				return false; //superfluous
	if (v0->GetPt() < 1.0)						return false;
	if (v0->GetPt() > 12.0)						return false;
	//if (v0->GetDCAPV() < 0.05)				return false; // dca of extrapolated v0 to pv? or daughters to pv?
	if (v0->GetDCAV0Daughters() > 0.5)			return false;
	if (v0->GetV0CosinePA() < 0.997)			return false; // :(
	if (v0->GetRadius() < 5.)					return false;
	if (v0->GetRadius() > 100.)					return false;

	if (!SelectV0daughter(v0->GetPosAnalysisTrack()))	return false;
	if (!SelectV0daughter(v0->GetNegAnalysisTrack()))	return false;

	return true;
}

bool SelectV0daughter(AliAnalysisPIDTrack* t) {

	if (fabs(t->GetEta()) > 0.8) return false;

	return true;
}

void readTree_V0(Int_t nEvents=10, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();

	if (!makeChain(inputFile)) printf("Couldn't create the chain! \n", );
	else printf("Chain created with %i entries \n", mChain->GetEntries());

	//*bTracks 	= new TClonesArray("AliAnalysisPIDTrack");
	//*bV0s 		= new TClonesArray("AliAnalysisPIDV0");
	mChain->SetBranchAddress("AnalysisTrack",&bTracks);
	mChain->SetBranchAddress("AnalysisV0Track",&bV0s);
	mChain->SetBranchAddress("AnalysisEvent",&mEvent);

	TH1F* hEventMonitor			= new TH1F("hEventMonitor","",10,-0.5,9.5);
	TH1F* hV0Monitor			= new TH1F("hV0Monitor","",10,-0.5,9.5);

	TH1F* hV0_IMK0s				= new TH1F("hV0_IMK0s","",2000,-1,1);
	TH1F* hV0_IML				= new TH1F("hV0_IML","",2000,-1,1);
	TH1F* hV0_IMAL				= new TH1F("hV0_IMAL","",2000,-1,1);

	nEvents = (nEvents < mChain->GetEntries()) ? nEvents : mChain->GetEntries();

	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		hEventMonitor->Fill(0);
		bTracks->Clear();
		mChain->GetEntry(iEv);
		if (!mEvent) continue;
		hEventMonitor->Fill(1);
		//printf("event vz is %f \n", mEvent->GetVertexZ());

		if (iEv==0) mEvent->PrintEventSelection();
		if (!SelectEvent(mEvent)) continue;
		hEventMonitor->Fill(2);

		Int_t nV0s = bV0s->GetEntriesFast();
		//printf("nV0s is %i \n", nV0s);
		for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
			
			hV0Monitor->Fill(0);
			AliAnalysisPIDV0* v0 = (AliAnalysisPIDV0*)bV0s->At(iV0);
			if (!v0) continue;
			hV0Monitor->Fill(1);

			if (!SelectV0(v0)) continue;
			hV0Monitor->Fill(2);

			hV0_IMK0s->Fill(v0->GetIMK0s());
			hV0_IML->Fill(v0->GetIML());
			hV0_IMAL->Fill(v0->GetIMAL());

			//printf("imko is %f and iml is %f \n",v0->GetIMK0s(),v0->GetIML());
		}

		Int_t nTracks = bTracks->GetEntriesFast();
		//printf("nTracks is %i \n", nTracks);
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
	//hEventMonitor->Draw();
	hV0_IMK0s->Draw();
	new TCanvas;
	hV0_IML->Draw();
	new TCanvas;
	hV0_IMAL->Draw();
}