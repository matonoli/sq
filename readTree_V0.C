// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>

using namespace std;
using namespace RooFit;

class BeforeMain {
	public:
		BeforeMain()	{ TH1::SetDefaultSumw2(); }	};
BeforeMain sumw2;

// GLOBALS
TChain* mChain;
AliAnalysisPIDEvent* mEvent;
TFile* mFout;
TClonesArray* bTracks = 0;
TClonesArray* bV0s = 0;
TClonesArray* bParticles = 0;
TNtuple* mMasses = new TNtuple("mMasses","v0 masses","MassDT:lPt:lPart");
Int_t flagMC = 0;


const Int_t PDG_IDS[] = {310, 3122, -3122};

TCanvas* cFits[6];
int canCounter = 0;
TString cNames[] = {"K^{0}_{s}", "K^{0}_{s} ub", "#Lambda", "#Lambda ub", "#bar{#Lambda}", "#bar{#Lambda} ub"};
const Int_t nPtBins = 35;
Double_t xBins[nPtBins+1] = { 0.00, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
	1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 
	3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 
	9.00, 10.00, 11.00, 12.00, 13.00, 14.00 };

// HISTOGRAMS
TH1F* hEventMonitor				= new TH1F("hEventMonitor","",10,-0.5,9.5);		// monitors
TH1F* hV0Monitor				= new TH1F("hV0Monitor","",10,-0.5,9.5);
TH1F* hTrackMonitor				= new TH1F("hTrackMonitor","",10,-0.5,9.5);
TH2F* hV0TrCounter				= new TH2F("hV0TrCounter","",100,0,100,100,0,100);

TH1F* hV0_IMK0s					= new TH1F("hV0_IMK0s","",2000,-0.2,0.2);		// v0 stats
TH1F* hV0_IML					= new TH1F("hV0_IML","",2000,-0.2,0.2);
TH1F* hV0_IMAL					= new TH1F("hV0_IMAL","",2000,-0.2,0.2);
TH1F* hV0_PtK0s					= new TH1F("hV0_PtK0s","",nPtBins,xBins);
TH1F* hV0_PtL					= new TH1F("hV0_PtL","",nPtBins,xBins);
TH1F* hV0_PtAL					= new TH1F("hV0_PtAL","",nPtBins,xBins);
TH2F* hV0_IMPtK0s				= new TH2F("hV0_IMPtK0s","",2000,-0.2,0.2,nPtBins,xBins);
TH2F* hV0_IMPtL					= new TH2F("hV0_IMPtL","",2000,-0.2,0.2,nPtBins,xBins);
TH2F* hV0_IMPtAL				= new TH2F("hV0_IMPtAL","",2000,-0.2,0.2,nPtBins,xBins);
TH1F* hYieldK0s					= new TH1F("hYieldK0s","",nPtBins,xBins);
TH1F* hYieldL					= new TH1F("hYieldL","",nPtBins,xBins);
TH1F* hYieldAL					= new TH1F("hYieldAL","",nPtBins,xBins);
TH1F* hYieldUBK0s				= new TH1F("hYieldUBK0s","",nPtBins,xBins);
TH1F* hYieldUBL					= new TH1F("hYieldUBL","",nPtBins,xBins);
TH1F* hYieldUBAL				= new TH1F("hYieldUBAL","",nPtBins,xBins);
TH1F* hV0_DPt 					= new TH1F("hV0_DPt","",200,0,10);
TH2F* hV0_Radiusvpt				= new TH2F("hV0_Radiusvpt","",300,0,15,400,0,100);
TH3F* hV0_PtPhiEta				= new TH3F("hV0_PtPhiEta","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH2F* hV0_AP 					= new TH2F("hV0_AP","",400,-1.,1.,400,0.,0.4);

TH1F* hSBinK0s					= new TH1F("hSBinK0s","",nPtBins,xBins);		// sideband stuff
TH1F* hSBinL					= new TH1F("hSBinL","",nPtBins,xBins);
TH1F* hSBinAL					= new TH1F("hSBinAL","",nPtBins,xBins);
TH1F* hSBoutK0s					= new TH1F("hSBoutK0s","",nPtBins,xBins);
TH1F* hSBoutL					= new TH1F("hSBoutL","",nPtBins,xBins);
TH1F* hSBoutAL					= new TH1F("hSBoutAL","",nPtBins,xBins);

TH1F* hV0_DHasTPC				= new TH1F("hV0_DHasTPC","",200,0,10);			// pid qa stuff
TH1F* hV0_DHasTOF				= new TH1F("hV0_DHasTOF","",200,0,10);
TH2F* hV0_DDTofPiPi				= new TH2F("hV0_DDTofPiPi","",300,-15,15,300,-15,15);
TH2F* hV0_DDTofPiP				= new TH2F("hV0_DDTofPiP","",300,-15,15,300,-15,15);
TH2F* hV0_DTofPivp				= new TH2F("hV0_DTofPivp","",100,0,10,300,-15,15);
TH2F* hV0_DDDedx				= new TH2F("hV0_DDDedx","",300,0,300,300,0,300);
TH2F* hV0_DDedxvp				= new TH2F("hV0_DDedxvp","",100,0,10,300,0,300);
TH3F* hV0_DTofBvpvr				= new TH3F("hV0_DTofBvpvr","",100,0,10,120,0.6,1.2,20,0,100);
TH2F* hV0_DTofminTpcvr			= new TH2F("hV0_DTofminTpcvr","",400,0,100,500,-5,5);
TH2F* hV0_DTofminTpcvrK0spi		= new TH2F("hV0_DTofminTpcvrK0spi","",400,0,100,500,-5,5);
TH2F* hV0_DTofminTpcvrLpi		= new TH2F("hV0_DTofminTpcvrLpi","",400,0,100,500,-5,5);
TH2F* hV0_DTofminTpcvrLpr		= new TH2F("hV0_DTofminTpcvrLpr","",400,0,100,500,-5,5);
TH2F* hV0_DTpcvpK0spi 			= new TH2F("hV0_DTpcvpK0spi","",450,0,15,400,-10.,10.);
TH2F* hV0_DTpcvpLpr 			= new TH2F("hV0_DTpcvpLpr","",450,0,15,400,-10.,10.);
TH2F* hV0_DTofvpK0spi 			= new TH2F("hV0_DTofvpK0spi","",450,0,15,400,-10.,10.);
TH2F* hV0_DTofvpLpr 			= new TH2F("hV0_DTofvpLpr","",450,0,15,400,-10.,10.);

TH1F* hTrDCA					= new TH1F("hTrDCA","",400,0,100);
TH3F* hTrPtPhiEta				= new TH3F("hTrPtPhiEta","",200,0,10,200,-0.5,6.5,200,-1.,1.);
TH1F* hTrTofminTpc				= new TH1F("hTrTofminTpc","",500,-5,5);

// MC HISTOGRAMS
TH1F* hMCV0Monitor				= new TH1F("hMCV0Monitor","",10,-0.5,9.5);
TH1F* hRCV0Counter				= new TH1F("hRCV0Counter","",100,0,100);
TH1F* hMCV0_PtK0s				= new TH1F("hMCV0_PtK0s","",nPtBins,xBins);
TH1F* hMCV0_PtL					= new TH1F("hMCV0_PtL","",nPtBins,xBins);
TH1F* hMCV0_PtAL				= new TH1F("hMCV0_PtAL","",nPtBins,xBins);
TH1F* hRCV0_PtK0s				= new TH1F("hRCV0_PtK0s","",nPtBins,xBins);
TH1F* hRCV0_PtL					= new TH1F("hRCV0_PtL","",nPtBins,xBins);
TH1F* hRCV0_PtAL				= new TH1F("hRCV0_PtAL","",nPtBins,xBins);
TH2F* hRCV0_APK0s 				= new TH2F("hRCV0_APK0s","",400,-1.,1.,400,0.,0.4);
TH2F* hRCV0_APL 				= new TH2F("hRCV0_APL","",400,-1.,1.,400,0.,0.4);
TH2F* hRCV0_APAL 				= new TH2F("hRCV0_APAL","",400,-1.,1.,400,0.,0.4);
TH3F* hRCV0_PtPhiEtaK0s 		= new TH3F("hRCV0_PtPhiEtaK0s","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hRCV0_PtPhiEtaL 			= new TH3F("hRCV0_PtPhiEtaL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hRCV0_PtPhiEtaAL 			= new TH3F("hRCV0_PtPhiEtaAL","",200,0,20,200,-0.5,6.5,200,-1.,1.);

// MC DAUGHTERS HISTOGRAMS
TH2F* hRCV0_DDedxvpK0s 			= new TH2F("hRCV0_DDedxvpK0s","",450,0,15,300,0,300);
TH2F* hRCV0_DDedxvpL 			= new TH2F("hRCV0_DDedxvpL","",450,0,15,300,0,300);
TH2F* hRCV0_DDedxvpAL 			= new TH2F("hRCV0_DDedxvpAL","",450,0,15,300,0,300);
TH1F* hRCV0_DTpcK0spi			= new TH1F("hRCV0_DTpcK0spi","",400,-10.,10.);
TH1F* hRCV0_DTpcLpi				= new TH1F("hRCV0_DTpcLpi","",400,-10.,10.);
TH1F* hRCV0_DTpcALpi			= new TH1F("hRCV0_DTpcALpi","",400,-10.,10.);
TH1F* hRCV0_DTpcLpr				= new TH1F("hRCV0_DTpcLpr","",400,-10.,10.);
TH1F* hRCV0_DTpcALpr			= new TH1F("hRCV0_DTpcALpr","",400,-10.,10.);
TH1F* hRCV0_DNDedxK0spi			= new TH1F("hRCV0_DNDedxK0spi","",200,0,200);
TH1F* hRCV0_DNDedxLpi			= new TH1F("hRCV0_DNDedxLpi","",200,0,200);
TH1F* hRCV0_DNDedxALpi			= new TH1F("hRCV0_DNDedxALpi","",200,0,200);
TH1F* hRCV0_DNDedxLpr			= new TH1F("hRCV0_DNDedxLpr","",200,0,200);
TH1F* hRCV0_DNDedxALpr			= new TH1F("hRCV0_DNDedxALpr","",200,0,200);
TH2F* hRCV0_DTofBetavpK0s 		= new TH2F("hRCV0_DTofBetavpK0s","",450,0,15,300,0,1.2);
TH2F* hRCV0_DTofBetavpL 		= new TH2F("hRCV0_DTofBetavpL","",450,0,15,300,0,1.2);
TH2F* hRCV0_DTofBetavpAL 		= new TH2F("hRCV0_DTofBetavpAL","",450,0,15,300,0,1.2);
TH1F* hRCV0_DTofK0spi			= new TH1F("hRCV0_DTofK0spi","",400,-10.,10.);
TH1F* hRCV0_DTofLpi				= new TH1F("hRCV0_DTofLpi","",400,-10.,10.);
TH1F* hRCV0_DTofALpi			= new TH1F("hRCV0_DTofALpi","",400,-10.,10.);
TH1F* hRCV0_DTofLpr				= new TH1F("hRCV0_DTofLpr","",400,-10.,10.);
TH1F* hRCV0_DTofALpr			= new TH1F("hRCV0_DTofALpr","",400,-10.,10.);
TH1F* hRCV0_DIPxyK0spi			= new TH1F("hRCV0_DIPxyK0spi","",400,0.,20.);
TH1F* hRCV0_DIPxyLpi			= new TH1F("hRCV0_DIPxyLpi","",400,0.,20.);
TH1F* hRCV0_DIPxyALpi			= new TH1F("hRCV0_DIPxyALpi","",400,0.,20.);
TH1F* hRCV0_DIPxyLpr			= new TH1F("hRCV0_DIPxyLpr","",400,0.,20.);
TH1F* hRCV0_DIPxyALpr			= new TH1F("hRCV0_DIPxyALpr","",400,0.,20.);
TH1F* hRCV0_DIPzK0spi			= new TH1F("hRCV0_DIPzK0spi","",400,0.,20.);
TH1F* hRCV0_DIPzLpi				= new TH1F("hRCV0_DIPzLpi","",400,0.,20.);
TH1F* hRCV0_DIPzALpi			= new TH1F("hRCV0_DIPzALpi","",400,0.,20.);
TH1F* hRCV0_DIPzLpr				= new TH1F("hRCV0_DIPzLpr","",400,0.,20.);
TH1F* hRCV0_DIPzALpr			= new TH1F("hRCV0_DIPzALpr","",400,0.,20.);
TH1F* hRCV0_DNclusK0spi 		= new TH1F("hRCV0_DNclusK0spi","",200,0,200);
TH1F* hRCV0_DNclusLpi 			= new TH1F("hRCV0_DNclusLpi","",200,0,200);
TH1F* hRCV0_DNclusALpi 			= new TH1F("hRCV0_DNclusALpi","",200,0,200);
TH1F* hRCV0_DNclusLpr 			= new TH1F("hRCV0_DNclusLpr","",200,0,200);
TH1F* hRCV0_DNclusALpr 			= new TH1F("hRCV0_DNclusALpr","",200,0,200);

TH2F* hMCvRC_PtK0s				= new TH2F("hMCvRC_PtK0s","",400,0,20,400,0,20);
TH2F* hMCvRC_PtL				= new TH2F("hMCvRC_PtL","",400,0,20,400,0,20);
TH2F* hMCvRC_PtAL				= new TH2F("hMCvRC_PtAL","",400,0,20,400,0,20);
TH2F* hMCvRC_PhiK0s				= new TH2F("hMCvRC_PhiK0s","",400,-0.5,6.5,400,-0.5,6.5);
TH2F* hMCvRC_PhiL				= new TH2F("hMCvRC_PhiL","",400,-0.5,6.5,400,-0.5,6.5);
TH2F* hMCvRC_PhiAL				= new TH2F("hMCvRC_PhiAL","",400,-0.5,6.5,400,-0.5,6.5);
TH2F* hMCvRC_EtaK0s				= new TH2F("hMCvRC_EtaK0s","",400,-1.,1.,400,-1.,1.);
TH2F* hMCvRC_EtaL				= new TH2F("hMCvRC_EtaL","",400,-1.,1.,400,-1.,1.);
TH2F* hMCvRC_EtaAL				= new TH2F("hMCvRC_EtaAL","",400,-1.,1.,400,-1.,1.);
TH2F* hMCvRC_DCAddK0s			= new TH2F("hMCvRC_DCAddK0s","",400,0,2,400,0,2);
TH2F* hMCvRC_DCAddL				= new TH2F("hMCvRC_DCAddL","",400,0,2,400,0,2);
TH2F* hMCvRC_DCAddAL			= new TH2F("hMCvRC_DCAddAL","",400,0,2,400,0,2);
TH2F* hMCvRC_CPAK0s				= new TH2F("hMCvRC_CPAK0s","",400,0.98,1.02,400,0.98,1.02);
TH2F* hMCvRC_CPAL				= new TH2F("hMCvRC_CPAL","",400,0.98,1.02,400,0.98,1.02);
TH2F* hMCvRC_CPAAL				= new TH2F("hMCvRC_CPAAL","",400,0.98,1.02,400,0.98,1.02);
TH2F* hMCvRC_RK0s				= new TH2F("hMCvRC_RK0s","",400,0,150,400,0,150);
TH2F* hMCvRC_RL					= new TH2F("hMCvRC_RL","",400,0,150,400,0,150);
TH2F* hMCvRC_RAL				= new TH2F("hMCvRC_RAL","",400,0,150,400,0,150);
TH2F* hMCvRC_APK0s				= new TH2F("hMCvRC_APK0s","",400,0,5,400,0,5);
TH2F* hMCvRC_APL				= new TH2F("hMCvRC_APL","",400,0,5,400,0,5);
TH2F* hMCvRC_APAL				= new TH2F("hMCvRC_APAL","",400,0,5,400,0,5);
TH2F* hMCvRC_DTofTimeK0spi 		= new TH2F("hMCvRC_DTofTimeK0spi","",400,1e4,35000,400,1e4,35000);
TH2F* hMCvRC_DTofLengthK0spi 	= new TH2F("hMCvRC_DTofLengthK0spi","",400,0,1000,400,0,1000);
TH2F* hMCvRC_DTofTimeLpi 		= new TH2F("hMCvRC_DTofTimeLpi","",400,1e4,35000,400,1e4,35000);
TH2F* hMCvRC_DTofLengthLpi 		= new TH2F("hMCvRC_DTofLengthLpi","",400,0,1000,400,0,1000);
TH2F* hMCvRC_DTofTimeALpi 		= new TH2F("hMCvRC_DTofTimeALpi","",400,1e4,35000,400,1e4,35000);
TH2F* hMCvRC_DTofLengthALpi 	= new TH2F("hMCvRC_DTofLengthALpi","",400,0,1000,400,0,1000);
TH2F* hMCvRC_DTofTimeLpr 		= new TH2F("hMCvRC_DTofTimeLpr","",400,1e4,35000,400,1e4,35000);
TH2F* hMCvRC_DTofLengthLpr 		= new TH2F("hMCvRC_DTofLengthLpr","",400,0,1000,400,0,1000);
TH2F* hMCvRC_DTofTimeALpr 		= new TH2F("hMCvRC_DTofTimeALpr","",400,1e4,35000,400,1e4,35000);
TH2F* hMCvRC_DTofLengthALpr 	= new TH2F("hMCvRC_DTofLengthALpr","",400,0,1000,400,0,1000);


bool MakeChain(const Char_t *inputFile="test.list") {

	if (!mChain) mChain = new TChain("PIDTree");
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
	if (v0->GetPt() > 14.0)						return false;
	//if (v0->GetDCAPV() < 0.05)				return false; // dca of extrapolated v0 to pv? or daughters to pv?
	if (v0->GetDCAV0Daughters() > 0.5)			return false;
	if (v0->GetV0CosinePA() < 0.997)			return false; // :(
	if (v0->GetRadius() < 5.)					return false;
	if (v0->GetRadius() > 100.)					return false;

	if (!SelectV0daughter(v0->GetPosAnalysisTrack()))	return false;
	if (!SelectV0daughter(v0->GetNegAnalysisTrack()))	return false;

	// IM photonic e+e- rejection cut?

	return true;
}

bool SelectV0_MC(AliAnalysisPIDParticle* p, Int_t part) {

	if (p->GetPdgCode() != PDG_IDS[part]) 				return false;
	if (fabs(p->GetEta()) > 0.8) 					return false;

	return true;
}

bool SelectV0daughter(AliAnalysisPIDTrack* t) {

	if (fabs(t->GetEta()) > 0.8) return false;
	if (fabs(t->GetImpactParameter(0)) < 0.05) return false;

	return true;
}

bool SelectTrack(AliAnalysisPIDTrack* t) {

	if (fabs(t->GetEta()) > 0.8) return false;

	return true;
}

bool IsK0s(AliAnalysisPIDV0* v0, Int_t cutFlag) {

	AliAnalysisPIDTrack* trP = v0->GetPosAnalysisTrack();
	AliAnalysisPIDTrack* trN = v0->GetNegAnalysisTrack();

	if ((cutFlag<1) && fabs(trP->GetNSigmaPionTOF())>3.) return false;
	if ((cutFlag<1) && fabs(trN->GetNSigmaPionTOF())>3.) return false;
	if ((cutFlag<2) && fabs(trP->GetNSigmaPionTPC())>3.) return false;
	if ((cutFlag<2) && fabs(trN->GetNSigmaPionTPC())>3.) return false;

	return true;
}


bool IsL(AliAnalysisPIDV0* v0, Int_t cutFlag) {

	AliAnalysisPIDTrack* trP = v0->GetPosAnalysisTrack();
	AliAnalysisPIDTrack* trN = v0->GetNegAnalysisTrack();

	if ((cutFlag<1) && fabs(trP->GetNSigmaProtonTOF())>3.) 	return false;
	if ((cutFlag<1) && fabs(trN->GetNSigmaPionTOF())>3.) 	return false;
	if ((cutFlag<2) && fabs(trP->GetNSigmaProtonTPC())>3.) 	return false;
	if ((cutFlag<2) && fabs(trN->GetNSigmaPionTPC())>3.) 	return false;

	return true;
}

bool IsAL(AliAnalysisPIDV0* v0, Int_t cutFlag) {

	AliAnalysisPIDTrack* trP = v0->GetPosAnalysisTrack();
	AliAnalysisPIDTrack* trN = v0->GetNegAnalysisTrack();

	if ((cutFlag<1) && fabs(trP->GetNSigmaPionTOF())>3.) 	return false;
	if ((cutFlag<1) && fabs(trN->GetNSigmaProtonTOF())>3.) 	return false;
	if ((cutFlag<2) && fabs(trP->GetNSigmaPionTPC())>3.) 	return false;
	if ((cutFlag<2) && fabs(trN->GetNSigmaProtonTPC())>3.) 	return false;

	return true;
}

void myLegendSetUp(TLegend *currentLegend=0, float currentTextSize=0.07, int columns=2)	{
	currentLegend->SetTextFont(42);
	currentLegend->SetBorderSize(0);
	currentLegend->SetFillStyle(0);
	currentLegend->SetFillColor(0);
	currentLegend->SetMargin(0.25);
	currentLegend->SetTextSize(currentTextSize);
	currentLegend->SetEntrySeparation(0.5);
	currentLegend->SetNColumns(columns);
	
	return;
}

Float_t* ExtractYield(TH1D* hist = 0, TTree* tree = 0, Int_t method = 0, Int_t part = 0) {	// extracting with RooFit, 0 is sideband, 1 is sideband bg, 2 is fit
																		// 3 is unbinned fit
																		// part: 0 for k0s, 1 for lambdas
	static Float_t val[2];
	val[0] = 0; val[1] = 0;
	Float_t fitMin = -0.03, fitMax = 0.03;

	switch (method) {
		case 0 :
			Float_t subRange = (!part) ? 0.03 : 0.015; 
			hist->GetXaxis()->SetRangeUser(-subRange,subRange);
			Float_t mean = hist->GetMean(); 
			Float_t rms = hist->GetRMS();
			printf("STATUS: part is %i rms is %f \n", part, rms);
			hist->GetXaxis()->SetRange();
			val[0] = hist->Integral(hist->FindBin(mean-3.*rms),hist->FindBin(mean+3.*rms));
			printf("STATUS: int inside is %f \n", val[0]);
			val[1] = sqrt(val[0]);
			break;

		case 1 :
			Float_t subRange = (!part) ? 0.03 : 0.015; 
			hist->GetXaxis()->SetRangeUser(-subRange,subRange); 
			Float_t mean = hist->GetMean(); 
			Float_t rms = hist->GetRMS();
			hist->GetXaxis()->SetRange();
			val[0] = hist->Integral(hist->FindBin(mean-6.*rms),hist->FindBin(mean-3.*rms));
			val[0] += hist->Integral(hist->FindBin(mean+3.*rms),hist->FindBin(mean+6.*rms));
			printf("STATUS: int outside is %f \n", val[0]);
			val[1] = sqrt(val[0]);
			break;

		case 2 :	// fallthrough
		case 3 :
			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			if (method==2) {	
				hist->Rebin(8);
				RooDataHist DT_set("DT_set","DT_hist",MassDT,Import(*hist)); }
			else RooDataSet DT_set("DT_set","DT_tree",MassDT,Import(*tree));
		
			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0,0.01);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e06);
		
			RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0,0.01);
			//RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus2A,pGaus2B); 
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus1A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e06);
		
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-200,200);
			RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);//RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e06);
		
			//RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
			RooFitResult* fR = fTotal.fitTo(DT_set,Save(),PrintLevel(-1));
		
			RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));
			//printf("Errors are %f and %f, total is %f or %f wrt to %f \n", nGaus1.getError(), nGaus2.getError(), nGaus1.getError()+nGaus2.getError(),sqrt(nGaus1.getError()*nGaus1.getError()+nGaus2.getError()*nGaus2.getError()),nGaus.getPropagatedError(*fR));
		
			cFits[canCounter%6]->cd(1+canCounter/6);
			RooPlot* plot1 = MassDT.frame(Title(" "));
			DT_set.plotOn(plot1,MarkerSize(0.4));
			fTotal.plotOn(plot1,LineWidth(1),LineColor(kRed));
			plot1->SetMinimum(1e-05);
			plot1->SetMaximum(1.35*plot1->GetMaximum());
			plot1->GetXaxis()->SetTitleSize(0.05);
			plot1->GetYaxis()->SetTitleSize(0.05);
			plot1->Draw();
			TLegend *leg1 = new TLegend(0.075,0.7,0.5,0.88);
			myLegendSetUp(leg1,0.065,1);
			leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[canCounter/6],xBins[1+canCounter/6])," ");
			leg1->AddEntry((TObject*)0,cNames[canCounter%6]+Form(" , #chi^{2}/ndf = %4.2f",plot1->chiSquare())," ");
			leg1->AddEntry((TObject*)0,Form("%4.1f #pm %4.1f",nGaus.getVal(),nGaus.getPropagatedError(*fR))," ");
			leg1->Draw();
	
			val[0] = nGaus.getVal();
			printf("STATUS: int from fit is %f \n", val[0]);
			val[1] = nGaus.getPropagatedError(*fR);
			canCounter++;
			break;

	}

	return val;
}

Float_t CalculateCPA(AliAnalysisPIDParticle* p, Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
  // calculates the pointing angle of the V0 wrt a reference point

  Double_t momV0[3]; //momentum of the V0
  TLorentzVector fLV;
  fLV.SetPtEtaPhiM(p->GetPt(),p->GetEta(),p->GetPhi(),p->GetMass());
  momV0[0] = fLV.X(); momV0[1] = fLV.Y(); momV0[2] = fLV.Z();

  Double_t deltaPos[3]; //vector between the reference point and the V0 vertex
  deltaPos[0] = fPos[0] - refPointX;
  deltaPos[1] = fPos[1] - refPointY;
  deltaPos[2] = fPos[2] - refPointZ;

  Double_t momV02    = momV0[0]*momV0[0] + momV0[1]*momV0[1] + momV0[2]*momV0[2];
  Double_t deltaPos2 = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];

  Double_t cosinePointingAngle = (deltaPos[0]*momV0[0] +
				  deltaPos[1]*momV0[1] +
				  deltaPos[2]*momV0[2] ) /
    TMath::Sqrt(momV02 * deltaPos2);
  
  return cosinePointingAngle;
}

Double_t* CalculateAP(AliAnalysisPIDV0* v0) { //[0]...alpha, [1]...pt

	static Double_t val[2]; 

	TVector3 momNeg, momPos, momTot;
	momNeg.SetPtEtaPhi(v0->GetNegAnalysisTrack()->GetPt(),v0->GetNegAnalysisTrack()->GetEta(),v0->GetNegAnalysisTrack()->GetPhi());
	momPos.SetPtEtaPhi(v0->GetPosAnalysisTrack()->GetPt(),v0->GetPosAnalysisTrack()->GetEta(),v0->GetPosAnalysisTrack()->GetPhi());
	//momTot.SetPtEtaPhi(v0->GetPt(),v0->GetEta(),v0->GetPhi());
	momTot = momNeg + momPos;

	Double_t lQlNeg = momNeg.Dot(momTot)/momTot.Mag();
	Double_t lQlPos = momPos.Dot(momTot)/momTot.Mag();

 	val[0] = (lQlPos - lQlNeg)/(lQlPos + lQlNeg);
	val[1] = momPos.Perp(momTot);

	return val;
}

Double_t CalculatePhi(AliAnalysisPIDV0* v0) {

	TVector3 momNeg, momPos, momTot;
	momNeg.SetPtEtaPhi(v0->GetNegAnalysisTrack()->GetPt(),v0->GetNegAnalysisTrack()->GetEta(),v0->GetNegAnalysisTrack()->GetPhi());
	momPos.SetPtEtaPhi(v0->GetPosAnalysisTrack()->GetPt(),v0->GetPosAnalysisTrack()->GetEta(),v0->GetPosAnalysisTrack()->GetPhi());
	momTot = momNeg + momPos;

	Double_t phi = (momTot.Phi() > 0) ? momTot.Phi() : (momTot.Phi() + 2.*TMath::Pi() );

	return phi;
}

//AliAnalysisPIDParticle

//cutflag: 0 is tpc+tof (no bg), 1 is tpc (+bg), 2 is bg
void readTree_V0(Int_t nEvents=10, Int_t cutFlag=0, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();
	TList* lHist = gDirectory->GetList();
	int iHist = 0;

	if (!MakeChain(inputFile)) printf("Couldn't create the chain! \n", );
	else printf("Chain created with %i entries \n", mChain->GetEntries());

	mChain->SetBranchAddress("AnalysisTrack",&bTracks);
	mChain->SetBranchAddress("AnalysisV0Track",&bV0s);
	mChain->SetBranchAddress("AnalysisEvent",&mEvent);
	if (flagMC) mChain->SetBranchAddress("AnalysisParticle",&bParticles);	

	nEvents = (nEvents < mChain->GetEntries()) ? nEvents : mChain->GetEntries();
	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		//if (iEv!=139254) continue;

		hEventMonitor->Fill(0);
		if (iEv%10000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		bTracks->Clear();
		mChain->GetEntry(iEv);
		if (!mEvent) continue;
		hEventMonitor->Fill(1);

		Float_t vz = mEvent->GetVertexZ();
		Int_t refmult = mEvent->GetReferenceMultiplicity();
		Float_t v0mult = mEvent->GetV0Mmultiplicity() ;
		//printf("EVENT #%i , VZ is %9.9f , rmult is %i , v0mult is %f \n", iEv, vz, refmult, v0mult);

		if (!iEv) mEvent->PrintEventSelection();
		if (!SelectEvent(mEvent)) continue;
		hEventMonitor->Fill(2);

		if (flagMC) //V0 MC analysis
		{
			Int_t nPs = bParticles->GetEntriesFast();
			for (int iP = 0; iP < nPs; ++iP)	{
				
				hMCV0Monitor->Fill(0);
				AliAnalysisPIDParticle* p = (AliAnalysisPIDParticle*)bParticles->At(iP);
				if (!p) continue;
				hMCV0Monitor->Fill(1);

				Int_t v0id;
				if (SelectV0_MC(p,0)) {				// select mc v0
					hMCV0Monitor->Fill(2);
					hMCV0_PtK0s->Fill(p->GetPt()); 
					v0id = 0;}
				else if (SelectV0_MC(p,1)) {
					hMCV0Monitor->Fill(2);
					hMCV0_PtL->Fill(p->GetPt()); 
					v0id = 1;}
				else if (SelectV0_MC(p,2)) {
					hMCV0Monitor->Fill(2);
					hMCV0_PtAL->Fill(p->GetPt()); 
					v0id = 2;}
				else continue;

				// find mc daughters
				//AliAnalysisPIDParticle* dP = FindPosDaughterMC(p,v0id);
				//AliAnalysisPIDParticle* dN = FindNegDaughterMC(p,v0id);

																// use vectors perhaps for matching?
				Int_t mcLabel = p->GetLabel();

				Int_t rcV0Count = 0;
				Int_t nV0s = bV0s->GetEntriesFast();
				//for (int iV0 = 0; (iV0 < nV0s && rcV0Count < 1); ++iV0)	{
				TString echo1;
				for (int iV0 = 0; iV0 < nV0s; ++iV0)	{

					AliAnalysisPIDV0* v0rc = (AliAnalysisPIDV0*)bV0s->At(iV0);
					if (!v0rc) continue;

					if (v0rc->GetMCPdgCode() != PDG_IDS[v0id]) continue;

					AliAnalysisPIDTrack* trPrc = v0rc->GetPosAnalysisTrack();
					AliAnalysisPIDTrack* trNrc = v0rc->GetNegAnalysisTrack();
					if (trPrc->GetMCMotherLabel() != mcLabel) continue;
					//if (trPrc->GetMCMotherPrimary() != 1) continue;
					
					rcV0Count++;
					if (rcV0Count>1) printf(echo1.Data());
					echo1 = Form("nEV %i nP %i v0 finds %i , v0id %i , mcpt %f , rcpt %f , mK0 %f , mL %f , tr+pt %f , tr-pt %f : \n", 
						iEv, iP, rcV0Count, v0id, p->GetPt(), v0rc->GetPt(), v0rc->GetIMK0s(), v0rc->GetIML(), trPrc->GetPt(), trNrc->GetPt());
					echo1 += Form("-------------------------- r %f , dcad %f , cos %f , eta %f , dcapv %f \n", 
						v0rc->GetRadius(), v0rc->GetDCAV0Daughters(), v0rc->GetV0CosinePA(), v0rc->GetEta(), v0rc->GetDCAPV());
					if (rcV0Count>1) printf(echo1.Data());

					Double_t* aprc = CalculateAP(v0rc);
					Double_t phirc = CalculatePhi(v0rc);

					if (v0id == 0) {
						hMCvRC_PtK0s->Fill(p->GetPt(),v0rc->GetPt());
						hMCvRC_PhiK0s->Fill(p->GetPhi(),phirc);
						hMCvRC_EtaK0s->Fill(p->GetEta(),v0rc->GetEta());
						hMCvRC_DCAddK0s->Fill(0,v0rc->GetDCAV0Daughters());		//trees don't store mc helices	//mc true should be 0
						hMCvRC_CPAK0s->Fill(1,v0rc->GetV0CosinePA());			//trees don't store PV_x,y 		//mc true should be 1
						hMCvRC_RK0s->Fill(1,v0rc->GetRadius());					//trees don't store mc helices
						if (fabs(*(aprc+0))>0.005) hMCvRC_APK0s->Fill(1,*(aprc+1)/fabs(*(aprc+0)));	//trees don't store v0 momentum vector
						hRCV0_APK0s->Fill(*(aprc+0),*(aprc+1));
						hRCV0_PtPhiEtaK0s->Fill(v0rc->GetPt(), phirc, v0rc->GetEta());
						if (trPrc->HasTPCPID()) {
							hRCV0_DDedxvpK0s->Fill(trPrc->GetP(),trPrc->GetTPCdEdx());
							hRCV0_DTpcK0spi->Fill(trPrc->GetNSigmaPionTPC());
							hRCV0_DNDedxK0spi->Fill(trPrc->GetTPCdEdxN()); 								}
						if (trNrc->HasTPCPID()) {
							hRCV0_DDedxvpK0s->Fill(trNrc->GetP(),trNrc->GetTPCdEdx());
							hRCV0_DTpcK0spi->Fill(trNrc->GetNSigmaPionTPC());
							hRCV0_DNDedxK0spi->Fill(trNrc->GetTPCdEdxN()); 								}
						if (trPrc->HasTOFPID()) {
							hRCV0_DTofBetavpK0s->Fill(trPrc->GetP(),trPrc->GetTOFExpBeta(AliPID::kPion));
							//hRCV0_DTofBetavpK0s->Fill(trPrc->GetP(),trPrc->GetTOFBeta());
							hRCV0_DTofK0spi->Fill(trPrc->GetNSigmaPionTOF());
							hMCvRC_DTofTimeK0spi->Fill(trPrc->GetMCTOFTime(),trPrc->GetTOFTime());
							hMCvRC_DTofLengthK0spi->Fill(trPrc->GetMCTOFLength(),trPrc->GetTOFLength());	}
						if (trNrc->HasTOFPID()) {
							hRCV0_DTofBetavpK0s->Fill(trNrc->GetP(),trNrc->GetTOFExpBeta(AliPID::kPion));
							hRCV0_DTofK0spi->Fill(trNrc->GetNSigmaPionTOF());
							hMCvRC_DTofTimeK0spi->Fill(trNrc->GetMCTOFTime(),trNrc->GetTOFTime());
							hMCvRC_DTofLengthK0spi->Fill(trNrc->GetMCTOFLength(),trNrc->GetTOFLength());	}
						hRCV0_DIPxyK0spi->Fill(trPrc->GetImpactParameter(0));
						hRCV0_DIPxyK0spi->Fill(trNrc->GetImpactParameter(0));
						hRCV0_DIPzK0spi->Fill(trPrc->GetImpactParameter(1));
						hRCV0_DIPzK0spi->Fill(trNrc->GetImpactParameter(1));
						hRCV0_DNclusK0spi->Fill(trPrc->GetTPCNcls());
						hRCV0_DNclusK0spi->Fill(trNrc->GetTPCNcls());
					}
					if (v0id == 1) {
						hMCvRC_PtL->Fill(p->GetPt(),v0rc->GetPt());
						hMCvRC_PhiL->Fill(p->GetPhi(),phirc);
						hMCvRC_EtaL->Fill(p->GetEta(),v0rc->GetEta());
						hMCvRC_DCAddL->Fill(0,v0rc->GetDCAV0Daughters());		//trees don't store mc helices
						hMCvRC_CPAL->Fill(1,v0rc->GetV0CosinePA());			//trees don't store PV_x,y
						hMCvRC_RL->Fill(1,v0rc->GetRadius());					//trees don't store mc helices
						if (fabs(*(aprc+0))>0.005) hMCvRC_APL->Fill(1,*(aprc+1)/fabs(*(aprc+0)));	//trees don't store v0 momentum vector
						hRCV0_APL->Fill(*(aprc+0),*(aprc+1));
						hRCV0_PtPhiEtaL->Fill(v0rc->GetPt(), phirc, v0rc->GetEta());
						if (trPrc->HasTPCPID()) {
							hRCV0_DDedxvpL->Fill(trPrc->GetP(),trPrc->GetTPCdEdx());
							hRCV0_DTpcLpr->Fill(trPrc->GetNSigmaProtonTPC());
							hRCV0_DNDedxLpr->Fill(trPrc->GetTPCdEdxN()); 								}
						if (trNrc->HasTPCPID()) {
							hRCV0_DDedxvpL->Fill(trNrc->GetP(),trNrc->GetTPCdEdx());
							hRCV0_DTpcLpi->Fill(trNrc->GetNSigmaPionTPC());
							hRCV0_DNDedxLpi->Fill(trNrc->GetTPCdEdxN()); 								}
						if (trPrc->HasTOFPID()) {
							hRCV0_DTofBetavpL->Fill(trPrc->GetP(),trPrc->GetTOFExpBeta(AliPID::kProton));
							hRCV0_DTofLpr->Fill(trPrc->GetNSigmaProtonTOF());
							hMCvRC_DTofTimeLpr->Fill(trPrc->GetMCTOFTime(),trPrc->GetTOFTime());
							hMCvRC_DTofLengthLpr->Fill(trPrc->GetMCTOFLength(),trPrc->GetTOFLength());	}
						if (trNrc->HasTOFPID()) {
							hRCV0_DTofBetavpL->Fill(trNrc->GetP(),trNrc->GetTOFExpBeta(AliPID::kPion));
							hRCV0_DTofLpi->Fill(trNrc->GetNSigmaPionTOF());
							hMCvRC_DTofTimeLpi->Fill(trNrc->GetMCTOFTime(),trNrc->GetTOFTime());
							hMCvRC_DTofLengthLpi->Fill(trNrc->GetMCTOFLength(),trNrc->GetTOFLength());	}
						hRCV0_DIPxyLpr->Fill(trPrc->GetImpactParameter(0));
						hRCV0_DIPxyLpi->Fill(trNrc->GetImpactParameter(0));
						hRCV0_DIPzLpr->Fill(trPrc->GetImpactParameter(1));
						hRCV0_DIPzLpi->Fill(trNrc->GetImpactParameter(1));
						hRCV0_DNclusLpr->Fill(trPrc->GetTPCNcls());
						hRCV0_DNclusLpi->Fill(trNrc->GetTPCNcls());
					}
					if (v0id == 2) {
						hMCvRC_PtAL->Fill(p->GetPt(),v0rc->GetPt());
						hMCvRC_PhiAL->Fill(p->GetPhi(),phirc);
						hMCvRC_EtaAL->Fill(p->GetEta(),v0rc->GetEta());
						hMCvRC_DCAddAL->Fill(0,v0rc->GetDCAV0Daughters());		//trees don't store mc helices
						hMCvRC_CPAAL->Fill(1,v0rc->GetV0CosinePA());			//trees don't store PV_x,y
						hMCvRC_RAL->Fill(1,v0rc->GetRadius());					//trees don't store mc helices
						if (fabs(*(aprc+0))>0.005) hMCvRC_APAL->Fill(1,*(aprc+1)/fabs(*(aprc+0)));	//trees don't store v0 momentum vector
						hRCV0_APAL->Fill(*(aprc+0),*(aprc+1));
						hRCV0_PtPhiEtaAL->Fill(v0rc->GetPt(), phirc, v0rc->GetEta());
						if (trPrc->HasTPCPID()) {
							hRCV0_DDedxvpAL->Fill(trPrc->GetP(),trPrc->GetTPCdEdx());
							hRCV0_DTpcALpi->Fill(trPrc->GetNSigmaPionTPC());
							hRCV0_DNDedxALpi->Fill(trPrc->GetTPCdEdxN()); 								}
						if (trNrc->HasTPCPID()) {
							hRCV0_DDedxvpAL->Fill(trNrc->GetP(),trNrc->GetTPCdEdx());
							hRCV0_DTpcALpr->Fill(trNrc->GetNSigmaProtonTPC());
							hRCV0_DNDedxALpr->Fill(trNrc->GetTPCdEdxN()); 								}
						if (trPrc->HasTOFPID()) {
							hRCV0_DTofBetavpAL->Fill(trPrc->GetP(),trPrc->GetTOFExpBeta(AliPID::kPion));
							hRCV0_DTofALpi->Fill(trPrc->GetNSigmaPionTOF());
							hMCvRC_DTofTimeALpi->Fill(trPrc->GetMCTOFTime(),trPrc->GetTOFTime());
							hMCvRC_DTofLengthALpi->Fill(trPrc->GetMCTOFLength(),trPrc->GetTOFLength());	}
						if (trNrc->HasTOFPID()) {
							hRCV0_DTofBetavpAL->Fill(trNrc->GetP(),trNrc->GetTOFExpBeta(AliPID::kProton));
							hRCV0_DTofALpr->Fill(trNrc->GetNSigmaProtonTOF());
							hMCvRC_DTofTimeALpr->Fill(trNrc->GetMCTOFTime(),trNrc->GetTOFTime());
							hMCvRC_DTofLengthALpr->Fill(trNrc->GetMCTOFLength(),trNrc->GetTOFLength());	}
						hRCV0_DIPxyALpi->Fill(trPrc->GetImpactParameter(0));
						hRCV0_DIPxyALpr->Fill(trNrc->GetImpactParameter(0));
						hRCV0_DIPzALpi->Fill(trPrc->GetImpactParameter(1));
						hRCV0_DIPzALpr->Fill(trNrc->GetImpactParameter(1));
						hRCV0_DNclusALpi->Fill(trPrc->GetTPCNcls());
						hRCV0_DNclusALpr->Fill(trNrc->GetTPCNcls());
					}

					if (!SelectV0(v0rc)) continue;
					if (v0id == 0 && IsK0s(v0rc,cutFlag)) 	{
						hRCV0_PtK0s->Fill(v0rc->GetPt());			}
					if (v0id == 1 && IsL(v0rc,cutFlag)) 		{
						hRCV0_PtL->Fill(v0rc->GetPt());				}
					if (v0id == 2 && IsAL(v0rc,cutFlag)) 		{
						hRCV0_PtAL->Fill(v0rc->GetPt());			}

					//break;
				}
				hRCV0Counter->Fill(rcV0Count);
				if (rcV0Count > 1) {
					printf("nEV %i nP %i counter %i Too many v0s reconstructed!!!\n", iEv, iP, rcV0Count);
					//return 0;
				}
			}
		}

		Int_t V0Count = 0; 			// VZERO ANALYSIS HERE
		Int_t nV0s = bV0s->GetEntriesFast();
		for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
			
			hV0Monitor->Fill(0);
			AliAnalysisPIDV0* v0 = (AliAnalysisPIDV0*)bV0s->At(iV0);
			if (!v0) continue;
			hV0Monitor->Fill(1);

			if (!SelectV0(v0)) continue;
			hV0Monitor->Fill(2);
			V0Count++;

			AliAnalysisPIDTrack* trP = v0->GetPosAnalysisTrack();
			AliAnalysisPIDTrack* trN = v0->GetNegAnalysisTrack();
			Double_t* ap = CalculateAP(v0rc);
			Double_t phi = CalculatePhi(v0rc);

			if (trP->HasTPCPID()) hV0_DHasTPC->Fill(trP->GetPt());
			if (trN->HasTPCPID()) hV0_DHasTPC->Fill(trN->GetPt());
			if (trP->HasTOFPID()) hV0_DHasTOF->Fill(trP->GetPt());
			if (trN->HasTOFPID()) hV0_DHasTOF->Fill(trN->GetPt());

			hV0_DPt->Fill(trP->GetPt());
			hV0_DPt->Fill(trN->GetPt());
			hV0_DDTofPiPi->Fill(trP->GetNSigmaPionTOF(),trN->GetNSigmaPionTOF());
			hV0_DDTofPiP->Fill(trP->GetNSigmaPionTOF(),trN->GetNSigmaProtonTOF());
			hV0_DTofPivp->Fill(trP->GetP(),trP->GetNSigmaPionTOF());
			hV0_DTofPivp->Fill(trN->GetP(),trN->GetNSigmaPionTOF());
			hV0_DDDedx->Fill(trP->GetTPCdEdx(),trN->GetTPCdEdx());
			hV0_DDedxvp->Fill(trP->GetP(),trP->GetTPCdEdx());
			hV0_DDedxvp->Fill(trN->GetP(),trN->GetTPCdEdx());
			hV0_DTofBvpvr->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kPion),v0->GetRadius());
			hV0_DTofBvpvr->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kPion),v0->GetRadius());
			hV0_DTofBvpvr->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kProton),v0->GetRadius());
			hV0_DTofBvpvr->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kProton),v0->GetRadius());
			hV0_Radiusvpt->Fill(v0->GetPt(),v0->GetRadius());
			hV0_PtPhiEta->Fill(v0->GetPt(),phi,v0->GetEta());
			hV0_AP->Fill(*(ap+0),*(ap+1));

			bool noCuts = 0; float masscut = 0.0075;
			if (noCuts || IsK0s(v0,cutFlag)) 	{
				hV0_IMK0s->Fill(v0->GetIMK0s());	
				hV0_IMPtK0s->Fill(v0->GetIMK0s(),v0->GetPt());
				mMasses->Fill(v0->GetIMK0s(),v0->GetPt(),0);
				if (fabs(v0->GetIMK0s())<2.*masscut) {
					hV0_PtK0s->Fill(v0->GetPt());
					Float_t delta = trP->GetNSigmaPionTOF() - trP->GetNSigmaPionTPC();	//TOF STUDY, cutflag should be 1 or 2
					hV0_DTofminTpcvrK0spi->Fill(v0->GetRadius(),delta);
					delta = trN->GetNSigmaPionTOF() - trN->GetNSigmaPionTPC();
					hV0_DTofminTpcvrK0spi->Fill(v0->GetRadius(),delta);
					hV0_DTpcvpK0spi->Fill(trP->GetP(),trP->GetNSigmaPionTPC());
					hV0_DTpcvpK0spi->Fill(trN->GetP(),trN->GetNSigmaPionTPC());
					hV0_DTofvpK0spi->Fill(trP->GetP(),trP->GetNSigmaPionTOF());
					hV0_DTofvpK0spi->Fill(trN->GetP(),trN->GetNSigmaPionTOF()); 		}
				}
			if (noCuts || IsL(v0,cutFlag)) 		{
				hV0_IML->Fill(v0->GetIML());
				hV0_IMPtL->Fill(v0->GetIML(),v0->GetPt());
				mMasses->Fill(v0->GetIML(),v0->GetPt(),1);		
				if (fabs(v0->GetIML())<masscut) {
					hV0_PtL->Fill(v0->GetPt());
					Float_t delta = trP->GetNSigmaProtonTOF() - trP->GetNSigmaProtonTPC();
					hV0_DTofminTpcvrLpr->Fill(v0->GetRadius(),delta);
					delta = trN->GetNSigmaPionTOF() - trN->GetNSigmaPionTPC();
					hV0_DTofminTpcvrLpi->Fill(v0->GetRadius(),delta);
					hV0_DTpcvpLpr->Fill(trP->GetP(),trP->GetNSigmaProtonTPC());
					hV0_DTofvpLpr->Fill(trP->GetP(),trP->GetNSigmaProtonTOF()); 			}
				}
			if (noCuts || IsAL(v0,cutFlag)) 	{
				hV0_IMAL->Fill(v0->GetIMAL());
				hV0_IMPtAL->Fill(v0->GetIMAL(),v0->GetPt());
				mMasses->Fill(v0->GetIMAL(),v0->GetPt(),2);
				if (fabs(v0->GetIMAL())<masscut) {
					hV0_PtAL->Fill(v0->GetPt());
					Float_t delta = trP->GetNSigmaPionTOF() - trP->GetNSigmaPionTPC();
					hV0_DTofminTpcvrLpi->Fill(v0->GetRadius(),delta);
					delta = trN->GetNSigmaProtonTOF() - trN->GetNSigmaProtonTPC();
					hV0_DTofminTpcvrLpr->Fill(v0->GetRadius(),delta);		}
				}


			// MC STUDY OF FAKE V0s
			

		}

		Int_t trCount = 0;		// TRACKS ANALYSIS HERE
		Int_t nTracks = bTracks->GetEntriesFast();
		for (int iTr = 0; iTr < nTracks; ++iTr)	{
			
			hTrackMonitor->Fill(0);
			AliAnalysisPIDTrack* track = (AliAnalysisPIDTrack*)bTracks->At(iTr);
			if (!track) continue;
			hTrackMonitor->Fill(1);

			hTrDCA->Fill(track->GetImpactParameter(0));
			hTrPtPhiEta->Fill(track->GetPt(),track->GetPhi(),track->GetEta());
			
			if (!SelectTrack(track)) continue;
			hTrackMonitor->Fill(2);
			trCount++;

			Float_t delta = track->GetNSigmaPionTOF() - track->GetNSigmaPionTPC();
			hTrTofminTpc->Fill(delta);
		}
		
		hV0TrCounter->Fill(V0Count,trCount);
	}	// EVENT LOOP FINISHED

	for (int iC = 0; iC < 6; ++iC)	{
		cFits[iC] = new TCanvas(Form("cFits%i",iC),"",2800,2000);
		cFits[iC]->Divide(7,5,0.0005,0.0005);	}

	gROOT->cd();	//TTree operations can't be done in read only files
	
	for (int iBin = 1; iBin < nPtBins+1; ++iBin)		// 0 is underflow
	{
		//if (iBin != 4) continue;
		Float_t* y;
		y = ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin),0,0,0);
		hSBinK0s->SetBinContent(iBin,*(y+0));	
		hSBinK0s->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin),0,1,0);
		hSBoutK0s->SetBinContent(iBin,*(y+0));	
		hSBoutK0s->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin),0,2,0);	// 0 is underflow bin
		hYieldK0s->SetBinContent(iBin,*(y+0));	
		hYieldK0s->SetBinError(iBin,*(y+1));
		y = ExtractYield(0,mMasses->CopyTree(Form("lPt>%f&&lPt<%f&&lPart==0",xBins[iBin-1],xBins[iBin])),3,0);
		hYieldUBK0s->SetBinContent(iBin,*(y+0));	
		hYieldUBK0s->SetBinError(iBin,*(y+1));
		//y = ExtractYield(0,mMasses->CopyTree("lPart==0"),3,0);

		y = ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin),0,0,1);
		hSBinL->SetBinContent(iBin,*(y+0));	
		hSBinL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin),0,1,1);
		hSBoutL->SetBinContent(iBin,*(y+0));	
		hSBoutL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin),0,2,1);
		hYieldL->SetBinContent(iBin,*(y+0));	
		hYieldL->SetBinError(iBin,*(y+1));
		y = ExtractYield(0,mMasses->CopyTree(Form("lPt>%f&&lPt<%f&&lPart==1",xBins[iBin-1],xBins[iBin])),3,1);
		hYieldUBL->SetBinContent(iBin,*(y+0));	
		hYieldUBL->SetBinError(iBin,*(y+1));

		y = ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin),0,0,1);
		hSBinAL->SetBinContent(iBin,*(y+0));	
		hSBinAL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin),0,1,1);
		hSBoutAL->SetBinContent(iBin,*(y+0));	
		hSBoutAL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin),0,2,1);
		hYieldAL->SetBinContent(iBin,*(y+0));	
		hYieldAL->SetBinError(iBin,*(y+1));
		y = ExtractYield(0,mMasses->CopyTree(Form("lPt>%f&&lPt<%f&&lPart==2",xBins[iBin-1],xBins[iBin])),3,1);
		hYieldUBAL->SetBinContent(iBin,*(y+0));	
		hYieldUBAL->SetBinError(iBin,*(y+1));
		//hYieldK0s->SetBinContent(iBin,*(ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin))+0));
		//hYieldL->SetBinContent(iBin,*(ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin))+0));
		//hYieldAL->SetBinContent(iBin,*(ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin))+0));
	}

	hV0_PtK0s->Scale(1,"width");
	hV0_PtL->Scale(1,"width");
	hV0_PtAL->Scale(1,"width");
	hYieldK0s->Scale(1,"width");
	hYieldL->Scale(1,"width");
	hYieldAL->Scale(1,"width");
	hYieldUBK0s->Scale(1,"width");
	hYieldUBL->Scale(1,"width");
	hYieldUBAL->Scale(1,"width");
	hV0_DHasTPC->Divide(hV0_DPt);
	hV0_DHasTOF->Divide(hV0_DPt);
	hSBinK0s->Add(hSBoutK0s,-1.);
	hSBinK0s->Scale(1,"width");
	hSBinL->Add(hSBoutL,-1.);
	hSBinL->Scale(1,"width");
	hSBinAL->Add(hSBoutAL,-1.);
	hSBinAL->Scale(1,"width");

	if (flagMC) {
		hMCV0_PtK0s->Scale(1,"width");
		hMCV0_PtL->Scale(1,"width");
		hMCV0_PtAL->Scale(1,"width");
		hRCV0_PtK0s->Scale(1,"width");
		hRCV0_PtL->Scale(1,"width");
		hRCV0_PtAL->Scale(1,"width");
	}

	hV0_DHasTPC->SetTitle("HasTPC PID; p_{T} (GeV/#it{c}); Efficiency");
	hV0_DHasTOF->SetTitle("HasTPC PID; p_{T} (GeV/#it{c}); Efficiency");
	hV0_PtK0s->SetTitle("; p_{T} (GeV/#it{c}); K_{0}^{s} yield");
	hV0_PtL->SetTitle("; p_{T} (GeV/#it{c}); #Lambda yield");
	hV0_PtAL->SetTitle("; p_{T} (GeV/#it{c}); #bar{#Lambda} yield");

	hV0_PtK0s->SetLineColor(kRed);
	hV0_PtL->SetLineColor(kRed);
	hV0_PtAL->SetLineColor(kRed);
	hSBinK0s->SetLineColor(kGreen+2);
	hSBinL->SetLineColor(kGreen+2);
	hSBinAL->SetLineColor(kGreen+2);
	hYieldUBK0s->SetLineColor(kMagenta);
	hYieldUBL->SetLineColor(kMagenta);
	hYieldUBAL->SetLineColor(kMagenta);

	hV0_PtK0s->SetMarkerStyle(20);
	hSBinK0s->SetMarkerStyle(20);
	hYieldK0s->SetMarkerStyle(20);
	hYieldUBK0s->SetMarkerStyle(20);
	hV0_PtK0s->SetMarkerColor(kRed);
	hSBinK0s->SetMarkerColor(kGreen+2);
	hYieldK0s->SetMarkerColor(kBlue);
	hYieldUBK0s->SetMarkerColor(kMagenta);

	hV0_PtL->SetMarkerStyle(20);
	hSBinL->SetMarkerStyle(20);
	hYieldL->SetMarkerStyle(20);
	hYieldUBL->SetMarkerStyle(20);
	hV0_PtL->SetMarkerColor(kRed);
	hSBinL->SetMarkerColor(kGreen+2);
	hYieldL->SetMarkerColor(kBlue);
	hYieldUBL->SetMarkerColor(kMagenta);

	hV0_PtAL->SetMarkerStyle(20);
	hSBinAL->SetMarkerStyle(20);
	hYieldAL->SetMarkerStyle(20);
	hYieldUBAL->SetMarkerStyle(20);
	hV0_PtAL->SetMarkerColor(kRed);
	hSBinAL->SetMarkerColor(kGreen+2);
	hYieldAL->SetMarkerColor(kBlue);
	hYieldUBAL->SetMarkerColor(kMagenta);

	TString path = Form("pics_%s/",outputFile);// ("$HOME/sq/pics/");//("$HOME/sq/pics/");
	path.ReplaceAll(".root","");
	gSystem->Exec(Form("mkdir %s", path.Data()));
	cFits[0]->SaveAs(path+"f_k0s.png");
	cFits[1]->SaveAs(path+"f_ubk0s.png");
	cFits[2]->SaveAs(path+"f_l.png");
	cFits[3]->SaveAs(path+"f_ubl.png");
	cFits[4]->SaveAs(path+"f_al.png");
	cFits[5]->SaveAs(path+"f_ubal.png");
	cFits[0]->SaveAs(path+"f_k0s.png");

	TCanvas* can1 = new TCanvas("can1","",1000,700);
	hV0_DHasTOF->Draw();
	can1->SaveAs(path+"eff_tof.png");
	hV0_DHasTPC->Draw();
	can1->SaveAs(path+"eff_tpc.png");

	hRCV0_PtK0s->Divide(hMCV0_PtK0s);
	hRCV0_PtK0s->SetTitle(";p_{T} (GeV/#it{c});K_{S}^{0} reconstruction efficiency");
	hRCV0_PtK0s->GetYaxis()->SetRangeUser(0.,1.1);
	hRCV0_PtK0s->SetLineColor(kRed);
	hRCV0_PtK0s->SetMarkerStyle(20);
	hRCV0_PtK0s->Draw();
	can1->SaveAs(path+"effrec_k0s.png");
	hRCV0_PtL->Divide(hMCV0_PtL);
	hRCV0_PtL->SetTitle(";p_{T} (GeV/#it{c});#Lambda reconstruction efficiency");
	hRCV0_PtL->GetYaxis()->SetRangeUser(0.,1.1);
	hRCV0_PtL->SetLineColor(kRed);
	hRCV0_PtL->SetMarkerStyle(20);
	hRCV0_PtL->Draw();
	can1->SaveAs(path+"effrec_l.png");
	hRCV0_PtAL->Divide(hMCV0_PtAL);
	hRCV0_PtAL->SetTitle(";p_{T} (GeV/#it{c});#bar{#Lambda} reconstruction efficiency");
	hRCV0_PtAL->GetYaxis()->SetRangeUser(0.,1.1);
	hRCV0_PtAL->SetLineColor(kRed);
	hRCV0_PtAL->SetMarkerStyle(20);
	hRCV0_PtAL->Draw();
	can1->SaveAs(path+"effrec_al.png");

	TLegend *legpt = new TLegend(0.6,0.55,0.9,0.75);
	myLegendSetUp(legpt,0.028,1);
	legpt->AddEntry(hV0_PtK0s,"bin count in fixed range","l");
	legpt->AddEntry(hYieldK0s,"gaus+gaus+pol1 fit","l");
	legpt->AddEntry(hYieldUBK0s,"gaus+gaus+pol1 ub fit","l");
	legpt->AddEntry(hSBinK0s,"sideband with 3 and 6 RMS","l");

	//make effi plots
			
	hV0_PtK0s->Draw();
	can1->SetLogy();
	hYieldK0s->Draw("same");
	hYieldUBK0s->Draw("same");
	hSBinK0s->Draw("same");
	legpt->Draw();
	can1->SaveAs(path+"pt_k0s.png");
	hV0_PtL->Draw();
	hYieldL->Draw("same");
	hYieldUBL->Draw("same");
	hSBinL->Draw("same");
	legpt->Draw();
	can1->SaveAs(path+"pt_l.png");
	hV0_PtAL->Draw();
	hYieldAL->Draw("same");
	hYieldUBAL->Draw("same");
	hSBinAL->Draw("same");
	legpt->Draw();
	can1->SaveAs(path+"pt_al.png");

	TH1D* hpy = hV0_DTpcvpK0spi->ProjectionY();
	hpy->SetTitle(";n#sigma_{TPC}^{#pi};Entries");
	hpy->Draw();
	can1->SaveAs(path+"eff_tpck0spi.png");
	hpy = hV0_DTofvpK0spi->ProjectionY();
	hpy->SetTitle(";n#sigma_{TOF}^{#pi};Entries");
	hpy->Draw();
	can1->SaveAs(path+"eff_tofk0spi.png");
	hpy = hV0_DTpcvpLpr->ProjectionY();
	hpy->SetTitle(";n#sigma_{TPC}^{p};Entries");
	hpy->Draw();
	can1->SaveAs(path+"eff_tpclpr.png");
	hpy = hV0_DTofvpLpr->ProjectionY();
	hpy->SetTitle(";n#sigma_{TOF}^{p};Entries");
	hpy->Draw();
	can1->SaveAs(path+"eff_toflpr.png");

	printf(" WHAT IS UP \n", );
	printf(" MC FLAG is %i \n", (int)flagMC);

	// WRITING OBJECTS TO OUTPUT FILE
	if (outputFile!="")	mFout = new TFile(outputFile,"RECREATE");
	iHist = 0; while (lHist->At(iHist)) {			// should use an iterator...
		TString objName(lHist->At(iHist)->GetName());
		if (objName.BeginsWith("h")) lHist->At(iHist)->Write();
		iHist++;
	}	// can be replaced with embed->GetHistList()->Write(); ?
	mMasses->Write();
}