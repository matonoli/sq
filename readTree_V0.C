// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TLegend.h>

/*#include <AliAnalysisPIDEvent.h>
#include <AliAnalysisPIDTrack.h>
#include <AliAnalysisPIDV0.h>
#include <AliAnalysisPIDParticle.h>*/

#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooPlot.h>

#include <iostream>
#include <fstream>

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

// DATA V0 QA
TH3F* hV0_PtPhiEtaK0s 			= new TH3F("hV0_PtPhiEtaK0s","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hV0_PtPhiEtaL 			= new TH3F("hV0_PtPhiEtaL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hV0_PtPhiEtaAL 			= new TH3F("hV0_PtPhiEtaAL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH1F* hV0_DCAddK0s 				= new TH1F("hV0_DCAddK0s","",400,0.,2.);
TH1F* hV0_DCAddL 				= new TH1F("hV0_DCAddL","",400,0.,2.);
TH1F* hV0_DCAddAL 				= new TH1F("hV0_DCAddAL","",400,0.,2.);
TH1F* hV0_CPAK0s 				= new TH1F("hV0_CPAK0s","",400,0.98,1.02);
TH1F* hV0_CPAL 					= new TH1F("hV0_CPAL","",400,0.98,1.02);
TH1F* hV0_CPAAL 				= new TH1F("hV0_CPAAL","",400,0.98,1.02);
TH1F* hV0_RK0s 					= new TH1F("hV0_RK0s","",400,0,150);
TH1F* hV0_RL 					= new TH1F("hV0_RL","",400,0,150);
TH1F* hV0_RAL 					= new TH1F("hV0_RAL","",400,0,150);
TH1F* hV0_APK0s 				= new TH1F("hV0_APK0s","",400,0,5);
TH1F* hV0_APL 					= new TH1F("hV0_APL","",400,0,5);
TH1F* hV0_APAL 					= new TH1F("hV0_APAL","",400,0,5);

// DATA DAUGHTERS QA
TH2F* hV0_DDedxvpK0s 			= new TH2F("hV0_DDedxvpK0s","",450,0,15,300,0,300);
TH2F* hV0_DDedxvpL 				= new TH2F("hV0_DDedxvpL","",450,0,15,300,0,300);
TH2F* hV0_DDedxvpAL 			= new TH2F("hV0_DDedxvpAL","",450,0,15,300,0,300);
TH1F* hV0_DTpcK0spi				= new TH1F("hV0_DTpcK0spi","",400,-10.,10.);
TH1F* hV0_DTpcLpi				= new TH1F("hV0_DTpcLpi","",400,-10.,10.);
TH1F* hV0_DTpcALpi				= new TH1F("hV0_DTpcALpi","",400,-10.,10.);
TH1F* hV0_DTpcLpr				= new TH1F("hV0_DTpcLpr","",400,-10.,10.);
TH1F* hV0_DTpcALpr				= new TH1F("hV0_DTpcALpr","",400,-10.,10.);
TH1F* hV0_DNDedxK0spi			= new TH1F("hV0_DNDedxK0spi","",200,0,200);
TH1F* hV0_DNDedxLpi				= new TH1F("hV0_DNDedxLpi","",200,0,200);
TH1F* hV0_DNDedxALpi			= new TH1F("hV0_DNDedxALpi","",200,0,200);
TH1F* hV0_DNDedxLpr				= new TH1F("hV0_DNDedxLpr","",200,0,200);
TH1F* hV0_DNDedxALpr			= new TH1F("hV0_DNDedxALpr","",200,0,200);
TH2F* hV0_DTofBetavpK0s 		= new TH2F("hV0_DTofBetavpK0s","",450,0,15,300,0,1.2);
TH2F* hV0_DTofBetavpL 			= new TH2F("hV0_DTofBetavpL","",450,0,15,300,0,1.2);
TH2F* hV0_DTofBetavpAL 			= new TH2F("hV0_DTofBetavpAL","",450,0,15,300,0,1.2);
TH1F* hV0_DTofK0spi				= new TH1F("hV0_DTofK0spi","",400,-10.,10.);
TH1F* hV0_DTofLpi				= new TH1F("hV0_DTofLpi","",400,-10.,10.);
TH1F* hV0_DTofALpi				= new TH1F("hV0_DTofALpi","",400,-10.,10.);
TH1F* hV0_DTofLpr				= new TH1F("hV0_DTofLpr","",400,-10.,10.);
TH1F* hV0_DTofALpr				= new TH1F("hV0_DTofALpr","",400,-10.,10.);
TH1F* hV0_DIPxyK0spi			= new TH1F("hV0_DIPxyK0spi","",400,0.,20.);
TH1F* hV0_DIPxyLpi				= new TH1F("hV0_DIPxyLpi","",400,0.,20.);
TH1F* hV0_DIPxyALpi				= new TH1F("hV0_DIPxyALpi","",400,0.,20.);
TH1F* hV0_DIPxyLpr				= new TH1F("hV0_DIPxyLpr","",400,0.,20.);
TH1F* hV0_DIPxyALpr				= new TH1F("hV0_DIPxyALpr","",400,0.,20.);
TH1F* hV0_DIPzK0spi				= new TH1F("hV0_DIPzK0spi","",400,0.,20.);
TH1F* hV0_DIPzLpi				= new TH1F("hV0_DIPzLpi","",400,0.,20.);
TH1F* hV0_DIPzALpi				= new TH1F("hV0_DIPzALpi","",400,0.,20.);
TH1F* hV0_DIPzLpr				= new TH1F("hV0_DIPzLpr","",400,0.,20.);
TH1F* hV0_DIPzALpr				= new TH1F("hV0_DIPzALpr","",400,0.,20.);
TH1F* hV0_DNclusK0spi 			= new TH1F("hV0_DNclusK0spi","",200,0,200);
TH1F* hV0_DNclusLpi 			= new TH1F("hV0_DNclusLpi","",200,0,200);
TH1F* hV0_DNclusALpi 			= new TH1F("hV0_DNclusALpi","",200,0,200);
TH1F* hV0_DNclusLpr 			= new TH1F("hV0_DNclusLpr","",200,0,200);
TH1F* hV0_DNclusALpr 			= new TH1F("hV0_DNclusALpr","",200,0,200);
TH2F* hV0_DTPCvphiK0spi			= new TH2F("hV0_DTPCvphiK0spi","",400,-0.5,6.5,400,-10.,10.);
TH2F* hV0_DTPCvphiLpr			= new TH2F("hV0_DTPCvphiLpr","",400,-0.5,6.5,400,-10.,10.);
TH2F* hV0_DTPCvphiALpr			= new TH2F("hV0_DTPCvphiALpr","",400,-0.5,6.5,400,-10.,10.);
TH2F* hV0_DTPCvetaK0spi			= new TH2F("hV0_DTPCvetaK0spi","",400,-1.,1.,400,-10.,10.);
TH2F* hV0_DTPCvetaLpr			= new TH2F("hV0_DTPCvetaLpr","",400,-1.,1.,400,-10.,10.);
TH2F* hV0_DTPCvetaALpr			= new TH2F("hV0_DTPCvetaALpr","",400,-1.,1.,400,-10.,10.);
TH3F* hV0_DPtPhiEtaK0spi		= new TH3F("hV0_DPtPhiEtaK0spi","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hV0_DPtPhiEtaLpr			= new TH3F("hV0_DPtPhiEtaLpr","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hV0_DPtPhiEtaLpi			= new TH3F("hV0_DPtPhiEtaLpi","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hV0_DPtPhiEtaALpr			= new TH3F("hV0_DPtPhiEtaALpr","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hV0_DPtPhiEtaALpi			= new TH3F("hV0_DPtPhiEtaALpi","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH1F* hV0_DHasTPCK0spi			= new TH1F("hV0_DHasTPCK0spi","",200,0,20);
TH1F* hV0_DHasTPCLpr			= new TH1F("hV0_DHasTPCLpr","",200,0,20);
TH1F* hV0_DHasTPCLpi			= new TH1F("hV0_DHasTPCLpi","",200,0,20);
TH1F* hV0_DHasTPCALpr			= new TH1F("hV0_DHasTPCALpr","",200,0,20);
TH1F* hV0_DHasTPCALpi			= new TH1F("hV0_DHasTPCALpi","",200,0,20);
TH1F* hV0_DHasTOFK0spi			= new TH1F("hV0_DHasTOFK0spi","",200,0,20);
TH1F* hV0_DHasTOFLpr			= new TH1F("hV0_DHasTOFLpr","",200,0,20);
TH1F* hV0_DHasTOFLpi			= new TH1F("hV0_DHasTOFLpi","",200,0,20);
TH1F* hV0_DHasTOFALpr			= new TH1F("hV0_DHasTOFALpr","",200,0,20);
TH1F* hV0_DHasTOFALpi			= new TH1F("hV0_DHasTOFALpi","",200,0,20);

// FAKE SIGNAL QA
TH1F* hFV0_PDGsK0s 				= new TH1F("hFV0_PDGsK0s","",20000,-10000,10000);
TH1F* hFV0_PDGsL 				= new TH1F("hFV0_PDGsL","",20000,-10000,10000);
TH1F* hFV0_PDGsAL 				= new TH1F("hFV0_PDGsAL","",20000,-10000,10000);
TH1F* hFV0_IMK0s				= new TH1F("hFV0_IMK0s","",2000,-0.2,0.2);		
TH1F* hFV0_IML					= new TH1F("hFV0_IML","",2000,-0.2,0.2);
TH1F* hFV0_IMAL					= new TH1F("hFV0_IMAL","",2000,-0.2,0.2);
TH3F* hFV0_PtPhiEtaK0s 			= new TH3F("hFV0_PtPhiEtaK0s","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hFV0_PtPhiEtaL 			= new TH3F("hFV0_PtPhiEtaL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hFV0_PtPhiEtaAL 			= new TH3F("hFV0_PtPhiEtaAL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH1F* hFV0_DCAddK0s 			= new TH1F("hFV0_DCAddK0s","",400,0.,2.);
TH1F* hFV0_DCAddL 				= new TH1F("hFV0_DCAddL","",400,0.,2.);
TH1F* hFV0_DCAddAL 				= new TH1F("hFV0_DCAddAL","",400,0.,2.);
TH1F* hFV0_CPAK0s 				= new TH1F("hFV0_CPAK0s","",400,0.98,1.02);
TH1F* hFV0_CPAL 				= new TH1F("hFV0_CPAL","",400,0.98,1.02);
TH1F* hFV0_CPAAL 				= new TH1F("hFV0_CPAAL","",400,0.98,1.02);
TH1F* hFV0_RK0s 				= new TH1F("hFV0_RK0s","",400,0,150);
TH1F* hFV0_RL 					= new TH1F("hFV0_RL","",400,0,150);
TH1F* hFV0_RAL 					= new TH1F("hFV0_RAL","",400,0,150);
TH1F* hFV0_APK0s 				= new TH1F("hFV0_APK0s","",400,0,5);
TH1F* hFV0_APL 					= new TH1F("hFV0_APL","",400,0,5);
TH1F* hFV0_APAL 				= new TH1F("hFV0_APAL","",400,0,5);
TH2F* hFV0_DDedxvpK0s 			= new TH2F("hFV0_DDedxvpK0s","",450,0,15,300,0,300);
TH2F* hFV0_DDedxvpL 			= new TH2F("hFV0_DDedxvpL","",450,0,15,300,0,300);
TH2F* hFV0_DDedxvpAL 			= new TH2F("hFV0_DDedxvpAL","",450,0,15,300,0,300);
TH1F* hFV0_DTpcK0spi			= new TH1F("hFV0_DTpcK0spi","",400,-10.,10.);
TH1F* hFV0_DTpcLpi				= new TH1F("hFV0_DTpcLpi","",400,-10.,10.);
TH1F* hFV0_DTpcALpi				= new TH1F("hFV0_DTpcALpi","",400,-10.,10.);
TH1F* hFV0_DTpcLpr				= new TH1F("hFV0_DTpcLpr","",400,-10.,10.);
TH1F* hFV0_DTpcALpr				= new TH1F("hFV0_DTpcALpr","",400,-10.,10.);
TH1F* hFV0_DNDedxK0spi			= new TH1F("hFV0_DNDedxK0spi","",200,0,200);
TH1F* hFV0_DNDedxLpi			= new TH1F("hFV0_DNDedxLpi","",200,0,200);
TH1F* hFV0_DNDedxALpi			= new TH1F("hFV0_DNDedxALpi","",200,0,200);
TH1F* hFV0_DNDedxLpr			= new TH1F("hFV0_DNDedxLpr","",200,0,200);
TH1F* hFV0_DNDedxALpr			= new TH1F("hFV0_DNDedxALpr","",200,0,200);
TH2F* hFV0_DTofBetavpK0s 		= new TH2F("hFV0_DTofBetavpK0s","",450,0,15,300,0,1.2);
TH2F* hFV0_DTofBetavpL 			= new TH2F("hFV0_DTofBetavpL","",450,0,15,300,0,1.2);
TH2F* hFV0_DTofBetavpAL 		= new TH2F("hFV0_DTofBetavpAL","",450,0,15,300,0,1.2);
TH1F* hFV0_DTofK0spi			= new TH1F("hFV0_DTofK0spi","",400,-10.,10.);
TH1F* hFV0_DTofLpi				= new TH1F("hFV0_DTofLpi","",400,-10.,10.);
TH1F* hFV0_DTofALpi				= new TH1F("hFV0_DTofALpi","",400,-10.,10.);
TH1F* hFV0_DTofLpr				= new TH1F("hFV0_DTofLpr","",400,-10.,10.);
TH1F* hFV0_DTofALpr				= new TH1F("hFV0_DTofALpr","",400,-10.,10.);
TH1F* hFV0_DIPxyK0spi			= new TH1F("hFV0_DIPxyK0spi","",400,0.,20.);
TH1F* hFV0_DIPxyLpi				= new TH1F("hFV0_DIPxyLpi","",400,0.,20.);
TH1F* hFV0_DIPxyALpi			= new TH1F("hFV0_DIPxyALpi","",400,0.,20.);
TH1F* hFV0_DIPxyLpr				= new TH1F("hFV0_DIPxyLpr","",400,0.,20.);
TH1F* hFV0_DIPxyALpr			= new TH1F("hFV0_DIPxyALpr","",400,0.,20.);
TH1F* hFV0_DIPzK0spi			= new TH1F("hFV0_DIPzK0spi","",400,0.,20.);
TH1F* hFV0_DIPzLpi				= new TH1F("hFV0_DIPzLpi","",400,0.,20.);
TH1F* hFV0_DIPzALpi				= new TH1F("hFV0_DIPzALpi","",400,0.,20.);
TH1F* hFV0_DIPzLpr				= new TH1F("hFV0_DIPzLpr","",400,0.,20.);
TH1F* hFV0_DIPzALpr				= new TH1F("hFV0_DIPzALpr","",400,0.,20.);
TH1F* hFV0_DNclusK0spi 			= new TH1F("hFV0_DNclusK0spi","",200,0,200);
TH1F* hFV0_DNclusLpi 			= new TH1F("hFV0_DNclusLpi","",200,0,200);
TH1F* hFV0_DNclusALpi 			= new TH1F("hFV0_DNclusALpi","",200,0,200);
TH1F* hFV0_DNclusLpr 			= new TH1F("hFV0_DNclusLpr","",200,0,200);
TH1F* hFV0_DNclusALpr 			= new TH1F("hFV0_DNclusALpr","",200,0,200);

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
TH3F* hMCV0_PtPhiEtaK0s 		= new TH3F("hMCV0_PtPhiEtaK0s","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hMCV0_PtPhiEtaL 			= new TH3F("hMCV0_PtPhiEtaL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hMCV0_PtPhiEtaAL 			= new TH3F("hMCV0_PtPhiEtaAL","",200,0,20,200,-0.5,6.5,200,-1.,1.);
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
TH3F* hRCV0_DPtPhiEtaK0spi		= new TH3F("hRCV0_DPtPhiEtaK0spi","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hRCV0_DPtPhiEtaLpr		= new TH3F("hRCV0_DPtPhiEtaLpr","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hRCV0_DPtPhiEtaLpi		= new TH3F("hRCV0_DPtPhiEtaLpi","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hRCV0_DPtPhiEtaALpr		= new TH3F("hRCV0_DPtPhiEtaALpr","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH3F* hRCV0_DPtPhiEtaALpi		= new TH3F("hRCV0_DPtPhiEtaALpi","",200,0,20,200,-0.5,6.5,200,-1.,1.);
TH1F* hRCV0_DHasTPCK0spi		= new TH1F("hRCV0_DHasTPCK0spi","",200,0,20);
TH1F* hRCV0_DHasTPCLpr			= new TH1F("hRCV0_DHasTPCLpr","",200,0,20);
TH1F* hRCV0_DHasTPCLpi			= new TH1F("hRCV0_DHasTPCLpi","",200,0,20);
TH1F* hRCV0_DHasTPCALpr			= new TH1F("hRCV0_DHasTPCALpr","",200,0,20);
TH1F* hRCV0_DHasTPCALpi			= new TH1F("hRCV0_DHasTPCALpi","",200,0,20);
TH1F* hRCV0_DHasTOFK0spi		= new TH1F("hRCV0_DHasTOFK0spi","",200,0,20);
TH1F* hRCV0_DHasTOFLpr			= new TH1F("hRCV0_DHasTOFLpr","",200,0,20);
TH1F* hRCV0_DHasTOFLpi			= new TH1F("hRCV0_DHasTOFLpi","",200,0,20);
TH1F* hRCV0_DHasTOFALpr			= new TH1F("hRCV0_DHasTOFALpr","",200,0,20);
TH1F* hRCV0_DHasTOFALpi			= new TH1F("hRCV0_DHasTOFALpi","",200,0,20);
TH2F* hRCV0_DTPCvphiK0spi		= new TH2F("hRCV0_DTPCvphiK0spi","",400,-0.5,6.5,400,-10.,10.);
TH2F* hRCV0_DTPCvphiLpr			= new TH2F("hRCV0_DTPCvphiLpr","",400,-0.5,6.5,400,-10.,10.);
TH2F* hRCV0_DTPCvphiALpr		= new TH2F("hRCV0_DTPCvphiALpr","",400,-0.5,6.5,400,-10.,10.);
TH2F* hRCV0_DTPCvetaK0spi		= new TH2F("hRCV0_DTPCvetaK0spi","",400,-1.,1.,400,-10.,10.);
TH2F* hRCV0_DTPCvetaLpr			= new TH2F("hRCV0_DTPCvetaLpr","",400,-1.,1.,400,-10.,10.);
TH2F* hRCV0_DTPCvetaALpr		= new TH2F("hRCV0_DTPCvetaALpr","",400,-1.,1.,400,-10.,10.);

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
	Double_t* ap = CalculateAP(v0);

	if (*(ap+1) < 0.2*fabs(*(ap+0)))	return false;
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
	Float_t bugcheck = -99; Float_t bugcheck2 = -99;
	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		//if (iEv!=139254) continue;

		hEventMonitor->Fill(0);
		if (iEv%10000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		bTracks->Clear();
		mChain->GetEntry(iEv);
		if (!mEvent) continue;

		AliAnalysisPIDV0* bugchecker = (AliAnalysisPIDV0*)bV0s->At(0);
		if (bugchecker) {
			if (fabs(bugchecker->GetRadius() - bugcheck) < 0.00001 &&
				fabs(bugchecker->GetPt() - bugcheck2) < 0.001)  {
				//printf("1 this %f past %f \n", bugchecker->GetRadius(), bugcheck);
				//printf("2 this %f past %f \n", bugchecker->GetPt(), bugcheck2);
				//printf("identified doubled event %i , skipping \n", iEv);
				continue; }
			bugcheck =  bugchecker->GetRadius();
			bugcheck2 =  bugchecker->GetPt();
		}
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
						hMCV0_PtPhiEtaK0s->Fill(p->GetPt(), p->GetPhi(), p->GetEta());
						hRCV0_PtPhiEtaK0s->Fill(v0rc->GetPt(), phirc, v0rc->GetEta());
						if (trPrc->HasTPCPID()) {
							hRCV0_DHasTPCK0spi->Fill(trPrc->GetPt());
							hRCV0_DDedxvpK0s->Fill(trPrc->GetP(),trPrc->GetTPCdEdx());
							hRCV0_DTpcK0spi->Fill(trPrc->GetNSigmaPionTPC());
							hRCV0_DNDedxK0spi->Fill(trPrc->GetTPCdEdxN());
							hRCV0_DTPCvphiK0spi->Fill(trPrc->GetPhi(),trPrc->GetNSigmaPionTPC()); 
							hRCV0_DTPCvetaK0spi->Fill(trPrc->GetEta(),trPrc->GetNSigmaPionTPC()); 			}
						if (trNrc->HasTPCPID()) {
							hRCV0_DHasTPCK0spi->Fill(trNrc->GetPt());
							hRCV0_DDedxvpK0s->Fill(trNrc->GetP(),trNrc->GetTPCdEdx());
							hRCV0_DTpcK0spi->Fill(trNrc->GetNSigmaPionTPC());
							hRCV0_DNDedxK0spi->Fill(trNrc->GetTPCdEdxN());
							hRCV0_DTPCvphiK0spi->Fill(trNrc->GetPhi(),trNrc->GetNSigmaPionTPC()); 
							hRCV0_DTPCvetaK0spi->Fill(trNrc->GetEta(),trNrc->GetNSigmaPionTPC()); 				}
						if (trPrc->HasTOFPID()) {
							hRCV0_DHasTOFK0spi->Fill(trPrc->GetPt());
							hRCV0_DTofBetavpK0s->Fill(trPrc->GetP(),trPrc->GetTOFExpBeta(AliPID::kPion));
							//hRCV0_DTofBetavpK0s->Fill(trPrc->GetP(),trPrc->GetTOFBeta());
							hRCV0_DTofK0spi->Fill(trPrc->GetNSigmaPionTOF());
							hMCvRC_DTofTimeK0spi->Fill(trPrc->GetMCTOFTime(),trPrc->GetTOFTime());
							hMCvRC_DTofLengthK0spi->Fill(trPrc->GetMCTOFLength(),trPrc->GetTOFLength());	}
						if (trNrc->HasTOFPID()) {
							hRCV0_DHasTOFK0spi->Fill(trNrc->GetPt());
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
						hRCV0_DPtPhiEtaK0spi->Fill(trPrc->GetPt(),trPrc->GetPhi(),trPrc->GetEta());
						hRCV0_DPtPhiEtaK0spi->Fill(trNrc->GetPt(),trNrc->GetPhi(),trNrc->GetEta());
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
						hMCV0_PtPhiEtaL->Fill(p->GetPt(), p->GetPhi(), p->GetEta());
						hRCV0_PtPhiEtaL->Fill(v0rc->GetPt(), phirc, v0rc->GetEta());
						if (trPrc->HasTPCPID()) {
							hRCV0_DHasTPCLpr->Fill(trPrc->GetPt());
							hRCV0_DDedxvpL->Fill(trPrc->GetP(),trPrc->GetTPCdEdx());
							hRCV0_DTpcLpr->Fill(trPrc->GetNSigmaProtonTPC());
							hRCV0_DNDedxLpr->Fill(trPrc->GetTPCdEdxN());
							hRCV0_DTPCvphiLpr->Fill(trPrc->GetPhi(),trPrc->GetNSigmaProtonTPC()); 
							hRCV0_DTPCvetaLpr->Fill(trPrc->GetEta(),trPrc->GetNSigmaProtonTPC()); 			}
						if (trNrc->HasTPCPID()) {
							hRCV0_DHasTPCLpi->Fill(trNrc->GetPt());
							hRCV0_DDedxvpL->Fill(trNrc->GetP(),trNrc->GetTPCdEdx());
							hRCV0_DTpcLpi->Fill(trNrc->GetNSigmaPionTPC());
							hRCV0_DNDedxLpi->Fill(trNrc->GetTPCdEdxN()); 								}
						if (trPrc->HasTOFPID()) {
							hRCV0_DHasTOFLpr->Fill(trPrc->GetPt());
							hRCV0_DTofBetavpL->Fill(trPrc->GetP(),trPrc->GetTOFExpBeta(AliPID::kProton));
							hRCV0_DTofLpr->Fill(trPrc->GetNSigmaProtonTOF());
							hMCvRC_DTofTimeLpr->Fill(trPrc->GetMCTOFTime(),trPrc->GetTOFTime());
							hMCvRC_DTofLengthLpr->Fill(trPrc->GetMCTOFLength(),trPrc->GetTOFLength());	}
						if (trNrc->HasTOFPID()) {
							hRCV0_DHasTOFLpi->Fill(trNrc->GetPt());
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
						hRCV0_DPtPhiEtaLpr->Fill(trPrc->GetPt(),trPrc->GetPhi(),trPrc->GetEta());
						hRCV0_DPtPhiEtaLpi->Fill(trNrc->GetPt(),trNrc->GetPhi(),trNrc->GetEta());
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
						hMCV0_PtPhiEtaAL->Fill(p->GetPt(), p->GetPhi(), p->GetEta());
						hRCV0_PtPhiEtaAL->Fill(v0rc->GetPt(), phirc, v0rc->GetEta());
						if (trPrc->HasTPCPID()) {
							hRCV0_DHasTPCALpi->Fill(trPrc->GetPt());
							hRCV0_DDedxvpAL->Fill(trPrc->GetP(),trPrc->GetTPCdEdx());
							hRCV0_DTpcALpi->Fill(trPrc->GetNSigmaPionTPC());
							hRCV0_DNDedxALpi->Fill(trPrc->GetTPCdEdxN()); 								}
						if (trNrc->HasTPCPID()) {
							hRCV0_DHasTPCALpr->Fill(trNrc->GetPt());
							hRCV0_DDedxvpAL->Fill(trNrc->GetP(),trNrc->GetTPCdEdx());
							hRCV0_DTpcALpr->Fill(trNrc->GetNSigmaProtonTPC());
							hRCV0_DNDedxALpr->Fill(trNrc->GetTPCdEdxN());
							hRCV0_DTPCvphiALpr->Fill(trNrc->GetPhi(),trNrc->GetNSigmaProtonTPC()); 
							hRCV0_DTPCvetaALpr->Fill(trNrc->GetEta(),trNrc->GetNSigmaProtonTPC()); 								}
						if (trPrc->HasTOFPID()) {
							hRCV0_DHasTOFALpi->Fill(trPrc->GetPt());
							hRCV0_DTofBetavpAL->Fill(trPrc->GetP(),trPrc->GetTOFExpBeta(AliPID::kPion));
							hRCV0_DTofALpi->Fill(trPrc->GetNSigmaPionTOF());
							hMCvRC_DTofTimeALpi->Fill(trPrc->GetMCTOFTime(),trPrc->GetTOFTime());
							hMCvRC_DTofLengthALpi->Fill(trPrc->GetMCTOFLength(),trPrc->GetTOFLength());	}
						if (trNrc->HasTOFPID()) {
							hRCV0_DHasTOFALpr->Fill(trNrc->GetPt());
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
						hRCV0_DPtPhiEtaALpi->Fill(trPrc->GetPt(),trPrc->GetPhi(),trPrc->GetEta());
						hRCV0_DPtPhiEtaALpr->Fill(trNrc->GetPt(),trNrc->GetPhi(),trNrc->GetEta());
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
			Double_t* ap = CalculateAP(v0);
			Double_t phi = CalculatePhi(v0);

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

			bool noCuts = 0; float masscut = 0.01;
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
					hV0_DTofvpK0spi->Fill(trN->GetP(),trN->GetNSigmaPionTOF());		

					hV0_PtPhiEtaK0s->Fill(v0->GetPt(), phi, v0->GetEta());
					hV0_DCAddK0s->Fill(v0->GetDCAV0Daughters());
					hV0_CPAK0s->Fill(v0->GetV0CosinePA());
					hV0_RK0s->Fill(v0->GetRadius());
					if (fabs(*(ap+0))>0.005) hV0_APK0s->Fill(*(ap+1)/fabs(*(ap+0)));
					if (trP->HasTPCPID()) {
							hV0_DHasTPCK0spi->Fill(trP->GetPt());
							hV0_DDedxvpK0s->Fill(trP->GetP(),trP->GetTPCdEdx());
							hV0_DTpcK0spi->Fill(trP->GetNSigmaPionTPC());
							hV0_DNDedxK0spi->Fill(trP->GetTPCdEdxN());
							hV0_DTPCvphiK0spi->Fill(trP->GetPhi(),trP->GetNSigmaPionTPC()); 
							hV0_DTPCvetaK0spi->Fill(trP->GetEta(),trP->GetNSigmaPionTPC());		}
						if (trN->HasTPCPID()) {
							hV0_DHasTPCK0spi->Fill(trN->GetPt());
							hV0_DDedxvpK0s->Fill(trN->GetP(),trN->GetTPCdEdx());
							hV0_DTpcK0spi->Fill(trN->GetNSigmaPionTPC());
							hV0_DNDedxK0spi->Fill(trN->GetTPCdEdxN());
							hV0_DTPCvphiK0spi->Fill(trN->GetPhi(),trN->GetNSigmaPionTPC()); 
							hV0_DTPCvetaK0spi->Fill(trN->GetEta(),trN->GetNSigmaPionTPC());		}
						if (trP->HasTOFPID()) {
							hV0_DHasTOFK0spi->Fill(trP->GetPt());
							hV0_DTofBetavpK0s->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kPion));
							//hV0_DTofBetavpK0s->Fill(trP->GetP(),trP->GetTOFBeta());
							hV0_DTofK0spi->Fill(trP->GetNSigmaPionTOF());	}
						if (trN->HasTOFPID()) {
							hV0_DHasTOFK0spi->Fill(trN->GetPt());
							hV0_DTofBetavpK0s->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kPion));
							hV0_DTofK0spi->Fill(trN->GetNSigmaPionTOF());	}
					hV0_DIPxyK0spi->Fill(trP->GetImpactParameter(0));
					hV0_DIPxyK0spi->Fill(trN->GetImpactParameter(0));
					hV0_DIPzK0spi->Fill(trP->GetImpactParameter(1));
					hV0_DIPzK0spi->Fill(trN->GetImpactParameter(1));
					hV0_DNclusK0spi->Fill(trP->GetTPCNcls());						
					hV0_DNclusK0spi->Fill(trN->GetTPCNcls());
					hV0_DPtPhiEtaK0spi->Fill(trP->GetPt(),trP->GetPhi(),trP->GetEta());
					hV0_DPtPhiEtaK0spi->Fill(trN->GetPt(),trN->GetPhi(),trN->GetEta());		}
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
					hV0_DTofvpLpr->Fill(trP->GetP(),trP->GetNSigmaProtonTOF()); 			

					hV0_PtPhiEtaL->Fill(v0->GetPt(), phi, v0->GetEta());
					hV0_DCAddL->Fill(v0->GetDCAV0Daughters());
					hV0_CPAL->Fill(v0->GetV0CosinePA());
					hV0_RL->Fill(v0->GetRadius());
					if (fabs(*(ap+0))>0.005) hV0_APL->Fill(*(ap+1)/fabs(*(ap+0)));
					if (trP->HasTPCPID()) {
							hV0_DHasTPCLpr->Fill(trP->GetPt());
							hV0_DDedxvpL->Fill(trP->GetP(),trP->GetTPCdEdx());
							hV0_DTpcLpr->Fill(trP->GetNSigmaProtonTPC());
							hV0_DNDedxLpr->Fill(trP->GetTPCdEdxN());
							hV0_DTPCvphiLpr->Fill(trP->GetPhi(),trP->GetNSigmaProtonTPC()); 
							hV0_DTPCvetaLpr->Fill(trP->GetEta(),trP->GetNSigmaProtonTPC());			}
						if (trN->HasTPCPID()) {
							hV0_DHasTPCLpi->Fill(trN->GetPt());
							hV0_DDedxvpL->Fill(trN->GetP(),trN->GetTPCdEdx());
							hV0_DTpcLpi->Fill(trN->GetNSigmaPionTPC());
							hV0_DNDedxLpi->Fill(trN->GetTPCdEdxN()); 								}
						if (trP->HasTOFPID()) {
							hV0_DHasTOFLpr->Fill(trP->GetPt());
							hV0_DTofBetavpL->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kProton));
							//hV0_DTofBetavpL->Fill(trP->GetP(),trP->GetTOFBeta());
							hV0_DTofLpr->Fill(trP->GetNSigmaProtonTOF());	}
						if (trN->HasTOFPID()) {
							hV0_DHasTOFLpi->Fill(trN->GetPt());
							hV0_DTofBetavpL->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kPion));
							hV0_DTofLpi->Fill(trN->GetNSigmaPionTOF());	}
					hV0_DIPxyLpr->Fill(trP->GetImpactParameter(0));
					hV0_DIPxyLpi->Fill(trN->GetImpactParameter(0));
					hV0_DIPzLpr->Fill(trP->GetImpactParameter(1));
					hV0_DIPzLpi->Fill(trN->GetImpactParameter(1));
					hV0_DNclusLpr->Fill(trP->GetTPCNcls());						
					hV0_DNclusLpi->Fill(trN->GetTPCNcls());
					hV0_DPtPhiEtaLpr->Fill(trP->GetPt(),trP->GetPhi(),trP->GetEta());
					hV0_DPtPhiEtaLpi->Fill(trN->GetPt(),trN->GetPhi(),trN->GetEta());		}
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
					hV0_DTofminTpcvrLpr->Fill(v0->GetRadius(),delta);

					hV0_PtPhiEtaAL->Fill(v0->GetPt(), phi, v0->GetEta());
					hV0_DCAddAL->Fill(v0->GetDCAV0Daughters());
					hV0_CPAAL->Fill(v0->GetV0CosinePA());
					hV0_RAL->Fill(v0->GetRadius());
					if (fabs(*(ap+0))>0.005) hV0_APAL->Fill(*(ap+1)/fabs(*(ap+0)));
					if (trP->HasTPCPID()) {
							hV0_DHasTPCALpi->Fill(trP->GetPt());
							hV0_DDedxvpAL->Fill(trP->GetP(),trP->GetTPCdEdx());
							hV0_DTpcALpi->Fill(trP->GetNSigmaPionTPC());
							hV0_DNDedxALpi->Fill(trP->GetTPCdEdxN()); 								}
						if (trN->HasTPCPID()) {
							hV0_DHasTPCALpr->Fill(trN->GetPt());
							hV0_DDedxvpAL->Fill(trN->GetP(),trN->GetTPCdEdx());
							hV0_DTpcALpr->Fill(trN->GetNSigmaProtonTPC());
							hV0_DNDedxALpr->Fill(trN->GetTPCdEdxN());
							hV0_DTPCvphiALpr->Fill(trN->GetPhi(),trN->GetNSigmaProtonTPC()); 
							hV0_DTPCvetaALpr->Fill(trN->GetEta(),trN->GetNSigmaProtonTPC());			}
						if (trP->HasTOFPID()) {
							hV0_DHasTOFALpi->Fill(trP->GetPt());
							hV0_DTofBetavpAL->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kPion));
							//hV0_DTofBetavpAL->Fill(trP->GetP(),trP->GetTOFBeta());
							hV0_DTofALpi->Fill(trP->GetNSigmaPionTOF());	}
						if (trN->HasTOFPID()) {
							hV0_DHasTOFALpr->Fill(trN->GetPt());
							hV0_DTofBetavpAL->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kProton));
							hV0_DTofALpr->Fill(trN->GetNSigmaProtonTOF());	}
					hV0_DIPxyALpi->Fill(trP->GetImpactParameter(0));
					hV0_DIPxyALpr->Fill(trN->GetImpactParameter(0));
					hV0_DIPzALpi->Fill(trP->GetImpactParameter(1));
					hV0_DIPzALpr->Fill(trN->GetImpactParameter(1));
					hV0_DNclusALpi->Fill(trP->GetTPCNcls());						
					hV0_DNclusALpr->Fill(trN->GetTPCNcls());
					hV0_DPtPhiEtaALpi->Fill(trP->GetPt(),trP->GetPhi(),trP->GetEta());
					hV0_DPtPhiEtaALpr->Fill(trN->GetPt(),trN->GetPhi(),trN->GetEta());	}
				}



			////////////////////////
			// MC STUDY OF FAKE V0s
			////////////////////////
			if (flagMC) {

				int mcPdg = v0->GetMCPdgCode();
				if (IsK0s(v0,cutFlag) && (mcPdg != PDG_IDS[0]) && (fabs(v0->GetIMK0s())<masscut)) {
					hFV0_PDGsK0s->Fill(mcPdg);
					hFV0_IMK0s->Fill(v0->GetIMK0s());

					hFV0_PtPhiEtaK0s->Fill(v0->GetPt(), phi, v0->GetEta());
					hFV0_DCAddK0s->Fill(v0->GetDCAV0Daughters());
					hFV0_CPAK0s->Fill(v0->GetV0CosinePA());
					hFV0_RK0s->Fill(v0->GetRadius());
					if (fabs(*(ap+0))>0.005) hFV0_APK0s->Fill(*(ap+1)/fabs(*(ap+0)));
					if (trP->HasTPCPID()) {
							hFV0_DDedxvpK0s->Fill(trP->GetP(),trP->GetTPCdEdx());
							hFV0_DTpcK0spi->Fill(trP->GetNSigmaPionTPC());
							hFV0_DNDedxK0spi->Fill(trP->GetTPCdEdxN()); 								}
						if (trN->HasTPCPID()) {
							hFV0_DDedxvpK0s->Fill(trN->GetP(),trN->GetTPCdEdx());
							hFV0_DTpcK0spi->Fill(trN->GetNSigmaPionTPC());
							hFV0_DNDedxK0spi->Fill(trN->GetTPCdEdxN()); 								}
						if (trP->HasTOFPID()) {
							hFV0_DTofBetavpK0s->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kPion));
							//hFV0_DTofBetavpK0s->Fill(trP->GetP(),trP->GetTOFBeta());
							hFV0_DTofK0spi->Fill(trP->GetNSigmaPionTOF());	}
						if (trN->HasTOFPID()) {
							hFV0_DTofBetavpK0s->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kPion));
							hFV0_DTofK0spi->Fill(trN->GetNSigmaPionTOF());	}
					hFV0_DIPxyK0spi->Fill(trP->GetImpactParameter(0));
					hFV0_DIPxyK0spi->Fill(trN->GetImpactParameter(0));
					hFV0_DIPzK0spi->Fill(trP->GetImpactParameter(1));
					hFV0_DIPzK0spi->Fill(trN->GetImpactParameter(1));
					hFV0_DNclusK0spi->Fill(trP->GetTPCNcls());						
					hFV0_DNclusK0spi->Fill(trN->GetTPCNcls());				}
				
				if (IsL(v0,cutFlag) && (mcPdg != PDG_IDS[1]) && (fabs(v0->GetIML())<masscut)) {
					hFV0_PDGsL->Fill(mcPdg);
					hFV0_IML->Fill(v0->GetIML());

					hFV0_PtPhiEtaL->Fill(v0->GetPt(), phi, v0->GetEta());
					hFV0_DCAddL->Fill(v0->GetDCAV0Daughters());
					hFV0_CPAL->Fill(v0->GetV0CosinePA());
					hFV0_RL->Fill(v0->GetRadius());
					if (fabs(*(ap+0))>0.005) hFV0_APL->Fill(*(ap+1)/fabs(*(ap+0)));
					if (trP->HasTPCPID()) {
							hFV0_DDedxvpL->Fill(trP->GetP(),trP->GetTPCdEdx());
							hFV0_DTpcLpr->Fill(trP->GetNSigmaProtonTPC());
							hFV0_DNDedxLpr->Fill(trP->GetTPCdEdxN()); 								}
						if (trN->HasTPCPID()) {
							hFV0_DDedxvpL->Fill(trN->GetP(),trN->GetTPCdEdx());
							hFV0_DTpcLpi->Fill(trN->GetNSigmaPionTPC());
							hFV0_DNDedxLpi->Fill(trN->GetTPCdEdxN()); 								}
						if (trP->HasTOFPID()) {
							hFV0_DTofBetavpL->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kProton));
							//hFV0_DTofBetavpL->Fill(trP->GetP(),trP->GetTOFBeta());
							hFV0_DTofLpr->Fill(trP->GetNSigmaProtonTOF());	}
						if (trN->HasTOFPID()) {
							hFV0_DTofBetavpL->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kPion));
							hFV0_DTofLpi->Fill(trN->GetNSigmaPionTOF());	}
					hFV0_DIPxyLpr->Fill(trP->GetImpactParameter(0));
					hFV0_DIPxyLpi->Fill(trN->GetImpactParameter(0));
					hFV0_DIPzLpr->Fill(trP->GetImpactParameter(1));
					hFV0_DIPzLpi->Fill(trN->GetImpactParameter(1));
					hFV0_DNclusLpr->Fill(trP->GetTPCNcls());						
					hFV0_DNclusLpi->Fill(trN->GetTPCNcls());				}

				if (IsAL(v0,cutFlag) && (mcPdg != PDG_IDS[2]) && (fabs(v0->GetIMAL())<masscut)) {
					hFV0_PDGsAL->Fill(mcPdg);
					hFV0_IMAL->Fill(v0->GetIMAL());

					hFV0_PtPhiEtaAL->Fill(v0->GetPt(), phi, v0->GetEta());
					hFV0_DCAddAL->Fill(v0->GetDCAV0Daughters());
					hFV0_CPAAL->Fill(v0->GetV0CosinePA());
					hFV0_RAL->Fill(v0->GetRadius());
					if (fabs(*(ap+0))>0.005) hFV0_APAL->Fill(*(ap+1)/fabs(*(ap+0)));
					if (trP->HasTPCPID()) {
							hFV0_DDedxvpAL->Fill(trP->GetP(),trP->GetTPCdEdx());
							hFV0_DTpcALpi->Fill(trP->GetNSigmaPionTPC());
							hFV0_DNDedxALpi->Fill(trP->GetTPCdEdxN()); 								}
						if (trN->HasTPCPID()) {
							hFV0_DDedxvpAL->Fill(trN->GetP(),trN->GetTPCdEdx());
							hFV0_DTpcALpr->Fill(trN->GetNSigmaProtonTPC());
							hFV0_DNDedxALpr->Fill(trN->GetTPCdEdxN()); 								}
						if (trP->HasTOFPID()) {
							hFV0_DTofBetavpAL->Fill(trP->GetP(),trP->GetTOFExpBeta(AliPID::kPion));
							//hFV0_DTofBetavpAL->Fill(trP->GetP(),trP->GetTOFBeta());
							hFV0_DTofALpi->Fill(trP->GetNSigmaPionTOF());	}
						if (trN->HasTOFPID()) {
							hFV0_DTofBetavpAL->Fill(trN->GetP(),trN->GetTOFExpBeta(AliPID::kProton));
							hFV0_DTofALpr->Fill(trN->GetNSigmaProtonTOF());	}
					hFV0_DIPxyALpi->Fill(trP->GetImpactParameter(0));
					hFV0_DIPxyALpr->Fill(trN->GetImpactParameter(0));
					hFV0_DIPzALpi->Fill(trP->GetImpactParameter(1));
					hFV0_DIPzALpr->Fill(trN->GetImpactParameter(1));
					hFV0_DNclusALpi->Fill(trP->GetTPCNcls());						
					hFV0_DNclusALpr->Fill(trN->GetTPCNcls());				}

			}	//fake mc study done			

		}	// v0 analysis done

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
	hV0_DHasTOF->SetTitle("HasTOF PID; p_{T} (GeV/#it{c}); Efficiency");
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

	// MAKE MC QA
	gSystem->Exec("mkdir qa");
	if (1 && flagMC) {
		gSystem->Exec("mkdir qa/mc");
		TCanvas* cq = new TCanvas("cq","",1000,1000);

		//////// k0s
		hMCV0_PtK0s->SetTitle("; MC pT (GeV/c); Entries");
		hMCV0_PtK0s->SetMarkerColor(2); hMCV0_PtK0s->SetMarkerStyle(20); hMCV0_PtK0s->SetLineColor(2);
		hMCV0_PtK0s->Draw();
		cq->SaveAs("qa/mc/mcPtK0s.png");
		hMCV0_PtK0s->Draw();
		cq->SaveAs("qa/mc/mcPtK0s.png");
		hRCV0_PtK0s->SetTitle("; RC pT (GeV/c); Entries");
		hRCV0_PtK0s->SetMarkerColor(1); hRCV0_PtK0s->SetMarkerStyle(20); hRCV0_PtK0s->SetLineColor(1);
		hRCV0_PtK0s->Draw();
		cq->SaveAs("qa/mc/rcPtK0s.png");

		TH2D* hPE = hMCV0_PtPhiEtaK0s->Project3D("yz");
		hPE->SetTitle(";MC phi;MC eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/mcPEK0s.png");

		hPE = (TH2D*)hRCV0_PtPhiEtaK0s->Project3D("yz");
		hPE->SetTitle(";RC phi;RC eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcPEK0s.png");

		TH1D* h1D = hMCvRC_DCAddK0s->ProjectionY();
		h1D->SetTitle(";DCA daughter-daughter (cm); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDCAddK0s.png");

		h1D = hMCvRC_CPAK0s->ProjectionY();
		h1D->SetTitle(";cosine of pointing angle; Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcCPAK0s.png");

		h1D = hMCvRC_RK0s->ProjectionY();
		h1D->SetTitle(";transverse radius (cm); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcRK0s.png");

		h1D = hRCV0_DPtPhiEtaK0spi->ProjectionX();
		h1D->SetTitle(";RC daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDPtK0spi.png");

		hPE = (TH2D*)hRCV0_DPtPhiEtaK0spi->Project3D("yz");
		hPE->SetTitle(";RC daughter phi;RC daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcDPEK0spi.png");

		hRCV0_DNclusK0spi->SetTitle("; number of TPC clusters; Entries");
		hRCV0_DNclusK0spi->SetMarkerColor(1); hRCV0_DNclusK0spi->SetMarkerStyle(20); hRCV0_DNclusK0spi->SetLineColor(1);
		hRCV0_DNclusK0spi->Draw();
		cq->SaveAs("qa/mc/rcDNclusK0spi.png");

		hRCV0_DIPxyK0spi->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hRCV0_DIPxyK0spi->SetMarkerColor(1); hRCV0_DIPxyK0spi->SetMarkerStyle(20); hRCV0_DIPxyK0spi->SetLineColor(1);
		hRCV0_DIPxyK0spi->Draw();
		cq->SaveAs("qa/mc/rcDIPxyK0spi.png");

		hRCV0_DIPzK0spi->SetTitle("; DCA_{z} daughter-PV; Entries");
		hRCV0_DIPzK0spi->SetMarkerColor(1); hRCV0_DIPzK0spi->SetMarkerStyle(20); hRCV0_DIPzK0spi->SetLineColor(1);
		hRCV0_DIPzK0spi->Draw();
		cq->SaveAs("qa/mc/rcDIPzK0spi.png");

		h1D = hRCV0_DPtPhiEtaK0spi->ProjectionX();
		hRCV0_DHasTPCK0spi->Divide(h1D);
		hRCV0_DHasTPCK0spi->SetTitle("; RC daughter pT (GeV/c); HasTPC efficiency");
		hRCV0_DHasTPCK0spi->SetMarkerColor(1); hRCV0_DHasTPCK0spi->SetMarkerStyle(20); hRCV0_DHasTPCK0spi->SetLineColor(1);
		hRCV0_DHasTPCK0spi->Draw();
		cq->SaveAs("qa/mc/rcDHasTPCK0spi.png");

		hRCV0_DDedxvpK0s->SetTitle(";RC daughter p (GeV/c); <dE/dx> (a.u.)");
		hRCV0_DDedxvpK0s->Draw("colz");
		cq->SaveAs("qa/mc/rcDDedxvpK0s.png");

		cq->SetLogy(1);
		hRCV0_DTpcK0spi->SetTitle("; n#sigma^{TPC}_{#pi}; Entries");
		hRCV0_DTpcK0spi->SetMarkerColor(1); hRCV0_DTpcK0spi->SetMarkerStyle(20); hRCV0_DTpcK0spi->SetLineColor(1);
		hRCV0_DTpcK0spi->Draw();
		cq->SaveAs("qa/mc/rcDTpcK0spi.png");
		cq->SetLogy(0);

		hRCV0_DTPCvphiK0spi->SetTitle("; RC daughter phi; n#sigma^{TPC}_{#pi}");
		hRCV0_DTPCvphiK0spi->Draw("colz");
		cq->SaveAs("qa/mc/rcDTPCvphiK0spi.png");

		hRCV0_DTPCvetaK0spi->SetTitle("; RC daughter eta; n#sigma^{TPC}_{#pi}");
		hRCV0_DTPCvetaK0spi->Draw("colz");
		cq->SaveAs("qa/mc/rcDTPCvetaK0spi.png");

		hRCV0_DNDedxK0spi->SetTitle("; number of dEdx points; Entries");
		hRCV0_DNDedxK0spi->SetMarkerColor(1); hRCV0_DNDedxK0spi->SetMarkerStyle(20); hRCV0_DNDedxK0spi->SetLineColor(1);
		hRCV0_DNDedxK0spi->Draw();
		cq->SaveAs("qa/mc/rcDNDedxK0spi.png");

		h1D = hRCV0_DPtPhiEtaK0spi->ProjectionX();
		hRCV0_DHasTOFK0spi->Divide(h1D);
		hRCV0_DHasTOFK0spi->SetTitle("; RC daughter pT (GeV/c); HasTOF efficiency");
		hRCV0_DHasTOFK0spi->SetMarkerColor(1); hRCV0_DHasTOFK0spi->SetMarkerStyle(20); hRCV0_DHasTOFK0spi->SetLineColor(1);
		hRCV0_DHasTOFK0spi->Draw();
		cq->SaveAs("qa/mc/rcDHasTOFK0spi.png");

		hRCV0_DTofBetavpK0s->SetTitle(";RC daughter p (GeV/c); TOF #beta");
		hRCV0_DTofBetavpK0s->Draw("colz");
		cq->SaveAs("qa/mc/rcDTofBetavpK0s.png");

		cq->SetLogy(1);
		hRCV0_DTofK0spi->SetTitle("; n#sigma^{TOF}_{#pi}; Entries");
		hRCV0_DTofK0spi->SetMarkerColor(1); hRCV0_DTofK0spi->SetMarkerStyle(20); hRCV0_DTofK0spi->SetLineColor(1);
		hRCV0_DTofK0spi->Draw();
		cq->SaveAs("qa/mc/rcDTofK0spi.png");
		cq->SetLogy(0);

		hMCvRC_DTofTimeK0spi->SetTitle(";MC TOF time; RC TOF time");
		hMCvRC_DTofTimeK0spi->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofTimeK0spi.png");

		hMCvRC_DTofLengthK0spi->SetTitle(";MC TOF length; RC TOF length");
		hMCvRC_DTofLengthK0spi->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofLengthK0spi.png");
		/////////////////////////

		//////// l
		hMCV0_PtL->SetTitle("; MC pT (GeV/c); Entries");
		hMCV0_PtL->SetMarkerColor(2); hMCV0_PtL->SetMarkerStyle(20); hMCV0_PtL->SetLineColor(2);
		hMCV0_PtL->Draw();
		cq->SaveAs("qa/mc/mcPtL.png");
		hRCV0_PtL->SetTitle("; RC pT (GeV/c); Entries");
		hRCV0_PtL->SetMarkerColor(1); hRCV0_PtL->SetMarkerStyle(20); hRCV0_PtL->SetLineColor(1);
		hRCV0_PtL->Draw();
		cq->SaveAs("qa/mc/rcPtL.png");

		hPE = (TH2D*) hMCV0_PtPhiEtaL->Project3D("yz");
		hPE->SetTitle(";MC phi;MC eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/mcPEL.png");

		hPE = (TH2D*)hRCV0_PtPhiEtaL->Project3D("yz");
		hPE->SetTitle(";RC phi;RC eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcPEL.png");

		h1D = hMCvRC_DCAddL->ProjectionY();
		h1D->SetTitle(";DCA daughter-daughter (cm); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDCAddL.png");

		h1D = hMCvRC_CPAL->ProjectionY();
		h1D->SetTitle(";cosine of pointing angle; Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcCPAL.png");

		h1D = hMCvRC_RL->ProjectionY();
		h1D->SetTitle(";transverse radius (cm); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcRL.png");

		h1D = hRCV0_DPtPhiEtaLpi->ProjectionX();
		h1D->SetTitle(";RC daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDPtLpi.png");

		h1D = hRCV0_DPtPhiEtaLpr->ProjectionX();
		h1D->SetTitle(";RC daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDPtLpr.png");

		hPE = (TH2D*)hRCV0_DPtPhiEtaLpi->Project3D("yz");
		hPE->SetTitle(";RC daughter phi;RC daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcDPELpi.png");

		hPE = (TH2D*)hRCV0_DPtPhiEtaLpr->Project3D("yz");
		hPE->SetTitle(";RC daughter phi;RC daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcDPELpr.png");

		hRCV0_DNclusLpi->SetTitle("; number of TPC clusters; Entries");
		hRCV0_DNclusLpi->SetMarkerColor(1); hRCV0_DNclusLpi->SetMarkerStyle(20); hRCV0_DNclusLpi->SetLineColor(1);
		hRCV0_DNclusLpi->Draw();
		cq->SaveAs("qa/mc/rcDNclusLpi.png");

		hRCV0_DNclusLpr->SetTitle("; number of TPC clusters; Entries");
		hRCV0_DNclusLpr->SetMarkerColor(1); hRCV0_DNclusLpr->SetMarkerStyle(20); hRCV0_DNclusLpr->SetLineColor(1);
		hRCV0_DNclusLpr->Draw();
		cq->SaveAs("qa/mc/rcDNclusLpr.png");

		hRCV0_DIPxyLpi->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hRCV0_DIPxyLpi->SetMarkerColor(1); hRCV0_DIPxyLpi->SetMarkerStyle(20); hRCV0_DIPxyLpi->SetLineColor(1);
		hRCV0_DIPxyLpi->Draw();
		cq->SaveAs("qa/mc/rcDIPxyLpi.png");

		hRCV0_DIPxyLpr->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hRCV0_DIPxyLpr->SetMarkerColor(1); hRCV0_DIPxyLpr->SetMarkerStyle(20); hRCV0_DIPxyLpr->SetLineColor(1);
		hRCV0_DIPxyLpr->Draw();
		cq->SaveAs("qa/mc/rcDIPxyLpr.png");

		hRCV0_DIPzLpi->SetTitle("; DCA_{z} daughter-PV; Entries");
		hRCV0_DIPzLpi->SetMarkerColor(1); hRCV0_DIPzLpi->SetMarkerStyle(20); hRCV0_DIPzLpi->SetLineColor(1);
		hRCV0_DIPzLpi->Draw();
		cq->SaveAs("qa/mc/rcDIPzLpi.png");

		hRCV0_DIPzLpr->SetTitle("; DCA_{z} daughter-PV; Entries");
		hRCV0_DIPzLpr->SetMarkerColor(1); hRCV0_DIPzLpr->SetMarkerStyle(20); hRCV0_DIPzLpr->SetLineColor(1);
		hRCV0_DIPzLpr->Draw();
		cq->SaveAs("qa/mc/rcDIPzLpr.png");

		h1D = hRCV0_DPtPhiEtaLpi->ProjectionX();
		hRCV0_DHasTPCLpi->Divide(h1D);
		hRCV0_DHasTPCLpi->SetTitle("; RC daughter pT (GeV/c); HasTPC efficiency");
		hRCV0_DHasTPCLpi->SetMarkerColor(1); hRCV0_DHasTPCLpi->SetMarkerStyle(20); hRCV0_DHasTPCLpi->SetLineColor(1);
		hRCV0_DHasTPCLpi->Draw();
		cq->SaveAs("qa/mc/rcDHasTPCLpi.png");

		h1D = hRCV0_DPtPhiEtaLpr->ProjectionX();
		hRCV0_DHasTPCLpr->Divide(h1D);
		hRCV0_DHasTPCLpr->SetTitle("; RC daughter pT (GeV/c); HasTPC efficiency");
		hRCV0_DHasTPCLpr->SetMarkerColor(1); hRCV0_DHasTPCLpr->SetMarkerStyle(20); hRCV0_DHasTPCLpr->SetLineColor(1);
		hRCV0_DHasTPCLpr->Draw();
		cq->SaveAs("qa/mc/rcDHasTPCLpr.png");

		hRCV0_DDedxvpL->SetTitle(";RC daughter p (GeV/c); <dE/dx> (a.u.)");
		hRCV0_DDedxvpL->Draw("colz");
		cq->SaveAs("qa/mc/rcDDedxvpL.png");

		cq->SetLogy(1);
		hRCV0_DTpcLpi->SetTitle("; n#sigma^{TPC}_{#pi}; Entries");
		hRCV0_DTpcLpi->SetMarkerColor(1); hRCV0_DTpcLpi->SetMarkerStyle(20); hRCV0_DTpcLpi->SetLineColor(1);
		hRCV0_DTpcLpi->Draw();
		cq->SaveAs("qa/mc/rcDTpcLpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hRCV0_DTpcLpr->SetTitle("; n#sigma^{TPC}_{p}; Entries");
		hRCV0_DTpcLpr->SetMarkerColor(1); hRCV0_DTpcLpr->SetMarkerStyle(20); hRCV0_DTpcLpr->SetLineColor(1);
		hRCV0_DTpcLpr->Draw();
		cq->SaveAs("qa/mc/rcDTpcLpr.png");
		cq->SetLogy(0);

		hRCV0_DTPCvphiLpr->SetTitle("; RC daughter phi; n#sigma^{TPC}_{p}");
		hRCV0_DTPCvphiLpr->Draw("colz");
		cq->SaveAs("qa/mc/rcDTPCvphiLpr.png");

		hRCV0_DTPCvetaLpr->SetTitle("; RC daughter eta; n#sigma^{TPC}_{p}");
		hRCV0_DTPCvetaLpr->Draw("colz");
		cq->SaveAs("qa/mc/rcDTPCvetaLpr.png");

		hRCV0_DNDedxLpi->SetTitle("; number of dEdx points; Entries");
		hRCV0_DNDedxLpi->SetMarkerColor(1); hRCV0_DNDedxLpi->SetMarkerStyle(20); hRCV0_DNDedxLpi->SetLineColor(1);
		hRCV0_DNDedxLpi->Draw();
		cq->SaveAs("qa/mc/rcDNDedxLpi.png");

		hRCV0_DNDedxLpr->SetTitle("; number of dEdx points; Entries");
		hRCV0_DNDedxLpr->SetMarkerColor(1); hRCV0_DNDedxLpr->SetMarkerStyle(20); hRCV0_DNDedxLpr->SetLineColor(1);
		hRCV0_DNDedxLpr->Draw();
		cq->SaveAs("qa/mc/rcDNDedxLpr.png");

		h1D = hRCV0_DPtPhiEtaLpi->ProjectionX();
		hRCV0_DHasTOFLpi->Divide(h1D);
		hRCV0_DHasTOFLpi->SetTitle("; RC daughter pT (GeV/c); HasTOF efficiency");
		hRCV0_DHasTOFLpi->SetMarkerColor(1); hRCV0_DHasTOFLpi->SetMarkerStyle(20); hRCV0_DHasTOFLpi->SetLineColor(1);
		hRCV0_DHasTOFLpi->Draw();
		cq->SaveAs("qa/mc/rcDHasTOFLpi.png");

		h1D = hRCV0_DPtPhiEtaLpr->ProjectionX();
		hRCV0_DHasTOFLpr->Divide(h1D);
		hRCV0_DHasTOFLpr->SetTitle("; RC daughter pT (GeV/c); HasTOF efficiency");
		hRCV0_DHasTOFLpr->SetMarkerColor(1); hRCV0_DHasTOFLpr->SetMarkerStyle(20); hRCV0_DHasTOFLpr->SetLineColor(1);
		hRCV0_DHasTOFLpr->Draw();
		cq->SaveAs("qa/mc/rcDHasTOFLpr.png");

		hRCV0_DTofBetavpL->SetTitle(";RC daughter p (GeV/c); TOF #beta");
		hRCV0_DTofBetavpL->Draw("colz");
		cq->SaveAs("qa/mc/rcDTofBetavpL.png");

		cq->SetLogy(1);
		hRCV0_DTofLpi->SetTitle("; n#sigma^{TOF}_{#pi}; Entries");
		hRCV0_DTofLpi->SetMarkerColor(1); hRCV0_DTofLpi->SetMarkerStyle(20); hRCV0_DTofLpi->SetLineColor(1);
		hRCV0_DTofLpi->Draw();
		cq->SaveAs("qa/mc/rcDTofLpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hRCV0_DTofLpr->SetTitle("; n#sigma^{TOF}_{#pr}; Entries");
		hRCV0_DTofLpr->SetMarkerColor(1); hRCV0_DTofLpr->SetMarkerStyle(20); hRCV0_DTofLpr->SetLineColor(1);
		hRCV0_DTofLpr->Draw();
		cq->SaveAs("qa/mc/rcDTofLpr.png");
		cq->SetLogy(0);

		hMCvRC_DTofTimeLpi->SetTitle(";MC TOF time; RC TOF time");
		hMCvRC_DTofTimeLpi->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofTimeLpi.png");

		hMCvRC_DTofTimeLpr->SetTitle(";MC TOF time; RC TOF time");
		hMCvRC_DTofTimeLpr->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofTimeLpr.png");

		hMCvRC_DTofLengthLpi->SetTitle(";MC TOF length; RC TOF length");
		hMCvRC_DTofLengthLpi->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofLengthLpi.png");

		hMCvRC_DTofLengthLpr->SetTitle(";MC TOF length; RC TOF length");
		hMCvRC_DTofLengthLpr->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofLengthLpr.png");
		/////////////////////////

		//////// al
		hMCV0_PtAL->SetTitle("; MC pT (GeV/c); Entries");
		hMCV0_PtAL->SetMarkerColor(2); hMCV0_PtAL->SetMarkerStyle(20); hMCV0_PtAL->SetLineColor(2);
		hMCV0_PtAL->Draw();
		cq->SaveAs("qa/mc/mcPtAL.png");
		hRCV0_PtAL->SetTitle("; RC pT (GeV/c); Entries");
		hRCV0_PtAL->SetMarkerColor(1); hRCV0_PtAL->SetMarkerStyle(20); hRCV0_PtAL->SetLineColor(1);
		hRCV0_PtAL->Draw();
		cq->SaveAs("qa/mc/rcPtAL.png");

		hPE = (TH2D*) hMCV0_PtPhiEtaAL->Project3D("yz");
		hPE->SetTitle(";MC phi;MC eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/mcPEAL.png");

		hPE = (TH2D*)hRCV0_PtPhiEtaAL->Project3D("yz");
		hPE->SetTitle(";RC phi;RC eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcPEAL.png");

		h1D = hMCvRC_DCAddAL->ProjectionY();
		h1D->SetTitle(";DCA daughter-daughter (cm); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDCAddAL.png");

		h1D = hMCvRC_CPAAL->ProjectionY();
		h1D->SetTitle(";cosine of pointing angle; Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcCPAAL.png");

		h1D = hMCvRC_RAL->ProjectionY();
		h1D->SetTitle(";transverse radius (cm); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcRAL.png");

		h1D = hRCV0_DPtPhiEtaALpi->ProjectionX();
		h1D->SetTitle(";RC daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDPtALpi.png");

		h1D = hRCV0_DPtPhiEtaALpr->ProjectionX();
		h1D->SetTitle(";RC daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(1); h1D->SetMarkerStyle(20); h1D->SetLineColor(1);
		h1D->Draw();
		cq->SaveAs("qa/mc/rcDPtALpr.png");

		hPE = (TH2D*)hRCV0_DPtPhiEtaALpi->Project3D("yz");
		hPE->SetTitle(";RC daughter phi;RC daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcDPEALpi.png");

		hPE = (TH2D*)hRCV0_DPtPhiEtaALpr->Project3D("yz");
		hPE->SetTitle(";RC daughter phi;RC daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/mc/rcDPEALpr.png");

		hRCV0_DNclusALpi->SetTitle("; number of TPC clusters; Entries");
		hRCV0_DNclusALpi->SetMarkerColor(1); hRCV0_DNclusALpi->SetMarkerStyle(20); hRCV0_DNclusALpi->SetLineColor(1);
		hRCV0_DNclusALpi->Draw();
		cq->SaveAs("qa/mc/rcDNclusALpi.png");

		hRCV0_DNclusALpr->SetTitle("; number of TPC clusters; Entries");
		hRCV0_DNclusALpr->SetMarkerColor(1); hRCV0_DNclusALpr->SetMarkerStyle(20); hRCV0_DNclusALpr->SetLineColor(1);
		hRCV0_DNclusALpr->Draw();
		cq->SaveAs("qa/mc/rcDNclusALpr.png");

		hRCV0_DIPxyALpi->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hRCV0_DIPxyALpi->SetMarkerColor(1); hRCV0_DIPxyALpi->SetMarkerStyle(20); hRCV0_DIPxyALpi->SetLineColor(1);
		hRCV0_DIPxyALpi->Draw();
		cq->SaveAs("qa/mc/rcDIPxyALpi.png");

		hRCV0_DIPxyALpr->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hRCV0_DIPxyALpr->SetMarkerColor(1); hRCV0_DIPxyALpr->SetMarkerStyle(20); hRCV0_DIPxyALpr->SetLineColor(1);
		hRCV0_DIPxyALpr->Draw();
		cq->SaveAs("qa/mc/rcDIPxyALpr.png");

		hRCV0_DIPzALpi->SetTitle("; DCA_{z} daughter-PV; Entries");
		hRCV0_DIPzALpi->SetMarkerColor(1); hRCV0_DIPzALpi->SetMarkerStyle(20); hRCV0_DIPzALpi->SetLineColor(1);
		hRCV0_DIPzALpi->Draw();
		cq->SaveAs("qa/mc/rcDIPzALpi.png");

		hRCV0_DIPzALpr->SetTitle("; DCA_{z} daughter-PV; Entries");
		hRCV0_DIPzALpr->SetMarkerColor(1); hRCV0_DIPzALpr->SetMarkerStyle(20); hRCV0_DIPzALpr->SetLineColor(1);
		hRCV0_DIPzALpr->Draw();
		cq->SaveAs("qa/mc/rcDIPzALpr.png");

		h1D = hRCV0_DPtPhiEtaALpi->ProjectionX();
		hRCV0_DHasTPCALpi->Divide(h1D);
		hRCV0_DHasTPCALpi->SetTitle("; RC daughter pT (GeV/c); HasTPC efficiency");
		hRCV0_DHasTPCALpi->SetMarkerColor(1); hRCV0_DHasTPCALpi->SetMarkerStyle(20); hRCV0_DHasTPCALpi->SetLineColor(1);
		hRCV0_DHasTPCALpi->Draw();
		cq->SaveAs("qa/mc/rcDHasTPCALpi.png");

		h1D = hRCV0_DPtPhiEtaALpr->ProjectionX();
		hRCV0_DHasTPCALpr->Divide(h1D);
		hRCV0_DHasTPCALpr->SetTitle("; RC daughter pT (GeV/c); HasTPC efficiency");
		hRCV0_DHasTPCALpr->SetMarkerColor(1); hRCV0_DHasTPCALpr->SetMarkerStyle(20); hRCV0_DHasTPCALpr->SetLineColor(1);
		hRCV0_DHasTPCALpr->Draw();
		cq->SaveAs("qa/mc/rcDHasTPCALpr.png");

		hRCV0_DDedxvpAL->SetTitle(";RC daughter p (GeV/c); <dE/dx> (a.u.)");
		hRCV0_DDedxvpAL->Draw("colz");
		cq->SaveAs("qa/mc/rcDDedxvpAL.png");

		cq->SetLogy(1);
		hRCV0_DTpcALpi->SetTitle("; n#sigma^{TPC}_{#pi}; Entries");
		hRCV0_DTpcALpi->SetMarkerColor(1); hRCV0_DTpcALpi->SetMarkerStyle(20); hRCV0_DTpcALpi->SetLineColor(1);
		hRCV0_DTpcALpi->Draw();
		cq->SaveAs("qa/mc/rcDTpcALpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hRCV0_DTpcALpr->SetTitle("; n#sigma^{TPC}_{p}; Entries");
		hRCV0_DTpcALpr->SetMarkerColor(1); hRCV0_DTpcALpr->SetMarkerStyle(20); hRCV0_DTpcALpr->SetLineColor(1);
		hRCV0_DTpcALpr->Draw();
		cq->SaveAs("qa/mc/rcDTpcALpr.png");
		cq->SetLogy(0);

		hRCV0_DTPCvphiALpr->SetTitle("; RC daughter phi; n#sigma^{TPC}_{p}");
		hRCV0_DTPCvphiALpr->Draw("colz");
		cq->SaveAs("qa/mc/rcDTPCvphiALpr.png");

		hRCV0_DTPCvetaALpr->SetTitle("; RC daughter eta; n#sigma^{TPC}_{p}");
		hRCV0_DTPCvetaALpr->Draw("colz");
		cq->SaveAs("qa/mc/rcDTPCvetaALpr.png");

		hRCV0_DNDedxALpi->SetTitle("; number of dEdx points; Entries");
		hRCV0_DNDedxALpi->SetMarkerColor(1); hRCV0_DNDedxALpi->SetMarkerStyle(20); hRCV0_DNDedxALpi->SetLineColor(1);
		hRCV0_DNDedxALpi->Draw();
		cq->SaveAs("qa/mc/rcDNDedxALpi.png");

		hRCV0_DNDedxALpr->SetTitle("; number of dEdx points; Entries");
		hRCV0_DNDedxALpr->SetMarkerColor(1); hRCV0_DNDedxALpr->SetMarkerStyle(20); hRCV0_DNDedxALpr->SetLineColor(1);
		hRCV0_DNDedxALpr->Draw();
		cq->SaveAs("qa/mc/rcDNDedxALpr.png");

		h1D = hRCV0_DPtPhiEtaALpi->ProjectionX();
		hRCV0_DHasTOFALpi->Divide(h1D);
		hRCV0_DHasTOFALpi->SetTitle("; RC daughter pT (GeV/c); HasTOF efficiency");
		hRCV0_DHasTOFALpi->SetMarkerColor(1); hRCV0_DHasTOFALpi->SetMarkerStyle(20); hRCV0_DHasTOFALpi->SetLineColor(1);
		hRCV0_DHasTOFALpi->Draw();
		cq->SaveAs("qa/mc/rcDHasTOFALpi.png");

		h1D = hRCV0_DPtPhiEtaALpr->ProjectionX();
		hRCV0_DHasTOFALpr->Divide(h1D);
		hRCV0_DHasTOFALpr->SetTitle("; RC daughter pT (GeV/c); HasTOF efficiency");
		hRCV0_DHasTOFALpr->SetMarkerColor(1); hRCV0_DHasTOFALpr->SetMarkerStyle(20); hRCV0_DHasTOFALpr->SetLineColor(1);
		hRCV0_DHasTOFALpr->Draw();
		cq->SaveAs("qa/mc/rcDHasTOFALpr.png");

		hRCV0_DTofBetavpAL->SetTitle(";RC daughter p (GeV/c); TOF #beta");
		hRCV0_DTofBetavpAL->Draw("colz");
		cq->SaveAs("qa/mc/rcDTofBetavpAL.png");

		cq->SetLogy(1);
		hRCV0_DTofALpi->SetTitle("; n#sigma^{TOF}_{#pi}; Entries");
		hRCV0_DTofALpi->SetMarkerColor(1); hRCV0_DTofALpi->SetMarkerStyle(20); hRCV0_DTofALpi->SetLineColor(1);
		hRCV0_DTofALpi->Draw();
		cq->SaveAs("qa/mc/rcDTofALpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hRCV0_DTofALpr->SetTitle("; n#sigma^{TOF}_{#pr}; Entries");
		hRCV0_DTofALpr->SetMarkerColor(1); hRCV0_DTofALpr->SetMarkerStyle(20); hRCV0_DTofALpr->SetLineColor(1);
		hRCV0_DTofALpr->Draw();
		cq->SaveAs("qa/mc/rcDTofALpr.png");
		cq->SetLogy(0);

		hMCvRC_DTofTimeALpi->SetTitle(";MC TOF time; RC TOF time");
		hMCvRC_DTofTimeALpi->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofTimeALpi.png");

		hMCvRC_DTofTimeALpr->SetTitle(";MC TOF time; RC TOF time");
		hMCvRC_DTofTimeALpr->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofTimeALpr.png");

		hMCvRC_DTofLengthALpi->SetTitle(";MC TOF length; RC TOF length");
		hMCvRC_DTofLengthALpi->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofLengthALpi.png");

		hMCvRC_DTofLengthALpr->SetTitle(";MC TOF length; RC TOF length");
		hMCvRC_DTofLengthALpr->Draw("colz");
		cq->SaveAs("qa/mc/mcrcDTofLengthALpr.png");
		/////////////////////////

		cq->SetLogz(1);
		hMCvRC_PtK0s->SetTitle("; MC pT (GeV/c); RC pT (GeV/c)");
		hMCvRC_PtK0s->Draw("colz");
		cq->SaveAs("qa/mc/mcrcPtK0s.png");
		hMCvRC_PtL->SetTitle("; MC pT (GeV/c); RC pT (GeV/c)");
		hMCvRC_PtL->Draw("colz");
		cq->SaveAs("qa/mc/mcrcPtL.png");
		hMCvRC_PtAL->SetTitle("; MC pT (GeV/c); RC pT (GeV/c)");
		hMCvRC_PtAL->Draw("colz");
		cq->SaveAs("qa/mc/mcrcPtAL.png");

		hMCvRC_PhiK0s->SetTitle("; MC Phi ; RC Phi ");
		hMCvRC_PhiK0s->Draw("colz");
		cq->SaveAs("qa/mc/mcrcPhiK0s.png");
		hMCvRC_PhiL->SetTitle("; MC Phi ; RC Phi ");
		hMCvRC_PhiL->Draw("colz");
		cq->SaveAs("qa/mc/mcrcPhiL.png");
		hMCvRC_PhiAL->SetTitle("; MC Phi ; RC Phi ");
		hMCvRC_PhiAL->Draw("colz");
		cq->SaveAs("qa/mc/mcrcPhiAL.png");
		
		hMCvRC_EtaK0s->SetTitle("; MC Eta ; RC Eta ");
		hMCvRC_EtaK0s->Draw("colz");
		cq->SaveAs("qa/mc/mcrcEtaK0s.png");
		hMCvRC_EtaL->SetTitle("; MC Eta ; RC Eta ");
		hMCvRC_EtaL->Draw("colz");
		cq->SaveAs("qa/mc/mcrcEtaL.png");
		hMCvRC_EtaAL->SetTitle("; MC Eta ; RC Eta ");
		hMCvRC_EtaAL->Draw("colz");
		cq->SaveAs("qa/mc/mcrcEtaAL.png");


	}



	// DATA QA
	if (1 && !flagMC) {
		gSystem->Exec("mkdir qa/d");
		TCanvas* cq = new TCanvas("cq","",1000,1000);

		//////// k0s
		TH2D* hPE = (TH2D*)hV0_PtPhiEtaK0s->Project3D("yz");
		hPE->SetTitle(";phi;eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dPEK0s.png");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dPEK0s.png");

		TH1D* h1D = (TH1D*)hV0_DCAddK0s;
		h1D->SetTitle(";DCA daughter-daughter (cm); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/DCAddK0s.png");

		h1D = (TH1D*)hV0_CPAK0s;
		h1D->SetTitle(";cosine of pointing angle; Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dCPAK0s.png");

		h1D = (TH1D*)hV0_RK0s;
		h1D->SetTitle(";transverse radius (cm); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dRK0s.png");

		h1D = hV0_DPtPhiEtaK0spi->ProjectionX();
		h1D->SetTitle(";daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDPtK0spi.png");

		hPE = (TH2D*)hV0_DPtPhiEtaK0spi->Project3D("yz");
		hPE->SetTitle(";daughter phi;RC daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dDPEK0spi.png");

		hV0_DNclusK0spi->SetTitle("; number of TPC clusters; Entries");
		hV0_DNclusK0spi->SetMarkerColor(4); hV0_DNclusK0spi->SetMarkerStyle(20); hV0_DNclusK0spi->SetLineColor(4);
		hV0_DNclusK0spi->Draw();
		cq->SaveAs("qa/d/dDNclusK0spi.png");

		hV0_DIPxyK0spi->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hV0_DIPxyK0spi->SetMarkerColor(4); hV0_DIPxyK0spi->SetMarkerStyle(20); hV0_DIPxyK0spi->SetLineColor(4);
		hV0_DIPxyK0spi->Draw();
		cq->SaveAs("qa/d/dDIPxyK0spi.png");

		hV0_DIPzK0spi->SetTitle("; DCA_{z} daughter-PV; Entries");
		hV0_DIPzK0spi->SetMarkerColor(4); hV0_DIPzK0spi->SetMarkerStyle(20); hV0_DIPzK0spi->SetLineColor(4);
		hV0_DIPzK0spi->Draw();
		cq->SaveAs("qa/d/dDIPzK0spi.png");

		h1D = hV0_DPtPhiEtaK0spi->ProjectionX();
		hV0_DHasTPCK0spi->Divide(h1D);
		hV0_DHasTPCK0spi->SetTitle("; daughter pT (GeV/c); HasTPC efficiency");
		hV0_DHasTPCK0spi->SetMarkerColor(4); hV0_DHasTPCK0spi->SetMarkerStyle(20); hV0_DHasTPCK0spi->SetLineColor(4);
		hV0_DHasTPCK0spi->Draw();
		cq->SaveAs("qa/d/dDHasTPCK0spi.png");

		hV0_DDedxvpK0s->SetTitle(";daughter p (GeV/c); <dE/dx> (a.u.)");
		hV0_DDedxvpK0s->Draw("colz");
		cq->SaveAs("qa/d/dDDedxvpK0s.png");

		cq->SetLogy(1);
		hV0_DTpcK0spi->SetTitle("; n#sigma^{TPC}_{#pi}; Entries");
		hV0_DTpcK0spi->SetMarkerColor(1); hV0_DTpcK0spi->SetMarkerStyle(20); hV0_DTpcK0spi->SetLineColor(1);
		hV0_DTpcK0spi->Draw();
		cq->SaveAs("qa/d/dDTpcK0spi.png");
		cq->SetLogy(0);

		hV0_DTPCvphiK0spi->SetTitle("; daughter phi; n#sigma^{TPC}_{#pi}");
		hV0_DTPCvphiK0spi->Draw("colz");
		cq->SaveAs("qa/d/dDTPCvphiK0spi.png");

		hV0_DTPCvetaK0spi->SetTitle("; daughter eta; n#sigma^{TPC}_{#pi}");
		hV0_DTPCvetaK0spi->Draw("colz");
		cq->SaveAs("qa/d/dDTPCvetaK0spi.png");

		hV0_DNDedxK0spi->SetTitle("; number of dEdx points; Entries");
		hV0_DNDedxK0spi->SetMarkerColor(4); hV0_DNDedxK0spi->SetMarkerStyle(20); hV0_DNDedxK0spi->SetLineColor(4);
		hV0_DNDedxK0spi->Draw();
		cq->SaveAs("qa/d/dDNDedxK0spi.png");

		h1D = hV0_DPtPhiEtaK0spi->ProjectionX();
		hV0_DHasTOFK0spi->Divide(h1D);
		hV0_DHasTOFK0spi->SetTitle("; daughter pT (GeV/c); HasTOF efficiency");
		hV0_DHasTOFK0spi->SetMarkerColor(4); hV0_DHasTOFK0spi->SetMarkerStyle(20); hV0_DHasTOFK0spi->SetLineColor(4);
		hV0_DHasTOFK0spi->Draw();
		cq->SaveAs("qa/d/dDHasTOFK0spi.png");

		hV0_DTofBetavpK0s->SetTitle(";daughter p (GeV/c); TOF #beta");
		hV0_DTofBetavpK0s->Draw("colz");
		cq->SaveAs("qa/d/dDTofBetavpK0s.png");

		cq->SetLogy(1);
		hV0_DTofK0spi->SetTitle("; n#sigma^{TOF}_{#pi}; Entries");
		hV0_DTofK0spi->SetMarkerColor(1); hV0_DTofK0spi->SetMarkerStyle(20); hV0_DTofK0spi->SetLineColor(1);
		hV0_DTofK0spi->Draw();
		cq->SaveAs("qa/d/dDTofK0spi.png");
		cq->SetLogy(0);
		/////////////////////////

		//////// l
		hPE = (TH2D*) hV0_PtPhiEtaL->Project3D("yz");
		hPE->SetTitle("; phi;eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dPEL.png");

		h1D = (TH1D*)hV0_DCAddL;
		h1D->SetTitle(";DCA daughter-daughter (cm); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDCAddL.png");

		h1D = (TH1D*)hV0_CPAL;
		h1D->SetTitle(";cosine of pointing angle; Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dCPAL.png");

		h1D = (TH1D*)hV0_RL;
		h1D->SetTitle(";transverse radius (cm); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dRL.png");

		h1D = (TH1D*)hV0_DPtPhiEtaLpi->ProjectionX();
		h1D->SetTitle(";daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDPtLpi.png");

		h1D = hV0_DPtPhiEtaLpr->ProjectionX();
		h1D->SetTitle(";daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDPtLpr.png");

		hPE = (TH2D*)hV0_DPtPhiEtaLpi->Project3D("yz");
		hPE->SetTitle(";daughter phi;daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dDPELpi.png");

		hPE = (TH2D*)hV0_DPtPhiEtaLpr->Project3D("yz");
		hPE->SetTitle(";daughter phi;daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dDPELpr.png");

		hV0_DNclusLpi->SetTitle("; number of TPC clusters; Entries");
		hV0_DNclusLpi->SetMarkerColor(4); hV0_DNclusLpi->SetMarkerStyle(20); hV0_DNclusLpi->SetLineColor(4);
		hV0_DNclusLpi->Draw();
		cq->SaveAs("qa/d/dDNclusLpi.png");

		hV0_DNclusLpr->SetTitle("; number of TPC clusters; Entries");
		hV0_DNclusLpr->SetMarkerColor(4); hV0_DNclusLpr->SetMarkerStyle(20); hV0_DNclusLpr->SetLineColor(4);
		hV0_DNclusLpr->Draw();
		cq->SaveAs("qa/d/dDNclusLpr.png");

		hV0_DIPxyLpi->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hV0_DIPxyLpi->SetMarkerColor(4); hV0_DIPxyLpi->SetMarkerStyle(20); hV0_DIPxyLpi->SetLineColor(4);
		hV0_DIPxyLpi->Draw();
		cq->SaveAs("qa/d/dDIPxyLpi.png");

		hV0_DIPxyLpr->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hV0_DIPxyLpr->SetMarkerColor(4); hV0_DIPxyLpr->SetMarkerStyle(20); hV0_DIPxyLpr->SetLineColor(4);
		hV0_DIPxyLpr->Draw();
		cq->SaveAs("qa/d/dDIPxyLpr.png");

		hV0_DIPzLpi->SetTitle("; DCA_{z} daughter-PV; Entries");
		hV0_DIPzLpi->SetMarkerColor(4); hV0_DIPzLpi->SetMarkerStyle(20); hV0_DIPzLpi->SetLineColor(4);
		hV0_DIPzLpi->Draw();
		cq->SaveAs("qa/d/dDIPzLpi.png");

		hV0_DIPzLpr->SetTitle("; DCA_{z} daughter-PV; Entries");
		hV0_DIPzLpr->SetMarkerColor(4); hV0_DIPzLpr->SetMarkerStyle(20); hV0_DIPzLpr->SetLineColor(4);
		hV0_DIPzLpr->Draw();
		cq->SaveAs("qa/d/dDIPzLpr.png");

		h1D = hV0_DPtPhiEtaLpi->ProjectionX();
		hV0_DHasTPCLpi->Divide(h1D);
		hV0_DHasTPCLpi->SetTitle(";  daughter pT (GeV/c); HasTPC efficiency");
		hV0_DHasTPCLpi->SetMarkerColor(4); hV0_DHasTPCLpi->SetMarkerStyle(20); hV0_DHasTPCLpi->SetLineColor(4);
		hV0_DHasTPCLpi->Draw();
		cq->SaveAs("qa/d/dDHasTPCLpi.png");

		h1D = hV0_DPtPhiEtaLpr->ProjectionX();
		hV0_DHasTPCLpr->Divide(h1D);
		hV0_DHasTPCLpr->SetTitle("; daughter pT (GeV/c); HasTPC efficiency");
		hV0_DHasTPCLpr->SetMarkerColor(4); hV0_DHasTPCLpr->SetMarkerStyle(20); hV0_DHasTPCLpr->SetLineColor(4);
		hV0_DHasTPCLpr->Draw();
		cq->SaveAs("qa/d/dDHasTPCLpr.png");

		hV0_DDedxvpL->SetTitle(";daughter p (GeV/c); <dE/dx> (a.u.)");
		hV0_DDedxvpL->Draw("colz");
		cq->SaveAs("qa/d/dDDedxvpL.png");

		cq->SetLogy(1);
		hV0_DTpcLpi->SetTitle("; n#sigma^{TPC}_{#pi}; Entries");
		hV0_DTpcLpi->SetMarkerColor(4); hV0_DTpcLpi->SetMarkerStyle(20); hV0_DTpcLpi->SetLineColor(4);
		hV0_DTpcLpi->Draw();
		cq->SaveAs("qa/d/dDTpcLpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hV0_DTpcLpr->SetTitle("; n#sigma^{TPC}_{p}; Entries");
		hV0_DTpcLpr->SetMarkerColor(4); hV0_DTpcLpr->SetMarkerStyle(20); hV0_DTpcLpr->SetLineColor(4);
		hV0_DTpcLpr->Draw();
		cq->SaveAs("qa/d/dDTpcLpr.png");
		cq->SetLogy(0);

		hV0_DTPCvphiLpr->SetTitle("; daughter phi; n#sigma^{TPC}_{p}");
		hV0_DTPCvphiLpr->Draw("colz");
		cq->SaveAs("qa/d/dDTPCvphiLpr.png");

		hV0_DTPCvetaLpr->SetTitle("; daughter eta; n#sigma^{TPC}_{p}");
		hV0_DTPCvetaLpr->Draw("colz");
		cq->SaveAs("qa/d/dDTPCvetaLpr.png");

		hV0_DNDedxLpi->SetTitle("; number of dEdx points; Entries");
		hV0_DNDedxLpi->SetMarkerColor(4); hV0_DNDedxLpi->SetMarkerStyle(20); hV0_DNDedxLpi->SetLineColor(4);
		hV0_DNDedxLpi->Draw();
		cq->SaveAs("qa/d/dDNDedxLpi.png");

		hV0_DNDedxLpr->SetTitle("; number of dEdx points; Entries");
		hV0_DNDedxLpr->SetMarkerColor(4); hV0_DNDedxLpr->SetMarkerStyle(20); hV0_DNDedxLpr->SetLineColor(4);
		hV0_DNDedxLpr->Draw();
		cq->SaveAs("qa/d/dDNDedxLpr.png");

		h1D = hV0_DPtPhiEtaLpi->ProjectionX();
		hV0_DHasTOFLpi->Divide(h1D);
		hV0_DHasTOFLpi->SetTitle("; daughter pT (GeV/c); HasTOF efficiency");
		hV0_DHasTOFLpi->SetMarkerColor(4); hV0_DHasTOFLpi->SetMarkerStyle(20); hV0_DHasTOFLpi->SetLineColor(4);
		hV0_DHasTOFLpi->Draw();
		cq->SaveAs("qa/d/dDHasTOFLpi.png");

		h1D = hV0_DPtPhiEtaLpr->ProjectionX();
		hV0_DHasTOFLpr->Divide(h1D);
		hV0_DHasTOFLpr->SetTitle("; daughter pT (GeV/c); HasTOF efficiency");
		hV0_DHasTOFLpr->SetMarkerColor(4); hV0_DHasTOFLpr->SetMarkerStyle(20); hV0_DHasTOFLpr->SetLineColor(4);
		hV0_DHasTOFLpr->Draw();
		cq->SaveAs("qa/d/dDHasTOFLpr.png");

		hV0_DTofBetavpL->SetTitle(";daughter p (GeV/c); TOF #beta");
		hV0_DTofBetavpL->Draw("colz");
		cq->SaveAs("qa/d/dDTofBetavpL.png");

		cq->SetLogy(1);
		hV0_DTofLpi->SetTitle("; n#sigma^{TOF}_{#pi}; Entries");
		hV0_DTofLpi->SetMarkerColor(4); hV0_DTofLpi->SetMarkerStyle(20); hV0_DTofLpi->SetLineColor(4);
		hV0_DTofLpi->Draw();
		cq->SaveAs("qa/d/dDTofLpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hV0_DTofLpr->SetTitle("; n#sigma^{TOF}_{#pr}; Entries");
		hV0_DTofLpr->SetMarkerColor(4); hV0_DTofLpr->SetMarkerStyle(20); hV0_DTofLpr->SetLineColor(4);
		hV0_DTofLpr->Draw();
		cq->SaveAs("qa/d/dDTofLpr.png");
		cq->SetLogy(0);
		/////////////////////////

		//////// al
		hPE = (TH2D*) hV0_PtPhiEtaAL->Project3D("yz");
		hPE->SetTitle("; phi;eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dPEAL.png");

		h1D = (TH1D*)hV0_DCAddAL;
		h1D->SetTitle(";DCA daughter-daughter (cm); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDCAddAL.png");

		h1D = (TH1D*)hV0_CPAAL;
		h1D->SetTitle(";cosine of pointing angle; Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dCPAAL.png");

		h1D = (TH1D*)hV0_RAL;
		h1D->SetTitle(";transverse radius (cm); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dRAL.png");

		h1D = (TH1D*)hV0_DPtPhiEtaALpi->ProjectionX();
		h1D->SetTitle(";daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDPtALpi.png");

		h1D = hV0_DPtPhiEtaALpr->ProjectionX();
		h1D->SetTitle(";daughter pT (GeV/c); Entries");
		h1D->SetMarkerColor(4); h1D->SetMarkerStyle(20); h1D->SetLineColor(4);
		h1D->Draw();
		cq->SaveAs("qa/d/dDPtALpr.png");

		hPE = (TH2D*)hV0_DPtPhiEtaALpi->Project3D("yz");
		hPE->SetTitle(";daughter phi;daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dDPEALpi.png");

		hPE = (TH2D*)hV0_DPtPhiEtaALpr->Project3D("yz");
		hPE->SetTitle(";daughter phi;daughter eta");
		hPE->Draw("colz");
		cq->SaveAs("qa/d/dDPEALpr.png");

		hV0_DNclusALpi->SetTitle("; number of TPC clusters; Entries");
		hV0_DNclusALpi->SetMarkerColor(4); hV0_DNclusALpi->SetMarkerStyle(20); hV0_DNclusALpi->SetLineColor(4);
		hV0_DNclusALpi->Draw();
		cq->SaveAs("qa/d/dDNclusALpi.png");

		hV0_DNclusALpr->SetTitle("; number of TPC clusters; Entries");
		hV0_DNclusALpr->SetMarkerColor(4); hV0_DNclusALpr->SetMarkerStyle(20); hV0_DNclusALpr->SetLineColor(4);
		hV0_DNclusALpr->Draw();
		cq->SaveAs("qa/d/dDNclusALpr.png");

		hV0_DIPxyALpi->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hV0_DIPxyALpi->SetMarkerColor(4); hV0_DIPxyALpi->SetMarkerStyle(20); hV0_DIPxyALpi->SetLineColor(4);
		hV0_DIPxyALpi->Draw();
		cq->SaveAs("qa/d/dDIPxyALpi.png");

		hV0_DIPxyALpr->SetTitle("; DCA_{xy} daughter-PV; Entries");
		hV0_DIPxyALpr->SetMarkerColor(4); hV0_DIPxyALpr->SetMarkerStyle(20); hV0_DIPxyALpr->SetLineColor(4);
		hV0_DIPxyALpr->Draw();
		cq->SaveAs("qa/d/dDIPxyALpr.png");

		hV0_DIPzALpi->SetTitle("; DCA_{z} daughter-PV; Entries");
		hV0_DIPzALpi->SetMarkerColor(4); hV0_DIPzALpi->SetMarkerStyle(20); hV0_DIPzALpi->SetLineColor(4);
		hV0_DIPzALpi->Draw();
		cq->SaveAs("qa/d/dDIPzALpi.png");

		hV0_DIPzALpr->SetTitle("; DCA_{z} daughter-PV; Entries");
		hV0_DIPzALpr->SetMarkerColor(4); hV0_DIPzALpr->SetMarkerStyle(20); hV0_DIPzALpr->SetLineColor(4);
		hV0_DIPzALpr->Draw();
		cq->SaveAs("qa/d/dDIPzALpr.png");

		h1D = hV0_DPtPhiEtaALpi->ProjectionX();
		hV0_DHasTPCALpi->Divide(h1D);
		hV0_DHasTPCALpi->SetTitle(";  daughter pT (GeV/c); HasTPC efficiency");
		hV0_DHasTPCALpi->SetMarkerColor(4); hV0_DHasTPCALpi->SetMarkerStyle(20); hV0_DHasTPCALpi->SetLineColor(4);
		hV0_DHasTPCALpi->Draw();
		cq->SaveAs("qa/d/dDHasTPCALpi.png");

		h1D = hV0_DPtPhiEtaALpr->ProjectionX();
		hV0_DHasTPCALpr->Divide(h1D);
		hV0_DHasTPCALpr->SetTitle("; daughter pT (GeV/c); HasTPC efficiency");
		hV0_DHasTPCALpr->SetMarkerColor(4); hV0_DHasTPCALpr->SetMarkerStyle(20); hV0_DHasTPCALpr->SetLineColor(4);
		hV0_DHasTPCALpr->Draw();
		cq->SaveAs("qa/d/dDHasTPCALpr.png");

		hV0_DDedxvpAL->SetTitle(";daughter p (GeV/c); <dE/dx> (a.u.)");
		hV0_DDedxvpAL->Draw("colz");
		cq->SaveAs("qa/d/dDDedxvpAL.png");

		cq->SetLogy(1);
		hV0_DTpcALpi->SetTitle("; n#sigma^{TPC}_{#pi}; Entries");
		hV0_DTpcALpi->SetMarkerColor(4); hV0_DTpcALpi->SetMarkerStyle(20); hV0_DTpcALpi->SetLineColor(4);
		hV0_DTpcALpi->Draw();
		cq->SaveAs("qa/d/dDTpcALpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hV0_DTpcALpr->SetTitle("; n#sigma^{TPC}_{p}; Entries");
		hV0_DTpcALpr->SetMarkerColor(4); hV0_DTpcALpr->SetMarkerStyle(20); hV0_DTpcALpr->SetLineColor(4);
		hV0_DTpcALpr->Draw();
		cq->SaveAs("qa/d/dDTpcALpr.png");
		cq->SetLogy(0);

		hV0_DTPCvphiALpr->SetTitle("; daughter phi; n#sigma^{TPC}_{p}");
		hV0_DTPCvphiALpr->Draw("colz");
		cq->SaveAs("qa/d/dDTPCvphiALpr.png");

		hV0_DTPCvetaALpr->SetTitle("; daughter eta; n#sigma^{TPC}_{p}");
		hV0_DTPCvetaALpr->Draw("colz");
		cq->SaveAs("qa/d/dDTPCvetaALpr.png");

		hV0_DNDedxALpi->SetTitle("; number of dEdx points; Entries");
		hV0_DNDedxALpi->SetMarkerColor(4); hV0_DNDedxALpi->SetMarkerStyle(20); hV0_DNDedxALpi->SetLineColor(4);
		hV0_DNDedxALpi->Draw();
		cq->SaveAs("qa/d/dDNDedxALpi.png");

		hV0_DNDedxALpr->SetTitle("; number of dEdx points; Entries");
		hV0_DNDedxALpr->SetMarkerColor(4); hV0_DNDedxALpr->SetMarkerStyle(20); hV0_DNDedxALpr->SetLineColor(4);
		hV0_DNDedxALpr->Draw();
		cq->SaveAs("qa/d/dDNDedxALpr.png");

		h1D = hV0_DPtPhiEtaALpi->ProjectionX();
		hV0_DHasTOFALpi->Divide(h1D);
		hV0_DHasTOFALpi->SetTitle("; daughter pT (GeV/c); HasTOF efficiency");
		hV0_DHasTOFALpi->SetMarkerColor(4); hV0_DHasTOFALpi->SetMarkerStyle(20); hV0_DHasTOFALpi->SetLineColor(4);
		hV0_DHasTOFALpi->Draw();
		cq->SaveAs("qa/d/dDHasTOFALpi.png");

		h1D = hV0_DPtPhiEtaALpr->ProjectionX();
		hV0_DHasTOFALpr->Divide(h1D);
		hV0_DHasTOFALpr->SetTitle("; daughter pT (GeV/c); HasTOF efficiency");
		hV0_DHasTOFALpr->SetMarkerColor(4); hV0_DHasTOFALpr->SetMarkerStyle(20); hV0_DHasTOFALpr->SetLineColor(4);
		hV0_DHasTOFALpr->Draw();
		cq->SaveAs("qa/d/dDHasTOFALpr.png");

		hV0_DTofBetavpAL->SetTitle(";daughter p (GeV/c); TOF #beta");
		hV0_DTofBetavpAL->Draw("colz");
		cq->SaveAs("qa/d/dDTofBetavpAL.png");

		cq->SetLogy(1);
		hV0_DTofALpi->SetTitle("; n#sigma^{TOF}_{#pi}; Entries");
		hV0_DTofALpi->SetMarkerColor(4); hV0_DTofALpi->SetMarkerStyle(20); hV0_DTofALpi->SetLineColor(4);
		hV0_DTofALpi->Draw();
		cq->SaveAs("qa/d/dDTofALpi.png");
		cq->SetLogy(0);

		cq->SetLogy(1);
		hV0_DTofALpr->SetTitle("; n#sigma^{TOF}_{#pr}; Entries");
		hV0_DTofALpr->SetMarkerColor(4); hV0_DTofALpr->SetMarkerStyle(20); hV0_DTofALpr->SetLineColor(4);
		hV0_DTofALpr->Draw();
		cq->SaveAs("qa/d/dDTofALpr.png");
		cq->SetLogy(0);
		/////////////////////////
	}



	TCanvas* can1 = new TCanvas("can1","",1000,700);
	hV0_DHasTOF->Draw();
	can1->SaveAs(path+"eff_tof.png");
	hV0_DHasTPC->Draw();
	can1->SaveAs(path+"eff_tpc.png");

	hRCV0_PtK0s->Divide(hMCV0_PtK0s);
	hRCV0_PtK0s->SetTitle(";p_{T} (GeV/#it{c});K_{S}^{0} reconstruction efficiency");
	hRCV0_PtK0s->GetYaxis()->SetRangeUser(0.,0.65);
	hRCV0_PtK0s->SetLineColor(kRed);
	hRCV0_PtK0s->SetMarkerStyle(20);
	hRCV0_PtK0s->Draw();
	can1->SaveAs(path+"effrec_k0s.png");
	hRCV0_PtL->Divide(hMCV0_PtL);
	hRCV0_PtL->SetTitle(";p_{T} (GeV/#it{c});#Lambda reconstruction efficiency");
	hRCV0_PtL->GetYaxis()->SetRangeUser(0.,0.65);
	hRCV0_PtL->SetLineColor(kRed);
	hRCV0_PtL->SetMarkerStyle(20);
	hRCV0_PtL->Draw();
	can1->SaveAs(path+"effrec_l.png");
	hRCV0_PtAL->Divide(hMCV0_PtAL);
	hRCV0_PtAL->SetTitle(";p_{T} (GeV/#it{c});#bar{#Lambda} reconstruction efficiency");
	hRCV0_PtAL->GetYaxis()->SetRangeUser(0.,0.65);
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

	printf(" WHAT IS UP \n" );
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