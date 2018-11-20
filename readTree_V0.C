// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
//#include <AliAnalysisPIDTrack.h>

#include <iostream>

using namespace std;
using namespace RooFit;

TChain* mChain;
AliAnalysisPIDEvent* mEvent;
TClonesArray* bTracks = 0;
TClonesArray* bV0s = 0;
TCanvas* cFits[3];
int canCounter = 0;

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

	// IM photonic e+e- rejection cut?

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

Float_t ExtractYield(TH1D* hist) {	// extracting with RooFit
	Float_t val = hist->Integral(hist->FindBin(-0.04),hist->FindBin(0.04));
	
	Float_t fitMin = -0.04, fitMax = 0.04;
	hist->Rebin(8);
	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
	RooDataHist DT_hist("DT_hist","DT_hist",MassDT,Import(*hist));

	RooRealVar pGaus1A("pGaus1A","Mean 1",-0.001,0.001);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0,0.2);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e06);

	RooRealVar pGaus2A("pGaus2A","Mean 2",-0.001,0.001);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0,0.2);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus2A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",1,0,1e06);

	RooRealVar pPolBgA("pPolBgA","Pol. par. A",-1,-100,100);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);//RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",1,0,1e06);

	//RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
	//fTotal.fitTo(DT_hist);


	//TCanvas* can1 = new TCanvas("can1","",700,700);
	//can1->cd();
	//cFits[canCounter%3]->cd(1+canCounter/3);
	printf("cc is %i and %i \n",canCounter%3, 1+canCounter/3);
	canCounter++;

	RooPlot* plot1 = MassDT.frame(Title(" "));
	DT_hist.plotOn(plot1);
	fTotal.plotOn(plot1);
	plot1->SetMinimum(1e-05);
	plot1->Draw();

	val = (nGaus1.getVal()+nGaus2.getVal());

	return val;
}

//cutflag: 0 is tpc+tof (no bg), 1 is tpc (+bg), 2 is bg
void readTree_V0(Int_t nEvents=10, Int_t cutFlag=0, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();

	if (!makeChain(inputFile)) printf("Couldn't create the chain! \n", );
	else printf("Chain created with %i entries \n", mChain->GetEntries());

	//*bTracks 	= new TClonesArray("AliAnalysisPIDTrack");
	//*bV0s 		= new TClonesArray("AliAnalysisPIDV0");
	mChain->SetBranchAddress("AnalysisTrack",&bTracks);
	mChain->SetBranchAddress("AnalysisV0Track",&bV0s);
	mChain->SetBranchAddress("AnalysisEvent",&mEvent);

	/*const Int_t nPtBins = 59;
  	Double_t xBins[nPtBins+1] = { 0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 
    	0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 
    	0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
    	1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 
    	3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 
    	9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00, 20.00 };*/
	const Int_t nPtBins = 35;
  	Double_t xBins[nPtBins+1] = { 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 
    	1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 
    	3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00, 
    	9.00, 10.00, 11.00, 12.00, 13.00, 14.00 }; 


	TH1F* hEventMonitor			= new TH1F("hEventMonitor","",10,-0.5,9.5);
	TH1F* hV0Monitor			= new TH1F("hV0Monitor","",10,-0.5,9.5);
	TH1F* hTrackMonitor			= new TH1F("hTrackMonitor","",10,-0.5,9.5);
	TH2F* hV0TrCounter			= new TH2F("hV0TrCounter","",100,0,100,100,0,100);

	TH1F* hV0_IMK0s				= new TH1F("hV0_IMK0s","",2000,-0.2,0.2);
	TH1F* hV0_PtK0s				= new TH1F("hV0_PtK0s","",nPtBins,xBins);
	TH2F* hV0_IMPtK0s			= new TH2F("hV0_IMPtK0s","",2000,-0.2,0.2,nPtBins,xBins);
	TH1F* hYieldK0s				= new TH1F("hYieldK0s","",nPtBins,xBins);
	TH1F* hV0_IML				= new TH1F("hV0_IML","",2000,-0.2,0.2);
	TH1F* hV0_PtL				= new TH1F("hV0_PtL","",nPtBins,xBins);
	TH2F* hV0_IMPtL				= new TH2F("hV0_IMPtL","",2000,-0.2,0.2,nPtBins,xBins);
	TH1F* hYieldL				= new TH1F("hYieldL","",nPtBins,xBins);
	TH1F* hV0_IMAL				= new TH1F("hV0_IMAL","",2000,-0.2,0.2);
	TH1F* hV0_PtAL				= new TH1F("hV0_PtAL","",nPtBins,xBins);
	TH2F* hV0_IMPtAL			= new TH2F("hV0_IMPtAL","",2000,-0.2,0.2,nPtBins,xBins);
	TH1F* hYieldAL				= new TH1F("hYieldAL","",nPtBins,xBins);
	hV0_IMK0s->Sumw2();
	hV0_PtK0s->Sumw2();
	hYieldK0s->Sumw2();
	hV0_IML->Sumw2();
	hV0_PtL->Sumw2();
	hYieldL->Sumw2();
	hV0_IMAL->Sumw2();
	hV0_PtAL->Sumw2();
	hYieldAL->Sumw2();

	TH1F* hV0_DHasTPC			= new TH1F("hV0_DHasTPC","",200,0,10);
	hV0_DHasTPC->Sumw2();
	TH1F* hV0_DHasTOF			= new TH1F("hV0_DHasTOF","",200,0,10);
	hV0_DHasTOF->Sumw2();
	TH1F* hV0_DPt 				= new TH1F("hV0_DPt","",200,0,10);
	TH2F* hV0_DDTofPiPi			= new TH2F("hV0_DDTofPiPi","",300,-15,15,300,-15,15);
	TH2F* hV0_DDTofPiP			= new TH2F("hV0_DDTofPiP","",300,-15,15,300,-15,15);
	TH2F* hV0_DTofPivp			= new TH2F("hV0_DTofPivp","",100,0,10,300,-15,15);
	TH2F* hV0_DDDedx			= new TH2F("hV0_DDDedx","",300,0,300,300,0,300);
	TH2F* hV0_DDedxvp			= new TH2F("hV0_DDedxvp","",100,0,10,300,0,300);

	nEvents = (nEvents < mChain->GetEntries()) ? nEvents : mChain->GetEntries();
	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		hEventMonitor->Fill(0);
		if (iEv%100000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		bTracks->Clear();
		mChain->GetEntry(iEv);
		if (!mEvent) continue;
		hEventMonitor->Fill(1);
		//printf("event vz is %f \n", mEvent->GetVertexZ());

		if (iEv==0) mEvent->PrintEventSelection();
		if (!SelectEvent(mEvent)) continue;
		hEventMonitor->Fill(2);

		Int_t V0Count = 0; 
		Int_t nV0s = bV0s->GetEntriesFast();
		for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
			
			hV0Monitor->Fill(0);
			AliAnalysisPIDV0* v0 = (AliAnalysisPIDV0*)bV0s->At(iV0);
			if (!v0) continue;
			hV0Monitor->Fill(1);

			if (!SelectV0(v0)) continue;
			hV0Monitor->Fill(2);
			V0Count++;

			if (v0->GetPosAnalysisTrack()->HasTPCPID()) hV0_DHasTPC->Fill(v0->GetPosAnalysisTrack()->GetPt());
			if (v0->GetNegAnalysisTrack()->HasTPCPID()) hV0_DHasTPC->Fill(v0->GetNegAnalysisTrack()->GetPt());
			if (v0->GetPosAnalysisTrack()->HasTOFPID()) hV0_DHasTOF->Fill(v0->GetPosAnalysisTrack()->GetPt());
			if (v0->GetNegAnalysisTrack()->HasTOFPID()) hV0_DHasTOF->Fill(v0->GetNegAnalysisTrack()->GetPt());
			hV0_DPt->Fill(v0->GetPosAnalysisTrack()->GetPt());
			hV0_DPt->Fill(v0->GetNegAnalysisTrack()->GetPt());
			hV0_DDTofPiPi->Fill(v0->GetPosAnalysisTrack()->GetNSigmaPionTOF(),v0->GetNegAnalysisTrack()->GetNSigmaPionTOF());
			hV0_DDTofPiP->Fill(v0->GetPosAnalysisTrack()->GetNSigmaPionTOF(),v0->GetNegAnalysisTrack()->GetNSigmaProtonTOF());
			hV0_DTofPivp->Fill(v0->GetPosAnalysisTrack()->GetP(),v0->GetPosAnalysisTrack()->GetNSigmaPionTOF());
			hV0_DTofPivp->Fill(v0->GetNegAnalysisTrack()->GetP(),v0->GetNegAnalysisTrack()->GetNSigmaPionTOF());
			hV0_DDDedx->Fill(v0->GetPosAnalysisTrack()->GetTPCdEdx(),v0->GetNegAnalysisTrack()->GetTPCdEdx());
			hV0_DDedxvp->Fill(v0->GetPosAnalysisTrack()->GetP(),v0->GetPosAnalysisTrack()->GetTPCdEdx());
			hV0_DDedxvp->Fill(v0->GetNegAnalysisTrack()->GetP(),v0->GetNegAnalysisTrack()->GetTPCdEdx());


			bool noCuts = 0;
			if (noCuts || IsK0s(v0,cutFlag)) 	{
				hV0_IMK0s->Fill(v0->GetIMK0s());	
				hV0_PtK0s->Fill(v0->GetPt());
				hV0_IMPtK0s->Fill(v0->GetIMK0s(),v0->GetPt()); 		}
			if (noCuts || IsL(v0,cutFlag)) 		{
				hV0_IML->Fill(v0->GetIML());		
				hV0_PtL->Fill(v0->GetPt());
				hV0_IMPtL->Fill(v0->GetIML(),v0->GetPt());			}
			if (noCuts || IsAL(v0,cutFlag)) 	{
				hV0_IMAL->Fill(v0->GetIMAL());
				hV0_PtAL->Fill(v0->GetPt());
				hV0_IMPtAL->Fill(v0->GetIMAL(),v0->GetPt());		}
		}

		Int_t trCount = 0;
		Int_t nTracks = bTracks->GetEntriesFast();
		for (int iTr = 0; iTr < nTracks; ++iTr)	{
			
			hTrackMonitor->Fill(0);
			AliAnalysisPIDTrack* track = (AliAnalysisPIDTrack*)bTracks->At(iTr);
			if (!track) continue;
			hTrackMonitor->Fill(1);
			
			if (!SelectTrack(track)) continue;
			hTrackMonitor->Fill(2);
			trCount++;
		}
		
		hV0TrCounter->Fill(V0Count,trCount);
	}

	TCanvas* cFits[3] = new TCanvas(,"",1400,1000);
	for (int iC = 0; iC < 3; ++iC)
	{
		cFits[iC] = new TCanvas(Form("cFits%i",iC),"",1400,1000);
		cFits[iC]->Divide(7,5);
	}

	//ExtractYield((TH1D*)hV0_IMK0s);
	for (int iBin = 0; iBin < nPtBins; ++iBin)
	{
		//if (iBin!= 279) continue;
		hYieldK0s->SetBinContent(iBin,ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin)));
		hYieldL->SetBinContent(iBin,ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin)));
		hYieldAL->SetBinContent(iBin,ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin)));
	}

	hV0_PtK0s->Scale(1,"width");
	hV0_PtL->Scale(1,"width");
	hV0_PtAL->Scale(1,"width");
	hYieldK0s->Scale(1,"width");
	hYieldL->Scale(1,"width");
	hYieldAL->Scale(1,"width");
	hV0_DHasTPC->Divide(hV0_DPt);
	hV0_DHasTOF->Divide(hV0_DPt);

	TString path("$HOME/sq/pics/");



	printf(" WHAT IS UP \n", );
	//hEventMonitor->Draw();

	new TCanvas;
}