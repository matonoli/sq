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
AliAnalysisPIDEvent* mEvent;
TClonesArray* bTracks = 0;
TClonesArray* bV0s = 0;
TClonesArray* bParticles = 0;
TFile* mFout;
const Int_t pdgIds[] = {310, 3122, -3122};

//cutflag: 0 is tpc+tof (no bg), 1 is tpc (+bg), 2 is bg
void readEsd(Int_t nEvents=10, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();
	TList* lHist = gDirectory->GetList();
	int iHist = 0;

	if (!makeChain(inputFile)) printf("Couldn't create the chain! \n", );
	else printf("Chain created with %i entries \n", mChain->GetEntries());

	mChain->SetBranchAddress("AnalysisTrack",&bTracks);
	mChain->SetBranchAddress("AnalysisV0Track",&bV0s);
	mChain->SetBranchAddress("AnalysisEvent",&mEvent);
	if (flagMC) mChain->SetBranchAddress("AnalysisParticle",&bParticles);	

	nEvents = (nEvents < mChain->GetEntries()) ? nEvents : mChain->GetEntries();
	for (int iEv = 0; iEv < nEvents; ++iEv)	{
		
		hEventMonitor->Fill(0);
		if (iEv%10000==0) printf("Processing: %i out of total %i events...\n", iEv, nEvents);
		bTracks->Clear();
		mChain->GetEntry(iEv);
		if (!mEvent) continue;
		hEventMonitor->Fill(1);

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
				if (SelectV0_MC(p,0)) {
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

				Int_t mcLabel = p->GetLabel();

				Int_t rcV0Count = 0;
				Int_t nV0s = bV0s->GetEntriesFast();
				for (int iV0 = 0; iV0 < nV0s; ++iV0)	{
					
					AliAnalysisPIDV0* v0rc = (AliAnalysisPIDV0*)bV0s->At(iV0);
					if (!v0rc) continue;

					if (v0rc->GetMCPdgCode() != pdgIds[v0id]) continue;

					AliAnalysisPIDTrack* trPrc = v0rc->GetPosAnalysisTrack();
					AliAnalysisPIDTrack* trNrc = v0rc->GetNegAnalysisTrack();
					if (trPrc->GetMCMotherLabel() != mcLabel) continue;
					//if (trPrc->GetMCMotherPrimary() != 1) continue;
					rcV0Count++;
					if (rcV0Count>0) {
						printf("nEV %i nP %i v0 finds %i , v0id %i , mcpt %f , rcpt %f , mK0 %f , mL %f , tr+pt %f , tr-pt %f : \n", iEv, iP, rcV0Count, v0id, p->GetPt(), v0rc->GetPt(), v0rc->GetIMK0s(), v0rc->GetIML(), trPrc->GetPt(), trNrc->GetPt());
						printf("-------------------------- r %f , dcad %f , cos %f , eta %f , dcapv %f \n", v0rc->GetRadius(), v0rc->GetDCAV0Daughters(), v0rc->GetV0CosinePA(), v0rc->GetEta(), v0rc->GetDCAPV());
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

			bool noCuts = 0; float masscut = 0.0075;
			if (noCuts || IsK0s(v0,cutFlag)) 	{
				hV0_IMK0s->Fill(v0->GetIMK0s());	
				hV0_IMPtK0s->Fill(v0->GetIMK0s(),v0->GetPt());
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
				if (fabs(v0->GetIMAL())<masscut) {
					hV0_PtAL->Fill(v0->GetPt());
					Float_t delta = trP->GetNSigmaPionTOF() - trP->GetNSigmaPionTPC();
					hV0_DTofminTpcvrLpi->Fill(v0->GetRadius(),delta);
					delta = trN->GetNSigmaProtonTOF() - trN->GetNSigmaProtonTPC();
					hV0_DTofminTpcvrLpr->Fill(v0->GetRadius(),delta);		}
				}
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

	for (int iC = 0; iC < 3; ++iC)	{
		cFits[iC] = new TCanvas(Form("cFits%i",iC),"",2800,2000);
		cFits[iC]->Divide(7,5,0.0005,0.0005);	}

	//ExtractYield((TH1D*)hV0_IMK0s);
	for (int iBin = 1; iBin > nPtBins+1; ++iBin)
	{
		Float_t* y;
		y = ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin),0,0);
		hSBinK0s->SetBinContent(iBin,*(y+0));	
		hSBinK0s->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin),1,0);
		hSBoutK0s->SetBinContent(iBin,*(y+0));	
		hSBoutK0s->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtK0s->ProjectionX("x",iBin,iBin),2,0);	// 0 is underflow bin
		hYieldK0s->SetBinContent(iBin,*(y+0));	
		hYieldK0s->SetBinError(iBin,*(y+1));

		y = ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin),0,1);
		hSBinL->SetBinContent(iBin,*(y+0));	
		hSBinL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin),1,1);
		hSBoutL->SetBinContent(iBin,*(y+0));	
		hSBoutL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtL->ProjectionX("x",iBin,iBin),2,1);
		hYieldL->SetBinContent(iBin,*(y+0));	
		hYieldL->SetBinError(iBin,*(y+1));

		y = ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin),0,1);
		hSBinAL->SetBinContent(iBin,*(y+0));	
		hSBinAL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin),1,1);
		hSBoutAL->SetBinContent(iBin,*(y+0));	
		hSBoutAL->SetBinError(iBin,*(y+1));
		y = ExtractYield(hV0_IMPtAL->ProjectionX("x",iBin,iBin),2,1);
		hYieldAL->SetBinContent(iBin,*(y+0));	
		hYieldAL->SetBinError(iBin,*(y+1));
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

	TString path = Form("pics_%s/",outputFile);// ("$HOME/sq/pics/");//("$HOME/sq/pics/");
	path.ReplaceAll(".root","");
	gSystem->Exec(Form("mkdir %s", path.Data()));
	cFits[0]->SaveAs(path+"f_k0s.png");
	cFits[1]->SaveAs(path+"f_l.png");
	cFits[2]->SaveAs(path+"f_al.png");
	cFits[0]->SaveAs(path+"f_k0s.png");

	TCanvas* can1 = new TCanvas("can1","",1000,700);
	hV0_DHasTOF->Draw();
	can1->SaveAs(path+"eff_tof.png");
	hV0_DHasTPC->Draw();
	can1->SaveAs(path+"eff_tpc.png");

	TLegend *legpt = new TLegend(0.6,0.55,0.9,0.75);
	myLegendSetUp(legpt,0.028,1);
	legpt->AddEntry(hV0_PtK0s,"bin count in fixed range","l");
	legpt->AddEntry(hYieldK0s,"gaus+gaus+pol1 fit","l");
	legpt->AddEntry(hSBinK0s,"sideband with 3 and 6 RMS","l");
			
	hV0_PtK0s->Draw();
	can1->SetLogy();
	hYieldK0s->Draw("same");
	hSBinK0s->Draw("same");
	legpt->Draw();
	can1->SaveAs(path+"pt_k0s.png");
	hV0_PtL->Draw();
	hYieldL->Draw("same");
	hSBinL->Draw("same");
	legpt->Draw();
	can1->SaveAs(path+"pt_l.png");
	hV0_PtAL->Draw();
	hYieldAL->Draw("same");
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
}