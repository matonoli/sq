// V0 reconstruction macro from local aurora V0 trees
// OliverM 2018 Lund

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>

using namespace std;
using namespace RooFit;

// generate peaks (\w different templates), test efficiency + purity as function of statistics
// try different fit methods, vary for S/B too, can test sideband too

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
			new TCanvas;
			hist->Draw();
			new TCanvas;
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
			new TCanvas;
			hist->Draw();
			new TCanvas;
			val[0] = hist->Integral(hist->FindBin(mean-6.*rms),hist->FindBin(mean-3.*rms));
			val[0] += hist->Integral(hist->FindBin(mean+3.*rms),hist->FindBin(mean+6.*rms));
			printf("STATUS: int outside is %f \n", val[0]);
			val[1] = sqrt(val[0]);
			break;

		case 2 :	// fallthrough
		case 3 :
			
			RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",fitMin,fitMax);
			if (method==2) {	
				//hist->Rebin(8);
				RooDataHist DT_set("DT_set","DT_hist",MassDT,Import(*hist)); }
			else RooDataSet DT_set("DT_set","DT_tree",MassDT,Import(*tree));
		
			RooRealVar pGaus1A("pGaus1A","Mean 1",-0.004,0.004);
			RooRealVar pGaus1B("pGaus1B","Sigma 1",0,0.01);
			RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
			RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e06);
		
			RooRealVar pGaus2A("pGaus2A","Mean 2",-0.004,0.004);
			RooRealVar pGaus2B("pGaus2B","Sigma 2",0,0.01);
			RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus2A,pGaus2B); 
			RooRealVar nGaus2("nGaus2","N_{Gaus2}",0,0,1e06);
		
			RooRealVar pPolBgA("pPolBgA","Pol. par. A",0,-200,200);
			RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);//RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC));
			RooRealVar nPolBg("nPolBg","N_{PolBg}",0,0,1e06);
		
			//RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
			RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
			RooFitResult* fR = fTotal.fitTo(DT_set,Save(),Silence(),PrintLevel(-1));
		
			RooFormulaVar nGaus("nGaus","nGaus1+nGaus2",RooArgList(nGaus1,nGaus2));
			//printf("Errors are %f and %f, total is %f or %f wrt to %f \n", nGaus1.getError(), nGaus2.getError(), nGaus1.getError()+nGaus2.getError(),sqrt(nGaus1.getError()*nGaus1.getError()+nGaus2.getError()*nGaus2.getError()),nGaus.getPropagatedError(*fR));
		
			//cFits[canCounter%6]->cd(1+canCounter/6);
			new TCanvas;
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
			//leg1->AddEntry((TObject*)0,Form("%4.2f < p_{T} < %4.2f (GeV/#it{c})",xBins[canCounter/6],xBins[1+canCounter/6])," ");
			//leg1->AddEntry((TObject*)0,cNames[canCounter%6]+Form(" , #chi^{2}/ndf = %4.2f",plot1->chiSquare())," ");
			leg1->Draw();
	
			val[0] = nGaus.getVal();
			//printf("STATUS: int from fit is %f \n", val[0]);
			val[1] = nGaus.getPropagatedError(*fR);
			//canCounter++;
			break;

	}

	return val;
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

void roofit_test() {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();

	//generatePeak();
	RooRealVar MassDT("MassDT","#Delta m_{inv} (GeV/#it{c}^{2})",-0.03,0.03);
	RooRealVar pGaus1A("pGaus1A","Mean 1",0.,-0.004,0.004);
	RooRealVar pGaus1B("pGaus1B","Sigma 1",0.002,0,0.01);
	RooGaussian fGaus1("fGaus1","fGaus1",MassDT,pGaus1A,pGaus1B); 
	RooRealVar nGaus1("nGaus1","N_{Gaus1}",1,0,1e06);

	RooRealVar pGaus2A("pGaus2A","Mean 2",0.,-0.004,0.004);
	RooRealVar pGaus2B("pGaus2B","Sigma 2",0.001,0,0.01);
	RooGaussian fGaus2("fGaus2","fGaus2",MassDT,pGaus2A,pGaus2B); 
	RooRealVar nGaus2("nGaus2","N_{Gaus2}",0,0,1e06);

	RooRealVar pPolBgA("pPolBgA","Pol. par. A",1,-200,200);
	RooChebychev fPolBg("fPolBg","fPolBg",MassDT,pPolBgA);//RooArgSet(CB_DT_ParA,CB_DT_ParB,CB_DT_ParC));
	RooRealVar nPolBg("nPolBg","N_{PolBg}",0,0,1e06);
	//RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2),RooArgList(nGaus1,nGaus2));
	RooAddPdf fTotal("fTotal","fTotal",RooArgList(fGaus1,fGaus2,fPolBg),RooArgList(nGaus1,nGaus2,nPolBg));
	
	//shold be done separately for bg and sig
	RooDataSet* DT_set = fTotal.generate(MassDT,1000);
	MassDT.setBins(100);
	RooDataHist* DT_hist = DT_set->binnedClone();


	/*TCanvas* can1 = new TCanvas("can1","",1000,700);
	RooPlot* plot1 = MassDT.frame(Title(" "));
	//DT_set->plotOn(plot1,MarkerColor(kBlue),MarkerSize(0.4));
	DT_hist->plotOn(plot1,MarkerSize(0.4));
	fTotal.plotOn(plot1,LineWidth(1),LineColor(kRed));
	fTotal.plotOn(plot1,Components(fPolBg),LineStyle(kDashed),LineWidth(1),LineColor(kRed));
	plot1->Draw();
	*/

	TH1D* hM = DT_hist->createHistogram("hM",MassDT);
	hM->Fit("gaus");
	//hM->Draw();

	Float_t* y;
	y = ExtractYield(hM,0,0,1);
	printf("Yield from method 0 is %f #pm %f \n", *(y+0), *(y+1));
	y = ExtractYield(hM,0,1,1);
	printf("Yield from method 1 is %f #pm %f \n", *(y+0), *(y+1));
	//y = ExtractYield(hM,0,2,0);
	printf("Yield from method 2 is %f #pm %f \n", *(y+0), *(y+1));



	}