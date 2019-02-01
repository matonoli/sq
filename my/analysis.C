#include <iostream>
#include <fstream>

using namespace std;
using namespace RooFit;

TChain* mChain;
AliAnalysisPIDEvent* mEvent;
TFile* mFout;
TClonesArray* bTracks = 0;
TClonesArray* bV0s = 0;
TClonesArray* bParticles = 0;

int main() {

	TFile* fin = new TFile("AnalysisResult.root","READ");
}