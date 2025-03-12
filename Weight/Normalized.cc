#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include "TBranch.h"
#include <string.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
	double  energy;
	double weight;
    TBranch *newBranch[6];
    TRandom *rand = new TRandom2(1);
    string file_name   = argv[1];
	TFile *file = new TFile(file_name.c_str(),"update");
	TTree *B4 = (TTree*)file->Get("B4");
	B4->SetBranchAddress("Energy"             , &energy);

	Long64_t nentries = B4->GetEntries();
    double total_weight = 0.0;
	newBranch[0] = B4->Branch("weight",&weight,"weight/D");
	for(Long64_t i=0;i<nentries; i++) 
	{
        B4->GetEntry(i);
        total_weight += pow(energy, -1.7);
	}
    // 计算归一化因子
    double norm_factor = nentries / total_weight;
    cout << "norm_factor" << norm_factor <<  endl;

    // 计算归一化权重
    for (Long64_t i = 0; i < nentries; i++) {
        B4->GetEntry(i);
        weight = pow(energy, -1.7) * norm_factor;
        newBranch[0]->Fill();
    }
	B4->Write("", TObject::kOverwrite); // save only the new version of the tree
	file->Close();
	exit(0);
}

