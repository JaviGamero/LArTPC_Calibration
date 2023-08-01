// Libraries
#include <string>
#include <iostream>
#include <unistd.h>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TPolyLine3D.h>


// Namespaces
using namespace std::filesystem;
using namespace std;

// Functions declaration
void PlotHist(Double_t X[], int dim, string title, string XTitle, string YTitle);
void plotSpatialdE_2D(Double_t X[], Double_t Y[], Double_t Z[], int dim);
void plotSpatialdE_3D(Double_t X[], Double_t Y[], Double_t Z[], int dim);
void NPhotonsHist_PerPD(Double_t NPhotons[], int dim);

// MAIN
void macro_VarGraphs()
{
	char dir[256];
	char skip[256] = "/Users/javigamero/MyMac/DS_Master/TFM/data/sample_particles_v2/.DS_Store";

  	getcwd(dir, 256); 
  	strcat(dir, "/data/sample_particles_v2"); // you must be in main folder of the project

  	for (const auto &folder: directory_iterator(dir)) // iterate above folders
  	{
  		if (folder.path() == skip) continue; // skip .DS_Store/ folder (only macOS)
	  	for (const auto &treePath: directory_iterator(folder)) // above trees in folder
	  	{
	  		// =================================================================
		    // Reading trees
		    // =================================================================
		    char *path = new char[treePath.path().string().length()+1];
		    strcpy(path, treePath.path().string().c_str());

		    // cout << path << endl;
		    TFile *fh = new TFile(path);		 
			TTree *tree = fh -> Get<TTree>("opanatree/OpAnaTree"); // real data

			
			// =================================================================
		    // Variables 
		    // =================================================================

			// example to read a branch:
			// vector<vector<double> >* fSimPhotonsLiteVUV=new vector<vector<double> >();
			// tree->SetBranchAddress("SimPhotonsLiteVIS",&fSimPhotonsLiteVIS);

			// True neutrino interaction
			vector<double> *nuvE = new vector<double>();
			tree -> SetBranchAddress("nuvE", &nuvE);

			// Deposited Energy
			vector<double> *dE = new vector<double>(); 
			vector<vector<double>> * energydep = new vector<vector<double>>();
			tree->SetBranchAddress("dE", &dE);
			tree -> SetBranchAddress("energydep", &energydep);


			// IDs
			vector<int> *trackID=new vector<int>; 
    		vector<int> *motherID=new vector<int>; 
    		vector<int> *PDGcode=new vector<int>; 
    		vector<string> *process=new vector<string>; 
    		unsigned int event; 
    		tree->SetBranchAddress("trackID",&trackID);
    		tree->SetBranchAddress("motherID",&motherID);
    		tree->SetBranchAddress("PDGcode",&PDGcode);
    		tree->SetBranchAddress("process",&process);
    		tree->SetBranchAddress("eventID",&event);

			// coordinates of the deposited Energy
			vector<vector<double>> *stepX=new vector<vector<double>>;
			vector<vector<double>> *stepY=new vector<vector<double>>;
			vector<vector<double>> *stepZ=new vector<vector<double>>;
			vector<vector<double>> *stepT=new vector<vector<double>>;
			tree->SetBranchAddress("stepX",&stepX);
			tree->SetBranchAddress("stepY",&stepY);
			tree->SetBranchAddress("stepZ",&stepZ);
			tree->SetBranchAddress("stepT",&stepT);

			vector<vector<double>> *energydepX=new vector<vector<double>>;
			vector<vector<double>> *energydepY=new vector<vector<double>>;
			vector<vector<double>> *energydepZ=new vector<vector<double>>;
			tree->SetBranchAddress("energydepX",&energydepX);
			tree->SetBranchAddress("energydepY",&energydepY);
			tree->SetBranchAddress("energydepZ",&energydepZ);

			// Number of photons arriving
			vector<double> *NPhotonsPerPD_VUV = new vector<double>;
			vector<double> *NPhotonsPerPD_VIS = new vector<double>;
			tree -> SetBranchAddress("SimPhotonsperOpChVUV", &NPhotonsPerPD_VUV);
			tree -> SetBranchAddress("SimPhotonsperOpChVIS", &NPhotonsPerPD_VIS);

			// StampTime --> time start of the waveforms 
			vector<double> *StampTime = new vector<double>();
			vector<double> *StampTimeDeco = new vector<double>();
			tree -> SetBranchAddress("StampTime", &StampTime);
			tree -> SetBranchAddress("StampTimeDeco", &StampTimeDeco);

			// =================================================================
		    // Graphs
		    // =================================================================

		    for (int i=0; i<tree->GetEntries(); i++)
		    {
		    	tree->GetEntry(i);
		    	int dim;

				// =============================================================
			    // Neutrino interaction
			    // =============================================================
			    // dim = nuvE -> size();
			    // Double_t E[dim];
			    // for (int j=0; j<dim; j++)
			    // {
			    // 	E[j] = nuvE -> at(j);
			    // }		    	

			    // PlotHist(E, dim, "Interaction Energy Histogram", 
			    // 		 "Energy", "Entries");


      			// =============================================================
			    // Spatial Coordinates of the deposited Energy
			    // =============================================================
		    	for(int j=0; j<stepX->size(); j++) 
			    {
					cout << "PDG: " << PDGcode -> at(j) << endl;
					if(PDGcode->at(j) != -11) continue; // check we have a Michel electron after decay

					// dE_electron = dE->at(j); // (MeV) 
					if(!(dE->at(j) > 0)) continue; // we need the electron in the active volume

					cout<<"PDG:"<<PDGcode->at(j)<<"  dE:" << dE->at(j) <<endl;	
					cout<<"Event: "<<event<<endl;
					
			    	dim = stepX->at(j).size();
			    	Double_t X[dim], Y[dim], Z[dim], T[dim];

			    	for (int k=0; k<stepX->at(j).size(); k++)
			    	{
			    		X[k] = stepX->at(j).at(k);
			    		Y[k] = stepY->at(j).at(k);
			    		Z[k] = stepZ->at(j).at(k);
			    		// T[k] = stepT->at(j).at(k);

			    		// X[k] = energydepX->at(j).at(k);
			    		// Y[k] = energydepY->at(j).at(k);
			    		// Z[k] = energydepZ->at(j).at(k);
			    	}
			    	
			    	// plotSpatialdE_2D(X, Y, Z, dim);
			    	plotSpatialdE_3D(X, Y, Z, dim);
      			}

      			// =============================================================
			    // N Scintillation Photons deposited 
			    // =============================================================
			    // dim = NPhotonsPerPD_VUV -> size();
			    // Double_t photons[dim];
			    // for (int j=0; j<dim; j++)
			    // 	// photons[j] = NPhotonsPerPD_VUV->at(j);
			    // 	photons[j] = NPhotonsPerPD_VIS->at(j);

      			// // NPhotonsHist_PerPD(photons, dim);
      			// PlotHist(photons, dim, "Number of photons histogram", 
      			// 		 "#Photons", "Entries");

      			// =============================================================
			    // Start time of the waveforms
			    // =============================================================
		    	// dim = StampTime -> size();
		    	// dim = StampTimeDeco -> size();
		    	// Double_t t[dim];
		    	// for (int j=0; j<dim; j++)
		    	// 	t[j] = StampTime -> at(j);
		    	// 	// t[j] = StampTimeDeco -> at(j);

		    	// PlotHist(t, dim, "StampTime", "Start Time (#mus)", "Entries");


		    } // loop over entries in one tree
		    delete tree;
	    	delete fh;

		} // loop over .roots in a folder
		

  	} // loop over folders

} // macro end










// FUNCTINONS
void PlotHist(Double_t X[], int dim, string title, string XTitle, string YTitle)
{
	char *cTitle = new char[title.length()+1]; strcpy(cTitle, title.c_str());
	char *cXTitle = new char[XTitle.length()+1]; strcpy(cXTitle, XTitle.c_str());
	char *cYTitle = new char[YTitle.length()+1]; strcpy(cYTitle, YTitle.c_str());

	TCanvas *canHist = new TCanvas("", cTitle, 200, 200, 800, 500);
	TH1F *Hist = new TH1F("", cTitle, 100, 
						  *min_element(X, X+dim), 
						  *max_element(X, X+dim));

	for (int i=0; i<dim; i++)
		Hist -> Fill(X[i]);

	Hist -> SetXTitle(cXTitle);
	Hist -> SetYTitle(cYTitle);
	Hist -> Draw();	
	
	canHist -> Update();
	canHist -> Modified();
	canHist -> WaitPrimitive();

	delete Hist;

	return;
}


void plotSpatialdE_2D(Double_t X[], Double_t Y[], Double_t Z[], int dim)
{
	TCanvas *canSpatialdE = new TCanvas("Spatial dE", "Spatial dE", 200,200,1800,600);
   
	canSpatialdE -> Divide(3,1);

	canSpatialdE -> cd(1);
	TGraph *XvsY = new TGraph(dim, X, Y);
	XvsY -> SetTitle("");
	XvsY -> GetXaxis() -> SetTitle("x (cm)");
	XvsY -> GetYaxis() -> SetTitle("y (cm)");
	XvsY -> SetLineColor(2);
	XvsY -> SetLineWidth(4);
	XvsY -> SetMarkerColor(4);
	XvsY -> SetMarkerStyle(21);
	XvsY -> Draw("al");

	canSpatialdE -> cd(2);
	TGraph *XvsZ = new TGraph(dim, X, Z);
	XvsZ -> SetTitle("Spatial Coordinates of dE");
	XvsZ -> GetXaxis() -> SetTitle("x (cm)");
	XvsZ -> GetYaxis() -> SetTitle("z (cm)");
	XvsZ -> SetLineColor(2);
	XvsZ -> SetLineWidth(4);
	XvsZ -> SetMarkerColor(4);
	XvsZ -> SetMarkerStyle(21);
	XvsZ -> Draw("al");

	canSpatialdE -> cd(3);
	TGraph *YvsZ = new TGraph(dim, Y, Z);
	YvsZ -> SetTitle("");
	YvsZ -> GetXaxis() -> SetTitle("y (cm)");
	YvsZ -> GetYaxis() -> SetTitle("z (cm)");
	YvsZ -> SetLineColor(2);
	YvsZ -> SetLineWidth(4);
	YvsZ -> SetMarkerColor(4);
	YvsZ -> SetMarkerStyle(21);
	YvsZ -> Draw("al");	

	canSpatialdE -> Update();
	canSpatialdE -> Modified();
	canSpatialdE -> WaitPrimitive();

	delete XvsY;
	delete XvsZ;
	delete YvsZ;

	return;
}

void plotSpatialdE_3D(Double_t X[], Double_t Y[], Double_t Z[], int dim)
{
	TCanvas *canSpatialdE = new TCanvas("Spatial dE", "Spatial dE", 200,200,700,700);
	TH3F* axis = new TH3F("h3", "Spatial dE", 
						20, *min_element(X, X+dim), *max_element(X, X+dim), 
						20, *min_element(Y, Y+dim), *max_element(Y, Y+dim), 
						20, *min_element(Z, Z+dim), *max_element(Z, Z+dim)); 
	axis -> SetStats(0);
	axis->Draw();
	
   	
	TPolyLine3D *graph = new TPolyLine3D(dim); // 3D line plot object
	for (int i=0; i<dim; i++)
	{
		graph -> SetPoint(i, X[i], Y[i], Z[i]); // set points to the objects
	}
	graph -> SetLineColor(kRed+1);
   	graph -> SetLineWidth(4);
   	// graph -> GetXaxis()->SetTitle("X (cm)"); // not this way
   	// graph -> GetYaxis()->SetTitle("Y (cm)");
   	// graph -> GetZaxis()->SetTitle("Z (cm)");
	graph -> Draw();

	canSpatialdE -> Update();
	canSpatialdE -> Modified();
	canSpatialdE -> WaitPrimitive();

	delete graph;

	return; 
}

void NPhotonsHist_PerPD(Double_t NPhotons[], int dim)
{
	TCanvas *canHist = new TCanvas("", "Photons per PD",200, 200, 800, 500);
	TGraph *graph = new TGraph(dim, NPhotons);

	graph -> SetTitle("Photons per PD");
	graph -> GetXaxis() -> SetTitle("PD"); 
	graph -> GetYaxis() -> SetTitle("Number of photons"); 
	graph -> Draw();

	canHist -> Update();
	canHist -> Modified();
	canHist -> WaitPrimitive();

	delete graph;	

	return;
}


