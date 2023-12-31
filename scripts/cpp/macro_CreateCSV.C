/*
The main purpose of this script is to create a csv with the temporary series 
that will be object of study by the several algorithms. 

Obviously, this way we will reduce the quantity of data in a huge quantity. 
*/

// Libraries
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>


// Namespaces
using namespace std::filesystem;
using namespace std::chrono;
using namespace std;

void writeCSV(int *v[], int narrays, int dim, int event, int ntree,
              string foldername, string headers);
void AddToCSV_T(int v[], int dim, string idx, ofstream &file);


void macro_CreateCSV()
{
    // Ids of the PDs
    // PMT IDs in the TPC with X<0
    int realisticPMT_IDs_TPC0[60] = {6, 8, 10, 12, 14, 16, 36, 38, 40, 60, 62, 64, 66, 68, 70, 84, 86, 88, 90, 92, 94, 114, 116, 118, 138, 140, 142, 144, 146, 148, 162, 164, 166, 168, 170, 172, 
    192, 194, 196, 216, 218, 220, 222, 224, 226, 240, 242, 244, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304};
    
    // PMT IDs in the TPC with X>0
    int realisticPMT_IDs_TPC1[60] = {7, 9, 11, 13, 15, 17, 37, 39, 41, 61, 63, 65, 67, 69, 71, 85, 87, 89, 91, 93, 95, 115, 117, 119, 139, 141, 143, 145, 147, 149, 163, 165, 167, 169, 171, 173, 
    193, 195, 197, 217, 219, 221, 223, 225, 227, 241, 243, 245, 247, 249, 251, 271, 273, 275, 295, 297, 299, 301, 303, 305};
    
    // All PMTs IDs
    int realisticPMT_IDs[120] = {6, 8, 10, 12, 14, 16, 36, 38, 40, 60, 62, 64, 66, 68, 70, 84, 86, 88, 90, 92, 94, 114, 116, 118, 138, 140, 142, 144, 146, 148, 162, 164, 166, 168, 170, 172, 192, 
    194, 196, 216, 218, 220, 222, 224, 226, 240, 242, 244, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304, 7, 9, 11, 13, 15, 17, 37, 39, 41, 61, 63, 65, 67, 69, 71, 85, 87, 89, 91, 
    93, 95, 115, 117, 119, 139, 141, 143, 145, 147, 149, 163, 165, 167, 169, 171, 173, 193, 195, 197, 217, 219, 221, 223, 225, 227, 241, 243, 245, 247, 249, 251, 271, 273, 275, 295, 297, 299, 
    301, 303, 305};

    vector<int> selected_pmtid; // array to select vectors
    vector<int>::iterator it;
    const double shiftStamp = 135.;//ns
  	const double fSamplingTime = 2.;//ns
  	const int fBaseline=8000;//ADC

    // Preparing Directory with data preprocessed
    char data_preproc[256];
    getcwd(data_preproc, 256);
    strcat(data_preproc, "/data_preproc/");

    if (!exists(data_preproc))
        if (create_directories(data_preproc))
            cout << "Data preprocessed directory created" << endl;

    // raw data
    char data_raw[256];
	char skip[256] = "/Users/javigamero/MyMac/DS_Master/TFM/data/sample_particles_v2/.DS_Store";

  	getcwd(data_raw, 256); 
  	strcat(data_raw, "/data/sample_particles_v2"); // you must be in main folder of the project

    auto start = high_resolution_clock::now();

    int ntree_number = 0;

    strcat(data_preproc, "??.csv");
    ofstream file; // file where all temporary series will be written
    file.open(data_preproc);
    string isclosed="False";

  	for (const auto &folder: directory_iterator(data_raw)) // iterate above folders
  	{
  		if (folder.path() == skip) continue; // skip .DS_Store/ folder (only macOS)
        
        string ntree = folder.path();
		ntree = ntree.substr(63, ntree.length()); // take folder name
		ntree = ntree.substr(ntree.find("_")+1, ntree.length()); // take tree number

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
		    // Reading Variables 
		    // =================================================================            
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

            // particles info variables
            vector<vector<double>> *stepX=new vector<vector<double>>;
            vector<double> *dE = new vector<double>(); 
			tree->SetBranchAddress("stepX",&stepX);
            tree->SetBranchAddress("dE", &dE);

            // ideal variables
            vector<vector<double>> *SimPhotonsLiteVUV = new vector<vector<double>>;
            vector<vector<double>> *SimPhotonsLiteVIS = new vector<vector<double>>;
            tree -> SetBranchAddress("SimPhotonsLiteVUV", &SimPhotonsLiteVUV);
            tree -> SetBranchAddress("SimPhotonsLiteVIS", &SimPhotonsLiteVIS);

            // digitalised variables
            vector<vector<double> >* fSignalsDigi=new vector<vector<double> >();
		    tree->SetBranchAddress("SignalsDigi",&fSignalsDigi);
		    vector<double> * fStampTime=new vector<double>; tree->SetBranchAddress("StampTime", &fStampTime);
		    vector<int> * fOpChDigi=new vector<int>; tree->SetBranchAddress("OpChDigi",&fOpChDigi);
		    vector<double> x_raw, y_raw;

            // Deconvolutioned signal variables
		    vector<vector<double> >* fSignalsDeco=new vector<vector<double> >();
		    tree->SetBranchAddress("SignalsDeco",&fSignalsDeco);
		    vector<double> * fStampTimeDeco=new vector<double>; tree->SetBranchAddress("StampTimeDeco",&fStampTimeDeco);
		    vector<int> * fOpChDeco=new vector<int>; tree->SetBranchAddress("OpChDeco",&fOpChDeco);
		    vector<double> x_deco, y_deco;


            // =================================================================
		    // MAIN LOOP
		    // =================================================================
            for (int i=0; i<tree->GetEntries(); i++)
            {
                tree->GetEntry(i);
                double dE_e = 0.;
                double x_e = 0.;
                
                for (int j=0; j<stepX->size(); j++)
                {
                    // check we use the Michel electron:
                    if (PDGcode->at(j)!=-11) continue; 

                    // energy and position
                    dE_e = dE->at(j);
                    x_e = stepX->at(j).at(0); // 0 --> start position
                }
                
                // check we have an electron with energy 
                if (!(dE_e > 0)) continue; 

                // Now, we select the realistics PMT depending on the start position 
                if(x_e <= 0) 
                    selected_pmtid.assign(realisticPMT_IDs_TPC0, realisticPMT_IDs_TPC0 + 60); 
                else 
                    selected_pmtid.assign(realisticPMT_IDs_TPC1, realisticPMT_IDs_TPC1 + 60);


                // =============================================================
		        // IDEAL SIGNAL (VUV, VIS)
		        // =============================================================
                TH1D *signalVUV = new TH1D("","",1000,0,10000);
                TH1D *signalVIS = new TH1D("","",1000,0,10000);

                // obtain data from realistic PMTs
                for (int k=0; k<SimPhotonsLiteVUV->size(); k++)
                {
                    it = find(selected_pmtid.begin(), selected_pmtid.end(), k);
                    if (it == selected_pmtid.end()) continue; // .end() points outside, dismiss it
                    
                    for(int j=0; j<SimPhotonsLiteVUV->at(k).size(); j++)
				 	    signalVUV->Fill(SimPhotonsLiteVUV->at(k).at(j));

                    for(int j=0; j<SimPhotonsLiteVIS->at(k).size(); j++)
				 	    signalVIS->Fill(SimPhotonsLiteVIS->at(k).at(j));	
                }

                int dim = signalVIS->GetNbinsX();
                TAxis *X = signalVIS->GetXaxis();

                int timeSerieVUV[dim];
                int timeSerieVIS[dim];
                int idealTimeSerie[dim];
                int x[dim];
                
                for (int j=0; j<=dim+1; j++) // finishes in Nbins+1, if there are N bins, there will be N+1 points
                {
                    timeSerieVUV[j] = signalVUV -> GetBinContent(j);
                    timeSerieVIS[j] = signalVIS -> GetBinContent(j);
                    idealTimeSerie[j] = timeSerieVIS[j] + timeSerieVUV[j];
                }   

                for (int j=0; j<dim; j++)
                    x[j] = X->GetBinCenter(j+1); // Bin indexes start at 1

                // to write one csv for each serie
                /* 
                string foldername = string(data_preproc) + "LightSignal/";
                string header = "NPhotonsVUV,NPhotonsVIS,x\n";
                int *arr[] = {timeSerieVUV,timeSerieVIS,x};
                writeCSV(arr, 3, dim, i, ntree, foldername, header);
                */

                // to write all series in one csv
                string idx = ntree + "_" + to_string(event);
                AddToCSV_T(idealTimeSerie, dim, idx, file);

                // if (ntree==0 && isclosed == "False") # this loop is for x axis
                // {
                //     AddToCSV_T(x, dim, "0", file);
                //     file.close();
                //     isclosed="True";
                // }                


                // =============================================================
		        // DIGITALISED SIGNAL
		        // =============================================================

                // for(int k=0; k<fOpChDigi->size(); k++) // for every photon detector
		      	// {
				// 	double max_digi=-1;
					
				// 	//selecting the PMTs in the array we want to use
				// 	it = find (selected_pmtid.begin(), selected_pmtid.end(), fOpChDigi->at(k));
				// 	if (it == selected_pmtid.end()) continue;
					
					// // This loop is for DIGITISED signal
					// x_raw.clear(); y_raw.clear();    
					// for(int j=0; j<fSignalsDigi->at(k).size(); j++)
					// {
			  		// 	double rawADC = fSignalsDigi->at(k).at(j);
			  		// 	rawADC-=fBaseline;

			  		// 	if(max_digi < -rawADC) max_digi = -rawADC; // to limit the plot
			  
			  		// 	double t = fSamplingTime*j+fStampTime->at(k)*1000-shiftStamp;
			  		// 	if(t>-1000 && t<11000)
		  			// 	{
		    		// 		x_raw.push_back(t);
		    		// 		y_raw.push_back(-rawADC);
		  			// 	}	  
					// }

                    // const int dim = y_raw.size();
                    // if(!(dim>0)) continue;

                    // int y[dim];
                    // for (int j=0; j<dim; j++)
                    //     y[j] = y_raw[j];

                    // string idx = to_string(ntree) + "_" + to_string(event);
                    // AddToCSV_T(y, dim, idx, file); 	
                // }


                // =============================================================
		        // DECONVOLVED SIGNAL
		        // =============================================================

                for(int k=0; k<fOpChDeco->size(); k++) // for every photon detector
		      	{
                    double max_deco=-1;

                    //selecting the PMTs in the array we want to use
					it = find (selected_pmtid.begin(), selected_pmtid.end(), fOpChDigi->at(k));
					if (it == selected_pmtid.end()) continue;

                    // This loop is for DECONVOLVED/RECONSTRUCTED signal
					x_deco.clear(); y_deco.clear();    				
					for(int j=0; j<fSignalsDeco->at(k).size(); j++)
					{
		  				double decoADC = fSignalsDeco->at(k).at(j)/500.; 
		  				double t = fSamplingTime*j+fStampTimeDeco->at(k)*1000-shiftStamp; 
		  
		  				if(max_deco < decoADC) max_deco = decoADC; // to limit the plot
		  
		  				if(t>-1000 && t<11000)
		  				{
		    				x_deco.push_back(t);
		    				y_deco.push_back(decoADC);
		  				}	    
					}	  

                    const int dim = y_deco.size();
                    if(!(dim>0)) continue;

                    int y[dim];
                    for (int j=0; j<dim; j++)
                        y[j] = y_deco[j];

                    string idx = to_string(ntree) + "_" + to_string(event);
                    // AddToCSV_T(y, dim, idx, file); 
	        
                }
            
            } // end of main loop 
            ntree_number++;
            cout << "Trees completed: " << ntree_number << endl;

            delete tree; 
            delete fh;
            

        }// end of loop above .roots in a folder
    } // end of loop above folders

    file.close();

    auto stop = high_resolution_clock::now();
    cout << "Duration of the script: " << duration_cast<seconds>(stop-start).count() <<
         "(s)"<< endl;
} // end of main 


void writeCSV(int *v[], int narrays, int dim, int event, int ntree,
              string foldername, string headers)
{
    /* This function creates one csv with the arrays pointed in v.
    
    Variables: 
    - v: array of pointers, each element points to the desired array 
    - narrays: length of v, i.e. number of arrays that will be written 
    - dim: dimension of the arrays that will be written, it must be equal for all
    - ntree: index to indicate a number to write it in the csv name
    - foldername: folder where the csv will be saved
    - headers: header of the csv, separated by ","
    */

    // first, check if the folder exists. If it does not, create it.
    if (!exists(foldername))
        create_directories(foldername);

    ofstream file; 
    string filename = foldername + "tree_" + to_string(ntree) + "_event_" + to_string(event) + ".csv";
    file.open(filename);
    
    // Headers
    file << headers;

    // v is a an array of pointers to arrays
    for (int j=0; j<dim; j++)
    {    
        for (int i=0; i<narrays-1; i++)      
        {  
            // cout << i << " " << j;
            file << to_string(v[i][j]) + ",";
        }
        // cout << endl;
        file << to_string(v[narrays-1][j]) << "\n";
    }

    file.close();

    return ;
}

void AddToCSV_T(int v[], int dim, string idx, ofstream &file)
{
    /* This function adds an array to a csv as a line. Hence the csv will have
    "dim" columns. 

    Variables: 
    - v: array to write. ONLY ONE
    - dim: dimension of the array, also the final number of columns in the csv
    - file: ofstream file, it MUST BE OPENED OUTSIDE the function
    */

    file << idx + ";";
    for (int i=0; i<dim-1; i++) 
        file << to_string(v[i]) + ";";

    file << to_string(v[dim-1]) + "\n";

    return;

}
