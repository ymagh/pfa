//----Script for computing inactive channels
//----from merged CRUZET DQM root files - Digi occupany plots
//----Root files available at: 
//eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/Calibration/CRUZET
//----Y.M - updated November 05, 2021

{	
	//-------Open files in order to merge digi occupancy plots
	//CRUZET runs
	TFile *f1 = new TFile("342810.root");
	TFile *f2 = new TFile("342966.root");
	TFile *f3 = new TFile("343034.root");
	TFile *f4 = new TFile("343082.root");
	TFile *f5 = new TFile("343171.root");
	TFile *f6 = new TFile("343266.root");
	TFile *f7 = new TFile("343387.root");
	TFile *f8 = new TFile("343496.root");
	TFile *f9 = new TFile("343498.root");
	TFile *f10 = new TFile("343621.root");
	TFile *f11 = new TFile("343642.root");
	TFile *f12 = new TFile("343677.root");
	
	TH2F *h1[36];
	TH2F *h2[36];
	TH2F *h3[36];
	TH2F *h4[36];
	TH2F *h5[36];
	TH2F *h6[36];
	TH2F *h7[36];
	TH2F *h8[36];
	TH2F *h9[36];
	TH2F *h10[36];
	TH2F *h11[36];
	TH2F *h12[36];

	char name[36];
	char name2[36];
	TList *list[36];
	TH1F *h[36];
	
	TDirectory* dir1 = (TDirectory*)f1->Get("DQMData/Run 342810/GEM/Run summary/digi");
	TDirectory* dir2 = (TDirectory*)f2->Get("DQMData/Run 342966/GEM/Run summary/digi");
	TDirectory* dir3 = (TDirectory*)f3->Get("DQMData/Run 343034/GEM/Run summary/digi");
	TDirectory* dir4 = (TDirectory*)f4->Get("DQMData/Run 343082/GEM/Run summary/digi");
	TDirectory* dir5 = (TDirectory*)f5->Get("DQMData/Run 343171/GEM/Run summary/digi");
	TDirectory* dir6 = (TDirectory*)f6->Get("DQMData/Run 343266/GEM/Run summary/digi");
	TDirectory* dir7 = (TDirectory*)f7->Get("DQMData/Run 343387/GEM/Run summary/digi");
	TDirectory* dir8 = (TDirectory*)f8->Get("DQMData/Run 343496/GEM/Run summary/digi");
	TDirectory* dir9 = (TDirectory*)f9->Get("DQMData/Run 343498/GEM/Run summary/digi");
	TDirectory* dir10 = (TDirectory*)f10->Get("DQMData/Run 343621/GEM/Run summary/digi");
	TDirectory* dir11 = (TDirectory*)f11->Get("DQMData/Run 343642/GEM/Run summary/digi");
	TDirectory* dir12 = (TDirectory*)f12->Get("DQMData/Run 343677/GEM/Run summary/digi");
	
	for (int i = 1; i <= 36; i++) {
		
		if (i<10) sprintf(name, "strip_occ_GE-11_L2_ch0%d",i); //choose +/- endcap, L1/2
		
		else if (i>=10) sprintf(name, "strip_occ_GE-11_L2_ch%d",i); //adapt the plots name in CRAFT
	
		h1[i] = (TH2F*)dir1->Get(name);
		h2[i] = (TH2F*)dir2->Get(name);
		h3[i] = (TH2F*)dir3->Get(name);
		h4[i] = (TH2F*)dir4->Get(name);
		h5[i] = (TH2F*)dir5->Get(name);
		h6[i] = (TH2F*)dir6->Get(name);
		h7[i] = (TH2F*)dir7->Get(name);
		h8[i] = (TH2F*)dir8->Get(name);
		h9[i] = (TH2F*)dir9->Get(name);
		h10[i] = (TH2F*)dir10->Get(name);
		h11[i] = (TH2F*)dir11->Get(name);
		h12[i] = (TH2F*)dir12->Get(name);
		
		//create lists
		list[i] = new TList;
		list[i]->Add(h1[i]);
		list[i]->Add(h2[i]);
		list[i]->Add(h3[i]);
		list[i]->Add(h4[i]);
		list[i]->Add(h5[i]);
		list[i]->Add(h6[i]);
		list[i]->Add(h7[i]);
		list[i]->Add(h8[i]);
		list[i]->Add(h9[i]);
		list[i]->Add(h10[i]);
		list[i]->Add(h11[i]);
		list[i]->Add(h12[i]);
		
		//clone and merge histograms
		sprintf(name2, "h%d",i);
		h[i] = (TH1F*)h1[i]->Clone("name2");
		h[i]->Reset();
		h[i]->Merge(list[i]);
		list[i]->Delete();
	}
	
//-------read scurve files and fill array of inactive channels in scurves	
//Dynamically allocate memory for 3D Array 
#define X 36
#define Y 24
#define Z 384
int*** A = new int**[X];
for (int i = 1; i <= X; i++)
{
	A[i] = new int*[Y];
	for (int j = 0; j < Y; j++) {
		A[i][j] = new int[Z];
	}
}	

// assign values to the allocated memory
for (int i = 1; i <= X; i++)
{
	for (int j = 0; j < Y; j++)
	{
		for (int k = 0; k < Z; k++) {
			A[i][j][k] = -1; //first assign all table values to zero
		}
	}
}

auto Ch=0;
string det;
auto vfat=0;
auto channel=0;
auto strip=0;
ifstream myfile("scurve-M-L2.txt");//(open scurve file - 4 scurve files for M/P endcap, L1/2)   
string line;

while(getline(myfile, line)){	
	stringstream ss(line);
	ss>>det>>Ch>>vfat>>channel>>strip;
	A[Ch][vfat][strip]=0; //only inactive channels in scurves will have 0 value, all others: -1
}

//-------find inactive channels in data
int n, s, v;//hits/strip/vfat
//counter of inactive strips
int t1 = 0;
int t2 = 0;
int t3 = 0;

//Open output file to write inactive strips in data
ofstream filep;
filep.open("test.csv"); // Create file

for (int ch = 1; ch <= 36; ch++)  { //loop on chambers
	for (int i = 0; i < 384; i++)  { //loop on strips
		for (int j = 1; j <= 8; j++){ //loop on eta partition
			
			n = h[ch]->GetBinContent(i,j); //hits in the strip
			s= i%128;                      //strip number in format 1-128
			v= (int(i/128)*8)+(8-j);       //vfat number

			//strips inactive in scurves but active in data
			if (A[ch][v][i]==0 && n!=0)
			{
				t1++;		
				cout<<"t2:"<<t2<<", chamber:"<<ch<<", eta:"<<j<<", vfat:"<<v<<" ,strip:"<<i<<endl;
				//filep<<ch<<","<<j<<","<<v<<","<<i<<endl;
			}
			//strips inactive in data but active in scurves
			if (A[ch][v][i]!=0 && n==0)
			{
				t2++;
				cout<<"t2:"<<t2<<", chamber:"<<ch<<", eta:"<<j<<", vfat:"<<v<<", strip:"<<i<<endl;
				//filep<<ch<<","<<j<<","<<v<<","<<i<<endl;
			}
			//strips inactive in both data and scurves
			if (A[ch][v][i]==0 && n==0)
			{
				t3++;
				cout<<"t3:"<<t3<<", chamber:"<<ch<<", eta:"<<j<<" ,vfat:"<<v<<" ,strip:"<<i<<endl;
				//filep<<ch<<","<<j<<","<<v<<","<<i<<endl;
			}
		}//end eta partition loop
	}//end strips loop
}//end chamber loop

//print percentage inactive channels
cout<<"inactive in scurve:"<<t1<<", %:"<<double(t1)/110592*100<<endl;
cout<<"inactive in data  :"<<t2<<", %:"<<double(t2)/110592*100<<endl;
cout<<"inactive in both:"<<t3<<", %:"<<double(t3)/110592*100<<endl;

//close output file
filep.close();

//deallocate memory
for (int i = 1; i <= X; i++)
{
	for (int j = 0; j < Y; j++) {
		delete[] A[i][j];
	}
	delete[] A[i];
}
delete[] A;
}
