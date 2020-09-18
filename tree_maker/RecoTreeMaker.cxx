// Main script
void RecoTreeMaker()
{

  // Set up output tree
  TFile tree_file("gd_water_reco.root", "recreate");
  TTree *tree = new TTree("tree", "");

  // Define branch variable
  std::string int_name;
  double mc_x, mc_y, mc_z;
  double reco_x, reco_y, reco_z;
  double mc_time, mc_energy;
  double inner_hits, veto_hits;
  double reco_pe, reco_time;
  double n9, good_pos, closest_pmt;
  double nevents;

  // Define branches
  tree->Branch("int_name", &int_name);  
  tree->Branch("mc_x", &mc_x);  
  tree->Branch("mc_y", &mc_y);  
  tree->Branch("mc_z", &mc_z);  
  tree->Branch("reco_x", &reco_x);  
  tree->Branch("reco_y", &reco_y);  
  tree->Branch("reco_z", &reco_z);  
  tree->Branch("mc_energy", &mc_energy);  
  tree->Branch("mc_time", &mc_time);  
  tree->Branch("inner_hits", &inner_hits);  
  tree->Branch("veto_hits", &veto_hits);  
  tree->Branch("reco_pe", &reco_pe);  
  tree->Branch("reco_time", &reco_time);  
  tree->Branch("n9", &n9);  
  tree->Branch("good_pos", &good_pos);  
  tree->Branch("closest_pmt", &closest_pmt);  
  tree->Branch("nevents", &nevents);  

  // List of all directories to include
  std::vector<TString> files = {
    "IBDPositronHeyshamSig_LIQUID_ibd_p_hs", 
    "IBDPositronHeyshamBkg_LIQUID_ibd_p_hb",
    "208Tl_LIQUID_CHAIN_232Th_NA", 
    "212Bi_LIQUID_CHAIN_232Th_NA", 
    "214Pb_LIQUID_CHAIN_222Rn_NA", 
    "40K_PMT_40K_NA", 
    "208Tl_PMT_CHAIN_232Th_NA", 
    "212Bi_PMT_CHAIN_232Th_NA", 
    "214Pb_PMT_CHAIN_238U_NA", 
    "IBDNeutron_LIQUID_ibd_n", 
    "210Bi_LIQUID_CHAIN_222Rn_NA", 
    "212Pb_LIQUID_CHAIN_232Th_NA", 
    "228Ac_LIQUID_CHAIN_232Th_NA", 
    "IBDPositron_LIQUID_ibd_p", 
    "210Bi_PMT_CHAIN_238U_NA", 
    "212Pb_PMT_CHAIN_232Th_NA", 
    "228Ac_PMT_CHAIN_232Th_NA", 
    "IBD_LIQUID_pn_ibd", 
    "210Tl_LIQUID_CHAIN_222Rn_NA", 
    "214Bi_LIQUID_CHAIN_222Rn_NA", 
    "234Pa_PMT_CHAIN_238U_NA", 
    "li9_LIQUID_A_Z",
    "210Tl_PMT_CHAIN_238U_NA", 
    "214Bi_PMT_CHAIN_238U_NA", 
    "40K_LIQUID_40K_NA",
    "n17_LIQUID_A_Z"};

  // Loop over the input directories
  for(int f = 0; f < files.size(); f++){

    int_name = files[f];
    // Parent directory
    TString directory = "/usr/WS2/brooks50/heysham_analysis/gd_water/bonsai_root_files_default/";
    TString fname = "merged_Watchman_"+files[f]+".root";
    std::cout<<fname<<"\n";

    // Get PMT info from the file
    TFile b_file(directory+fname); //path to your rat-pac files

    // Get the number of events
    TTreeReader sumreader("runSummary", &b_file);
    TTreeReaderValue<int> nEvents(sumreader, "nEvents");
    nevents = 0;
    while(sumreader.Next()){
      nevents += *nEvents;
    }

    // Set up a TTreeReader to access the tree information
    TTreeReader reader("data", &b_file);

    TTreeReaderValue<double> x(reader, "x");
    TTreeReaderValue<double> y(reader, "y");
    TTreeReaderValue<double> z(reader, "z");
    TTreeReaderValue<double> t(reader, "t");
    TTreeReaderValue<double> mcx(reader, "mcx");
    TTreeReaderValue<double> mcy(reader, "mcy");
    TTreeReaderValue<double> mcz(reader, "mcz");
    TTreeReaderValue<double> mcenergy(reader, "mc_energy");
    TTreeReaderValue<double> mct(reader, "mct");
    TTreeReaderValue<double> reco_n9(reader, "n9");
    TTreeReaderValue<double> pe(reader, "pe");
    TTreeReaderValue<int> subid(reader, "subid");
    TTreeReaderValue<int> inner_hit(reader, "inner_hit");
    TTreeReaderValue<int> veto_hit(reader, "veto_hit");
    TTreeReaderValue<double> goodpos(reader, "good_pos");
    TTreeReaderValue<double> closestPMT(reader, "closestPMT");

    // Loop through the events
    while(reader.Next()){

      // Only look at the first triggered event (skip retriggers)
      if(*subid != 0) continue;
      // Fill the tree
      mc_x = *mcx; 
      mc_y = *mcy;
      mc_z = *mcz;
      reco_x = *x; 
      reco_y = *y;
      reco_z = *z;
      mc_time = *mct;
      mc_energy = *mcenergy;
      n9 = *reco_n9;
      inner_hits = *inner_hit;
      veto_hits = *veto_hit;
      reco_pe = *pe;
      reco_time = *t;
      good_pos = *goodpos;
      closest_pmt = *closestPMT;

      tree->Fill();

    }

    std::cout<<"Number of events = "<<nevents<<"\n";

  }

  // Write the tree to file
  tree_file.cd();
  tree->Write();
  tree_file.Close();
}
