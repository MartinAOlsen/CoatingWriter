%------------------------ Analyzing the resulting optimum ------

%options_single_home.compile=0;
options_single_home.gravitation=1;
options_single_home.ncount=1*1e7;
options_single_home.mpi=1;
options_single_home.dir=[cpath '/' filename 'waveALL'];

%optimal
optimal=monitor(1).Data.Parameters;

% rest of options for cluster
if (select==2); load NUMCORES.DAT; warning off; else NUMCORES=0; end;


% general opptions
options_single_cluster.mpi=NUMCORES;
options_single_cluster.Display=[];
options_single_cluster.dir=[filename 'waveALL'];
options_single_cluster.gravitation=1;
options_single_cluster.compile=0;
options_single_cluster.ncount=1*5e8;

options_single={options_single_home options_single_cluster};

optimal.WaveMin='0.1';
optimal_ess.WaveMin='0.1';
optimal.WaveMax='8';
optimal_ess.WaveMax='8';
MaxWB=str2num(optimal.WaveMax);

% Run for all wavelengths
options_single{select}.dir=[filename 'waveLarge'];
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];

monitor_ALLW=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});

save([filename '_all.mat']);
