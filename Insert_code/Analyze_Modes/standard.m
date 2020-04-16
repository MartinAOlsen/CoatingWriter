%------------------------ Analyzing the resulting optimum ------

%options_single_home.compile=0;
options_single_home.gravitation=1;
options_single_home.ncount=1*1e7;
options_single_home.mpi=1;
options_single_home.dir=[cpath '/' filename 'waveALL'];


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

wavecenters=[ 3 3 3 3 3 ];
snapwidth=0.01*ones(1,5);
optimal=monitor(1).Data.Parameters;

optimal_ess = optimal;
names_ess=fieldnames(optimal_ess);
Foptimal = fieldnames(optimal);
Fp = fieldnames(p);
for i = 1:length(Fp)
    isField = 0;
    for j = 1:length(Foptimal)
        if strcmp(Foptimal{j},Fp{i})
            isField = 1;
        end
    end
   if isField==0
        eval(['optimal.' Fp{i} '= p.' Fp{i} ';'])
    end
end
names=fieldnames(optimal);

% Run for fom wavelengths
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat']
monitor_fom=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});

optimal.WaveMin='0.1';
optimal_ess.WaveMin='0.1';
optimal.WaveMax='6';
optimal_ess.WaveMax='6';
MaxWB=str2num(optimal.WaveMax);

% Run for all wavelengths
options_single{select}.dir=[filename 'waveLarge'];
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];
monitor_ALLW=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});

% Run for rough ESS source spectrum
options_single{select}.dir=[filename 'waveESS'];
monitor_ESSW=mcstas([instrument_name '_analyze_ess.instr'],optimal_ess,options_single{select});

% Run for robustness check

degraded=optimal;

for i=1:length(names)
tmpname=names{i}
 if strcmp(tmpname(1:1),'m') && length(tmpname)<4
   if ischar(optimal.(names{i}))
    degraded.(names{i})=0.8*str2num(optimal.(names{i}))
   else
    degraded.(names{i})=0.8*optimal.(names{i})
   end
 end
 if length(tmpname)>5
  if strcmp(tmpname(1:5),'alpha')
   if ischar(optimal.(names{i}))
     degraded.(names{i})=1.4*str2num(optimal.(names{i}))
   else
     degraded.(names{i})=1.4*optimal.(names{i})
   end
  end
 end
end
degraded_ess=optimal_ess;
for i=1:length(names_ess)
tmpname=names_ess{i}
 if strcmp(tmpname(1:1),'m') && length(tmpname)<4
   if ischar(optimal_ess.(names_ess{i}))
    degraded_ess.(names_ess{i})=0.8*str2num(optimal_ess.(names_ess{i}))
   else
    degraded_ess.(names_ess{i})=0.8*optimal_ess.(names_ess{i})
   end
 end
 if length(tmpname)>5
  if strcmp(tmpname(1:5),'alpha')
   if ischar(optimal_ess.(names_ess{i}))
     degraded_ess.(names_ess{i})=1.4*str2num(optimal_ess.(names_ess{i}))
   else
     degraded_ess.(names_ess{i})=1.4*optimal_ess.(names_ess{i})
   end
  end
 end
end

% Run for all wavelengths
options_single{select}.dir=[filename 'Large_degraded'];
degraded_visualizer = degraded;
degraded_visualizer.file_name = [filename '_geometry.dat'];
monitor_ALLW_degraded=mcstas([instrument_name '_analyze.instr'],degraded_visualizer,options_single{select});

% Run for rough ESS source spectrum
options_single{select}.dir=[filename 'ESS_degraded'];
monitor_ESSW_degraded=mcstas([instrument_name '_analyze_ess.instr'],degraded_ess,options_single{select});

% Run for individual wavelength snapshots
options_single_cluster.ncount=1*1e8;
options_single={options_single_home options_single_cluster};

options_single{select}.dir=[filename 'wave1'];
optimal.WaveMin=wavecenters(1) - snapwidth(1)*0.5;
optimal.WaveMax=wavecenters(1) + snapwidth(1)*0.5;
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];
monitor_W.wave1=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});
options_single{select}.dir=[filename 'wave_ess1'];
optimal_ess.WaveMin=wavecenters(1) - snapwidth(1)*0.5;
optimal_ess.WaveMax=wavecenters(1) + snapwidth(1)*0.5;
monitor_W_ess.wave1=mcstas([instrument_name '_analyze_ess.instr'],optimal_ess,options_single{select});


options_single{select}.dir=[filename 'wave2'];
optimal.WaveMin=wavecenters(2) - snapwidth(2)*0.5;
optimal.WaveMax=wavecenters(2) + snapwidth(2)*0.5;
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];
monitor_W.wave2=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});
options_single{select}.dir=[filename 'wave_ess2'];
optimal_ess.WaveMin=wavecenters(2) - snapwidth(2)*0.5;
optimal_ess.WaveMax=wavecenters(2) + snapwidth(2)*0.5;
monitor_W_ess.wave2=mcstas([instrument_name '_analyze_ess.instr'],optimal_ess,options_single{select});


options_single{select}.dir=[filename 'wave3'];
optimal.WaveMin=wavecenters(3) - snapwidth(3)*0.5;
optimal.WaveMax=wavecenters(3) + snapwidth(3)*0.5;
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];
monitor_W.wave3=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});
options_single{select}.dir=[filename 'wave_ess3'];
optimal_ess.WaveMin=wavecenters(3) - snapwidth(3)*0.5;
optimal_ess.WaveMax=wavecenters(3) + snapwidth(3)*0.5;
monitor_W_ess.wave3=mcstas([instrument_name '_analyze_ess.instr'],optimal_ess,options_single{select});


options_single{select}.dir=[filename 'wave4'];
optimal.WaveMin=wavecenters(4) - snapwidth(4)*0.5;
optimal.WaveMax=wavecenters(4) + snapwidth(4)*0.5;
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];
monitor_W.wave4=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});
options_single{select}.dir=[filename 'wave_ess4'];
optimal_ess.WaveMin=wavecenters(4) - snapwidth(4)*0.5;
optimal_ess.WaveMax=wavecenters(4) + snapwidth(4)*0.5;
monitor_W_ess.wave4=mcstas([instrument_name '_analyze_ess.instr'],optimal_ess,options_single{select});


options_single{select}.dir=[filename 'wave5'];
optimal.WaveMin=wavecenters(5) - snapwidth(5)*0.5;
optimal.WaveMax=wavecenters(5) + snapwidth(5)*0.5;
optimal_visualizer = optimal;
optimal_visualizer.file_name = [filename '_geometry.dat'];
monitor_W.wave5=mcstas([instrument_name '_analyze.instr'],optimal_visualizer,options_single{select});
options_single{select}.dir=[filename 'wave_ess5'];
optimal_ess.WaveMin=wavecenters(5) - snapwidth(5)*0.5;
optimal_ess.WaveMax=wavecenters(5) + snapwidth(5)*0.5;
monitor_W_ess.wave5=mcstas([instrument_name '_analyze_ess.instr'],optimal_ess,options_single{select});



save([filename '_all.mat']);
save([cpath '/../output/analysis/' instrument_name scanname '_all.mat']);
copyfile([filename '_geometry.dat'],[cpath '/../output/analysis/' filename '_geometry.dat'])

fid = fopen([cpath '/../output/analysis/analyze_all_ifit.m'], 'a');
fprintf(fid,['clear all;clc;close all;\n' instrument_name scanname '_ifit_analyse\n']);
fclose(fid);

string='';
for i=1:length(names)
string=[string  names{i} '=' num2str(optimal.(names{i})) '\n'];
end

fid = fopen([cpath '/' filename '.par'], 'w');
fprintf(fid,string);
fclose(fid);

string='';
for i=1:length(names_ess)
string=[string  names_ess{i} '=' num2str(optimal_ess.(names_ess{i})) '\n'];
end

fid = fopen([cpath '/' filename '_ess.par'], 'w');
fprintf(fid,string);
fclose(fid);

flux=monitor_fom(11).Data.values(1);

fid = fopen([cpath '/../master_record-done' scanname '.txt'],'a');
fprintf(fid,[num2str(flux) ' = ' filename ' - ' inputstring '\n'])
fclose(fid);
