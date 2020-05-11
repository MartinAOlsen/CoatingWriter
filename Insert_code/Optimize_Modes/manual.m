
options_home.plot=1;
options_home.seed=1;
options_home.populationsize=16;
options_home.weightrange = 4 ;
options_home.maxIter = 50 ;
options_home.optimizefile = [instrument_name '_optimize.instr'];
options_home.analyzefile = [instrument_name '_analyze.instr'];

options={options_home options_cluster};

options{select}


Best = CW_optimizer(p,options{select},{'intensity','price'},[instrument_name '_optimize.instr'])
