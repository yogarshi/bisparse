function fasta_biling_solver(mu_en, mu_fr, mu_xling, output_size, tol, iters, fr_vecs, en_vecs,  mat_file)
%path('./FASTA/solvers/', path)

% Load the english and french vectors
S_EN = importdata(en_vecs,' ');
S_FR = importdata(fr_vecs, ' ');

[v_en, dim1] = size(S_EN);
[v_fr, dim2] = size(S_FR);


% The two vector spaces should have the same dimensions
%assert(dim1 == dim2, 'The dimensions of the 2 vector spaces do not match')

%Randomly init X and Y for both langauges
X_EN0 = rand([v_en,output_size]);
Y_EN0 = rand([dim1,output_size]);
X_FR0 = rand([v_fr,output_size]);
Y_FR0 = rand([dim2,output_size]);



% Options for FASTA - look up FASTA documentation
opts = [];
opts.maxIters = iters;
opts.tol = tol;  		% Use super strict tolerance
%opts.recordObjective = true; 	% Record the objective function so we can plot it.
opts.verbose=2;			% Verbosity level	
opts.stringHeader='    ';      	% Append a tab to all text output from FISTA. Nicer formatting. 

% Load the alignment stats
load(mat_file);
AL_proper = AL(1:v_en,1:v_fr);
%Call Fasta
[ X_EN_sol, Y_EN_sol, X_FR_sol, Y_FR_sol, outs ] = fasta_biling(S_EN, X_EN0, Y_EN0, S_FR, X_FR0, Y_FR0, AL_proper, mu_en, mu_fr, mu_xling, opts);


%Write to file
dlmwrite(strcat('bisparse_en.txt'),X_EN_sol,'\t');
dlmwrite(strcat('bisparse_fr.txt'),X_FR_sol,'\t');


exit;
