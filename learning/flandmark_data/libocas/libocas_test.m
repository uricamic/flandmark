% This script tests functionality of SVMOCAS and MSVMOCAS solvers.
%
% It runs the solvers on example data and compares results to solutions
% stored in reference files. 
%  

% two-class problem 
BinaryTrnFile = './data/riply_trn.light';

% multi-class problem
MulticlassTrnFile = './data/example4_train.light';

% file to store/load reference solution
ReferenceFile = './data/refernce_solution';

% if 1 save results to reference files else compares the results to the
% reference solutions
CREATE_REFERNCE_FILES = 0;

% Solver options
opt.C = 1;
opt.Method = 1;
opt.TolRel = 0.01;
opt.TolAbs = 0;
opt.QPBound = 0;
opt.BufSize = 2000;
opt.MaxTime = inf;
opt.X0 = 1;
opt.verb = 0;

fprintf('Training binary SVM classifier by SVMOCAS...');
[bin.W,bin.W0,bin.stat] = svmocas(BinaryTrnFile,opt.X0,opt.C,opt.Method,opt.TolRel,...
                             opt.TolAbs,opt.QPBound,opt.BufSize,inf,opt.MaxTime,opt.verb);
fprintf('done.\n');

fprintf('Training multi-class SVM classifier by MSVMOCAS...');
[multi.W,multi.stat] = msvmocas(MulticlassTrnFile,opt.C,opt.Method,opt.TolRel,...
                             opt.TolAbs,opt.QPBound,opt.BufSize,inf,opt.MaxTime,opt.verb);
fprintf('done.\n');


if CREATE_REFERNCE_FILES == 1,    
    fprintf('Saving reference solutions to %s\n', ReferenceFile);
    save(ReferenceFile,'bin','multi');    
else
    ref = load(ReferenceFile);    
        
    test(1).dif = sum(abs(bin.W - ref.bin.W)+abs(bin.W0-ref.bin.W0));
    test(1).name = 'sum(|W-ref.W| + |W0-ref.W0])';
    test(2).dif = abs(bin.stat.Q_P - ref.bin.stat.Q_P);
    test(2).name = 'PrimalVal - ref.PrimalVal   ';
    test(3).dif = abs(bin.stat.Q_D - ref.bin.stat.Q_D);
    test(3).name = 'DualVal - ref.DualVal       ';
    
    fprintf('\nSVMOCAS (solver for binary classifiation problems):\n');
    for i=1:length(test)
        fprintf('   %s = %.20f ... ',test(i).name,test(i).dif);
        if test(i).dif == 0
            fprintf('SOLUTIONS EQUAL - OK\n');
        else
            fprintf('SOLUTION DIFFERS\n');
        end
    end        

    test(1).dif = sum(sum(abs(multi.W - ref.multi.W)));
    test(1).name = 'sum(|W-ref.W|)           ';
    test(2).dif = abs(multi.stat.Q_P - ref.multi.stat.Q_P);
    test(2).name = 'PrimalVal - ref.PrimalVal';
    test(3).dif = abs(multi.stat.Q_D - ref.multi.stat.Q_D);
    test(3).name = 'DualVal - ref.DualVal    ';
    
    fprintf('\nMSVMOCAS: (solver for multi-class problems\n');
    for i=1:length(test)
        fprintf('   %s = %.20f ... ',test(i).name, test(i).dif);
        if test(i).dif == 0
            fprintf('SOLUTIONS EQUAL - OK\n');
        else
            fprintf('SOLUTION DIFFERS\n');
        end
    end        
    
    
end


break;

for i=1:length(DataSets)  
   fprintf('\nDataset: %s\n', DataSets{i});
   
   [exp{i}.W,exp{i}.W0,exp{i}.stat] = svmocas(DataSets{i},opt.X0,opt.C,opt.Method,...
                opt.TolRel,opt.TolAbs,opt.QPBound,opt.BufSize,inf,opt.MaxTime);  
end

fprintf('\n\nRESULTS SUMMARY\n================================\n\n');
for i=1:length(DataSets)
   fprintf('\nDataset: %s\n--------------------------------\n', DataSets{i});
   
   % remove suffix
   sol_fname = DataSets{i};
   idx = findstr(sol_fname,'.');
   sol_fname = sol_fname(1:idx(end)-1);
         
   sol_fname = [sol_fname '_ocas_C' num2str(opt.C) '_solution.mat'];
   if SAVE_AS_REFERENCE,
      if exist(sol_fname)
         fprintf('Solution file %s already exists.\n', sol_fname);
         error('Erase the file or set SAVE_AS_REFERENCE = 0 and run the test again.');
      else
          fprintf('Saving solution to %s ...', sol_fname);
          ref_sol = exp{i};
          ref_opt = opt;
          save(sol_fname,'ref_sol','ref_opt');
          fprintf('done.\n');      
          ref_sol = [];
      end
   else
      if exist(sol_fname)
          load(sol_fname,'ref_sol','ref_opt');          

          fprintf('\nReference solution\n');
          fprintf('file: %s\n', sol_fname);
          fprintf(['settings: C: %f, Method: %d, TolRel: %f, TolAbs: %f, ' ...
                   'QPBound: %f, BufSize: %d, MaxTime: %f, X0: %f \n'], ...
                  ref_opt.C, ref_opt.Method, ref_opt.TolRel, ref_opt.TolAbs, ...
                  ref_opt.QPBound, ref_opt.BufSize, ref_opt.MaxTime, ref_opt.X0);
          fprintf(['solution: QP: %.10f, QD: %.10f, nIter: %d, nCutPlanes: %d, '...
                'exitflag: %d, ocas_time: %f, total_time: %f\n'],...
                  ref_sol.stat.Q_P, ref_sol.stat.Q_D, ref_sol.stat.nIter, ref_sol.stat.nCutPlanes, ...
                  ref_sol.stat.exitflag, ref_sol.stat.ocas_time, ref_sol.stat.total_time);
      end       
   end
   
   fprintf('\nCurrent solution \n');
   fprintf(['settings: C: %f, Method: %d, TolRel: %f, TolAbs: %f, ' ...
                   'QPBound: %f, BufSize: %d, MaxTime: %f, X0: %f \n'], ...
                  opt.C, opt.Method, opt.TolRel, opt.TolAbs, ...
                  opt.QPBound, opt.BufSize, opt.MaxTime, opt.X0);
   fprintf(['solution: QP: %.10f, QD: %.10f, nIter: %d, nCutPlanes: %d, '...
            'exitflag: %d, ocas_time: %f, total_time: %f\n'],...
           exp{i}.stat.Q_P, exp{i}.stat.Q_D, exp{i}.stat.nIter, exp{i}.stat.nCutPlanes, ...
           exp{i}.stat.exitflag, exp{i}.stat.ocas_time, exp{i}.stat.total_time);
   
   if ~isempty(ref_sol)
       opt_dif = any([opt.Method~=ref_opt.Method opt.C~=ref_opt.C opt.TolRel ~= ref_opt.TolRel ...
                      opt.TolAbs~=ref_opt.TolAbs opt.QPBound~=ref_opt.QPBound ...
                      opt.BufSize~=ref_opt.BufSize opt.MaxTime~=ref_opt.MaxTime opt.X0~=ref_opt.X0]);
       sol_dif = max(max(abs(exp{i}.W-ref_sol.W)),abs(exp{i}.W0-ref_sol.W0));
       fprintf('\nSolution difference: max(abs(ref.W-W)) = %f. Solutions are ', sol_dif); 
       if ~sol_dif 
           fprintf('EQUAL.\n'); 
       else
           fprintf('DIFFERENT.\n');
       end
       fprintf('Current and refrence settings are ');
       if ~opt_dif 
           fprintf('EQUAL.\n'); 
       else
           fprintf('DIFFERENT.\n');
       end
       fprintf('Comparison result: ');
       if sol_dif && ~opt_dif ,
           fprintf('THIS IS NOT GOOD ... BUT DON''T PANIC.\n');
       elseif ~sol_dif && ~opt_dif
           fprintf('OK.\n');
       elseif opt_dif
           fprintf('UNDECIDED since options are different.\n');
       end
       
   end
end
