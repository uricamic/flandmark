% Testing LIBQP_SPLX solver.
%
% It solves some benchmark problems and compares results with 
% reference solutions stored in file. 
%
% If switched on (QUADPROG=1) then it does the same with 
% Matlab QUADPROG for comparison.
%
% LIBQP_SPLX solves an instance of Quadratic Proramming task defined as:
%
%   min 0.5*x'*H*x + f'*x  
%    x
%
% subject to:   sum(x(find(I==k))) == b(k)  for all k in find(S==0)
%               sum(x(find(I==k))) <= b(k)  for all k in find(S==1)
%                             x(i) >= 0 for all i=1:n
%


help libqp_splx_test;

DataFolder = '../data/';
QPTask = {'splx_qp_random100x10','splx_qp_random300x30'};
%QPTask = {'splx_qp_random100x10','splx_qp_random300x30','splx_qp_delta1'};

QUADPROG = 1;  
%QUADPROG = 0;

SAVE_AS_REFERENCE = 0;  % if 1 then the results are saved and later they serve
                        % as reference solutions

for i=1:length(QPTask)

    qptask_fname = [DataFolder QPTask{i} '.mat'];
    libqp_ref_solution_fname = [DataFolder QPTask{i} '_libqp_sol.mat'];
    qprog_ref_solution_fname = [DataFolder QPTask{i} '_quadprog_sol.mat'];
    
    load(qptask_fname,'H','f','b','I','S');

    fprintf('\n--------------------------------------------------\n');
    fprintf('QP task: %s\n',qptask_fname);
    fprintf('n: %d, max(I): %d\n', length(I),max(I));
    fprintf('--------------------------------------------------\n');
        
    % CALL LIBQP_SPLX with dense H
    t0=clock;
    [libqp.x,libqp.stat] = libqp_splx(H,f,b,I,S);
        
    libqp.time = etime(clock, t0);
    libqp.nnz = nnz(libqp.x);

    % DISPLAY RESULTS
    fprintf(' (dense H)              prim_val      time      nnz    iter       dual_gap\n');
    fprintf('LIBQP_SPLX       %.10f     %.3f  %7d   %5d   %.10f\n', ...
            libqp.stat.QP, libqp.time, libqp.nnz, libqp.stat.nIter, libqp.stat.QP-libqp.stat.QD);

    % LOAD/SAVE reference solution
    if SAVE_AS_REFERENCE,
      if ~exist(libqp_ref_solution_fname)
          save(libqp_ref_solution_fname,'libqp');
      else
          error(sprintf('Reference solution %s already exists. Erase it manually and re-run the test.',...
                        libqp_ref_solution_fname));
      end
    else 
      if ~exist(libqp_ref_solution_fname)
          error(sprintf('Reference solution %s does not exist.',...
                        libqp_ref_solution_fname));
      else
          ref=load(libqp_ref_solution_fname);
          fprintf('LIBQP REF.       %.10f     %.3f  %7d   %5d   %.10f   %s\n', ...
                  ref.libqp.stat.QP, ref.libqp.time, ref.libqp.nnz, ...
                  ref.libqp.stat.nIter,...
                  ref.libqp.stat.QP - ref.libqp.stat.QD,libqp_ref_solution_fname);
          sol_dif = max(abs(ref.libqp.x-libqp.x));
          fprintf('Comparing to reference: max(abs(ref.x-x)) = %.10f.', sol_dif);
          if sol_dif == 0
              fprintf(' OK.\n\n');
          else
              fprintf(' THIS IS NOT GOOD ... BUT DON''T PANIC.\n\n');
          end
      end
    end 
    

    % CALL LIBQP_SPLX with sparse H
    t0=clock;
    [libqp.x,libqp.stat] = libqp_splx(sparse(H),f,b,I,S);
        
    libqp.time = etime(clock, t0);
    libqp.nnz = nnz(libqp.x);

    % DISPLAY RESULTS
    fprintf(' (sparse H)             prim_val      time      nnz    iter       dual_gap\n');
    fprintf('LIBQP_SPLX       %.10f     %.3f  %7d   %5d   %.10f\n', ...
            libqp.stat.QP, libqp.time, libqp.nnz, libqp.stat.nIter, libqp.stat.QP-libqp.stat.QD);

    % LOAD/SAVE reference solution
    if SAVE_AS_REFERENCE,
      if ~exist(libqp_ref_solution_fname)
          save(libqp_ref_solution_fname,'libqp');
      else
          error(sprintf('Reference solution %s already exists. Erase it manually and re-run the test.',...
                        libqp_ref_solution_fname));
      end
    else 
      if ~exist(libqp_ref_solution_fname)
          error(sprintf('Reference solution %s does not exist.',...
                        libqp_ref_solution_fname));
      else
          ref=load(libqp_ref_solution_fname);
          fprintf('LIBQP REF.       %.10f     %.3f  %7d   %5d   %.10f   %s\n', ...
                  ref.libqp.stat.QP, ref.libqp.time, ref.libqp.nnz, ...
                  ref.libqp.stat.nIter,...
                  ref.libqp.stat.QP - ref.libqp.stat.QD,libqp_ref_solution_fname);
          sol_dif = max(abs(ref.libqp.x-libqp.x));
          fprintf('Comparing to reference: max(abs(ref.x-x)) = %.10f.', sol_dif);
          if sol_dif == 0
              fprintf(' OK.\n\n');
          else
              fprintf(' THIS IS NOT GOOD ... BUT DON''T PANIC.\n\n');
          end
      end
    end 

    if QUADPROG
        % SETUP TASK 
        nEqu = sum(S==0);
        nInEqu = sum(S==1);
        nVar = length(I);
                
        Aeq = zeros(nEqu,nVar);
        beq = zeros(nEqu,1);
        cnt = 0;
        for i=find(S(:)'==0)
          cnt = cnt + 1;
          Aeq(cnt,find(I==i)) = 1;
          beq(cnt) = b(i);
        end

        Aineq = zeros(nInEqu,nVar);
        bineq = zeros(nInEqu,1);
        cnt = 0;
        for i=find(S(:)'==1)
          cnt = cnt + 1;
          Aineq(cnt,find(I==i)) = 1;
          bineq(cnt) = b(i);
        end

        % CALL QUADPROG
        opt=optimset('Display','off','Diagnostics','off');
        t0=clock;
        warning off;
        [qprog.x,qprog.QP,qprog.exitflag,qprog.output] = quadprog(H,f,Aineq,bineq,Aeq,beq,zeros(nVar,1),[],[],opt);
        warning on;
        qprog.time = etime(clock,t0);
        qprog.nnz = nnz(qprog.x);

        fprintf('QUADPROG         %.10f     %.3f  %7d\n', qprog.QP, qprog.time, qprog.nnz); 

        % LOAD/SAVE reference solution
        if SAVE_AS_REFERENCE,
            if ~exist(qprog_ref_solution_fname)
                save(qprog_ref_solution_fname,'qprog');
            else
                error(sprintf('Reference solution %s already exists. Erase it manually and re-run the test.',...
                              qprog_ref_solution_fname));
            end
        else
            if ~exist(qprog_ref_solution_fname)
                error(sprintf('Reference solution %s does not exist.',...
                              qprog_ref_solution_fname));
            else
                ref=load(qprog_ref_solution_fname);
                fprintf('QUADPROG REF.    %.10f     %.3f  %7d                  %s\n', ...
                        ref.qprog.QP, ref.qprog.time, ref.qprog.nnz, qprog_ref_solution_fname); 
                
                sol_dif = max(abs(ref.qprog.x-qprog.x));
                fprintf('Comparing to reference: max(abs(ref.x-x)) = %.10f.', sol_dif);
                if sol_dif == 0
                    fprintf(' OK.\n\n');
                else
                    fprintf(' THIS IS NOT GOOD ... BUT DON''T PANIC.\n\n');
                end
            end
        end 
        
        fprintf('Comparison between LIBQP_SSVM and QUADPROG: max(abs(libqp.x-quadprog.x)) = %.10f\n', ...
                max(abs(libqp.x - qprog.x)));
        fprintf('                                            libqp.prim_val-quadprog.prim_val = %.15f\n\n',...
                libqp.stat.QP-qprog.QP);
    end        
end    

% EOF