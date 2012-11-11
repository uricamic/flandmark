function [ w, stat ] = pp_pbmrm( data, risk, lambda, opts )
%PARBMRM Summary of this function goes here
%   Detailed explanation goes here
% 
% 03-01-12, timing statistics: rtime, subgradtime, qptime, hesstime, wtime
% 12-01-12, save euclidean distance of w_t and w_{t-1}

    if nargin<4, opts = struct(); end;
    if ~isfield(opts, 'TolRel'), opts.TolRel = 1e-3; end;
    if ~isfield(opts, 'TolAbs'), opts.TolAbs = 0; end;
    if ~isfield(opts, 'MaxIter'), opts.MaxIter = Inf; end;
    if ~isfield(opts, 'verb'), opts.verb = 0; end;
    if ~isfield(opts,'MaxMemory'), opts.MaxMemory = Inf; end
    if ~isfield(opts,'CleanBuffer'), opts.CleanBuffer = inf; end
%    if ~isfield(opts,'BufSize'), o         0pts.BufSize = 100; end
    
    if ~isfield(opts,'K'), opts.K = 0.4; end
    if ~isfield(opts, 'CleanICP'), opts.CleanICP = false; end;
    if ~isfield(opts, 'CleanAfter'), opts.CleanAfter = 10; end;
    if ~isfield(opts, 'W'), opts.W = []; end;
    if ~isfield(opts, 'Tmax'), opts.Tmax = 100; end;

    if ~iscell(data), data = {data}; end;
    nThreads = length(data);

    rtime = 0; 
    subgradtime = 0;
    qptime = 0; 
    hesstime = 0; 
    wtime = 0;
    
    clck0_ = clock;
    
    K = opts.K;
    alpha = 0;
    w = opts.W;
    Tmax = opts.Tmax;
    tuneAlpha = true;
    alphaChanged = false;
    
    % constants
    exitflag = -1;
    
    % get space dimension
    [R, SG] = risk(data{1});
    wDim  = length(SG);
   
    t = 0;
    Rt = zeros(1, nThreads);
    SGt = zeros(wDim, nThreads);
    
    %% Risk and subgradient computation
    if (isempty(w))
        w = zeros(wDim, 1, 'double');
    end;
    
    for p = 1:nThreads
        [Rt(p), SGt(:,p), ldata, risktim] = risk(data{p}, w);
        data{p} = ldata;
        rtime = rtime + risktim.rtime;
        subgradtime = subgradtime + risktim.subgradtime;
    end;
    
    prevw = w;
    prevw2 = prevw;
    
    hR = [];
    hSG = [];
    I = [];
    f = [];
    beta = [];
    nCP = 0;
        
    CPcounter = [];
    
    hist_Fd(1) = -inf;
    hist_Fp(1) = sum(Rt) + 0.5*lambda*norm(w)^2;
    hist_R(1) = R;
    clck1_ = clock;
    hist_runtime(1) = etime(clck1_, clck0_);
    
    % timing stats
    hist_rtime(1) = rtime;
    hist_subgradtime(1) = subgradtime;
    hist_qptime(1) = qptime;
    hist_hesstime(1) = hesstime;
    hist_wtime(1) = wtime;
    
    % wdists stats
    hist_wdists(1) = norm(prevw - w);
    
    fprintf('%4d: Starting BRMR in %d threads\n', t, nThreads);
    
    while exitflag == -1;
        %tic;
        clck0 = clock;    
        t = t + 1;
%         T = Tmax - t;
        T = Tmax;
        %% Hessian computation
        clock_s = clock;
        nCP = nCP + nThreads;
        hR = [hR Rt'];
        for i = 1:nThreads
            hSG{end+1} = SGt(:,i);
        end;
        f = [f, (w'*SGt - Rt)];
        if t == 1
            H = (SGt'*SGt);
        else
            tmp = zeros(nThreads, nCP);
            for i = 1:nCP
                tmp(:,i) = SGt'*hSG{i};
            end;
            H = [H tmp(:,1:end-nThreads)'];
            H = [H; tmp];
        end;
        clock_f = clock;
        hesstime = etime(clock_f, clock_s);
        
        %% QP solver
        % init temp variables
        beta_good = []; alpha_good = []; stat_good = []; I_good = [];
        % alpha setting cycle
        flag = 1;
        qptime = 0;
        beta_start = beta;
        I_start = I;
        qp_cnt = 0; % counter for QP-solver callss
        if (tuneAlpha)
%             prevw = prevw2;
            alpha_start = alpha; alpha = 0;
            %Solve QP for alpha = 0, if norm(wt-w) < K, then accept this solution, otherwise, set alpha back to alpha_start
            clock_s = clock;
            I = [I_start 1:nThreads];
            beta = [beta_start; zeros(nThreads, 1)];
            [beta, stat] = libqp_splx(H*(1/(lambda + 2*alpha)), f' - (((2*alpha)/(lambda + 2*alpha))*cell2mat(hSG)'*prevw), ones(1, nThreads), I, ones(1, nThreads), beta);
            clock_f = clock;
            qp_cnt = qp_cnt + 1;
            qptime = qptime + etime(clock_f, clock_s);
            Fd_alpha0 = -stat.QP;
            % reconstruct wt and check if norm(wt+1 - wt)^2 < K
            wt = (2*alpha*prevw)/(lambda + 2*alpha);
            for i = 1:nCP
                wt = wt - (hSG{i}*beta(i))/(lambda + 2*alpha);
            end;
            if (norm(prevw - wt) <= K)
                flag = 0;
                if (alpha ~= alpha_start)
                    alphaChanged = true;
                end;
            else
                alpha = alpha_start;
            end;
            while (flag)
                clock_s = clock;
                I = [I_start 1:nThreads];
                beta = [beta_start; zeros(nThreads, 1)];
                [beta, stat] = libqp_splx(H*(1/(lambda + 2*alpha)), f' - (((2*alpha)/(lambda + 2*alpha))*cell2mat(hSG)'*prevw), ones(1, nThreads), I, ones(1, nThreads), beta);
                clock_f = clock;
                qp_cnt = qp_cnt + 1;
                qptime = qptime + etime(clock_f, clock_s);
                % reconstruct wt and check if norm(wt+1 - wt)^2 < K
                wt = (2*alpha*prevw)/(lambda + 2*alpha);
                for i = 1:nCP
                    wt = wt - (hSG{i}*beta(i))/(lambda + 2*alpha);
                end;
                if (norm(prevw - wt) > K)
                    % if there exist a record of some good solution (i.e. we are adjusting alpha by division by 2)
                    if (~isempty(beta_good) && ~isempty(alpha_good) && ~isempty(stat_good) && ~isempty(I_good))
                        beta = beta_good;
                        alpha = alpha_good;
                        stat = stat_good;
                        I = I_good;
                        flag = 0;
                    else
                        if (alpha == 0)
                            alpha = 1;
                            alphaChanged = true;
                        else
                            alpha = alpha * 2;
                            alphaChanged = true;
                        end;
                    end;
                else 
                    if (alpha ~= 0)
                        % keep good solution and try for alpha /= 2 (or alpha = 0 if previous alpha is 1)
                        beta_good = beta; alpha_good = alpha; stat_good = stat; I_good = I;
                        if (alpha ~= 1)
                            alpha = alpha/2;
                            alphaChanged = true;
                        else 
                            alpha = 0;
                            alphaChanged = true;
                        end;
                    else 
                        flag = 0;
                    end;
                end;
            end;
        else 
%             % tune alpha only if || w_t-1 - w_t || > K
%             alphaChanged = false;
%             flag1 = true;
%             while(flag1)
%                 clock_s = clock;
%                 I = [I_start 1:nThreads];
%                 beta = [beta_start; zeros(nThreads, 1)];
%                 [beta, stat] = libqp_splx(H*(1/(lambda + 2*alpha)), f' - (((2*alpha)/(lambda + 2*alpha))*cell2mat(hSG)'*prevw), ones(1, nThreads), I, ones(1, nThreads), beta);
%                 clock_f = clock;
%                 qp_cnt = qp_cnt + 1;
%                 qptime = qptime + etime(clock_f, clock_s);
%                 % reconstruct wt and check if norm(wt+1 - wt)^2 < K
%                 wt = (2*alpha*prevw)/(lambda + 2*alpha);
%                 for i = 1:nCP
%                     wt = wt - (hSG{i}*beta(i))/(lambda + 2*alpha);
%                 end;
%                 if (norm(prevw - wt) > K)
%                     if (alpha == 0)
%                         alpha = 1;
%                         alphaChanged = true;
%                     else
%                         alpha = alpha * 2;
%                         alphaChanged = true;
%                     end;
%                 else 
%                     flag1 = false;
%                 end;
%             end;
            % relax the condition || w_t-1 - w_t || < K and solve just with current alpha
            alphaChanged = false;
            clock_s = clock;
            I = [I_start 1:nThreads];
            beta = [beta_start; zeros(nThreads, 1)];
            [beta, stat] = libqp_splx(H*(1/(lambda + 2*alpha)), f' - (((2*alpha)/(lambda + 2*alpha))*cell2mat(hSG)'*prevw), ones(1, nThreads), I, ones(1, nThreads), beta);
            clock_f = clock;
            qp_cnt = qp_cnt + 1;
            qptime = qptime + etime(clock_f, clock_s);
        end;
        % update CPCounter - add one unused CPs and reset used CPs
        CPcounter = [CPcounter+1; ones(nThreads, 1)];
        CPcounter(beta > 0) = 0;
        %% W update
        clock_s = clock;
        w = (2*alpha*prevw)/(lambda + 2*alpha);
        for i = 1:nCP
            w = w - (hSG{i}*beta(i))/(lambda + 2*alpha);
        end;
        clock_f = clock;
        wtime = etime(clock_f, clock_s);
        %% Risk and subgradient computation
        rtime = 0; subgradtime = 0;
        Rt = zeros(1, nThreads);
        SGt = zeros(wDim, nThreads);
        for p = 1:nThreads
            [Rt(p) SGt(:,p) ldata, risktim] = risk(data{p}, w);
            data{p} = ldata;
            rtime = rtime + risktim.rtime;
            subgradtime = subgradtime + risktim.subgradtime;
        end;
        %% Compute Fp and Fd
        Fp = 0.5*lambda*norm(w)^2 + sum(Rt) + alpha*norm(w-prevw)^2;
        Fd = -stat.QP + ((alpha*lambda)/(lambda + 2*alpha))*(prevw'*prevw);
        %% Gamma and tuneAlpha flag
        if (alphaChanged)
%             epsilon = 1 - Fd/Fp;
%             gamma = (hist_Fp(t)*(1-epsilon) - Fd_alpha0)/(T*(1-epsilon));
% % % % % % % % % 
%             epsilon = 1 - Fp/hist_Fp(t);
%             gamma = (hist_Fp(t)*(1-epsilon) - Fd_alpha0)/(T*(1-epsilon));
% % % % % % % % %
            epsilon = 1 - Fd/Fp;
            gamma = (hist_Fp(t)*(1-epsilon) - Fd_alpha0)/(T*(1-epsilon));
        end;
        
        if ((hist_Fp(t) - Fp) <= gamma)
            tuneAlpha = true;
        else 
            tuneAlpha = false;
        end;
        
        %% Stopping condition - set only if alpha == 0 or nIter > options.MaxIter
        if (alpha == 0 || t >= opts.MaxIter)
            if Fp-Fd<= opts.TolRel*abs(Fp)
                exitflag= 1;
            elseif Fp-Fd <= opts.TolAbs
                exitflag= 2;    
            elseif t >= opts.MaxIter
                exitflag= 0;
            end
        end;

        tmp = whos('H');
        mem = tmp.bytes;
        tmp = whos('hSG');
        mem = mem + tmp.bytes;
        tmp = whos('hR');
        mem = mem + tmp.bytes;
        
        %hist_runtime(t+1) = toc;
        clck1 = clock;
        hist_runtime(t+1) = etime(clck1, clck0);
        hist_Fp(t+1) = Fp;
        hist_Fd(t+1) = Fd;
        hist_R(t+1) = sum(Rt);
        
        % timing stats
        hist_rtime(t+1) = rtime;
        hist_subgradtime(t+1) = subgradtime;
        hist_qptime(t+1) = qptime;
        hist_hesstime(t+1) = hesstime;
        hist_wtime(t+1) = wtime;
        % wdists stats
        wdist = norm(prevw - w);
        hist_wdists(t+1) = wdist;
        
        prevw2 = prevw;
        prevw = w;
        
        nzA = length(find(beta > 0));
        if mod(t,opts.verb) == 0 || exitflag ~= -1
            if (alpha ~= 0)
            fprintf('%4d: tim=%.3f, Fp=%.2f, Fd=%.2f, (Fp-Fd)=%.3f, (Fp-Fd)/Fp=%f, R=%.2f, nCP=%d, nzA=%d, mem=%dMB, wdist=%.3f, alpha=%d, qp_cnt=%d, gamma=%.3f, tuneAlpha=%d\n', ...
                t, hist_runtime(t+1), Fp, Fd, Fp-Fd,(Fp-Fd)/Fp, sum(Rt), nCP, nzA, ceil(mem/1024/1024), wdist, alpha, qp_cnt, gamma, tuneAlpha);
            else
            fprintf('%4d: tim=%.3f, Fp=%.2f, Fd=%.2f, (Fp-Fd)=%.3f, (Fp-Fd)/Fp=%f, R=%.2f, nCP=%d, nzA=%d, mem=%dMB, wdist=%.3f, qp_cnt=%d, gamma=%.3f, tuneAlpha=%d\n', ...
                t, hist_runtime(t+1), Fp, Fd, Fp-Fd,(Fp-Fd)/Fp, sum(Rt), nCP, nzA, ceil(mem/1024/1024), wdist, qp_cnt, gamma, tuneAlpha);
            end;
        end;

        if opts.MaxMemory < (mem/1024/1024)
            exitflag = -2;
        end;
        
        if mod(t,opts.CleanBuffer) == 0 && nzA < t*nThreads
            old_nCP = nCP;
            idx = find(beta > 0);
            nCP = length(idx);

            beta = beta(idx);
            H(1:nCP,1:nCP) = H(idx,idx);
            hSG(:,1:nCP) = hSG(:,idx);
            hR(1:nCP) = hR(idx);

            fprintf('done. (old_nCP=%d, new_nCP=%d) \n',old_nCP, nCP);
        end
        
        % inactiveCP removal
        if opts.CleanICP
            iCP = find(CPcounter>=opts.CleanAfter);
            I(iCP) = [];
            beta(iCP) = [];
            CPcounter(iCP) = [];
            f(iCP) = [];
            hSG(iCP) = [];
            nCP = nCP - length(iCP);
            H(:, iCP) = []; H(iCP, :) = [];
        end;
        
    end;

    stat = [];
    stat.Fp = Fp;
    stat.Fd = Fd;
    stat.nIter = t;
    stat.hist.Fp = hist_Fp(1:t+1);
    stat.hist.Fd = hist_Fd(1:t+1);
    stat.hist.R = hist_R(1:t+1);
    stat.hist.runtime = hist_runtime(1:t+1);
    
    % timing stats
    stat.timings.r = hist_rtime(1:t+1);
    stat.timings.subgrad = hist_subgradtime(1:t+1);
    stat.timings.hessian = hist_hesstime(1:t+1);
    stat.timings.qp = hist_qptime(1:t+1);
    stat.timings.w = hist_wtime(1:t+1);
    stat.timings.overall = hist_runtime(1:t+1);
    
    % euclidean distance of consecutive w_{t-1} and w_t
    stat.hist.wdist = hist_wdists(1:t+1);
    
end

