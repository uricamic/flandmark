function [ PsiGS0, PsiGS1, PsiGS2 ] = getDisplacements( options )
%GETDISPLACEMENTS Computes 
%   Detailed explanation goes here
% 
% 10-08-10 Michal Uricar
% 11-07-11 Michal Uricar, corners dataset modification
% 04-08-11 Michal Uricar, added nose component (M = 8)
    
    %% S0
    % get size of S0 bbox
    S0_siz = [options.S(4, 1) - options.S(2, 1) + 1; options.S(3, 1) - options.S(1, 1) + 1];
    x0y0 = zeros(2, S0_siz(1)*S0_siz(2));
    S0_ofset = options.S(1:2, 1);
    
    % get all coordinates of S0
    for i = 1 : S0_siz(2)
        x0y0(1, (i-1)*S0_siz(1)+1:(i*S0_siz(1))) = i-1+S0_ofset(1);
        x0y0(2, (i-1)*S0_siz(1)+1:(i*S0_siz(1))) = S0_ofset(2):S0_ofset(2)+S0_siz(1)-1;
    end;
    
    q0_length = size(x0y0, 2);
    PsiGS0 = cell(q0_length, options.M - 3);
    
    for q = 1 : q0_length
        % fix S0 positoin
        S0 = x0y0(:, q);
        %for comp = 2 : options.M - 2
        % nose component displacement (i.e. g(s0, s7))
        for comp = [2:options.M-3 options.M]
            % get size of Si
            Si_siz = [options.S(4, comp) - options.S(2, comp) + 1; options.S(3, comp) - options.S(1, comp) + 1];
            dxdy = zeros(4, Si_siz(1)*Si_siz(2));
            Si_ofset = options.S(1:2, comp);
            
            % get all coordinates of Si
            for i = 1 : Si_siz(2)
                dxdy(1, (i-1)*Si_siz(1)+1:(i*Si_siz(1))) = i-1+Si_ofset(1) - S0(1);
                dxdy(2, (i-1)*Si_siz(1)+1:(i*Si_siz(1))) = Si_ofset(2)-S0(2):Si_ofset(2)+Si_siz(1)-1-S0(2);
            end;
            dxdy(3, :) = dxdy(1, :).^2;
            dxdy(4, :) = dxdy(2, :).^2;
            if (comp == options.M)
                PsiGS0{q, 5} = int32(dxdy);
            else 
                PsiGS0{q, comp - 1} = int32(dxdy);
            end;
        end;
    end;
      
    %% S1
    % get size of S1 bbox
    S1_siz = [options.S(4, 2) - options.S(2, 2) + 1; options.S(3, 2) - options.S(1, 2) + 1];
    x1y1 = zeros(2, S1_siz(1)*S1_siz(2));
    S1_ofset = options.S(1:2, 2);
    
    % get all coordinates of S1
    for i = 1 : S1_siz(2)
        x1y1(1, (i-1)*S1_siz(1)+1:(i*S1_siz(1))) = i-1+S1_ofset(1);
        x1y1(2, (i-1)*S1_siz(1)+1:(i*S1_siz(1))) = S1_ofset(2):S1_ofset(2)+S1_siz(1)-1;
    end;
    
    q1_length = size(x1y1, 2);
    PsiGS1 = cell(q1_length, 1);
    
    for q = 1 : q1_length
        % fix S1 position
        S1 = x1y1(:, q);
        %for comp = 2 : options.M - 2
        for comp = 6
            % get size of Si
            Si_siz = [options.S(4, comp) - options.S(2, comp) + 1; options.S(3, comp) - options.S(1, comp) + 1];
            dxdy1 = zeros(4, Si_siz(1)*Si_siz(2));
            Si_ofset = options.S(1:2, comp);
            
            % get all coordinates of Si
            for i = 1 : Si_siz(2)
                dxdy1(1, (i-1)*Si_siz(1)+1:(i*Si_siz(1))) = i-1+Si_ofset(1) - S1(1);
                dxdy1(2, (i-1)*Si_siz(1)+1:(i*Si_siz(1))) = Si_ofset(2)-S1(2):Si_ofset(2)+Si_siz(1)-1-S1(2);
            end;
            dxdy1(3, :) = dxdy1(1, :).^2;
            dxdy1(4, :) = dxdy1(2, :).^2;
%             PsiGS1{q, comp - 1} = dxdy1;
            PsiGS1{q, 1} = int32(dxdy1);
        end;
    end;
    
    %% S2
    % get size of S2 bbox
    S2_siz = [options.S(4, 3) - options.S(2, 3) + 1; options.S(3, 3) - options.S(1, 3) + 1];
    x2y2 = zeros(2, S2_siz(1)*S2_siz(2));
    S2_ofset = options.S(1:2, 3);
    
    % get all coordinates of S1
    for i = 1 : S2_siz(2)
        x2y2(1, (i-1)*S2_siz(1)+1:(i*S2_siz(1))) = i-1+S2_ofset(1);
        x2y2(2, (i-1)*S2_siz(1)+1:(i*S2_siz(1))) = S2_ofset(2):S2_ofset(2)+S2_siz(1)-1;
    end;
    
    q2_length = size(x2y2, 2);
    PsiGS2 = cell(q2_length, 1);
    
    for q = 1 : q2_length
        % fix S2 positoin
        S2 = x2y2(:, q);
        %for comp = 2 : options.M - 2
        for comp = 7
            % get size of Si
            Si_siz = [options.S(4, comp) - options.S(2, comp) + 1; options.S(3, comp) - options.S(1, comp) + 1];
            dxdy2 = zeros(4, Si_siz(1)*Si_siz(2));
            Si_ofset = options.S(1:2, comp);
            
            % get all coordinates of Si
            for i = 1 : Si_siz(2)
                dxdy2(1, (i-1)*Si_siz(1)+1:(i*Si_siz(1))) = i-1+Si_ofset(1) - S2(1);
                dxdy2(2, (i-1)*Si_siz(1)+1:(i*Si_siz(1))) = Si_ofset(2)-S2(2):Si_ofset(2)+Si_siz(1)-1-S2(2);
            end;
            dxdy2(3, :) = dxdy2(1, :).^2;
            dxdy2(4, :) = dxdy2(2, :).^2;
%             PsiGS2{q, comp - 1} = dxdy2;
            PsiGS2{q, 1} = int32(dxdy2);
        end;
    end;
    
end
