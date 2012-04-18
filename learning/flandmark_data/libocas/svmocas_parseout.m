function [iter,timing,trn_err] = svmocas_parseout(fname)
% SVMOCAS_PARSEOUT parses text output of SVMOCAS solver.
% 
% Synopsis:
%   [iter,timing,trn_err] = svmocas_parseout(fname)
%
% Inputs:
%  fname [string] file which contains OCAS output; e.g. created by running
%
%    ./svmocas -c 10 -b 1 riply_trn.light ocas.model > out.txt
%
% Outputs:
%  iter [1 x nIter (struct)] progress info printed by OCAS solver:
%   .time [1x1] runtime
%   .Q_P [1x1] primal objective value
%   .Q_D [1x1] dual (reduced) objective value
%   .trn_err [1x1] training error
%   .nza [1x1] number of active cutting planes
%   .qpf [1x1] return flag of inner QP solver
%
%  timing [struct] Timing statistics.
%  trn_err [1x1] Training error.
% 
% Example:  
%  iter = svmocas_parseout('out.txt');
%  figure; hold on;
%  plot([iter.time],[iter.Q_P]); 
%  plot([iter.time],[iter.Q_D],'r')
% 


lines = textread(fname,'%s','delimiter','\n','bufsize',100000);

from = strmatch('Starting optimization:',lines,'exact')+1;
to = strmatch('Stopping condition:',lines)-1;

for i=from:to
   cnt = i-from+1;
   l = lines{i};

   b = strfind(l,'tim=');
   e = min(strfind(l(b:end),','));
   iter(cnt).time = str2num(l(b+length('tim='):b+e-2));

   b = strfind(l,'Q_P=');
   e = min(strfind(l(b:end),','));
   iter(cnt).Q_P = str2num(l(b+length('Q_P='):b+e-2));

   b = strfind(l,'Q_D=');
   e = min(strfind(l(b:end),','));
   iter(cnt).Q_D = str2num(l(b+length('Q_D='):b+e-2));

   if i > from
       b = strfind(l,'nza');
       e = min(strfind(l(b:end),','));
       iter(cnt).nza = str2num(l(b+length('nza'):b+e-2));

       b = strfind(l,'err');
       e = min(strfind(l(b:end),','));
       iter(cnt).trn_err = str2num(l(b+length('err'):b+e-2));

       b = strfind(l,'qpf');
       iter(cnt).qpf = str2num(l(b+length('nza'):end));
   else
       iter(cnt).nza = [];
       iter(cnt).trn_err = [];
       iter(cnt).qpf = [];
   end
end

i = strmatch('load_time',lines);
if isempty(i)
    i = strmatch('init_time',lines);
end
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.load_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('qp_solver_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.qp_solver_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('sort_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.sort_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('output_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.output_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('add_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.add_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('w_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.w_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('ocas_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.ocas_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('total_time',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'[s]'));
timing.total_time = str2num(l(b+length(':'):b+e-2));

i = strmatch('Training error:',lines);
l = lines{i};
b = strfind(l,':');
e = min(strfind(l(b:end),'%'));
trn_err = str2num(l(b+length(':'):b+e-2))/100;



return;
% EOF