function [ tokens ] = strtokenize( str, delim )
%STRTOKENIZE Summary of this function goes here
%   Detailed explanation goes here

    tokens = cell(0);
    remain = str;

    while true
        [token remain] = strtok(remain, delim);
        if (isempty(token))
            break;
        end;
        tokens(end+1) = {token};
    end;

end

