    function [state,options,changed,str] = savebestineachgen(options,state,flag)
        file_name = 'SaveBest.mat';                                     % Name File
        if strcmp(flag,'init')
            Var = state.Population;
            save(file_name, 'Var')                                      % Write ‘Best Individual’ To File
        elseif strcmp(flag,'iter')
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);
            previous = load('SaveBest.mat');
            Var = [previous.Var; bestx];                                % Read Previous Results, Append New Value
            save(file_name, 'Var')                                      % Write ‘Best Individual’ To File
        end
        changed = true;                                                 % Necessary For Cide, Use  App
    end