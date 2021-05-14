function [HP] = H_Perspective(PW)
    %% Projects the homogenous coordinates of any degree
    % Input:  
    %         Homogenous coordinates PW
    % Return: 
    %         The homogenous coordinates projected for a hyperplane
    if(PW(end) ~= 0)
        HP = (PW(1:end-1))./PW(end);
    else
        HP = PW(1:end-1);
    end
end