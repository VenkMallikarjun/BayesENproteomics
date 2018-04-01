%% Empirical Bayes contrast function for comparing different treatments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output -
%       contrasted - cell array containg ProteinOutput with the specified
%       comparison.
%
%%  Input - 
%       ProteinOutput - Protein/PTM/Pathway-level output from
%       BayesENproteomics.m 
%       (ProteinOutput.Abds/ProteinOutput.PTMs/PathwayOutput.Abds,
%       respectively).
%
%       ctrl - Scalar value. Index of group that all others are to be
%       compared to. E.g. to compare to 2nd group listed in ProteinOutput,
%       ctrl = 2.
%
%       PTM - boolean (true or false). Denotes whether comparison is being
%       performed on a PTM output table (true) or not (false).
%
%       GroupNum - Scalar value. Denotes how many treatments are in the
%       dataset in total.
%
%       d0s0 - Optional nx2 double matrix. Where n = GroupNum. Contains
%       shape and scale parameters for chi^2 distrution fit to standard
%       errors for Empirical Bayes comparison. Defaults to value provided
%       in ProteinOutput{1,1}.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [contrasted] = DataProcessEBContrasts(ProteinOutput, ctrl, PTM, GroupNum, d0s0)
if nargin < 5, d0s0 = ProteinOutput{1,1}; end

isPTM = zeros(3,1);
switch PTM
    case true, 
        isPTM(1) = 6;
        isPTM(2) = 4;
        isPTM(3) = 0;
        DoFs = cell2mat(ProteinOutput(3:end,4));
    case false, 
        isPTM(1) = 3;
        isPTM(2) = 0;
        isPTM(3) = 3;
        DoFs = cell2mat(ProteinOutput(3:end,end-1));
end

contrasted = ProteinOutput;
pvals = ones(size(contrasted,1)-2,GroupNum);
for i = 0:1:GroupNum-1
    SEM = cell2mat(ProteinOutput(3:end,isPTM(1) + isPTM(2) + GroupNum + i));
    FC = abs(cell2mat(ProteinOutput(3:end,isPTM(1)+i))-cell2mat(ProteinOutput(3:end,isPTM(1)+ctrl-1)));
    contrasted(3:end,isPTM(1)+i) = num2cell(cell2mat(ProteinOutput(3:end,isPTM(1)+i))-cell2mat(ProteinOutput(3:end,isPTM(1)+ctrl-1)));
    t = FC./SEM;     
    pvals(:,i+1) = min(1,(1 - tcdf(t, DoFs + d0s0(i+1,1) - 1)) * 2 + 1e-15); %2-sided t-test, p-values cannot be 0
    contrasted(3:end,isPTM(1) + isPTM(2) + GroupNum*2 + i) = num2cell(pvals(:,i+1));
    contrasted(3:end,isPTM(1) + isPTM(2) + GroupNum*3 + isPTM(3) + i) = num2cell(mafdr(pvals(:,i+1), 'BHFDR',true));
end


end
