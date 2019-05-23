%Finds all unique modified sites within all peptides that were used in
%protein quantification. Only calculates fold changes unless using bayesian
%linear regression.
function [RMAList,ProtList] = CreatePTMList(normedpeplist, protabds,...
    normedlength, modlist, w, RA, idcol, uniprottable, upcol, reg, groups)
    
name = '';
switch reg
    case 'bayes'
        v = 3;
    otherwise
        v = 2;
end
for i = RA:RA-1+w
    if isempty(cell2mat(normedpeplist(v,i)))
        normedpeplist{v,i} = name;
    else
        name = normedpeplist{v,i};
    end
end

%Preload arrays to save time
%Groups = normedpeplist(1,RA:1:RA-1+w);
GroupList = groups;
GroupNum = numel(GroupList);
if strcmp(reg, 'bayes') || strcmp(reg,'ME')
    PTMList = cell(normedlength, 6 + GroupNum*4);
    w = GroupNum*4 + 4;
else
    PTMList = cell(normedlength, 6 + GroupNum);
end
    
ii = 1;
%fill mod array
for i = 2:normedlength
    ModFind = strfind(normedpeplist{i,10}, modlist);
    %Make array showing peptides with specified mod.
    if (size(ModFind, 2) > 0)
        PTMList(ii,1) = normedpeplist(i,idcol);  %ACCESSION/GENE ID
        PTMList(ii,2) = normedpeplist(i,10);     %PTMs
        PTMList(ii,3) = normedpeplist(i,9);      %PEPTIDE SEQUENCE
        PTMList(ii,4) = normedpeplist(i,8);      %SCORE
        PTMList(ii,5) = {i};                     %PEPTIDE NUMBER IN DATASET
        PTMList(ii,6) = normedpeplist(i,13);     %USE IN QUANTITATION
        if strcmp(reg, 'bayes') || strcmp(reg,'ME')
            PTMList(ii,7:end) = normedpeplist(i,RA:end-2);
        else PTMList(ii,7:6+w) = normedpeplist(i,RA:end);
        end
        ii = ii + 1;
    end
end
PTMLength = size(PTMList,1);
for index = 2:PTMLength
    if isempty(PTMList{index,3})
        PTMList = PTMList(1:(index - 1), :);
        break;
    end
end   
PTMLength = index - 1;
switch reg
    case 'bayes'; RMAList = cell(PTMLength, 10 + GroupNum*4);
    case 'ME'; RMAList = cell(PTMLength, 10 + GroupNum*4);
    otherwise; RMAList = cell(PTMLength, 9 + GroupNum);
end
ProtList = cell(PTMLength, 7 + GroupNum);
AllProtAbundances = zeros(1,GroupNum);
iiii = 3;
n = 3;

for i = 1:PTMLength
    ModPep = PTMList(i,1);
   
    %Skip over peptides whose relative mod abundances have already been
    %calculated.
    try if isequal(ModPep, RMAList(n-1,1))
            continue
        end
    catch
        continue;
    end
    
    %Find all modified and unmodified peptides corresponding to
    %that protein
    UnModProtFind = strcmp(normedpeplist(:,idcol), char(ModPep));
    
    %Skip over peptides whose proteins have been removed from 
    %protein abundance list.
    if sum(UnModProtFind) == 0
        continue
    end
    
    %Find position of modified protein total abundance values in
    %totalpeplist.
    UnModProtFind = strcmp(protabds(:,1),ModPep);
    
    %Get separate list of abundances for all protein that possess 'mod'.
    if ~strcmp(ProtList(iiii-1,1), ModPep)
        try ProtList(iiii,1) = ModPep;                  %IDENTIFIER
            ProtList(iiii,2) = protabds(UnModProtFind,2); %PEPTIDE NUMBER
            ii = 3:GroupNum+2;
            ProtList(iiii,ii) = protabds(UnModProtFind,ii);
            AllProtAbundances(1,:) = cell2mat(ProtList(iiii,ii));
            DF = protabds(UnModProtFind,end-1);
            iiii = iiii + 1; %Move ProtList to next entry.
        catch
            %warning('Skipping PTM for protein with no identifier')
            continue;
        end
    end
    
    %Get protein sequence, uniprot id, names and go ids associated with
    %ModPep
    ModPepnoDBID = char(ModPep);
    %if idcol == 11
    %    ModPepnoDBID = {ModPepnoDBID(1,4:end)};
    %end
    uniprotposition = strcmp(uniprottable(:,upcol), ModPepnoDBID);
    
    if sum(uniprotposition,1) > 0
        ProtSeq = uniprottable(uniprotposition,11);
        UniProtID = uniprottable(uniprotposition,1);
        ProteinName = uniprottable(uniprotposition,3);
        GOIDs = uniprottable(uniprotposition,9);
    else
        continue;
    end

    %Set up variables
    modfind = zeros(PTMLength,3); %Preallocate modfind
    peptidemods = char(zeros(PTMLength,10));
    ModPos = NaN(length(modfind),3,10);
    ModdedPepsPerProtein = cell(PTMLength,10 + w);
    jj = 1;
    
    %Find all modded peptides belonging to protein specificed in ModPep.
    for kk = 3:1:PTMLength
        if strcmp(PTMList(kk,1),char(ModPep))
            switch reg
                case 'bayes'; ModdedPepsPerProtein(jj,1:6 + GroupNum*4) = PTMList(kk,:);
                case 'ME'; ModdedPepsPerProtein(jj,1:6 + GroupNum*4) = PTMList(kk,:);
                otherwise; ModdedPepsPerProtein(jj,1:(6 + w)) = PTMList(kk,:);
            end
            jj = jj + 1;
        end
    end
    
    %Shorten ModdedPepsPerProtein array to remove empty cells from
    %preloading
    for index = 1:1:PTMLength
        if isempty(ModdedPepsPerProtein{index,1})
            ModdedPepsPerProtein = ModdedPepsPerProtein(1:(index - 1), :);
            break;
        end
    end
    
    MPPPLength = size(ModdedPepsPerProtein,1);
    %Get positions of mods within protein sequence
    for ll = 1:1:MPPPLength
        %positions in text of mods of peptides matching current peptide
        moddedpeps = strfind(ModdedPepsPerProtein{ll,2},modlist);  
        oo = 1:size(moddedpeps,2);
        modfind(ll,oo) = moddedpeps(1,oo);
        
        %List of mods for peptides matching current peptide
        PTMperPeplist = char(ModdedPepsPerProtein(ll,2));
        oo = 1:size(PTMperPeplist,2);
        peptidemods(ll,oo) = PTMperPeplist(1,oo);
        
        for m = 1:size(moddedpeps,2)
            if modfind(1,m) == 0
                continue
            end
            ModPos(ll,1,m) = str2double(peptidemods(ll,modfind(ll,m)-4))...
                * 10;
            ModPos(ll,2,m) = str2double(peptidemods(ll,modfind(ll,m)-3));
            if ~isnan(ModPos(ll,1,m))
                ModPos(ll,3,m) = ModPos(ll,1,m) + ModPos(ll,2,m);
            else ModPos(ll,3,m) = ModPos(ll,2,m);
            end
            
            %Catch for single digit mod numbers
            if isnan(ModPos(ll,2,m))
                ModPos(ll,3,m) = ModPos(ll,1,m) / 10;
            end
        end
        ModPositionInProt = strfind(ProtSeq,...
            char(ModdedPepsPerProtein(ll,3)));
        t = zeros(size(ModPositionInProt,2),size(ModPositionInProt,1));
        for p = 1:size(ModPositionInProt,1)
            t(1:size(ModPositionInProt{p},2),p) = ModPositionInProt{p};
            if isempty(ModPositionInProt{p})
                t(1,p) = NaN;
            end
        end
        p = 1:size(ModPositionInProt,1);
        t = unique(t(t~=0),'stable');
        %Check which PTMs belong to which isoform
        %t_member = zeros(size(t,1),size(p,2)); 
        for ii = p
            t_member(:,ii) = ismember(t,cell2mat(ModPositionInProt(ii)))';
        end
        t_member = max(t_member,[],1);
        PTMPOSITIONS = permute(bsxfun(@plus, ModPos(ll,3,:), t)-1,[2,3,1]);
        PTMPOSITIONS = PTMPOSITIONS(~isnan(PTMPOSITIONS));
        %Get info for correct isoform with PTM
        try ModdedPepsPerProtein(ll,7 + w) = ProtSeq(p(t_member > 0),1);
            ModdedPepsPerProtein(ll,9 + w) = UniProtID(p(t_member > 0),1);
            ModdedPepsPerProtein(ll,10 + w) = ...
                ProteinName(p(t_member > 0),1);
        catch %If several isoforms can contain PTM just pick the first one.
            if ~isempty(ProtSeq(p(t_member > 0),1))
                temp = ProtSeq(p(t_member > 0),1);
                ModdedPepsPerProtein(ll,7 + w) = temp(1,1);
                temp = UniProtID(p(t_member > 0),1);
                ModdedPepsPerProtein(ll,9 + w) = temp(1,1);
                temp = ProteinName(p(t_member > 0),1);
                ModdedPepsPerProtein(ll,10 + w) = temp(1,1);
            end
        end
        if ~isempty(PTMPOSITIONS)
            ModdedPepsPerProtein(ll,8 + w) = {PTMPOSITIONS(1,:)};
        end
        clearvars('t_member');
    end
    
    %Get ranges for mod search whereby abundance values for peptides
    %containing mods at abundances positions are gathered.
    numarray = zeros(MPPPLength,2);
    infoarray = cell(MPPPLength,3);
    ccc = 1;
    for bb = 1:MPPPLength
        if ~isequal(ModdedPepsPerProtein(bb,8 + w),{[]})
            cc = 1:size(cell2mat(ModdedPepsPerProtein(bb,8 + w)),2);
            numarray(ccc,cc) = cell2mat(ModdedPepsPerProtein(bb,8 + w));
            infoarray(ccc,cc) = ...
                {{ModdedPepsPerProtein(bb,9 + w),...%Uniprot ID
                ModdedPepsPerProtein(bb,10 + w),... %Name
                ModdedPepsPerProtein(bb,7 + w)}};  %Sequence
            %infoarray{ccc,4} = bb; %PTM ID - called on later
            ccc = ccc + 1;
        end
    end
    [numarray,ia] = unique(numarray(numarray ~= 0), 'stable');
    infoarray = infoarray(:);
    infoarray = infoarray(~cellfun('isempty',infoarray));
    infoarray = infoarray(ia,:);
   
    for cc = 1:size(numarray,1)
        %In cases where no abundances (i.e. where proteins have been
        %removed from list) -> skip.
        if isempty(numarray)
            continue
        end
        
        %Get abundances for all peptides that have mod in position cc.
        dd = 1;
        ModProtAbundances = cell(MPPPLength,w);
        if strcmp(reg, 'bayes') || strcmp(reg, 'ME')
            ModProtAbundances = cell(MPPPLength,GroupNum);
            ModProtStats = cell(MPPPLength,GroupNum*3);
        end
        peps = ' [';
        score = NaN;
        for bb = 1:MPPPLength 
            switch reg
                case 'bayes', index = cell2mat(ModdedPepsPerProtein(bb,end-2));
                case 'ME', index = cell2mat(ModdedPepsPerProtein(bb,end-2));
                otherwise, index = cell2mat(ModdedPepsPerProtein(bb,8 + w));
            end
            if ismember(numarray(cc,1),index)
                if strcmp(reg, 'bayes') || strcmp(reg, 'ME')
                    ModProtAbundances(dd,:) = ModdedPepsPerProtein(bb,7:6 + GroupNum);
                    ModProtStats(dd,:) = ModdedPepsPerProtein(bb, 7+GroupNum:6+GroupNum*4);
                else
                    ModProtAbundances(dd,:) = ModdedPepsPerProtein(bb,7:6 + w);
                end
                dd = dd + 1;
                
                peps = cat(2,peps,...
                    num2str(cell2mat(ModdedPepsPerProtein(bb,5))),'] ',...
                    char(ModdedPepsPerProtein(bb,3)),'; [');
                
                if isequal({cell2mat(ModdedPepsPerProtein(bb,4))},{'---'})
                    score(end+1,1) = 0;
                else
                    score(end+1,1) = str2double(cell2mat(ModdedPepsPerProtein(bb,4)));
                end
            end
        end
        score = nanmean(score,1);
        
        %Calculate relative mod abundances per sample
        if ~isempty(ModProtAbundances)
            for index = 1:size(ModProtAbundances,1)
                if isempty(ModProtAbundances{index,1})
                    ModProtAbundances = ModProtAbundances(1:(index - 1),:);
                    break;
                end
            end
            Abds = NaN(size(ModProtAbundances));
            for ii = 1:1:size(ModProtAbundances,2)
                Abds(:,ii) = cell2mat(ModProtAbundances(:,ii));
                Abds(:,ii) = sort(Abds(:,ii), 1, 'descend');
            end
            if strcmp(reg, 'bayes') || strcmp(reg, 'ME')
                Stats = nanmedian(cell2mat(ModProtStats),1);
            end
            %Calculate relative mod abundance for the current sample 
            %based on method stipulated by calcmethod. Collapse all
            %peptides with PTM at same position via medians. Take mean
            %accross all replicates and normalise to fold change
            %observed at protein level (ProteinAbundance).                
            Abds = nanmedian(Abds,1);
            if ~strcmp(reg, 'bayes') && ~strcmp(reg, 'ME')
                PTMAbds_final = zeros(1,GroupNum);
                for zz = 1:GroupNum
                    PTMAbds_final(1,zz) = nanmean(Abds(1,strcmp(GroupList,...
                        GroupList(zz))),2);
                end
                ctrl = repmat(PTMAbds_final(:,1), 1, GroupNum);
                PTMAbds_final = PTMAbds_final - ctrl; %Get log2(fold changes)
                %Normalise PTM fold changes to protein abundance fold changes
                PTMAbds_final = PTMAbds_final - AllProtAbundances(1,:);
            else
                PTMAbds_final = Abds;% - AllProtAbundances(1,:);
                if strcmp(reg, 'bayes') || strcmp(reg, 'ME')
                    RMAList(n,10+GroupNum:9+GroupNum+size(Stats,2)) = num2cell(Stats);
                end
            end
              
            RMAList(n,6:5+GroupNum) = num2cell(PTMAbds_final);
                
            num = numarray(cc,1);
            protein = infoarray{cc,1};
            ProtSeqChar = char(protein{1,3});
            title = cat(2,ProtSeqChar(1,num),num2str(num),'; ',...
            char(protein{1,2}));
            RMAList(n,1) = ModPep;   %GENE NAME
            RMAList(n,2) = {peps};   %List of peptides used in relative mod
                                     %abundance calculation.
            RMAList(n,3) = {num}; %Position of mod in protein sequence
            RMAList(n,4) = DF;   %Unmodified degrees of freedom.
            RMAList(n,5) = {ProtSeqChar(1,num)}; %Modified residue.
            RMAList(n,6+GroupNum:9+GroupNum) = [protein{1,1},...
                protein{1,2}, {score}, {title}];
            n = n + 1;
        end       
    end
end

RMAList{2,1} = 'IDENTIFIER';
RMAList{2,2} = 'PEPTIDES USED IN CACLULATION';
RMAList{2,3} = 'MOD POSITION IN PROTEIN';
RMAList{2,4} = 'UNMODIFIED DEGREES OF FREEDOM';
RMAList{2,5} = 'MODIFIED RESIDUE';
RMAList(2,6:5+GroupNum) = GroupList;
RMAList{1,6} = 'RELATIVE PTM LOG2(FOLD-CHANGE)';
RMAList{1,6+GroupNum} = 'UNIPROT ID';
RMAList{1,7+GroupNum} = 'PROTEIN NAME';
RMAList{1,8+GroupNum} = 'MASCOT SCORE';
RMAList{1,9+GroupNum} = 'SITE TITLE FOR GRAPHS';
if strcmp(reg,'bayes')  || strcmp(reg, 'ME')
    RMAList(2,10+GroupNum:9+GroupNum*4) = repmat(GroupList,1,3);
    RMAList(1,10+GroupNum) = {'SE'};
    RMAList(1,10+GroupNum*2) = {'P-VALUE'};
    RMAList(1,10+GroupNum*3) = {'FDR-ADJUSTED P-VALUE'};
end

RMALength = size(RMAList,1);
for index = 3:RMALength
    if isempty(RMAList{index,3})
        RMAList = RMAList(1:(index - 1), :);
        break;
    end
end   
end