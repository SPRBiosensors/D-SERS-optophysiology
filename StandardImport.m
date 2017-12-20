function [StandardSpectra,StandardShift,nbref,RefName] = StandardImport(SERSDatab,varargin)

if strcmp(varargin,'all') == 1; nbref = size(SERSDatab,1);   
    else nbref = nargin - 1;
end

StandardSpectra = zeros(790,nbref); StandardShift = zeros(1024,nbref);
Row = zeros(size(SERSDatab,1),1); RefName = cell(nbref,1);

    for k = 1:nbref;
        if strcmp(varargin,'all') == 1; StandardSpectra(:,k) = SERSDatab{k,4};
           StandardShift(:,k) = SERSDatab{k,5}; RefName{k,:} = SERSDatab{k,1};
        else R = strcmp(SERSDatab,varargin{k});
        for i = 1:size(R,1);
            if R(i,1) == 1; Row(i,1) = i;
            else Row(i,1) = 0;
            end
        end
        Row(Row == 0) = []; StandardSpectra(:,k) = SERSDatab{Row,4};
        StandardShift(:,k) = SERSDatab{Row,5}; RefName{k,:} = SERSDatab{Row,1};
        end
    end
end

