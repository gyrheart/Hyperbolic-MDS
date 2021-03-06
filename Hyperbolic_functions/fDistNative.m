function M = fDistNative(X,varargin)

paramNames = {'Ref'};
paramDflts = {[]};
[RefM] = internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
if length(RefM) ==  length(X)                     
    RefN = 0;
else
    RefN = size(RefM,1);
end
range = RefN+1:length(X);

[n,p] = size(X);
cos_angle = cos(X(range,p)-X(:,p)');
for count_angle = p-1:-1:2
   cos_angle = sin(X(range,count_angle))*sin(X(:,count_angle)').*cos_angle + cos(X(range,count_angle))*cos(X(:,count_angle)');
end
M(range,:) = real(acosh(cosh(X(range,1))*cosh(X(:,1)')-sinh(X(range,1))*sinh(X(:,1)').*cos_angle));
M(:,range) = M(range,:)';
M(1:RefN,1:RefN) = RefM;
for k = 1:n
    M(k,k) = 0;
end

