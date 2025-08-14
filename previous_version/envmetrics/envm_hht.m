function [imf,w,t] = envm_hht(x,Ts)

%Uses emd functions from 
%plot_hht.m from Alan Tan (http://www.mathworks.com/matlabcentral/fileexchange/19681-hilbert-huang-transform/content/plot_hht.m) 
%but modified in one case (see below)

%% Hilbert-Huang transform:
imf = emd(x); 
for k = 1:length(imf)
   th = angle(hilbert(imf{k}));
   th = unwrap(th); %ST added
   w{k} = diff(th)/Ts/(2*pi);           %instantaneous frequencies
end

%%
N = length(x);
t = linspace(0,(N-1)*Ts,N);

end


function imf = emd(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% imf = emd(x)
% Func : findpeaks

x   = transpose(x(:));
imf = [];
while ~ismonotonic(x)
   x1 = x;
   sd = Inf;
   cnt = 0; %ST added
   while (sd > 0.1) || ~isimf(x1)
      s1 = getspline(x1);
      s2 = -getspline(-x1);
      x2 = x1-(s1+s2)/2;
      
      sd = sum((x1-x2).^2)/sum(x1.^2);
      x1 = x2;
      cnt=cnt+1; if cnt>10, break; end %ST
      
   end
   
   imf{end+1} = x1;
   x          = x-x1;
end
imf{end+1} = x;
end

% FUNCTIONS

function u = ismonotonic(x)

u1 = length(findpeaks(x))*length(findpeaks(-x));
if u1 > 0, u = 0;
else      u = 1; end
end

function u = isimf(x)

N  = length(x);
%u1 = sum(x(1:N-1).*x(2:N) < 0); %original
%---------------ST modified because 0-valued samples make this not work
u1 = sum(abs(diff(sign(x)))/2);
u2 = length(findpeaks(x))+length(findpeaks(-x));
if abs(u1-u2) > 1, u = 0;
else              u = 1; end
end

function s = getspline(x)

N = length(x);
p = findpeaks(x);
s = spline([0 p N+1],[0 x(p) 0],1:N);
end

function n = findpeaks(x)
% Find peaks.
% n = findpeaks(x)

n    = find(diff(diff(x) > 0) < 0);
u    = find(x(n+1) > x(n));
n(u) = n(u)+1;
end