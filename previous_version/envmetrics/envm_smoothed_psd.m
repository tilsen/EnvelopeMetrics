function [smpsd,freq] = envm_smoothed_psd(x,par)

% signal: an amplitude envelope (from envm_band_energy.m)
% par:additional parameters (see below)

dbstop if error;


%default parameters if not supplied
if ~isfield(par,'nfft'), par.nfft = 2048; end
if ~isfield(par,'sm_Hz'), par.sm_Hz = 1; end
if ~isfield(par,'L')
    par.L = fix(par.nfft*par.sm_Hz/par.Fs);
    if ~mod(par.L,2), par.L = par.L+1; end
end

if length(x)>par.nfft
    fprintf('WARNING: signal longer than nfft\n');    
end

Fs      = par.Fs;
nfft    = par.nfft;
L       = par.L;
x       = [x(:)' zeros(1,nfft-length(x))]';
N       = length(x);                                    %length of padded signal
X       = fft(x,nfft);                                  %fast Fourier transform
P       = (abs(X).^2)/N;                                %power normalized (/N)
psd     = 2*P(1:N/2);                            %powsd=P(1:N/2); multiply by 2 b/c want normalized power
freq    = Fs*(0:N/2-1)/N;                               %frequencies

psdsym  = [flipud(psd); psd; flipud(psd)];              %treat spectrum as symmetric about 0 and N (cf. Chatfield)

smpsd = filter((1/L)*ones(1,L),1,psdsym);               %moving average (can use other smoothing filters here)
smpsd = smpsd((N/2)+1:(N/2)+(N/2));

end




