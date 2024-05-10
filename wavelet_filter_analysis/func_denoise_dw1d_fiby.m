function sigDEN = func_denoise_dw1d(SIG)
% FUNC_DENOISE_DW1D Saved Denoising Process.
%   SIG: vector of data
%   -------------------
%   sigDEN: vector of denoised data

%  Auto-generated by Wavelet Toolbox on 25-May-2019 16:16:27

% Analysis parameters.
%---------------------
wname = 'db3';
level = 4;

% Denoising parameters.
%----------------------
% meth = 'rigrsure';
% scal_or_alfa = one;
sorh = 's';    % Specified soft or hard thresholding
thrSettings =  [...
    0.095388131516469 ; ...
    0.044538744413803 ; ...
    0.139077307216529 ; ...
    0.186060288708523   ...
    ];

% Denoise using CMDDENOISE.
%--------------------------
sigDEN = cmddenoise(SIG,wname,level,sorh,NaN,thrSettings);
