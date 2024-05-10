function sigDEN = func_denoise_dw1d_wibx(SIG)
% FUNC_DENOISE_DW1D Saved Denoising Process.
%   SIG: vector of data
%   -------------------
%   sigDEN: vector of denoised data

%  Auto-generated by Wavelet Toolbox on 25-May-2019 15:59:33

% Analysis parameters.
%---------------------
wname = 'db3';
level = 8;

% Denoising parameters.
%----------------------
% meth = 'rigrsure';
% scal_or_alfa = one;
sorh = 's';    % Specified soft or hard thresholding
thrSettings =  [...
    0.000010909454841 ; ...
    0.000011732943091 ; ...
    0.000011525098660 ; ...
    0.000010528126270 ; ...
    0.000011519742595 ; ...
    0.000008892886480 ; ...
    0.000010415756810 ; ...
    0.000007477494662   ...
    ];

% Denoise using CMDDENOISE.
%--------------------------
sigDEN = cmddenoise(SIG,wname,level,sorh,NaN,thrSettings);