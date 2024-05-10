function Fib_den=Fib_Filter(Fib_b_nois)


FibX_nois=Fib_b_nois(:,1);
FibY_nois=Fib_b_nois(:,2);
FibZ_nois=Fib_b_nois(:,3);

%% Denoising Step 01  => Wavelet


FibX_den=func_denoise_dw1d_fibx(FibX_nois);
FibY_den=func_denoise_dw1d_fiby(FibY_nois);
FibZ_den=func_denoise_dw1d_fibz(FibZ_nois);


Fib_den=[FibX_den' FibY_den' FibZ_den'];
