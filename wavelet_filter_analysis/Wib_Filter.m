function Wib_den=Wib_Filter(Wib_b_nois)

WibX_nois=Wib_b_nois(:,1);
WibY_nois=Wib_b_nois(:,2);
WibZ_nois=Wib_b_nois(:,3);


%% Denoising Step 01  => Wavelet

WibX_den=func_denoise_dw1d_wibx(WibX_nois);
WibY_den=func_denoise_dw1d_wibx(WibY_nois);
WibZ_den=func_denoise_dw1d_wibz(WibZ_nois);


Wib_den=[WibX_den' WibY_den' WibZ_den'];
