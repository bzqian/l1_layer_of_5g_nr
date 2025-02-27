function gen_time_offset_table()

global frame_cfg;
global PuschCom;

frame_cfg.N_fft = 4096;
frame_cfg.N_sc = 273*12;
frame_cfg.N_cp_first = 19; % scs is 30KHz, normal cyclic prefix

for fft_idx=0:frame_cfg.N_fft-1
    tempTO = (2*pi*(fft_idx))/frame_cfg.N_fft;
    pTemp(fft_idx+1) = cos(tempTO) + 1i*sin(tempTO);
end

for sc_idx=0:frame_cfg.N_sc-1
    if sc_idx < frame_cfg.N_sc/2
        RealN(sc_idx+1) = sc_idx + (frame_cfg.N_fft - frame_cfg.N_sc/2);
    elseif sc_idx >= frame_cfg.N_sc/2
        RealN(sc_idx+1) = sc_idx - frame_cfg.N_sc/2;
    end
end

nSeq = 0;
for nTimeOffset = -frame_cfg.N_cp_first:frame_cfg.N_cp_first
    for sc_idx = 0:frame_cfg.N_sc-1
        a = RealN(sc_idx+1);
        b = mod(a*nTimeOffset, frame_cfg.N_fft);
        PuschCom.TimeOffsetTable(nSeq+1) = pTemp(b+1);
        nSeq = nSeq+1;
    end
end

end