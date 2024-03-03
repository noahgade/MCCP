function CPs = SBKC(data, arl)
    [L, ~] = size(data);
    ref_size = 160;
    ref_data = data(1:ref_size, :);
    test_data = data((ref_size + 1):L, :);
    M = 31;
    omega_B = 2:2:M;
    N = 5;
    flagScanB = 0;
    nptsScanB = 0;
    flagKCUSUM = 0;
    nptsKCUSUM = 0;
    cpScanB = 0;
    cpKCUSUM = 0;
    r = 1;
    bandw = r * bandw1(ref_data);
    ARL = arl;
    threshold = find_thre(M, ARL);
    stat_seq0 = online_kernel_cusum(ref_data, test_data, omega_B, N, bandw);
    stat_seqScanB = stat_seq0;
    while(flagScanB < 1)
        if(size(find(stat_seqScanB(:, 2) > threshold, 1, "first"), 1) > 0)
            tempcpScanB = find(stat_seqScanB(:, 2) > threshold, 1, "first") + ref_size + cpScanB(1, 1 + nptsScanB);
            cpScanB = cat(2, cpScanB, tempcpScanB);
            nptsScanB = nptsScanB + 1;
               if(tempcpScanB > L - ref_size - M)
                   flagScanB = 1;
               else
                   ref_data = data((tempcpScanB + 1):(tempcpScanB + ref_size), :);
                   test_data = data((tempcpScanB + ref_size + 1):L, :);
                   bandw = r * bandw1(ref_data);
                   stat_seqScanB = online_kernel_cusum(ref_data, test_data, omega_B, N, bandw);
               end
        else
            cpScanB = cat(2, cpScanB, L * ones(1, 1));
            flagScanB = 1;
        end
    end
    stat_seqKCUSUM = stat_seq0;
    while(flagKCUSUM < 1)
        if(size(find(stat_seqKCUSUM(:, 1) > threshold, 1, "first"), 1) > 0)
            tempcpKCUSUM = find(stat_seqKCUSUM(:, 1) > threshold, 1, "first") + ref_size + cpKCUSUM(1, 1 + nptsKCUSUM);
            cpKCUSUM = cat(2, cpKCUSUM, tempcpKCUSUM);
            nptsKCUSUM = nptsKCUSUM + 1;
               if(tempcpKCUSUM > L - ref_size - M)
                   flagKCUSUM = 1;
               else
                   ref_data = data((tempcpKCUSUM + 1):(tempcpKCUSUM + ref_size), :);
                   test_data = data((tempcpKCUSUM + ref_size + 1):L, :);
                   bandw = r * bandw1(ref_data);
                   stat_seqKCUSUM = online_kernel_cusum(ref_data, test_data, omega_B, N, bandw);
               end
        else
            cpKCUSUM = cat(2, cpKCUSUM, L * ones(1, 1));
            flagKCUSUM = 1;
        end
    end
    CPs = {cpScanB, cpKCUSUM};
end
