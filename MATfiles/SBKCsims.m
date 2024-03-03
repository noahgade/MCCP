clear
load("~/Desktop/Research/Rodu/MCCP/Data/tar0.mat")
load("~/Desktop/Research/Rodu/MCCP/Data/tar1.mat")
load("~/Desktop/Research/Rodu/MCCP/Data/tar2.mat")
load("~/Desktop/Research/Rodu/MCCP/Data/gp0.mat")
load("~/Desktop/Research/Rodu/MCCP/Data/gp1.mat")
load("~/Desktop/Research/Rodu/MCCP/Data/gp2.mat")

for arl = [1, 2, 3, 4, 5, 6]
    tSB0 = cell(100, 1); 
    gSB0 = cell(100, 1); 
    tSB1 = cell(100, 1); 
    gSB1 = cell(100, 1);
    tSB2 = cell(100, 1); 
    gSB2 = cell(100, 1);
    tKC0 = cell(100, 1); 
    gKC0 = cell(100, 1);
    tKC1 = cell(100, 1); 
    gKC1 = cell(100, 1);
    tKC2 = cell(100, 1); 
    gKC2 = cell(100, 1);

    for i = 1:100
        tSBKC0 = SBKC(transpose(tar0(:, :, i)), 10 ^ arl);
        tSB0{i, 1} = tSBKC0{1, 1};
        tKC0{i, 1} = tSBKC0{1, 2};

        tSBKC1 = SBKC(transpose(tar1(:, :, i)), 10 ^ arl);
        tSB1{i, 1} = tSBKC1{1, 1};
        tKC1{i, 1} = tSBKC1{1, 2};

        tSBKC2 = SBKC(transpose(tar2(:, :, i)), 10 ^ arl);
        tSB2{i, 1} = tSBKC2{1, 1};
        tKC2{i, 1} = tSBKC2{1, 2};

        gSBKC0 = SBKC(transpose(gp0(:, :, i)), 10 ^ arl);
        gSB0{i, 1} = gSBKC0{1, 1};
        gKC0{i, 1} = gSBKC0{1, 2};

        gSBKC1 = SBKC(transpose(gp1(:, :, i)), 10 ^ arl);
        gSB1{i, 1} = gSBKC1{1, 1};
        gKC1{i, 1} = gSBKC1{1, 2};

        gSBKC2 = SBKC(transpose(gp2(:, :, i)), 10 ^ arl);
        gSB2{i, 1} = gSBKC2{1, 1};
        gKC2{i, 1} = gSBKC2{1, 2};
        disp(i + "/" + arl);
    end

    writecell(tSB0, "~/Desktop/Research/Rodu/MCCP/SB/tSB0_" + arl + ".csv")
    writecell(tSB1, "~/Desktop/Research/Rodu/MCCP/SB/tSB1_" + arl + ".csv")
    writecell(tSB2, "~/Desktop/Research/Rodu/MCCP/SB/tSB2_" + arl + ".csv")
    writecell(tKC0, "~/Desktop/Research/Rodu/MCCP/KC/tKC0_" + arl + ".csv")
    writecell(tKC1, "~/Desktop/Research/Rodu/MCCP/KC/tKC1_" + arl + ".csv")
    writecell(tKC2, "~/Desktop/Research/Rodu/MCCP/KC/tKC2_" + arl + ".csv")

    writecell(gSB0, "~/Desktop/Research/Rodu/MCCP/SB/gSB0_" + arl + ".csv")
    writecell(gSB1, "~/Desktop/Research/Rodu/MCCP/SB/gSB1_" + arl + ".csv")
    writecell(gSB2, "~/Desktop/Research/Rodu/MCCP/SB/gSB2_" + arl + ".csv")
    writecell(gKC0, "~/Desktop/Research/Rodu/MCCP/KC/gKC0_" + arl + ".csv")
    writecell(gKC1, "~/Desktop/Research/Rodu/MCCP/KC/gKC1_" + arl + ".csv")
    writecell(gKC2, "~/Desktop/Research/Rodu/MCCP/KC/gKC2_" + arl + ".csv")
end


       