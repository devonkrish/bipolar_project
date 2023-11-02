%For high density grids only
freq_table = '/Users/devonkrish/Desktop/IED/_BipolarReref/baseline-high-density-data/tableoutputfile.xlsx';
freqs = [];
bpds = [];
powers = [];

T = readtable(freq_table);
pts_ids = T{1:9,1};
for i = 1:height(T)
    freqs = [freqs, T{i,2}];
    bpds = [bpds, T{i,3}];
    powers = [powers, T{i,4}];
end

delta_index = [];
theta_index = [];
alpha_index = [];
beta_index = [];
gamma_index = [];
highgamma_index = [];

for j = 1:length(freqs)

    if strcmp(freqs{j},'delta')
        delta_index = [delta_index, j];
    elseif strcmp(freqs{j},'theta')
        theta_index = [theta_index, j];
    elseif strcmp(freqs{j},'alpha')
        alpha_index = [alpha_index, j];
    elseif strcmp(freqs{j},'beta')
        beta_index = [beta_index, j];
    elseif strcmp(freqs{j},'gamma')
        gamma_index = [gamma_index, j];
    else
        highgamma_index = [highgamma_index, j];
    end
end

indexarr = delta_index;
freqbandarr = [];
hold = [];
j = indexarr(1);
initj = 1;

for i = 1:9
    while j<=indexarr(54)
        hold = [hold; powers(j)];
        j = j+9;
    freqbandarr = [freqbandarr; hold];
    hold = [];

    initj = initj + 1;
    j = initj + indexarr(1);
    end
end











