function [startint,endint,state_tally]=ECG_Statistics(ECG_input,Accel_input,Stats_Matrix, ECG_Emissions_File, Accel_Emissions_File,subject_ID,frequency)
    subj_ecg_data=readtable(ECG_input,'delimiter',','); %read in subject data file
    subj_ecg=subj_ecg_data(:,2); %read in ecg data
    
    subj_accel_data=readtable(Accel_input,'delimiter',','); %read in subject ECG data file
    subj_accel=subj_accel_data(:,2:4); %read in accel data
    
    ecg_emissions_matrix=table2array(readtable(ECG_Emissions_File,'delimiter','\t')); %read in accel emission matrix
    accel_emissions_matrix=table2array(readtable(Accel_Emissions_File,'delimiter','\t')); %read in accel emissions matrix
    stats_matrix=readtable(Stats_Matrix,'delimiter','\t'); %read in statistics file for mean and SD 
    transition_matrix=[0.21207 0.03185 0.02766; 0.03772 0.24476 0.06621; 0.02347 0.08382 0.27242]; %make transition matrix
    duration_params=[3.4554 1.0628;3.2273 0.7843;3.7889 1.0853]; %mean and SD for lognormal duration probabilities
    %duration_params=[3.3156 1.0250;3.2273 0.7843;3.7889 1.0853];
    
    timestamps=subj_ecg(:,1); %get time stamp information
    lint=size(timestamps,1)-(250*300); %ignore indices after this one
    %startint=randi([0 lint],1,1);%randomly choose start index for testing period
    startint=664345;
    endint=startint+(250*300)-1; %calculate end index for testing period
    startint_accel=ceil((startint/250)*31.25);
    endint_accel=startint_accel+9359;
    test_ecg_data=table2array(subj_ecg(startint:endint,:)); %subset ecg data
    
    %Filter,calculate magnitude, and standardize accel data
    hd=design_filter;
    s_filt=filter(hd,table2array(subj_accel));
    accel_magnitude=sqrt(sum((s_filt.^2),2)); %get magnitude of accel data
    accel_mean=mean(accel_magnitude);
    accel_sd=std(accel_magnitude);
    accel_subset=accel_magnitude(startint_accel:endint_accel,:); %subset accel data
    standardized_accel=(accel_subset-accel_mean)./accel_sd;
   
    %generate ecg cell array for viterbi
    num_mat_ecg=(endint+1-startint)/(250*5); %figure out how many submatrics there will need to be
    sec_mat_ecg=(250*5)* ones(1,num_mat_ecg);
    splitmats_ecg=mat2cell(test_ecg_data,sec_mat_ecg,1); %split ecg into parts of 1250
    
    %generate accel cell array for viterbi
    num_mat_accel=9360/156;
    sec_mat_accel=156*ones(1,num_mat_accel);
    splitmats_accel=mat2cell(standardized_accel,sec_mat_accel,1); %split accel into parts of 156
    
    %generate vector for storing labels from viterbi
    state_tally=zeros(length(splitmats_ecg),1);
    
    %generate edges vector for histogram
    X=(1:100)';
    X=[0;X]; 
    X=X/100;
    
    %iterate through splitmats and generate labels vector using viterbi
    init_dur=5; %initiate state duration counter
    init_prob=1; %initial previous state probability
    init_prev=randi([1 3],1,1);%randomly select a state to transition from
    
    %perform viterbi and generate most probable state labels:
    %1=rest,2=aroused,3=active
    
    for ii=1:length(splitmats_ecg)
        set(0,'DefaultFigureVisible','off');
        ecg_matrix=splitmats_ecg{ii};
        accel_matrix=splitmats_accel{ii};
        f=@pan_tompkin; %perform pan tompkins algorithm, calculate RR interval, standardize ECG
        [qrs_amp_raw,qrs_i_raw,delay]=f(ecg_matrix,frequency,0);
        R_indices=qrs_i_raw';
        R_diff=diff(R_indices); %get the differences between R interval indices
        RR_interval=R_diff/250; %convert to RR interval per second
        ecg_matrix=(RR_interval-stats_matrix.Mean(subject_ID))./stats_matrix.SD(subject_ID);
        scaled_ecg_matrix=(ecg_matrix-min(ecg_matrix))/(max(ecg_matrix)-min(ecg_matrix));
        scaled_accel_matrix=(accel_matrix-min(accel_matrix))/(max(accel_matrix)-min(accel_matrix));
        [N1,edges1,bin1]=histcounts(scaled_ecg_matrix,X); 
        [N2,edges2,bin2]=histcounts(scaled_accel_matrix,X);
        rd_prob=pdf('lognormal',init_dur,duration_params(1,1),duration_params(1,2)); %calculate state duration probability
        rd_ECG_IDX=find(N1);
        rd_Accel_IDX=find(N2);
        rd_ecg_emissions=sum(ecg_emissions_matrix(rd_ECG_IDX,1));
        rd_accel_emissions=sum(accel_emissions_matrix(rd_Accel_IDX,1));
        total_rd_emissions=rd_ecg_emissions*rd_accel_emissions;
        max_rd=total_rd_emissions*init_prob*rd_prob*transition_matrix(init_prev,1); %calculate total probability of rest state given params
        ard_prob=pdf('lognormal',init_dur,duration_params(2,1),duration_params(2,2));
        ard_ecg_emissions=sum(ecg_emissions_matrix(rd_ECG_IDX,2));
        ard_accel_emissions=sum(accel_emissions_matrix(rd_Accel_IDX,2));
        total_rd_emissions=ard_ecg_emissions*ard_accel_emissions;
        max_ard=total_rd_emissions*init_prob*ard_prob*transition_matrix(init_prev,2); %calculate total probability of aroused state given params
        at_prob=pdf('lognormal',init_dur,duration_params(3,1),duration_params(3,2));
        at_ecg_emissions=sum(ecg_emissions_matrix(rd_ECG_IDX,3));
        at_accel_emissions=sum(accel_emissions_matrix(rd_Accel_IDX,3));
        total_at_emissions=at_ecg_emissions*at_accel_emissions;
        max_at=total_at_emissions*init_prob*at_prob*transition_matrix(init_prev,3); %calculate total probability of active state given params
        state_total=[max_rd;max_ard;max_at];
        [~,idx]=max(state_total);
        init_prev=idx; %update previous state
        init_prob=max(state_total); %update previous state probability
        %if there is a state change reset to five
        state_tally(ii,:)=idx; %append most probable state label to vector
        if idx ~= 1  
            init_dur=5;
        elseif idx~=1 && idx~=state_tally(ii-1,:)
            init_dur=5;
        elseif idx~=1 && idx==state_tally(ii-1,:)
            init_dur=init_dur+5; %update state duration
        end 
    end
    writetable(array2table(state_tally),'Predicted_Labels_Subj6.txt','Delimiter','\t');
end 