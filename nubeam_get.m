function nubeam_struct=nubeam_get(transp,t)
    %
    % gets information about nubeam

    nubeam_struct={};
    %% collect dimensions
    time = transp.coords.TIME.data;
    rho = transp.coords.X.data; %time and index dependant
    if isempty(t) | t<min(time) | t>max(time)
        fprintf('Set time to 0, timerange = %s - %s \n',num2str(min(time)), num2str(max(time)));
        t=0;
    end
    [~,ind]  = min(abs(t-time));
    
    nubeam_struct.time=time;
    nubeam_struct.rho=rho(:,ind);
    nubeam_struct.ind=ind;
    nubeam_struct.id = transp.id;
    %% collect 1D data
    p_inj = transp.allvars.PINJ.data*1e-6; % total injected power
    p_ST = transp.allvars.BPSHI.data*1e-6; % shine-through
    p_OL = transp.allvars.BPLIM.data*1e-6; % orbit losses
    p_e = transp.allvars.BPTE.data*1e-6;
    p_i = transp.allvars.BPTI.data*1e-6;
    p_th = transp.allvars.BPTH.data*1e-6;
    p_CX = p_inj-p_ST-p_OL-p_e-p_i-p_th;

    neut = transp.allvars.NEUTT.data;
    neut_DD = transp.allvars.NEUTX_DD.data;
    neut_thnuc = transp.allvars.NEUTX.data;
    nubeam_struct.d1 = {};
    names_1d={'p_inj', 'p_ST','p_OL', 'p_CX', 'p_e', 'p_i', 'p_th', 'neut'};
    for el=1:length(names_1d)
        fn = sprintf('%s', names_1d{el});
        nubeam_struct.d1.(fn)=eval(fn);
    end
    %% collect 2D data
    area = transp.allvars.DAREA.data*1e-4; %differential area (area of one flux tube)
    vol = transp.allvars.DVOL.data*1e-6;%differential volume (volume of one flux tube)
    j_beam = transp.allvars.CURB.data*10.; % kA/m^2
    pe_beam = transp.allvars.PBE.data; %MW/m3
    pi_beam = transp.allvars.PBI.data;%MW/m3
    n_beam = transp.allvars.BDENS.data*1e6; %1/m3
    pr_beam = transp.allvars.PMHDF_IN.data;
    nubeam_struct.d2 = {};

    for el={'area', 'vol', 'j_beam','pe_beam', 'pi_beam', 'n_beam', 'pr_beam'}
        fn = sprintf('%s', el{1});
        vv=eval(fn); vv=vv(:, ind);
        nubeam_struct.d2.(fn)=vv;
    end   
    %% integrate 2D data
   % I_beam = dot(j_beam/10.,area)*1e-3; %kA

return