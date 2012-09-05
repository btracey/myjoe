function ExtractProfile()
    clear all; clc; close all;

    turbmod = menu('turbulence model being used?','k-e','k-w','SST','SA');
    % import data file according to turbulence model
    switch turbmod
        case 1
            [XY,scalars,kine,eps,muLam] = importKEprofile();
        case 2
            [XY,scalars,kine,omega,muLam] = importKWprofile();
        case 3
            [XY,scalars,kine,omega,muLam,dudy,muT] = importSSTprofile();
        case 4
            [XY,scalars,nusa,muLam] = importSAprofile();
    end
    
    % from saved scalars compute rho, vel, press
    [rho,vel,press] = scalar_decomp(scalars);

    % Summarize Data about current flow to make sure it is what the user
    % desired: Mach, Density, Pressure, Temp of free stream
    % average all values from cell row next to top boundary as a measure of
    % freestream conditions - may be some effects from neumann condition of
    % boundary layer is too close to top
    srho = mean(rho(:,end-1));
    spress = mean(press(:,end-1));
    stemp = spress / srho / 287.06;
    smach = mean(vel(:,end-1,1)) / sqrt(1.4*spress/srho);
    disp('Freestream conditions are:')
    disp(['rho=',num2str(srho)])
    disp(['press=',num2str(spress)])
    disp(['temp=',num2str(stemp)])
    disp(['Mach=',num2str(smach)])
    
    
    %Calculate delta99, ReTheta, displacement thickness, momentum thickness
    [delta,dispth,theta,Retheta] = calc_height(XY,vel,rho,muLam);
    
    varcut = menu('What variable would you like to select with?','delta 99%','disp. thick.'...
        ,'mom. thick.','Re_theta');
    
    switch varcut
        case 1
            des_var = delta;
        case 2
            des_var = dispth;
        case 3
            des_var = theta;
        case 4
            des_var = Retheta;
    end
    
    %ask for Re_theta from user
    desired = input('What value would you like to extract?\n','s');
    while size(str2num(desired)) == 0
        desired = input('Input not a real number, try again.\nWhat value would you like to extract?\n','s');
    end
    while (str2num(desired) < min(des_var)) || (str2num(desired) > max(des_var))
        s = sprintf('Input not in range of data (%f - %f).',min(des_var),max(des_var));
        disp(s);
        desired = input('What value would you like to extract?\n','s');
        while size(str2num(desired)) == 0
            desired = input('Input not a real number, try again.\nWhat value would you like to extract?\n','s');
        end
    end
    desired = str2num(desired);
    
    % sweep through data and cut profile
    for i=(length(des_var)):-1:2
       if desired <= des_var(i) && desired > des_var(i-1)
          imatch = i;
          fak = (desired - des_var(i-1))/(des_var(i) - des_var(i-1));
          break;
       end
    end

    % output profile
    switch turbmod
        case 1
            
        case 2
            
        case 3
            gen_profile_SST(XY,vel,rho,press,kine,omega,imatch,fak);
        case 4
            
    end

    % Plot Figuers of Interest
    close all;
    
    figure(1)
    subplot(2,2,1);
    plot(XY(:,1,1),delta(:,1,1),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[delta(imatch,1,1)*.8 delta(imatch,1,1)*1.2],'--r');
    title('BL height (99%)');xlabel('X (m)');ylabel('\delta_{99%} (m)');

    subplot(2,2,2);
    plot(XY(:,1,1),dispth(:,1,1),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[dispth(imatch,1,1)*.95 dispth(imatch,1,1)*1.05],'--r');
    title('displacement thickness');xlabel('X (m)');ylabel('\delta^* (m)');

    subplot(2,2,3);
    plot(XY(:,1,1),theta(:,1,1),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[theta(imatch,1,1)*.8 theta(imatch,1,1)*1.2],'--r');
    title('displacement momentum');xlabel('X (m)');ylabel('\theta (m)');

    subplot(2,2,4);
    plot(XY(:,1,1),Retheta(:,1,1),'-o',[XY(imatch,1,1) XY(imatch,1,1)],[Retheta(imatch,1,1)*.8 Retheta(imatch,1,1)*1.2],'--r');
    title('Re_theta');xlabel('X (m)');ylabel('Re_{\theta}');
    
    
end

function [XY,scalars,kine,omega,muLam,dudy,muT] = importSSTprofile()
    % read in unstructured output profile, sort by X and Y location to get
    % into structured (i,j) format

    file = 'KOMSST_profile.txt'; %created when last restart is run on 1 processor
    Data = importdata(file);
    
    %spatial locations
    XY = Data(:,1:2);
    
    %count num columns in i
    sortedXY = sortrows(XY,[1 2]);
    imax = 1;
    val = sortedXY(1,2);
    while sortedXY(imax+1,2) > val
        imax = imax+1;
    end
  
    jmax = size(Data,1) / imax;
    
    %sort in X
    Data = sortrows(Data,[1]);
    
    %spatial locations
    XY = zeros(jmax,imax,2);
    
    %Regular Scalars - rho,Xvel,Yvel,press
    scalars = zeros(jmax,imax,4);
    
    %KOMSST specific scalars
    kine = zeros(jmax,imax);
    omega = zeros(jmax,imax);
    dudy = zeros(jmax,imax);
    muT = zeros(jmax,imax);
    
    muLam = zeros(jmax,imax);
    
    %temporary holder
    temp = zeros(imax,size(Data,2));
    
    % sort into structured (i,j) format
    for j = 1:jmax
        temp(1:imax,:) = Data(((j-1)*imax+1):(j*imax),:);
        temp = sortrows(temp,[2]);
        
        for i =1:imax
            XY(j,i,1) = temp(i,1);
            XY(j,i,2) = temp(i,2);
           
            scalars(j,i,1) = temp(i,3);
            scalars(j,i,2) = temp(i,4)./scalars(j,i,1);
            scalars(j,i,3) = temp(i,5)./scalars(j,i,1);            
            scalars(j,i,4) = temp(i,6);  

            kine(j,i,1) = temp(i,7);
            omega(j,i,1) = temp(i,8);
            muLam(j,i,1) = temp(i,9);
            dudy(j,i,1) = temp(i,10);
            muT(j,i,1) = temp(i,11);
        end
    end
    
    clear Data file
end

function [XY,scalars,kine,omega,muLam] = importKWprofile()
    % read in unstructured output profile, sort by X and Y location to get
    % into structured (i,j) format

    file = 'KOM_profile.txt'; %created when last restart is run on 1 processor
    Data = importdata(file);
    
    %spatial locations
    XY = Data(:,1:2);
    
    %count num columns in i
    sortedXY = sortrows(XY,[1 2]);
    imax = 1;
    val = sortedXY(1,2);
    while sortedXY(imax+1,2) > val
        imax = imax+1;
    end
  
    jmax = size(Data,1) / imax;
    
    %sort in X
    Data = sortrows(Data,[1]);
    
    %spatial locations
    XY = zeros(jmax,imax,2);
    
    %Regular Scalars - rho,Xvel,Yvel,press
    scalars = zeros(jmax,imax,4);
    
    %KOM specific scalars
    kine = zeros(jmax,imax);
    omega = zeros(jmax,imax);
    muLam = zeros(jmax,imax);
    
    %temporary holder
    temp = zeros(imax,size(Data,2));
    
    % sort into structured (i,j) format
    for j = 1:jmax
        temp(1:imax,:) = Data(((j-1)*imax+1):(j*imax),:);
        temp = sortrows(temp,[2]);
        
        for i =1:imax
            XY(j,i,1) = temp(i,1);
            XY(j,i,2) = temp(i,2);
           
            scalars(j,i,1) = temp(i,3);
            scalars(j,i,2) = temp(i,4)./scalars(j,i,1);
            scalars(j,i,3) = temp(i,5)./scalars(j,i,1);            
            scalars(j,i,4) = temp(i,6);  

            kine(j,i,1) = temp(i,7);
            omega(j,i,1) = temp(i,8);
            muLam(j,i,1) = temp(i,9);
        end
    end
    
    clear Data file
end

function [XY,scalars,kine,eps,muLam] = importKEprofile()
    % read in unstructured output profile, sort by X and Y location to get
    % into structured (i,j) format

    file = 'KEps_profile.txt'; %created when last restart is run on 1 processor
    Data = importdata(file);
    
    %spatial locations
    XY = Data(:,1:2);
    
    %count num columns in i
    sortedXY = sortrows(XY,[1 2]);
    imax = 1;
    val = sortedXY(1,2);
    while sortedXY(imax+1,2) > val
        imax = imax+1;
    end
  
    jmax = size(Data,1) / imax;
    
    %sort in X
    Data = sortrows(Data,[1]);
    
    %spatial locations
    XY = zeros(jmax,imax,2);
    
    %Regular Scalars - rho,Xvel,Yvel,press
    scalars = zeros(jmax,imax,4);
    
    %KOM specific scalars
    kine = zeros(jmax,imax);
    eps = zeros(jmax,imax);
    muLam = zeros(jmax,imax);
    
    %temporary holder
    temp = zeros(imax,size(Data,2));
    
    % sort into structured (i,j) format
    for j = 1:jmax
        temp(1:imax,:) = Data(((j-1)*imax+1):(j*imax),:);
        temp = sortrows(temp,[2]);
        
        for i =1:imax
            XY(j,i,1) = temp(i,1);
            XY(j,i,2) = temp(i,2);
           
            scalars(j,i,1) = temp(i,3);
            scalars(j,i,2) = temp(i,4)./scalars(j,i,1);
            scalars(j,i,3) = temp(i,5)./scalars(j,i,1);            
            scalars(j,i,4) = temp(i,6);  

            kine(j,i,1) = temp(i,7);
            eps(j,i,1) = temp(i,8);
            muLam(j,i,1) = temp(i,9);
        end
    end
    
    clear Data file
end

function [XY,scalars,nusa,muLam] = importSAprofile()
    % read in unstructured output profile, sort by X and Y location to get
    % into structured (i,j) format

    file = 'SA_profile.txt'; %created when last restart is run on 1 processor
    Data = importdata(file);
    
    %spatial locations
    XY = Data(:,1:2);
    
    %count num columns in i
    sortedXY = sortrows(XY,[1 2]);
    imax = 1;
    val = sortedXY(1,2);
    while sortedXY(imax+1,2) > val
        imax = imax+1;
    end
  
    jmax = size(Data,1) / imax;
    
    %sort in X
    Data = sortrows(Data,[1]);
    
    %spatial locations
    XY = zeros(jmax,imax,2);
    
    %Regular Scalars - rho,Xvel,Yvel,press
    scalars = zeros(jmax,imax,4);
    
    %KOM specific scalars
    nusa = zeros(jmax,imax);
    muLam = zeros(jmax,imax);
    
    %temporary holder
    temp = zeros(imax,size(Data,2));
    
    % sort into structured (i,j) format
    for j = 1:jmax
        temp(1:imax,:) = Data(((j-1)*imax+1):(j*imax),:);
        temp = sortrows(temp,[2]);
        
        for i =1:imax
            XY(j,i,1) = temp(i,1);
            XY(j,i,2) = temp(i,2);
           
            scalars(j,i,1) = temp(i,3);
            scalars(j,i,2) = temp(i,4)./scalars(j,i,1);
            scalars(j,i,3) = temp(i,5)./scalars(j,i,1);            
            scalars(j,i,4) = temp(i,6);  

            nusa(j,i,1) = temp(i,7);
            muLam(j,i,1) = temp(i,8);
        end
    end
    
    clear Data file
end

function [rho,vel,press] = scalar_decomp(scalars)
    % decompose matrix scalars into individual values: rho,u,v,P
    rho = scalars(:,:,1);
    vel(:,:,1) = scalars(:,:,2);
    vel(:,:,2) = scalars(:,:,3);
    press = scalars(:,:,4);
end

function [delta dispth,theta,Retheta] = calc_height(XY,vel,rho,muLam)
    % determine boundary layer thickness (delta 99%)
    delta = zeros(size(XY,1),1);
    for i=1:size(XY,1)
        u0 = vel(i,end,1);
        for j=1:size(XY,2)
           if vel(i,j,1) >= .99*u0
               delta(i,1) = XY(i,j,2);
               break
           end
        end
        clear u0
    end

    % for each line of i generate ReTheta based on compressible
    % formulation:
    % ReTheta = integral( rho(y)/rho0 * u(y)/u0 *(1 - u(y)/u0) dy)
    
    Retheta = zeros(size(XY,1),1);
    
    % calculate for each line i
    for i=1:size(XY,1)
        
        y = XY(i,:,2); %y values
        u0 = vel(i,end,1); %freestream velocity
        inv_u0 = 1 / u0;
        rho0 = rho(i,end); %free stream density
        inv_rho0 = 1 / rho0;
        inv_mu0 = 1 / 8e-6;
        %inv_mu0 = 1 / muLam(i,end);
        
        u_u0 = vel(i,:,1) * inv_u0;
        rho_rho0 = rho(i,:) * inv_rho0;
        
        mom = rho_rho0 .* u_u0 .* (1 - u_u0);
        %displacement thickness
        dispth(i,1) = trapz(y,u_u0);
        %momentum thickness
        theta(i,1) = trapz(y,mom);
        %Retheta
        Retheta(i,1) = theta(i,1) * u0 * rho0 * inv_mu0;
    end

end

function gen_profile_KE(XY,vel,rho,press,kine,eps,imatch,fak)

    y(:,1) = XY(imatch,:,2) + fak*(XY(imatch+1,:,2) - XY(imatch,:,2));
    
    u(:,1) = vel(imatch,:,1) + fak*(vel(imatch+1,:,1) - vel(imatch,:,1));
    v(:,1) = vel(imatch,:,2) + fak*(vel(imatch+1,:,2) - vel(imatch,:,2));
    
    temp = press ./ rho ./ 287.06;
    T(:,1) = temp(imatch,:) + fak*(temp(imatch+1,:) - temp(imatch,:));    
    P(:,1) = press(imatch,:) + fak*(press(imatch+1,:) - press(imatch,:));
    k(:,1) = kine(imatch,:) + fak*(kine(imatch+1,:) - kine(imatch,:));
    ep(:,1) = eps(imatch,:) + fak*(eps(imatch+1,:) - eps(imatch,:));
    
    %output to file here
    prof = fopen('KE_turb_inlet.txt','w');
    fprintf(prof,'n=%d\td=%d\n',length(y),7);
    for j = 1:length(y)
       fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',y(j),u(j),v(j),P(j),T(j),k(j),ep(j)); 
    end

end

function gen_profile_KW(XY,vel,rho,press,kine,omega,imatch,fak)

    y(:,1) = XY(imatch,:,2) + fak*(XY(imatch+1,:,2) - XY(imatch,:,2));
    
    u(:,1) = vel(imatch,:,1) + fak*(vel(imatch+1,:,1) - vel(imatch,:,1));
    v(:,1) = vel(imatch,:,2) + fak*(vel(imatch+1,:,2) - vel(imatch,:,2));
    
    temp = press ./ rho ./ 287.06;
    T(:,1) = temp(imatch,:) + fak*(temp(imatch+1,:) - temp(imatch,:));    
    P(:,1) = press(imatch,:) + fak*(press(imatch+1,:) - press(imatch,:));
    k(:,1) = kine(imatch,:) + fak*(kine(imatch+1,:) - kine(imatch,:));
    om(:,1) = omega(imatch,:) + fak*(omega(imatch+1,:) - omega(imatch,:));
    
    %output to file here
    prof = fopen('KW_turb_inlet.txt','w');
    fprintf(prof,'n=%d\td=%d\n',length(y),7);
    for j = 1:length(y)
       fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',y(j),u(j),v(j),P(j),T(j),k(j),om(j)); 
    end

end

function gen_profile_SST(XY,vel,rho,press,kine,omega,imatch,fak)

    y(:,1) = XY(imatch,:,2) + fak*(XY(imatch+1,:,2) - XY(imatch,:,2));
    
    u(:,1) = vel(imatch,:,1) + fak*(vel(imatch+1,:,1) - vel(imatch,:,1));
    v(:,1) = vel(imatch,:,2) + fak*(vel(imatch+1,:,2) - vel(imatch,:,2));
    
    temp = press ./ rho ./ 287.06;
    T(:,1) = temp(imatch,:) + fak*(temp(imatch+1,:) - temp(imatch,:));    
    P(:,1) = press(imatch,:) + fak*(press(imatch+1,:) - press(imatch,:));
    k(:,1) = kine(imatch,:) + fak*(kine(imatch+1,:) - kine(imatch,:));
    om(:,1) = omega(imatch,:) + fak*(omega(imatch+1,:) - omega(imatch,:));
    
    %output to file here
    prof = fopen('SST_turb_inlet.txt','w');
    fprintf(prof,'n=%d\td=%d\n',length(y),7);
    for j = 1:length(y)
       fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',y(j),u(j),v(j),P(j),T(j),k(j),om(j)); 
    end

end

function gen_profile_SA(XY,vel,rho,press,nusa,imatch,fak)

    y(:,1) = XY(imatch,:,2) + fak*(XY(imatch+1,:,2) - XY(imatch,:,2));
    
    u(:,1) = vel(imatch,:,1) + fak*(vel(imatch+1,:,1) - vel(imatch,:,1));
    v(:,1) = vel(imatch,:,2) + fak*(vel(imatch+1,:,2) - vel(imatch,:,2));
    
    temp = press ./ rho ./ 287.06;
    T(:,1) = temp(imatch,:) + fak*(temp(imatch+1,:) - temp(imatch,:));    
    P(:,1) = press(imatch,:) + fak*(press(imatch+1,:) - press(imatch,:));
    nsa(:,1) = nusa(imatch,:) + fak*(nusa(imatch+1,:) - nusa(imatch,:));
    
    %output to file here
    prof = fopen('SA_turb_inlet.txt','w');
    fprintf(prof,'n=%d\td=%d\n',length(y),6);
    for j = 1:length(y)
       fprintf(prof,'%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n',y(j),u(j),v(j),P(j),T(j),nsa(j)); 
    end

end
