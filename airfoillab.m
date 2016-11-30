clear all
close all
clc

addpath('WindTunnelData');

for k = 1:9
    filename = ['AirFoilPressure_S011_G0',sprintf('%d', 3), '.csv'];
    data = csvread(filename,1);
    Patm1 = data(:,1);
    Tatm1 = data(:,2);
    Rhoatm1 = data(:,3);
    Airspeed = data(:,4);
    PitotDynP = data(:,5);
    AuxDynP = data(:,6);
    PortPres = data(:,7:22);
    P1=data(:,7);
    P2=data(:,8);
    P3=data(:,9);
    P4=data(:,10);
    P5=data(:,11);
    P6=data(:,12);
    P7=data(:,13);
    P8=data(:,14);
    P9=data(:,15);
    P10=data(:,16);
    P11=data(:,17);
    P12=data(:,18);
    P13=data(:,19);
    P14=data(:,20);
    P15=data(:,21);
    P16=data(:,22);
    AngleofA=data(:,23);
    NormalF = data(:,24);
    AxialF = data(:,25);
    PitchMom = data(:,26);
    ELDX = data(:,27);
    ELDY = data(:,28);
    
    Patm=mean(Patm1); %taking the average of the atm pressure
    Rhoatm=mean(Rhoatm1);%take the average of the atm density
    
    airfoilX =0.025*[ 0 .175 .35 .7 1.05 1.4 1.75 2.1 2.8 2.8 2.1 1.4 1.05 .7 .35 .175];
    airfoilY =0.025*[0.14665 .33075 .4018 .476 .49 .4774 .4403 .38325  .21875 0 0 0 0 .0014 .0175 .03885];
    cord = sqrt(airfoilY(1)^2+3.5^2);
    cord = cord *0.025;
    Aoavg = [];
     for i=1:9
        % There are 20 repeat values
        count = 20;
        j = i;
        if count > 1
            % Determine average of the number of repeat values
            av_Vel = 0;
            for m = i:i+count-1
                av_Vel = av_Vel + Airspeed(m);
            end
            av_Vel = av_Vel/count;
            % Set the first of the repeat values as the average
            Airspeed(i) = av_Vel;
            while count > 1
                % Get rid of repeat values
                n = i+1;
                Airspeed(n) = [];
                count = count - 1;
                AngleofA(n) = [];
            end
        end
    end
    PressAvg = [];
    count = 1;
    for i = 1:20:160
        for j = 1:16
            avg = mean(PortPres(i:(i+20),j));
            PressAvg(count,j) = avg;
        end
        count=count+1;
    end
    for i = 1:16
        avg = mean(PortPres(161:180,i));
        PressAvg(9,i) = avg;
    end
    %PressAvg
    
    Cn=0;
    Ca=0;
    sumDx=0;
    sumDy=0;
    Cn=zeros(9,16);
    Ca =zeros(9,16);

    %AOA -6
    for i=1:9
        for j = 1:16
            if j == 16
                Cp1=((PressAvg(i,j))/(.5*Rhoatm*Airspeed(i).^2));
                Cp2=((PressAvg(i,1))/(.5*Rhoatm*Airspeed(i).^2));
                if k == 2
                    CpVec(i,j) = Cp1;
                end
                Dx=(airfoilX(1)-airfoilX(j))/cord;
                Dy=(airfoilY(1)-airfoilY(j))/cord;
                SumDx=sumDx+Dx;
                sumDy=sumDy+Dy;
                Cn(i,j)=(Cn(i,j)+.5*(Cp1+Cp2)*Dx);
                Ca(i,j)= Ca(i,j)+.5*(Cp1+Cp2)*Dy;
            else
                Cp1=((PressAvg(i,j))/(.5*Rhoatm*Airspeed(i).^2));
                Cp2=((PressAvg(i,j+1))/(.5*Rhoatm*Airspeed(i).^2));
                if k == 2
                    
                    CpVec(i,j) = Cp1;
                end
                if i == 9
                    plot(airfoilX(j)/cord,Cp1,'*');hold on
                end
                Dx=(airfoilX(j+1)-airfoilX(j))/cord;
                Dy=(airfoilY(j+1)-airfoilY(j))/cord;
                SumDx=sumDx+Dx;
                sumDy=sumDy+Dy;
                Cn(i,j)=Cn(i,j)+.5*(Cp1+Cp2)*Dx;
                Ca(i,j)= Ca(i,j)+.5*(Cp1+Cp2)*Dy;
            end
            
        end
    end
    
    Cn=-1*Cn;
    for i = 1:9
        AngleofAVec(k,i) = AngleofA(i);
        if AngleofA(i) == 0
            Cl = Cn(i,:);
            Cd = Ca(i,:);
        else
            
            Cl=Cn(i,:)*cos(deg2rad(AngleofA(i)))-Ca(i,:)*sin(deg2rad(AngleofA(i)));
            Cd=Cn(i,:)*sin(deg2rad(AngleofA(i)))+Ca(i,:)*cos(deg2rad(AngleofA(i)));
        end
        a = sum(Cl(1,:));
       % b = sum(Cl(2,:));
       % c = sum(Cl(3,:));
        ClVec(k,i) = a;
        a = sum(Cd(1,:));
        %b = sum(Cd(2,:));
        %c = sum(Cd(3,:))';
        CdVec(k,i) = a;
    end
end
for i = 1:3
end

hold on
for i = 1:3
    if i == 1 
        tmpAoa = [AngleofAVec(9:-1:1,i); AngleofAVec(9:-1:1,i+3); AngleofAVec(9:-1:1, i+6);];
        tmpCl = [ClVec(9:-1:1,i); ClVec(9:-1:1,i+3); ClVec(9:-1:1,i+6);];
        plot(tmpAoa,tmpCl,'r*-')
    elseif i == 2 
        tmpAoa = [AngleofAVec(9:-1:1,i); AngleofAVec(9:-1:1,i+3); AngleofAVec(9:-1:1, i+6);];
        tmpCl = [ClVec(9:-1:1,i); ClVec(9:-1:1,i+3); ClVec(9:-1:1,i+6);];
        plot(tmpAoa,tmpCl,'b*-')
    else
       tmpAoa = [AngleofAVec(9:-1:1,i); AngleofAVec(9:-1:1,i+3); AngleofAVec(9:-1:1, i+6);];
        tmpCl = [ClVec(9:-1:1,i); ClVec(9:-1:1,i+3); ClVec(9:-1:1,i+6);];
        plot(tmpAoa,tmpCl,'g*-')
    end
end
title('Lift Coeffecient vs Angle of Attack')
xlabel('Angle of Attack (^o)')
ylabel('Cl')
legend('10 m/s', '20 m/s', '30 m/s', 'Location', 'northwest')
saveas(gcf, 'Cl_a.png');
figure 
hold on
for i = 1:3
    if i == 1 
        tmpAoa = [AngleofAVec(9:-1:1,i); AngleofAVec(9:-1:1,i+3); AngleofAVec(9:-1:1, i+6);];
        tmpCl = [CdVec(9:-1:1,i); CdVec(9:-1:1,i+3); CdVec(9:-1:1,i+6);];
        plot(tmpAoa,tmpCl,'r*-')
    elseif i == 2 
        tmpAoa = [AngleofAVec(9:-1:1,i); AngleofAVec(9:-1:1,i+3); AngleofAVec(9:-1:1, i+6);];
        tmpCl = [CdVec(9:-1:1,i); CdVec(9:-1:1,i+3); CdVec(9:-1:1,i+6);];
        plot(tmpAoa,tmpCl,'b*-')
    else
       tmpAoa = [AngleofAVec(9:-1:1,i); AngleofAVec(9:-1:1,i+3); AngleofAVec(9:-1:1, i+6);];
        tmpCl = [CdVec(9:-1:1,i); CdVec(9:-1:1,i+3); CdVec(9:-1:1,i+6);];
        plot(tmpAoa,tmpCl,'g*-')
    end
end
title('Drag Coefficient vs Angle of Attack')
xlabel('Angle of Attack (^o)')
ylabel('Cd')
legend('10 m/s', '20 m/s', '30 m/s', 'Location', 'northwest')
saveas(gcf, 'Cd_a.png');
topLine = -(CpVec(:,9)-CpVec(:,8))/((airfoilX(9)-airfoilX(8))/cord);
botLine = -(CpVec(:,11) - CpVec(:,10))/((airfoilX(11)-airfoilX(10))/cord);
theo = [];
theo = topLine*((airfoilX(9)/cord)-1)+CpVec(:,9);
theo1 = [];
theo2 = botLine*((airfoilX(11)/cord)-1)+CpVec(:,11);
avg = (theo+theo2)/2;
CpVec = [CpVec(:,1:9) avg CpVec(:,10:end)];
airfoilX = [airfoilX(1:9) 3.5*.025 airfoilX(10:end)];


for i = 1:9
    figure;
    hold on
    plot(airfoilX/cord, CpVec(i,:), ':r*')
    aoa = [-6 -6 -6 4 4 4 14 14 14];
    airVel = [10 20 30 10 20 30 10 20 30];
    tittle = ['Velocity ',sprintf('%d',airVel(i)),' Angle of Attack ',sprintf('%d',aoa(i))];
    title(tittle)
    ylabel('Cp')
    xlabel('Normalized Chord Length')
    set(gca,'Ydir','reverse')
    hold off
end
