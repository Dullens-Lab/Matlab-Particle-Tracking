function [MSD] = Msd(xytID)

%input=xyztID
%output=[tau MSD]

lmax=max(xytID(:,3)); 
IDs=unique(xytID(:,4)); 

MSD=zeros(lmax+1,1);
MSD_count=zeros(lmax+1,1);

for i=1:length(IDs)   
    thisID=IDs(i);    
    thispart=xytID(xytID(:,4)==thisID,:);   %Look at each particle in turn

    for dt=0:size(thispart,1)
        for start=1:size(thispart,1)-dt            
        t1=start;
        t2=start+dt;

        realdt=thispart(t2,3)-thispart(t1,3);   %calculate real dt (incase miss frame)  

        v2=[thispart(t2,1:2)]; %position in frame 2
        v1=[thispart(t1,1:2)]; %position in frame 1
            
        MSD_count(realdt+1)=MSD_count(realdt+1)+1;
        MSD(realdt+1)=MSD(realdt+1)+((v2(1)-v1(1)).^2+(v2(2)-v1(2)).^2); 
        end
    end
end

MSD=[0:lmax MSD./MSD_count]; 
loglog([0:lmax],MSD,'b')
hold on


 end


