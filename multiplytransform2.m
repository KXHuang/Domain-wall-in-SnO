clear,close,clc 
fidin=fopen('initial.txt');                                           
fidout=fopen('mkmatlab.txt','w');                      
while ~feof(fidin)                                                   
    tline=fgetl(fidin);                                
    if double(tline(1))==32     
       fprintf(fidout,'%s\n\n',tline);                  
       continue                                        
    end
end
fclose(fidout);
format long
MK=importdata('MKMATLAB.txt'); 

atom1=MK(4,1);
atom2=MK(4,2);
pos=MK(5:end,:);   %设置初始参数，可能加atom，其余不改
sizepos=size(pos);
a=repmat(MK(1,:),sizepos(1),1);
b=repmat(MK(2,:),sizepos(1),1);
c=repmat(MK(3,:),sizepos(1),1);

ma=100;
mb=100;
mc=1;   %设置扩展倍数

ad=[];
for i=0:ma-1
   for j=0:mb-1
       for k=0:mc-1
           ad=[ad;i*a+j*b+k*c];
       end
   end
end

position=repmat(pos,ma*mb*mc,1)+ad;

posa=[];
posb=[];
for m=0:ma*mb*mc-1
    posa=[posa;position((1+m*sizepos(1)):(atom1+m*sizepos(1)),:)];
    posb=[posb;position((atom1+m*sizepos(1)+1):(atom1+atom2+m*sizepos(1)),:)];
end

%transformation
A=MK(1,1)*ma;
B=MK(2,2)*mb;
T=[A/sqrt(A^2+B^2),-B/sqrt(A^2+B^2),0;B/sqrt(A^2+B^2),A/sqrt(A^2+B^2),0;0,0,1];
X=MK(1,1:3)*ma*T;
Y=MK(2,1:3)*mb*T;
posa=posa*T;
posb=posb*T;

for i = 1:length(posa)
    if posa(i,2)<0
       posa(i,1)=-posa(i,1)+sqrt(A^2+B^2)+2.6865529402297097;
    end
end

for i = 1:length(posb)
    if posb(i,2)<0
       posb(i,1)=-posb(i,1)+sqrt(A^2+B^2);
    end
end

%take out the sandwich
wid=1;
len=ma/2;
x2=2.6865529402297097;
y2=1.342202276323422594;%MK(8,1)*MK(8,2)/sqrt(MK(8,1)^2+MK(8,2)^2)+0.0000001;
X=2*wid*[x2,0,0];
Y=2*len*[0,y2,0];
posa2=[];
posb2=[];
newa=0;
newb=0;
for i =1:length(posa)
    if posa(i,1)<=sqrt(A^2+B^2)/2+wid*x2 && posa(i,1)>sqrt(A^2+B^2)/2-wid*x2  && posa(i,2)<len*y2 && posa(i,2)>=-len*y2
        posa2=[posa2;[posa(i,1)-sqrt(A^2+B^2)/2+wid*x2,posa(i,2)+len*y2,posa(i,3)]];
        newa=newa+1;
    end
end

for i =1:length(posb)
    if posb(i,1)<=sqrt(A^2+B^2)/2+wid*x2 && posb(i,1)>sqrt(A^2+B^2)/2-wid*x2 && posb(i,2)<len*y2 && posb(i,2)>=-len*y2
        posb2=[posb2;[posb(i,1)-sqrt(A^2+B^2)/2+wid*x2,posb(i,2)+len*y2,posb(i,3)]];
        newb=newb+1;
    end
end

posa=sortrows(posa2,2);

posb=sortrows(posb2,2);

%change to direct
cordinate=[X;Y;mc*MK(3,1:3)];
posa = ((cordinate*cordinate')\(cordinate*posa'))';
posb = ((cordinate*cordinate')\(cordinate*posb'))';

%change to direct


fout=fopen('example.txt','wt');
fprintf(fout,'%s\n','SnO'); %设置晶体名称
fprintf(fout,'%s\n','1.0');
fprintf(fout,'%s',' ' );

fprintf(fout,'%14.10f\t',X(1:2));
fprintf(fout,'%14.10f\n',X(3));
fprintf(fout,'%s',' ' );

fprintf(fout,'%14.10f\t',Y(1:2));
fprintf(fout,'%14.10f\n',Y(3));
fprintf(fout,'%s',' ' );

fprintf(fout,'%14.10f\t',mc*MK(3,1:2));
fprintf(fout,'%14.10f\n',mc*MK(3,3));

fprintf(fout,'%s','   ' );
fprintf(fout,'%s','Sn'); %设置原子
fprintf(fout,'%s','   ' );
fprintf(fout,'%s\n','O');
fprintf(fout,'%s','   ' );

fprintf(fout,'%d',newa);
fprintf(fout,'%s','   ' );
fprintf(fout,'%d\n',newb);
%此处可类似加原子
fprintf(fout,'%s\n','Direct');

[m,n]=size(posa);
for i=1:m
    for j=1:n
        if j==n
            fprintf(fout,'%14.10f\n',posa(i,j));
            
        else
            fprintf(fout,'%14.10f\t',posa(i,j));
            
        end
    end
end
[m,n]=size(posb);
for i=1:m
    for j=1:n
        if j==n
            fprintf(fout,'%14.10f\n',posb(i,j));
            
        else
            fprintf(fout,'%14.10f\t',posb(i,j));
          
        end
    end
end

fclose(fout);


