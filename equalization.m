function  [GN_stand MIG_Max_cellj0]= GN_stand(data, PPL, MX, T, T0)

%% 2.1 Determination of particle size standard GN_stand
MIG_Max=[];NNUM=[];
for i5=1:size(data,2)
    
    MIG_max=cell2mat(PPL(i5));
    MM=MIG_max(:,end);
    for i6=1:length(MM)
        NUM(i6)=sum(MM([1:i6],:))-(i6-1);
    end
    
    NUM=NUM(:,1:length(MM));
    
    TT=[0 T];
    TTT=unique([T NUM]);
    
    IGG_max=[];
    for i8=1:length(TT)-1
        label=find(NUM>TT(i8)&NUM<=TT(i8+1));
        
        if length(label)==0
            label1=find(NUM>TT(i8+1));
            label=find(NUM==NUM(label1(1)));
        elseif NUM(label(end))~=TT(i8+1)
            label=[label label(end)+1];
        else
            label=label;
        end
        
        G_max=MIG_max(label,:); 
        
        D_i8=find(TTT>TT(i8)&TTT<=TT(i8+1));
        
        if i8==1
            for i9=1:length(D_i8)
                if i9==1
                    Label_i9(i9)=TTT(D_i8(i9));
                else
                    Label_i9(i9)=TTT(D_i8(i9))-TTT(D_i8(i9)-1)+1;
                end
            end
            
        else
            
            for i9=1:length(D_i8)
                Label_i9(i9)=TTT(D_i8(i9))-TTT(D_i8(i9)-1)+1;
            end
        end
        
        GG_max=[G_max Label_i9'];
        GG_max(:,end-1)=[];
        Gg_max=mat2cell(GG_max,size(GG_max,1),4);
        IGG_max=[IGG_max Gg_max];
        clearvars -except D lambda_max y  n data Obje X FIG MObje_max MX Alpha_opt  MIG_opt T0 T MIG_max NUM TT TTT IGG_max MIG_Max NNUM Att Att1 PPL
    end
    MIG_Max=[MIG_Max IGG_max']; 
    NNUM=[NNUM mat2cell(NUM,1,size(NUM,2))];
end

Yy_i11=[]; 
for i11=1:size(MIG_Max,2) 
    y_i11=MX(:,i11);
    
    Yy_i12=[];
    for i12=1:size(MIG_Max,1) 
        if i12==1
            y_i12=y_i11([1:T0]);
        else
            y_i12=y_i11([(i12-1)*T0:i12*T0]);
        end
        MIG_Max_12=MIG_Max{i12,i11};
        
        Yy_i13=[];
        for i13=1:size(MIG_Max_12,1) 
            %%实际值
            if i13==1
                yy_i13=y_i12([1:MIG_Max_12(i13,4)]);
            else
                label_begain=sum(MIG_Max_12([1:i13-1],4))-length([1:i13-1])+1;
                label_end=label_begain+(MIG_Max_12(i13,4)-1);
                yy_i13=y_i12([label_begain:label_end]);
            end
            Yy_i13=[Yy_i13 mat2cell(yy_i13,size(yy_i13,1),size(yy_i13,2))];
        end
        Yy_i12=[Yy_i12 mat2cell(Yy_i13,size(Yy_i13,1),size(Yy_i13,2)) ];
    end
    Yy_i11=[Yy_i11 Yy_i12'];
end


MIG_Max_14=[];
for i14=1:size(Yy_i11,2) 
    MIG_Max_15=[];
    for i15=1:size(Yy_i11,1) 
        MIG_Max_45=MIG_Max{i15,i14};
        K_i16=[]; B_i16=[];
        for i16=1:size(MIG_Max_45,1) 
            k_i16=(Yy_i11{i15,i14}{1,i16}(end)-Yy_i11{i15,i14}{1,i16}(1))/(length(Yy_i11{i15,i14}{1,i16})-1);
            b_i16=Yy_i11{i15,i14}{1,i16}(1)-k_i16;
            K_i16=[K_i16 k_i16]; B_i16=[B_i16 b_i16];
        end
        MIG_Max_45(:,2)=B_i16';
        MIG_Max_15=[MIG_Max_15 mat2cell(MIG_Max_45, size(MIG_Max_45,1),size(MIG_Max_45,2))];
    end
    MIG_Max_14=[MIG_Max_14 MIG_Max_15'];
end


MIG_Max=MIG_Max_14;

%%%Determine granularity criteria
for i10=1:size(MIG_Max,1)
    for i11=1:size(MIG_Max,2)
        GN(i10,i11)=size(cell2mat(MIG_Max(i10,i11)),1);
    end
end

GN_stand=max(max(GN));  

%% 2.2 Granularity reference
G_ref_i12=[];
for i12=1:size(data,2)
    G_ref=[];
    for i13=2:size(data,1)
        k_i12=data(i13,i12)-data(i13-1,i12);
        b_i12=data(i13,i12)-2*k_i12;
        %At this point, two neighboring points are used for linear fitting, so the deviation is always 0 and the time span is all 2
        g_ref=[k_i12  b_i12 0 2];
        G_ref=[G_ref g_ref];
    end
    G_ref=(reshape(G_ref,[4,size(data,1)-1]))';
    
    GG_ref=[];
    for i14=1:length(T)
        if i14==1
            G_ref14=G_ref(1:T(i14)-1,:);
        else
            G_ref14=G_ref(T(i14-1):T(i14)-1,:);
        end
        G_refi=mat2cell(G_ref14,size(G_ref14,1),size(G_ref14,2));
        GG_ref=[GG_ref G_refi];
    end
    G_ref_i12=[G_ref_i12 GG_ref']; 
end

%% 2.3 Granularity equalization (granularity in all cycles)）

MIG_Max_cellj0=[];
for j0=1:size(data,2) 
    
    y_j0=data(:,j0);
    
    MIG_Max_cellj1=[];
    for j1=1:size(MIG_Max,1) 
        
        if j1==1
            y_j1=y_j0([1:T0]);
        else
            y_j1=y_j0([(j1-1)*T0:j1*T0]);
        end
        
        MIG_Max_01=MIG_Max{j1,j0};
        
        while size(MIG_Max_01,1)<GN_stand
            Mat_j2=[]; 
            Delta_j2=[]; 
            YY_j1=[];
            for j2=1:size(MIG_Max_01,1) 
                
                if j2==1
                    yy_j1=y_j1([1:MIG_Max_01(j2,4)]);
                else
                    label_begain=sum(MIG_Max_01([1:j2-1],4))-length([1:j2-1])+1;
                    label_end=label_begain+(MIG_Max_01(j2,4)-1);
                    yy_j1=y_j1([label_begain:label_end]);
                end
                
                %%fitted value
                fit_value=MIG_Max_01(j2,1).*[1:MIG_Max_01(j2,4)]+MIG_Max_01(j2,2);
                
                %%Updated standard deviationdelta
                delta_j2=sqrt(sum((yy_j1'-fit_value).^2)/length(fit_value));
                
                %%match value
                mat=sum(exp(-(yy_j1'-fit_value).^2))/delta_j2; 
                
                Delta_j2=[Delta_j2 delta_j2];
                Mat_j2=[Mat_j2 mat];
                YY_j1=[YY_j1 mat2cell(yy_j1,size(yy_j1,1),size(yy_j1,2))];
            end
            MIG_Max_01(:,3)=Delta_j2'; 
            
            [px,v]=sort(Mat_j2); 
            
            lab_01=MIG_Max_01(:,end);
            
            while lab_01(v(1))<=2
                loc=find(lab_01(v)>2)
                v(1)=v(loc(1));
            end
             
            Distance=[];MIG_reshape=[];
            for j3=2:lab_01(v(1))-1
                point_first=MIG_Max_01(v(1),1)+MIG_Max_01(v(1),2); 
                point_end=MIG_Max_01(v(1),1)*lab_01(v(1))+MIG_Max_01(v(1),2); 
                 
%                 
                point_turn=YY_j1{1,v(1)}(j3,:); 
                
                k1=(point_turn-point_first)/(j3-1); 
                b1=point_first-k1;
                L1=j3;
                deta1= sqrt(sum((YY_j1{1,v(1)}(1:j3,:)-(k1.*[1:L1]+b1)').^2)/L1);
                G1=[k1 b1 deta1 L1];
                
                %%%%%%%%%%
                
                k2=(point_end-point_turn)/(lab_01(v(1))-j3); 
                b2=point_turn-k2;
                L2=lab_01(v(1))-j3+1;
                deta2= sqrt(sum((YY_j1{1,v(1)}(j3:end,:)-(k2.*[1:L2]+b2)').^2)/L2);
                G2=[k2 b2 deta2 L2];
                
                
                MIG_Max_01j3=[(MIG_Max_01(1:v(1)-1,:))' G1' G2'  (MIG_Max_01(v(1)+1:end,:))']'; 
                
                %%%%%%%%%%% 
                MIG_ref=G_ref_i12{j1,j0};
                
                for j4=1:size(MIG_ref,1)
                    for j5=1:size(MIG_Max_01j3,1)                 
                        d(j4,j5)=sum(MIG_ref(j4,:).*MIG_Max_01j3(j5,:))/(sqrt(sum((MIG_ref(j4,:)).^2))*sqrt(sum((MIG_Max_01j3(j5,:)).^2)));
                    end 
                end
                dd=1-d;
               
                %Accumulate the distances to get the DP matrix
                DP=zeros(size(dd,1),size(dd,2));
                DP(1,1)=dd(1,1);
                for j6=2:size(dd,1)
                    DP(j6,1)=dd(j6,1)+DP(j6-1,1);
                end
                for j7=2:size(dd,2)
                    DP(1,j7)=dd(1,j7)+DP(1,j7-1);
                end
                for j8=2:size(dd,1)
                    for j9=2:size(dd,2)
                        DP(j8,j9)=dd(j8,j9)+GetMin(DP(j8-1,j9),DP(j8,j9-1),DP(j8-1,j9-1));
                    end
                end
                Distance=[Distance DP(end,end)];
                MIG_reshape=[MIG_reshape mat2cell(MIG_Max_01j3,size(MIG_Max_01j3,1),size(MIG_Max_01j3,2))];
            end
            MIG_Max_01=MIG_reshape{find(Distance==min(Distance))};
        end
        
        MIG_Max_01
        
        MIG_Max_cellj1=[MIG_Max_cellj1 mat2cell(MIG_Max_01,size(MIG_Max_01,1),size(MIG_Max_01,2)) ];  
    end
    
    MIG_Max_cellj0=[MIG_Max_cellj0 MIG_Max_cellj1']; 
end
