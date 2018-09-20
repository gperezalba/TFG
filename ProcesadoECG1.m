function ProcesadoECG1(Name)
 
%Ruido de la señal
s1=Name';
s2=smooth(s1,10000);
ecgsmooth=s1-s2;
 
%Filtrado de la señal
[C,L]=wavedec(ecgsmooth,8,'db4');
[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
[thr,sorh,keepapp]=ddencmp('den','wv',ecgsmooth);
cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);
 
%Obtener señal ruido
ruido = Name' - cleanecg - s2;
 
%Localizar picos R en la señal
[ondasRECG, nR] = findpeaks(ecgsmooth,'MinpeakDistance',150); 
nR=nR(2:(length(nR)-1));
 
%Localizar picos Q en la señal
nQ=zeros(1,length(nR));
for k = 1:(length(nR)) 
    n1 = nR(k) - 50; 
    n2 = nR(k);
    ecgPart = cleanecg(n1:n2);
    A = ecgPart;
    [fila, columnaQ] = findpeaks (-A,'MinpeakDistance',50); 
    if length(columnaQ)>1
        columnaQ=columnaQ(length(columnaQ));
    end
    if length(columnaQ)==0
        columnaQ=40;
    end
    nQ(k) = [nR(k) - length(A) + columnaQ]; 
end
nQ = nQ';
 
%Localizar picos S en la señal
nS=zeros(1,length(nR));
for j = 1:length(nR)
    m1 = nR(j) + 50;
    m2 = nR(j);
    ecgPart = cleanecg(m2:m1);
    B = ecgPart;
    [fila, columnaS] = findpeaks (-B,'MinpeakDistance',50);
    if length(columnaS)>1
        columnaS=columnaS(length(columnaS));
    end
    if length(columnaS)==0
        columnaS=10;
    end
    nS(j) = [nR(j) + columnaS];
end
nS = nS';
 
%Localizar picos T en la señal
nT=zeros(1,length(nR));
for i = 1:length(nR)
    n1 = nR(i) + 100 + 20;
    if n1>length(cleanecg) 
        n1=length(cleanecg); 
    end
    n2 = nR(i)+20;
    ecgPart = cleanecg(n2:n1);
    B = ecgPart;
    [fila, columnaT] = findpeaks (B,'MinpeakDistance',75);
    if length(columnaT)>1
        columnaT=columnaT(length(columnaT));
    end
    if length(columnaT)==0
        columnaT=31;
    end
    nT(i) = [nR(i) + 20 + columnaT];
end
nT = nT';
 
%Localizar picos P en la señal
nP=zeros(1,length(nR));
for l = 1:length(nR)
    n1 = nR(l) - 75 - 20;
    n2 = nR(l)-20;
    ecgPart = cleanecg(n1:n2);
    A = ecgPart;
    [fila, columnaP] = findpeaks (A,'MinpeakDistance',75);
    if length(columnaP)>1
        columnaP=columnaP(length(columnaP));
    end
    if length(columnaP)==0
        columnaP=10;
    end
    nP(l) = [nR(l) - 20 - length(A) + columnaP];
end
nP = nP';
 
%GENERAR SECUENCIA RR RQ RS RP RT ETC
    rr = zeros(1, (length(nR)-1));
    rq = zeros(1, (length(nR)-1));
    rs = zeros(1, (length(nR)-1));
    rp = zeros(1, (length(nR)-1));
    rt = zeros(1, (length(nR)-1));
 
     
for z = 1:(length(nR)-1)
    rr(z) = nR(z+1) - nR(z) + cleanecg(nR(z));
    rq(z) = nR(z) - nQ(z) + cleanecg(nQ(z));
    rs(z) = nS(z) - nR(z) + cleanecg(nS(z));
    rp(z) = nR(z) - nP(z) + cleanecg(nP(z));
    rt(z) = nT(z) - nR(z) + cleanecg(nT(z));
 
end
 
 
%ELIMINAR MEDIA POR TRAMOS DE CADA SEÑAL Y ACOTAR EN AMPLITUD
rr_nm=rr;
l0=length(rr_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = rr_nm(i:i+l);
            rr_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
rr_nm(l0-l:l0)=rr_nm(l0-l:l0)-mean(rr_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if rr_nm(k)>0
        positivos(n1) = rr_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = rr_nm(k);
        n2 = n2+1;
    end
end
mean(positivos);
mean(negativos);
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
k=1;
for i=1:length(rr_nm)
    if abs(rr_nm(i))<corte
        rr_nmN(k)=rr_nm(i);
        k=k+1;
    end
end
 
rq_nm=rq;
l0=length(rq_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = rq_nm(i:i+l);
            rq_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
rq_nm(l0-l:l0)=rq_nm(l0-l:l0)-mean(rq_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if rq_nm(k)>0
        positivos(n1) = rq_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = rq_nm(k);
        n2 = n2+1;
    end
end
mean(positivos);
mean(negativos);
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
k=1;
for i=1:length(rq_nm)
    if abs(rq_nm(i))<corte
        rq_nmN(k)=rq_nm(i);
        k=k+1;
    end
end
 
 
rs_nm=rs;
l0=length(rs_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = rs_nm(i:i+l);
            rs_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
rs_nm(l0-l:l0)=rs_nm(l0-l:l0)-mean(rs_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if rs_nm(k)>0
        positivos(n1) = rs_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = rs_nm(k);
        n2 = n2+1;
    end
end
mean(positivos);
mean(negativos);
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
k=1;
for i=1:length(rs_nm)
    if abs(rs_nm(i))<corte
        rs_nmN(k)=rs_nm(i);
        k=k+1;
    end
end
 
 
rp_nm=rp;
l0=length(rp_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = rp_nm(i:i+l);
            rp_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
rp_nm(l0-l:l0)=rp_nm(l0-l:l0)-mean(rp_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if rp_nm(k)>0
        positivos(n1) = rp_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = rp_nm(k);
        n2 = n2+1;
    end
end
mean(positivos);
mean(negativos);
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
k=1;
for i=1:length(rp_nm)
    if abs(rp_nm(i))<corte
        rp_nmN(k)=rp_nm(i);
        k=k+1;
    end
end
 
 
rt_nm=rt;
l0=length(rt_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = rt_nm(i:i+l);
            rt_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
rt_nm(l0-l:l0)=rt_nm(l0-l:l0)-mean(rt_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if rt_nm(k)>0
        positivos(n1) = rt_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = rt_nm(k);
        n2 = n2+1;
    end
end
mean(positivos);
mean(negativos);
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
k=1;
for i=1:length(rt_nm)
    if abs(rt_nm(i))<corte
        rt_nmN(k)=rt_nm(i);
        k=k+1;
    end
end
 
%NORMALIZAR SEÑALES
 
rpNor=rp_nmN/max(abs(rp_nmN));
rp256=127.5*rpNor+127.5;
 
rrNor=rr_nmN/max(abs(rr_nmN));
rr256=127.5*rrNor+127.5;
 
rqNor=rq_nmN/max(abs(rq_nmN));
rq256=127.5*rqNor+127.5;
 
rsNor=rs_nmN/max(abs(rs_nmN));
rs256=127.5*rsNor+127.5;
 
rtNor=rt_nmN/max(abs(rt_nmN));
rt256=127.5*rtNor+127.5;
 
 
%OBTENER RUIDO ENTRE 1-128 Y 1-5
 
k=1;
for i=1:length(ruido)
    if abs(ruido(i))<12
        ruidoN(k)=ruido(i);
        k=k+1;
    end
end
 
ruido1=ruidoN/max(abs(ruidoN));
ruido128=63*ruido1+64;
ruido128=fix(ruido128);
 
ruido09=2.5*ruido1+3.5;
ruido09=fix(ruido09);
 
 
%CONSTRUIR SEÑAL SALIDA COMBINANDO RP RQ RR RS RT
 
salida=[];
c=1;
 
while(length(salida)<1000000)
    
    inicio=ruido128(c);
    fin=inicio+ruido128(c+1);
    
    if(ruido09(c)==1)&&(fin<(length(rp256)))
         trozo=rp256(inicio:fin);
    end
    if(ruido09(c)==2)&&(fin<(length(rq256)))
         trozo=rq256(inicio:fin);
    end
    if(ruido09(c)==3)&&(fin<(length(rr256)))
         trozo=rr256(inicio:fin);
    end
    if(ruido09(c)==4)&&(fin<(length(rs256)))
         trozo=rs256(inicio:fin);
    end
    if(ruido09(c)==5)&&(fin<(length(rt256)))
         trozo=rt256(inicio:fin);
    end
    
    salida=[salida trozo];
    c=c+1;
    
end
 
%ELIMINAR MEDIA POR TRAMOS Y ACOTAR EN AMPLITUD A LA SEÑAL DE SALIDA.
salida_nm=salida;
l0=length(salida_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = salida_nm(i:i+l);
            salida_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
salida_nm(l0-l:l0)=salida_nm(l0-l:l0)-mean(salida_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if salida_nm(k)>0
        positivos(n1) = salida_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = salida_nm(k);
        n2 = n2+1;
    end
end
mean(positivos);
mean(negativos);
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
k=1;
for i=1:length(salida_nm)
    if abs(salida_nm(i))<corte
        salida_nmN(k)=salida_nm(i);
        k=k+1;
    end
end
 
 
%NORMALIZAR Y PASAR A BINARIO
outNor=salida_nmN/max(abs(salida_nmN));
out256=127.5*outNor+127.5;
outDec=dec2bin(out256,8);
 
%GENERAR SECUENCIA DE SALIDA DE BITS
%BUCLE EQUILIBRAR 1 Y 0
cont_zeros=0;
cont_unos=0;
ctramo1=0;
ctramo0=0;
insertados=0;
 
insert=fix(ruido128(1)*2);
cr=1;
 
for i=1:size(outDec,1)
    
  for j=1:4
    if outDec(i,j+4)=='1'
        cont_unos=cont_unos+1;
        ctramo1=ctramo1+1;
        insert=insert-1;
        outFinal(4*(i-1)+j+insertados)='1';
    else
        cont_zeros=cont_zeros+1;
        ctramo0=ctramo0+1;
        insert=insert-1;
        outFinal(4*(i-1)+j+insertados)='0';
        
    end
    if insert==0
            if ctramo1>ctramo0
                outFinal(4*(i-1)+j+1+insertados)='0';
                insertados=insertados+1;
                cont_zeros=cont_zeros+1;
            else
                outFinal(4*(i-1)+j+1+insertados)='1';
                insertados=insertados+1;
                cont_unos=cont_unos+1;
            end
        cr=cr+1;
        insert=fix(ruido128(cr)*2);
        ctramo0=0;
        ctramo1=0;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%ESCRIBIR SALIDA EN FICHERO
fileID=fopen('filename.e','w');
fprintf(fileID,outFinal);
fclose(fileID);
 
end