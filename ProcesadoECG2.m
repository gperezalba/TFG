function ProcesadoECG2(Name)
 
%Eliminar ruido de la señal
s1=Name';
s2=smooth(s1,10000);
ecgsmooth=s1-s2;
 
%Filtrar la señal
[C,L]=wavedec(ecgsmooth,8,'db4');
[d1,d2,d3,d4,d5,d6,d7,d8]=detcoef(C,L,[1,2,3,4,5,6,7,8]);
[thr,sorh,keepapp]=ddencmp('den','wv',ecgsmooth);
cleanecg=wdencmp('gbl',C,L,'db4',8,thr,sorh,keepapp);
 
%Obtener el ruido restandole la señal limpia
ruido = Name' - cleanecg - s2;
 
%ELIMINO MEDIA POR TRAMOS DE LA SEÑAL RUIDO
 
ruido_nm=ruido;
l0=length(ruido_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = ruido_nm(i:i+l);
            ruido_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if ruido_nm(k)>0
        positivos(n1) = ruido_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = ruido_nm(k);
        n2 = n2+1;
    end
end
mean(positivos)
mean(negativos)
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
%ACOTAR EN AMPLITUD LA SEÑAL RUIDO
k=1;
for i=1:length(ruido_nm)
    if abs(ruido_nm(i))<corte
        ruido_nmN(k)=ruido_nm(i);
        k=k+1;
    end
end
 
%OBTENER RUIDO ENTRE 1-128 
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
 
%Normalizar y pasar a binario
outNor=ruido_nmN/max(abs(ruido_nmN));
out256=127.5*outNor+127.5;
outDec=dec2bin(out256,8);
 
%Generar secuencia de bits
cont_zeros=0;
cont_unos=0;
for i=1:size(outDec,1)
  for j=1:4
    if outDec(i,j+4)=='1'
        cont_unos=cont_unos+1
        outFinal(4*(i-1)+j)='1';
    else
        cont_zeros=cont_zeros+1
        outFinal(4*(i-1)+j)='0';
    end
  end
end
 
 
%ESCRIBIR SALIDA EN FICHERO
fileID=fopen('filename.e','w');
fprintf(fileID,outFinal);
fclose(fileID);
 
end