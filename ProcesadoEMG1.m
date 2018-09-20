function ProcesadoEMG1(signal)
 
%ELIMINAR MEDIA POR TRAMOS DE LA SEÑAL
signal_nm=signal;
l0=length(signal_nm);
l=l0-1;
 
while (l>100)
    
    for i=1:l:l0
        if(i+l)<l0
            tramo = signal_nm(i:i+l);
            signal_nm(i:i+l) = tramo - mean(tramo);
        end
    end
    l = l/2;
end
 
l=l*2;
signal_nm(l0-l:l0)=signal_nm(l0-l:l0)-mean(signal_nm(l0-l:l0));
 
n1=1;
n2=1;
positivos=0;
negativos=0;
for k=1:l0
    if signal_nm(k)>0
        positivos(n1) = signal_nm(k);
        n1 = n1+1;
    else
        negativos(n2) = signal_nm(k);
        n2 = n2+1;
    end
end
mean(positivos)
mean(negativos)
corte = (mean(positivos)+abs(mean(negativos)))/2;
 
%ACOTAR EN AMPLITUD LA SEÑAL 
k=1;
for i=1:length(signal_nm)
    if abs(signal_nm(i))<corte
        signal_nmN(k)=signal_nm(i);
        k=k+1;
    end
end
 
%Normalizar y paso a binario
outNor=signal_nmN/max(abs(signal_nmN));
out256=127.5*outNor+127.5;
outDec=dec2bin(out256,8);
 
%GENERAR SECUENCIA DE BITS
cont_zeros=0;
cont_unos=0;
for i=1:size(outDec,1)
  for j=1:4
    if outDec(i,j+4)=='1'
        cont_unos=cont_unos+1;
        outFinal(4*(i-1)+j)='1';
    else
        cont_zeros=cont_zeros+1;
        outFinal(4*(i-1)+j)='0';
    end
  end
end
 
%ESCRIBO SALIDA EN FICHERO
fileID=fopen('filename.e','w');
fprintf(fileID,outFinal);
fclose(fileID);
 
end
