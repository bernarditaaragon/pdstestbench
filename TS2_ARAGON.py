# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 21:41:34 2026

@author: Aragon
"""

#importacion de modulos
import numpy as np
import matplotlib.pyplot as plt 

#definiciones
N= 1000 #cantidad de muestras
fs=1000 #frecuencia de muestreo
ts= 1/fs #tiempo entre muestras
B= 4 #cantidad de bits
vf=2 #rango
f0= fs/N #frecuencia de la senoide
P=1 #potencia de la senoide
vmax= np.sqrt(2*P)
dc = 0
ph = 0
kn = 1  # factor de escala del ruido analogico


#funciones
def mi_funcion_sen( vmax, dc, ff, ph , N, fs):
    
    tt= np.arange(0,N) / fs
    xx= dc + vmax*np.sin(2 * np.pi * ff * tt + ph)

    return(tt, xx)

#INCISO A
tt, xx = mi_funcion_sen(vmax, dc, f0, ph, N, fs)
q = (2*vf)/(2**B) #paso de cuantizacion
Pq= (q**2)/12 #potencia de cuantizacion
Pr= kn*Pq #potencia de ruido
desvio = np.sqrt(Pr)
ruido = np.random.normal(0, desvio, N)
#senal
ruido_xx = xx + ruido
ruido_xx_q = q * np.round(ruido_xx / q) #cuantizacion


#grafico
plt.figure(figsize=(12 , 4))
plt.plot(tt, ruido_xx_q, label="sQ =QB,Vf{sR} (ADC out)")
plt.plot(tt, ruido_xx, 'g:o', markersize=2, linewidth=0.8, label= "sR = s + n (ADC in)")
plt.plot(tt, xx, ':',linewidth=1.2, label= "s (analog)")
plt.title(f'Señal muestreada por un ADC de {B} bits - ±Vr = {vf:.1f} V - q = {q:.3f}V')
plt.xlabel("tiempo [segundos]")
plt.ylabel("Amplitud [V]")
plt.legend()
plt.show()

#vector de frecuencias
ff_vec = np.fft.fftfreq(N, 1/fs)
ff_vec = ff_vec[:N//2] #me quedo con las frecuencias positivas hasta Nyquist

#espectros
#senal
XX = np.fft.fft(xx)/N
XX= XX[:N//2] #mismo tamano que el vector de frecuencias
PXX_mod = 2 * (np.abs(XX)**2)
PXX_mod_db = 10 * np.log10(PXX_mod + 1e-12)
#senal con ruido
R_XX = np.fft.fft(ruido_xx)/N
R_XX= R_XX[:N//2] #mismo tamano que el vector de frecuencias
PR_XX_mod = 2 * (np.abs(R_XX)**2)
PR_XX_mod_db = 10 * np.log10(PR_XX_mod + 1e-12)
#senal con ruido cuantizada
R_XX_Q = np.fft.fft(ruido_xx_q)/N
R_XX_Q= R_XX_Q[:N//2] #mismo tamano que el vector de frecuencias
PR_XX_Q_mod = 2 * (np.abs(R_XX_Q)**2)
PR_XX_Q_mod_db = 10 * np.log10(PR_XX_Q_mod + 1e-12)

#pisos de ruido
piso_analogico = 10 * np.log10(Pr/(N/2) + 1e-12)
piso_digital = 10 * np.log10(Pq/(N/2) + 1e-12)

print("Piso analogico:", piso_analogico, "dB")
print("Piso digital:", piso_digital, "dB")

#graficos
plt.figure(figsize=(12,5))

plt.plot(ff_vec, PR_XX_Q_mod_db, color='b', linewidth=1.2,
         label="sQ = QB,Vf{sR} (ADC out)")

plt.plot(ff_vec, PXX_mod_db, color='orange', linestyle=':', linewidth=1.2,
         label="s (analog)")

plt.plot(ff_vec, PR_XX_mod_db, color='g', linestyle=':',linewidth=1.2,
         label="sR = s + n (ADC in)")

plt.axhline(piso_analogico, color='r', linestyle='--',linewidth=2.5,
            label=f"Piso analógico = {piso_analogico:.1f} dB")

plt.axhline(piso_digital, color='c', linestyle=':', linewidth=1.5,
            label=f"Piso digital = {piso_digital:.1f} dB")

plt.title(f'Señal muestreada por un ADC de {B} bits - ±Vr = {vf:.1f} V - q = {q:.3f} V')
plt.xlabel("Frecuencia [Hz]")
plt.ylabel("Densidad de Potencia [dB]")
plt.xlim(0, fs/2)
plt.ylim(-85, 10)
plt.legend()
plt.show()


# histograma
eq = ruido_xx_q - ruido_xx #error ruido de cuantizacion

#grafico
plt.figure(figsize=(12,5))

plt.hist(eq, bins=20, density=True,
         color='b', alpha=0.6)

#distribucion uniforme teorica
plt.plot([-q/2, -q/2, q/2, q/2],
         [0, 1/q, 1/q, 0],
         color='r', linestyle='--', linewidth=2)

plt.title(f'Ruido de cuantización para {B} bits - ±Vr = {vf:.1f} V - q = {q:.3f} V')
plt.xlim(-q/2 - 0.02, q/2 + 0.02)
plt.show()


#INCISO B
#para no tener que copiar el codigo repetidas veces voy a hacer un for que pase por todos los valores
valores_B = [4, 8, 16]
valores_kn = [0.1, 1, 10]

for B in valores_B:
    for kn in valores_kn:
        if B == 4 and kn == 1: #este caso ya lo hice en el inciso A
           continue
        print("B =", B, "- kn =", kn)
        
        tt, xx = mi_funcion_sen(vmax, dc, f0, ph, N, fs)
        q = (2*vf)/(2**B) #paso de cuantizacion
        Pq= (q**2)/12 #potencia de cuantizacion
        Pr= kn*Pq #potencia de ruido
        desvio = np.sqrt(Pr)
        ruido = np.random.normal(0, desvio, N)
        #senal
        ruido_xx = xx + ruido
        ruido_xx_q = q * np.round(ruido_xx / q) #cuantizacion
        #vector de frecuencias
        ff_vec = np.fft.fftfreq(N, 1/fs)
        ff_vec = ff_vec[:N//2] #me quedo con las frecuencias positivas hasta Nyquist
        #espectros
        #senal
        XX = np.fft.fft(xx)/N
        XX= XX[:N//2] #mismo tamano que el vector de frecuencias
        PXX_mod = 2 * (np.abs(XX)**2)
        PXX_mod_db = 10 * np.log10(PXX_mod + 1e-30)
        #senal con ruido
        R_XX = np.fft.fft(ruido_xx)/N
        R_XX= R_XX[:N//2] #mismo tamano que el vector de frecuencias
        PR_XX_mod = 2 * (np.abs(R_XX)**2)
        PR_XX_mod_db = 10 * np.log10(PR_XX_mod + 1e-30)
        #senal con ruido cuantizada
        R_XX_Q = np.fft.fft(ruido_xx_q)/N
        R_XX_Q= R_XX_Q[:N//2] #mismo tamano que el vector de frecuencias
        PR_XX_Q_mod = 2 * (np.abs(R_XX_Q)**2)
        PR_XX_Q_mod_db = 10 * np.log10(PR_XX_Q_mod + 1e-30)

        #pisos de ruido
        piso_analogico = 10 * np.log10(Pr/(N/2))
        piso_digital = 10 * np.log10(Pq/(N/2))

        print("Piso analogico:", piso_analogico, "dB")
        print("Piso digital:", piso_digital, "dB")
        
        eq = ruido_xx_q - ruido_xx #error ruido de cuantizacion

        #graficos
        plt.figure(figsize=(12,14))
        
        #senal
        plt.subplot(3,1,1)
        plt.plot(tt, ruido_xx_q, label="sQ =QB,Vf{sR} (ADC out)")
        plt.plot(tt, ruido_xx, 'g:o', markersize=2, linewidth=0.8, label= "sR = s + n (ADC in)")
        plt.plot(tt, xx, ':',linewidth=1.2, label= "s (analog)")
        plt.title(f'Señal muestreada por un ADC de {B} bits - kn = {kn} - ±Vr = {vf:.1f} V - q = {q:.6f} V')
        plt.xlabel("tiempo [segundos]")
        plt.ylabel("Amplitud [V]")
        plt.legend()
        
        # espectros
        plt.subplot(3,1,2)
        plt.plot(ff_vec, PR_XX_Q_mod_db, color='b', linewidth=1.2,
                 label="sQ = QB,Vf{sR} (ADC out)")

        plt.plot(ff_vec, PXX_mod_db, color='orange', linestyle=':', linewidth=1.2,
                 label="s (analog)")

        plt.plot(ff_vec, PR_XX_mod_db, color='g', linestyle=':',linewidth=1.2,
                 label="sR = s + n (ADC in)")

        plt.axhline(piso_analogico, color='r', linestyle='--',linewidth=2.5,
                    label=f"Piso analógico = {piso_analogico:.1f} dB")

        plt.axhline(piso_digital, color='c', linestyle=':', linewidth=1.5,
                    label=f"Piso digital = {piso_digital:.1f} dB")

        plt.title(f'Señal muestreada por un ADC de {B} bits - kn = {kn} - ±Vr = {vf:.1f} V - q = {q:.6f} V')
        plt.xlabel("Frecuencia [Hz]")
        plt.ylabel("Densidad de Potencia [dB]")
        plt.xlim(0, fs/2)
        plt.ylim(min(piso_analogico, piso_digital)-15, 10)
        plt.legend()
        
        # histograma
        plt.subplot(3,1,3)
        plt.hist(eq, bins=20, density=True,
                 color='b', alpha=0.6)

        #distribucion uniforme teorica
        plt.plot([-q/2, -q/2, q/2, q/2],
                 [0, 1/q, 1/q, 0],
                 color='r', linestyle='--', linewidth=2)

        plt.title(f'Ruido de cuantización para {B} bits - kn = {kn} - ±Vr = {vf:.1f} V - q = {q:.6f} V')
        plt.xlim(-0.6*q, 0.6*q)

        plt.tight_layout()
        plt.show()
        