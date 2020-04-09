from scipy import signal
import csv
import matplotlib.pyplot as plt
import numpy as np

with open('/home/pruthvirg/Downloads/h.csv', 'r') as f:
  reader = csv.reader(f)
  your_list = list(reader)

w, h = signal.freqz(your_list)

fig, ax1 = plt.subplots()
ax1.set_title('Digital filter frequency response')
ax1.plot(w, abs(h), 'b')
ax1.set_ylabel('Amplitude [dB]', color='b')
ax1.set_xlabel('Frequency [rad/sample]')
ax2 = ax1.twinx()
angles = np.unwrap(np.angle(h))
ax2.plot(w, angles, 'g')
ax2.set_ylabel('Angle (radians)', color='g')
ax2.grid()
ax2.axis('tight')
plt.show()


x = np.linspace(1,2**10,2**10)
y = np.cos(0.2*(np.pi)*x) + np.cos(0.85*(np.pi)*x)
plt.plot(abs(y))
plt.title('Plot of input signal')
plt.xlabel('n')
plt.ylabel('x[n]')
plt.grid(True)
plt.show()

your_list_float = [float(i[0]) for i in your_list]



plt.plot(your_list_float)
plt.title('Plot of h[n]')
plt.xlabel('n')
plt.ylabel('|h[n]|')
plt.grid(True)
plt.show()



z = np.convolve(your_list_float,y)
plt.plot(z)
plt.title('Output of simple convolution')
plt.xlabel('n')
plt.ylabel('|y_1[n]|')
plt.grid(True)
plt.show()

a_1 = np.concatenate([ your_list_float,np.zeros(len(x)-len(your_list_float)) ])
y1=np.fft.ifft(np.fft.fft(y) * np.fft.fft(a_1))
plt.plot(y1)
plt.title('Output of convolution with zero padding')
plt.xlabel('n')
plt.ylabel('|y_2[n]|')
plt.grid(True)
plt.show()



N = len(x)+len(your_list_float)-1
fil = np.concatenate([your_list_float,np.zeros(N-len(your_list_float))])
y_new = np.concatenate([y,np.zeros(N-len(y))])
y2=np.fft.ifft(np.fft.fft(y_new) * np.fft.fft(fil))
plt.plot(y2)
plt.title('Output of convolution with padding and circular convolution')
plt.xlabel('n')
plt.ylabel('|y_3[n]|')
plt.grid(True)
plt.show()



with open('/home/pruthvirg/Downloads/x1.csv', 'r') as f:
  reader = csv.reader(f)
  z_seq = list(reader)

z_seq_shifted = z_seq[-5:]+ z_seq[:-5]
z_seq_shifted = [complex((i[0].replace("i", "j")).strip()) for i in z_seq_shifted]
z_seq = [complex((i[0].replace("i", "j")).strip()) for i in z_seq]
z_1 = np.correlate(z_seq,z_seq_shifted,"full")
z_1_mag = [abs(i) for i in z_1] 
plt.plot(z_1_mag)
plt.title('Cross-correlation output')
plt.xlabel('n')
plt.ylabel('|y_3[n]|')
plt.grid(True)
plt.show()



