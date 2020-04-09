from sympy import *      
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal as sp

init_session      
s = symbols('s')  



def lowpass(R1,R2,C1,C2,G,Vi):
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0], [0,-G,G,1],[-(1/R1)-(1/R2)-(s*C1),1/R2,0,s*C1]])
    b=Matrix([0,0,0,-Vi/R1])
    V=A.inv()*b
    return (A,b,V)


def highpass(R1,R3,C1,C2,G,Vi):
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],[-(s*R3*C2)/(1+s*R3*C2),1,0,0], [0,-G,G,1],[-(s*C1)-(s*C2)-(1/R1),s*C2,0,1/R1]])
    b=Matrix([0,0,0,-Vi*s*C1])
    V=A.inv()*b
    return (A,b,V)


def get_response_freq_domain(h,filterType):
    if (filterType == 'l'):              
        A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,h)
    elif (filterType == 'h'):           
        A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,h)
    Vo=V[3]                              

    
    w=np.logspace(0,8,801)
    ss=1j*w
    hf=lambdify(s,Vo,'numpy')
    v=hf(ss)
    return Vo,w,abs(v),np.angle(v)       



def get_numpy_array_from_Poly(num,den):   
    isFloat = False
    try:
        num = Poly(num).all_coeffs()      
    except GeneratorsNeeded:
        num = num                        
        isFloat = True
    den = Poly(den).all_coeffs()          
    
    
    den2 = []
    num2 = []
    for i in den:
        den2.append(float(i))
    den2 = np.array(den2)
    
    if (isFloat): 
        num2 = num
    else:
        for i in num:
            num2.append(float(i))
        num2 = np.array(num2)
    return num2,den2


def get_output_time_domain(Y,t,steps):
    simplY = simplify(Y)     
    num = fraction(simplY)[0]   
    den = fraction(simplY)[1]   
    num2,den2 = get_numpy_array_from_Poly(num,den)  
    
    num2 = np.poly1d(num2) 
    den2 = np.poly1d(den2)
    
   
    Y = sp.lti(num2,den2)
    t = np.linspace(0.0,t,steps)
    t,y=sp.impulse(Y,None,t)
    return t,y



def get_output_with_lsim(H,x,t): 

 
    simplH = simplify(H)
    num = fraction(simplH)[0]
    den = fraction(simplH)[1]
    num2,den2 = get_numpy_array_from_Poly(num,den)
    
    num2 = np.poly1d(num2)
    den2 = np.poly1d(den2)
    
    H = sp.lti(num2,den2)
    t,y,sec=sp.lsim(H,x,t)   
    return t,y

Vi = 1
H_l,w,v,ph = get_response_freq_domain(Vi,'l')


fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex = True)
plt.suptitle('Bode plots of low pass transfer function')
axes[0].loglog(w,v,lw=2)
axes[0].grid()
axes[0].set_ylabel('Magnitude in log scale')
axes[0].set_title('Magnitude Plot')

axes[1].semilogx(w,ph,lw=2)
axes[1].grid()
axes[1].set_xlabel('Frequency in log scale')
axes[1].set_ylabel('Phase in radians')
axes[1].set_title('Phase Plot')
plt.show()




Vi = 1/s
Y,w,v,ph = get_response_freq_domain(Vi,'l')
t,y = get_output_time_domain(Y,4e-3,10001)


plt.loglog(w,v,lw=2)
plt.grid()
plt.title('Step Response of the low pass filter')
plt.xlabel('Frequency in log scale')
plt.ylabel('Magnitude in log scale')
plt.show()

 
plt.plot(t,y)
plt.grid()
plt.title('Time domain output to step input')
plt.xlabel('Time')
plt.ylabel('y(t)')
plt.ylim(-1,1)
plt.show()


t = np.linspace(0.0,4e-3,100001)     
x = np.sin(2000*math.pi*t) + np.cos(2*(10**6)*math.pi*t)   
t,y = get_output_with_lsim(H_l,x,t)  

fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey = True)
axes[0].plot(t,y)
axes[0].grid()
axes[0].set_title('Output of lowpass filter to given input')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('y(t)')

axes[1].plot(t,x)
axes[1].grid()
axes[1].set_title('Input to lowpass filter')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('y(t)')
plt.show()



Vi = 1
H_h,w,v,ph = get_response_freq_domain(Vi,'h') 


fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex = True)
plt.suptitle('Bode plots of high pass transfer function')
axes[0].loglog(w,v,lw=2)
axes[0].grid()
axes[0].set_ylabel('Magnitude in log scale')
axes[0].set_title('Magnitude Plot')

axes[1].semilogx(w,ph,lw=2)
axes[1].grid()
axes[1].set_xlabel('Frequency in log scale')
axes[1].set_ylabel('Phase in radians')
axes[1].set_title('Phase Plot')
plt.show()




Vi = 1/s
Y,w,v,ph = get_response_freq_domain(Vi,'h')
t,y = get_output_time_domain(Y,4e-3,10001)


plt.loglog(w,v,lw=2)
plt.grid()
plt.title('Step Response of the high pass filter')
plt.xlabel('Frequency in log scale')
plt.ylabel('Magnitude in log scale')
plt.show()

 
plt.plot(t,y)
plt.grid()
plt.title('Time domain output to step input')
plt.xlabel('Time')
plt.ylabel('y(t)')
plt.ylim(-1,1)
plt.show()



t = np.linspace(0.0,4e-3,100001)     
x = np.sin(2000*math.pi*t) + np.cos(2*(10**6)*math.pi*t)   
t,y = get_output_with_lsim(H_h,x,t)   

fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey = True)
axes[0].plot(t,y)
axes[0].grid()
axes[0].set_title('Output of highpass filter to given input')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('y(t)')

axes[1].plot(t,x)
axes[1].grid()
axes[1].set_title('Input to highpass filter')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('y(t)')
plt.show()




t = np.linspace(0.0,4e-3,100001)     
x = (np.sin(2000*math.pi*t))*np.exp((-10**3)*t)
t,y = get_output_with_lsim(H_h,x,t)   

fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey = True)
axes[0].plot(t,y)
axes[0].grid()
axes[0].set_title('Output of highpass filter to given input')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('y(t)')

axes[1].plot(t,x)
axes[1].grid()
axes[1].set_title('Input to highpass filter')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('y(t)')
plt.show()





t = np.linspace(0.0,3e-5,100001)     
x = (np.sin(2*(10**6)*math.pi*t))*np.exp((-10**5)*t)
t,y = get_output_with_lsim(H_h,x,t)   

fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey = True)
axes[0].plot(t,y)
axes[0].grid()
axes[0].set_title('Output of highpass filter to given input')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('y(t)')

axes[1].plot(t,x)
axes[1].grid()
axes[1].set_title('Input to highpass filter')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('y(t)')
plt.show()

