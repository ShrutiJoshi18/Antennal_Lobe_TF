import numpy as np
import matplotlib.pyplot as plt
import json

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

n_n = int(input("Number of Neurons (Default = 120) : ") or "120")
trials = int(input("Number of Trials (Default = 100) : ") or "100")
#time = float(input("Time in ms (Default = 1000) : ") or "1000")
#eps = float(input("Time Resolution in ms (Default = 0.01) : ") or "0.01")
eps = 0.01
pulse_width = float(input("Pulse width (Default = 100) : ") or "100")
start_pulse = float(input("Inter-pulse interval (Default = 50) : ") or "50")

with open('parameters.json', 'r') as read_file:
    params = json.load(read_file)

init_start = 1000    
time_stop = trials*(start_pulse+pulse_width) + start_pulse + init_start 
time = np.arange(start = 0, stop = time_stop, step = eps)
print(len(time)*eps)

p_n = int(3*n_n/4)
l_n = int(n_n/4)



#p_n = 90
#l_n = 60

#n_n = p_n + l_n



I_pn = np.zeros(int(p_n))
I_ln = np.zeros(int(l_n))

peak_pn = 3.5
peak_ln = 3
p_input = 1/3
l_input = 1
current_input = np.zeros((n_n,int(len(time))))

#current_index = np.zeros((int(time_stop/eps),n_n))
#current_index = np.arange(int(start_pulse/eps) , int(start_pulse/eps)+int(pulse_width/eps))

np.random.seed(0)
# I_pn[np.random.choice(int(p_n), size = int(p_input*p_n), replace = False)] = peak_pn*gaussian(np.arange(int(p_input*p_n)), p_input*p_n/2, p_input*p_n*5)

I_pn[np.arange(int(p_input*p_n))] = peak_pn*gaussian(np.arange(int(p_input*p_n)), p_input*p_n/2, p_input*p_n*5)

#I_pn[np.arange(int(p_input*p_n))] = peak_pn

I_ln[np.random.choice(int(l_n), size = int(l_input*l_n), replace = False)] = peak_ln*gaussian(np.arange(int(l_input*l_n)), l_input*l_n/2, 5*l_n*l_input)
# I_ln[np.arange(int(l_input*l_n))] = peak_ln*gaussian(np.arange(int(l_input*l_n)), l_input*l_n/2, 5*l_n*l_input)

#I_ln[np.arange(int(l_n))] = peak_ln

#Generate trials with given pulse widths  

for y in range(n_n):
    if y < p_n:
        
        current_index = np.arange(int(start_pulse/eps), int(start_pulse/eps)+int(pulse_width/eps)) #+ int(np.random.randint(-5,5)/eps)
        
        for trial in range(1,trials):
            current_index = np.concatenate((current_index, np.arange(int(start_pulse/eps) + trial*int((start_pulse+pulse_width)/eps), int(start_pulse/eps)+int(pulse_width/eps) + trial*int((start_pulse+pulse_width)/eps))))
        
        current_index = current_index + int(init_start/eps)
        current_input[y, current_index] = I_pn[y-1]
        
    
    else:
        
        current_index = np.arange(int(start_pulse/eps) , int(start_pulse/eps)+int(pulse_width/eps))#+ int(np.random.randint(-10,10)/eps)
        
        for trial in range(1,trials):
            current_index = np.concatenate((current_index, np.arange(int(start_pulse/eps) + trial*int((start_pulse+pulse_width)/eps)  , int(start_pulse/eps)+int(pulse_width/eps) + trial*int((start_pulse+pulse_width)/eps))))
            
        current_index = current_index + int(init_start/eps)
        current_input[y, current_index] = I_ln[y-int(p_n)-1]

        
#np.random.seed(5)        
current_input[:p_n,int(init_start/eps):] = current_input[:p_n,int(init_start/eps):] + 2.5#np.transpose([1*np.random.random(p_n) + 3]*len(time))    # PN DC
current_input[p_n:,int(init_start/eps):] = current_input[p_n:,int(init_start/eps):] + 2.5#+ np.transpose([0.5*np.random.random(l_n) + 1]*len(time))  # LN DC

current_input[:p_n,:int(init_start/eps)] = current_input[:p_n,:int(init_start/eps)] + 3
current_input[p_n:,:int(init_start/eps)] = current_input[p_n:,:int(init_start/eps)] + 3

# Additive Gaussian Noise
np.random.seed(35)
for y in range(n_n):
    if y < p_n:
        if I_pn[y-1] == 0:
            current_input[y,:] = current_input[y,:] + (np.random.normal(size = current_input[y,:].shape, scale = current_input[y,0]/50))
        else:
            current_input[y,:] = current_input[y,:] + (np.random.normal(size = current_input[y,:].shape, scale = I_pn[y-1]/50)) # Additive Gaussian Noise
    else: 
        if I_ln[y-p_n-1] == 0:
            current_input[y,:] = current_input[y,:] + (np.random.normal(size = current_input[y,:].shape, scale = current_input[y,0]/50))
        else:
            current_input[y,:] = current_input[y,:] + (np.random.normal(size = current_input[y,:].shape, scale = I_ln[y-p_n-1]/50)) 

#plt.plot(current_input[50,::100])
np.save("current",current_input)


    
total_time_dict = {'total_time': time_stop}
trials_dict = {'trials' : trials}
pulse_width_dict = {'pulse_width' : pulse_width}
start_pulse_dict = {'start_pulse' : start_pulse}

params.update(trials_dict)
params.update(pulse_width_dict)
params.update(start_pulse_dict)
params.update(total_time_dict)
   
print(params)

with open('parameters_used.json', 'w') as write_file:
    params = json.dump(params, write_file)