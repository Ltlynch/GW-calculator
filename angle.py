
# coding: utf-8

# In 3c66b paper, equation for total residual (1). Try setting phi and theta_n equal to 0, and orbital inclination angle equal to 90', for edge-on. Try to find maximum amplitude? (Also need to calculate amplitudes A (9) and B (10). For now set eccentricity to 0, eventually add it in cuz life's not always circular.

# In[137]:


import numpy as np
import scipy as sci
import astropy.constants as const


# In[138]:


map = np.transpose(np.genfromtxt("11yr_skymap_v4.txt", skip_header=1))

#raj = input("Raj: ")
#dec = input("Dec: ")

ul_col = map[3]
theta_col = map[1]
phi_col = map[2]


# In[139]:


"""
theta_dict = {}
for i in range(len(phi_col)):
    #print(theta_col[i])
    if not theta_col[i] in theta_dict:
        theta_dict[theta_col[i]] = {}
    theta_dict[theta_col[i]][phi_col[i]] = ul_col[i]
    #print(theta_dict[theta_col[i]][phi_col[i]])
print(theta_dict)
"""
theta_dict = {}
#for i in range(len(map[4])):
for i in range(10):
    theta_dict.setdefault(theta_col[i], {})[phi_col[i]] = ul_col[i]
print(theta_dict)


# In[140]:


def convert_angles(raj, dec):
    # For now
    return raj, dec


# In[141]:


def binary_search(array, value):
    mid = (len(array))/2
    if value == array[mid]:
        return mid
    elif len(array) < 2:
        return array
    elif value > array[mid]:
        return binary_search(array[mid:], value)
    elif value < array[mid]:
        return binary_search(array[:mid], value)



# In[142]:


def find_ul(theta, phi):
    target_theta = theta
    target_phi = phi

    theta_array = np.array(list(theta_dict.keys()))
    #print(theta_array)
    
    theta_arg = np.argmin(abs(theta_array - target_theta))
    best_theta = theta_array[theta_arg]
    print("The closest theta value to {0}, is {1}".format(theta, best_theta))
    #print(theta_dict)
    #print(theta_dict[best_theta])
    
    phi_array = np.array(list(theta_dict[best_theta].keys()))
    phi_arg = np.argmin(abs(phi_array - target_phi))
    best_phi = phi_array[phi_arg]
    print("The closest phi value to {0}, given theta ~ {1}, is {2}".format(target_phi, best_theta, best_phi))
    
    print("Given the best values of theta and phi, the continuous source upper limit in the region given is: {0}".format(theta_dict[best_theta][best_phi]))
    #for theta, phi in theta_dict.items():
        #print(phi)
        #print(np.array(list(phi.keys())))
        #print(np.array(list(phi.values())))
    
    #for theta in theta_dict.items
    #nested = list(theta_dict.values())
    #print(nested)
    
    #phi_array = nested.keys()
    
    #ul_array = phi_array(list(phi_array.values()))
    #dist_theta = abs(ul_array - target_phi)
    #arg = np.argmin(dist_theta)
    #answer(k[arg])
    #return binary_search(theta_row, theta)


# In[143]:


def strain(Distance, Radius, Frequency, Mass):
    return (32*(np.pi**2)*const.G.value)/(Distance*const.c.value**4)* Mass*(Radius**2)*(Frequency**2)


# In[144]:


print(strain(4.5e20, 20, 1/(1.05*3.154e+7), 6e30))


# In[145]:


find_ul(.2, .3)

