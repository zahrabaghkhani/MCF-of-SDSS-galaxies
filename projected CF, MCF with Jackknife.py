import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo

#This code computes projected correlation function(w_p) and projected marked correlation function(m_p) of observational galaxies(ra,dec,z)!
#Errorbars are calculated using Jackknife resampling!
#You need to have (ra,dec,z) of galaxies as the positions and a propertiy e.g. stellar mass as the mark! Also, you should have random catalogue!
#For SDSS, see random_sdss.py
#To use this code, you should first convert ra,dec,z into X,Y,Z! 

h=cosmo.H(0)/100

def sky_to_cart(ra,dec,z):
    ra_rad = ra*np.pi/180
    dec_rad = dec*np.pi/180
    r = cosmo.comoving_distance(z).value*h
    x = r*np.cos(dec_rad)*np.cos(ra_rad)
    y = r*np.cos(dec_rad)*np.sin(ra_rad)
    z = r*np.sin(dec_rad)
    return x, y, z


def distance_calculator(X1,X2): #between a point and an array

        
    pi = np.abs(((X1[:,0]+X2[0])*(X1[:,0]-X2[0]) + (X1[:,1]+X2[1])*(X1[:,1]-X2[1]) + (X1[:,2]+X2[2])*(X1[:,2]-X2[2]))/np.sqrt((X1[:,0]+X2[0])**2 + (X1[:,1]+X2[1])**2 + (X1[:,2]+X2[2])**2))
    rp = np.sqrt((X1[:,0]-X2[0])**2 + (X1[:,1]-X2[1])**2 + (X1[:,2]-X2[2])**2 - pi**2)
    
    return pi,rp

def MCF(data,marks,rbins,pi_max,random_catalog):#data = XYZ positions of galaxies, marks, bins, maximum distance parallel to LOS, random catalogue
    
    
    r_max = np.amax(rbins)
    pibins = np.arange(0,pi_max,1)
    N_random =len(random_catalog)
    s_max = np.sqrt(r_max**2 + pi_max**2)
     
    
    
    distance = []
    distance_random = []
    distance_DR = []
    DD = []
    RR = []
    DR= []
    mark_array=[]
    mark_array_WR = []


    WW = []
    WR = []
    
    #To make the code faster, I considered a box with the size of 2*s_max around each galaxy.
    # I just computed separarions for galaxies in this box
    for i in range(len(data)): 
        
 
        
        condition_1 = (data[:,0]>data[i,0]-s_max)&(data[:,0]<data[i,0]+s_max)
        condition_2 = (data[:,1]>data[i,1]-s_max)&(data[:,1]<data[i,1]+s_max)
        condition_3 = (data[:,2]>data[i,2]-s_max)&(data[:,2]<data[i,2]+s_max)
       
        distance.append(distance_calculator(data[i+1:,:][(condition_1&condition_2&condition_3)[i+1:]],data[i]))
        
        
        mark_array.append(marks[i]*marks[i+1:][(condition_1&condition_2&condition_3)[i+1:]]) #m_i*m_j/m^2
        if (i%100000==0):
            print(i)
    print('DD and WW')
    del condition_1
    del condition_2
    del condition_3

    DD = np.zeros((len(pibins)-1,len(rbins)-1))
    WW = np.zeros((len(pibins)-1,len(rbins)-1))
    distance = np.array(distance)



    pi_DD = [i for sublist in distance[:,0] for i in sublist]
    r_DD = [i for sublist in distance[:,1] for i in sublist]

    del distance 

    mark_array = [i for sublist in mark_array for i in sublist]
    mark_array = np.array(mark_array)

    pi_DD = np.array(pi_DD)
    r_DD = np.array(r_DD)
    
    for i in range(len(pibins)-1):
        for j in range(len(rbins)-1):
            condition = (r_DD<rbins[j+1])&(r_DD>rbins[j])&(pi_DD<pibins[i+1])&(pi_DD>pibins[i])
 
            
            DD[i,j] = np.count_nonzero(condition)
       
        
            WW[i,j] = mark_array[condition].sum()
    
    
    print('DD is done!')

    #To reduce memory used for this code variables are removed when we do not need them anymore!
    del condition
     
    del r_DD
    del pi_DD
    del mark_array


    
    for i in range(len(random_catalog)):
      
       
     
        condition_1 = (random_catalog[:,0]>random_catalog[i,0]-s_max)&(random_catalog[:,0]<random_catalog[i,0]+s_max)
        condition_2 = (random_catalog[:,1]>random_catalog[i,1]-s_max)&(random_catalog[:,1]<random_catalog[i,1]+s_max)
        condition_3 = (random_catalog[:,2]>random_catalog[i,2]-s_max)&(random_catalog[:,2]<random_catalog[i,2]+s_max)
     
        distance_random.append(distance_calculator(random_catalog[i+1:,:][(condition_1&condition_2&condition_3)[i+1:]],random_catalog[i]))
        if (i%100000==0):
            print(i)
    print('RR')

    del condition_1
    del condition_2
    del condition_3
    RR = np.zeros((len(pibins)-1,len(rbins)-1))

    distance_random = np.array(distance_random)

    pi_RR = [i for sublist in distance_random[:,0] for i in sublist]
    r_RR = [i for sublist in distance_random[:,1] for i in sublist]
    del distance_random
    pi_RR = np.array(pi_RR)
    r_RR = np.array(r_RR)
    
    for i in range(len(pibins)-1):
        for j in range(len(rbins)-1):
        
            condition = (r_RR<rbins[j+1])&(r_RR>rbins[j])&(pi_RR<pibins[i+1])&(pi_RR>pibins[i])
            RR[i,j] = np.count_nonzero(condition)
            
    del r_RR
    del pi_RR
    del condition
    print('RR is done!')



   
    for i in range(len(data)): 
       
        condition_1 = (random_catalog[:,0]>data[i,0]-s_max)&(random_catalog[:,0]<data[i,0]+s_max)
        condition_2 = (random_catalog[:,1]>data[i,1]-s_max)&(random_catalog[:,1]<data[i,1]+s_max)
        condition_3 = (random_catalog[:,2]>data[i,2]-s_max)&(random_catalog[:,2]<data[i,2]+s_max)
       
        distance_DR.append(distance_calculator(random_catalog[(condition_1&condition_2&condition_3)],data[i]))
        marks_WR = [marks[i]]*len(random_catalog[(condition_1&condition_2&condition_3)])
        mark_array_WR.append(marks_WR)
        if (i%100000==0):
            print(i/100000)
    print('DR and WR')
    
 
    del condition_1
    del condition_2
    del condition_3
 
    
    
    DR = np.zeros((len(pibins)-1,len(rbins)-1))
    
    WR = np.zeros((len(pibins)-1,len(rbins)-1))

    
   
    distance_DR = np.array(distance_DR)

    pi_DR = [i for sublist in distance_DR[:,0] for i in sublist]
    r_DR = [i for sublist in distance_DR[:,1] for i in sublist]
    del distance_DR
    pi_DR = np.array(pi_DR)
    r_DR = np.array(r_DR)
    
    
    
    mark_array_WR = [i for sublist in mark_array_WR for i in sublist]
    
    
    mark_array_WR=np.array(mark_array_WR)
    
    
    
    

    
 
    
    

    for i in range(len(pibins)-1):
        for j in range(len(rbins)-1):
            condition = (r_DR<rbins[j+1])&(r_DR>rbins[j])&(pi_DR<pibins[i+1])&(pi_DR>pibins[i])
            DR[i,j] = np.count_nonzero(condition)
            WR[i,j] = mark_array_WR[condition].sum()

    del r_DR
    del pi_DR
    del condition 
    del mark_array_WR
    print('DR is done!')
    



   

  
 

   

  
        

    xi = ((np.array(DD)/(len(data)*(len(data)-1)/2)) - (2*np.array(DR)/(len(data)*N_random)) +(np.array(RR)/(N_random*(N_random-1)/2)))/(np.array(RR)/(N_random*(N_random-1)/2)) 
    
 

    

        
  
        
   

    
    

        
        
        
        
  

        
    W = ((np.array(WW)/(len(data)*(len(data)-1)/2)) - (2*np.array(WR)/(len(data)*N_random)) +(np.array(RR)/(N_random*(N_random-1)/2)))/(np.array(RR)/(N_random*(N_random-1)/2)) 

   
    xi_p = 2*xi.sum(axis=0)
    W_p = 2*W.sum(axis=0)
    bincenters = 0.5*(rbins[1:]+rbins[:-1])


    M = (1+(W_p/bincenters))/(1+(xi_p/bincenters))
        
    
    
    return xi_p,W_p,M

def jackknife(data,marks,random,rbins,pi_max,n_jackknife): #data['ra','dec','X','Y','Z'],marks,random catalogue, bins, max value of distance parallel to LOS, # Jackknife samples
    M_subsample = [] 
    C_ij = np.zeros((len(rbins)-1,len(rbins)-1))
    Error = []
    
    ra_max = np.amax(data[:,0])
    ra_min = np.amin(data[:,0])
    dec_max = np.amax(data[:,1])
    dec_min = np.amin(data[:,1])

    l_ra = ra_max - ra_min
    l_dec = dec_max - dec_min
    for i in range(n_jackknife):
        for j in range(n_jackknife):
            
            condition_ra = (data[:,0]<(ra_min + ((i+1)*l_ra/n_jackknife)))&(data[:,0]>(ra_min + (i*l_ra/n_jackknife)))
            condition_dec = (data[:,1]<(dec_min + ((j+1)*l_dec/n_jackknife)))&(data[:,1]>(dec_min + (j*l_dec/n_jackknife)))
            
            sub_sample = data[((condition_ra&condition_dec)==False)]
            mark_sample = marks[((condition_ra&condition_dec)==False)]
                
            condition_ra = (random[:,0]<(ra_min + ((i+1)*l_ra/n_jackknife)))&(random[:,0]>(ra_min + (i*l_ra/n_jackknife)))
            condition_dec = (random[:,1]<(dec_min + ((j+1)*l_dec/n_jackknife)))&(random[:,1]>(dec_min + (j*l_dec/n_jackknife)))
            
            sub_sample_random = random[((condition_ra&condition_dec)==False)]
            M_subsample.append(MCF(sub_sample[:,2:5],mark_sample,rbins,pi_max,sub_sample_random[:,2:5]))
            print('n=',j)
                
    M_subsample = np.array(M_subsample)
    M_tot = MCF(data[:,2:5],marks,rbins,pi_max,random[:,2:5])
    print('M_tot has been done!')
    
    n_sample = n_jackknife**2  
    for i in range(len(rbins)-1):
        for j in range(len(rbins)-1):
            C_ij[i,j] = ((n_sample-1)/n_sample)*np.sum((M_subsample[:,i]-np.mean(M_subsample[:,i]))*(M_subsample[:,j]-np.mean(M_subsample[:,j])))
            if (i == j):
                Error.append(np.sqrt(C_ij[i,j]))
    Error = np.array(Error)
    M_mean = np.mean(M_subsample,axis=0)
    
    return M_mean, Error, C_ij   
        
    
    

