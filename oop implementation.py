#########################
#LIBRARY IMPORTS
import streamlit as st
import numpy as np
import pandas as pd
import vegas
import time
from scipy.optimize import curve_fit
from statistics import mean
from numpy import arange
from pandas import read_csv
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import logging

from numeric import num_diff
from constant import constants
from collections import namedtuple

###############################

logging.basicConfig(filename='log_information.log', encoding='utf-8', level=logging.DEBUG)

st.title('Second Virial Coefficient Calculator for A2B2 Molecule')

uploaded_file = st.file_uploader("upload a file")


potential = st.selectbox(
    'What potential energy do you want to use?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential', 'both'))


step = st.selectbox(
    'What is the Temperature step that you want to use?',
    (100, 50, 25, 200, 300))

gas_type = st.selectbox(
    'What gas structure are you studying?',
    ('A2B2', 'AB',))
if potential == 'both':
    uploaded_file_2 = st.file_uploader("upload the rydberg file")

if gas_type == 'A2B2':

    monomero1 = st.selectbox(
        'What is your first monomer?',
        ('H2', 'F2', 'Cl2', 'Br2'))
    monomero2 = st.selectbox(
        'What is your second monomer?',
        ('H2', 'F2', 'Cl2', 'Br2'))
    gas = monomero1 + monomero2


st.write('The calculus will be made using ', potential,
          'using a temperature step of', step)

ILJ_result = namedtuple('ILJresult', ["Reqs", "LC_lst", "harm_esf_lst", "Tstr", "B_clas", "De_lim"])


class ILJ:
    ...

class Data_manage:
    def __init__(self, uploaded_file, st):
        self.uploaded_file = uploaded_file
        self.st = st  # You'll need to provide the streamlit object reference here
        self.T = None
        self.data = None
    def process_data(self):

        data = pd.read_csv(self.uploaded_file, sep="\s+", header=None)
        data.columns = ["alpha", "beta", "mp", "De", "Req"]
        self.st.subheader('DataFrame')
        self.st.write(data)
        self.data = data

        

        # Extracting data


class LC:
    def __init__(self, data):
        self.data = data

    def UH(self, r):
        n = self.data.loc[0, 'beta'] + self.data.loc[0, 'alpha'] * (r / self.data.loc[0, 'Req']) ** 2
        #st.write('uh good')
        return self.data.loc[0, 'De'] * ((self.data.loc[0, 'mp'] / (n - self.data.loc[0, 'mp']) * (self.data.loc[0, 'Req'] / r) ** n) - (n / (n - self.data.loc[0, 'mp']) * (self.data.loc[0, 'Req'] / r) ** self.data.loc[0, 'mp']))
    
    def UX(self, r):
        n = self.data.loc[1, 'beta'] + self.data.loc[1, 'alpha'] * (r / self.data.loc[1, 'Req']) ** 2
        return self.data.loc[1, 'De'] * ((self.data.loc[1, 'mp'] / (n - self.data.loc[1, 'mp']) * (self.data.loc[1, 'Req'] / r) ** n) - (n / (n - self.data.loc[1, 'mp']) * (self.data.loc[1, 'Req'] / r) ** self.data.loc[1, 'mp']))
    
    def UZ(self, r):
        n = self.data.loc[2, 'beta'] + self.data.loc[2, 'alpha'] * (r / self.data.loc[2, 'Req']) ** 2
        return self.data.loc[2, 'De'] * ((self.data.loc[2, 'mp'] / (n - self.data.loc[2, 'mp']) * (self.data.loc[2, 'Req'] / r) ** n) - (n / (n - self.data.loc[2, 'mp']) * (self.data.loc[2, 'Req'] / r) ** self.data.loc[2, 'mp']))
    
    def UTa(self, r):
        n = self.data.loc[3, 'beta'] + self.data.loc[3, 'alpha'] * (r / self.data.loc[3, 'Req']) ** 2
        return self.data.loc[3, 'De'] * ((self.data.loc[3, 'mp'] / (n - self.data.loc[3, 'mp']) * (self.data.loc[3, 'Req'] / r) ** n) - (n / (n - self.data.loc[3, 'mp']) * (self.data.loc[3, 'Req'] / r) ** self.data.loc[3, 'mp']))
    
    def UTb(self, r):
        n = self.data.loc[4, 'beta'] + self.data.loc[4, 'alpha'] * (r / self.data.loc[4, 'Req']) ** 2
        return self.data.loc[4, 'De'] * ((self.data.loc[4, 'mp'] / (n - self.data.loc[4, 'mp']) * (self.data.loc[4, 'Req'] / r) ** n) - (n / (n - self.data.loc[4, 'mp']) * (self.data.loc[4, 'Req'] / r) ** self.data.loc[4, 'mp']))
    
    def UL(self, r):
        n = self.data.loc[5, 'beta'] + self.data.loc[5, 'alpha'] * (r / self.data.loc[5, 'Req']) ** 2
        return self.data.loc[5, 'De'] * ((self.data.loc[5, 'mp'] / (n - self.data.loc[5, 'mp']) * (self.data.loc[5, 'Req'] / r) ** n) - (n / (n - self.data.loc[5, 'mp']) * (self.data.loc[5, 'Req'] / r) ** self.data.loc[5, 'mp']))
class Momentum:
    def __init__(self, lc_instance):
        self.lc_instance = lc_instance
    
    def UM_000(self, r):
        #st.write('um000 momentum good')
        return (2 * self.lc_instance.UH(r) + self.lc_instance.UL(r) + 2 *
                (self.lc_instance.UTa(r) + self.lc_instance.UTb(r) + self.lc_instance.UX(r))) / 9

    def UM_202(self, r):
        return 2 * (self.lc_instance.UH(r) - self.lc_instance.UL(r) + self.lc_instance.UTa(r) -
                    2 * self.lc_instance.UTb(r) + self.lc_instance.UX(r)) / (9 * (5**(1 / 2)))

    def UM_022(self, r):
        return 2 * (self.lc_instance.UH(r) - self.lc_instance.UL(r) - 2 * self.lc_instance.UTa(r) +
                    self.lc_instance.UTb(r) + self.lc_instance.UX(r)) / (9 * (5**(1 / 2)))

    def UM_220(self, r):
        return 2 * (4 * self.lc_instance.UH(r) - self.lc_instance.UL(r) - 5 * (self.lc_instance.UTa(r) +
                    self.lc_instance.UTb(r) + self.lc_instance.UX(r)) + 12 * self.lc_instance.UZ(r)) / (45 * (5**(1 / 2)))

    def UM_222(self, r):
        return ((2 / 7)**(1 / 2)) * (13 * self.lc_instance.UH(r) - self.lc_instance.UL(r) + 7 *
                                    (self.lc_instance.UTa(r) + self.lc_instance.UTb(r) - 2 * self.lc_instance.UX(r)) - 12 * self.lc_instance.UZ(r)) / 45

    def UM_224(self, r):
        return ((2 / 35)**(1 / 2) * 8 * self.lc_instance.UH(r) + self.lc_instance.UL(r) + 2 * self.lc_instance.UZ(r)) / 15
    
    '''
    def call_all(self, r):
        UM_000 = self.momentum_instance.UM_000(r)
        UM_202 = self.momentum_instance.UM_202(r)
        UM_022 = self.momentum_instance.UM_022(r)
        UM_220 = self.momentum_instance.UM_220(r)
        UM_222 = self.momentum_instance.UM_222(r)
        UM_224 = self.momentum_instance.UM_224(r)
        return UM_000, UM_202, UM_022
            UM_000, UM_202 = momentum.momentos()

    '''    

class UM_FINAL:
    def __init__(self, momentum_instance):
        self.momentum_instance = momentum_instance

    def calculate(self, r, th_a, th_b, phi):
        #st.write('calculate um final good')
        UM_000 = self.momentum_instance.UM_000(r)
        UM_202 = self.momentum_instance.UM_202(r)
        UM_022 = self.momentum_instance.UM_022(r)
        UM_220 = self.momentum_instance.UM_220(r)
        UM_222 = self.momentum_instance.UM_222(r)
        UM_224 = self.momentum_instance.UM_224(r)

        UM_FINAL = (
            UM_000 +
            (5 ** (1 / 2)) / 4 * (3 * (np.cos(2 * th_a)) + 1) * UM_202 +
            (5 ** (1 / 2)) / 4 * (3 * (np.cos(2 * th_b)) + 1) * UM_022 +
            (5 ** (1 / 2)) / 16 * UM_220 * (
                (3 * (np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) +
                12 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) +
                3 * (1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi)
            ) -
            (14 ** (1 / 2)) * 5 / 112 * UM_222 * (
                (3 * np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) +
            6 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) -
            3 * (1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi) +
            ((3 * (70) ** (1 / 2)) / 112) * UM_224 * (
                (3 * (np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) -
                8 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) +
                ((1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi)) / 2
            )
        )
        return UM_FINAL



class Integrand:
    def __init__(self, um_final_instance, T):
        self.um_final_instance = um_final_instance
        self.T = temperature_instance

    def integrand_vegas(self, x):
        #st.write('integrand good')
        r, th_a, th_b, phi = x
        # Pass LC as an argument to UM_FINAL
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * constants.rk / self.T
        if F < -1:
            F = -1

        ff = constants.N_A / 4 * np.sin(th_a) * np.sin(th_b) * (r ** 2) * (1 - np.exp(-F))

        return ff
class Mean_req:
    def __init__(self, data):
        self.data = data

    def calculate_mean(self):
        Reqs = [self.data.loc[0, 'Req'], self.data.loc[1, 'Req'], self.data.loc[2, 'Req'], 
                self.data.loc[3, 'Req'], self.data.loc[0, 'Req'], self.data.loc[0, 'Req']]

        return mean(Reqs)

class Calculate_virial:
    def __init__(self, integrand_instance, mean_req_instance, temperature_instance):
        self.integrand_instance = integrand_instance
        self.mean_req_instance = mean_req_instance
        self.temperature_instance = temperature_instance
    
    def compute_results(self, step):

        #st.write('outisde look good')
        # Integrator and definition of integration limits
        integ = vegas.Integrator(
        [[0.2 * self.mean_req_instance.calculate_mean(), 2 * self.mean_req_instance.calculate_mean()],
        [0, np.pi], [0, np.pi], [0, 2 * np.pi]])
        B_clas = []
        Tstr = []
        B_main = []

        for temp in range(50, 1000, step):
            logging.info(f"Starting integration for temperature {temp}")
            self.temperature_instance.set_temperature(temp)
            result = integ(self.integrand_instance.integrand_vegas, nitn=10, neval=10000)
            B_clas.append(result.mean)
            Tstr.append(temp)
            B_main.append(result.mean)
            logging.info(f"Integration for temperature {temp} completed")

    
        return B_clas, Tstr, B_main
            

if st.button('Calculate'):
    if uploaded_file is not None:
        
        if potential == 'Improved Leonard-Jonnes Potential':
            imported = Data_manage(uploaded_file, st)
            imported.process_data()
            if imported.data is not None:  # Check if data has been processed
                lc_instance = LC(imported.data)

                momentum_instance = Momentum(lc_instance)

                um_final_instance = UM_FINAL(momentum_instance)

                temperature_instance = TemperatureSetter()

                integrand_instance = Integrand(um_final_instance, temperature_instance)

                x = [1, 2*np.pi, 2*np.pi, np.pi]
                testing = integrand_instance.integrand_vegas(x)

                st.write(testing)

                mean_req_instance = Mean_req(imported.data)
                st.write(mean_req_instance.calculate_mean())

                virial_calculation_instance = Calculate_virial(integrand_instance, mean_req_instance,temperature_instance)
                B_values, T_values, B_main_values = virial_calculation_instance.compute_results(step)

                st.write("B values:", B_values)
                st.write("T values:", T_values)
                st.write("Main B values:", B_main_values)
        elif potential == 'Rydberg Potential':
            pass
        elif potential == "both":

            pass


       
    else:
        st.info('☝️ Upload a .dat file')
