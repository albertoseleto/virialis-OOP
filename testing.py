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
    def __init__(self, uploaded_file, st):
        self.uploaded_file = uploaded_file
        self.st = st  # You'll need to provide the streamlit object reference here
        self.T = None
    def process_data(self):

        data = pd.read_csv(self.uploaded_file, sep="\s+", header=None)
        data.columns = ["alpha", "beta", "mp", "De", "Req"]
        self.st.subheader('DataFrame')
        self.st.write(data)

        # Extracting data
        self.extract_data(data)

    def extract_data(self, data):
        self.h_alpha = data.loc[0, 'alpha']
        self.h_beta = data.loc[0, 'beta']
        self.h_mp = data.loc[0, 'mp']
        self.h_De = data.loc[0, 'De']
        self.h_Req = data.loc[0, 'Req']

        self.x_alpha = data.loc[1, 'alpha']
        self.x_beta = data.loc[1, 'beta']
        self.x_mp = data.loc[1, 'mp']
        self.x_De = data.loc[1, 'De'] 
        self.x_Req = data.loc[1, 'Req']
            

        self.z_alpha = data.loc[2, 'alpha']
        self.z_beta = data.loc[2, 'beta']
        self.z_mp = data.loc[2, 'mp']
        self.z_De = data.loc[2, 'De'] 
        self.z_Req = data.loc[2, 'Req']
            

        self.ta_alpha = data.loc[3, 'alpha']
        self.ta_beta = data.loc[3, 'beta']
        self.ta_mp = data.loc[3, 'mp']
        self.ta_De = data.loc[3, 'De'] 
        self.ta_Req = data.loc[3, 'Req']
        

        self.tb_alpha = data.loc[4, 'alpha']
        self.tb_beta = data.loc[4, 'beta']
        self.tb_mp = data.loc[4, 'mp']
        self.tb_De = data.loc[4, 'De'] 
        self.tb_Req = data.loc[4, 'Req']
            

        self.l_alpha = data.loc[5, 'alpha']
        self.l_beta = data.loc[5, 'beta']
        self.l_mp = data.loc[5, 'mp']
        self.l_De = data.loc[5, 'De'] 
        self.l_Req = data.loc[5, 'Req']

        st.write('all good in extract data')

    class LC:
        def UH(self, r):
            n = self.h_beta + self.h_alpha * (r / self.h_Req) ** 2
            return  self.h_De * ((self.h_mp/(n - self.h_mp) * (self.h_Req/r) ** n) - (n/(n - self.h_mp) * (self.h_Req/r) ** self.h_mp))
        
        def UX(self, r):
            n = self.x_beta + self.x_alpha * (r / self.x_Req) ** 2
            return self.x_De * ((self.x_mp/(n - self.x_mp) * (self.x_Req/r) ** n) - (n/(n - self.x_mp) * (self.x_Req/r) ** self.x_mp))
        
        def UZ(self, r):
            n = self.z_beta + self.z_alpha * (r / self.z_Req) ** 2
            return self.z_De * ((self.z_mp/(n - self.z_mp) * (self.z_Req/r) ** n) - (n/(n - self.z_mp) * (self.z_Req/r) ** self.z_mp))
        
        def UTa(self, r):
            n = self.ta_beta + self.ta_alpha * (r / self.ta_Req) ** 2
            return self.ta_De * ((self.ta_mp/(n - self.ta_mp) * (self.ta_Req/r) ** n) - (n/(n - self.ta_mp) * (self.ta_Req/r) ** self.ta_mp))
        
        def UTb(self, r):
            n = self.tb_beta + self.tb_alpha * (r / self.tb_Req) ** 2
            return self.tb_De * ((self.tb_mp/(n - self.tb_mp) * (self.tb_Req/r) ** n) - (n/(n - self.tb_mp) * (self.tb_Req/r) ** self.tb_mp))
        
        def UL(self, r):
            n = self.l_beta + self.l_alpha * (r / self.l_Req) ** 2
            return self.l_De * ((self.l_mp/(n - self.l_mp) * (self.l_Req/r) ** n) - (n/(n - self.l_mp) * (self.l_Req/r) ** self.l_mp))
        st.write('all good in LC')
        '''
        def USa(r):
            n = Sa.beta + Sa.alpha * (r / Sa.Req) ** 2
            return Sa.De * ((Sa.mp/(n - Sa.mp) * (Sa.Req/r) ** n) - (n/(n - Sa.mp) * (Sa.Req/r) ** Sa.mp))

        def USb(r):
            n = Sb.beta + Sb.alpha * (r / Sb.Req) ** 2
            return Sb.De * ((Sb.mp/(n - Sb.mp) * (Sb.Req/r) ** n) - (n/(n - Sb.mp) * (Sb.Req/r) ** Sb.mp))
        '''
    class momentum:
        def UM_000(self, r, LC):
            st.write('working fine in UM_000')
            return (2 * LC.UH(r) + LC.UL(r) + 2 *
                    (LC.UTa(r) + LC.UTb(r) + LC.UX(r))) / 9

        def UM_202(self, r, LC):
            return 2 * (LC.UH(r) - LC.UL(r) + LC.UTa(r) -
                        2 * LC.UTb(r) + LC.UX(r)) / (9 * (5**(1 / 2)))

        def UM_022(self, r, LC):
            return 2 * (LC.UH(r) - LC.UL(r) - 2 * LC.UTa(r) +
                        LC.UTb(r) + LC.UX(r)) / (9 * (5**(1 / 2)))

        def UM_220(self, r, LC):
            return 2 * (4 * LC.UH(r) - LC.UL(r) - 5 * (LC.UTa(r) +
                        LC.UTb(r) + LC.UX(r)) + 12 * LC.UZ(r)) / (45 * (5**(1 / 2)))

        def UM_222(self, r, LC):
            return ((2 / 7)**(1 / 2)) * (13 * LC.UH(r) - LC.UL(r) + 7 *
                                        (LC.UTa(r) + LC.UTb(r) - 2 * LC.UX(r)) - 12 * LC.UZ(r)) / 45

        def UM_224(self, r, LC):
            return ((2 / 35)**(1 / 2) * 8 * LC.UH(r) + LC.UL(r) + 2 * LC.UZ(r)) / 15

    def UM_FINAL(self, r, th_a, th_b, phi, LC):
        UM_000 = self.UM_000(r, LC)
        UM_202 = self.UM_202(r, LC)
        UM_022 = self.UM_022(r, LC)
        UM_220 = self.UM_220(r, LC)
        UM_222 = self.UM_222(r, LC)
        UM_224 = self.UM_224(r, LC)

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
        st.write('working fine in UM_FINAL')

        return UM_FINAL


    last_F_container = [-1e100]

    def set_temperature(self, T):
        self.T = T
    
    def integrand_vegas(self, x, LC):
        last_F = ILJ.last_F_container[0]
        r = x[0]
        th_a = x[1]
        th_b = x[2]
        phi = x[3]
        # Pass LC as an argument to UM_FINAL
        F = self.UM_FINAL(r, th_a, th_b, phi, LC) * constants.rk / self.T
        if F < -1:
            F = -1

        ff = constants.N_A / 4 * np.sin(th_a) * np.sin(th_b) * (r ** 2) * (1 - np.exp(-F))

        if ff < -1e4:
            if F >= last_F:
                logging.debug(f'ff value: {ff}, F value: {F}, r: {r}, th_a: {th_a}, th_b: {th_b}, phi: {phi}')
                ILJ.last_F_container[0] = F
            return ff
        else:
            return ff

        

    def compute_results(self, step, LC):
        Reqs = [self.h_Req, self.x_Req, self.z_Req, self.ta_Req, self.tb_Req, self.l_Req]
        media_Reqs = mean(Reqs)
        self.st.write(media_Reqs)

        # Integrator and definition of integration limits
        integ = vegas.Integrator(
            [[0.2 * media_Reqs, 2 * media_Reqs], [0, np.pi], [0, np.pi], [0, 2 * np.pi]])
        B_clas = []
        Tstr = []
        B_main = []

        for T in range(50, 1000, step):
            self.set_temperature(T)  # Set the temperature for the current iteration
            # Pass LC as an argument to the integrand_vegas method
            result = integ(lambda x: self.integrand_vegas(x, LC), nitn=10, neval=10000)
            self.st.write('result of classic virial = ', result, 'for T = ', T)
            B_clas.append(result.mean)
            Tstr.append(T)

        self.st.write('result of final virial =', result.mean, 'for T = ', T)
        B_main.append(result.mean)

    

if st.button('Calculate'):
    if uploaded_file is not None:
        
        if potential == 'Improved Leonard-Jonnes Potential':
            ilj_instance = ILJ(uploaded_file, st)  # Pass the uploaded file to the class constructor
            ilj_instance.process_data()
        
            ilj_instance.compute_results(step, ILJ.LC)
        elif potential == 'Rydberg Potential':
            pass
        elif potential == "both":

            pass


       
    else:
        st.info('☝️ Upload a .dat file')
