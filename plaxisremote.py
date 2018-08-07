__author__ = 'Dazhong Li'
__version__ = '1.0'

import os
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt

#plaxis related module

import imp
path = r'C:\Program Files\Plaxis\PLAXIS 2D'
found_module = imp.find_module('plxscripting', [path])
plx = imp.load_module('plxscripting', *found_module)
from plxscripting.easy import *

def find_hoop_force (anchor_result,y_loc,R=15.597):
    '''
    This function calculate the hoop force given the anchor result. 
    Parameter:
        anchor_result: a data frame containing the anchor result in a form of 'xa','ya','xb','yb','F','Fmax'
        y_loc: a string indicating the location to be calculated, either 'ya' or 'yb'
        R: radius of the circular structure, the default value is what we've been using for the project at BCF
    Return:
        A data from contaning the hoop force in a form of y|force
    '''
    hoop_force = pd.DataFrame(dict(y=[],force=[]))
    hoop_force.set_index('y')
    for y in sorted(set(anchor_result[y_loc]),reverse=True):
        point_result = anchor_result[abs(anchor_result[y_loc] -y) <0.1]
        H = 0
        R = 15.597
        N = 0
        for xa, ya, xb, yb, Nmax in zip(point_result.xa,
                                  point_result.ya,
                                  point_result.xb,
                                  point_result.yb,point_result.F):
            H += 0.5*abs(yb-ya)
            theta = np.arctan((yb - ya)/(xb - xa))
            N += Nmax*np.cos(theta)
        hoop_force.loc[y,'force']=N*R/H
        hoop_force.loc[y,'y']=y
    return hoop_force

class Input:

    def __init__(self,port):
        self.s_i = new_server('localhost', port)[0]
        self.g_i = new_server('localhost', port)[1]

class Output:

    def __init__(self,port):
        '''
        This function initialize the plaxis remote with path and port
        Parameter:
            path - path to where plaxis is installed, e.g. 'C:\Program Files\Plaxis\PLAXIS 2D'
            port - port number to be fitt in
        return:
        '''
        self.s_o = new_server('localhost', port)[0]
        self.g_o = new_server('localhost', port)[1]

    def open_file(self,filename):
        self.s_o.open(filename)

    def get_result_epwp(self,Phase):
        '''
        This function return the excess pore water perssure at a certain phase
        Parameter:
        g_o: Plaxis output object
        Phase: Phase of interested,example 'g_o.Phase_6'
        '''
        g_o = self.g_o
        Phase = self.get_phase_by_name(Phase)
        x = np.array(g_o.getresults(Phase, g_o.ResultTypes.Soil.X, "node"))
        y = np.array(g_o.getresults(Phase, g_o.ResultTypes.Soil.Y, "node"))
        epwp =np.array(g_o.getresults(Phase, g_o.ResultTypes.Soil.PExcess, "node"))
        df = pd.DataFrame(dict(x=[],y=[],epwp=[]))
        df.x = x; df.y = y; df.epwp = epwp
        return df

    def get_phase_by_name(self,phase_name):

        #it's strange that the name show in phase ID is different from the internal name by Plaxis, which is a bug, 
        #Instead, we have to surcharge ID rather than Name
        #phase =  [x for x in self.g_o.Phases[:] if x.Name==phase_name]
        phase =  [x for x in self.g_o.Phases[:] if phase_name in x.Identification[:]]
        assert(len(phase) ==1)
        return phase[0]

    def get_anchor_force(self,phase):
        '''
        This function get the anchor force in a certain phase
        '''
        phase = self.get_phase_by_name(phase)
        F = np.array(self.g_o.getresults(phase,
                                         self.g_o.ResultTypes.NodeToNodeAnchor.AnchorForce2D,'node') )
        Fmax = np.array(self.g_o.getresults(phase,
                                         self.g_o.ResultTypes.NodeToNodeAnchor.AnchorForceMax2D,'node') )
        x = np.array(self.g_o.getresults(phase, self.g_o.ResultTypes.NodeToNodeAnchor.X,'node') )
        y = np.array(self.g_o.getresults(phase, self.g_o.ResultTypes.NodeToNodeAnchor.Y,'node') )
        df = pd.DataFrame(dict(x=[],y=[],F=[],Fmax=[]))
        df.x = x; df.y = y; df.F =F; df.Fmax =Fmax
        # put the two end points at the same row 
        df_start = df.iloc[0:-1:2,:]
        df_end = df.iloc[1:-1:2,:]
        ix = 0
        df_anchor = pd.DataFrame(dict(xa=[], ya=[],xb=[],yb=[],F=[], Fmax=[]))
        for xa, ya, xb, yb, fa, fb, fmax in zip (df_start.x,df_start.y, df_end.x, df_end.y, df_start.F,df_end.F, df_start.Fmax):
            #import pdb; pdb.set_trace()
            #we always put the left point as xa
            if xa > xb:
                xa, xb = xb, xa
                ya, yb = yb, ya
            assert(abs(fa-fb)<1e-6) # we want make sure the two points form the same anchor by checking the force. 
            df_anchor.loc[ix] =[xa, ya,xb, yb, fa,fmax]
            ix += 1
        return df_anchor.sort_values(by ='ya',ascending=False)
