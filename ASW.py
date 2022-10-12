
import mat_properties as prop
import numpy as n
import pandas as pd



class Accum():
    def __init__(self, **kwargs):
        
        # инициализация на случай если не введут set_construct
#         self._V = 1

        if 'water_streams_eto_tyt' in kwargs.keys():
            self.water_streams = kwargs['water_streams_eto_tyt']
        if 'water' in kwargs.keys():
            self._water = kwargs['water']
        self._D = 1
        self._F = 1
        self._H = 1
        self._P_accum = 1e-1
        self._T_accum = 95
        
        self._D = 1
        self._kolichestvo = 1
        self._V = n.pi*self._D**3/4
        self._F = 1.5*n.pi*self._D**2

        self._khi = 1
        self._lambda_min_vata = 0.045
        self.delta_min_vata = 0.01
        self._T_nar_vozd = 15
#         self._water = water    
#         self.water_streams=water_streams

        if 'stream11' in kwargs.keys():
            self._stream11 = kwargs['stream11']
        if 'stream12' in kwargs.keys():
            self._stream12 = kwargs['stream12']
        if 'stream_obratnoi_setevoi_vody' in kwargs.keys():
            self._stream_obratnoi_setevoi_vody = kwargs['stream_obratnoi_setevoi_vody']
        if 'stream_pryamoi_setevoi_vody' in kwargs.keys():
            self._stream_pryamoi_setevoi_vody = kwargs['stream_pryamoi_setevoi_vody']
        if 'T_nar_vozd' in kwargs.keys():
            self._T_nar_vozd = kwargs['T_nar_vozd']
        
        # параметры обратной сетевой воды
 
        self._T_obr_set_voda = self.water_streams.at[self._stream_obratnoi_setevoi_vody,'T']
        self._P_obr_set_voda = self.water_streams.at[self._stream_obratnoi_setevoi_vody,'P']
        #         self._h_obr_set_voda = self.water_streams.at["SWIN-OD",'H'] #в экселе пусто
        self._h_obr_set_voda = self._water.p_t(self._P_obr_set_voda,self._T_obr_set_voda)['h']
        self._G_obr_set_voda = self.water_streams.at[self._stream_obratnoi_setevoi_vody,'G'] # пока не использовал
        # параметры аккумулирующей воды 
        self.water_streams.at[self._stream11,'T'] = self.water_streams.at[self._stream_pryamoi_setevoi_vody,'T']
        
        self.water_streams.at[self._stream11,'P'] = self.water_streams.at[self._stream_pryamoi_setevoi_vody,'P']
        
        #Конструкция
        

        
    def set_construct(self,**kwargs):

        if 'D' in kwargs.keys():
            self._D = kwargs['D']
        if 'd' in kwargs.keys():
            self._D = kwargs['d']
        if 'Diametr' in kwargs.keys():
            self._D = kwargs['Diametr']
        if 'kolichestvo' in kwargs.keys():
            self._kolichestvo = kwargs['kolichestvo']
        if 'Visota' in kwargs.keys():
            self._H = kwargs['Visota']    
        if 'lambda_min_vata' in kwargs.keys():
            self._lambda_min_vata = kwargs['lambda_min_vata']
        if 'delta_min_vata' in kwargs.keys():
            self._delta_min_vata = kwargs['delta_min_vata'] 
            
        self._V = n.pi*self._D**2/4 * self._H
        self._F = n.pi*self._D*self._H + n.pi*self._D**2/2
         
        pass

    def zaryadka(self, tau):
        
        self._T_accum = self.water_streams.at[self._stream_pryamoi_setevoi_vody,'T'] # тут уточнить
        self._h_accum = self.water_streams.at[self._stream_pryamoi_setevoi_vody,'H']# тут уточнить
        self._P_accum = self.water_streams.at[self._stream_pryamoi_setevoi_vody,'P']# тут уточнить

        self.water_streams.at[self._stream11,'T'] = self._T_accum
        self.water_streams.at[self._stream11,'H'] = self._h_accum 
        self.water_streams.at[self._stream11,'P'] = self._P_accum
        
#         self._Q = self._kolichestvo*self._V * (self._h_accum - self._h_obr_set_voda)*self._water.p_t(self._P_accum, self._T_accum)['rho']#kJ

        self._G = self._kolichestvo*self._V *self._water.p_t(self._P_accum, self._T_accum)['rho']/(tau*3600)
        self.water_streams.at[self._stream11,'G'] = self._G
        
        self._Mass = self._kolichestvo*self._V *self._water.p_t(self._P_accum, self._T_accum)['rho']
        self._Q = self._Mass* (self._h_accum - self._h_obr_set_voda)#kJ

        print(self._Q,'OLD Q')
        print(self._h_accum,'OLD _h_accum')
        return {'T_accum': self._T_accum,'Q': self._Q}
     
    def razryadka(self,tau):
        
        self.water_streams.at[self._stream12,'T'] = self._T_accum
        self.water_streams.at[self._stream12,'H'] = self._h_accum 
        self.water_streams.at[self._stream12,'P'] = self._P_accum

  
        self._G = self._Mass/(tau*3600)    
        self.water_streams.at[self._stream12,'G'] = self._G                     
        return {'T_accum': self._T_accum,'h_accum': self._h_accum,'P_accum': self._P_accum,'G': self._G,} 
            
    
    def jdat(self,tau):
#         print("+")
        self._poteri = self._khi*(self._lambda_min_vata/self.delta_min_vata)*self._F*self._kolichestvo *(self._T_accum-self._T_nar_vozd)*tau*3600/1000 #kJ
#         print(self._F,'m2')
        print(self._poteri,"poteri")
        self._Q = self._Q - self._poteri
        print(self._Q,"new Q")
        self._h_accum = self._Q/self._Mass + self._water.p_t(self._P_obr_set_voda, self._T_obr_set_voda)['h']
        self._T_accum = self._water.p_h(self._P_accum, self._h_accum)['T']
        
        
        self.water_streams.at[self._stream12,'T'] = self._T_accum
        self.water_streams.at[self._stream12,'H'] = self._h_accum 
        self.water_streams.at[self._stream12,'P'] = self._P_accum
        
        
        return {'T_accum': self._T_accum, 'poteri': self._poteri, 'Q': self._Q}
    
    
    def chetamostalos_po_Q(self):
        print(self._Q, 'kWt - Мощность, которую можно передать ')
        pass
    
    def chetamostalos_po_T(self):
        print(self._T_accum, 'C - температура воды в баке')
        pass