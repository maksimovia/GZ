import mat_properties as prop
import pandas as pd

class HTS:
    def __init__(self,stream11,table,Tout, KPD):
        self.KPD=KPD
        self.table=table
        self.stream11=stream11
        self.Tout=Tout

 
    def calc(self):
        #Параметры газа на входе
        syngas_mix = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac=  list(self.table.loc[self.stream11,'N2':'CO'])
        tsyngas= self.table.loc[self.stream11,'T']
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")

        psyngas= self.table.loc[self.stream11,'P']
        gas_KU = prop.Materials_prop(
        syngas_mix,
        SGfrac,
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
        )
        h_sg1=gas_KU.p_t(psyngas,tsyngas)['h']

        #Параметры газа на выходе
        # SGfrac1=  [0.0, 0, 0.06356, 0.0, 0.536615, 0.0500421, 0.325899, 0.0238773000]
                gas_KU1 = prop.Materials_prop(
        syngas_mix,
        SGfrac1,
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
        )
        Tout = self.Tout
        h_sg2 = gas_KU1.p_t(psyngas,self.Tout)['h']

        #Разность энтальпий
        Hc_syngas = h_sg1 - h_sg2
        KPD=self.KPD
        #Тепловая мощность ТО
        gsyngas=  self.table.loc[self.stream11,'G']

        Qc_syngas = Hc_syngas * gsyngas
        Qc_sg_w = Hc_syngas * gsyngas* self.KPD
        res = {"T":Tout,"Qc_syngas":Qc_syngas,"Qc_sg_w":Qc_sg_w,"psyngas":psyngas,"Hc_syngas":Hc_syngas,"gsyngas":gsyngas,"h_sg1":h_sg1,"h_sg2":h_sg2,"KPD":KPD,"SGfrac":SGfrac,"KPD":KPD,"SGfrac":SGfrac}
        
        return res
    
    
        

