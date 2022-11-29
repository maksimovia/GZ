from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
import CoolProp.CoolProp as CP
os.environ["RPPREFIX"] = r"C:/Program Files (x86)/REFPROP"

def init_REFPROP(path_to_refplot):
    try:
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        #RP = REFPROPFunctionLibrary('/'.join([path_to_refplot, 'REFPROP.dll']))
    except ValueError:
        print('Не удалось загрузить библиотеку REFPROP, проверьте путь к папке')
        return None
    else:
        RP.SETPATHdll(path_to_refplot)
        return RP

def REFPROP_h_s(h, s, gas,fraction, RP):
    RP.PREOSdll(0)
    prop = RP.REFPROPdll(gas, 'HS', 'P;T;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, h, s, fraction)
    res = dict()
    res['p'] = prop.Output[0]
    res['T'] = prop.Output[1]-273.15
    res['rho'] = prop.Output[2]
    res['cv'] = prop.Output[3]
    res['cp'] = prop.Output[4]
    fraction_local=list(fraction)
    if fraction_local[3]>0.05:
        fraction_local[2]=fraction_local[2]+fraction_local[3]-0.05
        fraction_local[3]=0.05
    prop1 = RP.REFPROPdll(gas, 'HS', 'P;T;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, h, s, fraction_local)
    res['nu'] = prop1.Output[5] / 100
    res['Prandtl'] = prop1.Output[6] 
    res['L'] = prop1.Output[7]
    res['Q'] = prop1.Output[9]
    return res

def REFPROP_p_t(p, t, gas,fraction, RP):
    if fraction[0]==1 or fraction[0]==0: #or gas.split("*")[0]!='Nitrogen':   
        RP.PREOSdll(0)
    else:
        RP.PREOSdll(2)
    prop = RP.REFPROPdll(gas, 'PT', 'H;S;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, p, t, fraction)
    res = dict()
    res['h'] = prop.Output[0]/1000
    res['s'] = prop.Output[1]/1000
    res['rho'] = prop.Output[2]
    res['cv'] = prop.Output[3]
    res['cp'] = prop.Output[4]
    fraction_local=list(fraction)
    if fraction_local[3]>0.05 and fraction_local[0]>0.5:
        fraction_local[2]=fraction_local[2]+fraction_local[3]-0.05
        fraction_local[3]=0.05
    prop1 = RP.REFPROPdll(gas, 'PT', 'H;S;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, p, t, fraction_local)
    res['nu'] = prop1.Output[5] / 100
    res['Prandtl'] = prop1.Output[6] 
    res['L'] = prop1.Output[7]
    res['Q'] = prop1.Output[9]
    if res['L']<0 and fraction_local[0]>0.5 and gas.split("*")[0]=='Nitrogen':
        print("Ошибка в расчете по p и t")
        # print(f"p: {p}, t: {t}, gas: {gas},fraction: {fraction}, RP: {RP}")
    return res

def REFPROP_p_h(p, h, gas,fraction, RP):
    if h<200000 and gas == "water":
        res = dict()
        res['T'] = CP.PropsSI('T','P', p,'H',h,gas)-273.15
        res['s'] = CP.PropsSI('S','P', p,'H',h,gas)/1000
        res['rho'] = CP.PropsSI('D','P', p,'H',h,gas)
        res['cv'] = CP.PropsSI('C','P', p,'H',h,gas)
        res['cp'] = CP.PropsSI('C','P', p,'H',h,gas)
        res['Q'] = CP.PropsSI('Q','P', p,'H',h,gas)
    else:    
        if fraction[0]==1 or gas.split("*")[0]!='Nitrogen':   
            RP.PREOSdll(0)
        else:
            RP.PREOSdll(2)
        prop = RP.REFPROPdll(gas, 'PH', 'T;S;D;CV;CP;KV;Q;Prandtl;TCX;VIS;QMass', 21, 0, 0, p, h, fraction)
        res = dict()
        res['T'] = prop.Output[0]-273.15
        res['s'] = prop.Output[1]/1000
        res['rho'] = prop.Output[2]
        res['cv'] = prop.Output[3]
        res['cp'] = prop.Output[4]
        k = prop.Output[5]
        res['q']  = prop.Output[6]
        fraction_local=list(fraction)
        if fraction_local[3]>0.05:
            fraction_local[2]=fraction_local[2]+fraction_local[3]-0.05
            fraction_local[3]=0.05
        prop1 = RP.REFPROPdll(gas, 'PH', 'T;S;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, p, h, fraction_local)
        res['nu'] = prop1.Output[5] / 100
        res['Prandtl'] = prop1.Output[6] 
        res['L'] = prop1.Output[7] 
        res['Q'] = prop1.Output[9]
    # if res['T']<-999 or res['s']<-1000:
    #     print("Ошибка в расчете по p и h, переход на расчет по Coolprop")
    #     # print(f"p: {p}, q: {q}, gas: {gas},fraction: {fraction}, RP: {RP}")
    #     res['T'] = CP.PropsSI('T','P', p,'H',h,gas)-273.15
    #     res['rho'] = CP.PropsSI('D','P', p,'H',h,gas)
    #     res['s'] = CP.PropsSI('H','P', p,'H',s,gas)/1000
    #     print(f"t: {res['T']}, p: {p},h: {h}, gas: {gas},fraction: {fraction}, RP: {RP}")       
    return res

def REFPROP_p_s(p, s, gas,fraction, RP):
    RP.PREOSdll(0)
    prop = RP.REFPROPdll(gas, 'PS', 'T;H;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, p, s, fraction)
    res = dict()
    res['T'] = prop.Output[0]-273.15
    res['h'] = prop.Output[1]/1000
    res['rho'] = prop.Output[2]
    res['cv'] = prop.Output[3]
    res['cp'] = prop.Output[4]
    fraction_local=list(fraction)
    if fraction_local[3]>0.05:
        fraction_local[2]=fraction_local[2]+fraction_local[3]-0.05
        fraction_local[3]=0.05
    prop1 = RP.REFPROPdll(gas, 'PS', 'T;H;D;CV;CP;KV;Prandtl;TCX;VIS;Qmass', 21, 0, 0, p, s, fraction_local)
    res['nu'] = prop1.Output[5] / 100
    res['Prandtl'] = prop1.Output[6] 
    res['L'] = prop1.Output[7] 
    res['Q'] = prop1.Output[9]
    if res['h']<0:
        print("Ошибка в расчете по p и s, переход на расчет по Coolprop")
        # print(f"p: {p}, q: {q}, gas: {gas},fraction: {fraction}, RP: {RP}")
        res['T'] = CP.PropsSI('T','P', p,'S',s,gas)-273.15
        res['rho'] = CP.PropsSI('D','P', p,'S',s,gas)
        res['h'] = CP.PropsSI('H','P', p,'S',s,gas)/1000
        print(f"h: {res['h']}, s: {s}, gas: {gas},fraction: {fraction}, RP: {RP}")    
    return res

def REFPROP_p_q(p, q, gas,fraction, RP):
    RP.PREOSdll(0)
    prop = RP.REFPROPdll(gas, 'PQ', 'T;D;H;S', 21, 0, 0, p, q,fraction)
    res = dict()
    res['T'] = prop.Output[0]-273.15
    res['rho'] = prop.Output[1]
    res['h'] = prop.Output[2]/1000
    res['s'] = prop.Output[3]/1000
    if res['h']<0:
        print("Ошибка в расчете по p и q, переход на расчет по Coolprop")
        # print(f"p: {p}, q: {q}, gas: {gas},fraction: {fraction}, RP: {RP}")
        res['T'] = CP.PropsSI('T','P', p,'Q',q,gas)-273.15
        res['rho'] = CP.PropsSI('D','P', p,'Q',q,gas)
        res['h'] = CP.PropsSI('H','P', p,'Q',q,gas)/1000
        res['s'] = CP.PropsSI('S','P', p,'Q',q,gas)/1000
        print(f"h: {res['h']}, q: {q}, gas: {gas},fraction: {fraction}, RP: {RP}")
    return res

def REFPROP_s_q(s, q, gas,fraction, RP):
    RP.PREOSdll(0)
    prop = RP.REFPROPdll(gas, 'SQ', 'T;D;H;P', 21, 0, 0, s, q,fraction)
    res = dict()
    res['T'] = prop.Output[0]-273.15
    res['rho'] = prop.Output[1]
    res['h'] = prop.Output[2]/1000
    res['P'] = prop.Output[3]/1000000
    return res

def REFPROP_p_rho(p, rho, gas,fraction, RP):
    RP.PREOSdll(0)
    prop = RP.REFPROPdll(gas, 'PD', 'T;H;S', 21, 0, 0, p, rho,fraction)
    res = dict()
    res['T'] = prop.Output[0]-273.15
    res['h'] = prop.Output[1]/1000
    res['s'] = prop.Output[2]/1000
    return res

def REFPROP_t_q(t, q, gas,fraction, RP):
    RP.PREOSdll(0)
    prop = RP.REFPROPdll(gas, 'TQ', 'P;D;H;S', 21, 0, 0, t, q,fraction)
    res = dict()
    res['P'] = prop.Output[0]/1e6
    res['rho'] = prop.Output[1]
    res['h'] = prop.Output[2]/1000
    res['s'] = prop.Output[3]/1000

    return res

class Materials_prop:
    def __init__(self, mat_name, fraction, hs_func, pt_func, ph_func, ps_func,pq_func, tq_func, prho_func,sq_func, *args, **kwargs):
        self.__mat_name = mat_name
        self.__fraction = fraction
        self.__hs_func = hs_func
        self.__pt_func = pt_func
        self.__ph_func = ph_func
        self.__ps_func = ps_func
        self.__pq_func = pq_func
        self.__tq_func = tq_func
        self.__prho_func = prho_func
        self.__sq_func = sq_func
        self.__params = dict()
        for key, val in kwargs.items():
            self.__params[key] = val

    def h_s(self, h, s):
        return self.__hs_func(h=h*1000, s=s*1000, gas=self.__mat_name,fraction=self.__fraction,  **self.__params)

    def p_t(self, p, t):
        return self.__pt_func(p=p*1e6, t=t+273.15, gas=self.__mat_name,fraction=self.__fraction, **self.__params)

    def p_h(self, p, h):
        return self.__ph_func(p=p*1e6, h=h*1000, gas=self.__mat_name,fraction=self.__fraction, **self.__params)

    def p_s(self, p, s):
        return self.__ps_func(p=p*1e6, s=s*1000, gas=self.__mat_name,fraction=self.__fraction, **self.__params)
    
    def p_q(self, p, q):
        return self.__pq_func(p=p*1e6, q=q, gas=self.__mat_name,fraction=self.__fraction, **self.__params)
    
    def t_q(self, t, q):
        return self.__tq_func(t=t+273.15, q=q, gas=self.__mat_name,fraction=self.__fraction, **self.__params)
    
    def p_rho(self, p, rho):
        return self.__prho_func(p=p*1e6, rho=rho, gas=self.__mat_name,fraction=self.__fraction, **self.__params)
    
    def s_q(self, s, q):
        return self.__sq_func(s=s*1e3, q=q, gas=self.__mat_name,fraction=self.__fraction, **self.__params)