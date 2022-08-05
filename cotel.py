import numpy as n
from scipy.optimize import root

class heatex:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,calctolerance,gas0,gas1,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.gas=gas1
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P021 = water_streams0.at[stream21,'P']
        self.P022 = water_streams0.at[stream22,'P']
        self.H11  = gas_streams.at[stream11,'H']
        self.H21  = water_streams.at[stream21,'H']
        self.G1   = gas_streams.at[stream11,'G']
        self.G2   = water_streams.at[stream21,'G']
        self.P1   = gas_streams.at[stream11,'P']
        self.P21  = water_streams.at[stream21,'P']
    def calc(self):
        Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P021,self.H021)['T']
        T022   = self.water.p_h(self.P022,self.H022)['T']
        dTmin0 = min(T011-T022,T012-T021)
        dTmax0 = max(T011-T022,T012-T021)
        LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        T02av  = (T021+T022)/2
        P02av  = (self.P021+self.P022)/2
        lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        ro01av = self.gas0.p_t(self.P01,T01av)['rho']
        ro02av = self.water.p_t(P02av,T02av)['rho']
        ro021  = self.water.p_q(self.P021,1)['rho']        
        ro21   = self.water.p_q(self.P21,1)['rho']
        ddp    = (ro21/ro021)*((self.G02/self.G2)**2)
        P22    = self.P21 - ((self.P021-self.P022)/ddp)
        P12 = self.P1
        def T12sved(T12):           
            T11 = self.gas.p_h(self.P1,self.H11)['T']
            T21 = self.water.p_h(self.P21,self.H21)['T']
            H12 = self.gas.p_t(self.P1,T12)['h']
            Q = self.G1*(self.H11-H12)*self.KPD
            H22 = self.H21 + (Q/self.G2)
            T22 = self.water.p_h(P22,H22)['T']
            dTmin= min(T11-T22,T12-T21)
            dTmax= max(T11-T22,T12-T21)
            if dTmin<0 or dTmax<0 or dTmin==dTmax:
                LMTD = (dTmax+dTmin)/2
            else:
                LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))            
            dt = LMTD0/LMTD
            T1av = (T11+T12)/2
            lambda1av = self.gas.p_t(self.P1,T1av)['L']
            Pr1av = self.gas.p_t(self.P1,T1av)['Prandtl']
            nu1av = self.gas.p_t(self.P1,T1av)['nu']
            ro1av = self.gas.p_t(self.P1,T1av)['rho']
            kk = (lambda01av/lambda1av)*((Pr01av/Pr1av)**0.33)*(((self.G01/self.G1)*(ro1av/ro01av)*(nu1av/nu01av))**0.685)
            Qq = Q0 / (kk*dt)
            return ((Q-Qq)/Q)
        sol = root(T12sved, T012, method=self.calcmethod, tol=self.calctolerance)
        T12=float(sol.x)
        H12 = self.gas.p_t(self.P1,T12)['h']
        Q = self.G1*(self.H11-H12)*self.KPD
        H22 = self.H21 + (Q/self.G2)
        T22 = self.water.p_h(P22,H22)['T']
        return [T12,P12,H12,self.G1,T22,P22,H22,self.G2,Q]
class vapor:
    def __init__(self, stream11, stream12,stream21,stream22,KPD,calctolerance,gas,gas0,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.gas=gas
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P02  = water_streams0.at[stream21,'P']
        self.H11  = gas_streams.at[stream11,'H']
        self.H21  = water_streams.at[stream21,'H']
        self.G1   = gas_streams.at[stream11,'G']
        self.P1   = gas_streams.at[stream11,'P']
        self.P2   = water_streams.at[stream21,'P']
    def calc(self):
        Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P02,self.H021)['T']
        T022   = self.water.p_h(self.P02,self.H022)['T']
        T0np   = self.water.p_q(self.P2,1)['T']        
        dTmin0 = min(T011-T0np,T012-T0np)
        dTmax0 = max(T011-T0np,T012-T0np)
        LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        ro01av = self.gas0.p_t(self.P01,T01av)['rho']
                
        def T12sved(T12):
            T11 = self.gas.p_h(self.P1,self.H11)['T']
            T21 = self.water.p_h(self.P2,self.H21)['T']
            H12 = self.gas.p_t(self.P1,T12)['h']
            Q = self.G1*(self.H11-H12)*self.KPD
            H22 = self.water.p_q(self.P2,1)['h']  
            T22 = self.water.p_h(self.P2,H22)['T']
            G2 = Q/(H22-self.H21)
            dTmin= min(T11-T22,T12-T21)
            dTmax= max(T11-T22,T12-T21)
            if dTmin<0 or dTmax<0 or dTmin==dTmax:
                LMTD = (dTmax+dTmin)/2
            else:
                LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))            
            dt = LMTD0/LMTD
            T1av = (T11+T12)/2
            lambda1av = self.gas.p_t(self.P1,T1av)['L']
            Pr1av = self.gas.p_t(self.P1,T1av)['Prandtl']
            nu1av = self.gas.p_t(self.P1,T1av)['nu']
            ro1av = self.gas.p_t(self.P1,T1av)['rho']
            kk = (lambda01av/lambda1av)*((Pr01av/Pr1av)**0.33)*(((self.G01/self.G1)*(ro1av/ro01av)*(nu1av/nu01av))**0.685)
            Qq = Q0 / (kk*dt)
            return ((Q-Qq)/Q)
        sol = root(T12sved, T012, method=self.calcmethod, tol=self.calctolerance)
        T12=float(sol.x)
        H12 = self.gas.p_t(self.P1,T12)['h']
        Q = self.G1*(self.H11-H12)*self.KPD
        H22 = self.water.p_q(self.P2,1)['h']
        T22 = self.water.p_h(self.P2,H22)['T']
        G2 = Q/(H22-self.H21)
        return [T12,self.P1,H12,self.G1,T22,self.P2,H22,G2,Q]         
class vaporND:
    def __init__(self, stream11, stream12,stream21,stream22,KPD,calctolerance,gas,gas0,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.gas=gas
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P02  = water_streams0.at[stream21,'P']
        self.H11  = gas_streams.at[stream11,'H']
        self.H21  = water_streams.at[stream21,'H']
        self.G1   = gas_streams.at[stream11,'G']
        self.P1   = gas_streams.at[stream11,'P']
        self.P2   = water_streams.at[stream21,'P']
    def calc(self):
        Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P02,self.H021)['T']
        T022   = self.water.p_h(self.P02,self.H022)['T']
        T0np   = self.water.p_q(self.P2,1)['T']        
        dTmin0 = min(T011-T0np,T012-T0np)
        dTmax0 = max(T011-T0np,T012-T0np)
        LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        ro01av = self.gas0.p_t(self.P01,T01av)['rho']
        def T12sved(T12):
            T11 = self.gas.p_h(self.P1,self.H11)['T']
            T21 = self.water.p_h(self.P2,self.H21)['T']
            H12 = self.gas.p_t(self.P1,T12)['h']
            Q = self.G1*(self.H11-H12)*self.KPD
            H22 = self.water.p_q(self.P2,1)['h']  
            T22 = self.water.p_h(self.P2,H22)['T']
            Dvd = self.water_streams.at['PEVD-DROS','G']
            Hvd =self.water.p_q(self.P2,0)['h']
            G2 = (Q-Dvd*(Hvd-self.H21))/(H22-self.H21)
            dTmin= min(T11-T22,T12-T21)
            dTmax= max(T11-T22,T12-T21)
            if dTmin<0 or dTmax<0 or dTmin==dTmax:
                LMTD = (dTmax+dTmin)/2
            else:
                LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))            
            dt = LMTD0/LMTD
            T1av = (T11+T12)/2
            lambda1av = self.gas.p_t(self.P1,T1av)['L']
            Pr1av = self.gas.p_t(self.P1,T1av)['Prandtl']
            nu1av = self.gas.p_t(self.P1,T1av)['nu']
            ro1av = self.gas.p_t(self.P1,T1av)['rho']
            kk = (lambda01av/lambda1av)*((Pr01av/Pr1av)**0.33)*(((self.G01/self.G1)*(ro1av/ro01av)*(nu1av/nu01av))**0.685)
            Qq = Q0 / (kk*dt)
            return ((Q-Qq)/Q)
        sol = root(T12sved, T012, method=self.calcmethod, tol=self.calctolerance)
        T12=float(sol.x)
        H12 = self.gas.p_t(self.P1,T12)['h']
        Q = self.G1*(self.H11-H12)*self.KPD
        H22 = self.water.p_q(self.P2,1)['h']
        T22 = self.water.p_h(self.P2,H22)['T']
        Dvd = self.water_streams.at['PEVD-DROS','G']
        Hvd =self.water.p_q(self.P2,0)['h']
        G2 = (Q- Dvd*(Hvd-self.H21))/(H22-self.H21)
        Hvd = self.water.p_q(self.P2,0)['h']
        Tvd = self.water.p_q(self.P2,0)['T']
        Pvd = self.P2
        Dvd = self.water_streams.at['PEVD-DROS','G']
        
        return [T12,self.P1,H12,self.G1,T22,self.P2,H22,G2,Q,Tvd,Pvd,Hvd,Dvd]