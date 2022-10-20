import numpy as n
from scipy.optimize import root
import time
import nasos


class heatex:
    def __init__(self, stream11, stream12, stream21, stream22, KPD, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.heaters = heaters
        self.gas = gas1
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P021 = water_streams0.at[stream21, 'P']
        self.P022 = water_streams0.at[stream22, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_h(self.P021, self.H021)['T']
        T022 = self.water.p_h(self.P022, self.H022)['T']
        dTmin0 = min(T011-T022, T012-T021)
        dTmax0 = max(T011-T022, T012-T021)
        
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        T02av = (T021+T022)/2
        P02av = (self.P021+self.P022)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']
        self.ro02av = self.water.p_t(P02av, T02av)['rho']
        self.ro021 = self.water.p_q(self.P021, 1)['rho']

    def calc(self, calctolerance):
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        G2 = self.water_streams.at[self.stream21, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P22 = self.water_streams.at[self.stream22, 'P']

        def dPsved(P21):
            if P21 < P22 or P21 > 2*P22:
                return 10**9
            else:
                ro21 = self.water.p_q(P21, 1)['rho']
                ddp = (ro21/self.ro021)*((self.G02/G2)**2)
                P22new = P21 - ((self.P021-self.P022)/ddp)
                return P22-P22new
        sol = root(dPsved, 1.1*P22, tol=10**-8)
        P21 = float(sol.x)
        ro21 = self.water.p_q(P21, 1)['rho']
        P12 = P1
        T21 = self.water.p_h(P21, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']

        def T12sved(T12):
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = H21 + (Q/G2)
                T22 = self.water.p_h(P22, H22)['T']
                dTmin = min(T11-T22, T12-T21)
                dTmax = max(T11-T22, T12-T21)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
        sol = root(T12sved, max(Tfirst, T21+5),
                   method=self.calcmethod, tol=calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        Qg = G1*(H11-H12)
        H22 = H21 + (Q/G2)
        T22 = self.water.p_h(P22, H22)['T']
        return {'Tg': T12, 'Pg': P12, 'Hg': H12, 'Gg': G1, 'Qg': Qg, 'Tw': T22, 'Pw1': P21, 'Pw2': P22, 'Hw': H22, 'Gw': G2, 'Qw': Q, 'KPD': self.KPD}


class heatexPEND:
    def __init__(self, stream11, stream12, stream21, stream22, KPD,  gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.heaters = heaters
        self.gas = gas1
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P021 = water_streams0.at[stream21, 'P']
        self.P022 = water_streams0.at[stream22, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_h(self.P021, self.H021)['T']
        T022 = self.water.p_h(self.P022, self.H022)['T']
        dTmin0 = min(T011-T022, T012-T021)
        dTmax0 = max(T011-T022, T012-T021)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        T02av = (T021+T022)/2
        P02av = (self.P021+self.P022)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']
        self.ro02av = self.water.p_t(P02av, T02av)['rho']
        self.ro021 = self.water.p_q(self.P021, 1)['rho']

    def calc(self, calctolerance):
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        G2 = self.water_streams.at[self.stream21, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P22 = self.water_streams.at[self.stream22, 'P']

        def dPsved(P21):
            if P21 < P22 or P21 > 2*P22:
                return 10**9
            else:
                ro21 = self.water.p_q(P21, 1)['rho']
                ddp = (ro21/self.ro021)*((self.G02/G2)**2)
                P22new = P21 - ((self.P021-self.P022)/ddp)
                return P22-P22new
        sol = root(dPsved, 1.1*P22, tol=10**-8)
        P21 = float(sol.x)
        ro21 = self.water.p_q(P21, 1)['rho']
        P12 = P1
        T21 = self.water.p_h(P21, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']

        def T12sved(T12):
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = H21 + (Q/G2)
                T22 = self.water.p_h(P22, H22)['T']
                dTmin = min(T11-T22, T12-T21)
                dTmax = max(T11-T22, T12-T21)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
        sol = root(T12sved, T11*0.99,
                   method=self.calcmethod, tol=calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        Qg = G1*(H11-H12)
        H22 = H21 + (Q/G2)
        T22 = self.water.p_h(P22, H22)['T']
        return {'Tg': T12, 'Pg': P12, 'Hg': H12, 'Gg': G1, 'Qg': Qg, 'Tw': T22, 'Pw1': P21, 'Pw2': P22, 'Hw': H22, 'Gw': G2, 'Qw': Q, 'KPD': self.KPD}


class evaporVD:
    def __init__(self, stream11, stream12, stream21, stream22, KPD, gas, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.heaters = heaters
        self.gas = gas
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P02 = water_streams0.at[stream21, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_q(self.P02, 0)['T']
        T022 = self.water.p_h(self.P02, self.H022)['T']
        T0np = self.water.p_q(self.P02, 1)['T']
        dTmin0 = min(T011-T0np, T012-T0np)
        dTmax0 = max(T011-T0np, T012-T0np)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']

    def calc(self, calctolerance):
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P2 = self.water_streams.at[self.stream21, 'P']
        T21 = self.water.p_h(P2, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']
        T11 = self.gas.p_h(P1, H11)['T']
        T21 = self.water.p_q(P2, 0)['T']

        def T12sved(T12):
            T12 = float(T12)
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = self.water.p_q(P2, 1)['h']
                T22 = self.water.p_h(P2, H22)['T']
                G2 = Q/(H22-H21)
                dTmin = min(T11-T22, T12-T21)
                dTmax = max(T11-T22, T12-T21)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
        sol = root(T12sved, max(Tfirst, T21+5),
                   method=self.calcmethod, tol=calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        Qg = G1*(H11-H12)
        H22 = self.water.p_q(P2, 1)['h']
        T22 = self.water.p_h(P2, H22)['T']
        G2 = Q/(H22-H21)
        return {'Tg': T12, 'Pg': P1, 'Hg': H12, 'Gg': G1, 'Qg': Qg, 'Tw': T22, 'Pw': P2, 'Hw': H22, 'Gw': G2, 'Qw': Q, 'KPD': self.KPD}


class evaporND:
    def __init__(self, stream11, stream12, stream21, stream22, streamVD, KPD, gas, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.streamVD = streamVD
        self.heaters = heaters
        self.gas = gas
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P02 = water_streams0.at[stream21, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_q(self.P02, 0)['T']
        T022 = self.water.p_h(self.P02, self.H022)['T']
        T0np = self.water.p_q(self.P02, 1)['T']
        dTmin0 = min(T011-T0np, T012-T0np)
        dTmax0 = max(T011-T0np, T012-T0np)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']

    def calc(self, calctolerance):
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P2 = self.water_streams.at[self.stream21, 'P']
        T21 = self.water.p_q(P2, 0)['T']
        T11 = self.gas.p_h(P1, H11)['T']
        Dvd = self.water_streams.at[self.streamVD, 'G']

        def T12sved(T12):
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = self.water.p_q(P2, 1)['h']
                T22 = self.water.p_h(P2, H22)['T']
                Hvd = self.water.p_q(P2, 0)['h']
                G2 = (Q-Dvd*(Hvd-H21))/(H22-H21)
                dTmin = min(T11-T22, T12-T22)
                dTmax = max(T11-T22, T12-T22)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
        sol = root(T12sved, max(Tfirst, T21+5),
                   method=self.calcmethod, tol=calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        Qg = G1*(H11-H12)
        H22 = self.water.p_q(P2, 1)['h']
        T22 = self.water.p_h(P2, H22)['T']
        Hvd = self.water.p_q(P2, 0)['h']
        G2 = (Q - Dvd*(Hvd-H21))/(H22-H21)
        Tvd = self.water.p_q(P2, 0)['T']
        Pvd = P2
        return {'Tg': T12, 'Pg': P1, 'Hg': H12, 'Gg': G1, 'Qg': Qg, 'Tw': T22, 'Pw': P2, 'Hw': H22, 'Gw': G2, 'Qw': Q, 'KPD': self.KPD, 'Tvd': Tvd, 'Pvd': Pvd, 'Hvd': Hvd, 'Dvd': Dvd}


class cotel_all:
    def __init__(self, KPD, KPDnasos,  gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters, electric):
        self.PEVD_obj = heatex('GTU-PEVD', 'PEVD-IVD', 'IVD-PEVD', 'PEVD-DROSVD',
                               KPD,  gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters)

        self.IVD_obj = evaporVD('PEVD-IVD', 'IVD-EVD', 'EVD-IVD', 'IVD-PEVD',
                                KPD,  gas1, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters)

        self.EVD_obj = heatex('IVD-EVD', 'EVD-PPND', 'PEN-EVD', 'EVD-IVD',
                              KPD,  gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters)

        self.PEN_obj = nasos.nasos(
            'BND-PEN', 'PEN-EVD', water, KPDnasos, water_streams, water_streams0)

        self.PPND_obj = heatexPEND('EVD-PPND', 'PPND-IND', 'IND-PPND', 'PPND-DROSND',
                                   KPD,  gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters)
        
        self.IND_obj = evaporND('PPND-IND', 'IND-GPK', 'GPK-IND', 'IND-PPND',  'PEVD-DROSVD',
                                KPD, gas1, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters)

        self.GPK_obj = heatex('IND-GPK', 'GPK-out', 'REC-GPK', 'GPK-REC',
                              KPD,  gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams, heaters)
        self.KPD = KPD
        self.water_streams = water_streams
        self.gas_streams = gas_streams
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas1 = gas1
        self.water = water
        self.heaters = heaters
        self.electric = electric

    def calc(self, calctolerance, maxiterations=50):
        
        it = maxiterations
        start_time = time.time()
        
        
        self.gas_streams.loc['GTU-PEVD', 'H'] = self.gas1.p_t(
            self.gas_streams.loc['GTU-PEVD', 'P'], self.gas_streams.loc['GTU-PEVD', 'T'])['h']
        for k in range(it):
            if k>it-it/2:
                calctolerance_new=calctolerance/10
                print('Повышена точность расчета котла для увеличения сходимости')

            else:
                calctolerance_new=calctolerance
            # Связвка высокого давления
            for j in range(it):

                # Расчёт ПЕВД
                PEVD = self.PEVD_obj.calc(calctolerance_new)
                self.gas_streams.loc['PEVD-IVD', 'T':'G'] = [PEVD['Tg'],
                                                             PEVD['Pg'], PEVD['Hg'], PEVD['Gg']]
                self.water_streams.loc['PEVD-DROSVD', 'T':'G'] = [
                    PEVD['Tw'], PEVD['Pw2'], PEVD['Hw'], PEVD['Gw']]
                self.heaters.loc['PEVD', 'Qw':'KPD'] = [
                    PEVD['Qw'], PEVD['Qg'], PEVD['KPD']]
                self.water_streams.loc["IVD-PEVD":"PEN-EVD", "P"] = PEVD['Pw1']

                # Расчёт ИВД
                IVD = self.IVD_obj.calc(calctolerance_new)
                self.gas_streams.loc['IVD-EVD', 'T':'G'] = [IVD['Tg'],
                                                            IVD['Pg'], IVD['Hg'], IVD['Gg']]
                self.water_streams.loc['IVD-PEVD', 'T':'G'] = [
                    IVD['Tw'], IVD['Pw'], IVD['Hw'], IVD['Gw']]
                self.heaters.loc['IVD', 'Qw':'KPD'] = [
                    IVD['Qw'], IVD['Qg'], IVD['KPD']]

                # Переопределение расхода в ВД
                self.water_streams.loc['PEVD-DROSVD':'PEN-EVD',
                                       'G'] = IVD['Gw']
                self.water_streams.loc['BND-PEN', 'G'] = IVD['Gw']

                # Расчёт ЭВД
                EVD = self.EVD_obj.calc(calctolerance_new)

                self.gas_streams.loc['EVD-PPND', 'T':'G'] = [EVD['Tg'],
                                                             EVD['Pg'], EVD['Hg'], EVD['Gg']]
                self.water_streams.loc['EVD-IVD', 'T':'G'] = [
                    EVD['Tw'], EVD['Pw2'], EVD['Hw'], EVD['Gw']]
                self.heaters.loc['EVD', 'Qw':'KPD'] = [
                    EVD['Qw'], EVD['Qg'], EVD['KPD']]

                # Баланс ПЕВ+ИВД+ЭВД
                Qgas1VD = self.KPD*self.gas_streams.at['GTU-PEVD', 'G'] * \
                    (self.gas_streams.at['GTU-PEVD', 'H'] -
                     self.gas_streams.at['IVD-EVD', 'H'])
                Qwat1VD = self.water_streams.at['PEVD-DROSVD', 'G']*(
                    self.water_streams.at['PEVD-DROSVD', 'H']-self.water_streams.at['EVD-IVD', 'H'])
                ErrorVD = (Qgas1VD-Qwat1VD)/Qgas1VD*100
                if k > it/2:
                    print('dQ/Q ПЕВД+ИВД+ЭВД', ErrorVD)
                if abs(ErrorVD) < calctolerance:
                    break
                if j == it - 1:
                    print(
                        "Достигнуто максимальное количество итераций контура высокого давления")
            # Для сходимости
            if k == 0:
                self.gas_streams.loc['PPND-IND',
                                     'T'] = self.gas_streams.loc['EVD-PPND', 'T'] - 3
                self.gas_streams.loc['PPND-IND', 'H'] = self.gas1.p_t(
                    self.gas_streams.loc['PPND-IND', 'P'], self.gas_streams.loc['PPND-IND', 'T'])['h']

            # Связка низкого давления
            for j in range(it):
                # Расчёт ППНД
                PPND = self.PPND_obj.calc(calctolerance_new)
                self.gas_streams.loc['PPND-IND', 'T':'G'] = [PPND['Tg'],
                                                             PPND['Pg'], PPND['Hg'], PPND['Gg']]
                self.water_streams.loc['PPND-DROSND', 'T':'G'] = [
                    PPND['Tw'], PPND['Pw2'], PPND['Hw'], PPND['Gw']]
                self.heaters.loc['PPND', 'Qw':'KPD'] = [
                    PPND['Qw'], PPND['Qg'], PPND['KPD']]
                self.water_streams.loc["IND-PPND":"SMESH-GPK",
                                       "P"] = PPND['Pw1']
                self.water_streams.loc["BND-PEN", "P"] = PPND['Pw1']

                # Расчёт ИНД
                IND = self.IND_obj.calc(calctolerance_new)
                self.gas_streams.loc['IND-GPK', 'T':'G'] = [IND['Tg'],
                                                            IND['Pg'], IND['Hg'], IND['Gg']]
                self.water_streams.loc['IND-PPND',
                                       'T':'G'] = [IND['Tw'], IND['Pw'], IND['Hw'], IND['Gw']]
                self.heaters.loc['IND', 'Qw':'KPD'] = [
                    IND['Qw'], IND['Qg'], IND['KPD']]

                # Переопределение расхода в НД
                self.water_streams.loc['PPND-DROSND':'IND-PPND',
                                       'G'] = IND['Gw']

                # ПЭН
                self.water_streams.loc['BND-PEN', 'T':'G'] = [IND['Tvd'],
                                                              IND['Pvd'], IND['Hvd'], IND['Dvd']]
                PEN = self.PEN_obj.calc()
                self.water_streams.loc['PEN-EVD',
                                       'T':'G'] = [PEN['T2'], PEN['P2'], PEN['h2real'], PEN['G1']]
                self.electric.loc['PEN', 'Ni':'KPD'] = [
                    PEN['Ni'], PEN['Ngm'], PEN['KPDm'], PEN['KPD']]

                # Баланс ППНД+ИНД
                Qgas = self.KPD*self.gas_streams.at['EVD-PPND', 'G'] * \
                    (self.gas_streams.at['EVD-PPND', 'H'] -
                     self.gas_streams.at['IND-GPK', 'H'])
                Qwat = self.water_streams.at['IND-PPND', 'G']*(self.water_streams.at['IND-PPND', 'H']-self.water_streams.at['GPK-IND', 'H']) +\
                    self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H'] -
                                                               self.water_streams.at['IND-PPND', 'H']) +\
                    self.water_streams.at['BND-PEN', 'G']*(
                        self.water_streams.at['BND-PEN', 'H']-self.water_streams.at['GPK-IND', 'H'])

                # Расчет ГПК
                for i in range(it):

                    # Расчёт ГПК
                    GPK = self.GPK_obj.calc(calctolerance_new)
                    self.gas_streams.loc['GPK-out', 'T':'G'] = [GPK['Tg'],
                                                                GPK['Pg'], GPK['Hg'], GPK['Gg']]
                    self.water_streams.loc['GPK-REC', 'T':'G'] = [
                        GPK['Tw'], GPK['Pw2'], GPK['Hw'], GPK['Gw']]
                    Qw_gpk1 = self.water_streams.at['GPK-IND', 'G']*(
                        self.water_streams.at['GPK-IND', 'H']-self.water_streams.at['SMESH-GPK', 'H'])
                    Qw_gpk2 = self.water_streams.at['GPK-REC', 'G']*(
                        self.water_streams.at['GPK-REC', 'H']-self.water_streams.at['REC-GPK', 'H'])
                    Error_gpk = (Qw_gpk1-Qw_gpk2)/Qw_gpk1*100

                    # Расчёт расхода в ГПК (рециркуляция)
                    tgpk_in = self.water_streams.at['SMESH-GPK', 'T']
                    tgpk_in0 = self.water_streams0.at['REC-GPK', 'T']
                    p_gpk = self.water_streams.at['REC-GPK', 'P']
                    h_gpk_in_rec = self.water_streams.at['SMESH-GPK', 'H']
                    h_gpk_in_60 = self.water.p_t(p_gpk, tgpk_in0)['h']
                    h_gpk_out = self.water_streams.at['GPK-REC', 'H']
                    G_all = self.water_streams.at['PPND-DROSND',
                                                  'G']+self.water_streams.at['PEVD-DROSVD', 'G']
                    if tgpk_in < tgpk_in0:
                        self.water_streams.at['REC-GPK', 'H'] = h_gpk_in_60
                        self.water_streams.at['REC-GPK', 'T'] = tgpk_in0
                        G_rec = G_all*(h_gpk_in_60-h_gpk_in_rec) / \
                            (h_gpk_out-h_gpk_in_60)
                    else:
                        G_rec = 0
                        self.water_streams.loc['REC-GPK',
                                              'T':'H'] = self.water_streams.loc['SMESH-GPK', 'T':'H']
                    G_gpk = G_all+G_rec
                    self.water_streams.at['REC-GPK', 'G'] = G_gpk
                    self.water_streams.loc['GPK-IND',
                                           'T':'H'] = self.water_streams.loc['GPK-REC', 'T':'H']
                    self.water_streams.at['GPK-IND', 'G'] = G_all
                    self.heaters.loc['GPK', 'Qw':'KPD'] = [
                        GPK['Qw'], GPK['Qg'], GPK['KPD']]
                    if abs(Error_gpk) < calctolerance_new:
                        break
                    if i == it - 1:
                        print(
                            "Достигнуто максимальное количество итераций контура ГПК")

                # Баланс ППНД+ИНД+ГПК
                Qgas1ND = self.KPD*self.gas_streams.at['EVD-PPND', 'G'] * \
                    (self.gas_streams.at['EVD-PPND', 'H'] -
                     self.gas_streams.at['GPK-out', 'H'])
                Qwat1ND = self.water_streams.at['GPK-IND', 'G']*(self.water_streams.at['GPK-IND', 'H']-self.water_streams.at['SMESH-GPK', 'H']) +\
                    self.water_streams.at['IND-PPND', 'G']*(self.water_streams.at['IND-PPND', 'H']-self.water_streams.at['GPK-IND', 'H']) +\
                    self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H']-self.water_streams.at['IND-PPND', 'H']) +\
                    self.water_streams.at['BND-PEN', 'G'] * \
                    (self.water_streams.at['BND-PEN', 'H'] -
                     self.water_streams.at['GPK-IND', 'H'])
                Qwat2ND = self.water_streams.at['GPK-REC', 'G']*(self.water_streams.at['GPK-REC', 'H']-self.water_streams.at['REC-GPK', 'H']) +\
                    self.water_streams.at['IND-PPND', 'G']*(self.water_streams.at['IND-PPND', 'H']-self.water_streams.at['GPK-IND', 'H']) +\
                    self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H']-self.water_streams.at['IND-PPND', 'H']) +\
                    self.water_streams.at['BND-PEN', 'G'] * \
                    (self.water_streams.at['BND-PEN', 'H'] -
                     self.water_streams.at['GPK-IND', 'H'])
                ErrorND = (Qgas1ND-Qwat1ND)/Qgas1ND*100
                ErrorND2 = (Qgas1ND-Qwat2ND)/Qgas1ND*100
                if k > it/2:
                    print('dQ/Q ППНД+ИНД+ГПК', ErrorND)
                if abs(ErrorND) < calctolerance and abs(ErrorND2) < calctolerance:
                    break
                if j == it - 1:
                    print(
                        "Достигнуто максимальное количество итераций контура низкого давления")

            # Баланс общий
            Qgasall = self.KPD*self.gas_streams.at['GTU-PEVD', 'G'] * \
                (self.gas_streams.at['GTU-PEVD', 'H'] -
                 self.gas_streams.at['GPK-out', 'H'])
            Qwatall = self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H']-self.water_streams.at['SMESH-GPK', 'H'])+self.water_streams.at['PEVD-DROSVD', 'G']*(
                self.water_streams.at['PEVD-DROSVD', 'H']-self.water_streams.at['SMESH-GPK', 'H'])-self.water_streams.at['BND-PEN', 'G']*(self.water_streams.at['PEN-EVD', 'H']-self.water_streams.at['BND-PEN', 'H'])
            ErrorALL = (Qgasall-Qwatall)/Qgasall*100
            # print('dQ/Qsumm', ErrorALL)
            if abs((Qgasall-Qwatall)/Qgasall*100) < calctolerance:
                print("Fin котел-утилизатора:--- %s сек. ---" %
                      round((time.time() - start_time), 2))
                print('dQ/Qsumm', ErrorALL)
                # print('dQ/Qvd', ErrorVD)
                # print('dQ/Qnd', ErrorND)
                break
            if k == it - 1:
                print("Достигнуто максимальное количество итераций котла-утилизатора")
