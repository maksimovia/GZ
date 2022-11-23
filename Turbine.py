import mat_properties as prop
import numpy as n


class turbine:
    def __init__(
        self,
        water,
        stream1,
        stream2,
        stream3,
        stream4,
        stream5,
        stream6,
        stream7,
        stream8,
        water_streams0,
        water_streams,
    ):

        # - для расчета Теплоперепадов
        self.water = water
        self.stream1 = stream1
        self.stream2 = stream2
        self.stream3 = stream3
        self.stream4 = stream4
        self.stream5 = stream5
        self.stream6 = stream6
        self.stream7 = stream7
        self.stream8 = stream8
        self.water_streams0 = water_streams0
        self.water_streams = water_streams

        # Номинальные параметры в точках
        self.Pvd0 = water_streams0.at[stream1, "P"]
        self.Hvd0 = water_streams0.at[stream1, "H"]
        self.Gvd0 = water_streams0.at[stream1, "G"]
        Vd0 = water.p_h(self.Pvd0, self.Hvd0)
        self.Svd0 = Vd0["s"]
        self.Vvd0 = 1 / Vd0["rho"]

        self.Pvd_out0 = water_streams0.at[stream2, "P"]
        self.Hvd_out0 = water_streams0.at[stream2, "H"]
        self.Gvd_out0 = water_streams0.at[stream2, "G"]
        Vd_out0 = water.p_h(self.Pvd_out0, self.Hvd_out0)
        self.Svd_out0 = Vd_out0["s"]
        self.Vvd_out0 = 1 / Vd_out0["rho"]

        self.Pnd0 = water_streams0.at[stream3, "P"]
        self.Hnd0 = water_streams0.at[stream3, "H"]
        self.Gnd0 = water_streams0.at[stream3, "G"]
        Nd0 = water.p_h(self.Pnd0, self.Hnd0)
        self.Snd0 = Nd0["s"]
        self.Vnd0 = 1 / Nd0["rho"]

        self.Psmesh0 = water_streams0.at[stream4, "P"]
        self.Gsmesh0 = self.Gnd0 + self.Gvd0
        self.Hsmesh0 = (
            self.Hnd0 * self.Gnd0 + self.Hvd_out0 * self.Gvd0
        ) / self.Gsmesh0
        Smesh0 = water.p_h(self.Psmesh0, self.Hsmesh0)
        self.Ssmesh0 = Smesh0["s"]
        self.Vsmesh0 = 1 / Smesh0["rho"]
        self.delta_Pnd_smesh0 = self.Pnd0 - self.Psmesh0

        self.Potb20 = water_streams0.at[stream5, "P"]
        self.Hotb20 = water_streams0.at[stream5, "H"]
        self.Gotb20 = water_streams0.at[stream5, "G"]
        Otb20 = water.p_h(self.Potb20, self.Hotb20)
        self.Sotb20 = Otb20["s"]
        self.Votb20 = 1 / Otb20["rho"]

        self.Potb10 = water_streams0.at[stream6, "P"]
        self.Hotb10 = water_streams0.at[stream6, "H"]
        self.Gotb10 = water_streams0.at[stream6, "G"]
        Otb10 = water.p_h(self.Potb10, self.Hotb10)
        self.Sotb10 = Otb10["s"]
        self.Votb10 = 1 / Otb10["rho"]

        self.Pin_cnd0 = water_streams0.at[stream7, "P"]
        self.Hin_cnd0 = water_streams0.at[stream7, "H"]
        self.Gin_cnd0 = water_streams0.at[stream7, "G"]
        In_cnd0 = water.p_h(self.Pin_cnd0, self.Hin_cnd0)
        self.Sin_cnd0 = In_cnd0["s"]
        self.Vin_cnd0 = 1 / In_cnd0["rho"]

        self.Pin_kond0 = water_streams0.at[stream8, "P"]
        self.Hin_kond0 = water_streams0.at[stream8, "H"]
        self.Gin_kond0 = water_streams0.at[stream8, "G"]
        In_kond0 = water.p_h(self.Pin_kond0, self.Hin_kond0)
        self.Sin_kond0 = In_kond0["s"]
        self.Vin_kond0 = 1 / In_kond0["rho"]

        Hvd_outt0 = water.p_s(self.Pvd_out0, self.Svd0)["h"]
        self.effiency0_ots1 = (self.Hvd0 - self.Hvd_out0) / \
            (self.Hvd0 - Hvd_outt0)

        Hotb2t0 = water.p_s(self.Potb20, self.Ssmesh0)["h"]
        self.effiency0_ots2 = (self.Hsmesh0 - self.Hotb20) / \
            (self.Hsmesh0 - Hotb2t0)

        Hotb1t0 = water.p_s(self.Potb10, self.Sotb20)["h"]
        self.effiency0_ots3 = (self.Hotb20 - self.Hotb10) / \
            (self.Hotb20 - Hotb1t0)

        Hin_kondt0 = water.p_s(self.Pin_kond0, self.Sin_cnd0)["h"]
        self.effiency0_ots4 = (self.Hin_cnd0 - self.Hin_kond0) / (
            self.Hin_cnd0 - Hin_kondt0
        )
        self.Tair = self.water_streams.at['AIR', "T"]
        Hkond_out0 = water_streams0.at['KOND-KN', "H"]
        self.Q_kond0 = self.Gin_kond0 * (self.Hin_kond0-Hkond_out0)*(10**-3)

    def expansion(self, p1, h1, p2, eff):
        s1 = self.water.p_h(p1, h1)["s"]
        s2t = s1
        h2t = self.water.p_s(p2, s1)["h"]
        h2 = h1 - (h1 - h2t) * eff
        end = self.water.p_h(p2, h2)
        s2 = end["s"]
        t2 = end["T"]
        res = {"h": h2, "s": s2, "t": t2}
        return res

    def stodola_flugel(self, D0, D1, pn0, pv0, Vv0, pn1, pv1, Vv1):
        under_square=pn1**2 + ((D1 / D0) ** 2) * (pv0**2 - pn0**2) *(pv1 * Vv1 / (pv0 * Vv0))
        if under_square<0:
            pv1=pv0
            print("При расчете Стодола-флюгеля давление в верхнем отборе оказалось меньше 0, а именно корень из:",under_square)
        else:
            pv1 = n.sqrt(pn1**2 + ((D1 / D0) ** 2) * (pv0**2 - pn0**2) *(pv1 * Vv1 / (pv0 * Vv0)))
        if pv1<=pn1:
            print(f"Давление по Стодола-Флюгеля не правильно рассчитано:pv1={pv1}< pn1={pn1}")
            delta=0.001
            print(f"Принимаем давление в верхнем отборе на delta={delta} больше чем в нижнем")
            pv1=pn1+delta
        return pv1

    def calculate_cond(self, tair, G):
        T1v = 21.215*n.exp(0.0123*tair)
        # Проверка по номиналу
        # T1v = 15
        Q = self.Gin_kond * \
            (self.Hin_kond-self.water.p_q(self.Pin_kond, 0)['h'])*(10**-3)
        q = Q/self.Q_kond0
        Pin_kond = max(0.001, ((-0.0174000000+0.0169740000*q+0.0036920000*T1v-0.0001400000*(T1v**2)+0.0000022900*(
            T1v**3))/(1-0.5925300000*q+0.1835860000*(q**2)-0.0173900000*T1v+0.0002330000*(T1v**2)))*0.09806650124809235)
        return Pin_kond

    def off_design_relative_efficiency(self, G0, Vin0, Vout0, G1, Vin1, Vout1):
        Vave0 = (Vin0+Vout0)/2
        Q0 = G0*Vave0
        Vave1 = (Vin1+Vout1)/2
        Q1 = G1*Vave1
        eff_0 = 0.962
        q = Q1/Q0
        eff_1 = (-0.2366*q**2+0.475*q+0.7236)/eff_0
        if eff_1<0:
            # print("Удельная эффективность в отсеке меньше нуля ",eff_1)
            eff_1=0
    
        return eff_1

    def off_design_relative_efficiency_CND(self, pin, Hin, pout, Hout, G0, Vin0, Vout0, G1, Vin1, Vout1):
        betta_vl = 0.15
        state0 = self.water.p_h(pin, Hin)
        s0 = state0["s"]
        x0 = state0["Q"]
        Hteor_vl = self.water.s_q(s0, 1)['h']
        Hteor_out = self.water.p_s(pout, s0)['h']
        H0 = Hin-Hteor_out
        if x0 < 0:
            y0 = 0
        else:
            y0 = (1-x0)*100
        xz = self.water.p_h(pout, Hout)["Q"]
        if xz < 0:
            yz = 0
        else:
            yz = (1-xz)*100

        H_vl = 0
        # все влажное
        if yz > 0 and y0 > 0:
            H_vl = Hin-Hteor_out
        # сухое на входе
        if y0 == 0 and yz > 0:
            H_vl = Hteor_vl-Hteor_out
        # все сухое
        if y0 == 0 and yz == 0:
            H_vl = 0

        K_vl = 1-0.4*(1-betta_vl)*(y0-yz)/100*H_vl/H0
        Hvs = 40*(Vout1*G1/G0/Vout0)**2
        KPD0i = 0.87*(1+(H0-400)/10000)*K_vl-Hvs/H0

        # Vave0 = (Vin0+Vout0)/2
        # Q0 = G0*Vave0
        # Vave1 = (Vin1+Vout1)/2
        # Q1 = G1*Vave1
        # q = Q1/Q0
        # Delta_eff = -1.0702*q**2+1.7951*q-0.6597
        # Efficiency_out = KPD0i-Delta_eff

        q = G1/G0
        Eff_massflow = 2.2965*q**3-5.9155*q**2+4.8421*q-0.2163
        if q > 1:
            Eff_massflow = 1
        Efficiency_out = KPD0i*Eff_massflow
#         print(Efficiency_out,'check')
        Efficiency_out = max(Efficiency_out, -0.2)
        Efficiency_out = min(Efficiency_out, 1)  # Опарин внес

        return Efficiency_out

    def retrive_values(self):
        # Параметры в точках из таблицы
        self.Pvd = self.water_streams.at[self.stream1, "P"]
        self.Hvd = self.water_streams.at[self.stream1, "H"]
        self.Gvd = self.water_streams.at[self.stream1, "G"]

        self.Pvd_out = self.water_streams.at[self.stream2, "P"]
        self.Hvd_out = self.water_streams.at[self.stream2, "H"]
        self.Gvd_out = self.water_streams.at[self.stream2, "G"]

        self.Pnd = self.water_streams.at[self.stream3, "P"]
        self.Hnd = self.water_streams.at[self.stream3, "H"]
        self.Gnd = self.water_streams.at[self.stream3, "G"]

        self.Psmesh = self.water_streams.at[self.stream4, "P"]
        self.Gsmesh = self.Gnd + self.Gvd
        self.Hsmesh = (self.Hnd * self.Gnd +
                       self.Hvd_out * self.Gvd) / self.Gsmesh

        self.Potb2 = self.water_streams.at[self.stream5, "P"]
        self.Hotb2 = self.water_streams.at[self.stream5, "H"]
        self.Gotb2 = self.water_streams.at[self.stream5, "G"]

        self.Potb1 = self.water_streams.at[self.stream6, "P"]
        self.Hotb1 = self.water_streams.at[self.stream6, "H"]
        self.Gotb1 = self.water_streams.at[self.stream6, "G"]

        self.Pin_cnd = self.water_streams.at[self.stream7, "P"]
        self.Hin_cnd = self.water_streams.at[self.stream7, "H"]
        self.Gin_cnd = self.water_streams.at[self.stream7, "G"]

        self.Pin_kond = self.water_streams.at[self.stream8, "P"]
        self.Hin_kond = self.water_streams.at[self.stream8, "H"]
        self.Gin_kond = self.water_streams.at[self.stream8, "G"]

    def calculate_turbine_expansion(self):

        # Уточняем давления в смешении, на выходе из турбины и НД
        self.Pnd = self.Psmesh
        self.Pvd_out = self.Psmesh

        # отсек 1
        self.Gvd_out = self.Gvd
        self.Vvd = 1 / self.water.p_h(self.Pvd, self.Hvd)["rho"]
        self.Vvd_out = 1 / self.water.p_h(self.Pvd_out, self.Hvd_out)["rho"]
        self.KPD_ots1 = min(1,self.effiency0_ots1*self.off_design_relative_efficiency(
            self.Gvd0, self.Vvd0, self.Vvd_out0, self.Gvd_out, self.Vvd, self.Vvd_out))
        ots1_out = self.expansion(
            self.Pvd, self.Hvd, self.Pvd_out, self.KPD_ots1)
        self.Hvd_out = ots1_out["h"]

        # расчет смешения
        self.Gsmesh = self.Gnd + self.Gvd
        self.Hsmesh = (self.Hnd * self.Gnd +
                       self.Hvd_out * self.Gvd) / self.Gsmesh
        if self.Hsmesh<0:
            print("self.Hsmesh<0 ")
            print("self.Hsmesh",self.Hsmesh)
            print("self.Hvd_out",self.Hvd_out)
            print("self.Hnd",self.Hnd)
            print("self.Gvd",self.Gvd)
            print("self.Gnd",self.Gnd)
            self.Hsmesh=max(self.Hvd_out,self.Hnd)
        

        # отсек 2
        self.Gotb2 = self.Gsmesh
        self.Vsmesh = 1 / self.water.p_h(self.Psmesh, self.Hsmesh)["rho"]
        self.Votb2 = 1 / self.water.p_h(self.Potb2, self.Hotb2)["rho"]
        self.KPD_ots2 = min(1,self.effiency0_ots2*self.off_design_relative_efficiency(
            self.Gsmesh0, self.Vsmesh0, self.Votb20, self.Gsmesh, self.Vsmesh, self.Votb2))
        ots2_out = self.expansion(
            self.Psmesh, self.Hsmesh, self.Potb2, self.KPD_ots2)
        self.Hotb2 = ots2_out["h"]

        # отсек 3
        self.Votb2 = 1 / self.water.p_h(self.Potb2, self.Hotb2)["rho"]
        self.Votb1 = 1 / self.water.p_h(self.Potb1, self.Hotb1)["rho"]
        self.KPD_ots3 = min(1,self.effiency0_ots3*self.off_design_relative_efficiency(
            self.Gotb10, self.Votb20, self.Votb10, self.Gotb1, self.Votb2, self.Votb1))
        if self.KPD_ots3<0:
            print("Эффективность 3 отсека меньше нуля ",0)
            self.KPD_ots3=0
        ots3_out = self.expansion(
            self.Potb2, self.Hotb2, self.Potb1, self.KPD_ots3)
        self.Hotb1 = ots3_out["h"]

        # диафрагма
        self.Pin_cnd = self.Potb1 - self.diafragma
        if self.Pin_cnd < 0:
            self.Pin_cnd = abs(self.Pin_cnd)
            print('Давление на входе в ЦНД меньше нуля при температуре воздуха:',
                  self.water_streams.at['AIR', 'T'])

        self.Hin_cnd = self.Hotb1

        # расчет конденсатора
        self.Pin_kond = self.calculate_cond(self.Tair, self.Gin_kond)

        # отсек 4
        self.Vin_cnd = 1 / self.water.p_h(self.Pin_cnd, self.Hin_cnd)["rho"]
        self.Vin_kond = 1 / self.water.p_h(self.Pin_kond, self.Hin_kond)["rho"]
        self.KPD_ots4 = self.off_design_relative_efficiency_CND(
            self.Pin_cnd, self.Hin_cnd, self.Pin_kond, self.Hin_kond, self.Gin_cnd0, self.Vin_cnd0, self.Vin_kond0, self.Gin_cnd, self.Vin_cnd, self.Vin_kond)
        # if self.KPD_ots4 < 0:
        #     print("Pin_cnd",self.Pin_cnd, "Hin_cnd",self.Hin_cnd, "Pin_kond",self.Pin_kond, "KPD_ots4",self.KPD_ots4)
        ots4_out = self.expansion(
            self.Pin_cnd, self.Hin_cnd, self.Pin_kond, self.KPD_ots4)
        self.Hin_kond = ots4_out["h"]

        return

    def calculate_turbine_stodol_flugel(self):

        # отсек 4
        # print('G_ots4',self.Gin_cnd)
        self.Vin_cnd = 1 / self.water.p_h(self.Pin_cnd, self.Hin_cnd)["rho"]
        self.Pin_cnd = self.stodola_flugel(
            self.Gin_cnd0,
            self.Gin_cnd,
            self.Pin_kond0,
            self.Pin_cnd0,
            self.Vin_cnd0,
            self.Pin_kond,
            self.Pin_cnd,
            self.Vin_cnd,
        )

        # диафрагма
        self.Potb1 = self.Pin_cnd + self.diafragma

        # отсек 3
        # print('G_ots3',self.Gotb1)
        self.Votb2 = 1 / self.water.p_h(self.Potb2, self.Hotb2)["rho"]
        self.Potb2 = self.stodola_flugel(
            self.Gotb10,
            self.Gotb1,
            self.Potb10,
            self.Potb20,
            self.Votb20,
            self.Potb1,
            self.Potb2,
            self.Votb2,
        )

        # D0, D1, pn0, pv0, Vv0, pn1, pv1, Vv1
        # отсек 2
        # print('G_ots2',self.Gsmesh)
        self.Vsmesh = 1 / self.water.p_h(self.Psmesh, self.Hsmesh)["rho"]
        self.Psmesh = self.stodola_flugel(
            self.Gsmesh0,
            self.Gsmesh,
            self.Potb20,
            self.Psmesh0,
            self.Vsmesh0,
            self.Potb2,
            self.Psmesh,
            self.Vsmesh,
        )

        # Уточняем давления в смешении, на выходе из турбины и НД
        self.Vnd = 1 / self.water.p_h(self.Pnd, self.Hnd)["rho"]
        self.Pnd = self.Psmesh + self.delta_Pnd_smesh0 * \
            (self.Gnd/self.Gnd0)**2
        # self.Pnd = self.Psmesh + self.delta_Pnd_smesh0 * (self.Gnd*self.Vnd/self.Gnd0*self.Vnd0)**2*self.Vnd/self.Vnd0
        self.Pvd_out = self.Psmesh

        # отсек 1
        self.Vvd = 1 / self.water.p_h(self.Pvd, self.Hvd)["rho"]
        self.Pvd = self.stodola_flugel(
            self.Gvd0,
            self.Gvd,
            self.Pvd_out0,
            self.Pvd0,
            self.Vvd0,
            self.Pvd_out,
            self.Pvd,
            self.Vvd,
        )

        # stodola_flugel(self, D0, D1, pn0, pv0, Vv0, pn1, pv1, Vv1):
        return

    def calculate(self, diafragma, maxiterations, calctolerance):

        # Передаем параметры дифрагмы (потерю давления например)
        self.diafragma = diafragma
        self.retrive_values()
        for i in range(maxiterations):
            self.calculate_turbine_expansion()
            H_out = [
                self.Hvd,
                self.Hvd_out,
                self.Hnd,
                self.Hsmesh,
                self.Hotb2,
                self.Hotb1,
                self.Hin_cnd,
                self.Hin_kond]
            P_out = [
                self.Pvd,
                self.Pvd_out,
                self.Pnd,
                self.Psmesh,
                self.Potb2,
                self.Potb1,
                self.Pin_cnd,
                self.Pin_kond]
            P_out1 = P_out.copy()
            Eff_0 = [
                self.effiency0_ots1,
                self.effiency0_ots2,
                self.effiency0_ots3,
                self.effiency0_ots4]
            Eff_out = [
                self.KPD_ots1,
                self.KPD_ots2,
                self.KPD_ots3,
                self.KPD_ots4]
            # print("P_out_exp",P_out)
            # print("H_out",H_out)
            # print("Eff_out",Eff_out)
            # print("Eff_0",Eff_0)
            self.calculate_turbine_stodol_flugel()
            

            P_out = [
                self.Pvd,
                self.Pvd_out,
                self.Pnd,
                self.Psmesh,
                self.Potb2,
                self.Potb1,
                self.Pin_cnd,
                self.Pin_kond]
            # print("P_out_stodol",P_out)
            
            P_out2 = P_out.copy()
            Errors = list(map(lambda x, y: abs((x - y) / x * 100), P_out1, P_out2))
            
            Max_error = max(Errors)
            if abs(Max_error) < calctolerance and i > 1:
                print("Максимальная погрешность определения давления в отборах", Max_error)
                if min(H_out) < 0:
                    print(Eff_out)
                    print(H_out)
                break
            # if i == maxiterations-1:
                # print('Достигнуто максимальное количество итераций турбины')
                # print("P_out2",P_out2)
        G_out = [self.Gvd,
                 self.Gvd_out,
                 self.Gnd,
                 self.Gsmesh,
                 self.Gotb2,
                 self.Gotb1,
                 self.Gin_cnd,
                 self.Gin_kond, ]
        self.water_streams.loc[self.stream1:self.stream8, "G"] = G_out
        self.water_streams.loc[self.stream1:self.stream8, "P"] = P_out2
        self.water_streams.loc[self.stream1:self.stream8, "H"] = H_out
        Temperatures = list(
            map(lambda p, h: self.water.p_h(p, h)["T"], P_out2, H_out))
        Entropies = list(
            map(lambda p, h: self.water.p_h(p, h)["s"], P_out2, H_out))
        # Humidities=list(map( lambda p,h: round((1-self.water.p_h(p, h)["Q"])*100,3),P_out2,H_out))
        Humidities = list(
            map(lambda p, h: self.water.p_h(p, h)["Q"], P_out2, H_out))
        Humidities = list(map(lambda x: 100 if x < 0 else x*100, Humidities))
        self.water_streams.loc[self.stream1:self.stream8, "T"] = Temperatures
        self.water_streams.loc[self.stream1:self.stream8, "S"] = Entropies
        self.water_streams.loc[self.stream1:self.stream8, "X"] = Humidities

        return Eff_out

    def calculate_power(self):
        self.retrive_values()

        # отсек 1
        N1 = self.Gvd*(self.Hvd-self.Hvd_out)/1000

        # отсек 2
        N2 = self.Gsmesh*(self.Hsmesh-self.Hotb2)/1000

        # отсек 3
        N3 = self.Gotb1*(self.Hotb2-self.Hotb1)/1000

        # отсек 4
        N4 = self.Gin_cnd*(self.Hin_cnd-self.Hin_kond)/1000

        Ni_all = N1+N2+N3+N4

        return {'Ni': Ni_all, 'Ni1': N1, 'Ni2': N2, 'Ni3': N3, 'Ni4': N4}
