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
        self.effiency0_ots1 = (self.Hvd0 - self.Hvd_out0) / (self.Hvd0 - Hvd_outt0)

        Hotb2t0 = water.p_s(self.Potb20, self.Ssmesh0)["h"]
        self.effiency0_ots2 = (self.Hsmesh0 - self.Hotb20) / (self.Hsmesh0 - Hotb2t0)

        Hotb1t0 = water.p_s(self.Potb10, self.Sotb20)["h"]
        self.effiency0_ots3 = (self.Hotb20 - self.Hotb10) / (self.Hotb20 - Hotb1t0)

        Hin_kondt0 = water.p_s(self.Pin_kond0, self.Sin_cnd0)["h"]
        self.effiency0_ots4 = (self.Hin_cnd0 - self.Hin_kond0) / (
            self.Hin_cnd0 - Hin_kondt0
        )

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
        pv1 = n.sqrt(
            pn1**2 + ((D1 / D0) ** 2) * (pv0**2 - pn0**2) * pv1 * Vv1 / pv0 / Vv0
        )
        return pv1

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
        self.Hsmesh = (self.Hnd * self.Gnd + self.Hvd_out * self.Gvd) / self.Gsmesh

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
        KPD_ots1 = self.effiency0_ots1
        ots1_out = self.expansion(self.Pvd, self.Hvd, self.Pvd_out, KPD_ots1)
        self.Hvd_out = ots1_out["h"]

        # расчет смешения
        self.Gsmesh = self.Gnd + self.Gvd
        self.Hsmesh = (self.Hnd * self.Gnd + self.Hvd_out * self.Gvd) / self.Gsmesh

        # отсек 2
        KPD_ots2 = self.effiency0_ots2
        ots2_out = self.expansion(self.Psmesh, self.Hsmesh, self.Potb2, KPD_ots2)
        self.Hotb2 = ots2_out["h"]

        # отсек 3
        KPD_ots3 = self.effiency0_ots3
        ots3_out = self.expansion(self.Potb2, self.Hotb2, self.Potb1, KPD_ots3)
        self.Hotb1 = ots3_out["h"]

        # диафрагма
        self.Pin_cnd = self.Potb1 - self.diafragma
        self.Hin_cnd = self.Hotb1

        # отсек 4
        KPD_ots4 = self.effiency0_ots4
        ots4_out = self.expansion(self.Pin_cnd, self.Hin_cnd, self.Pin_kond, KPD_ots4)
        self.Hin_kond = ots4_out["h"]
        return

    def calculate_turbine_stodol_flugel(self):

        # отсек 4
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
        self.Votb2 = 1 / self.water.p_h(self.Potb2, self.Hotb2)["rho"]
        self.Potb2 = self.stodola_flugel(
            self.Gotb20,
            self.Gotb20,
            self.Potb10,
            self.Potb20,
            self.Votb20,
            self.Potb1,
            self.Potb2,
            self.Votb2,
        )

        # отсек 2
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
        self.Pnd = self.Psmesh + self.delta_Pnd_smesh0 * 1
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
                self.Hin_kond,
            ]
            P_out = [
                self.Pvd,
                self.Pvd_out,
                self.Pnd,
                self.Psmesh,
                self.Potb2,
                self.Potb1,
                self.Pin_cnd,
                self.Pin_kond,
            ]
            P_out1=P_out.copy()
            Eff_out = [
                self.effiency0_ots1,
                self.effiency0_ots2,
                self.effiency0_ots3,
                self.effiency0_ots4,
            ]
            self.calculate_turbine_stodol_flugel()
            P_out = [
                self.Pvd,
                self.Pvd_out,
                self.Pnd,
                self.Psmesh,
                self.Potb2,
                self.Potb1,
                self.Pin_cnd,
                self.Pin_kond,
            ]
            P_out2=P_out.copy()
            Errors = list(map(lambda x, y: (x - y) / x * 100, P_out1, P_out2))
            Max_error = max(Errors)
            if abs(Max_error) < calctolerance:
                print("Максимальная погрешность определения давления в отборах", Max_error)
                break
        self.water_streams.loc[self.stream1:self.stream8, "P"]=P_out2
        self.water_streams.loc[self.stream1:self.stream8, "H"]=H_out
        Temperatures=list(map( lambda p,h: self.water.p_h(p, h)["T"],P_out2,H_out))
        Entropies=list(map( lambda p,h: self.water.p_h(p, h)["s"],P_out2,H_out))
        # Humidities=list(map( lambda p,h: round((1-self.water.p_h(p, h)["Q"])*100,3),P_out2,H_out))
        Humidities=list(map( lambda p,h: self.water.p_h(p, h)["Q"],P_out2,H_out))
        Humidities=list(map( lambda x: 100 if x<0 else x*100,Humidities))
        self.water_streams.loc[self.stream1:self.stream8, "T"]=Temperatures
        self.water_streams.loc[self.stream1:self.stream8, "S"]=Entropies
        self.water_streams.loc[self.stream1:self.stream8, "X"]=Humidities

        
        

        return {'Max_error':Max_error, 'P_out':P_out2, 'Eff_out':Eff_out}