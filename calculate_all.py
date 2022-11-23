

def calculate_CCGT(args):
    
    import os
    import time
    import GTU
    import KU_TU
    import mat_properties as prop
    import numpy as n
    import pandas as pd
    import SP
    import ASW
    from scipy.optimize import root
#     import calculate_all as calc
    from numpy import linalg as LA
    gas_streams0, water_streams0,GTU_ISO,GTU_input,gas_streams,water_streams, heaters, electric, Tnv, KPD_PN, KPD_KN, KPD_to, KPD_SP, Calcmethod,Calctolerance,Maxiterations_KU_TU,Maxiterations_cotel,Maxiterations_turbine,Сalculate_minimum,Teplo,steamVD_fraction_to_turbine,gasmix,RP,fractiongas0,gas0,water,fractionwaterMethane,waterMethanemix,waterMethane,Methane,gas_KU_PKM,PKM_zaryad,Temperatures,Pressure,Enthalpies,Gas_turbine,Diametr,kolichestvo,Visota,lambda_min_vata,delta_min_vata,ASWbul,accumulation,vremya,vremyaojidaniya,ASWatm = args

    # Расчет ГТУ
    Gas_turbine_res = Gas_turbine.calc(GTU_input)

    # Запись данных об электричестве
    electric.at["GTU", "N"] = Gas_turbine_res["N"]
    electric.at["GTU", "KPD"] = Gas_turbine_res["eff"]
    electric.at["DK", "N"] = Gas_turbine_res["Ndk"]

    # Запись данных о газе на выходе из ГТУ
    gas_streams = pd.read_excel("streams.xlsx", sheet_name="gas", index_col=0)
    gas_streams.at["GTU-KU", "T"] = Gas_turbine_res["T"]
    gas_streams.at["GTU-KU", "G"] = Gas_turbine_res["G"]
    gas_streams.at["GTU-KU", "P"] = 0.1
    gas_streams.at["GTU-KU", "H"] = gas0.p_t(
        gas_streams.at["GTU-KU", "P"], gas_streams.at["GTU-KU", "T"]
    )["h"]
    Gas_turbine_composition = pd.read_excel(
        "input.xlsx", sheet_name="Gas_composition0", index_col=0
    )
    
    #####################пкм
    #####Максимов Игорь#######
    if PKM_zaryad == True:
        
        steamVD_to_turbine = 0.25*water_streams0.at["PEVD-DROSVD", "G"]   ##!!!максимов
        from PKM import steam_transformer
        # Пар в паротрансформатор
        water_streams.loc["DROSVD-ST", "T":"H"] = water_streams.loc["PEVD-DROSVD", "T":"H"]
        water_streams.loc["DROSVD-ST", "G"] = water_streams.at["PEVD-DROSVD", "G"] - steamVD_to_turbine  ##!!!!
        # паротрансформатор
        ST = steam_transformer(
            stream11="DROSVD-ST",
            water=water,
            water_streams=water_streams,
            heaters=heaters,
            Pdr1=2,Pdr2=0.8,P2=2,dT=15,dTmin=5,Tdec=10)
        steam_trans = ST.calc()
        
        # Ввод в табл выходов из паротрансформатора
        water_streams.loc["ST-GPK", "T":"G"] = [steam_trans["T17"],steam_trans["P17"],steam_trans["H17"],steam_trans["G1"]]
        water_streams.loc["ST-PKM", "T":"G"] = [steam_trans["T24"],steam_trans["P2"],steam_trans["H24"],steam_trans["G2"]]
        heaters.at["Strans", "Qw"] = steam_trans["Q"]
        heaters.at["Strans_Qcool", "Qw"] = steam_trans['Qcool80']

        print(heaters.at["Strans_Qcool", "Qw"])
        # реформер
        from PKM import reformer

        ref = reformer(
            stream11="ST-PKM",
            water=water,
            gas_KU=gas_KU_PKM,
            Methane=Methane,
            waterMethane=waterMethane,
            water_streams=water_streams,
            heaters=heaters,
            Tref=700,
            Pref=2,
            T1gas=1968.58395330148,
            T2gas=800,
        )
        reform = ref.calc()
        # Газы реформера
        gas_streams.loc["AIR-REF", "T":"G"] = [15, 0.1, 414.38, reform["Gair"]]
        gas_streams.loc["CH4-REF", "T":"G"] = [15, 0.7, 881.50, reform["Gch4"]]
        gas_streams.loc["REF-SMESH", "T":"G"] = [800, 0.1, reform["H2gas"], reform["Ggas"]]
        gas_streams.loc["REF-SMESH", "N2":"Ar"] = list(reform["Gasfrac"].values())

        # Смешение
        gas_streams.loc["GTU-PEVD", "G"] = (gas_streams.at["REF-SMESH", "G"] + gas_streams.at["GTU-KU", "G"])
        gas_streams.loc["GTU-PEVD", "H"] = (gas_streams.at["REF-SMESH", "G"] * gas_streams.at["REF-SMESH", "H"]
                                            + gas_streams.at["GTU-KU", "G"] * gas_streams.at["GTU-KU", "H"]) / gas_streams.loc["GTU-PEVD", "G"]
        gas_streams.loc["GTU-PEVD", "P"] = 0.1
        
        from PKM import mixing_gases_molar
        
        mixing_gases_molar("GTU-KU", "REF-SMESH", "GTU-PEVD", gas_streams)
        for stream in gas_streams.index[4:10]:
            gas_streams.loc[stream, "N2":"Ar"] = gas_streams.loc["GTU-PEVD", "N2":"Ar"]
    else:
        gas_streams.loc["GTU-PEVD", "T":"Ar"] = gas_streams.loc["GTU-KU", "T":"Ar"]
        water_streams.loc["ST-GPK", "T":"G"] = [80,2,320,0]
        steamVD_to_turbine=0
# ----------------------------- конец Максимова Игоря


    # Параметры газа на входе в КУ
    gas_streams.loc["GTU-PEVD", "T":"P"] = gas_streams.loc["GTU-KU", "T":"P"]
    gas_streams.at["GTU-PEVD", "G"] = gas_streams.loc["GTU-KU", "G"]
    gas_streams.loc["GTU-PEVD", "N2":"Ar"] = Gas_turbine_composition.loc["Fraction", "N2":"Ar"]

#     Состав газов при частичной нагрузке
    fractiongas = list(gas_streams.loc["GTU-PEVD", "N2":"Ar"])

    gas1 = prop.Materials_prop(
        gasmix,
        fractiongas,
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


 # ----------------------------- Начало Максима Опарина

    ASW = ASW.Accum(water,water_streams,accumulation,ASWatm,
                         stream12 = "ASW-WOUT",
                         stream11 = "SP2-ASW",
                         stream_obratnoi_setevoi_vody = "SWIN-TURB",
                         stream_pryamoi_setevoi_vody = "SP2-WOUT",
                         T_nar_vozd = water_streams.at["AIR", "T"])
    
    ASW.set_construct(Diametr=Diametr,
                           kolichestvo = kolichestvo,
                           Visota = Visota,
                           lambda_min_vata = lambda_min_vata,
                           delta_min_vata = delta_min_vata)
    
    if ASWatm == False:
      # # если Зарядка   
        if ASWbul == 1:
            G_ASW_zarydka = ASW.zaryadka(vremya)['G']
            water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka
            water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka
                ## если разрядка
        if ASWbul == 2:
            ASW.zaryadka(vremya)
            ASW.jdat_n(vremyaojidaniya)
            G_ASW_razryadka = ASW.razryadka(vremya)['G']
            water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
            water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
            water_streams.at["SWOUT", "H"] = water.p_t(water_streams.at["SWOUT","P"],water_streams.at["SWOUT", "T"])['h']
            water_streams.at["SP2-WOUT", "H"] = (water_streams.at["SWOUT", "H"]*water_streams.at["SWOUT", "G"] -
                                                water_streams.at["ASW-WOUT", "H"]*G_ASW_razryadka)/(water_streams.at["SWOUT", "G"]-G_ASW_razryadka)
            water_streams.at["SP2-WOUT", "T"] = water.p_h(water_streams.at["SP2-WOUT", "P"],water_streams.at["SP2-WOUT", "H"])['T']
    else:

        # # если Зарядка   
        if ASWbul == 1:
            G_ASW_zarydka = ASW.zaryadka(vremya)['G']

            a = [[1, 1],
            [1*water_streams.at["SP2-WOUT", "H"], 1*water_streams.at["SWIN-TURB", "H"]]]
            b = [G_ASW_zarydka, G_ASW_zarydka*water.p_t(0.1,min(95,water_streams.at["SP2-WOUT", "T"]))['h']]
            x = LA.solve(a, b)
            G1 = x[0]
            G2 = x[1]
            water_streams.at["SWIN-ASWatm", "G"] = G2 = x[1]
            water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka - G2 
            water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka - G2
                ## если разрядка
        if ASWbul == 2:
            ASW.zaryadka(vremya)
            ASW.jdat_n(vremyaojidaniya)
            G_ASW_razryadka = ASW.razryadka(vremya)['G']
            water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
            water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
            water_streams.at["SWOUT", "H"] = water.p_t(water_streams.at["SWOUT","P"],water_streams.at["SWOUT", "T"])['h']
            water_streams.at["SP2-WOUT", "H"] = (water_streams.at["SWOUT", "H"]*water_streams.at["SWOUT", "G"] -
                                                water_streams.at["ASW-WOUT", "H"]*G_ASW_razryadka)/(water_streams.at["SWOUT", "G"]-G_ASW_razryadka)
            water_streams.at["SP2-WOUT", "T"] = water.p_h(water_streams.at["SP2-WOUT", "P"],water_streams.at["SP2-WOUT", "H"])['T']

# ----------------------------- Конец Максима Опарина

      # Инициализаця KU+TU, она здесь потому что нжно менять состав газа на входе в КУ
        
    
    KU_and_TU = KU_TU.ku_tu(
        gas0,
        gas1,
        water,
        gas_streams0,
        gas_streams,
        water_streams0,
        water_streams,
        heaters,
        electric,
        "PEVD-DROSVD",
        "PPND-DROSND",
        "DROSVD-TURBVD",
        "DROSND-TURBND",
        Calcmethod,
        KPD_SP,
        KPD_KN,
        KPD_to,
        KPD_PN,
        steamVD_fraction_to_turbine,
        steamVD_to_turbine,
    )
    start_time = time.time()
    # Расчет КУ и ТУ
    print(water_streams)
    KU_and_TU.calculate(
        Teplo,
        Calctolerance,
        Maxiterations_KU_TU,
        Maxiterations_cotel,
        Maxiterations_turbine,
    )

    print(f"fin КУ и ТУ:--- {round((time.time() - start_time), 1)} сек. ---")


