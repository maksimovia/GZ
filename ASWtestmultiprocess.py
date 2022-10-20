
# Diametr=1
# kolichestvo = 2
# Visota = 1
# ASWbul = 1
# lambda_min_vata = 0.045
# delta_min_vata = 0.01
# nagr = 1
# vremya=1

def test(args):

        # Импорт библиотек
    import os
    import time
    import ASW
    import GTU
    import KU_TU
    import mat_properties as prop
    import numpy as n
    import pandas as pd
    import SP
    from scipy.optimize import root

    Diametr,kolichestvo, Visota,lambda_min_vata,delta_min_vata,ASWbul,nagr,vremya,vremyaojidaniya, t_air_list = args
    # таблица номинального режима
    gas_streams0 = pd.read_excel("streams0.xlsx", sheet_name="gas", index_col=0)
    water_streams0 = pd.read_excel("streams0.xlsx", sheet_name="water", index_col=0)
    GTU_ISO = pd.read_excel("input.xlsx", sheet_name="ISO", index_col=0)
    GTU_input = pd.read_excel("input.xlsx", sheet_name="GTU_input", index_col=0)
    # рабочая таблица (=номинал в 1 итерации)
    gas_streams = pd.read_excel("streams.xlsx", sheet_name="gas", index_col=0)
    water_streams = pd.read_excel("streams.xlsx", sheet_name="water", index_col=0)
    # рабочая таблица показателей блоков
    heaters = pd.read_excel("blocks.xlsx", sheet_name="heaters", index_col=0)
    electric = pd.read_excel("blocks.xlsx", sheet_name="electric", index_col=0)
    accumulation = pd.read_excel("blocks.xlsx", sheet_name="accumulation", index_col=0)
    # Состав газов в номинале
    gasmix = "Nitrogen*Oxygen*CO2*Water*Argon"
    # Считывание рефпропа
    RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
    fractiongas0 = (
        gas_streams0.at["GTU-PEVD", "N2"],
        gas_streams0.at["GTU-PEVD", "O2"],
        gas_streams0.at["GTU-PEVD", "CO2"],
        gas_streams0.at["GTU-PEVD", "H2O"],
        gas_streams0.at["GTU-PEVD", "Ar"],
    )

    gas0 = prop.Materials_prop(
        gasmix,
        fractiongas0,
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
    water = prop.Materials_prop(
        "water",
        [1.0, 0, 0, 0, 0],
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

    # Задание энтальпий газа в номинальном режиме
    Temperatures = gas_streams0.loc["GTU-KU":"GPK-out", "T"]
    Pressure = gas_streams0.loc["GTU-KU", "P"]
    Enthalpies = list(map(lambda x: gas0.p_t(Pressure, x)["h"], Temperatures))
    gas_streams0.loc["GTU-KU":"GPK-out", "H"] = Enthalpies



    # Состав газов при частичной нагрузке
    fractiongas = (
        gas_streams.at["GTU-PEVD", "N2"],
        gas_streams.at["GTU-PEVD", "O2"],
        gas_streams.at["GTU-PEVD", "CO2"],
        gas_streams.at["GTU-PEVD", "H2O"],
        gas_streams.at["GTU-PEVD", "Ar"],
    )
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

    # Основные константы
    Calcmethod = "hybr"
    KPD_PN = 0.8074
    KPD_KN = 0.75
    KPD_to = 0.99
    KPD_SP = 0.99

    Calctolerance = 10**-2
    Teplo = 1
    Maxiterations_KU_TU = 20
    Maxiterations_cotel = 4
    Maxiterations_turbine = 30
    steamVD_fraction_to_turbine=1
    steamVD_to_turbine=0

#     n = name.n

    # Задаем нагрузку
    GTU_input.at["n", 1] = nagr
    GTU_input.at["tair", 1] = t_air_list
#     vremya = 4
     # # Выбор происходит ли зарядка : 1 - заряжается, 2 - разряжается, любое другое число - не участвует в расчете
#     ASWbul = name.ASWbul
    #----------------------------------------------------------
    # print(GTU_input)
    ############################################################
    # Теплосеть
    gas_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    water_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    Tnv = gas_streams.at["AIR", "T"]
    water_streams.at["SWIN-TURB", "T"] = SP.Tset(Tnv)[1]
    water_streams.at["SP2-WOUT", "T"] = SP.Tset(Tnv)[0]

    water_streams.at["SWOUT", "T"] = SP.Tset(Tnv)[0]
    water_streams.at["SWIN", "T"] = SP.Tset(Tnv)[1]
    # print(water_streams)
    ############################################################
    #--------------------------------Расчет ГТУ
    Gas_turbine = GTU.gtu(GTU_ISO, GTU_input, "GTU-KU")
    Gas_turbine_res = Gas_turbine.calc()
    electric.at["GTU", "N"] = Gas_turbine_res["N"]
    electric.at["GTU", "KPD"] = Gas_turbine_res["eff"]
    electric.at["DK", "N"] = Gas_turbine_res["Ndk"]
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

    # Параметры газа на выходе в КУ
    gas_streams.loc["GTU-PEVD", "T":"P"] = gas_streams.loc["GTU-KU", "T":"P"]
    gas_streams.at["GTU-PEVD", "G"] = gas_streams.loc["GTU-KU", "G"]
    gas_streams.loc["GTU-PEVD", "N2":"Ar"] = Gas_turbine_composition.loc[
        "Fraction", "N2":"Ar"
    ]
    #--------------------------------------------------------------------------
    # Class KU+TU
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

    ASW = ASW.Accum(water,water_streams,accumulation,
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

    start_time = time.time()



        # # если Зарядка 
    if ASWbul == 1:
        G_ASW_zarydka = ASW.zaryadka(vremya)['G']
        water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka
        water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka
            ## если разрядка
    if ASWbul == 2:
        ASW.zaryadka(vremya)
        ASW.jdat(vremyaojidaniya)
        G_ASW_razryadka = ASW.razryadka(vremya)['G']
        water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
        water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
        water_streams.at["SWOUT", "H"] = water.p_t(water_streams.at["SWOUT","P"],water_streams.at["SWOUT", "T"])['h']
        water_streams.at["SP2-WOUT", "H"] = (water_streams.at["SWOUT", "H"]*water_streams.at["SWOUT", "G"] -
                                            water_streams.at["ASW-WOUT", "H"]*G_ASW_razryadka)/(water_streams.at["SWOUT", "G"]-G_ASW_razryadka)
        water_streams.at["SP2-WOUT", "T"] = water.p_h(water_streams.at["SP2-WOUT", "P"],water_streams.at["SP2-WOUT", "H"])['T']

    KU_and_TU.calculate(
        Teplo,
        Calctolerance,
        Maxiterations_KU_TU,
        Maxiterations_cotel,
        Maxiterations_turbine,
    )

    print(
        "Степень сухости пара в ЭВД: ",
        water.p_h(water_streams.at["EVD-IVD", "P"], water_streams.at["EVD-IVD", "H"])["Q"],
    )
    print(
        "Степень сухости пара в ГПК: ",
        water.p_h(water_streams.at["GPK-IND", "P"], water_streams.at["GPK-IND", "H"])["Q"],
    )

    print(f"fin КУ и ТУ:--- {round((time.time() - start_time), 1)} сек. ---")

    return {'GTU':electric.at("GTU","N"), "Turbine":electric.at("Turbine","Ni"),"KN":electric.at("KN","Ni"),"DK":electric.at("DK","N"),"PEN":electric.at("PEN","Ni")}
    

