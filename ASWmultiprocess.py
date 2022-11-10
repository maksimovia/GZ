def ParallelCompute(args):
    Diametr,kolichestvo, Visota,lambda_min_vata,delta_min_vata,ASWbul,nagr,vremya,vremyaojidaniya, t_air_list, PKM_zaryad,Сalculate_minimum,KPD_PN,KPD_KN,KPD_to,KPD_SP,Calcmethod,Calctolerance,Maxiterations_KU_TU,Maxiterations_cotel,Maxiterations_turbine,Сalculate_minimum,Teplo,steamVD_fraction_to_turbine,PKM_zaryad,ASWatm= args
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
    import calculate_all

    # таблица номинального режима
    gas_streams0 = pd.read_excel("streams0.xlsx", sheet_name="gas", index_col=0)
    water_streams0 = pd.read_excel("streams0.xlsx", sheet_name="water", index_col=0)
    GTU_ISO = pd.read_excel("input.xlsx", sheet_name="ISO", index_col=0)
    GTU_input = pd.read_excel("input.xlsx", sheet_name="GTU_input", index_col=0)
    GTU_input.at["n", 1] = nagr
    GTU_input.at["tair", 1] = t_air_list
    
    # рабочая таблица (=номинал в 1 итерации)
    gas_streams = pd.read_excel("streams.xlsx", sheet_name="gas", index_col=0)
    water_streams = pd.read_excel("streams.xlsx", sheet_name="water", index_col=0)
    # рабочая таблица показателей блоков
    heaters = pd.read_excel("blocks.xlsx", sheet_name="heaters", index_col=0)
    electric = pd.read_excel("blocks.xlsx", sheet_name="electric", index_col=0)
    accumulation = pd.read_excel("blocks.xlsx", sheet_name="accumulation", index_col=0)
    ############################################################
    # Теплосеть и перекидка температуры воздуха
    gas_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    water_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    Tnv = gas_streams.at["AIR", "T"]
    water_streams.at["SWINB", "T"] = SP.Tset(Tnv)[1]
    water_streams.at["SWOUT", "T"] = SP.Tset(Tnv)[0]
    water_streams.at["SWIN-TURB", "T"] = water_streams.at["SWINB", "T"]
    water_streams.at["SP2-WOUT", "T"] = water_streams.at["SWOUT", "T"]
    ############################################################
    # Состав газов в номинале в ГТУ
    gasmix = "Nitrogen*Oxygen*CO2*Water*Argon"
    # Считывание рефпропа
    RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
    fractiongas0 = list(gas_streams0.loc["GTU-PEVD", "N2":"Ar"])
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

    #####################
    ######Максимов#######
    fractionwaterMethane = (0.833372660622383, 0.166627339377617, 0, 0, 0)
    waterMethanemix = "Water*METHANE"

    waterMethane = prop.Materials_prop(
        waterMethanemix,
        fractionwaterMethane,
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
    Methane = prop.Materials_prop(
        "METHANE",
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
    gas_KU_PKM = prop.Materials_prop(
        gasmix,
        (
        0.710320591016015,
        0.00996710270335893,
        0.090538556815177,
        0.180531273012258,
        0.00864247645319178,
    ),
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

    ############################################################
    # Задание ГТУ
    Gas_turbine = GTU.gtu(GTU_ISO, "GTU-KU")

    args = [gas_streams0,
            water_streams0,
            GTU_ISO,
            GTU_input,
            gas_streams,
            water_streams,
            heaters,
            electric,
            Tnv,
            KPD_PN,
            KPD_KN,
            KPD_to,
            KPD_SP,
            Calcmethod,
           Calctolerance,
           Maxiterations_KU_TU,
           Maxiterations_cotel,
           Maxiterations_turbine,
           Сalculate_minimum,
           Teplo,
           steamVD_fraction_to_turbine,
           gasmix,
           RP,
           fractiongas0,
           gas0,
           water,
           fractionwaterMethane,
           waterMethanemix,
           waterMethane,
           Methane,
           gas_KU_PKM,
           PKM_zaryad,
           Temperatures,
           Pressure,
           Enthalpies,
           Gas_turbine,
           Diametr,
           kolichestvo,
           Visota,
           lambda_min_vata,
           delta_min_vata,
           ASWbul,
           accumulation,
           vremya,
           vremyaojidaniya,
           ASWatm,
           ]

    calculate_all.calculate_CCGT(args)
    # steamVD_to_turbine = 0
#     GTU_input = pd.read_excel("input.xlsx", sheet_name="GTU_input", index_col=0)
#     GTU_input.at["n", 1] = nagr
#     GTU_input.at["tair", 1] = t_air_list

    Max_iterations_minimum = 10
    if Teplo == 1:
        n_GTU = GTU_input.at["n", 1]
        Delt_Gcnd = (water_streams.at["INKOND", "G"] - 4.44) / 4.44
        Delt_Nturb = (electric.at["Turbine", "Ni"] - 17.6) / 17.6
        Delt_Gvd = (water_streams.at["PEVD-DROSVD", "G"] / water_streams0.at["PEVD-DROSVD", "G"]- 0.25)/ 0.25
        Delt_Gnd = (water_streams.at["PPND-DROSND", "G"] / water_streams0.at["PPND-DROSND", "G"]- 0.5) / 0.5
        Delta_min = min(Delt_Gcnd, Delt_Nturb, Delt_Gvd, Delt_Gnd)
        print(Delta_min,'check')
        if n_GTU == 1 and Delta_min < 0:
            print("Мощность ГТУ 100% и расход пара все еще слишком мал")

    if Сalculate_minimum == True:
        n_GTU = GTU_input.at["n", 1]
        start_time = time.time()
        for i in range(Max_iterations_minimum):
            print(f"Началась {i+1} итерация расчета ПГУ")
            if i < round(Max_iterations_minimum/2,1):
                New_iterations_KU_TU, New_iterations_cotel, New_iterations_turbine = (
                    2,
                    2,
                    15,
                )
                args = [gas_streams0,
                        water_streams0,
                        GTU_ISO,
                        GTU_input,
                        gas_streams,
                        water_streams,
                        heaters,
                        electric,
                        Tnv,
                        KPD_PN,
                        KPD_KN,
                        KPD_to,
                        KPD_SP,
                        Calcmethod,
                       Calctolerance,
                       2,#Maxiterations_turbine
                       2,#Maxiterations_turbine
                       15,#Maxiterations_turbine
                       Сalculate_minimum,
                       Teplo,
                       steamVD_fraction_to_turbine,
                       gasmix,
                       RP,
                       fractiongas0,
                       gas0,
                       water,
                       fractionwaterMethane,
                       waterMethanemix,
                       waterMethane,
                       Methane,
                       gas_KU_PKM,
                       PKM_zaryad,
                       Temperatures,
                       Pressure,
                       Enthalpies,
                       Gas_turbine,
                       Diametr,
                       kolichestvo,
                       Visota,
                       lambda_min_vata,
                       delta_min_vata,
                       ASWbul,
                       accumulation,
                       vremya,
                       vremyaojidaniya,
                       ASWatm,
                       ]
            else:
                New_iterations_KU_TU, New_iterations_cotel, New_iterations_turbine = (
                    Maxiterations_KU_TU,
                    Maxiterations_cotel,
                    Maxiterations_turbine,
                )


            Delt_Gcnd = 100
            Delt_Nturb = 100
            if Teplo == 1:
                Delt_Gcnd = (water_streams.at["INKOND", "G"] - 4.44) / 4.44
                Delt_Nturb = (electric.at["Turbine", "Ni"] - 17.6) / 17.6
            Delt_Gvd = (
                water_streams.at["PEVD-DROSVD", "G"] / water_streams0.at["PEVD-DROSVD", "G"]
                - 0.25
            ) / 0.25
            Delt_Gnd = (
                water_streams.at["PPND-DROSND", "G"] / water_streams0.at["PPND-DROSND", "G"]
                - 0.5
            ) / 0.5
            Delta_min = min(Delt_Gcnd, Delt_Nturb, Delt_Gvd, Delt_Gnd)
            if n_GTU == 1 and Delta_min < 0:
                print("Мощность ГТУ 100% и расход пара все еще слишком мал")
            n_GTU = n_GTU - Delta_min / 10
            GTU_input.at["n", 1] = n_GTU

            calculate_all.calculate_CCGT(args)
            print(f"Отклонение от ограничения минимальное равно {Delta_min}")
            if abs(Delta_min) < Calctolerance:
                calculate_all.calculate_CCGT(args)
                print(f"Отклонение от ограничения минимальное равно {Delta_min}")
                print(f"Относительная мощность ГТУ равна {n_GTU}")
                print(
                    f"fin минимальная мощность ПГУ:--- {round((time.time() - start_time), 1)} сек. ---"
                )
                break
            if i == Max_iterations_minimum - 1:
                print("Достигнуто максимальное количество итераций минимального расхода в ПГУ", i+1)

    result = {
    "T_air":round(t_air_list,2),
    "n_GTU":round(GTU_input.at["n", 1],2),
    "GTU": round(electric.at["GTU", "N"], 4),
    "GTU_KPD": round(electric.at["GTU", "KPD"], 4),
    "Turbine": round(electric.at["Turbine", "Ni"], 4),
    "KN": round(electric.at["KN", "Ni"], 4),
    "DK": round(electric.at["DK", "N"], 4),
    "PEN": round(electric.at["PEN", "Ni"], 4),
    "Turbine_Qt":round(heaters.at["SP2", "Qw"]+heaters.at["SP1", "Qw"]+heaters.at["OD", "Qw"], 4),
    "ASW_Qt":round(accumulation.at["ASW", "Qw"]/(vremya*3600), 4),
    "ASW_bull":ASWbul,
    "Delta_P_Diafragma":round(water_streams.at["INCND", "P"]-water_streams.at["DOOTB1", "P"],4),
    "INKOND": round(water_streams.at["INKOND", "G"],4),
    
#     "Delta_P_Diafragma":round(Delta_P_Diafragma, 4)
    }
    return result