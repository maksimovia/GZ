import os
import time

import GTU
import KU_TU
import mat_properties as prop
import numpy as n
import pandas as pd
import SP
from scipy.optimize import root


def calculate_CCGT_PKM(arguments_all):
    Iterations_KU_TU, Iterations_cotel, Iterations_turbine, gas_streams0, water_streams0, GTU_ISO, GTU_input, gas_streams, water_streams, heaters, electric, Gas_turbine, gas0,    water, PKM_zaryad, PKM_razryad,    syngas_streams, Calcmethod,    Calctolerance, KPD_PN, KPD_KN, KPD_to, KPD_SP,    steamVD_fraction_to_turbine, accumulation, time_ac = arguments_all

    gasmix = "Nitrogen*Oxygen*CO2*Water*Argon"
    RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")

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

    #####################Максимов#####################
    from PKM import accum
    if PKM_zaryad:
        Accumulator = accum.zaryad(time_ac, accumulation, gas_streams,
                                   syngas_streams, water_streams, water_streams0, heaters, electric)
        steamVD_to_turbine = Accumulator['steamVD_to_turbine']
        Teplo = Accumulator['Teplo']
        # print("Zaryad")

    elif PKM_razryad:
        Accumulator = accum.razryad(time_ac, accumulation, gas_streams,
                                    syngas_streams, water_streams, water_streams0, heaters, electric)
        Teplo = Accumulator['Teplo']
        steamVD_to_turbine = Accumulator['steamVD_to_turbine']

    else:
        gas_streams.loc["GTU-PEVD",
                        "T":"Ar"] = gas_streams.loc["GTU-KU", "T":"Ar"]
        water_streams.loc["ST-GPK", "T":"G"] = [0, 0, 0, 0]
        steamVD_to_turbine = water_streams.at["PEVD-DROSVD", "G"]
        Teplo = 1
    ##################################################

    # Состав газов при частичной нагрузке
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

    # Инициализаця KU+TU, она здесь потому что нужно менять состав газа на входе в КУ

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

    # print('Teplo',Teplo)
    # Расчет КУ и ТУ
    KU_and_TU.calculate(
        Teplo,
        Calctolerance,
        Iterations_KU_TU,
        Iterations_cotel,
        Iterations_turbine,
    )
#     print(gas_streams)
    print(f"fin КУ и ТУ:--- {round((time.time() - start_time), 1)} сек. ---")
    return gas_streams

def Calculate_CCGT_PKM_iter(arguments_all_it,Iter_pkm,pkm_pgu_tol):
    Maxiterations_KU_TU,    Maxiterations_cotel = arguments_all_it[0], arguments_all_it[1]
    start_time = time.time()
    water_streams0=arguments_all_it[4]
    water_streams=arguments_all_it[8]
    print('Gst',water_streams.at["PEVD-DROSVD", "G"],'Gst',round(water_streams.at["SMESH-GPK", "G"]))

    Gst = [max([water_streams0.at["DROSVD-ST", "G"]],round(water_streams.at["PEVD-DROSVD", "G"], 2))]
    Ggpk = [max([water_streams0.at["SMESH-GPK", "G"]],round(water_streams.at["SMESH-GPK", "G"], 2))]

    for i in range(Iter_pkm):
        if i < 6:
            Maxiterations_KU_TU_new = 2
            Maxiterations_cotel_new = 2
            # print(Maxiterations_KU_TU_new,Maxiterations_cotel_new)

        else:
            Maxiterations_KU_TU_new = Maxiterations_KU_TU
            Maxiterations_cotel_new = Maxiterations_cotel

        # from calculate_CCGT_PKM import calculate_CCGT_PKM

        arguments_all = arguments_all_it.copy()
        arguments_all[0], arguments_all[1] = [
            Maxiterations_KU_TU_new,
            Maxiterations_cotel_new,
        ]

        gas_streams = calculate_CCGT_PKM(arguments_all)
        ####################################
        print(
            f"Время {i+1} итерации расчета КУ+ТУ с ПКМ: --- {round((time.time() - start_time), 1)} сек. ---"
        )
        Gst.append(round(water_streams.at["PEVD-DROSVD", "G"], 2))
        Ggpk.append(round(water_streams.at["SMESH-GPK", "G"], 2))



        Err1 = abs((Gst[i] - Gst[i - 1]) / (Gst[i]) * 100)
        Err3 = abs((Ggpk[i] - Ggpk[i - 1]) / (Ggpk[i]) * 100)

        if i == Iter_pkm - 1:
            print("Достигнуто максимальное количество итераций КУ+ПТУ+ПКМ")
            print("Gst", Gst)
            print("Ggpk", Ggpk)

        if Err1 < pkm_pgu_tol and Err3 < pkm_pgu_tol:

            print(
                f"Расчет КУ+ПТУ+ПКМ закончен:--- %s сек. ---{round((time.time() - start_time), 1)}"
            )
            print("Gst", Gst)
            print("Ggpk", Ggpk)
            # print(Err1, Err3)

            break
    return gas_streams
