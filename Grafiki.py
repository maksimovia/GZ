import matplotlib.pyplot as plt


def Q_t_diagram(gas_streams, water_streams, heaters, xsize=15, ysize=10):
    fig = plt.figure(figsize=(xsize, ysize),dpi=500)
    Qg = [0,
          heaters.loc['PEVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'] +
          heaters.loc['EVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'] +
          heaters.loc['EVD', 'Qg']+heaters.loc['PPND', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg']+heaters.loc['EVD',
                                                                         'Qg']+heaters.loc['PPND', 'Qg']+heaters.loc['IND', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg']+heaters.loc['EVD', 'Qg'] +
          heaters.loc['PPND', 'Qg'] +
          heaters.loc['IND', 'Qg']+heaters.loc['GPK', 'Qg']
          ]
    Qw = [0,
          heaters.loc['PEVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'] +
          heaters.loc['EVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'] +
          heaters.loc['EVD', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg'] +
          heaters.loc['EVD', 'Qg']+heaters.loc['PPND', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg']+heaters.loc['EVD',
                                                                         'Qg']+heaters.loc['PPND', 'Qg']+heaters.loc['IND', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg']+heaters.loc['EVD',
                                                                         'Qg']+heaters.loc['PPND', 'Qg']+heaters.loc['IND', 'Qg'],
          heaters.loc['PEVD', 'Qg']+heaters.loc['IVD', 'Qg']+heaters.loc['EVD', 'Qg'] +
          heaters.loc['PPND', 'Qg'] +
          heaters.loc['IND', 'Qg']+heaters.loc['GPK', 'Qg']
          ]
    Tg = gas_streams.loc['GTU-PEVD':'GPK-out', 'T']
    Tw = [water_streams.loc['PEVD-DROSVD', 'T'],
          water_streams.loc['IVD-PEVD', 'T'],
          water_streams.loc['IVD-PEVD', 'T'],
          water_streams.loc['EVD-IVD', 'T'],
          water_streams.loc['PEN-EVD', 'T'],
          water_streams.loc['PPND-DROSND', 'T'],
          water_streams.loc['IND-PPND', 'T'],
          water_streams.loc['IND-PPND', 'T'],
          water_streams.loc['GPK-REC', 'T'],
          water_streams.loc['REC-GPK', 'T'],
          ]
    plt.plot(Qg, Tg, color='red')
    plt.plot(Qw[0:2], Tw[0:2], color='blue')
    plt.plot(Qw[1:3], Tw[1:3], color='blue')
    plt.plot(Qw[3:5], Tw[3:5], color='blue')
    plt.plot(Qw[5:7], Tw[5:7], color='blue')
    plt.plot(Qw[6:8], Tw[6:8], color='blue')
    plt.plot(Qw[8:10], Tw[8:10], color='blue')
    plt.xlabel('Q')
    plt.ylabel('T')
    plt.legend(['gas', 'water'])
    # plt.show()
    return fig

    
def H_S_diagram(water, water_streams, xsize=15, ysize=10):
    water_streams.at['PEVD-DROSVD', 'S'] = water.p_h(water_streams.at['PEVD-DROSVD', 'P'],water_streams.at['PEVD-DROSVD', 'H'])['s']
    water_streams.at['PPND-DROSND', 'S'] = water.p_h(water_streams.at['PPND-DROSND', 'P'],water_streams.at['PPND-DROSND', 'H'])['s']

    fig=plt.figure(figsize=(xsize,ysize),dpi=500)
    Hvd = [water_streams.at['PEVD-DROSVD', 'H'],
           water_streams.at['DROSVD-TURBVD', 'H'],
           water_streams.at['ENDOFVD', 'H'],
           water_streams.at['SMESHEND', 'H'],
    ]
    Svd = [water_streams.at['PEVD-DROSVD', 'S'],
           water_streams.at['DROSVD-TURBVD', 'S'],
           water_streams.at['ENDOFVD', 'S'],
           water_streams.at['SMESHEND', 'S'],
    ]
    Hsm = [water_streams.at['PPND-DROSND', 'H'],
           water_streams.at['DROSND-TURBND', 'H'],
           water_streams.at['SMESHEND', 'H'],
    ]
    Ssm = [water_streams.at['PPND-DROSND', 'S'],
           water_streams.at['DROSND-TURBND', 'S'],
           water_streams.at['SMESHEND', 'S'],
    ]
    Hnd = [water_streams.at['SMESHEND', 'H'],
           water_streams.at['DOOTB2', 'H'],
           water_streams.at['DOOTB1', 'H'],
           water_streams.at['INCND', 'H'],
           water_streams.at['INKOND', 'H'],
    ]
    Snd = [water_streams.at['SMESHEND', 'S'],
           water_streams.at['DOOTB2', 'S'],
           water_streams.at['DOOTB1', 'S'],
           water_streams.at['INCND', 'S'],
           water_streams.at['INKOND', 'S'],
    ]

    #Давления
    # stream = 'PEVD-DROSVD'
    # H1=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    # S1=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    # stream = 'DROSVD-TURBVD'
    # H2=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    # S2=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'ENDOFVD'
    H3=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S3=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'PEVD-DROSVD'
    H4=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S4=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'DROSVD-TURBVD'
    H5=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S5=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'PPND-DROSND'
    H6=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S6=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'DROSND-TURBND'
    H7=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S7=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]

    stream = 'DOOTB2'
    H8=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S8=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'DOOTB1'
    H9=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S9=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'INCND'
    H10=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S10=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]
    stream = 'INKOND'
    H11=[water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']-0.05)['h'],water_streams.at[stream,'H'],water.p_s(water_streams.at[stream, 'P'],water_streams.at[stream, 'S']+0.05)['h']]
    S11=[water_streams.at[stream, 'S']-0.05,water_streams.at[stream, 'S'],water_streams.at[stream, 'S']+0.05]


    plt.plot(Svd,Hvd)
    plt.plot(Ssm,Hsm)
    plt.plot(Snd,Hnd)
    plt.plot(S3,H3,S4,H4,S5,H5,S6,H6,S7,H7,S8,H8,S9,H9,S10,H10,S11,H11, color = "gray", alpha=0.3)
    plt.xlabel('S')
    plt.ylabel('H')
    # plt.show()
    return fig