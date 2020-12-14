from math import sqrt, log10, log2, pow, erfc, e
from scipy.special import erfcinv
from matplotlib import pyplot as plt

import argparse

#Parâmetros de referência: Link Budget de Uplink
uplink_pot_tx = 26 #dBm
uplink_Gtx = 1.2
uplink_tx_loss = 0
uplink_SNR_req = 0
uplink_sens_req_rx = -101.5 #dBm
uplink_Grx = 19 #dBi
uplink_rx_loss = 3 #dB
uplink_multipath = 3 #dB
uplink_shadow_margin = 4 #dB

#Service mix
#0 -> average bit rate
#1 -> Service mix
service_mix = {
    "VoLTE": [22000, 0.22],
    "Video Call": [384000, 0.08],
    "Video Stream": [2000000, 0.28],
    "Music Stream": [196000, 0.2],
    "Web Browse": [2000000, 0.1],
    "File Share": [2000000, 0.08],
    "E-mail": [1000000, 0.04]}

def start_mcs():
    # 0 -> bit rate
    # 1 -> code rate

    mcs = {
        "QPSK 1/3":[2, 0.333333, 2342010, 14.0051, -0.577897],
        "16-QAM 1/2":[4, 0.5, 47613.1, 0.0926275, -0.29583],
        "64-QAM 3/4":[6, 0.75, 26405.8, 0.0220186, -0.24491]}

    return mcs

def count_hex(area, radius):

    hex_area_m = 1.5*sqrt(3)*pow(radius, 2)
    hex_area_km = hex_area_m/1000000
    total_hex = area//hex_area_km

    return hex_area_km, total_hex

def link_budget(pot_tx, pot_rx, Ms, Gt, Pt, decrease):
    Lmax = pot_tx - pot_rx + Gt - Pt - Ms - decrease
    return Lmax

def max_radius(max_loss, freq):
    alfa = 4
    beta = 10.2
    gama = 2.36

    log_R_max = (max_loss - beta - 10*gama*log10(freq))/(10*alfa)
    R_max = pow(10, log_R_max)
    #print(modulation, "log:", log_R_max, "normal:", R_max)

    return R_max

def loss_model(dist, freq):
    alfa = 4
    beta = 10.2
    gama = 2.36

    loss = 10*alfa*log10(dist) + beta + 10*gama*log10(freq)

    return loss

def outdoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs):
    #prob_cobertura_borda = float()
    Ms = 4*n/sigma - 3
    Q = 1 - (erfc(Ms/(sigma*sqrt(2))))/2
    #prob_cobertura_borda = 1 - Q

    sigma_sir = sqrt(2*sigma*(1-ro))
    Q_inv = sqrt(2)*erfcinv(2*Q)
    M_in = -1*Q_inv*sigma_sir

    m_in = pow(10, M_in/10)
    D_in = 10*log10(m_in*gama*eta)

    setor_gain = 4.77
    rate_rb = float()
    rate_total = float()
    sinr_mod = [1, 6, 14]

    for n, modulation in enumerate(mcs):
        sens = sinr_mod[n] + figura + 10*log10(bw*5*180000) - 174 + D_in
        loss = link_budget(uplink_pot_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss,
                            0)
        mcs[modulation].append(max_radius(loss, 2.6))

    return mcs, m_in

def indoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs):
    sigma = sqrt(64 + pow(sigma, 2))

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2))))/2
    prob_cobertura_borda = 1 - Q

    sigma_sir = sqrt(2*sigma*(1-ro))
    Q_inv = sqrt(2)*erfcinv(2*prob_cobertura_borda)
    M_in = -1*Q_inv*sigma_sir

    m_in = pow(10, M_in/10)
    D_in = 10*log10(m_in*gama*eta)

    rate_rb = float()
    rate_total = float()
    temp_rate = float()
    temp_sinr = float()
    sinr_mod = [1, 6, 14]

    for n, modulation in enumerate(mcs):
        sens = sinr_mod[n] + figura + 10*log10(bw*5*180000) - 174 + D_in
        loss = link_budget(uplink_pot_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss,
                            20)
        mcs[modulation].append(max_radius(loss, 2.6))

    return mcs

def incar_radius(sigma, n, ro, gama, eta, bw, figura, mcs):
    sigma = sqrt(9 + pow(sigma, 2))

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2))))/2
    prob_cobertura_borda = 1 - Q

    sigma_sir = sqrt(2*sigma*(1-ro))
    Q_inv = sqrt(2)*erfcinv(2*prob_cobertura_borda)
    M_in = -1*Q_inv*sigma_sir

    m_in = pow(10, M_in/10)
    D_in = 10*log10(m_in*gama*eta)

    rate_rb = float()
    rate_total = float()
    temp_rate = float()
    temp_sinr = float()
    sinr_mod = [1, 6, 14]

    for n, modulation in enumerate(mcs):
        sens = sinr_mod[n] + figura + 10*log10(bw*5*180000) - 174 + D_in
        loss = link_budget(uplink_pot_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss,
                            10)
        mcs[modulation].append(max_radius(loss, 2.6))

    return mcs

def active_users(users, area, prob_outdoor, prob_indoor, prob_incar, mcs, cell_area):
    active_users_total = users*0.15*0.1 #Penetration Ratio de 15% e Usage Ratio de 10%
    user_density = active_users_total/area

    active_users_per_cell = user_density*cell_area

    #Raios iniciais para cada modulação
    radius_QPSK = prob_outdoor*mcs["QPSK 1/3"][5] + prob_indoor*mcs["QPSK 1/3"][6] + prob_incar*mcs["QPSK 1/3"][7]
    radius_16QAM = prob_outdoor*mcs["16-QAM 1/2"][5] + prob_indoor*mcs["16-QAM 1/2"][6] + prob_incar*mcs["16-QAM 1/2"][7]
    radius_64QAM = prob_outdoor*mcs["64-QAM 3/4"][5] + prob_indoor*mcs["64-QAM 3/4"][6] + prob_incar*mcs["64-QAM 3/4"][7]

    #Usuários ativos por modulação por célula
    N_users_QPSK_numerador = pow(0.95*radius_QPSK,2) - pow(radius_16QAM,2)
    N_users_16QAM_numerador = pow(radius_16QAM,2) - pow(radius_64QAM,2)
    N_users_64QAM_numerador = pow(radius_64QAM,2)

    N_users_QPSK = (N_users_QPSK_numerador / pow(radius_QPSK,2)) * active_users_per_cell
    N_users_16QAM = (N_users_16QAM_numerador / pow(radius_QPSK,2)) * active_users_per_cell
    N_users_64QAM = (N_users_64QAM_numerador / pow(radius_QPSK,2)) * active_users_per_cell

    active_users = {
        "QPSK 1/3": [N_users_QPSK],
        "16-QAM 1/2": [N_users_16QAM],
        "64-QAM 3/4": [N_users_64QAM]}

    print("Average users using QPSK:", active_users["QPSK 1/3"][0])
    print("Average users using 16QAM:", active_users["16-QAM 1/2"][0])
    print("Average users using 64QAM:", active_users["64-QAM 3/4"][0])
    print("Average active users per cell:", N_users_QPSK + N_users_16QAM + N_users_64QAM)
    return active_users

def required_throughput(active_users, service_mix):
    #per modulation per cell
    for modulation in active_users:
        total_throughput = 0
        for service in service_mix:
            service_bitrate = active_users[modulation][0]*service_mix[service][0]*service_mix[service][1]
            total_throughput = total_throughput + service_bitrate
        print("Required Total Throughput for (", modulation, "): ", round(total_throughput,2))
        active_users[modulation].append(total_throughput)

    return active_users

def check_sinr(active_users, mcs, m_in, gama):
    min_sinr = int()
    sinr_list = [1, 6, 14]

    for n, modulation in enumerate(active_users):
        for sinr in reversed(range(sinr_list[n], 34)):
            rate_rb = mcs[modulation][2]/(mcs[modulation][3] + pow(e, mcs[modulation][4]*sinr))
            rate_total = rate_rb*(5*bw)

            if rate_total > active_users[modulation][1]:
                min_sinr = sinr
            elif rate_total == active_users[modulation][1]:
                min_sinr = sinr
                break
            else:
                break
        
        mean_sir = m_in*gama*min_sinr

        for k in range(1, 14):
            mean_sir = mean_sir + 10*log10(pow(k, 2))

            if mean_sir <= min_sinr:
                mean_sir = mean_sir + 4.77

                if mean_sir <= min_sinr:
                    mean_sir = mean_sir + 3

                    if mean_sir <= min_sinr:
                        mean_sir = mean_sir - 7.77 - 10*log10(pow(k, 2))
                    else:
                        print("Mean SIR of", mean_sir, "dB using six sectors per cell and reuse factor", k)
                        break
                else:    
                    print("Mean SIR of", mean_sir, "dB using three sectors per cell and reuse factor", k)
                    break
            else:
                print("Mean SIR of", mean_sir, "dB using no sectorization and reuse factor", k)
                break

def var_to_rad(sigma, n, ro, gama, eta, bw, figura, mcs):
    mcs, m_in = outdoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs)
    mcs = indoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs)
    mcs = incar_radius(sigma, n, ro, gama, eta, bw, figura, mcs)

    radius = 0.3*mcs["QPSK 1/3"][5] + 0.4*mcs["QPSK 1/3"][6] + 0.3*mcs["QPSK 1/3"][7]

    return radius

def print_graph(rad_list, type_var):
    variables = list()
    
    for x in range(1, 100):
        y = x/100

        if type_var == "eta":
            y = y + 1

        variables.append(y)
    
    fig_title = "Variação do raio de acordo com a variação de " + type_var
    fig_name = "Raio x " + type_var + ".png"

    plt.title(fig_title)
    plt.xlabel(type_var)
    plt.ylabel("raio [m]")
    plt.grid(True)
    plt.plot(variables, rad_list)
    plt.savefig(fig_name, format="png")
    plt.clf()

if __name__ == "__main__":

    inputs = {
        "bandwidth": 0,
        "area": 0,
        "ro": 0,
        "gama": 0,
        "eta": 0,
        "prob_outdoor": 0,
        "prob_indoor": 0,
        "prob_incar": 0,
        "figura": 0,
        "users": 0
    }

    #inputs
    inputs["bandwidth"] = float(input("Please select Bandwidth in MHz[Default:20]: ") or "20")
    inputs["area"] = float(input("Select desired area to be covered in km2[Default:5]: ") or "5")
    inputs["ro"] = float(input("Select neighbour cells signals correlation ro[Default:0.8]: ") or "0.8")
    inputs["gama"] = float(input("Select charge factor gama[Default:0.8]: ") or "0.8")
    inputs["eta"] = float(input("Select mean SNR degradation due to interference eta[Default:1.5]: ") or "1.5")
    inputs["prob_outdoor"] = float(input("Select probability of user being outdoors[Default:0.3]: ") or "0.3")
    inputs["prob_indoor"] = float(input("Select probability of user being indoors[Default:0.4]: ") or "0.4")
    inputs["prob_incar"] = float(input("Select probability of user being in a car[Default: 0.3]: ") or "0.3")
    inputs["figura"] = float(input("Select noise figure[Default:8]: ") or "8")
    inputs["users"] = int(input("Select amount of users in the area[Default:4000]: ") or "4000")

    bw = inputs["bandwidth"]
    area = inputs["area"]
    ro = inputs["ro"]
    gama = inputs["gama"]
    eta = inputs["eta"]
    prob_outdoor = inputs["prob_outdoor"]
    prob_indoor = inputs["prob_indoor"]
    prob_incar = inputs["prob_incar"]
    figura = inputs["figura"]
    users = inputs["users"]
    mcs = start_mcs()
    n = 4
    sigma = 7.6
    prob_cobertura_celula = 0.9

    mcs_copy = mcs

    mcs, m_in = outdoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs)
    mcs = indoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs)
    mcs = incar_radius(sigma, n, ro, gama, eta, bw, figura, mcs)

    rate_list = list()

    for modulation in mcs:
        radius = prob_outdoor*mcs[modulation][5] + prob_indoor*mcs[modulation][6] + prob_incar*mcs[modulation][7]
        radius = round(radius,4)
        print(modulation, "radius:", radius, "m")
        mcs[modulation].append(radius)
        
    cell_area, cell_quant = count_hex(area, mcs["QPSK 1/3"][8]*0.95)
    print("Cell area:", round(cell_area,2), "km²")
    print("Cell quantity:", cell_quant)

    active_users = active_users(users, area, prob_outdoor, prob_indoor, prob_incar, mcs, cell_area)
    active_users = required_throughput(active_users, service_mix) #active_users[modulation][1] is now required throughput
    check_sinr(active_users, mcs, m_in, gama)

    print("Generating graphs...")

    radius_gama = list()
    radius_eta = list()
    radius_ro = list()

    for i in range(1, 100):
        var01 = i/100
        var12 = (i/10)

        radius_gama.append(var_to_rad(sigma, n, ro, var01, eta, bw, figura, mcs_copy))
        radius_eta.append(var_to_rad(sigma, n, ro, gama, var12, bw, figura, mcs_copy))
        radius_ro.append(var_to_rad(sigma, n, var01, gama, eta, bw, figura, mcs_copy))
    
    print_graph(radius_gama, "gama")
    print_graph(radius_eta, "eta")
    print_graph(radius_ro, "ro")