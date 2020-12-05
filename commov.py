from math import sqrt, log10, log2, pow
from matplotlib import pyplot as plt

import argparse

#Parâmetros de referência: Link Budget de Uplink
uplink_pot_tx = 23 #dBm
uplink_Gtx = 0
uplink_tx_loss = 0
uplink_SNR_req = 0
uplink_sens_req_rx = -101.5 #dBm
uplink_Grx = 18 #dBi
uplink_rx_loss = 3 #dB
uplink_multipath = 3 #dB
uplink_shadow_margin = 4 #dB

def start_mcs():
    # 0 -> bit rate
    # 1 -> code rate

    mcs = {
        "QPSK 1/2":[2, 0.5],
        "QPSK 3/4":[2, 0.75],
        "16-QAM 1/2":[4, 0.5],
        "16-QAM 3/4":[4, 0.75],
        "64-QAM 1/2":[6, 0.5],
        "64-QAM 2/3":[6, 0.666666],
        "64-QAM 3/4":[6, 0.75],
        "64-QAM 5/6":[6, 0.833333]}

    return mcs

def count_hex(area, radius):

    tri_area = radius*(radius*sqrt(3)/2)/2
    hex_area = 6*tri_area
    total_hex = area/hex_area

    return hex_area, total_hex

def link_budget(pot_tx, pot_rx, snr, Md, Gt, Pt):
    Lmax = pot_tx - pot_rx + Gt - Pt - snr - Md
    return Lmax

#Raio máximo é determinado a partir do modelo de Okumura-Hata, resolvido para a distância ao utilizar-se a maior perda aceitável (Link Budget)
def max_radius(max_loss, freq, hb, hm):
    b = 44.9 - 6.55*log10(hb)
    a = (1.11*log10(freq) - 0.7)*hm - 1.56*log10(freq) - 0.8
    log_R_max = (max_loss - 69.55 - 26.16*log10(freq) + 13.82*log10(hb) + a)/b
    R_max = pow(10, log_R_max)
    R_max = R_max*pow(10, 3)

    return R_max

def max_throughput(bw, bpsimb, code_rate):
    t_simb = 0.000071367
    subcarries = bw*60
    throughput = (1/t_simb)*subcarries*bpsimb*code_rate/pow(10, 6)

    return throughput

def efficiency(throughput, bw):
    se = throughput/bw

    return se

def snr(bw, c):
    snr = pow(2, c/(bw)) - 1

    return snr

def sir_triple_setor(hb):
    q = sqrt(3)
    gama = (44.9 - 6.55*log10(hb))/10

    sir = pow(q, gama)/2

    return sir

def loss_model(dist, freq, hb, hm):
    b = 44.9 - 6.55*log10(hb)
    a = (1.11*log10(freq) - 0.7)*hm - 1.56*log10(freq) - 0.8
    loss = 69.55 + 26.16*log10(freq) - 13.82*log10(hb) + b*log10(dist) - a

    return loss

def snr_link_budget(pot_tx, pot_rx, loss, Md, Gt, Pt):
    snr = pot_tx - pot_rx + Gt - Pt - loss - Md

    return snr

def rate_Mbps(bw, snr):
    try:
        r_Mbps = bw*log2(1 + snr)
    except:
        r_Mbps = 0
        print(snr)

    return r_Mbps

def plot_rate(mcs):
    colors = ["red", "blue", "green", "pink", "magenta", "yellow", "black", "orange"]

    plt.title("Throughput da rede de acordo com a distância do UE")
    plt.ylabel("Throughput [Mbps]")
    plt.xlabel("Distância [m]")
    plt.grid(True)
    
    for n, key in enumerate(mcs):
        plt.plot(range(1, int(mcs[key][5])), mcs[key][6], color=colors[n], label=key)
        print(key, " SNR mínimo:", mcs[key][4])
    
    plt.legend()
    plt.savefig("Throughput x Distância.png", format="png")

def signal_rate(rate, signal):
    noise = signal/rate
    return noise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Input arguments")
    parser.add_argument('--bandwidth', '-b',
                        type=float,
                        help="Channel bandwidth to work on. It must be 1.4, 3, 5, 10, 15 or 20 [MHz]. Default = 20.",
                        default=20)
    parser.add_argument('--area', '-a',
                        type=float,
                        help="Area to cover with LTE [km]. Required argument, without default value.",
                        required=True)
    parser.add_argument('--throughput', '-t',
                        type=float,
                        help="Mean throughput desired for the LTE network [Mbps]. Default = 25.",
                        default=25)

    inputs = parser.parse_args()

    bw = inputs.bandwidth
    area = inputs.area
    desired_throughput = inputs.throughput

    mcs = start_mcs()
    chosen_modulation = str()

    for modulation in mcs:
        loss_rate = float()
        snr_rate = float()
        noise = float()
        interference = float()
        
        rate = list()
        sinr = list()

        bit_rate = mcs[modulation][0]
        code_rate = mcs[modulation][1]
        
        throughput = max_throughput(bw, bit_rate, code_rate)
        mcs[modulation].append(throughput)

        if desired_throughput >= throughput:
            chosen_modulation = modulation

        snr_min = snr(bw, throughput)
        mcs[modulation].append(snr_min)
        
        sir = sir_triple_setor(30)
        mcs[modulation].append(sir)
        interference = signal_rate(sir, uplink_pot_tx + uplink_Gtx)

        loss_max = link_budget(uplink_pot_tx, 
                                uplink_sens_req_rx,
                                snr_min,
                                uplink_shadow_margin,
                                uplink_Grx + uplink_Gtx + uplink_multipath,
                                uplink_rx_loss + uplink_tx_loss)
        radius = max_radius(loss_max, 900, 30, 1.5)
        mcs[modulation].append(radius)

        for r in range(1, int(radius)):
            loss_rate = loss_model(r/1000, 900, 30, 1.5)
            snr_rate = snr_link_budget(uplink_pot_tx, 
                                        uplink_sens_req_rx,
                                        loss_rate,
                                        uplink_shadow_margin,
                                        uplink_Grx + uplink_Gtx + uplink_multipath,
                                        uplink_rx_loss + uplink_tx_loss)
            noise = signal_rate(snr_rate, uplink_pot_tx + uplink_Gtx)
            
            rate.append(rate_Mbps(bw, snr_rate))
            sinr.append((uplink_pot_tx + uplink_Gtx)/(interference + noise))

        mcs[modulation].append(rate)
    
    enb_area, enb_quantity = count_hex(area, mcs[chosen_modulation][5])

    plot_rate(mcs)

    #print(mcs)