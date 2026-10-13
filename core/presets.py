"""Presets das listas da 2ª unidade.

Os presets apenas preenchem os campos normais da interface; todos os valores
continuam editáveis antes da execução.
"""
from .specs import DesignSpecs


EXERCISES={
    "Lista U2.1 — Questão 1 — PD (Mp e ts 5%)":{
        "list":"U2.1","question":"1",
        "controller":"PD",
        "num_g":[4,16],"den_g":[1,4,4,0],
        "num_h":[1],"den_h":[1],
        "specs":lambda:DesignSpecs.from_mp_ts(10,4,0.05),
        "description":"G(s)=4(s+4)/(s³+4s²+4s), H(s)=1; Mp≤10%, ts(5%)<4 s.",
    },
    "Lista U2.1 — Questão 2 — PD (ξ e ωn)":{
        "list":"U2.1","question":"2",
        "controller":"PD",
        "num_g":[1],"den_g":[10000,0,-11772],
        "num_h":[1],"den_h":[1],
        "specs":lambda:DesignSpecs.from_zeta_wn(0.7,0.5),
        "description":"G(s)=1/[10000(s²−1,1772)], H(s)=1; ξ=0,7, ωn=0,5 rad/s.",
    },
    "Lista U2.1 — Questão 3 — PI (polos desejados)":{
        "list":"U2.1","question":"3",
        "controller":"PI",
        "num_g":[5,25,20],"den_g":[1,4,4],
        "num_h":[0.2],"den_h":[1,1],
        "specs":lambda:DesignSpecs.from_pole(-4,4),
        "description":"G(s)=5(s²+5s+4)/(s²+4s+4), H(s)=0,2/(s+1); polos −4±4j.",
    },
    "Lista U2.1 — Questão 4 — PID (zeros iguais)":{
        "list":"U2.1","question":"4",
        "controller":"PID",
        "controller_config":{"pid_zero_mode":"equal"},
        "num_g":[5],"den_g":[1,12,22,20],
        "num_h":[0.4],"den_h":[1],
        "specs":lambda:DesignSpecs.from_mp_ts(20,5,0.02),
        "description":"G(s)=5/(s³+12s²+22s+20), H(s)=0,4; zeros iguais, Mp≤20%, ts(2%)<5 s.",
    },

    "Lista U2.2 — Questão 1 — PID + Tustin":{
        "list":"U2.2","question":"1",
        "controller":"PID",
        "controller_config":{"pid_zero_mode":"equal"},
        "num_g":[5,15],"den_g":[1,4,0],
        "num_h":[1],"den_h":[1,1],
        "specs":lambda:DesignSpecs.from_mp_ts(10,3,0.05),
        "discretize_after_design":True,
        "discretization_target":"Controlador Gc(s)",
        "discretization_method":"Tustin",
        "sample_time":2.0,
        "description":(
            "G(s)=5(s+3)/[s(s+4)], H(s)=1/(s+1). "
            "Projetar PID com z1=z2 para Mp≤10% e ts(5%)<3 s; "
            "depois discretizar o controlador por Tustin com T=2 s."
        ),
    },
    "Lista U2.2 — Questão 2 — Genérico a/(s+b) + Euler":{
        "list":"U2.2","question":"2",
        "controller":"Genérico",
        "controller_config":{
            "fixed_zeros":[],
            "fixed_poles":[],
            "unknown_kind":"pole",
            "unknown_count":1,
            "require_lhp":True,
        },
        "num_g":[2,2],"den_g":[1,2,2],
        "num_h":[1,3],"den_h":[1,5],
        "specs":lambda:DesignSpecs.from_pole(-2.5,2.0),
        "discretize_after_design":True,
        "discretization_target":"G(s)H(s)Gc(s)",
        "discretization_method":"Euler / Forward",
        "sample_time":1.0,
        "description":(
            "G(s)=2(s+1)/(s²+2s+2), H(s)=(s+3)/(s+5), "
            "Gc(s)=a/(s+b). Determinar a,b para polos −2,5±2j; "
            "depois discretizar G(s)H(s)Gc(s) por Euler com T=1 s."
        ),
    },
}


# Aliases antigos para não quebrar testes/código que usava os nomes da Lista U2.1.
LEGACY_NAMES={
    "Questão 1 — PD (Mp e ts 5%)":"Lista U2.1 — Questão 1 — PD (Mp e ts 5%)",
    "Questão 2 — PD (ξ e ωn)":"Lista U2.1 — Questão 2 — PD (ξ e ωn)",
    "Questão 3 — PI (polos desejados)":"Lista U2.1 — Questão 3 — PI (polos desejados)",
    "Questão 4 — PID (zeros iguais)":"Lista U2.1 — Questão 4 — PID (zeros iguais)",
}


def get_exercise(name):
    return EXERCISES[LEGACY_NAMES.get(name,name)]
