"""Presets da 1ª lista da 2ª unidade."""
from .specs import DesignSpecs

EXERCISES={
"Questão 1 — PD (Mp e ts 5%)":{
"controller":"PD","num_g":[4,16],"den_g":[1,4,4,0],"num_h":[1],"den_h":[1],
"specs":lambda:DesignSpecs.from_mp_ts(10,4,0.05),
"description":"G(s)=4(s+4)/(s³+4s²+4s), H(s)=1; Mp≤10%, ts(5%)<4 s."},
"Questão 2 — PD (ξ e ωn)":{
"controller":"PD","num_g":[1],"den_g":[10000,0,-11772],"num_h":[1],"den_h":[1],
"specs":lambda:DesignSpecs.from_zeta_wn(0.7,0.5),
"description":"G(s)=1/[10000(s²−1,1772)], H(s)=1; ξ=0,7, ωn=0,5 rad/s."},
"Questão 3 — PI (polos desejados)":{
"controller":"PI","num_g":[5,25,20],"den_g":[1,4,4],"num_h":[0.2],"den_h":[1,1],
"specs":lambda:DesignSpecs.from_pole(-4,4),
"description":"G(s)=5(s²+5s+4)/(s²+4s+4), H(s)=0,2/(s+1); polos −4±4j."},
"Questão 4 — PID (zeros iguais)":{
"controller":"PID","num_g":[5],"den_g":[1,12,22,20],"num_h":[0.4],"den_h":[1],
"specs":lambda:DesignSpecs.from_mp_ts(20,5,0.02),
"description":"G(s)=5/(s³+12s²+22s+20), H(s)=0,4; zeros iguais, Mp≤20%, ts(2%)<5 s."}
}
