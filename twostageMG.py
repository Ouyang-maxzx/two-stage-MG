# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 19:09:18 2023

@author: wyx
"""
# 导入
from gurobipy import *
import numpy as np

# 定义常量
dT = 1
PG_max = 800  # DG出力上限
PG_min = 80  # DG出力下限
a = 0.67  # 成本系数
b = 0  # 成本系数
K_S = 0.38  # 单位充放电成本
PS_max = 500  # 充放电功率上限
ES_max = 1800  # 荷电状态上限
ES_min = 400  # 荷电状态下限
ES0 = 1000  # 初始荷电状态
eta = 0.95  # 储能充放电效率
K_DR = 0.32  # 单位调度成本
D_DR = 2940  # 总用电需求
D_DR_max = 200  # 最大用电需求
D_DR_min = 50  # 最小用电需求
P_DR0=[80,70,60,50,70,70,90,100,120,150,170,200,140,100,100,120,140,150,190,200,200,190,100,80]
PM_max = 1500  # 最大交互功率
lambda1 = [0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.9, 1.35, 1.35, 1.35,
           0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1.35, 1.35, 1.35, 1.35, 1.35, 0.48]

# 不确定
P_PV0=np.array([0,0,0,0,0,0.0465,0.1466,0.3135,0.4756,0.5213,0.6563,1,0.7422,0.6817,0.4972,0.4629,
       0.2808,0.0948,0.0109,0,0,0,0,0])*1500
P_L0=np.array([0.4658,0.4601,0.5574,0.5325,0.5744,0.6061,0.6106,0.6636,0.741,0.708,0.7598,0.8766,0.7646,
    0.7511,0.6721,0.5869,0.6159,0.6378,0.6142,0.6752,0.6397,0.5974,0.5432,0.4803])*1200
delta_u_PV=0.15
delta_u_L=0.1
tau_L = 12
tau_PV = 6
# 定义变量
T_set = np.arange(24)
model = Model('MGs')
# 第一阶段
Us_t = model.addVars(T_set, vtype=GRB.BINARY, name='Us_t')
Um_t = model.addVars(T_set, vtype=GRB.BINARY, name='Um_t')
# 第二阶段
Pg_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Pg_t')
Ps_ch_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Ps_ch_t')
Ps_dis_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Ps_dis_t')
Pdr_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Pdr_t')
Pdr1_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Pdr1_t')
Pdr2_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Pdr2_t')
Pbuy_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Pbuy_t')
Psell_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Psell_t')
Ppv_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Ppv_t')
Pl_t = model.addVars(T_set,lb=-GRB.INFINITY, name='Pl_t')
# 不确定
B_pv_t = model.addVars(T_set,vtype=GRB.BINARY,name='B_pv_t')
B_l_t = model.addVars(T_set,vtype=GRB.BINARY,name='B_l_t')

model.addConstrs((Pg_t[t]>=PG_min for t in T_set),name='min-Pg_t')
model.addConstrs((Pg_t[t]<=PG_max for t in T_set),name='max-Pg_t')
model.addConstrs((Pdr_t[t]<=D_DR_max for t in T_set),name='max-Pdr_t')
model.addConstrs((Pdr_t[t]>=D_DR_min for t in T_set),name='min-Pdr_t')
model.addConstrs((Pdr1_t[t]>=0 for t in T_set),name='Pdr1_t')
model.addConstrs((Pdr2_t[t]>=0 for t in T_set),name='Pdr2_t')
model.addConstrs((Pbuy_t[t]>=0 for t in T_set),name='Pbuy_t')
model.addConstrs((Psell_t[t]>=0 for t in T_set),name='Psell_t')
model.addConstrs((Ppv_t[t]>=0 for t in T_set),name='Ppv_t')
model.addConstrs((Pl_t[t]>=0 for t in T_set),name='Pl_t')
model.addConstrs((Ps_dis_t[t]>=0 for t in T_set),name='Ppv_t')
model.addConstrs((Ps_ch_t[t]>=0 for t in T_set),name='Pl_t')


# 添加约束
# 可控分布式电源
obj1 = quicksum((a*Pg_t[t]+b) for t in T_set)
# 储能
obj2 = quicksum(K_S*(Ps_dis_t[t]/eta+Ps_ch_t[t]*eta) for t in T_set)
model.addConstrs(
    (Ps_dis_t[t] <= Us_t[t]*PS_max for t in T_set), name='Ps_disminmax')
model.addConstrs((Ps_ch_t[t] <= (1-Us_t[t]) *
                 PS_max for t in T_set), name='Ps_disminmax')

model.addConstr(eta*quicksum(Ps_ch_t[t] for t in T_set) -
                1/eta*quicksum(Ps_dis_t[t] for t in T_set) == 0, name='capacity')

model.addConstrs((ES0+eta*quicksum(Ps_ch_t[t]*dT for t in range(0, t_hat))-1/eta*quicksum(
    Ps_dis_t[t]*dT for t in range(0, t_hat)) >= ES_min for t_hat in T_set), name='SOCmin')
model.addConstrs((ES0+eta*quicksum(Ps_ch_t[t]*dT for t in range(0, t_hat))-1/eta*quicksum(
    Ps_dis_t[t]*dT for t in range(0, t_hat)) <= ES_max for t_hat in T_set), name='SOCmax')
# 需求响应负荷
obj3 = quicksum(K_DR*(Pdr1_t[t]+Pdr2_t[t]) for t in T_set)

model.addConstr(quicksum(Pdr_t[t] for t in T_set)==D_DR)
model.addConstrs((Pdr_t[t]-P_DR0[t]+Pdr1_t[t]-Pdr2_t[t]==0 for t in T_set),name='PDRauc')

# 电网交互功率
obj4 = quicksum(lambda1[t]*(Pbuy_t[t]-Psell_t[t]) for t in T_set)
model.addConstrs((Pbuy_t[t]-Psell_t[t]==Ps_ch_t[t]+Pdr_t[t]+Pl_t[t]-Pg_t[t]-Ps_dis_t[t]-Ppv_t[t] for t in T_set),name='Pbalance')
# Importance!!! 因为Pl和Ppv为不确定参数，所以相当于变量
model.addConstrs((Pbuy_t[t]<=Um_t[t]*PM_max for t in T_set),name='PMbuy')
model.addConstrs((Psell_t[t]<=(1-Um_t[t])*PM_max for t in T_set),name='PMsell')

# 不确定参数
model.addConstrs((Ppv_t[t]==P_PV0[t]-B_pv_t[t]*delta_u_PV for t in T_set),name='uncertainty-1')
model.addConstrs((Pl_t[t]==P_L0[t]+B_l_t[t]*delta_u_L for t in T_set),name='uncertainty-2')
model.addConstr(quicksum(B_pv_t[t]for t in T_set)<=tau_PV,name='B-1')
model.addConstr(quicksum(B_l_t[t]for t in T_set)<=tau_L,name='B-2')

# 设置目标函数
model.setObjective(obj1+obj2+obj3+obj4,GRB.MINIMIZE)

# 优化
model.optimize()
# model.update()


