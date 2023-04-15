# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 19:43:04 2023

@author: wyx
"""
from twostageMG import *

A = model.getA().toarray()  #
b = np.array(model.RHS)
sense = np.array(model.sense)
f = np.array(model.obj)

# 2. 加工系数矩阵，将A分为不等式约束Aineq和等式约束Aeq两块
Aeq = A[sense == '=', :]  # Submatrix corresponding to equality constraints
beq = b[sense == '=']  # RHS for equality constraints

# Submatrix corresponding to less-or-equal constraints
Ale = A[sense == '<', :]
ble = b[sense == '<']  # RHS for LE constraints

# Submatrix corresponding to greater-or-equal constraints
Age = A[sense == '>', :]
bge = b[sense == '>']  # RHS for GE constraints

Aineq = np.vstack((-Ale, Age))  # 把所有的<=和>=组合在一起
bineq = np.append(-ble, bge)  # 这里用append使bineq为一个一维矩阵，而不是2行1列的二维矩阵，避免后面的运行错误
# 提取第一阶段变量
num_1 = len(Us_t)+len(Um_t)
num_2 = len(Pg_t)+len(Ps_ch_t)+len(Ps_dis_t)+len(Pdr_t)+len(Pdr1_t)+len(Pdr2_t)+len(Pbuy_t)+len(Psell_t)+len(Ppv_t)+len(Pl_t)
num_3 = len(B_pv_t)+len(B_l_t)

G, M,  E, h = [], [], [], []
G1, M1,  E1, h1 = [], [], [], []
for row in np.arange(Aineq.shape[0]):
    # 第一阶段变量
    if not np.all(Aineq[row][num_1:] == 0): #第一阶段独有约束
        E.append(Aineq[row][:num_1]) #x
        G.append(Aineq[row][num_1:num_1+num_2]) #y
        M.append(Aineq[row][num_1+num_2:]) #u
        h.append(bineq[row]) #b
for row in np.arange(Aeq.shape[0]):
    # 第一阶段变量
    if not np.all(Aeq[row][num_1:] == 0): #第一阶段独有约束
        E1.append(Aeq[row][:num_1]) #x
        G1.append(Aeq[row][num_1:num_1+num_2]) #y
        M1.append(Aeq[row][num_1+num_2:]) #u
        h1.append(beq[row]) #b
       
G, M,  E, h = np.array(G),np.array(M),np.array(E),np.array(h)
G1, M1,  E1, h1 = np.array(G1),np.array(M1),np.array(E1),np.array(h1)
PVL = np.append(P_PV0, P_L0)
# 建立矩阵模型
m = Model('specification')
x = m.addMVar((48,),vtype=GRB.BINARY,name='x')
y = m.addMVar((240,),lb=-GRB.INFINITY,name='y')
u = m.addMVar((48,),vtype=GRB.BINARY,name='u')
f = np.array(model.obj)
c = f[num_1:num_1+num_2]
# 构建确定性模型
m.setObjective(c@y,GRB.MINIMIZE)
m.addConstr(G@y>=h-E@x-M@u)
m.addConstr(G1@y==h1-E1@x-M1@u)
m.optimize()

if abs(m.objVal-model.objVal)<=0.000001:
    print("矩阵形式的模型与元素形式的模型一致")
else:
    print("模型构建不一致，需要检查")







